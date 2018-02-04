/****************************************************************************
 *
 *   Copyright (c) 2015 Estimation and Control Library (ECL). All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name ECL nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file covariance.cpp
 * Contains functions for initialising, predicting and updating the state
 * covariance matrix
 *
 * @author Roman Bast <bastroman@gmail.com>
 *
 */

#include "../ecl.h"
#include "ekf.h"
#include <math.h>
#include "mathlib.h"

void Ekf::initialiseCovariance()
{
	memset(&P_UKF, 0, sizeof(P_UKF));
	memset(&Q_UKF, 0, sizeof(Q_UKF));
	memset(&SPA_UKF, 0, sizeof(SPA_UKF));

	// calculate average prediction time step in sec
	float dt = 0.001f * (float)FILTER_UPDATE_PERIOD_MS;

	// define the initial angle uncertainty as variances for a rotation vector
	P_UKF(2,2) = P_UKF(1,1) = P_UKF(0,0) = sq(_params.initial_tilt_err);

	// velocity
	P_UKF(4,4) = P_UKF(3,3) = sq(fmaxf(_params.gps_vel_noise, 0.01f));
	P_UKF(5,5) = sq(1.5f) * P_UKF(4,4);

	// position
	P_UKF(7,7) = P_UKF(6,6) = sq(fmaxf(_params.gps_pos_noise, 0.01f));
	if (_control_status.flags.rng_hgt) {
		P_UKF(8,8) = sq(fmaxf(_params.range_noise, 0.01f));
	} else if (_control_status.flags.gps_hgt) {
		float lower_limit = fmaxf(_params.gps_pos_noise, 0.01f);
		float upper_limit = fmaxf(_params.pos_noaid_noise, lower_limit);
		P_UKF(8,8) = sq(1.5f * math::constrain(_gps_sample_delayed.vacc, lower_limit, upper_limit));
	} else {
		P_UKF(8,8) = sq(fmaxf(_params.baro_noise, 0.01f));
	}

	// gyro bias
	P_UKF(11,11) = P_UKF(10,10) = P_UKF(9,9) = sq(_params.switch_on_gyro_bias * dt);

	// accel bias
	_prev_dvel_bias_var(0) = P_UKF(12,12) = sq(_params.switch_on_accel_bias * dt);
	_prev_dvel_bias_var(1) = P_UKF(13,13) = P_UKF(12,12);
	_prev_dvel_bias_var(2) = P_UKF(14,14) = P_UKF(12,12);

	// record IMU bias state covariance reset time - used to prevent resets being performed too often
	_last_imu_bias_cov_reset_us = _imu_sample_delayed.time_us;

	// variances for optional states

	// earth frame and body frame magnetic field
	// set to observation variance
	for (uint8_t index = 15; index <= 20; index ++) {
		P_UKF(index,index) = sq(_params.mag_noise);
	}

	// wind
	P_UKF(22,22) = P_UKF(21,21) = sq(_params.initial_wind_uncertainty);

	// set the control input noise variances
	Q_UKF(2,2) = Q_UKF(1,1) = Q_UKF(0,0) = sq(_params.gyro_noise * dt);
	Q_UKF(5,5) = Q_UKF(4,4) = Q_UKF(3,3) = sq(_params.accel_noise * dt);

	// calculate cholesky decomposition matrices
	SP_UKF = matrix::cholesky(P_UKF);
	SQ_UKF = matrix::cholesky(Q_UKF);

}

void Ekf::get_pos_var(Vector3f &pos_var)
{
	pos_var(0) = P_UKF(6,6);
	pos_var(1) = P_UKF(7,7);
	pos_var(2) = P_UKF(8,8);
}

void Ekf::get_vel_var(Vector3f &vel_var)
{
	vel_var(0) = P_UKF(3,3);
	vel_var(1) = P_UKF(4,4);
	vel_var(2) = P_UKF(5,5);
}

void Ekf::prediction()
{

	float dt = math::constrain(_imu_sample_delayed.delta_ang_dt, 0.0005f * FILTER_UPDATE_PERIOD_MS, 0.002f * FILTER_UPDATE_PERIOD_MS);

	// compute noise variance for stationary processes
	float process_noise_variance[14] = {};

	// convert rate of change of rate gyro bias (rad/s**2) as specified by the parameter to an expected change in delta angle (rad) since the last update
	float d_ang_bias_sig = dt * dt * math::constrain(_params.gyro_bias_p_noise, 0.0f, 1.0f);

	// convert rate of change of accelerometer bias (m/s**3) as specified by the parameter to an expected change in delta velocity (m/s) since the last update
	float d_vel_bias_sig = dt * dt * math::constrain(_params.accel_bias_p_noise, 0.0f, 1.0f);

	// inhibit learning of imu acccel bias if the manoeuvre levels are too high to protect against the effect of sensor nonlinearities or bad accel data is detected
	float alpha = 1.0f - math::constrain((dt / _params.acc_bias_learn_tc), 0.0f, 1.0f);
	_ang_rate_mag_filt = fmaxf(_imu_sample_delayed.delta_ang.norm(), alpha * _ang_rate_mag_filt);
	_accel_mag_filt = fmaxf(_imu_sample_delayed.delta_vel.norm(), alpha * _accel_mag_filt);
	if (_ang_rate_mag_filt > dt * _params.acc_bias_learn_gyr_lim
			|| _accel_mag_filt > dt * _params.acc_bias_learn_acc_lim
			|| _bad_vert_accel_detected) {
		// store the bias state variances to be reinstated later
		if (!_accel_bias_inhibit) {
			_prev_dvel_bias_var(0) = P_UKF(12,12);
			_prev_dvel_bias_var(1) = P_UKF(13,13);
			_prev_dvel_bias_var(2) = P_UKF(14,14);
		}
		_accel_bias_inhibit = true;
	} else {
		if (_accel_bias_inhibit) {
			// reinstate the bias state variances
			P_UKF(12,12) = _prev_dvel_bias_var(0);
			P_UKF(13,13) = _prev_dvel_bias_var(1);
			P_UKF(14,14) = _prev_dvel_bias_var(2);
		} else {
			// store the bias state variances to be reinstated later
			_prev_dvel_bias_var(0) = P_UKF(12,12);
			_prev_dvel_bias_var(1) = P_UKF(13,13);
			_prev_dvel_bias_var(2) = P_UKF(14,14);
		}
		_accel_bias_inhibit = false;
	}

	// Don't continue to grow the earth field variances if they are becoming too large or we are not doing 3-axis fusion as this can make the covariance matrix badly conditioned
	float mag_I_sig;
	if (_control_status.flags.mag_3D && (P_UKF(15,15) + P_UKF(16,16) + P_UKF(17,17)) < 0.1f) {
		mag_I_sig = dt * math::constrain(_params.mage_p_noise, 0.0f, 1.0f);

	} else {
		mag_I_sig = 0.0f;
	}

	// Don't continue to grow the body field variances if they is becoming too large or we are not doing 3-axis fusion as this can make the covariance matrix badly conditioned
	float mag_B_sig;
	if (_control_status.flags.mag_3D && (P_UKF(18,18) + P_UKF(19,19) + P_UKF(20,20)) < 0.1f) {
		mag_B_sig = dt * math::constrain(_params.magb_p_noise, 0.0f, 1.0f);

	} else {
		mag_B_sig = 0.0f;
	}

	float wind_vel_sig;

	// Don't continue to grow wind velocity state variances if they are becoming too large or we are not using wind velocity states as this can make the covariance matrix badly conditioned
	if (_control_status.flags.wind && (P_UKF(21,21) + P_UKF(22,22)) < 2.0f * sq(_params.initial_wind_uncertainty)) {
		wind_vel_sig = dt * math::constrain(_params.wind_vel_p_noise, 0.0f, 1.0f);

	} else {
		wind_vel_sig = 0.0f;
	}

	// Construct the process noise variance diagonal for those states with a stationary process model
	// These are kinematic states and their error growth is controlled separately by the IMU noise variances
	// delta angle bias states
	process_noise_variance[2] = process_noise_variance[1] = process_noise_variance[0] = sq(d_ang_bias_sig);
	// delta_velocity bias states
	process_noise_variance[5] = process_noise_variance[4] = process_noise_variance[3] = sq(d_vel_bias_sig);
	// earth frame magnetic field states
	process_noise_variance[8] = process_noise_variance[7] = process_noise_variance[6] = sq(mag_I_sig);
	// body frame magnetic field states
	process_noise_variance[11] = process_noise_variance[10] = process_noise_variance[9] = sq(mag_B_sig);
	// wind velocity states
	process_noise_variance[13] = process_noise_variance[12] = sq(wind_vel_sig);

	// Add variances to covariance matrix
	for (unsigned i = 0; i < 14; i++) {
		P_UKF(i,i) += process_noise_variance[i];
	}

	// calculate an array of sigma points for the augmented state vector
	CalcSigmaPoints();

	// convert the attitude error vector sigma points to equivalent delta quaternions
	// Don't calculate the first column - it is zero by definition and not used
	Quatf dq_sigma[UKF_N_SIGMA];
	float normsigmaX2;
	for (uint8_t s=1; s<(2*_ukf_L); s++) {
		normsigmaX2 = sq(_sigma_x_a(0,s)) + sq(_sigma_x_a(1,s)) + sq(_sigma_x_a(2,s));
		dq_sigma[s](0) = (-_grp_a * normsigmaX2 + _grp_f*sqrtf( sq(_grp_f) + (1.0f - sq(_grp_a)) * normsigmaX2)) / (sq(_grp_f) + normsigmaX2);
		for (uint8_t i=0; i<3; i++) {
			dq_sigma[s](i+1) = (_grp_a + dq_sigma[s](0)) * _sigma_x_a(i,s) / _grp_f;
		}
	}

	// Apply the delta quaternions to the previous estimate to calculate the
	// quaternion sigma points. Although we could propagate the delta angles
	// through the vehicle state prediction, it is more accurate to use
	// quaternions and convert back to a set of GRP attitude error vectors
	// when the covariance information needs to be extracted.
	 _sigma_quat[0] = _ukf_states.data.quat;
	 for (uint8_t s=1; s<(2*_ukf_L); s++) {
	     _sigma_quat[s] = dq_sigma[s] * _sigma_quat[0];
	 }

	// Propagate each sigma point forward using the INS equations
	predictSigmaPoints();

	// Use the sigma points to calculate the mean state vector and the covariances

	// The zero perturbation sigma points which represent the nominal expected value are propagated as the state vector
	for (uint8_t i=0; i<UKF_N_AUG_STATES; i++) {
		_ukf_states.vector(i) = _sigma_x_a(i,0);
	}
	_ukf_states.data.quat = _sigma_quat[0];

	// Convert propagated quaternions to delta quaternions around the
	// nominal expected value
	Quatf sigma_dq[UKF_N_SIGMA];
	Quatf qm_inverse = _sigma_quat[0].inversed();
	sigma_dq[0](0) = 1.0f;
	sigma_dq[0](3) = sigma_dq[0](2) = sigma_dq[0](1) = 0.0f;
	for (uint8_t s=1; s<UKF_N_SIGMA; s++) {
	    sigma_dq[s] = _sigma_quat[s] * qm_inverse;
	}

	// Convert error quaternions to attitude error vector
	// By definition the first column is the expected value
	_sigma_x_a(0,0) = _sigma_x_a(0,1) = _sigma_x_a(0,2) = 0.0f;
	for (uint8_t s=1; s<UKF_N_SIGMA; s++) {
		for (uint8_t i=0; i<3; i++) {
			_sigma_x_a(i,s) = _grp_f * sigma_dq[s](i+1)/(_grp_a + sigma_dq[s](0));
		}
	}

	// Calculate mean of predicted states from the sigma points
	for (uint8_t i=0; i<UKF_N_STATES; i++) {
		_ukf_states_mean.vector(i) = 0.0f;
		for (uint8_t s=0; s<UKF_N_SIGMA; s++) {
			_ukf_states_mean.vector(i) += _sigma_x_a(i,s) * _ukf_wm[s];
		}
	}

	// Calculate covariance of predicted vehicle states from the sigma points using
	// the unscented transform
	memset(&P_UKF, 0, sizeof(P_UKF));
	for (uint8_t s=0; s<UKF_N_SIGMA; s++) {
		// P = P + param.ukf.wc(i)*(sigma_x_a(1:param.ukf.nP,s) - x_m)*(sigma_x_a(1:param.ukf.nP,s) - x_m)';
		float state_delta[UKF_N_STATES];
		for (uint8_t i=0; i<UKF_N_STATES; i++) {
			state_delta[i] = _sigma_x_a(i,s) - _ukf_states_mean.vector(i);
		}
		for (uint8_t i=0; i<UKF_N_STATES; i++) {
			for (uint8_t j=0; j<UKF_N_STATES; j++) {
				P_UKF(i,j) = _ukf_wc[s] * state_delta[i] * state_delta[j];
			}
		}
	}

	// fix gross errors in the covariance matrix
	fixCovarianceErrors();

}

void Ekf::fixCovarianceErrors()
{
	// NOTE: This limiting is a last resort and should not be relied on
	// TODO: Split covariance prediction into separate F*P*transpose(F) and Q contributions
	// and set corresponding entries in Q to zero when states exceed 50% of the limit
	// Covariance diagonal limits. Use same values for states which
	// belong to the same group (e.g. vel_x, vel_y, vel_z)
	float P_lim[8] = {};
	P_lim[0] = 0.5f;		// attitude max var
	P_lim[1] = 1e6f;		// velocity max var
	P_lim[2] = 1e6f;		// positiion max var
	P_lim[3] = 1.0f;		// gyro bias max var
	P_lim[4] = 1.0f;		// delta velocity z bias max var
	P_lim[5] = 1.0f;		// earth mag field max var
	P_lim[6] = 1.0f;		// body mag field max var
	P_lim[7] = 1e6f;		// wind max var

	for (int i = 0; i <= 2; i++) {
		// attitude vector states
		P_UKF(i,i) = math::constrain(P_UKF(i,i), 1e-9f, P_lim[0]);
	}
	for (int i = 3; i <= 5; i++) {
		// NED velocity states
		P_UKF(i,i) = math::constrain(P_UKF(i,i), 1e-9f, P_lim[1]);
	}
	for (int i = 6; i <= 8; i++) {
		// NED position states
		P_UKF(i,i) = math::constrain(P_UKF(i,i), 1e-9f, P_lim[2]);
	}
	for (int i = 9; i <= 11; i++) {
		P_UKF(i,i) = math::constrain(P_UKF(i,i), 1e-9f, P_lim[3]);
	}
	for (int i = 12; i <= 14; i++) {
		P_UKF(i,i) = math::constrain(P_UKF(i,i), 1e-9f, P_lim[4]);
	}
	for (int i = 15; i <= 17; i++) {
		P_UKF(i,i) = math::constrain(P_UKF(i,i), 1e-9f, P_lim[5]);
	}
	for (int i = 18; i <= 20; i++) {
		P_UKF(i,i) = math::constrain(P_UKF(i,i), 1e-9f, P_lim[6]);
	}
	for (int i = 21; i <= 22; i++) {
		P_UKF(i,i) = math::constrain(P_UKF(i,i), 1e-9f, P_lim[7]);
	}
}

void Ekf::resetMagCovariance()
{
	// set the attitude covariance terms to zero
	zeroCovMat(0,2);

	// set the magnetic field covariance terms to zero
	zeroCovMat(15,20);

	// set the field state variance to the observation variance
	for (uint8_t rc_index=15; rc_index <= 20; rc_index ++) {
		P_UKF(rc_index,rc_index) = sq(_params.mag_noise);
	}
}

void Ekf::resetWindCovariance()
{
	// set the wind  covariance terms to zero
	zeroCovMat(21,22);

	if (_tas_data_ready && (_imu_sample_delayed.time_us - _airspeed_sample_delayed.time_us < (uint64_t)5e5)) {
		// Use airspeed and zer sideslip assumption to set initial covariance values for wind states

		// calculate the wind speed and bearing
		float spd = sqrtf(sq(_ukf_states.data.wind_vel(0))+sq(_ukf_states.data.wind_vel(1)));
		float yaw = atan2f(_ukf_states.data.wind_vel(1),_ukf_states.data.wind_vel(0));

		// calculate the uncertainty in wind speed and direction using the uncertainty in airspeed and sideslip angle
		// used to calculate the initial wind speed
		float R_spd = sq(math::constrain(_params.eas_noise, 0.5f, 5.0f) * math::constrain(_airspeed_sample_delayed.eas2tas, 0.9f, 10.0f));
		float R_yaw = sq(0.1745f);

		// calculate the variance and covariance terms for the wind states
		float cos_yaw = cosf(yaw);
		float sin_yaw = sinf(yaw);
		float cos_yaw_2 = sq(cos_yaw);
		float sin_yaw_2 = sq(sin_yaw);
		float sin_cos_yaw = sin_yaw*cos_yaw;
		float spd_2 = sq(spd);
		P_UKF(21,21) = R_yaw*spd_2*sin_yaw_2 + R_spd*cos_yaw_2;
		P_UKF(21,22) = - R_yaw*sin_cos_yaw*spd_2 + R_spd*sin_cos_yaw;
		P_UKF(22,21) = P_UKF(21,22);
		P_UKF(22,22) = R_yaw*spd_2*cos_yaw_2 + R_spd*sin_yaw_2;

		// Now add the variance due to uncertainty in vehicle velocity that was used to calculate the initial wind speed
		P_UKF(21,21) += P_UKF(3,3);
		P_UKF(22,22) += P_UKF(4,4);

	} else {
		// without airspeed, start with a small initial uncertainty to improve the initial estimate
		P_UKF(21,21) = sq(1.0f);
		P_UKF(22,22) = sq(1.0f);

	}
}

void Ekf::CalcSigmaPoints()
{
	// Calculate the lower diagonal Cholesky decomposition for the vehicle
	// state covariance matrix. This requires UKF_N_STATES^3 operations
	SP_UKF = matrix::cholesky(P_UKF);

	// Assemble the augmented covariance matrix
	for (int i = 0; i < UKF_N_STATES; i++) {
		for (int j = 0; j < UKF_N_STATES; j++) {
			SPA_UKF(i,j) = SP_UKF(i,j);
		}
	}
	for (int i = 0; i < UKF_N_Q; i++) {
		for (int j = 0; j < UKF_N_Q; j++) {
			SPA_UKF(i+UKF_N_STATES,j+UKF_N_STATES) = SQ_UKF(i,j);
		}
	}

	// Expected value of augmented state vector from previous frame
	matrix::Vector<float, UKF_N_AUG_STATES> x_a_prev;
	for (int i = 0; i < UKF_N_STATES; i++) {
		// vehicle states
		x_a_prev(i) = _ukf_states.vector(i);
	}
	for (int i = UKF_N_STATES; i < UKF_N_AUG_STATES; i++) {
		// IMU noise has zero mean expected value
		x_a_prev(i) = 0.0f;
	}

	// Generate sigma points for the augmented state vector

	// first column represents expected value (zero delta)
	for (int i = 0; i < UKF_N_AUG_STATES; i++) {
		// vehicle states
		_sigma_x_a(i,0) = x_a_prev(i);
	}

	// remaining columns
	float temp_var1 = sqrtf((float)_ukf_L + _ukf_lambda);
	for (int j = 0; j < UKF_N_AUG_STATES; j++) {
		for (int i = 0; i < UKF_N_AUG_STATES; i++) {
			float temp_var2 = temp_var1 * SPA_UKF(i,j);
			_sigma_x_a(i,j+1) = x_a_prev(i) + temp_var2;
			_sigma_x_a(i,j+1+UKF_N_AUG_STATES) = x_a_prev(i) - temp_var2;
		}
	}

	_sigma_points_are_stale = false;
}
