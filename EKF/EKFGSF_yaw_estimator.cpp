/****************************************************************************
 *
 *   Copyright (c) 2020 Estimation and Control Library (ECL). All rights reserved.
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

/*
 * @file EKFGSF_yaw_estimator.cpp
 *
 * Functions used to provide a backup estimate of yaw angle that doesn't use a
 * magnetometer. A bank of AHRS solutions, each starting from a different initial
 * yaw hoypothesis are used to provide an atitude reference using only IMU and
 * optional true airspeed data. Each AHRS output is then usd to convert IMU data
 * into a front acceleration, right acceleration and yaw rate as control inputs
 * used by a bank of 3-state EKF's whose states are velocity N,E (m/s) and yaw
 * angle (rad). The EKF's use the GPS horizontal NE velocity vector as observations
 * which enables the yaw angle to be estimated provided there is a change in
 * horizontal velocity vector. The output from each EKF is combined using a
 * Gaussian Sum Filter.
 *
 * This code has been developed from previous work by Paul Riseborough and
 * Rudaba Khan which includes the following Matlab code:
 *
 * https://github.com/priseborough/3_state_filter
 *
 * @author Paul Riseborough <gncsolns@gmail.com>
 *
 */

#include "ekf.h"
#include <ecl.h>
#include <mathlib/mathlib.h>
#include <cstdlib>

void Ekf::ahrsPredictEKFGSF(const uint8_t model_index)
{
	// generate attitude solution using simple complementary filter for the selected model

	// Accelerometer correction
	// Project 'k' unit vector of earth frame to body frame
	// Vector3f k = quaterion.conjugate_inversed(Vector3f(0.0f, 0.0f, 1.0f));
	// Optimized version with dropped zeros
	const Vector3f k(_ahrs_ekf_gsf[model_index].R(2,0), _ahrs_ekf_gsf[model_index].R(2,1), _ahrs_ekf_gsf[model_index].R(2,2));

	// Perform angular rate correction using accel data and reduce correction as accel magnitude moves away from 1 g (reduces drift when vehicle picked up and moved).
	// During fixed wing flight, compensate for centripetal acceleration assuming coordinated turns and X axis forward
	Vector3f correction = {};
	if (_ahrs_accel_fusion_gain > 0.0f) {

		Vector3f accel = _ahrs_accel;

		if (_ahrs_turn_comp_enabled) {
			// turn rate is component of gyro rate about vertical (down) axis
			const float turn_rate = _ahrs_ekf_gsf[model_index].R(2,0) * _ang_rate_delayed_raw(0)
					  + _ahrs_ekf_gsf[model_index].R(2,1) * _ang_rate_delayed_raw(1)
					  + _ahrs_ekf_gsf[model_index].R(2,2) * _ang_rate_delayed_raw(2);

			// use measured airspeed to calculate centripetal acceleration if available
			float centripetal_accel;
			if (_imu_sample_delayed.time_us - _airspeed_sample_delayed.time_us < 1000000) {
				centripetal_accel = _airspeed_sample_delayed.true_airspeed * turn_rate;
			} else {
				centripetal_accel = _params.EKFGSF_tas_default * turn_rate;
			}

			// project Y body axis onto horizontal and multiply by centripetal acceleration to give estimated
			// centripetal acceleration vector in earth frame due to coordinated turn
			Vector3f centripetal_accel_vec_ef = {_ahrs_ekf_gsf[model_index].R(0,1), _ahrs_ekf_gsf[model_index].R(1,1), 0.0f};
			if (_ahrs_ekf_gsf[model_index].R(2,2) > 0.0f) {
				// vehicle is upright
				centripetal_accel_vec_ef *= centripetal_accel;
			} else {
				// vehicle is inverted
				centripetal_accel_vec_ef *= - centripetal_accel;
			}

			// rotate into body frame
			Vector3f centripetal_accel_vec_bf = _ahrs_ekf_gsf[model_index].R.transpose() * centripetal_accel_vec_ef;

			// correct measured accel for centripetal acceleration
			accel -= centripetal_accel_vec_bf;
		}

		correction = (k % accel) * _ahrs_accel_fusion_gain / _ahrs_accel_norm;

	}

	// Gyro bias estimation
	const float gyro_bias_limit = 0.05f;
	const float spinRate = _ang_rate_delayed_raw.length();
	if (spinRate < 0.175f) {
		_ahrs_ekf_gsf[model_index].gyro_bias -= correction * (_params.EKFGSF_gyro_bias_gain * _imu_sample_delayed.delta_ang_dt);

		for (int i = 0; i < 3; i++) {
			_ahrs_ekf_gsf[model_index].gyro_bias(i) = math::constrain(_ahrs_ekf_gsf[model_index].gyro_bias(i), -gyro_bias_limit, gyro_bias_limit);
		}
	}

	const Vector3f rates = _ang_rate_delayed_raw - _ahrs_ekf_gsf[model_index].gyro_bias;

	// delta angle from previous to current frame
	Vector3f delta_angle = (correction + rates) * _imu_sample_delayed.delta_ang_dt;

	// Apply delta angle to rotation matrix
	_ahrs_ekf_gsf[model_index].R = updateRotMatEKFGSF(_ahrs_ekf_gsf[model_index].R, delta_angle);

}

void Ekf::ahrsAlignTiltEKFGSF()
{
	// Rotation matrix is constructed directly from acceleration measurement and will be the same for
	// all models so only need to calculate it once. Assumptions are:
	// 1) Yaw angle is zero - yaw is aligned later for each model when velocity fusion commences.
	// 2) The vehicle is not accelerating so all of the measured acceleration is due to gravity.

	// Calculate earth frame Down axis unit vector rotated into body frame
	Vector3f down_in_bf = -_imu_sample_delayed.delta_vel;
	down_in_bf.normalize();

	// Calculate earth frame North axis unit vector rotated into body frame, orthogonal to 'down_in_bf'
	// * operator is overloaded to provide a dot product
	const Vector3f i_vec_bf(1.0f,0.0f,0.0f);
	Vector3f north_in_bf = i_vec_bf - down_in_bf * (i_vec_bf * down_in_bf);
	north_in_bf.normalize();

	// Calculate earth frame East axis unit vector rotated into body frame, orthogonal to 'down_in_bf' and 'north_in_bf'
	// % operator is overloaded to provide a cross product
	const Vector3f east_in_bf = down_in_bf % north_in_bf;

	// Each column in a rotation matrix from earth frame to body frame represents the projection of the
	// corresponding earth frame unit vector rotated into the body frame, eg 'north_in_bf' would be the first column.
	// We need the rotation matrix from body frame to earth frame so the earth frame unit vectors rotated into body
	// frame are copied into corresponding rows instead.
	Dcmf R;
	R.setRow(0, north_in_bf);
	R.setRow(1, east_in_bf);
	R.setRow(2, down_in_bf);

}

void Ekf::ahrsAlignYawEKFGSF()
{
	// Align yaw angle for each model
	for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++) {
		if (fabsf(_ahrs_ekf_gsf[model_index].R(2, 0)) < fabsf(_ahrs_ekf_gsf[model_index].R(2, 1))) {
			// get the roll, pitch, yaw estimates from the rotation matrix using a  321 Tait-Bryan rotation sequence
			Eulerf euler_init(_ahrs_ekf_gsf[model_index].R);

			// set the yaw angle
			euler_init(2) = wrap_pi(_ekf_gsf[model_index].X[2]);

			// update the rotation matrix
			_ahrs_ekf_gsf[model_index].R = Dcmf(euler_init);

		} else {
			// Calculate the 312 Tait-Bryan rotation sequence that rotates from earth to body frame
			Vector3f rot312;
			rot312(0) = wrap_pi(_ekf_gsf[model_index].X[2]); // first rotation (yaw) taken from EKF model state
			rot312(1) = asinf(_ahrs_ekf_gsf[model_index].R(2, 1)); // second rotation (roll)
			rot312(2) = atan2f(-_ahrs_ekf_gsf[model_index].R(2, 0), _ahrs_ekf_gsf[model_index].R(2, 2));  // third rotation (pitch)

			// Calculate the body to earth frame rotation matrix
			_ahrs_ekf_gsf[model_index].R = taitBryan312ToRotMat(rot312);

		}
		_ahrs_ekf_gsf[model_index].aligned = true;
	}
}

// predict states and covariance for specified model index
void Ekf::statePredictEKFGSF(const uint8_t model_index)
{
	// generate an attitude reference using IMU data
	ahrsPredictEKFGSF(model_index);

	// we don't start running the EKF part of the algorithm until there are regular velocity observations
	if (!_ekf_gsf_vel_fuse_started) {
		return;
	}

	// Calculate the yaw state using a projection onto the horizontal that avoids gimbal lock
	if (fabsf(_ahrs_ekf_gsf[model_index].R(2, 0)) < fabsf(_ahrs_ekf_gsf[model_index].R(2, 1))) {
		// use 321 Tait-Bryan rotation to define yaw state
		_ekf_gsf[model_index].X[2] = atan2f(_ahrs_ekf_gsf[model_index].R(1, 0), _ahrs_ekf_gsf[model_index].R(0, 0));
	} else {
		// use 312 Tait-Bryan rotation to define yaw state
		_ekf_gsf[model_index].X[2] = atan2f(-_ahrs_ekf_gsf[model_index].R(0, 1), _ahrs_ekf_gsf[model_index].R(1, 1)); // first rotation (yaw)
	}

	// calculate delta velocity in a horizontal front-right frame
	const Vector3f del_vel_NED = _ahrs_ekf_gsf[model_index].R * _imu_sample_delayed.delta_vel;
	const float dvx =   del_vel_NED(0) * cosf(_ekf_gsf[model_index].X[2]) + del_vel_NED(1) * sinf(_ekf_gsf[model_index].X[2]);
	const float dvy = - del_vel_NED(0) * sinf(_ekf_gsf[model_index].X[2]) + del_vel_NED(1) * cosf(_ekf_gsf[model_index].X[2]);

	// sum delta velocities in earth frame:
	_ekf_gsf[model_index].X[0] += del_vel_NED(0);
	_ekf_gsf[model_index].X[1] += del_vel_NED(1);

	// predict covariance - autocode from https://github.com/priseborough/3_state_filter/blob/flightLogReplay-wip/calcPupdate.txt

	// Local short variable name copies required for readability
	// Compiler might be smart enough to optimise these out
	const float P00 = _ekf_gsf[model_index].P[0][0];
	const float P01 = _ekf_gsf[model_index].P[0][1];
	const float P02 = _ekf_gsf[model_index].P[0][2];
	const float P10 = _ekf_gsf[model_index].P[1][0];
	const float P11 = _ekf_gsf[model_index].P[1][1];
	const float P12 = _ekf_gsf[model_index].P[1][2];
	const float P20 = _ekf_gsf[model_index].P[2][0];
	const float P21 = _ekf_gsf[model_index].P[2][1];
	const float P22 = _ekf_gsf[model_index].P[2][2];

	// Use fixed values for delta velocity and delta angle process noise variances
	const float dvxVar = sq(_params.EKFGSF_accel_noise * _imu_sample_delayed.delta_vel_dt); // variance of forward delta velocity - (m/s)^2
	const float dvyVar = dvxVar; // variance of right delta velocity - (m/s)^2
	const float dazVar = sq(_params.EKFGSF_gyro_noise * _imu_sample_delayed.delta_ang_dt); // variance of yaw delta angle - rad^2

	const float t2 = sinf(_ekf_gsf[model_index].X[2]);
	const float t3 = cosf(_ekf_gsf[model_index].X[2]);
	const float t4 = dvy*t3;
	const float t5 = dvx*t2;
	const float t6 = t4+t5;
	const float t8 = P22*t6;
	const float t7 = P02-t8;
	const float t9 = dvx*t3;
	const float t11 = dvy*t2;
	const float t10 = t9-t11;
	const float t12 = dvxVar*t2*t3;
	const float t13 = t2*t2;
	const float t14 = t3*t3;
	const float t15 = P22*t10;
	const float t16 = P12+t15;

	const float min_var = 1e-6f;
	_ekf_gsf[model_index].P[0][0] = fmaxf(P00-P20*t6+dvxVar*t14+dvyVar*t13-t6*t7 , min_var);
	_ekf_gsf[model_index].P[0][1] = P01+t12-P21*t6+t7*t10-dvyVar*t2*t3;
	_ekf_gsf[model_index].P[0][2] = t7;
	_ekf_gsf[model_index].P[1][0] = P10+t12+P20*t10-t6*t16-dvyVar*t2*t3;
	_ekf_gsf[model_index].P[1][1] = fmaxf(P11+P21*t10+dvxVar*t13+dvyVar*t14+t10*t16 , min_var);
	_ekf_gsf[model_index].P[1][2] = t16;
	_ekf_gsf[model_index].P[2][0] = P20-t8;
	_ekf_gsf[model_index].P[2][1] = P21+t15;
	_ekf_gsf[model_index].P[2][2] = fmaxf(P22+dazVar , min_var);

	// force symmetry
	makeCovSymEKFGSF(model_index);
}

// Update EKF states and covariance for specified model index using velocity measurement
bool Ekf::stateUpdateEKFGSF(const uint8_t model_index)
{
	// set observation variance from accuracy estimate supplied by GPS and apply a sanity check minimum
	const float velObsVar = sq(fmaxf(_gps_sample_delayed.sacc, _params.gps_vel_noise));

	// calculate velocity observation innovations
	_ekf_gsf[model_index].innov[0] = _ekf_gsf[model_index].X[0] - _gps_sample_delayed.vel(0);
	_ekf_gsf[model_index].innov[1] = _ekf_gsf[model_index].X[1] - _gps_sample_delayed.vel(1);

	// copy covariance matrix to temporary variables
	const float P00 = _ekf_gsf[model_index].P[0][0];
	const float P01 = _ekf_gsf[model_index].P[0][1];
	const float P02 = _ekf_gsf[model_index].P[0][2];
	const float P10 = _ekf_gsf[model_index].P[1][0];
	const float P11 = _ekf_gsf[model_index].P[1][1];
	const float P12 = _ekf_gsf[model_index].P[1][2];
	const float P20 = _ekf_gsf[model_index].P[2][0];
	const float P21 = _ekf_gsf[model_index].P[2][1];
	const float P22 = _ekf_gsf[model_index].P[2][2];

	// calculate innovation variance
	_ekf_gsf[model_index].S[0][0] = P00 + velObsVar;
	_ekf_gsf[model_index].S[1][1] = P11 + velObsVar;
	_ekf_gsf[model_index].S[0][1] = P01;
	_ekf_gsf[model_index].S[1][0] = P10;

	// Perform a chi-square innovation consistency test and calculate a compression scale factor that limits the magnitude of innovations to 5-sigma
	float S_det_inv = (_ekf_gsf[model_index].S[0][0]*_ekf_gsf[model_index].S[1][1] - _ekf_gsf[model_index].S[0][1]*_ekf_gsf[model_index].S[1][0]);
	float innov_comp_scale_factor = 1.0f;
	if (fabsf(S_det_inv) > 1E-6f) {
		// Calculate elements for innovation covariance inverse matrix assuming symmetry
		S_det_inv = 1.0f / S_det_inv;
		const float S_inv_NN = _ekf_gsf[model_index].S[1][1] * S_det_inv;
		const float S_inv_EE = _ekf_gsf[model_index].S[0][0] * S_det_inv;
		const float S_inv_NE = _ekf_gsf[model_index].S[0][1] * S_det_inv;

		// The following expression was derived symbolically from test ratio = transpose(innovation) * inverse(innovation variance) * innovation = [1x2] * [2,2] * [2,1] = [1,1]
		const float test_ratio = _ekf_gsf[model_index].innov[0]*(_ekf_gsf[model_index].innov[0]*S_inv_NN + _ekf_gsf[model_index].innov[1]*S_inv_NE) + _ekf_gsf[model_index].innov[1]*(_ekf_gsf[model_index].innov[0]*S_inv_NE + _ekf_gsf[model_index].innov[1]*S_inv_EE);

		// If the test ratio is greater than 25 (5 Sigma) then reduce the length of the innovation vector to clip it at 5-Sigma
		// This protects from large measurement spikes
		if (test_ratio > 25.0f) {
			innov_comp_scale_factor = sqrtf(25.0f / test_ratio);
		}
	} else {
		// skip this fusion step because calculation is badly conditioned
		return false;
	}

	// calculate Kalman gain K and covariance matrix P
	// autocode from https://github.com/priseborough/3_state_filter/blob/flightLogReplay-wip/calcK.txt
	// and https://github.com/priseborough/3_state_filter/blob/flightLogReplay-wip/calcPmat.txt
	const float t2 = P00*velObsVar;
 	const float t3 = P11*velObsVar;
	const float t4 = velObsVar*velObsVar;
	const float t5 = P00*P11;
	const float t9 = P01*P10;
	const float t6 = t2+t3+t4+t5-t9;
	float t7;
	if (fabsf(t6) > 1e-6f) {
		t7 = 1.0f/t6;
	} else {
		// skip this fusion step
		return false;
	}
	const float t8 = P11+velObsVar;
	const float t10 = P00+velObsVar;
	float K[3][2];

 	K[0][0] = -P01*P10*t7+P00*t7*t8;
	K[0][1] = -P00*P01*t7+P01*t7*t10;
	K[1][0] = -P10*P11*t7+P10*t7*t8;
	K[1][1] = -P01*P10*t7+P11*t7*t10;
	K[2][0] = -P10*P21*t7+P20*t7*t8;
	K[2][1] = -P01*P20*t7+P21*t7*t10;

	const float t11 = P00*P01*t7;
	const float t15 = P01*t7*t10;
	const float t12 = t11-t15;
	const float t13 = P01*P10*t7;
	const float t16 = P00*t7*t8;
	const float t14 = t13-t16;
	const float t17 = t8*t12;
	const float t18 = P01*t14;
	const float t19 = t17+t18;
	const float t20 = t10*t14;
	const float t21 = P10*t12;
	const float t22 = t20+t21;
	const float t27 = P11*t7*t10;
	const float t23 = t13-t27;
	const float t24 = P10*P11*t7;
	const float t26 = P10*t7*t8;
	const float t25 = t24-t26;
	const float t28 = t8*t23;
	const float t29 = P01*t25;
	const float t30 = t28+t29;
	const float t31 = t10*t25;
	const float t32 = P10*t23;
	const float t33 = t31+t32;
	const float t34 = P01*P20*t7;
	const float t38 = P21*t7*t10;
	const float t35 = t34-t38;
	const float t36 = P10*P21*t7;
	const float t39 = P20*t7*t8;
	const float t37 = t36-t39;
	const float t40 = t8*t35;
	const float t41 = P01*t37;
	const float t42 = t40+t41;
	const float t43 = t10*t37;
	const float t44 = P10*t35;
	const float t45 = t43+t44;

	const float min_var = 1e-6f;
	_ekf_gsf[model_index].P[0][0] = fmaxf(P00-t12*t19-t14*t22 , min_var);
	_ekf_gsf[model_index].P[0][1] = P01-t19*t23-t22*t25;
	_ekf_gsf[model_index].P[0][2] = P02-t19*t35-t22*t37;
	_ekf_gsf[model_index].P[1][0] = P10-t12*t30-t14*t33;
	_ekf_gsf[model_index].P[1][1] = fmaxf(P11-t23*t30-t25*t33 , min_var);
	_ekf_gsf[model_index].P[1][2] = P12-t30*t35-t33*t37;
	_ekf_gsf[model_index].P[2][0] = P20-t12*t42-t14*t45;
	_ekf_gsf[model_index].P[2][1] = P21-t23*t42-t25*t45;
	_ekf_gsf[model_index].P[2][2] = fmaxf(P22-t35*t42-t37*t45 , min_var);

	// force symmetry
	makeCovSymEKFGSF(model_index);

	float oldYaw = _ekf_gsf[model_index].X[2];
	for (uint8_t obs_index = 0; obs_index < 2; obs_index++) {
		// apply the state corrections including the compression scale factor
		for (unsigned row = 0; row < 3; row++) {
			_ekf_gsf[model_index].X[row] -= K[row][obs_index] * _ekf_gsf[model_index].innov[obs_index] * innov_comp_scale_factor;
		}
	}
	float yawDelta = _ekf_gsf[model_index].X[2] - oldYaw;

	// apply the change in yaw angle to the AHRS
	// take advantage of sparseness in the yaw rotation matrix
	const float cosYaw = cosf(yawDelta);
	const float sinYaw = sinf(yawDelta);
	float  R_prev[2][3];
	memcpy(&R_prev, &_ahrs_ekf_gsf[model_index].R, sizeof(R_prev));
	_ahrs_ekf_gsf[model_index].R(0,0) = R_prev[0][0] * cosYaw - R_prev[1][0] * sinYaw;
	_ahrs_ekf_gsf[model_index].R(0,1) = R_prev[0][1] * cosYaw - R_prev[1][1] * sinYaw;
	_ahrs_ekf_gsf[model_index].R(0,2) = R_prev[0][2] * cosYaw - R_prev[1][2] * sinYaw;
	_ahrs_ekf_gsf[model_index].R(1,0) = R_prev[0][0] * sinYaw + R_prev[1][0] * cosYaw;
	_ahrs_ekf_gsf[model_index].R(1,1) = R_prev[0][1] * sinYaw + R_prev[1][1] * cosYaw;
	_ahrs_ekf_gsf[model_index].R(1,2) = R_prev[0][2] * sinYaw + R_prev[1][2] * cosYaw;

	return true;
}

void Ekf::initialiseEKFGSF()
{
	memset(&X_GSF, 0, sizeof(X_GSF));
	_ekf_gsf_vel_fuse_started = false;
	_emergency_yaw_reset_time = 0;
	_time_last_on_ground_us = 0;
	_yaw_extreset_counter = 0;
	_do_emergency_yaw_reset = false;
	_ekf_gsf_yaw_variance = sq(M_PI_2_F);
	memset(&_ekf_gsf, 0, sizeof(_ekf_gsf));
	const float yaw_increment = 2.0f * M_PI_F / (float)N_MODELS_EKFGSF;
	for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++) {
		// evenly space initial yaw estimates in the region between +-Pi
		_ekf_gsf[model_index].X[2] = -M_PI_F + (0.5f * yaw_increment) + ((float)model_index * yaw_increment);

		// All filter models start with the same weight
		_ekf_gsf[model_index].W = 1.0f / (float)N_MODELS_EKFGSF;

		if (_control_status.flags.gps) {
			_ekf_gsf[model_index].X[0] = _gps_sample_delayed.vel(0);
			_ekf_gsf[model_index].X[1] = _gps_sample_delayed.vel(1);
			_ekf_gsf[model_index].P[0][0] = sq(_gps_sample_delayed.sacc);
			_ekf_gsf[model_index].P[1][1] = _ekf_gsf[model_index].P[0][0];
		} else {
			// no GPS so assume we are in a semi-static condition and moving no more that 0.5 m/s
			_ekf_gsf[model_index].P[0][0] = sq(0.5f);
			_ekf_gsf[model_index].P[1][1] = sq(0.5f);
		}

		// use half yaw interval for yaw uncertainty
		_ekf_gsf[model_index].P[2][2] = sq(0.5f * yaw_increment);
	}
}

void Ekf::runEKFGSF()
{
	// Iniitialise states first time
	_ahrs_accel = _imu_sample_delayed.delta_vel / fmaxf(_imu_sample_delayed.delta_vel_dt, FILTER_UPDATE_PERIOD_S / 4);
	if (!_ahrs_ekf_gsf_tilt_aligned) {
		// check for excessive acceleration to reduce likelihood of large inital roll/pitch errors
		// due to vehicle movement
		const float accel_norm_sq = _ahrs_accel.norm_squared();
		const float upper_accel_limit = CONSTANTS_ONE_G * 1.1f;
		const float lower_accel_limit = CONSTANTS_ONE_G * 0.9f;
		const bool ok_to_align = ((accel_norm_sq > lower_accel_limit * lower_accel_limit &&
			  accel_norm_sq < upper_accel_limit * upper_accel_limit));
		if (ok_to_align) {
			initialiseEKFGSF();
			ahrsAlignTiltEKFGSF();
			_ahrs_ekf_gsf_tilt_aligned = true;
		}
		return;
	}

	// AHRS prediction cycle for each model - this always runs
	calcAccelGainEKFGSF();
	for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
		statePredictEKFGSF(model_index);
	}

	// The 3-state EKF models only run when flying to avoid corrupted estimates due to operator handling and GPS interference
	if (_control_status.flags.gps && _gps_data_ready && _control_status.flags.in_air) {
		if (!_ekf_gsf_vel_fuse_started) {
			for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
				// use the first measurement to set the velocities and corresponding covariances
				_ekf_gsf[model_index].X[0] = _gps_sample_delayed.vel(0);
				_ekf_gsf[model_index].X[1] = _gps_sample_delayed.vel(1);
				_ekf_gsf[model_index].P[0][0] = sq(_gps_sample_delayed.sacc);
				_ekf_gsf[model_index].P[1][1] = _ekf_gsf[model_index].P[0][0];
			}
			ahrsAlignYawEKFGSF();
			_ekf_gsf_vel_fuse_started = true;
		} else {
			float total_w = 0.0f;
			float newWeight[N_MODELS_EKFGSF];
			bool weight_update_inhibit = false;
			for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
				// subsequent measurements are fused as direct state observations
				if (!stateUpdateEKFGSF(model_index)) {
					weight_update_inhibit = true;
				}
			}

			for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
				if (weight_update_inhibit) {
					// keep old weights
					newWeight[model_index] = _ekf_gsf[model_index].W;
				} else {
					// calculate weighting for each model assuming a normal distribution
					newWeight[model_index] = fmaxf(gaussianDensityEKFGSF(model_index) * _ekf_gsf[model_index].W, 0.0f);
				}
				total_w += newWeight[model_index];
			}

			// normalise the weighting function
			if (_ekf_gsf_vel_fuse_started && total_w > 1e-6f) {
				float total_w_inv = 1.0f / total_w;
				for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
					_ekf_gsf[model_index].W  = newWeight[model_index] * total_w_inv;
				}
			}

			// Enforce a minimum weighting value. This was added during initial development but has not been needed
			// subsequently, so this block of code and the corresponding EKFGSF_weight_min can be removed if we get
			// through testing without any weighting function issues.
			if (_params.EKFGSF_weight_min > FLT_EPSILON) {
				float correction_sum = 0.0f; // amount the sum of weights has been increased by application of the limit
				bool change_mask[N_MODELS_EKFGSF] = {}; // true when the weighting for that model has been increased
				float unmodified_weights_sum = 0.0f; // sum of unmodified weights
				for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
					if (_ekf_gsf[model_index].W < _params.EKFGSF_weight_min) {
						correction_sum += _params.EKFGSF_weight_min - _ekf_gsf[model_index].W;
						_ekf_gsf[model_index].W = _params.EKFGSF_weight_min;
						change_mask[model_index] = true;
					} else {
						unmodified_weights_sum += _ekf_gsf[model_index].W;
					}
				}

				// rescale the unmodified weights to make the total sum unity
				const float scale_factor = (unmodified_weights_sum - correction_sum - _params.EKFGSF_weight_min) / (unmodified_weights_sum - _params.EKFGSF_weight_min);
				for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
					if (!change_mask[model_index]) {
						_ekf_gsf[model_index].W = _params.EKFGSF_weight_min + scale_factor * (_ekf_gsf[model_index].W - _params.EKFGSF_weight_min);
					}
				}
			}
		}
	} else if (_ekf_gsf_vel_fuse_started && !_control_status.flags.in_air) {
		// reset EKF states and wait to fly again
		initialiseEKFGSF();
		_ekf_gsf_vel_fuse_started = false;
	}

	// Calculate a composite state vector as a weighted average of the states for each model.
	// To avoid issues with angle wrapping, the yaw state is converted to a vector with legnth
	// equal to the weighting value before it is summed.
	memset(&X_GSF, 0, sizeof(X_GSF));
	Vector2f yaw_vector = {};
	for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
		for (uint8_t state_index = 0; state_index < 2; state_index++) {
			X_GSF[state_index] += _ekf_gsf[model_index].X[state_index] * _ekf_gsf[model_index].W;
		}
		yaw_vector(0) += _ekf_gsf[model_index].W * cosf(_ekf_gsf[model_index].X[2]);
		yaw_vector(1) += _ekf_gsf[model_index].W * sinf(_ekf_gsf[model_index].X[2]);
	}
	X_GSF[2] = atan2f(yaw_vector(1),yaw_vector(0));

	/* Example of how a full covaraince matrix can be generated for the GSF state vector
	// calculate a composite covariance matrix from a weighted average of the covariance for each model
	// models with larger innovations are weighted less
	memset(&P_GSF, 0, sizeof(P_GSF));
	for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
		float Xdelta[3];
		for (uint8_t row = 0; row < 3; row++) {
			Xdelta[row] = _ekf_gsf[model_index].X[row] - X_GSF[row];
		}
		for (uint8_t row = 0; row < 3; row++) {
			for (uint8_t col = 0; col < 3; col++) {
				P_GSF[row][col] +=  _ekf_gsf[model_index].W * (_ekf_gsf[model_index].P[row][col] + Xdelta[row] * Xdelta[col]);
			}
		}
	}
	*/

	// calculate a composite variance for the yaw state from a weighted average of the variance for each model
	// models with larger innovations are weighted less
	_ekf_gsf_yaw_variance = 0.0f;
	for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
		float yaw_delta = wrap_pi(_ekf_gsf[model_index].X[2] - X_GSF[2]);
		_ekf_gsf_yaw_variance += _ekf_gsf[model_index].W * (_ekf_gsf[model_index].P[2][2] + sq(yaw_delta));
	}
}

float Ekf::gaussianDensityEKFGSF(const uint8_t model_index) const
{
	const float t2 = _ekf_gsf[model_index].S[0][0] * _ekf_gsf[model_index].S[1][1];
	const float t5 = _ekf_gsf[model_index].S[0][1] * _ekf_gsf[model_index].S[1][0];
	const float t3 = t2 - t5; // determinant
	float t4; // determinant inverse
	if (fabsf(t3) > 1e-6f) {
		t4 = 1.0f / t3;
	} else {
		t4 = 1.0f / fmaxf(t2, 1e-6f);
	}

	// inv(S)
	float invMat[2][2];
	invMat[0][0] =   t4 * _ekf_gsf[model_index].S[1][1];
	invMat[1][1] =   t4 * _ekf_gsf[model_index].S[0][0];
	invMat[0][1] = - t4 * _ekf_gsf[model_index].S[0][1];
	invMat[1][0] = - t4 * _ekf_gsf[model_index].S[1][0];

 	// inv(S) * innovation
	float tempVec[2];
	tempVec[0] = invMat[0][0] * _ekf_gsf[model_index].innov[0] + invMat[0][1] * _ekf_gsf[model_index].innov[1];
	tempVec[1] = invMat[1][0] * _ekf_gsf[model_index].innov[0] + invMat[1][1] * _ekf_gsf[model_index].innov[1];

	// transpose(innovation) * inv(S) * innovation
	float normDist = tempVec[0] * _ekf_gsf[model_index].innov[0] + tempVec[1] * _ekf_gsf[model_index].innov[1];

	return M_TWOPI_INV * sqrtf(t4) * expf(-0.5f * normDist);
}

bool Ekf::getDataEKFGSF(float *yaw_composite, float *yaw_variance, float yaw[N_MODELS_EKFGSF], float innov_VN[N_MODELS_EKFGSF], float innov_VE[N_MODELS_EKFGSF], float weight[N_MODELS_EKFGSF]) const
{
	if (_ekf_gsf_vel_fuse_started) {
		memcpy(yaw_composite, &X_GSF[2], sizeof(X_GSF[2]));
		memcpy(yaw_variance, &_ekf_gsf_yaw_variance, sizeof(_ekf_gsf_yaw_variance));
		for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++) {
			yaw[model_index] = _ekf_gsf[model_index].X[2];
			innov_VN[model_index] = _ekf_gsf[model_index].innov[0];
			innov_VE[model_index] = _ekf_gsf[model_index].innov[1];
			weight[model_index] = _ekf_gsf[model_index].W;
		}
		return true;
	}
	return false;
}

void Ekf::makeCovSymEKFGSF(const uint8_t model_index)
{
	float P01 = 0.5f * (_ekf_gsf[model_index].P[0][1] + _ekf_gsf[model_index].P[1][0]);
	float P02 = 0.5f * (_ekf_gsf[model_index].P[0][2] + _ekf_gsf[model_index].P[2][0]);
	float P12 = 0.5f * (_ekf_gsf[model_index].P[1][2] + _ekf_gsf[model_index].P[2][1]);
	_ekf_gsf[model_index].P[0][1] = _ekf_gsf[model_index].P[1][0] = P01;
	_ekf_gsf[model_index].P[0][2] = _ekf_gsf[model_index].P[2][0] = P02;
	_ekf_gsf[model_index].P[1][2] = _ekf_gsf[model_index].P[2][1] = P12;
}

// Reset main nav filter yaw to value from EKF-GSF and reset velocity and position to GPS
bool Ekf::resetYawToEKFGSF()
{
	// don't allow reset using the EKF-GSF estimate until the filter has started fusing velocity
	// data and the yaw estimate has converged
	if (_ekf_gsf_vel_fuse_started && _ekf_gsf_yaw_variance < sq(_params.EKFGSF_yaw_err_max)) {
		// save a copy of the quaternion state for later use in calculating the amount of reset change
		Quatf quat_before_reset = _state.quat_nominal;
		Quatf quat_after_reset = _state.quat_nominal;

		// update transformation matrix from body to world frame using the current estimate
		_R_to_earth = Dcmf(_state.quat_nominal);

		// calculate the initial quaternion
		// determine if a 321 or 312 Euler sequence is best
		if (fabsf(_R_to_earth(2, 0)) < fabsf(_R_to_earth(2, 1))) {
			// use a 321 sequence

			Eulerf euler321(_state.quat_nominal);
			euler321(2) = X_GSF[2];
			quat_after_reset = Quatf(euler321);

		} else {
			// Calculate the 312 Tait-Bryan rotation sequence that rotates from earth to body frame
			// We use a 312 sequence as an alternate when there is more pitch tilt than roll tilt
			// to avoid gimbal lock
			Vector3f rot312;
			rot312(0) = X_GSF[2]; // first rotation about Z is taken from EKF-GSF
			rot312(1) = asinf(_R_to_earth(2, 1)); // second rotation about X is taken from main EKF
			rot312(2) = atan2f(-_R_to_earth(2, 0), _R_to_earth(2, 2));  // third rotation about Y is taken from main EKF

			// Calculate the body to earth frame rotation matrix
			const Dcmf R = taitBryan312ToRotMat(rot312);

			// calculate initial quaternion states for the main EKF
			// we don't change the output attitude to avoid jumps
			quat_after_reset = Quatf(R);
		}

		// record the time for the magnetic field alignment event
		_flt_mag_align_start_time = _imu_sample_delayed.time_us;

		// calculate the amount that the quaternion has changed by
		Quatf q_error =  quat_after_reset * quat_before_reset.inversed();
		q_error.normalize();

		// update quaternion states
		_state.quat_nominal = quat_after_reset;
		uncorrelateQuatStates();

		// record the state change
		_state_reset_status.quat_change = q_error;

		// update transformation matrix from body to world frame using the current estimate
		_R_to_earth = Dcmf(_state.quat_nominal);

		// reset the rotation from the EV to EKF frame of reference if it is being used
		if ((_params.fusion_mode & MASK_ROTATE_EV) && !_control_status.flags.ev_yaw) {
			resetExtVisRotMat();
		}

		// update the yaw angle variance using half the nominal yaw separation between models
		increaseQuatYawErrVariance(sq(fmaxf(M_PI_F / (float)N_MODELS_EKFGSF, 1.0e-2f)));

		// add the reset amount to the output observer buffered data
		for (uint8_t i = 0; i < _output_buffer.get_length(); i++) {
			_output_buffer[i].quat_nominal = _state_reset_status.quat_change * _output_buffer[i].quat_nominal;
		}

		// apply the change in attitude quaternion to our newest quaternion estimate
		// which was already taken out from the output buffer
		_output_new.quat_nominal = _state_reset_status.quat_change * _output_new.quat_nominal;

		// capture the reset event
		_state_reset_status.quat_counter++;

		// reset velocity and position states to GPS - if yaw is fixed then the filter should start to operate correctly
		resetVelocity();
		resetPosition();

		// stop using the magnetometer in the main EKF otherwise it's fusion could drag the yaw around
		// and cause another navigation failure
		_control_status.flags.mag_fault = true;

		ECL_INFO_TIMESTAMPED("Emergency yaw reset - mag use stopped");

		return true;
	}

	return false;

}

void Ekf::requestEmergencyNavReset(uint8_t counter)
{
	if (counter > _yaw_extreset_counter) {
		_yaw_extreset_counter = counter;
		_do_emergency_yaw_reset = true;
	}
}

// converts Tait-Bryan 312 sequence of rotations from frame 1 to frame 2
// to the corresponding rotation matrix that rotates from frame 2 to frame 1
// rot312(0) - First rotation is a RH rotation about the Z axis (rad)
// rot312(1) - Second rotation is a RH rotation about the X axis (rad)
// rot312(2) - Third rotation is a RH rotation about the Y axis (rad)
// See http://www.atacolorado.com/eulersequences.doc
Dcmf Ekf::taitBryan312ToRotMat(Vector3f &rot312)
{
		// Calculate the frame2 to frame 1 rotation matrix from a 312 rotation sequence
		const float c2 = cosf(rot312(2));
		const float s2 = sinf(rot312(2));
		const float s1 = sinf(rot312(1));
		const float c1 = cosf(rot312(1));
		const float s0 = sinf(rot312(0));
		const float c0 = cosf(rot312(0));

		Dcmf R;
		R(0, 0) = c0 * c2 - s0 * s1 * s2;
		R(1, 1) = c0 * c1;
		R(2, 2) = c2 * c1;
		R(0, 1) = -c1 * s0;
		R(0, 2) = s2 * c0 + c2 * s1 * s0;
		R(1, 0) = c2 * s0 + s2 * s1 * c0;
		R(1, 2) = s0 * s2 - s1 * c0 * c2;
		R(2, 0) = -s2 * c1;
		R(2, 1) = s1;

		return R;
}

void Ekf::calcAccelGainEKFGSF()
{
	// calculate common values used by the AHRS complementary filter models
	_ahrs_accel_norm = _ahrs_accel.norm();
	_ahrs_turn_comp_enabled = _control_status.flags.fixed_wing && _params.EKFGSF_tas_default > FLT_EPSILON;

	// Calculate the acceleration fusion gain using a continuous function that is unity at 1g and zero
	// at the min and mas g value. Allow for more acceleration when flying as a fixed wing vehicle using centripetal
	// acceleration correction as higher and more sustained g will be experienced.
	// Use a quadratic finstead of linear unction to prevent vibration around 1g reducing the tilt correction effectiveness.
	float accel_g = _ahrs_accel_norm / CONSTANTS_ONE_G;
	if (accel_g > 1.0f) {
		if (_ahrs_turn_comp_enabled && accel_g < 2.0f) {
			_ahrs_accel_fusion_gain = _params.EKFGSF_tilt_gain * sq(2.0f - accel_g);
		} else if (accel_g < 1.5f) {
			_ahrs_accel_fusion_gain = _params.EKFGSF_tilt_gain * sq(3.0f - 2.0f * accel_g);
		} else {
			_ahrs_accel_fusion_gain = 0.0f;
		}
	} else if (accel_g > 0.5f) {
		_ahrs_accel_fusion_gain = _params.EKFGSF_tilt_gain * sq(2.0f * accel_g - 1.0f);
	} else {
		_ahrs_accel_fusion_gain = 0.0f;
	}
}

// Efficient propagation of a delta angle in body frame applied to the body to earth frame rotation matrix
Matrix3f Ekf::updateRotMatEKFGSF(Matrix3f &R, Vector3f &g)
{
	Matrix3f ret = R;
	ret(0,0) += R(0,1) * g(2) - R(0,2) * g(1);
	ret(0,1) += R(0,2) * g(0) - R(0,0) * g(2);
	ret(0,2) += R(0,0) * g(1) - R(0,1) * g(0);
	ret(1,0) += R(1,1) * g(2) - R(1,2) * g(1);
	ret(1,1) += R(1,2) * g(0) - R(1,0) * g(2);
	ret(1,2) += R(1,0) * g(1) - R(1,1) * g(0),
        ret(2,0) += R(2,1) * g(2) - R(2,2) * g(1);
	ret(2,1) += R(2,2) * g(0) - R(2,0) * g(2);
	ret(2,2) += R(2,0) * g(1) - R(2,1) * g(0);

	// Renormalise rows
	float rowLengthSq;
	for (uint8_t r = 0; r < 3; r++) {
		rowLengthSq = ret(r,0) * ret(r,0) + ret(r,1) * ret(r,1) + ret(r,2) * ret(r,2);
		if (rowLengthSq > FLT_EPSILON) {
			// Use linear approximation for inverse sqrt taking advantage of the row length being close to 1.0
			const float rowLengthInv = 1.5f - 0.5f * rowLengthSq;
			ret(r,0) *= rowLengthInv;
			ret(r,1) *= rowLengthInv;
			ret(r,2) *= rowLengthInv;
		}
        }

	return ret;
}
