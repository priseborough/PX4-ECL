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
 * @file mag_calibration.cpp
 * Magnetometer calibration methods.
 *
 * @author Paul Riseborough <p_riseborough@live.com.au>
 *
 */

#include "ekf.h"
#include <ecl.h>
#include <mathlib/mathlib.h>

void Ekf::runMagCal()
{
	// reset all class variables required to restart the sampling and processing steps
	if (_restart_sampling) {
		_mag_cal_sampling_complete = false;
		_mag_cal_complete = false;
		_mag_sample_index = 0;
		_mag_cal_iteration_index = 0;
		_mag_cal_sample_time_us = 0;
		_mag_cal_direction = 0;
		_retry_count = 0;
		_restart_sampling = false;
	}

	// calculate the AHRS solution for the mag calibrator using IMU data only
	// we need to keep running this so a calibration can be restarted
	calcQuatMagCal();

	if (!_mag_cal_complete) {
		sampleDataMagCal();

		if (_mag_cal_sampling_complete) {
			processDataMagCal();
		}
	}

	// print results if calibration completes successfully
	if (_mag_cal_complete) {

		printf("completed at %5.1f sec after %u iterations\n",1E-6*(double)_imu_sample_new.time_us,(unsigned)_mag_cal_iteration_index);
		printf("bias = %5.3f,%5.3f,%5.3f yaw = %5.3f,\n",
		(double)_mag_cal_states.mag_bias(0), (double)_mag_cal_states.mag_bias(1), (double)_mag_cal_states.mag_bias(2),
		(double)_mag_cal_states.yaw_offset);
		printf("residuals = %5.3f,%5.3f,%5.3f\n\n",(double)_mag_cal_residual[0],(double)_mag_cal_residual[1],(double)_mag_cal_residual[2]);
		_mag_cal_complete = false;
	}
}

void Ekf::calcQuatMagCal()
{
	// generate attitude solution using simple complementary filter

	// check for excessive acceleration.
	Vector3f accel = _imu_sample_delayed.delta_vel / _imu_sample_delayed.delta_vel_dt;
	const float accel_norm_sq = accel.norm_squared();
	const float upper_accel_limit = CONSTANTS_ONE_G * 1.1f;
	const float lower_accel_limit = CONSTANTS_ONE_G * 0.9f;

	// Angular rate of correction
	bool ok_to_align = ((accel_norm_sq > lower_accel_limit * lower_accel_limit &&
			  accel_norm_sq < upper_accel_limit * upper_accel_limit));

	// Iniitialise quaternion first time
	if (!_mag_cal_quat_initialised) {
		if (ok_to_align) {
			// Rotation matrix can be easily constructed from acceleration
			// assuming zero yaw
			// 'k' is Earth Z axis (Down) unit vector in body frame
			Vector3f k_init = -_imu_sample_delayed.delta_vel;
			k_init.normalize();

			// 'i' is Earth X axis (North) unit vector in body frame, orthogonal with 'k'
			Vector3f temp(1.0f,0.0f,0.0f);
			Vector3f i_init = (temp - k_init * (temp * k_init));
			i_init.normalize();

			// 'j' is Earth Y axis (East) unit vector in body frame, orthogonal with 'k' and 'i'
			Vector3f j_init = k_init % i_init;

			// Fill rotation matrix
			Dcmf R;
			R.setRow(0, i_init);
			R.setRow(1, j_init);
			R.setRow(2, k_init);

			// Convert to quaternion
			_mag_cal_quat = R;

			_mag_cal_quat_initialised = true;
// debug print to compare between this and main EKF
Eulerf temp1(_mag_cal_quat);
Eulerf temp2(_state.quat_nominal);
printf("roll=%5.1f,%5.1f\n",(double)math::degrees(temp1(0)),(double)math::degrees(temp2(0)));
printf("pitch=%5.1f,%5.1f\n",(double)math::degrees(temp1(1)),(double)math::degrees(temp2(1)));
printf("yaw=%5.1f,%5.1f\n",(double)math::degrees(temp1(2)),(double)math::degrees(temp2(2)));

		}
	} else {
		// Accelerometer correction
		// Project 'k' unit vector of earth frame to body frame
		// Vector3f k = quaterion.conjugate_inversed(Vector3f(0.0f, 0.0f, 1.0f));
		// Optimized version with dropped zeros
		Vector3f k(
			2.0f * (_mag_cal_quat(1) * _mag_cal_quat(3) - _mag_cal_quat(0) * _mag_cal_quat(2)),
			2.0f * (_mag_cal_quat(2) * _mag_cal_quat(3) + _mag_cal_quat(0) * _mag_cal_quat(1)),
			(_mag_cal_quat(0) * _mag_cal_quat(0) - _mag_cal_quat(1) * _mag_cal_quat(1) - _mag_cal_quat(2) * _mag_cal_quat(2) + _mag_cal_quat(3) * _mag_cal_quat(3))
		);

		// fuse accel data only if its norm is close to 1 g (reduces drift when vehicle picked up and moved).
		// Angular rate correction
		const float accel_fusion_gain = 0.4f;
		Vector3f corr;
		if (ok_to_align) {
			corr = (k % accel.normalized()) * accel_fusion_gain;
		}

		// Gyro bias estimation
		const float gyro_bias_gain = 0.0f;
		const float gyro_bias_limit = 0.05f;
		Vector3f gyro = _imu_sample_delayed.delta_ang / _imu_sample_delayed.delta_ang_dt;
		float spinRate = gyro.length();
		if (spinRate < 0.175f) {
			_mag_cal_gyro_bias -= corr * (gyro_bias_gain * _imu_sample_delayed.delta_ang_dt);

			for (int i = 0; i < 3; i++) {
				_mag_cal_gyro_bias(i) = math::constrain(_mag_cal_gyro_bias(i), -gyro_bias_limit, gyro_bias_limit);
			}
		}

		Vector3f rates = gyro - _mag_cal_gyro_bias;

		// Feed forward gyro
		corr += rates;

		// Apply correction to state
		_mag_cal_quat += _mag_cal_quat.derivative1(corr) * _imu_sample_delayed.delta_ang_dt;

		// Normalize quaternion
		_mag_cal_quat.normalize();

	}
}

void Ekf::sampleDataMagCal()
{
	// don't run unless mag data is fresh
	if (_mag_sample_index < _mag_cal_nsamples) {
		// save field, quaternion and time step to struct until we have sufficient coverage.

		// check if yaw rate and tilt is sufficient to perform calibration
		float yaw_rate;
		if (_imu_sample_delayed.delta_ang_dt > 0.0001f) {
		// apply imu bias corrections to sensor data
			Vector3f corrected_delta_ang = _imu_sample_delayed.delta_ang - _mag_cal_gyro_bias * _imu_sample_delayed.delta_ang_dt;
			Dcmf Teb(_mag_cal_quat);
			yaw_rate = Teb(2,0) * corrected_delta_ang(0)
					+ Teb(2,1) * corrected_delta_ang(1)
					+ Teb(2,2) * corrected_delta_ang(2);
			yaw_rate = yaw_rate / _imu_sample_delayed.delta_ang_dt;

			bool tilt_ok = Teb(2,2) > cosf(math::radians(45.0f));

			if (!_mag_cal_sampling_active && fabsf(yaw_rate) > math::radians(10.0f) && tilt_ok) {
				_mag_cal_sampling_active = true;
			} else if (_mag_cal_sampling_active && (fabsf(yaw_rate) < math::radians(5.0f) || !tilt_ok)) {
				_mag_cal_sampling_active = false;
			}
		} else {
			// invalid dt so can't proceed
			return;

		}

		// don't run if excessively tilted of if not rotating quickly enough
		if (!_mag_cal_sampling_active) {
			return;
		}

		// limit to run once per 8 degrees of yaw rotation and check for reversal of rotation
		Eulerf euler321(_mag_cal_quat);
		float yaw_delta = euler321(2) - _mag_bias_ekf_yaw_last;
		if (yaw_delta > M_PI_F) {
			yaw_delta -= M_TWOPI_F;
		} else if (yaw_delta < -M_PI_F) {
			yaw_delta += M_TWOPI_F;
		}
		if (fabsf(yaw_delta) < math::radians(8.0f) || (yaw_delta * (float)_mag_cal_direction) < -0.001f) {
			return;
		}
		_mag_bias_ekf_yaw_last = euler321(2);

		// calculate the time lapsed since the last sample
		float time_delta_sec =  1E-6f * (float)(_imu_sample_delayed.time_us - _mag_cal_sample_time_us);

		// reset the calibrator first time or if data sampling is interrupted for more than 10 seconds
		if (_mag_cal_sample_time_us == 0 || time_delta_sec > 10 ) {
			// reset covariance matrix
			memset(_mag_cov_mat, 0, sizeof(_mag_cov_mat));
			_mag_cov_mat[0][0] = sq(0.2f);
			_mag_cov_mat[1][1] = sq(0.2f);
			_mag_cov_mat[2][2] = sq(0.2f);
			_mag_cov_mat[3][3] = sq(0.2f);

			// reset states to zero
			_mag_cal_states.mag_bias(0) = 0.0f;
			_mag_cal_states.mag_bias(1) = 0.0f;
			_mag_cal_states.mag_bias(2) = 0.0f;
			_mag_cal_states.yaw_offset = 0.0f;

			// reset all counters
			_mag_cal_iteration_index = 0;
			_mag_sample_index = 0;
			_mag_cal_direction = 0;

			// record time to prevent reset repeating
			_mag_cal_sample_time_us = _imu_sample_delayed.time_us;

			return;
		}

		/*
		When storing the first sample in the data set we need to:
		- store a earth field vector for use by the EKF that is consistent with an initial magnetic heading error of zero
		- record the direction of rotation to enable detection  of any rotation reversal
		*/
		if (_mag_sample_index == 0) {
			// the time lapsed since previous sample is set to zero for the first sample in a data set
			time_delta_sec = 0.0f;

			// rotate the magnetometer measurements into earth frame
			Vector3f mag_earth_meas = _R_to_earth * _mag_sample_delayed.mag;

			// calulate how much we need to rotate the earth field to start with zero yaw error
			float decl_offset = atan2f(mag_earth_meas(1), mag_earth_meas(0)) - getMagDeclination();

			// get earth field from tables
			_mag_field_EF = getGeoMagNED();

			// rotate to match local yaw
			float mn = _mag_field_EF(0) * cosf(decl_offset) - _mag_field_EF(1) * sinf(decl_offset);
			float me = _mag_field_EF(1) * cosf(decl_offset) + _mag_field_EF(0) * sinf(decl_offset);
			_mag_field_EF(0) = mn;
			_mag_field_EF(1) = me;

			// record direction of yaw sample
			if (yaw_delta > 0.0f) {
				_mag_cal_direction = 1; // CW
			} else {
				_mag_cal_direction = -1; // CCW
			}
			printf("\nstarting sampling at %5.1f sec, yaw = %5.3f, dirn = %i\n",1.E-6*(double)_imu_sample_new.time_us,(double)euler321(2),_mag_cal_direction);
		}

		_mag_cal_fit_data[_mag_sample_index].mag_data = _mag_sample_delayed.mag;
		_mag_cal_fit_data[_mag_sample_index].quaternion = _mag_cal_quat;
		_mag_cal_fit_data[_mag_sample_index].time_step = time_delta_sec;
		_mag_cal_sample_time_us = _imu_sample_delayed.time_us;
		_mag_sample_index++;

	} else {
		_mag_cal_sampling_complete = true;
	}
}

void Ekf::processDataMagCal()
{

	// Process stored measurements

	// XYZ Measurement noise.
	float R_MAG = fmaxf(_params.mag_noise, 0.0f);
	R_MAG = R_MAG * R_MAG;

	// copy to variable names used by autocode
	// TODO remove
	float mn = _mag_field_EF(0);
	float me = _mag_field_EF(1);
	float md = _mag_field_EF(2);

	// Observation jacobian and Kalman gain vectors
	float H_MAG[4];
	float Kfusion[4];
	float innovation[3];

	float rss_innov[3] = {};

	for (uint8_t replay_index = 0; replay_index < 36; replay_index++) {

		// Apply process noise of 0.5 deg/sec to yaw state variance
		float yaw_process_noise_variance = _mag_cal_fit_data[replay_index].time_step * math::radians(0.5f);
		yaw_process_noise_variance *= yaw_process_noise_variance;
		_mag_cov_mat[3][3] += yaw_process_noise_variance;

		// update the states and covariance using sequential fusion of the magnetometer components
		for (uint8_t index = 0; index < 3; index++) {
			// rotate the quaternions by the yaw offset state
			Quatf quat_relative;
			quat_relative(0) = cosf(_mag_cal_states.yaw_offset);
			quat_relative(1) = 0.0f;
			quat_relative(2) = 0.0f;
			quat_relative(3) = sinf(_mag_cal_states.yaw_offset);
			quat_relative =  quat_relative * _mag_cal_fit_data[replay_index].quaternion;

			// get equivalent rotation matrix
			Matrix3f Teb = matrix::Dcmf(quat_relative).transpose();

			// rotate earth field into body frame and add bias states to get predicted measurement
			Vector3f mag_obs_predicted = Teb * _mag_field_EF + _mag_cal_states.mag_bias;

			// copy to variable names used by autocode
			// TODO remove
			float q0 = quat_relative(0);
			float q1 = quat_relative(1);
			float q2 = quat_relative(2);

			// intermediate variables from algebraic optimisation
			float t2 = cosf(_mag_cal_states.yaw_offset);
			float t3 = sinf(_mag_cal_states.yaw_offset);
			float t4 = q1*t2;
			float t5 = q0*t3;
			float t6 = t4+t5;
			float t7 = q2*t2;
			float t8 = q1*t3;
			float t9 = t7+t8;
			float t10 = q0*t2;
			float t15 = q2*t3;
			float t11 = t4-t15;
			float t12 = q0*t2*t9*2.0f;
			float t13 = t8-t10;
			float t14 = q0*t3*t13*2.0f;
			float t16 = q0*t3*t11*2.0f;
			float t17 = q0*q0;
			float t18 = t9*t11*2.0f;
			float t19 = q0*t2*t6*2.0f;
			float t20 = t6*t11*2.0f;
			float t21 = t2*t3*t17*4.0f;
			float t22 = t6*t13*2.0f;
			float t23 = t18+t22;
			float t24 = t12-t14+t16+t19;

			// Calculate observation jacobians and innovation
			if (index == 0) {
				// Calculate X axis observation jacobians
				memset(H_MAG, 0, sizeof(H_MAG));
				H_MAG[0] = 1.0f;
				H_MAG[1] = 0.0f;
				H_MAG[2] = 0.0f;
				H_MAG[3] = -md*(t12+t14+t16-q0*t2*t6*2.0f)-me*t24+mn*(t18+t6*(t10-q1*t3)*2.0f-t2*t3*t17*4.0f);

			} else if (index == 1) {
				// Calculate Y axis observation jacobians
				memset(H_MAG, 0, sizeof(H_MAG));
				H_MAG[0] = 0.0f;
				H_MAG[1] = 1.0f;
				H_MAG[2] = 0.0f;
				H_MAG[3] = me*t23-md*(t20+t21-t9*t13*2.0f)+mn*(t12+t14+t16-t19);

			} else if (index == 2) {
				// calculate Z axis observation jacobians
				memset(H_MAG, 0, sizeof(H_MAG));
				H_MAG[0] = 0.0f;
				H_MAG[1] = 0.0f;
				H_MAG[2] = 1.0f;
				H_MAG[3] = -md*t23-mn*t24+me*(-t20+t21+t9*t13*2.0f);

			}
			innovation[index] = mag_obs_predicted(index) - _mag_cal_fit_data[replay_index].mag_data(index);
			rss_innov[index] += sq(innovation[index]);

			// calculate the innovation variance
			float PH[4];
			float mag_innov_var = R_MAG;
			for (unsigned row = 0; row < 4; row++) {
				PH[row] = 0.0f;

				for (uint8_t col = 0; col < 4; col++) {
					PH[row] += _mag_cov_mat[row][col] * H_MAG[col];
				}

				mag_innov_var += H_MAG[row] * PH[row];
			}

			float mag_innov_var_inv;

			// check if the innovation variance calculation is badly conditioned
			if (mag_innov_var >= R_MAG) {
				// the innovation variance contribution from the state covariances is not negative, no fault
				mag_innov_var_inv = 1.0f / mag_innov_var;

			} else {
				// we reinitialise the covariance matrix and abort this fusion step
				memset(_mag_cov_mat, 0, sizeof(_mag_cov_mat));
				_mag_cov_mat[0][0] = 0.25f;
				_mag_cov_mat[1][1] = 0.25f;
				_mag_cov_mat[2][2] = 0.25f;
				_mag_cov_mat[3][3] = 1.0f;

				ECL_ERR("EKF mag bias cal fusion numerical error - covariance reset");

				return;

			}

			// calculate the Kalman gains
			for (uint8_t row = 0; row < 4; row++) {
				Kfusion[row] = 0.0f;

				for (uint8_t col = 0; col < 4; col++) {
					Kfusion[row] += _mag_cov_mat[row][col] * H_MAG[col];
				}

				Kfusion[row] *= mag_innov_var_inv;

			}

			// apply covariance correction via P_new = (I -K*H)*P
			// first calculate expression for KHP
			// then calculate P - KHP
			float KHP[4][4];
			float KH[4];

			for (unsigned row = 0; row < 4; row++) {

				KH[0] = Kfusion[row] * H_MAG[0];
				KH[1] = Kfusion[row] * H_MAG[1];
				KH[2] = Kfusion[row] * H_MAG[2];
				KH[3] = Kfusion[row] * H_MAG[3];

				for (unsigned col = 0; col < 4; col++) {
					float tmp = KH[0] * _mag_cov_mat[0][col];
					tmp += KH[1] * _mag_cov_mat[1][col];
					tmp += KH[2] * _mag_cov_mat[2][col];
					tmp += KH[3] * _mag_cov_mat[3][col];
					KHP[row][col] = tmp;
				}
			}

			// apply the covariance corrections
			for (unsigned row = 0; row < 4; row++) {
				for (unsigned col = 0; col < 4; col++) {
					_mag_cov_mat[row][col] = _mag_cov_mat[row][col] - KHP[row][col];
				}
			}

			// correct the covariance matrix for gross errors
			// force symmetry
			for (unsigned row = 0; row < 4; row++) {
				for (unsigned col = 0; col < row; col++) {
					float tmp = (_mag_cov_mat[row][col] + _mag_cov_mat[col][row]) / 2;
					_mag_cov_mat[row][col] = tmp;
					_mag_cov_mat[col][row] = tmp;
				}
			}
			// force positive variances
			for (unsigned col = 0; col < 2; col++) {
				_mag_cov_mat[col][col] = fmaxf(_mag_cov_mat[col][col], sq(0.005f));
			}
			_mag_cov_mat[3][3] = fmaxf(_mag_cov_mat[3][3], sq(0.005f));

			// apply the state corrections
			const float innov_limit = 0.5f;
			innovation[index] = math::constrain(innovation[index], -innov_limit, innov_limit);
			_mag_cal_states.mag_bias(0) -= Kfusion[0] * innovation[index];
			_mag_cal_states.mag_bias(1) -= Kfusion[1] * innovation[index];
			_mag_cal_states.mag_bias(2) -= Kfusion[2] * innovation[index];
			_mag_cal_states.yaw_offset -= Kfusion[3] * innovation[index];

			// Constrain state estimates
			const float bias_limit = 0.5f;
			_mag_cal_states.mag_bias(0) = math::constrain(_mag_cal_states.mag_bias(0), -bias_limit, bias_limit);
			_mag_cal_states.mag_bias(1) = math::constrain(_mag_cal_states.mag_bias(1), -bias_limit, bias_limit);
			_mag_cal_states.mag_bias(2) = math::constrain(_mag_cal_states.mag_bias(2), -bias_limit, bias_limit);
			_mag_cal_states.yaw_offset = math::constrain(_mag_cal_states.yaw_offset, -math::radians(180.0f), math::radians(180.0f));
		}
	}
	// get the RSS residiual for each sensor axis
	const float k1 = 1.0f / (float)_mag_cal_nsamples;
	rss_innov[0] = k1 * rss_innov[0];
	rss_innov[1] = k1 * rss_innov[1];
	rss_innov[2] = k1 * rss_innov[2];

	// Simple convergence test that runs a minimum of 5 and a maximum of 100 iterations
	// Assume converged if  worst single axis residual is less than 50 mGauss and the reduction in residual for all axes from the prrevious iteration is less than 0.1%
	// If convergence fails then rotate the earth field by 45 degrees and try again for a maximum of 7 retries.
	// This will give a maximum of 8 yaw starting points 45 degrees apart so that a bad initial yaw cguess aused by a large mag bias does not casue convergene failure
	bool solution_converged = rss_innov[0] / _mag_cal_residual[0] > 0.999f &&
					rss_innov[1] / _mag_cal_residual[1] > 0.999f &&
					rss_innov[2] / _mag_cal_residual[2] > 0.999f &&
					fmaxf(fmaxf(rss_innov[0],rss_innov[1]),rss_innov[2]) < 0.05f;

	if (_mag_cal_iteration_index > 4 && solution_converged) {
		_mag_cal_complete = true;
		_restart_sampling = true;

	} else if (_mag_cal_iteration_index >= 99) {
		if (fmaxf(fmaxf(rss_innov[0],rss_innov[1]),rss_innov[2]) >= 0.05f) {
			// rotate mag field by 45 degrees and try again
			if (_retry_count < 7) {
				printf("trying new earth field yaw offset\n");
				float mn_temp = _mag_field_EF(0) * cosf(math::radians(45.0f)) - _mag_field_EF(1) * sinf(math::radians(45.0f));
				float me_temp = _mag_field_EF(1) * cosf(math::radians(45.0f)) + _mag_field_EF(0) * sinf(math::radians(45.0f));
				_mag_field_EF(0) = mn_temp;
				_mag_field_EF(1) = me_temp;

				// restart iteration
				_mag_cal_iteration_index = 0;
				_retry_count ++;
				_mag_cal_complete = false;

				// reset covariance matrix
				memset(_mag_cov_mat, 0, sizeof(_mag_cov_mat));
				_mag_cov_mat[0][0] = sq(0.2f);
				_mag_cov_mat[1][1] = sq(0.2f);
				_mag_cov_mat[2][2] = sq(0.2f);
				_mag_cov_mat[3][3] = sq(0.2f);

				// reset states to zero
				_mag_cal_states.mag_bias(0) = 0.0f;
				_mag_cal_states.mag_bias(1) = 0.0f;
				_mag_cal_states.mag_bias(2) = 0.0f;
				_mag_cal_states.yaw_offset = 0.0f;

			} else {
				printf("convergence failed, residuals = %5.3f,%5.3f,%5.3f\n",(double)rss_innov[0],(double)rss_innov[1],(double)rss_innov[2]);
				_mag_cal_complete = false;
				_restart_sampling = true;

			}
		} else {
			printf("convergence was slow\n");
			_mag_cal_complete = true;
			_restart_sampling = true;

		}

	}

	// store residuals for use by convergence check next iteration
	for (uint8_t index = 0; index < 3; index++) {
		_mag_cal_residual[index] = rss_innov[index];
	}
	_mag_cal_iteration_index++;

}

// Return the magnetic field in Gauss to be used by  alignment and fusion processing
Vector3f Ekf::getGeoMagNED()
{
	// use parameter value until GPS is available, then use value returned by geo library
	if (_NED_origin_initialised) {
		// predicted earth field vector
		float mag_h = _mag_strength_gps * cosf(-_mag_inclination_gps);
		Vector3f mag_EF{cosf(_mag_declination_gps) * mag_h, sinf(_mag_declination_gps) * mag_h, _mag_strength_gps * sinf(-_mag_inclination_gps) };
		return mag_EF;

	} else {
		float mag_h = _params.mag_strength_gauss * cosf(math::radians(-_params.mag_inclination_deg));
		Vector3f mag_EF{cosf(math::radians(_params.mag_declination_deg)) * mag_h, sinf(math::radians(_params.mag_declination_deg)) * mag_h, _params.mag_strength_gauss * sinf(math::radians(-_params.mag_inclination_deg))};
		return mag_EF;
	}
}
