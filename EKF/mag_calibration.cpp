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

void Ekf::fuseMagCal()
{
	if (_mag_cal_complete) {
		return;
	}

	if (_mag_sample_index < 36 && fabsf(_mag_cal_yaw_delta_sum) < math::radians(360.0f)) {
		// save field, quaternion and time step to struct until we have 360 deg coverage.

		// check if yaw rate and tilt is sufficient to perform calibration
		float yaw_rate;
		if (_imu_sample_delayed.delta_ang_dt > 0.0001f) {
		// apply imu bias corrections to sensor data
			Vector3f corrected_delta_ang = _imu_sample_delayed.delta_ang - _state.gyro_bias;

			yaw_rate = _R_to_earth(2,0) * corrected_delta_ang(0)
					+ _R_to_earth(2,1) * corrected_delta_ang(1)
					+ _R_to_earth(2,2) * corrected_delta_ang(2);
			yaw_rate = yaw_rate / _imu_sample_delayed.delta_ang_dt;

			bool tilt_ok = _R_to_earth(2,2) > cosf(math::radians(45.0f));

			if (!_mag_cal_sampling_active && fabsf(yaw_rate) > math::radians(10.0f) && tilt_ok) {
				_mag_cal_sampling_active = true;
			} else if (_mag_cal_sampling_active && (fabsf(yaw_rate) < math::radians(5.0f) || !tilt_ok)) {
				_mag_cal_sampling_active = false;
			}

		} else {
			// invalid dt so can't proceed
			return;

		}

		// don't run if main filter is using the magnetomer or if excessively tilted of if not rotating quickly enough
		if (!_mag_use_inhibit || !_mag_cal_sampling_active) {
			return;
		}

		// limit to run once per 8 degrees of yaw rotation and check for reversal of rotation
		Eulerf euler321(_state.quat_nominal);
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
		_mag_cal_yaw_delta_sum += yaw_delta;

		// reset the calibrator first time or if data sampling is interrupted for more than 10 seconds
		float time_delta_sec =  1E-6f * (float)(_imu_sample_delayed.time_us - _mag_cal_sample_time_us);
		if (_mag_cal_sample_time_us == 0 || time_delta_sec > 10 ) {
			// reset covariance matrix
			memset(_mag_cov_mat, 0, sizeof(_mag_cov_mat));
			_mag_cov_mat[0][0] = sq(0.05f);
			_mag_cov_mat[1][1] = sq(0.05f);
			_mag_cov_mat[2][2] = sq(0.05f);
			_mag_cov_mat[3][3] = sq(0.1f);

			// reset states to zero
			_mag_cal_states.mag_bias(0) = 0.0f;
			_mag_cal_states.mag_bias(1) = 0.0f;
			_mag_cal_states.mag_bias(2) = 0.0f;
			_mag_cal_states.yaw_offset = 0.0f;

			// reset all counters
			_mag_cal_iteration_index = 0;
			_mag_sample_index = 0;
			_mag_cal_direction = 0;
			_mag_cal_yaw_delta_sum = 0.0f;

			// record time to prevent reset repeating
			_mag_cal_sample_time_us = _imu_sample_delayed.time_us;

			return;

		}

		if (_mag_sample_index == 0) {
			time_delta_sec = 0.0f;

			// rotate the magnetometer measurements into earth frame assuming a zero yaw angle
			Vector3f mag_earth_meas = _R_to_earth * _mag_sample_delayed.mag;

			// calulate how much we need to rotate the earth field to start with zero yaw error
			_mag_cal_decl_offset = atan2f(mag_earth_meas(1), mag_earth_meas(0)) - getMagDeclination();

			// get earth field from tables
			_mag_field_EF = getGeoMagNED();

			// rotate to match local yaw
			float mn = _mag_field_EF(0) * cosf(_mag_cal_decl_offset) - _mag_field_EF(1) * sinf(_mag_cal_decl_offset);
			float me = _mag_field_EF(1) * cosf(_mag_cal_decl_offset) + _mag_field_EF(0) * sinf(_mag_cal_decl_offset);
			_mag_field_EF(0) = mn;
			_mag_field_EF(1) = me;

			// record direction of yaw sample
			if (yaw_delta > 0.0f) {
				_mag_cal_direction = 1;
			} else {
				_mag_cal_direction = -1;
			}
			printf("start sampling at %5.1f sec, yaw = %5.3f, dirn = %i\n",1.E-6*(double)_imu_sample_new.time_us,(double)euler321(2),_mag_cal_direction);
		}

		_mag_cal_fit_data[_mag_sample_index].mag_data = _mag_sample_delayed.mag;
		_mag_cal_fit_data[_mag_sample_index].quaternion = _state.quat_nominal;
		_mag_cal_fit_data[_mag_sample_index].time_step = time_delta_sec;
		_mag_cal_sample_time_us = _imu_sample_delayed.time_us;
		_mag_sample_index++;

	} else {
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
		const float k1 = 1.0f / (float)_mag_sample_index;
		rss_innov[0] = k1 * rss_innov[0];
		rss_innov[1] = k1 * rss_innov[1];
		rss_innov[2] = k1 * rss_innov[2];

		// Simple convergence test that runs a minimum of 5 iterations and completes when reduction in residuals has fallen below 0.1%
		// up to a maximum of 50
		if (_mag_cal_iteration_index > 4) {
			_mag_cal_complete = rss_innov[0] / _mag_cal_residual[0] > 0.999f &&
						rss_innov[1] / _mag_cal_residual[1] > 0.999f &&
						rss_innov[2] / _mag_cal_residual[2] > 0.999f;
		} else if (_mag_cal_iteration_index > 49) {
			_mag_cal_complete = true;
			printf("calibration failed to converge\n");
		}
		for (uint8_t index = 0; index < 3; index++) {
			_mag_cal_residual[index] = rss_innov[index];
		}
		_mag_cal_iteration_index++;

		// replay debug print
		if (_mag_cal_complete) {
			printf("\ncompleted at %5.1f sec after %u iterations\n",1E-6*(double)_imu_sample_new.time_us,(unsigned)_mag_cal_iteration_index);
			printf("bias = %5.3f,%5.3f,%5.3f yaw = %5.3f,\n",
			(double)_mag_cal_states.mag_bias(0), (double)_mag_cal_states.mag_bias(1), (double)_mag_cal_states.mag_bias(2),
			(double)_mag_cal_states.yaw_offset);
			printf("residuals = %5.3f,%5.3f,%5.3f\n\n",(double)rss_innov[0],(double)rss_innov[1],(double)rss_innov[2]);
			_mag_cal_complete = false;
			_mag_cal_iteration_index = 0;
			_mag_sample_index = 0;
			_mag_cal_sample_time_us = 0;
			_mag_cal_direction = 0;
		}

		//Matrix3f R = matrix::Dcmf(_mag_cal_fit_data[replay_index].quaternion).transpose();
		//mag_obs_predicted = R * mag_EF;
		//printf("EF = %5.3f,%5.3f,%5.3f\n",(double)mag_EF(0),(double)mag_EF(1),(double)mag_EF(2));
		//printf("pred = %5.3f,%5.3f,%5.3f\n",(double)mag_obs_predicted(0),(double)mag_obs_predicted(1),(double)mag_obs_predicted(2));
		//printf("meas = %5.3f,%5.3f,%5.3f\n\n",(double)_mag_cal_fit_data[replay_index].mag_data(0),(double)_mag_cal_fit_data[replay_index].mag_data(1),(double)_mag_cal_fit_data[replay_index].mag_data(2));

	}
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
