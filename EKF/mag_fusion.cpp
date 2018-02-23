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
 * @file heading_fusion.cpp
 * Magnetometer fusion methods.
 *
 * @author Roman Bast <bapstroman@gmail.com>
 * @author Paul Riseborough <p_riseborough@live.com.au>
 *
 */
#include "../ecl.h"
#include "ekf.h"
#include "mathlib.h"

void Ekf::fuseMag()
{
	// TODO
}

void Ekf::fuseHeading()
{
	float R_YAW = 1.0f;
	Vector3f mag_earth_pred;

	// Calculate the observation variance
	if (_control_status.flags.mag_hdg) {
		// using magnetic heading tuning parameter
		R_YAW = sq(fmaxf(_params.mag_heading_noise, 1.0e-2f));
	} else if (_control_status.flags.ev_yaw) {
		// using error estimate from external vision data
		R_YAW = sq(fmaxf(_ev_sample_delayed.angErr, 1.0e-2f));
	} else {
		// there is no yaw observation
		return;
	}

	if (_sigma_points_are_stale) {
		CalcSigmaPoints();
	}

	float sigma_y_pred[UKF_N_SIGMA];
	float sigma_y_meas[UKF_N_SIGMA];
	// determine if a 321 or 312 Euler sequence is best
	if (fabsf(_R_to_earth(2, 0)) < fabsf(_R_to_earth(2, 1))) {

		for (uint8_t s=0; s<UKF_N_SIGMA; s++) {
			// rotate the magnetometer measurement into earth frame
			Eulerf euler321(_sigma_quat[s]);
			sigma_y_pred[s] = euler321(2); // we will need the predicted heading to calculate the innovation

			// calculate the observed yaw angle
			if (_control_status.flags.mag_hdg) {
				// Set the yaw angle to zero and rotate the measurements into earth frame using the zero yaw angle
				euler321(2) = 0.0f;
				Dcmf R_to_earth(euler321);

				// rotate the magnetometer measurements into earth frame using a zero yaw angle
				if (_control_status.flags.mag_3D) {
					// don't apply bias corrections if we are learning them
					mag_earth_pred = R_to_earth * _mag_sample_delayed.mag;
				} else {
					Vector3f mag_body;
					mag_body(0) = _mag_sample_delayed.mag(0) - _sigma_x_a(18,s);
					mag_body(1) = _mag_sample_delayed.mag(1) - _sigma_x_a(19,s);
					mag_body(2) = _mag_sample_delayed.mag(2) - _sigma_x_a(20,s);
					mag_earth_pred = R_to_earth * mag_body;
				}

				// the angle of the projection onto the horizontal gives the yaw angle
				sigma_y_meas[s] = -atan2f(mag_earth_pred(1), mag_earth_pred(0)) + _mag_declination;

			} else if (_control_status.flags.ev_yaw) {
				// calculate the yaw angle for a 321 sequence
				// Expressions obtained from yaw_input_321.c produced by https://github.com/PX4/ecl/blob/master/matlab/scripts/Inertial%20Nav%20EKF/quat2yaw321.m
				float Tbn_1_0 = 2.0f*(_ev_sample_delayed.quat(0)*_ev_sample_delayed.quat(3)+_ev_sample_delayed.quat(1)*_ev_sample_delayed.quat(2));
				float Tbn_0_0 = sq(_ev_sample_delayed.quat(0))+sq(_ev_sample_delayed.quat(1))-sq(_ev_sample_delayed.quat(2))-sq(_ev_sample_delayed.quat(3));
				sigma_y_meas[s] = atan2f(Tbn_1_0,Tbn_0_0);

			} else {
				// there is no yaw observation
				return;
			}
		}
	} else {

		for (uint8_t s=0; s<UKF_N_SIGMA; s++) {
			Dcmf R_to_earth = quat_to_invrotmat(_sigma_quat[s]);
			/* Calculate the 312 sequence euler angles that rotate from earth to body frame
			 * Derived from https://github.com/PX4/ecl/blob/master/matlab/scripts/Inertial%20Nav%20EKF/quat2yaw312.m
			 * Body to nav frame transformation using a yaw-roll-pitch rotation sequence is given by:
			 *
			[ cos(pitch)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw), -cos(roll)*sin(yaw), cos(yaw)*sin(pitch) + cos(pitch)*sin(roll)*sin(yaw)]
			[ cos(pitch)*sin(yaw) + cos(yaw)*sin(pitch)*sin(roll),  cos(roll)*cos(yaw), sin(pitch)*sin(yaw) - cos(pitch)*cos(yaw)*sin(roll)]
			[                               -cos(roll)*sin(pitch),           sin(roll),                                cos(pitch)*cos(roll)]
			*/
			float yaw = atan2f(-R_to_earth(0, 1), R_to_earth(1, 1)); // first rotation (yaw)
			float roll = asinf(R_to_earth(2, 1)); // second rotation (roll)
			float pitch = atan2f(-R_to_earth(2, 0), R_to_earth(2, 2)); // third rotation (pitch)

			sigma_y_pred[s] = yaw; // we will need the predicted heading to calculate the innovation

			// calculate the observed yaw angle
			if (_control_status.flags.mag_hdg) {
				// Set the first rotation (yaw) to zero and rotate the measurements into earth frame
				yaw = 0.0f;

				// Calculate the body to earth frame rotation matrix from the euler angles using a 312 rotation sequence
				// Equations from Tbn_312.c produced by https://github.com/PX4/ecl/blob/master/matlab/scripts/Inertial%20Nav%20EKF/quat2yaw312.m
				float sy = sinf(yaw);
				float cy = cosf(yaw);
				float sp = sinf(pitch);
				float cp = cosf(pitch);
				float sr = sinf(roll);
				float cr = cosf(roll);
				R_to_earth(0,0) = cy*cp-sy*sp*sr;
				R_to_earth(0,1) = -sy*cr;
				R_to_earth(0,2) = cy*sp+sy*cp*sr;
				R_to_earth(1,0) = sy*cp+cy*sp*sr;
				R_to_earth(1,1) = cy*cr;
				R_to_earth(1,2) = sy*sp-cy*cp*sr;
				R_to_earth(2,0) = -sp*cr;
				R_to_earth(2,1) = sr;
				R_to_earth(2,2) = cp*cr;

				// rotate the magnetometer measurements into earth frame using a zero yaw angle
				if (_control_status.flags.mag_3D) {
					// don't apply bias corrections if we are learning them
					mag_earth_pred = R_to_earth * _mag_sample_delayed.mag;
				} else {
					Vector3f mag_body;
					mag_body(0) = _mag_sample_delayed.mag(0) - _sigma_x_a(18,s);
					mag_body(1) = _mag_sample_delayed.mag(1) - _sigma_x_a(19,s);
					mag_body(2) = _mag_sample_delayed.mag(2) - _sigma_x_a(20,s);
					mag_earth_pred = R_to_earth * mag_body;
				}

				// the angle of the projection onto the horizontal gives the yaw angle
				sigma_y_meas[s] = -atan2f(mag_earth_pred(1), mag_earth_pred(0)) + _mag_declination;

			} else if (_control_status.flags.ev_yaw) {
				// calculate the yaw angle for a 312 sequence
				// Values from yaw_input_312.c file produced by https://github.com/PX4/ecl/blob/master/matlab/scripts/Inertial%20Nav%20EKF/quat2yaw312.m
				float Tbn_0_1_neg = 2.0f*(_ev_sample_delayed.quat(0)*_ev_sample_delayed.quat(3)-_ev_sample_delayed.quat(1)*_ev_sample_delayed.quat(2));
				float Tbn_1_1 = sq(_ev_sample_delayed.quat(0))-sq(_ev_sample_delayed.quat(1))+sq(_ev_sample_delayed.quat(2))-sq(_ev_sample_delayed.quat(3));
				sigma_y_meas[s] = atan2f(Tbn_0_1_neg,Tbn_1_1);

			} else {
				// there is no yaw observation
				return;

			}
		}
	}

	// Calculate covariance of predicted output and cross-covariance between state and output
	float Pyy;
	float Pxy[UKF_N_STATES] = {};
	Pyy = R_YAW;
	for (int sigma_index=0; sigma_index<UKF_N_SIGMA; sigma_index++) { // lopo through sigma points
		//Pyy +=  param.ukf.wc(s)*(psi_m(:,s) - y_m)*(psi_m(:,s) - y_m)';
		Pyy += _ukf_wc[sigma_index] * sq(sigma_y_pred[sigma_index] - sigma_y_pred[0]);
		for (int state_index=0; state_index<UKF_N_STATES; state_index++) { // loop through states
			//Pxy += param.ukf.wc(s)*(sigma_x_a(1:param.ukf.nP,si) - x_m)*(psi_m(:,s) - y_m)';
			Pxy[state_index] += _ukf_wc[sigma_index] * (_sigma_x_a(state_index,sigma_index) - _sigma_x_a(state_index,0)) * (sigma_y_pred[sigma_index] - sigma_y_pred[0]);
		}
	}

	// Check innovation variance
	_heading_innov_var = Pyy;
	float heading_innov_var_inv;
	// check if the innovation variance calculation is badly conditioned
	if (_heading_innov_var >= R_YAW) {
		// the innovation variance contribution from the state covariances is not negative, no fault
		_fault_status.flags.bad_mag_hdg = false;
		heading_innov_var_inv = 1.0f / _heading_innov_var;

	} else {
		// the innovation variance contribution from the state covariances is negative which means the covariance matrix is badly conditioned
		_fault_status.flags.bad_mag_hdg = true;

		// we reinitialise the covariance matrix and abort this fusion step
		initialiseCovariance();
		ECL_ERR("EKF mag yaw fusion numerical error - covariance reset");
		return;
	}

	// calculate the Kalman gains
	// only calculate gains for states we are using - don't modify magnetic field states
	float Kfusion[UKF_N_STATES] = {};
	for (uint8_t row = 0; row <= 14; row++) {
		Kfusion[row] = Pxy[row] / Pyy;
	}
	if (_control_status.flags.wind) {
		for (uint8_t row = 21; row <= 22; row++) {
			Kfusion[row] = Pxy[row] / Pyy;
		}
	}

	// wrap the heading to the interval between +-pi
	sigma_y_meas[0] = wrap_pi(sigma_y_meas[0]);

	// calculate the innovation
	_heading_innov = sigma_y_pred[0] - sigma_y_meas[0];

	// wrap the innovation to the interval between +-pi
	_heading_innov = wrap_pi(_heading_innov);

	// innovation test ratio
	_yaw_test_ratio = sq(_heading_innov) / (sq(math::max(_params.heading_innov_gate, 1.0f)) * _heading_innov_var);

	// we are no longer using 3-axis fusion so set the reported test levels to zero
	memset(_mag_test_ratio, 0, sizeof(_mag_test_ratio));

	// set the magnetometer unhealthy if the test fails
	if (_yaw_test_ratio > 1.0f) {
		_innov_check_fail_status.flags.reject_yaw = true;

		// if we are in air we don't want to fuse the measurement
		// we allow to use it when on the ground because the large innovation could be caused
		// by interference or a large initial gyro bias
		if (_control_status.flags.in_air) {
			return;

		} else {
			// constrain the innovation to the maximum set by the gate
			float gate_limit = sqrtf((sq(math::max(_params.heading_innov_gate, 1.0f)) * _heading_innov_var));
			_heading_innov = math::constrain(_heading_innov, -gate_limit, gate_limit);
		}

	} else {
		_innov_check_fail_status.flags.reject_yaw = false;
	}

	// update covariance matrix via P = P - K*Pyy*K'
	float KPK[UKF_N_STATES][UKF_N_STATES];
	for (unsigned row = 0; row < UKF_N_STATES; row++) {
		for (unsigned column = 0; column < UKF_N_STATES; column++) {
			KPK[row][column] = Kfusion[row] * Pyy * Kfusion[column];
		}
	}

	// if the covariance correction will result in a negative variance, then
	// the covariance marix is unhealthy and must be corrected
	bool healthy = true;
	_fault_status.flags.bad_mag_hdg = false;
	for (int i = 0; i < UKF_N_STATES; i++) {
		if (P_UKF(i,i) < KPK[i][i]) {
			// zero rows and columns
			zeroCovMat(i,i);

			//flag as unhealthy
			healthy = false;

			// update individual measurement health status
			_fault_status.flags.bad_mag_hdg = true;

		}
	}

	// only apply covariance and state corrrections if healthy
	if (healthy) {
		// apply the covariance corrections
		for (unsigned row = 0; row < UKF_N_STATES; row++) {
			for (unsigned column = 0; column < UKF_N_STATES; column++) {
				P_UKF(row,column) -= KPK[row][column];
			}
		}
		_sigma_points_are_stale = true;

		// correct the covariance marix for gross errors
		fixCovarianceErrors();

		// apply the state corrections
		fuse(Kfusion, _heading_innov);

	}
}

void Ekf::fuseDeclination()
{
	// TODO
}
