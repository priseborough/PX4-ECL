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
 * @file vel_pos_fusion.cpp
 * Function for fusion of direct observations of  velocity and position states/
 *
 * @author Paul Riseborough <p_riseborough@live.com.au>
 *
 */

#include "ukf.h"
#include "mathlib.h"

void Ukf::fuseVelPosHeight()
{
	bool fuse_map[6] = {}; // map of booleans true when [VN,VE,VD,PN,PE,PD] observations are available
	bool innov_check_pass_map[6] = {}; // true when innovations consistency checks pass for [VN,VE,VD,PN,PE,PD] observations
	float R[6] = {}; // observation variances for [VN,VE,VD,PN,PE,PD]
	float gate_size[6] = {}; // innovation consistency check gate sizes for [VN,VE,VD,PN,PE,PD] observations
	float Kfusion[UKF_N_STATES] = {}; // Kalman gain vector for any single observation - sequential fusion is used
	float innovation[6]; // local copy of innovations for  [VN,VE,VD,PN,PE,PD]
	memcpy(innovation, _vel_pos_innov, sizeof(_vel_pos_innov));

	// calculate innovations, innovations gate sizes and observation variances
	if (_fuse_hor_vel || _fuse_hor_vel_aux) {
		// enable fusion for NE velocity axes
		fuse_map[0] = fuse_map[1] = true;

		// handle special case where we are getting velocity observations from an auxiliary source
		if (!_fuse_hor_vel) {
			innovation[0] = _aux_vel_innov[0];
			innovation[1] = _aux_vel_innov[1];
		}

		// Set observation noise variance and innovation consistency check gate size for the NE position observations
		R[0] = _velObsVarNE(0);
		R[1] = _velObsVarNE(1);
		gate_size[1] = gate_size[0] = _hvelInnovGate;

	}

	if (_fuse_vert_vel) {
		fuse_map[2] = true;
		// observation variance - use receiver reported accuracy with parameter setting the minimum value
		R[2] = fmaxf(_params.gps_vel_noise, 0.01f);
		// use scaled horizontal speed accuracy assuming typical ratio of VDOP/HDOP
		R[2] = 1.5f * fmaxf(R[2], _gps_sample_delayed.sacc);
		R[2] = R[2] * R[2];
		// innovation gate size
		gate_size[2] = fmaxf(_params.vel_innov_gate, 1.0f);
	}

	if (_fuse_pos) {
		// enable fusion for the NE position axes
		fuse_map[3] = fuse_map[4] = true;

		// Set observation noise variance and innovation consistency check gate size for the NE position observations
		R[4] = R[3] = sq(_posObsNoiseNE);
		gate_size[4] = gate_size[3] = _posInnovGateNE;

	}

	if (_fuse_height) {
		if (_control_status.flags.baro_hgt) {
			fuse_map[5] = true;
			// vertical position innovation - baro measurement has opposite sign to earth z axis
			innovation[5] = _ukf_states.data.pos(2) + _baro_sample_delayed.hgt - _baro_hgt_offset - _hgt_sensor_offset;
			// observation variance - user parameter defined
			R[5] = fmaxf(_params.baro_noise, 0.01f);
			R[5] = R[5] * R[5];
			// innovation gate size
			gate_size[5] = fmaxf(_params.baro_innov_gate, 1.0f);

			// Compensate for positive static pressure transients (negative vertical position innovations)
			// casued by rotor wash ground interaction by applying a temporary deadzone to baro innovations.
			float deadzone_start = 0.25f * _params.baro_noise;
			float deadzone_end = deadzone_start + _params.gnd_effect_deadzone;
			if (_control_status.flags.gnd_effect) {
				if (innovation[5] < -deadzone_start) {
					if (innovation[5] <= -deadzone_end) {
						innovation[5] += deadzone_end;
					} else {
						innovation[5] = -deadzone_start;
					}
				}
			}

		} else if (_control_status.flags.gps_hgt) {
			fuse_map[5] = true;
			// vertical position innovation - gps measurement has opposite sign to earth z axis
			innovation[5] = _ukf_states.data.pos(2) + _gps_sample_delayed.hgt - _gps_alt_ref - _hgt_sensor_offset;
			// observation variance - receiver defined and parameter limited
			// use scaled horizontal position accuracy assuming typical ratio of VDOP/HDOP
			float lower_limit = fmaxf(_params.gps_pos_noise, 0.01f);
			float upper_limit = fmaxf(_params.pos_noaid_noise, lower_limit);
			R[5] = 1.5f * math::constrain(_gps_sample_delayed.vacc, lower_limit, upper_limit);
			R[5] = R[5] * R[5];
			// innovation gate size
			gate_size[5] = fmaxf(_params.baro_innov_gate, 1.0f);

		} else if (_control_status.flags.rng_hgt && (_R_rng_to_earth_2_2 > _params.range_cos_max_tilt)) {
			fuse_map[5] = true;
			// use range finder with tilt correction
			innovation[5] = _ukf_states.data.pos(2) - (-math::max(_range_sample_delayed.rng * _R_rng_to_earth_2_2,
							     _params.rng_gnd_clearance)) - _hgt_sensor_offset;
			// observation variance - user parameter defined
			R[5] = fmaxf((sq(_params.range_noise) + sq(_params.range_noise_scaler * _range_sample_delayed.rng)) * sq(_R_rng_to_earth_2_2), 0.01f);
			// innovation gate size
			gate_size[5] = fmaxf(_params.range_innov_gate, 1.0f);
		} else if (_control_status.flags.ev_hgt) {
			fuse_map[5] = true;
			// calculate the innovation assuming the external vision observaton is in local NED frame
			innovation[5] = _ukf_states.data.pos(2) - _ev_sample_delayed.posNED(2);
			// observation variance - defined externally
			R[5] = fmaxf(_ev_sample_delayed.posErr, 0.01f);
			R[5] = R[5] * R[5];
			// innovation gate size
			gate_size[5] = fmaxf(_params.ev_innov_gate, 1.0f);
		}

		// update innovation class variable for logging purposes
		_vel_pos_innov[5] = innovation[5];
	}

	if (_sigma_points_are_stale) {
		CalcSigmaPoints();
	}

	// Calculate covariance of predicted output and cross-covariance between state and output taking advantae of direct state observation
	float Pyy[6];
	float Pxy[6][UKF_N_STATES] = {};
	for (int obs_index=0; obs_index<6; obs_index++) { // loop through observations
		if (!fuse_map[obs_index]) {
			continue;
		}
		Pyy[obs_index] = R[obs_index];
		for (int sigma_index=0; sigma_index<UKF_N_SIGMA; sigma_index++) { // lopo through sigma points
			//Pyy +=  param.ukf.wc(s)*(psi_m(:,s) - y_m)*(psi_m(:,s) - y_m)';
			Pyy[obs_index] += _ukf_wc[sigma_index] * sq(_sigma_x_a(obs_index+3,sigma_index) - _sigma_x_a(obs_index+3,0));
			for (int state_index=0; state_index<UKF_N_STATES; state_index++) { // loop through states
				//Pxy += param.ukf.wc(s)*(sigma_x_a(1:param.ukf.nP,si) - x_m)*(psi_m(:,s) - y_m)';
				Pxy[obs_index][state_index] += _ukf_wc[sigma_index] * (_sigma_x_a(state_index,sigma_index) - _sigma_x_a(state_index,0)) * (_sigma_x_a(obs_index+3,sigma_index) - _sigma_x_a(obs_index+3,0));
			}
		}
	}

	// calculate innovation test ratios
	for (unsigned obs_index = 0; obs_index < 6; obs_index++) {
		if (fuse_map[obs_index]) {
			// compute the innovation variance SK = HPH + R
			unsigned state_index = obs_index + 3;	// we start with vx and this is the 3. state
			_vel_pos_innov_var[obs_index] = P_UKF(state_index,state_index) + R[obs_index];
			// Compute the ratio of innovation to gate size
			_vel_pos_test_ratio[obs_index] = sq(innovation[obs_index]) / (sq(gate_size[obs_index]) * Pyy[obs_index]);
		}
	}

	// check position, velocity and height innovations
	// treat 3D velocity, 2D position and height as separate sensors
	// always pass position checks if using synthetic position measurements or yet to complete tilt alignment
	// always pass height checks if yet to complete tilt alignment
	bool vel_check_pass = (_vel_pos_test_ratio[0] <= 1.0f) && (_vel_pos_test_ratio[1] <= 1.0f)
			      && (_vel_pos_test_ratio[2] <= 1.0f);
	innov_check_pass_map[2] = innov_check_pass_map[1] = innov_check_pass_map[0] = vel_check_pass;
	bool pos_check_pass = ((_vel_pos_test_ratio[3] <= 1.0f) && (_vel_pos_test_ratio[4] <= 1.0f)) || !_control_status.flags.tilt_align;
	innov_check_pass_map[4] = innov_check_pass_map[3] = pos_check_pass;
	innov_check_pass_map[5] = (_vel_pos_test_ratio[5] <= 1.0f) || !_control_status.flags.tilt_align;

	// record the successful velocity fusion event
	if ((_fuse_hor_vel || _fuse_hor_vel_aux) && vel_check_pass) {
		_time_last_vel_fuse = _time_last_imu;
		_innov_check_fail_status.flags.reject_vel_NED = false;
	} else if (!vel_check_pass) {
		_innov_check_fail_status.flags.reject_vel_NED = true;
	}
	_fuse_hor_vel = _fuse_hor_vel_aux = false;

	// record the successful position fusion event
	if (pos_check_pass && _fuse_pos) {
		if (!_fuse_hpos_as_odom) {
			_time_last_pos_fuse = _time_last_imu;
		} else {
			_time_last_delpos_fuse = _time_last_imu;
		}
		_innov_check_fail_status.flags.reject_pos_NE = false;
	} else if (!pos_check_pass) {
		_innov_check_fail_status.flags.reject_pos_NE = true;
	}
	_fuse_pos = false;

	// record the successful height fusion event
	if (innov_check_pass_map[5] && _fuse_height) {
		_time_last_hgt_fuse = _time_last_imu;
		_innov_check_fail_status.flags.reject_pos_D = false;
	} else if (!innov_check_pass_map[5]) {
		_innov_check_fail_status.flags.reject_pos_D = true;
	}
	_fuse_height = false;

	for (unsigned obs_index = 0; obs_index < 6; obs_index++) {
		// skip fusion if not requested or checks have failed
		if (!fuse_map[obs_index] || !innov_check_pass_map[obs_index]) {
			continue;
		}

		// calculate kalman gain
		for (int row = 0; row < UKF_N_STATES; row++) {
			Kfusion[row] = Pxy[obs_index][row] / Pyy[obs_index];
		}

		// update covariance matrix via P = P - K*Pyy*K'
		float KPK[UKF_N_STATES][UKF_N_STATES];
		for (unsigned row = 0; row < UKF_N_STATES; row++) {
			for (unsigned column = 0; column < UKF_N_STATES; column++) {
				KPK[row][column] = Kfusion[row] * Pyy[obs_index] * Kfusion[column];
			}
		}

		// if the covariance correction will result in a negative variance, then
		// the covariance matrix is unhealthy and must be corrected
		bool healthy = true;
		for (int i = 0; i < UKF_N_STATES; i++) {
			if (P_UKF(i,i) < KPK[i][i]) {
				// zero rows and columns
				zeroCovMat(i,i);

				//flag as unhealthy
				healthy = false;

				// update individual measurement health status
				if (obs_index == 0) {
					_fault_status.flags.bad_vel_N = true;
				} else if (obs_index == 1) {
					_fault_status.flags.bad_vel_E = true;
				} else if (obs_index == 2) {
					_fault_status.flags.bad_vel_D = true;
				} else if (obs_index == 3) {
					_fault_status.flags.bad_pos_N = true;
				} else if (obs_index == 4) {
					_fault_status.flags.bad_pos_E = true;
				} else if (obs_index == 5) {
					_fault_status.flags.bad_pos_D = true;
				}
			} else {
				// update individual measurement health status
				if (obs_index == 0) {
					_fault_status.flags.bad_vel_N = false;
				} else if (obs_index == 1) {
					_fault_status.flags.bad_vel_E = false;
				} else if (obs_index == 2) {
					_fault_status.flags.bad_vel_D = false;
				} else if (obs_index == 3) {
					_fault_status.flags.bad_pos_N = false;
				} else if (obs_index == 4) {
					_fault_status.flags.bad_pos_E = false;
				} else if (obs_index == 5) {
					_fault_status.flags.bad_pos_D = false;
				}
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
			fuse(Kfusion, innovation[obs_index]);
		}
	}
}

void Ukf::fuseHeight()
{
	float R_obs = 1.0f; // observation variances for PD
	float gate_size = 0.0f; // innovation consistency check gate size
	float Kfusion[UKF_N_STATES] = {}; // Kalman gain vector for any single observation - sequential fusion is used
	float innovation = 0.0f; // local copy of innovations for PD

	// calculate innovations, innovations gate sizes and observation variances
	if (_fuse_height) {
		if (_control_status.flags.baro_hgt) {
			// vertical position innovation - baro measurement has opposite sign to earth z axis
			innovation = _ukf_states.data.pos(2) + _baro_sample_delayed.hgt - _baro_hgt_offset - _hgt_sensor_offset;
			// observation variance - user parameter defined
			R_obs = fmaxf(_params.baro_noise, 0.01f);
			R_obs = R_obs * R_obs;
			// innovation gate size
			gate_size = fmaxf(_params.baro_innov_gate, 1.0f);

			// Compensate for positive static pressure transients (negative vertical position innovations)
			// casued by rotor wash ground interaction by applying a temporary deadzone to baro innovations.
			float deadzone_start = 0.25f * _params.baro_noise;
			float deadzone_end = deadzone_start + _params.gnd_effect_deadzone;
			if (_control_status.flags.gnd_effect) {
				if (innovation < -deadzone_start) {
					if (innovation <= -deadzone_end) {
						innovation += deadzone_end;
					} else {
						innovation = -deadzone_start;
					}
				}
			}

		} else if (_control_status.flags.gps_hgt) {
			// vertical position innovation - gps measurement has opposite sign to earth z axis
			innovation = _ukf_states.data.pos(2) + _gps_sample_delayed.hgt - _gps_alt_ref - _hgt_sensor_offset;
			// observation variance - receiver defined and parameter limited
			// use scaled horizontal position accuracy assuming typical ratio of VDOP/HDOP
			float lower_limit = fmaxf(_params.gps_pos_noise, 0.01f);
			float upper_limit = fmaxf(_params.pos_noaid_noise, lower_limit);
			R_obs = 1.5f * math::constrain(_gps_sample_delayed.vacc, lower_limit, upper_limit);
			R_obs = R_obs * R_obs;
			// innovation gate size
			gate_size = fmaxf(_params.baro_innov_gate, 1.0f);

		} else if (_control_status.flags.rng_hgt && (_R_rng_to_earth_2_2 > _params.range_cos_max_tilt)) {
			// use range finder with tilt correction
			innovation = _ukf_states.data.pos(2) - (-math::max(_range_sample_delayed.rng * _R_rng_to_earth_2_2,
							     _params.rng_gnd_clearance)) - _hgt_sensor_offset;
			// observation variance - user parameter defined
			R_obs = fmaxf((sq(_params.range_noise) + sq(_params.range_noise_scaler * _range_sample_delayed.rng)) * sq(_R_rng_to_earth_2_2), 0.01f);
			// innovation gate size
			gate_size = fmaxf(_params.range_innov_gate, 1.0f);
		} else if (_control_status.flags.ev_hgt) {
			// calculate the innovation assuming the external vision observaton is in local NED frame
			innovation = _ukf_states.data.pos(2) - _ev_sample_delayed.posNED(2);
			// observation variance - defined externally
			R_obs = fmaxf(_ev_sample_delayed.posErr, 0.01f);
			R_obs = R_obs * R_obs;
			// innovation gate size
			gate_size = fmaxf(_params.ev_innov_gate, 1.0f);
		}

		// update innovation class variable for logging purposes
		_vel_pos_innov[5] = innovation;
	}

	if (_sigma_points_are_stale) {
		CalcSigmaPoints();
	}

	// Calculate covariance of predicted output and cross-covariance between state and output taking advantage of direct state observation
	float Pyy;
	float Pxy[UKF_N_STATES] = {};
	Pyy = R_obs;
	for (int sigma_index=0; sigma_index<UKF_N_SIGMA; sigma_index++) { // lopo through sigma points
		//Pyy +=  param.ukf.wc(s)*(psi_m(:,s) - y_m)*(psi_m(:,s) - y_m)';
		Pyy += _ukf_wc[sigma_index] * sq(_sigma_x_a(8,sigma_index) - _sigma_x_a(8,0));
		for (int state_index=0; state_index<UKF_N_STATES; state_index++) { // loop through states
			//Pxy += param.ukf.wc(s)*(sigma_x_a(1:param.ukf.nP,si) - x_m)*(psi_m(:,s) - y_m)';
			Pxy[state_index] += _ukf_wc[sigma_index] * (_sigma_x_a(state_index,sigma_index) - _sigma_x_a(state_index,0)) * (_sigma_x_a(8,sigma_index) - _sigma_x_a(8,0));
		}
	}

	// calculate innovation test ratio
	// compute the innovation variance SK = HPH + R
	_vel_pos_innov_var[5] = P_UKF(8,8) + R_obs;
	// Compute the ratio of innovation to gate size
	_vel_pos_test_ratio[5] = sq(innovation) / (sq(gate_size) * Pyy);

	// check vertical position innovations
	// always pass height checks if yet to complete tilt alignment
	bool innov_check_pass = (_vel_pos_test_ratio[5] <= 1.0f) || !_control_status.flags.tilt_align;

	// record the successful height fusion event
	if (innov_check_pass && _fuse_height) {
		_time_last_hgt_fuse = _time_last_imu;
		_innov_check_fail_status.flags.reject_pos_D = false;
	} else if (!innov_check_pass) {
		_innov_check_fail_status.flags.reject_pos_D = true;
	}
	_fuse_height = false;

	// exit  if innovation checks have failed
	if (!innov_check_pass) {
		return;
	}

	// calculate kalman gain
	for (int row = 0; row < UKF_N_STATES; row++) {
		Kfusion[row] = Pxy[row] / Pyy;
	}

	// update covariance matrix via P = P - K*Pyy*K'
	float KPK[UKF_N_STATES][UKF_N_STATES];
	for (unsigned row = 0; row < UKF_N_STATES; row++) {
		for (unsigned column = 0; column < UKF_N_STATES; column++) {
			KPK[row][column] = Kfusion[row] * Pyy * Kfusion[column];
		}
	}

	// if the covariance correction will result in a negative variance, then
	// the covariance matrix is unhealthy and must be corrected
	bool healthy = true;
	for (int i = 0; i < UKF_N_STATES; i++) {
		if (P_UKF(i,i) < KPK[i][i]) {
			// zero rows and columns, record health status and exit
			zeroCovMat(i,i);
			healthy = false;
			_fault_status.flags.bad_pos_D = true;
			return;
		} else {
			_fault_status.flags.bad_pos_D = false;
		}
	}

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
	fuse(Kfusion, innovation);
}
