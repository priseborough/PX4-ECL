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

#include "ekf.h"

#include <ecl.h>
#include <math.h>
#include <mathlib/mathlib.h>

void Ekf::initialiseCovariance()
{
	for (unsigned i = 0; i < _k_num_states; i++) {
		for (unsigned j = 0; j < _k_num_states; j++) {
			P[i][j] = 0.0f;
		}
	}

	// calculate average prediction time step in sec
	float dt = FILTER_UPDATE_PERIOD_S;

	// define the initial angle uncertainty as variances for a rotation vector
	Vector3f rot_vec_var;
	rot_vec_var(2) = rot_vec_var(1) = rot_vec_var(0) = sq(_params.initial_tilt_err);

	// update the quaternion state covariances
	initialiseQuatCovariances(rot_vec_var);

	// velocity
	P[4][4] = sq(fmaxf(_params.gps_vel_noise, 0.01f));
	P[5][5] = P[4][4];
	P[6][6] = sq(1.5f) * P[4][4];

	// position
	P[7][7] = sq(fmaxf(_params.gps_pos_noise, 0.01f));
	P[8][8] = P[7][7];

	if (_control_status.flags.rng_hgt) {
		P[9][9] = sq(fmaxf(_params.range_noise, 0.01f));

	} else if (_control_status.flags.gps_hgt) {
		float lower_limit = fmaxf(_params.gps_pos_noise, 0.01f);
		float upper_limit = fmaxf(_params.pos_noaid_noise, lower_limit);
		P[9][9] = sq(1.5f * math::constrain(_gps_sample_delayed.vacc, lower_limit, upper_limit));

	} else {
		P[9][9] = sq(fmaxf(_params.baro_noise, 0.01f));
	}

	// gyro bias
	P[10][10] = sq(_params.switch_on_gyro_bias * dt);
	P[11][11] = P[10][10];
	P[12][12] = P[10][10];

	// accel bias
	_prev_dvel_bias_var(0) = P[13][13] = sq(_params.switch_on_accel_bias * dt);
	_prev_dvel_bias_var(1) = P[14][14] = P[13][13];
	_prev_dvel_bias_var(2) = P[15][15] = P[13][13];

	// record IMU bias state covariance reset time - used to prevent resets being performed too often
	_last_imu_bias_cov_reset_us = _imu_sample_delayed.time_us;

	// variances for optional states

	// earth frame and body frame magnetic field
	// set to observation variance
	for (uint8_t index = 16; index <= 21; index ++) {
		P[index][index] = sq(_params.mag_noise);
	}

	// save covariance data for re-use when auto-switching between heading and 3-axis fusion
	save_mag_cov_data();

	// wind
	P[22][22] = sq(_params.initial_wind_uncertainty);
	P[23][23] = sq(_params.initial_wind_uncertainty);

}

void Ekf::get_pos_var(Vector3f &pos_var)
{
	pos_var(0) = P[7][7];
	pos_var(1) = P[8][8];
	pos_var(2) = P[9][9];
}

void Ekf::get_vel_var(Vector3f &vel_var)
{
	vel_var(0) = P[4][4];
	vel_var(1) = P[5][5];
	vel_var(2) = P[6][6];
}

void Ekf::predictCovariance()
{
	// assign intermediate state variables
	float q0 = _state.quat_nominal(0);
	float q1 = _state.quat_nominal(1);
	float q2 = _state.quat_nominal(2);
	float q3 = _state.quat_nominal(3);

	float dax = _imu_sample_delayed.delta_ang(0);
	float day = _imu_sample_delayed.delta_ang(1);
	float daz = _imu_sample_delayed.delta_ang(2);

	float dvx = _imu_sample_delayed.delta_vel(0);
	float dvy = _imu_sample_delayed.delta_vel(1);
	float dvz = _imu_sample_delayed.delta_vel(2);

	float dax_b = _state.gyro_bias(0);
	float day_b = _state.gyro_bias(1);
	float daz_b = _state.gyro_bias(2);

	float dvx_b = _state.accel_bias(0);
	float dvy_b = _state.accel_bias(1);
	float dvz_b = _state.accel_bias(2);

	float dt = math::constrain(_imu_sample_delayed.delta_ang_dt, 0.5f * FILTER_UPDATE_PERIOD_S, 2.0f * FILTER_UPDATE_PERIOD_S);
	float dt_inv = 1.0f / dt;

	// compute noise variance for stationary processes
	float process_noise[_k_num_states] = {};

	// convert rate of change of rate gyro bias (rad/s**2) as specified by the parameter to an expected change in delta angle (rad) since the last update
	float d_ang_bias_sig = dt * dt * math::constrain(_params.gyro_bias_p_noise, 0.0f, 1.0f);

	// convert rate of change of accelerometer bias (m/s**3) as specified by the parameter to an expected change in delta velocity (m/s) since the last update
	float d_vel_bias_sig = dt * dt * math::constrain(_params.accel_bias_p_noise, 0.0f, 1.0f);

	// inhibit learning of imu accel bias if the manoeuvre levels are too high to protect against the effect of sensor nonlinearities or bad accel data is detected
	float alpha = math::constrain((dt / _params.acc_bias_learn_tc), 0.0f, 1.0f);
	float beta = 1.0f - alpha;
	_ang_rate_mag_filt = fmaxf(dt_inv * _imu_sample_delayed.delta_ang.norm(), beta * _ang_rate_mag_filt);
	_accel_mag_filt = fmaxf(dt_inv * _imu_sample_delayed.delta_vel.norm(), beta * _accel_mag_filt);
	_accel_vec_filt(0) = alpha * dt_inv * _imu_sample_delayed.delta_vel(0) + beta * _accel_vec_filt(0);
	_accel_vec_filt(1) = alpha * dt_inv * _imu_sample_delayed.delta_vel(1) + beta * _accel_vec_filt(1);
	_accel_vec_filt(2) = alpha * dt_inv * _imu_sample_delayed.delta_vel(2) + beta * _accel_vec_filt(2);

	if (_ang_rate_mag_filt > _params.acc_bias_learn_gyr_lim
	    || _accel_mag_filt > _params.acc_bias_learn_acc_lim
	    || _bad_vert_accel_detected) {

		// store the bias state variances to be reinstated later
		if (!_accel_bias_inhibit) {
			_prev_dvel_bias_var(0) = P[13][13];
			_prev_dvel_bias_var(1) = P[14][14];
			_prev_dvel_bias_var(2) = P[15][15];
		}

		_accel_bias_inhibit = true;

	} else {
		if (_accel_bias_inhibit) {
			// reinstate the bias state variances
			P[13][13] = _prev_dvel_bias_var(0);
			P[14][14] = _prev_dvel_bias_var(1);
			P[15][15] = _prev_dvel_bias_var(2);

		} else {
			// store the bias state variances to be reinstated later
			_prev_dvel_bias_var(0) = P[13][13];
			_prev_dvel_bias_var(1) = P[14][14];
			_prev_dvel_bias_var(2) = P[15][15];
		}

		_accel_bias_inhibit = false;
	}

	// Don't continue to grow the earth field variances if they are becoming too large or we are not doing 3-axis fusion as this can make the covariance matrix badly conditioned
	float mag_I_sig;

	if (_control_status.flags.mag_3D && (P[16][16] + P[17][17] + P[18][18]) < 0.1f) {
		mag_I_sig = dt * math::constrain(_params.mage_p_noise, 0.0f, 1.0f);

	} else {
		mag_I_sig = 0.0f;
	}

	// Don't continue to grow the body field variances if they is becoming too large or we are not doing 3-axis fusion as this can make the covariance matrix badly conditioned
	float mag_B_sig;

	if (_control_status.flags.mag_3D && (P[19][19] + P[20][20] + P[21][21]) < 0.1f) {
		mag_B_sig = dt * math::constrain(_params.magb_p_noise, 0.0f, 1.0f);

	} else {
		mag_B_sig = 0.0f;
	}

	float wind_vel_sig;

	// Calculate low pass filtered height rate
	float alpha_height_rate_lpf = 0.1f * dt; // 10 seconds time constant
	_height_rate_lpf = _height_rate_lpf * (1.0f - alpha_height_rate_lpf) + _state.vel(2) * alpha_height_rate_lpf;

	// Don't continue to grow wind velocity state variances if they are becoming too large or we are not using wind velocity states as this can make the covariance matrix badly conditioned
	if (_control_status.flags.wind && (P[22][22] + P[23][23]) < sq(_params.initial_wind_uncertainty)) {
		wind_vel_sig = dt * math::constrain(_params.wind_vel_p_noise, 0.0f, 1.0f) * (1.0f + _params.wind_vel_p_noise_scaler * fabsf(_height_rate_lpf));

	} else {
		wind_vel_sig = 0.0f;
	}

	// Construct the process noise variance diagonal for those states with a stationary process model
	// These are kinematic states and their error growth is controlled separately by the IMU noise variances
	for (unsigned i = 0; i <= 9; i++) {
		process_noise[i] = 0.0f;
	}

	// delta angle bias states
	process_noise[12] = process_noise[11] = process_noise[10] = sq(d_ang_bias_sig);
	// delta_velocity bias states
	process_noise[15] = process_noise[14] = process_noise[13] = sq(d_vel_bias_sig);
	// earth frame magnetic field states
	process_noise[18] = process_noise[17] = process_noise[16] = sq(mag_I_sig);
	// body frame magnetic field states
	process_noise[21] = process_noise[20] = process_noise[19] = sq(mag_B_sig);
	// wind velocity states
	process_noise[23] = process_noise[22] = sq(wind_vel_sig);

	// assign IMU noise variances
	// inputs to the system are 3 delta angles and 3 delta velocities
	float daxVar, dayVar, dazVar;
	float dvxVar, dvyVar, dvzVar;
	float gyro_noise = math::constrain(_params.gyro_noise, 0.0f, 1.0f);
	daxVar = dayVar = dazVar = sq(dt * gyro_noise);
	float accel_noise = math::constrain(_params.accel_noise, 0.0f, 1.0f);

	if (_bad_vert_accel_detected) {
		// Increase accelerometer process noise if bad accel data is detected. Measurement errors due to
		// vibration induced clipping commonly reach an equivalent 0.5g offset.
		accel_noise = BADACC_BIAS_PNOISE;
	}

	dvxVar = dvyVar = dvzVar = sq(dt * accel_noise);

	// predict the covariance

	// intermediate calculations
	const float PS0 = ecl::powf(q1, 2);
	const float PS1 = 0.25F*daxVar;
	const float PS2 = ecl::powf(q2, 2);
	const float PS3 = 0.25F*dayVar;
	const float PS4 = ecl::powf(q3, 2);
	const float PS5 = 0.25F*dazVar;
	const float PS6 = 0.5F*q1;
	const float PS7 = 0.5F*q2;
	const float PS8 = PS7*P[10][11];
	const float PS9 = 0.5F*q3;
	const float PS10 = PS9*P[10][12];
	const float PS11 = 0.5F*dax - 0.5F*dax_b;
	const float PS12 = 0.5F*day - 0.5F*day_b;
	const float PS13 = 0.5F*daz - 0.5F*daz_b;
	const float PS14 = PS10 - PS11*P[1][10] - PS12*P[2][10] - PS13*P[3][10] + PS6*P[10][10] + PS8 + P[0][10];
	const float PS15 = PS6*P[10][11];
	const float PS16 = PS9*P[11][12];
	const float PS17 = -PS11*P[1][11] - PS12*P[2][11] - PS13*P[3][11] + PS15 + PS16 + PS7*P[11][11] + P[0][11];
	const float PS18 = PS6*P[10][12];
	const float PS19 = PS7*P[11][12];
	const float PS20 = -PS11*P[1][12] - PS12*P[2][12] - PS13*P[3][12] + PS18 + PS19 + PS9*P[12][12] + P[0][12];
	const float PS21 = PS12*P[1][2];
	const float PS22 = -PS13*P[1][3];
	const float PS23 = -PS11*P[1][1] - PS21 + PS22 + PS6*P[1][10] + PS7*P[1][11] + PS9*P[1][12] + P[0][1];
	const float PS24 = -PS11*P[1][2];
	const float PS25 = PS13*P[2][3];
	const float PS26 = -PS12*P[2][2] + PS24 - PS25 + PS6*P[2][10] + PS7*P[2][11] + PS9*P[2][12] + P[0][2];
	const float PS27 = PS11*P[1][3];
	const float PS28 = -PS12*P[2][3];
	const float PS29 = -PS13*P[3][3] - PS27 + PS28 + PS6*P[3][10] + PS7*P[3][11] + PS9*P[3][12] + P[0][3];
	const float PS30 = PS11*P[0][1];
	const float PS31 = PS12*P[0][2];
	const float PS32 = PS13*P[0][3];
	const float PS33 = -PS30 - PS31 - PS32 + PS6*P[0][10] + PS7*P[0][11] + PS9*P[0][12] + P[0][0];
	const float PS34 = 0.5F*q0;
	const float PS35 = q2*q3;
	const float PS36 = q0*q1;
	const float PS37 = q1*q3;
	const float PS38 = q0*q2;
	const float PS39 = q1*q2;
	const float PS40 = q0*q3;
	const float PS41 = -PS2;
	const float PS42 = ecl::powf(q0, 2);
	const float PS43 = -PS4 + PS42;
	const float PS44 = PS0 + PS41 + PS43;
	const float PS45 = -PS11*P[1][13] - PS12*P[2][13] - PS13*P[3][13] + PS6*P[10][13] + PS7*P[11][13] + PS9*P[12][13] + P[0][13];
	const float PS46 = PS37 + PS38;
	const float PS47 = -PS11*P[1][15] - PS12*P[2][15] - PS13*P[3][15] + PS6*P[10][15] + PS7*P[11][15] + PS9*P[12][15] + P[0][15];
	const float PS48 = 2*PS47;
	const float PS49 = dvy - dvy_b;
	const float PS50 = dvx - dvx_b;
	const float PS51 = dvz - dvz_b;
	const float PS52 = PS49*q0 + PS50*q3 - PS51*q1;
	const float PS53 = 2*PS29;
	const float PS54 = -PS39 + PS40;
	const float PS55 = -PS11*P[1][14] - PS12*P[2][14] - PS13*P[3][14] + PS6*P[10][14] + PS7*P[11][14] + PS9*P[12][14] + P[0][14];
	const float PS56 = 2*PS55;
	const float PS57 = -PS49*q3 + PS50*q0 + PS51*q2;
	const float PS58 = 2*PS33;
	const float PS59 = PS49*q1 - PS50*q2 + PS51*q0;
	const float PS60 = 2*PS59;
	const float PS61 = PS49*q2 + PS50*q1 + PS51*q3;
	const float PS62 = 2*PS61;
	const float PS63 = -PS11*P[1][4] - PS12*P[2][4] - PS13*P[3][4] + PS6*P[4][10] + PS7*P[4][11] + PS9*P[4][12] + P[0][4];
	const float PS64 = -PS0;
	const float PS65 = PS2 + PS43 + PS64;
	const float PS66 = PS39 + PS40;
	const float PS67 = 2*PS45;
	const float PS68 = -PS35 + PS36;
	const float PS69 = -PS11*P[1][5] - PS12*P[2][5] - PS13*P[3][5] + PS6*P[5][10] + PS7*P[5][11] + PS9*P[5][12] + P[0][5];
	const float PS70 = PS4 + PS41 + PS42 + PS64;
	const float PS71 = PS35 + PS36;
	const float PS72 = 2*PS57;
	const float PS73 = -PS37 + PS38;
	const float PS74 = 2*PS52;
	const float PS75 = -PS11*P[1][6] - PS12*P[2][6] - PS13*P[3][6] + PS6*P[6][10] + PS7*P[6][11] + PS9*P[6][12] + P[0][6];
	const float PS76 = -PS34*P[10][11];
	const float PS77 = PS11*P[0][11] - PS12*P[3][11] + PS13*P[2][11] - PS19 + PS76 + PS9*P[11][11] + P[1][11];
	const float PS78 = PS13*P[0][2];
	const float PS79 = PS12*P[0][3];
	const float PS80 = PS11*P[0][0] - PS34*P[0][10] - PS7*P[0][12] + PS78 - PS79 + PS9*P[0][11] + P[0][1];
	const float PS81 = PS11*P[0][2];
	const float PS82 = PS13*P[2][2] + PS28 - PS34*P[2][10] - PS7*P[2][12] + PS81 + PS9*P[2][11] + P[1][2];
	const float PS83 = PS9*P[10][11];
	const float PS84 = PS7*P[10][12];
	const float PS85 = PS11*P[0][10] - PS12*P[3][10] + PS13*P[2][10] - PS34*P[10][10] + PS83 - PS84 + P[1][10];
	const float PS86 = -PS34*P[10][12];
	const float PS87 = PS11*P[0][12] - PS12*P[3][12] + PS13*P[2][12] + PS16 - PS7*P[12][12] + PS86 + P[1][12];
	const float PS88 = PS11*P[0][3];
	const float PS89 = -PS12*P[3][3] + PS25 - PS34*P[3][10] - PS7*P[3][12] + PS88 + PS9*P[3][11] + P[1][3];
	const float PS90 = PS13*P[1][2];
	const float PS91 = PS12*P[1][3];
	const float PS92 = PS30 - PS34*P[1][10] - PS7*P[1][12] + PS9*P[1][11] + PS90 - PS91 + P[1][1];
	const float PS93 = PS11*P[0][13] - PS12*P[3][13] + PS13*P[2][13] - PS34*P[10][13] - PS7*P[12][13] + PS9*P[11][13] + P[1][13];
	const float PS94 = PS11*P[0][15] - PS12*P[3][15] + PS13*P[2][15] - PS34*P[10][15] - PS7*P[12][15] + PS9*P[11][15] + P[1][15];
	const float PS95 = 2*PS94;
	const float PS96 = PS11*P[0][14] - PS12*P[3][14] + PS13*P[2][14] - PS34*P[10][14] - PS7*P[12][14] + PS9*P[11][14] + P[1][14];
	const float PS97 = 2*PS96;
	const float PS98 = PS11*P[0][4] - PS12*P[3][4] + PS13*P[2][4] - PS34*P[4][10] - PS7*P[4][12] + PS9*P[4][11] + P[1][4];
	const float PS99 = 2*PS93;
	const float PS100 = PS11*P[0][5] - PS12*P[3][5] + PS13*P[2][5] - PS34*P[5][10] - PS7*P[5][12] + PS9*P[5][11] + P[1][5];
	const float PS101 = PS11*P[0][6] - PS12*P[3][6] + PS13*P[2][6] - PS34*P[6][10] - PS7*P[6][12] + PS9*P[6][11] + P[1][6];
	const float PS102 = -PS34*P[11][12];
	const float PS103 = -PS10 + PS102 + PS11*P[3][12] + PS12*P[0][12] - PS13*P[1][12] + PS6*P[12][12] + P[2][12];
	const float PS104 = PS11*P[3][3] + PS22 - PS34*P[3][11] + PS6*P[3][12] + PS79 - PS9*P[3][10] + P[2][3];
	const float PS105 = PS13*P[0][1];
	const float PS106 = -PS105 + PS12*P[0][0] - PS34*P[0][11] + PS6*P[0][12] + PS88 - PS9*P[0][10] + P[0][2];
	const float PS107 = PS6*P[11][12];
	const float PS108 = PS107 + PS11*P[3][11] + PS12*P[0][11] - PS13*P[1][11] - PS34*P[11][11] - PS83 + P[2][11];
	const float PS109 = PS11*P[3][10] + PS12*P[0][10] - PS13*P[1][10] + PS18 + PS76 - PS9*P[10][10] + P[2][10];
	const float PS110 = PS12*P[0][1];
	const float PS111 = PS110 - PS13*P[1][1] + PS27 - PS34*P[1][11] + PS6*P[1][12] - PS9*P[1][10] + P[1][2];
	const float PS112 = PS11*P[2][3];
	const float PS113 = PS112 + PS31 - PS34*P[2][11] + PS6*P[2][12] - PS9*P[2][10] - PS90 + P[2][2];
	const float PS114 = PS11*P[3][13] + PS12*P[0][13] - PS13*P[1][13] - PS34*P[11][13] + PS6*P[12][13] - PS9*P[10][13] + P[2][13];
	const float PS115 = PS11*P[3][15] + PS12*P[0][15] - PS13*P[1][15] - PS34*P[11][15] + PS6*P[12][15] - PS9*P[10][15] + P[2][15];
	const float PS116 = 2*PS115;
	const float PS117 = PS11*P[3][14] + PS12*P[0][14] - PS13*P[1][14] - PS34*P[11][14] + PS6*P[12][14] - PS9*P[10][14] + P[2][14];
	const float PS118 = 2*PS117;
	const float PS119 = PS11*P[3][4] + PS12*P[0][4] - PS13*P[1][4] - PS34*P[4][11] + PS6*P[4][12] - PS9*P[4][10] + P[2][4];
	const float PS120 = 2*PS114;
	const float PS121 = PS11*P[3][5] + PS12*P[0][5] - PS13*P[1][5] - PS34*P[5][11] + PS6*P[5][12] - PS9*P[5][10] + P[2][5];
	const float PS122 = PS11*P[3][6] + PS12*P[0][6] - PS13*P[1][6] - PS34*P[6][11] + PS6*P[6][12] - PS9*P[6][10] + P[2][6];
	const float PS123 = -PS11*P[2][10] + PS12*P[1][10] + PS13*P[0][10] - PS15 + PS7*P[10][10] + PS86 + P[3][10];
	const float PS124 = PS105 + PS12*P[1][1] + PS24 - PS34*P[1][12] - PS6*P[1][11] + PS7*P[1][10] + P[1][3];
	const float PS125 = PS110 + PS13*P[0][0] - PS34*P[0][12] - PS6*P[0][11] + PS7*P[0][10] - PS81 + P[0][3];
	const float PS126 = -PS107 - PS11*P[2][12] + PS12*P[1][12] + PS13*P[0][12] - PS34*P[12][12] + PS84 + P[3][12];
	const float PS127 = PS102 - PS11*P[2][11] + PS12*P[1][11] + PS13*P[0][11] - PS6*P[11][11] + PS8 + P[3][11];
	const float PS128 = -PS11*P[2][2] + PS21 - PS34*P[2][12] - PS6*P[2][11] + PS7*P[2][10] + PS78 + P[2][3];
	const float PS129 = -PS112 + PS32 - PS34*P[3][12] - PS6*P[3][11] + PS7*P[3][10] + PS91 + P[3][3];
	const float PS130 = -PS11*P[2][13] + PS12*P[1][13] + PS13*P[0][13] - PS34*P[12][13] - PS6*P[11][13] + PS7*P[10][13] + P[3][13];
	const float PS131 = -PS11*P[2][15] + PS12*P[1][15] + PS13*P[0][15] - PS34*P[12][15] - PS6*P[11][15] + PS7*P[10][15] + P[3][15];
	const float PS132 = 2*PS131;
	const float PS133 = -PS11*P[2][14] + PS12*P[1][14] + PS13*P[0][14] - PS34*P[12][14] - PS6*P[11][14] + PS7*P[10][14] + P[3][14];
	const float PS134 = 2*PS133;
	const float PS135 = -PS11*P[2][4] + PS12*P[1][4] + PS13*P[0][4] - PS34*P[4][12] - PS6*P[4][11] + PS7*P[4][10] + P[3][4];
	const float PS136 = 2*PS130;
	const float PS137 = -PS11*P[2][5] + PS12*P[1][5] + PS13*P[0][5] - PS34*P[5][12] - PS6*P[5][11] + PS7*P[5][10] + P[3][5];
	const float PS138 = -PS11*P[2][6] + PS12*P[1][6] + PS13*P[0][6] - PS34*P[6][12] - PS6*P[6][11] + PS7*P[6][10] + P[3][6];
	const float PS139 = 2*PS46;
	const float PS140 = 2*PS54;
	const float PS141 = -PS139*P[13][15] + PS140*P[13][14] - PS44*P[13][13] + PS60*P[2][13] + PS62*P[1][13] + PS72*P[0][13] - PS74*P[3][13] + P[4][13];
	const float PS142 = -PS139*P[15][15] + PS140*P[14][15] - PS44*P[13][15] + PS60*P[2][15] + PS62*P[1][15] + PS72*P[0][15] - PS74*P[3][15] + P[4][15];
	const float PS143 = PS62*P[1][3];
	const float PS144 = PS72*P[0][3];
	const float PS145 = -PS139*P[3][15] + PS140*P[3][14] + PS143 + PS144 - PS44*P[3][13] + PS60*P[2][3] - PS74*P[3][3] + P[3][4];
	const float PS146 = -PS139*P[14][15] + PS140*P[14][14] - PS44*P[13][14] + PS60*P[2][14] + PS62*P[1][14] + PS72*P[0][14] - PS74*P[3][14] + P[4][14];
	const float PS147 = PS60*P[0][2];
	const float PS148 = PS74*P[0][3];
	const float PS149 = -PS139*P[0][15] + PS140*P[0][14] + PS147 - PS148 - PS44*P[0][13] + PS62*P[0][1] + PS72*P[0][0] + P[0][4];
	const float PS150 = PS62*P[1][2];
	const float PS151 = PS72*P[0][2];
	const float PS152 = -PS139*P[2][15] + PS140*P[2][14] + PS150 + PS151 - PS44*P[2][13] + PS60*P[2][2] - PS74*P[2][3] + P[2][4];
	const float PS153 = PS60*P[1][2];
	const float PS154 = PS74*P[1][3];
	const float PS155 = -PS139*P[1][15] + PS140*P[1][14] + PS153 - PS154 - PS44*P[1][13] + PS62*P[1][1] + PS72*P[0][1] + P[1][4];
	const float PS156 = 4*dvyVar;
	const float PS157 = 4*dvzVar;
	const float PS158 = -PS139*P[4][15] + PS140*P[4][14] - PS44*P[4][13] + PS60*P[2][4] + PS62*P[1][4] + PS72*P[0][4] - PS74*P[3][4] + P[4][4];
	const float PS159 = 2*PS141;
	const float PS160 = 2*PS68;
	const float PS161 = PS65*dvyVar;
	const float PS162 = 2*PS66;
	const float PS163 = PS44*dvxVar;
	const float PS164 = -PS139*P[5][15] + PS140*P[5][14] - PS44*P[5][13] + PS60*P[2][5] + PS62*P[1][5] + PS72*P[0][5] - PS74*P[3][5] + P[4][5];
	const float PS165 = 2*PS71;
	const float PS166 = 2*PS73;
	const float PS167 = PS70*dvzVar;
	const float PS168 = -PS139*P[6][15] + PS140*P[6][14] - PS44*P[6][13] + PS60*P[2][6] + PS62*P[1][6] + PS72*P[0][6] - PS74*P[3][6] + P[4][6];
	const float PS169 = PS160*P[14][15] - PS162*P[13][14] - PS60*P[1][14] + PS62*P[2][14] - PS65*P[14][14] + PS72*P[3][14] + PS74*P[0][14] + P[5][14];
	const float PS170 = PS160*P[13][15] - PS162*P[13][13] - PS60*P[1][13] + PS62*P[2][13] - PS65*P[13][14] + PS72*P[3][13] + PS74*P[0][13] + P[5][13];
	const float PS171 = PS74*P[0][1];
	const float PS172 = PS150 + PS160*P[1][15] - PS162*P[1][13] + PS171 - PS60*P[1][1] - PS65*P[1][14] + PS72*P[1][3] + P[1][5];
	const float PS173 = PS160*P[15][15] - PS162*P[13][15] - PS60*P[1][15] + PS62*P[2][15] - PS65*P[14][15] + PS72*P[3][15] + PS74*P[0][15] + P[5][15];
	const float PS174 = PS62*P[2][3];
	const float PS175 = PS148 + PS160*P[3][15] - PS162*P[3][13] + PS174 - PS60*P[1][3] - PS65*P[3][14] + PS72*P[3][3] + P[3][5];
	const float PS176 = PS60*P[0][1];
	const float PS177 = PS144 + PS160*P[0][15] - PS162*P[0][13] - PS176 + PS62*P[0][2] - PS65*P[0][14] + PS74*P[0][0] + P[0][5];
	const float PS178 = PS72*P[2][3];
	const float PS179 = -PS153 + PS160*P[2][15] - PS162*P[2][13] + PS178 + PS62*P[2][2] - PS65*P[2][14] + PS74*P[0][2] + P[2][5];
	const float PS180 = 4*dvxVar;
	const float PS181 = PS160*P[5][15] - PS162*P[5][13] - PS60*P[1][5] + PS62*P[2][5] - PS65*P[5][14] + PS72*P[3][5] + PS74*P[0][5] + P[5][5];
	const float PS182 = PS160*P[6][15] - PS162*P[6][13] - PS60*P[1][6] + PS62*P[2][6] - PS65*P[6][14] + PS72*P[3][6] + PS74*P[0][6] + P[5][6];
	const float PS183 = -PS165*P[14][15] + PS166*P[13][15] + PS60*P[0][15] + PS62*P[3][15] - PS70*P[15][15] - PS72*P[2][15] + PS74*P[1][15] + P[6][15];
	const float PS184 = -PS165*P[14][14] + PS166*P[13][14] + PS60*P[0][14] + PS62*P[3][14] - PS70*P[14][15] - PS72*P[2][14] + PS74*P[1][14] + P[6][14];
	const float PS185 = -PS165*P[13][14] + PS166*P[13][13] + PS60*P[0][13] + PS62*P[3][13] - PS70*P[13][15] - PS72*P[2][13] + PS74*P[1][13] + P[6][13];
	const float PS186 = -PS165*P[6][14] + PS166*P[6][13] + PS60*P[0][6] + PS62*P[3][6] - PS70*P[6][15] - PS72*P[2][6] + PS74*P[1][6] + P[6][6];

	// covariance update
	float nextP[24][24];

	// calculate variances and upper diagonal covariances for quaternion, velocity, position and gyro bias states
	nextP[0][0] = PS0*PS1 - PS11*PS23 - PS12*PS26 - PS13*PS29 + PS14*PS6 + PS17*PS7 + PS2*PS3 + PS20*PS9 + PS33 + PS4*PS5;
	nextP[0][1] = -PS1*PS36 + PS11*PS33 - PS12*PS29 + PS13*PS26 - PS14*PS34 + PS17*PS9 - PS20*PS7 + PS23 + PS3*PS35 - PS35*PS5;
	nextP[1][1] = PS1*PS42 + PS11*PS80 - PS12*PS89 + PS13*PS82 + PS2*PS5 + PS3*PS4 - PS34*PS85 - PS7*PS87 + PS77*PS9 + PS92;
	nextP[0][2] = -PS1*PS37 + PS11*PS29 + PS12*PS33 - PS13*PS23 - PS14*PS9 - PS17*PS34 + PS20*PS6 + PS26 - PS3*PS38 + PS37*PS5;
	nextP[1][2] = PS1*PS40 + PS11*PS89 + PS12*PS80 - PS13*PS92 - PS3*PS40 - PS34*PS77 - PS39*PS5 + PS6*PS87 + PS82 - PS85*PS9;
	nextP[2][2] = PS0*PS5 + PS1*PS4 + PS103*PS6 + PS104*PS11 + PS106*PS12 - PS108*PS34 - PS109*PS9 - PS111*PS13 + PS113 + PS3*PS42;
	nextP[0][3] = PS1*PS39 - PS11*PS26 + PS12*PS23 + PS13*PS33 + PS14*PS7 - PS17*PS6 - PS20*PS34 + PS29 - PS3*PS39 - PS40*PS5;
	nextP[1][3] = -PS1*PS38 - PS11*PS82 + PS12*PS92 + PS13*PS80 - PS3*PS37 - PS34*PS87 + PS38*PS5 - PS6*PS77 + PS7*PS85 + PS89;
	nextP[2][3] = -PS1*PS35 - PS103*PS34 + PS104 + PS106*PS13 - PS108*PS6 + PS109*PS7 - PS11*PS113 + PS111*PS12 + PS3*PS36 - PS36*PS5;
	nextP[3][3] = PS0*PS3 + PS1*PS2 - PS11*PS128 + PS12*PS124 + PS123*PS7 + PS125*PS13 - PS126*PS34 - PS127*PS6 + PS129 + PS42*PS5;
	nextP[0][4] = PS23*PS62 + PS26*PS60 - PS44*PS45 - PS46*PS48 - PS52*PS53 + PS54*PS56 + PS57*PS58 + PS63;
	nextP[1][4] = -PS44*PS93 - PS46*PS95 + PS54*PS97 + PS60*PS82 + PS62*PS92 + PS72*PS80 - PS74*PS89 + PS98;
	nextP[2][4] = -PS104*PS74 + PS106*PS72 + PS111*PS62 + PS113*PS60 - PS114*PS44 - PS116*PS46 + PS118*PS54 + PS119;
	nextP[3][4] = PS124*PS62 + PS125*PS72 + PS128*PS60 - PS129*PS74 - PS130*PS44 - PS132*PS46 + PS134*PS54 + PS135;
	nextP[4][4] = -PS139*PS142 + PS140*PS146 - PS141*PS44 - PS145*PS74 + PS149*PS72 + PS152*PS60 + PS155*PS62 + PS156*ecl::powf(PS54, 2) + PS157*ecl::powf(PS46, 2) + PS158 + ecl::powf(PS44, 2)*dvxVar;
	nextP[0][5] = -PS23*PS60 + PS26*PS62 + PS48*PS68 + PS52*PS58 + PS53*PS57 - PS55*PS65 - PS66*PS67 + PS69;
	nextP[1][5] = PS100 - PS60*PS92 + PS62*PS82 - PS65*PS96 - PS66*PS99 + PS68*PS95 + PS72*PS89 + PS74*PS80;
	nextP[2][5] = PS104*PS72 + PS106*PS74 - PS111*PS60 + PS113*PS62 + PS116*PS68 - PS117*PS65 - PS120*PS66 + PS121;
	nextP[3][5] = -PS124*PS60 + PS125*PS74 + PS128*PS62 + PS129*PS72 + PS132*PS68 - PS133*PS65 - PS136*PS66 + PS137;
	nextP[4][5] = -PS140*PS161 + PS142*PS160 + PS145*PS72 - PS146*PS65 + PS149*PS74 + PS152*PS62 - PS155*PS60 - PS157*PS46*PS68 - PS159*PS66 + PS162*PS163 + PS164;
	nextP[5][5] = PS157*ecl::powf(PS68, 2) + PS160*PS173 - PS162*PS170 - PS169*PS65 - PS172*PS60 + PS175*PS72 + PS177*PS74 + PS179*PS62 + PS180*ecl::powf(PS66, 2) + PS181 + ecl::powf(PS65, 2)*dvyVar;
	nextP[0][6] = PS23*PS74 - PS26*PS72 - PS47*PS70 + PS53*PS61 - PS56*PS71 + PS58*PS59 + PS67*PS73 + PS75;
	nextP[1][6] = PS101 + PS60*PS80 + PS62*PS89 - PS70*PS94 - PS71*PS97 - PS72*PS82 + PS73*PS99 + PS74*PS92;
	nextP[2][6] = PS104*PS62 + PS106*PS60 + PS111*PS74 - PS113*PS72 - PS115*PS70 - PS118*PS71 + PS120*PS73 + PS122;
	nextP[3][6] = PS124*PS74 + PS125*PS60 - PS128*PS72 + PS129*PS62 - PS131*PS70 - PS134*PS71 + PS136*PS73 + PS138;
	nextP[4][6] = PS139*PS167 - PS142*PS70 + PS145*PS62 - PS146*PS165 + PS149*PS60 - PS152*PS72 + PS155*PS74 - PS156*PS54*PS71 + PS159*PS73 - PS163*PS166 + PS168;
	nextP[5][6] = -PS160*PS167 + PS161*PS165 - PS165*PS169 + PS166*PS170 + PS172*PS74 - PS173*PS70 + PS175*PS62 + PS177*PS60 - PS179*PS72 - PS180*PS66*PS73 + PS182;
	nextP[6][6] = PS156*ecl::powf(PS71, 2) - PS165*PS184 + PS166*PS185 + PS180*ecl::powf(PS73, 2) - PS183*PS70 + PS186 + PS60*(-PS151 - PS165*P[0][14] + PS166*P[0][13] + PS171 + PS60*P[0][0] + PS62*P[0][3] - PS70*P[0][15] + P[0][6]) + PS62*(PS154 - PS165*P[3][14] + PS166*P[3][13] - PS178 + PS60*P[0][3] + PS62*P[3][3] - PS70*P[3][15] + P[3][6]) + ecl::powf(PS70, 2)*dvzVar - PS72*(PS147 - PS165*P[2][14] + PS166*P[2][13] + PS174 - PS70*P[2][15] - PS72*P[2][2] + PS74*P[1][2] + P[2][6]) + PS74*(PS143 - PS165*P[1][14] + PS166*P[1][13] + PS176 - PS70*P[1][15] - PS72*P[1][2] + PS74*P[1][1] + P[1][6]);
	nextP[0][7] = -PS11*P[1][7] - PS12*P[2][7] - PS13*P[3][7] + PS6*P[7][10] + PS63*dt + PS7*P[7][11] + PS9*P[7][12] + P[0][7];
	nextP[1][7] = PS11*P[0][7] - PS12*P[3][7] + PS13*P[2][7] - PS34*P[7][10] - PS7*P[7][12] + PS9*P[7][11] + PS98*dt + P[1][7];
	nextP[2][7] = PS11*P[3][7] + PS119*dt + PS12*P[0][7] - PS13*P[1][7] - PS34*P[7][11] + PS6*P[7][12] - PS9*P[7][10] + P[2][7];
	nextP[3][7] = -PS11*P[2][7] + PS12*P[1][7] + PS13*P[0][7] + PS135*dt - PS34*P[7][12] - PS6*P[7][11] + PS7*P[7][10] + P[3][7];
	nextP[4][7] = -PS139*P[7][15] + PS140*P[7][14] + PS158*dt - PS44*P[7][13] + PS60*P[2][7] + PS62*P[1][7] + PS72*P[0][7] - PS74*P[3][7] + P[4][7];
	nextP[5][7] = PS160*P[7][15] - PS162*P[7][13] - PS60*P[1][7] + PS62*P[2][7] - PS65*P[7][14] + PS72*P[3][7] + PS74*P[0][7] + P[5][7] + dt*(PS160*P[4][15] - PS162*P[4][13] - PS60*P[1][4] + PS62*P[2][4] - PS65*P[4][14] + PS72*P[3][4] + PS74*P[0][4] + P[4][5]);
	nextP[6][7] = -PS165*P[7][14] + PS166*P[7][13] + PS60*P[0][7] + PS62*P[3][7] - PS70*P[7][15] - PS72*P[2][7] + PS74*P[1][7] + P[6][7] + dt*(-PS165*P[4][14] + PS166*P[4][13] + PS60*P[0][4] + PS62*P[3][4] - PS70*P[4][15] - PS72*P[2][4] + PS74*P[1][4] + P[4][6]);
	nextP[7][7] = P[4][7]*dt + P[7][7] + dt*(P[4][4]*dt + P[4][7]);
	nextP[0][8] = -PS11*P[1][8] - PS12*P[2][8] - PS13*P[3][8] + PS6*P[8][10] + PS69*dt + PS7*P[8][11] + PS9*P[8][12] + P[0][8];
	nextP[1][8] = PS100*dt + PS11*P[0][8] - PS12*P[3][8] + PS13*P[2][8] - PS34*P[8][10] - PS7*P[8][12] + PS9*P[8][11] + P[1][8];
	nextP[2][8] = PS11*P[3][8] + PS12*P[0][8] + PS121*dt - PS13*P[1][8] - PS34*P[8][11] + PS6*P[8][12] - PS9*P[8][10] + P[2][8];
	nextP[3][8] = -PS11*P[2][8] + PS12*P[1][8] + PS13*P[0][8] + PS137*dt - PS34*P[8][12] - PS6*P[8][11] + PS7*P[8][10] + P[3][8];
	nextP[4][8] = -PS139*P[8][15] + PS140*P[8][14] + PS164*dt - PS44*P[8][13] + PS60*P[2][8] + PS62*P[1][8] + PS72*P[0][8] - PS74*P[3][8] + P[4][8];
	nextP[5][8] = PS160*P[8][15] - PS162*P[8][13] + PS181*dt - PS60*P[1][8] + PS62*P[2][8] - PS65*P[8][14] + PS72*P[3][8] + PS74*P[0][8] + P[5][8];
	nextP[6][8] = -PS165*P[8][14] + PS166*P[8][13] + PS60*P[0][8] + PS62*P[3][8] - PS70*P[8][15] - PS72*P[2][8] + PS74*P[1][8] + P[6][8] + dt*(-PS165*P[5][14] + PS166*P[5][13] + PS60*P[0][5] + PS62*P[3][5] - PS70*P[5][15] - PS72*P[2][5] + PS74*P[1][5] + P[5][6]);
	nextP[7][8] = P[4][8]*dt + P[7][8] + dt*(P[4][5]*dt + P[5][7]);
	nextP[8][8] = P[5][8]*dt + P[8][8] + dt*(P[5][5]*dt + P[5][8]);
	nextP[0][9] = -PS11*P[1][9] - PS12*P[2][9] - PS13*P[3][9] + PS6*P[9][10] + PS7*P[9][11] + PS75*dt + PS9*P[9][12] + P[0][9];
	nextP[1][9] = PS101*dt + PS11*P[0][9] - PS12*P[3][9] + PS13*P[2][9] - PS34*P[9][10] - PS7*P[9][12] + PS9*P[9][11] + P[1][9];
	nextP[2][9] = PS11*P[3][9] + PS12*P[0][9] + PS122*dt - PS13*P[1][9] - PS34*P[9][11] + PS6*P[9][12] - PS9*P[9][10] + P[2][9];
	nextP[3][9] = -PS11*P[2][9] + PS12*P[1][9] + PS13*P[0][9] + PS138*dt - PS34*P[9][12] - PS6*P[9][11] + PS7*P[9][10] + P[3][9];
	nextP[4][9] = -PS139*P[9][15] + PS140*P[9][14] + PS168*dt - PS44*P[9][13] + PS60*P[2][9] + PS62*P[1][9] + PS72*P[0][9] - PS74*P[3][9] + P[4][9];
	nextP[5][9] = PS160*P[9][15] - PS162*P[9][13] + PS182*dt - PS60*P[1][9] + PS62*P[2][9] - PS65*P[9][14] + PS72*P[3][9] + PS74*P[0][9] + P[5][9];
	nextP[6][9] = -PS165*P[9][14] + PS166*P[9][13] + PS186*dt + PS60*P[0][9] + PS62*P[3][9] - PS70*P[9][15] - PS72*P[2][9] + PS74*P[1][9] + P[6][9];
	nextP[7][9] = P[4][9]*dt + P[7][9] + dt*(P[4][6]*dt + P[6][7]);
	nextP[8][9] = P[5][9]*dt + P[8][9] + dt*(P[5][6]*dt + P[6][8]);
	nextP[9][9] = P[6][9]*dt + P[9][9] + dt*(P[6][6]*dt + P[6][9]);
	nextP[0][10] = PS14;
	nextP[1][10] = PS85;
	nextP[2][10] = PS109;
	nextP[3][10] = PS123;
	nextP[4][10] = -PS139*P[10][15] + PS140*P[10][14] - PS44*P[10][13] + PS60*P[2][10] + PS62*P[1][10] + PS72*P[0][10] - PS74*P[3][10] + P[4][10];
	nextP[5][10] = PS160*P[10][15] - PS162*P[10][13] - PS60*P[1][10] + PS62*P[2][10] - PS65*P[10][14] + PS72*P[3][10] + PS74*P[0][10] + P[5][10];
	nextP[6][10] = -PS165*P[10][14] + PS166*P[10][13] + PS60*P[0][10] + PS62*P[3][10] - PS70*P[10][15] - PS72*P[2][10] + PS74*P[1][10] + P[6][10];
	nextP[7][10] = P[4][10]*dt + P[7][10];
	nextP[8][10] = P[5][10]*dt + P[8][10];
	nextP[9][10] = P[6][10]*dt + P[9][10];
	nextP[10][10] = P[10][10];
	nextP[0][11] = PS17;
	nextP[1][11] = PS77;
	nextP[2][11] = PS108;
	nextP[3][11] = PS127;
	nextP[4][11] = -PS139*P[11][15] + PS140*P[11][14] - PS44*P[11][13] + PS60*P[2][11] + PS62*P[1][11] + PS72*P[0][11] - PS74*P[3][11] + P[4][11];
	nextP[5][11] = PS160*P[11][15] - PS162*P[11][13] - PS60*P[1][11] + PS62*P[2][11] - PS65*P[11][14] + PS72*P[3][11] + PS74*P[0][11] + P[5][11];
	nextP[6][11] = -PS165*P[11][14] + PS166*P[11][13] + PS60*P[0][11] + PS62*P[3][11] - PS70*P[11][15] - PS72*P[2][11] + PS74*P[1][11] + P[6][11];
	nextP[7][11] = P[4][11]*dt + P[7][11];
	nextP[8][11] = P[5][11]*dt + P[8][11];
	nextP[9][11] = P[6][11]*dt + P[9][11];
	nextP[10][11] = P[10][11];
	nextP[11][11] = P[11][11];
	nextP[0][12] = PS20;
	nextP[1][12] = PS87;
	nextP[2][12] = PS103;
	nextP[3][12] = PS126;
	nextP[4][12] = -PS139*P[12][15] + PS140*P[12][14] - PS44*P[12][13] + PS60*P[2][12] + PS62*P[1][12] + PS72*P[0][12] - PS74*P[3][12] + P[4][12];
	nextP[5][12] = PS160*P[12][15] - PS162*P[12][13] - PS60*P[1][12] + PS62*P[2][12] - PS65*P[12][14] + PS72*P[3][12] + PS74*P[0][12] + P[5][12];
	nextP[6][12] = -PS165*P[12][14] + PS166*P[12][13] + PS60*P[0][12] + PS62*P[3][12] - PS70*P[12][15] - PS72*P[2][12] + PS74*P[1][12] + P[6][12];
	nextP[7][12] = P[4][12]*dt + P[7][12];
	nextP[8][12] = P[5][12]*dt + P[8][12];
	nextP[9][12] = P[6][12]*dt + P[9][12];
	nextP[10][12] = P[10][12];
	nextP[11][12] = P[11][12];
	nextP[12][12] = P[12][12];

	// add process noise that is not from the IMU
	for (unsigned i = 0; i <= 12; i++) {
		nextP[i][i] += process_noise[i];
	}

	// Don't calculate these covariance terms if IMU delta velocity bias estimation is inhibited
	if (!(_params.fusion_mode & MASK_INHIBIT_ACC_BIAS) && !_accel_bias_inhibit) {

		// calculate variances and upper diagonal covariances for IMU delta velocity bias states
		nextP[0][13] = PS45;
		nextP[1][13] = PS93;
		nextP[2][13] = PS114;
		nextP[3][13] = PS130;
		nextP[4][13] = PS141;
		nextP[5][13] = PS170;
		nextP[6][13] = PS185;
		nextP[7][13] = P[4][13]*dt + P[7][13];
		nextP[8][13] = P[5][13]*dt + P[8][13];
		nextP[9][13] = P[6][13]*dt + P[9][13];
		nextP[10][13] = P[10][13];
		nextP[11][13] = P[11][13];
		nextP[12][13] = P[12][13];
		nextP[13][13] = P[13][13];
		nextP[0][14] = PS55;
		nextP[1][14] = PS96;
		nextP[2][14] = PS117;
		nextP[3][14] = PS133;
		nextP[4][14] = PS146;
		nextP[5][14] = PS169;
		nextP[6][14] = PS184;
		nextP[7][14] = P[4][14]*dt + P[7][14];
		nextP[8][14] = P[5][14]*dt + P[8][14];
		nextP[9][14] = P[6][14]*dt + P[9][14];
		nextP[10][14] = P[10][14];
		nextP[11][14] = P[11][14];
		nextP[12][14] = P[12][14];
		nextP[13][14] = P[13][14];
		nextP[14][14] = P[14][14];
		nextP[0][15] = PS47;
		nextP[1][15] = PS94;
		nextP[2][15] = PS115;
		nextP[3][15] = PS131;
		nextP[4][15] = PS142;
		nextP[5][15] = PS173;
		nextP[6][15] = PS183;
		nextP[7][15] = P[4][15]*dt + P[7][15];
		nextP[8][15] = P[5][15]*dt + P[8][15];
		nextP[9][15] = P[6][15]*dt + P[9][15];
		nextP[10][15] = P[10][15];
		nextP[11][15] = P[11][15];
		nextP[12][15] = P[12][15];
		nextP[13][15] = P[13][15];
		nextP[14][15] = P[14][15];
		nextP[15][15] = P[15][15];

		// add process noise that is not from the IMU
		for (unsigned i = 13; i <= 15; i++) {
			nextP[i][i] += process_noise[i];
		}

	} else {
		// Inhibit delta velocity bias learning by zeroing the covariance terms
		zeroRows(nextP, 13, 15);
		zeroCols(nextP, 13, 15);
	}

	// Don't do covariance prediction on magnetic field states unless we are using 3-axis fusion
	if (_control_status.flags.mag_3D) {
		// calculate variances and upper diagonal covariances for earth and body magnetic field states
		nextP[0][16] = -PS11*P[1][16] - PS12*P[2][16] - PS13*P[3][16] + PS6*P[10][16] + PS7*P[11][16] + PS9*P[12][16] + P[0][16];
		nextP[1][16] = PS11*P[0][16] - PS12*P[3][16] + PS13*P[2][16] - PS34*P[10][16] - PS7*P[12][16] + PS9*P[11][16] + P[1][16];
		nextP[2][16] = PS11*P[3][16] + PS12*P[0][16] - PS13*P[1][16] - PS34*P[11][16] + PS6*P[12][16] - PS9*P[10][16] + P[2][16];
		nextP[3][16] = -PS11*P[2][16] + PS12*P[1][16] + PS13*P[0][16] - PS34*P[12][16] - PS6*P[11][16] + PS7*P[10][16] + P[3][16];
		nextP[4][16] = -PS139*P[15][16] + PS140*P[14][16] - PS44*P[13][16] + PS60*P[2][16] + PS62*P[1][16] + PS72*P[0][16] - PS74*P[3][16] + P[4][16];
		nextP[5][16] = PS160*P[15][16] - PS162*P[13][16] - PS60*P[1][16] + PS62*P[2][16] - PS65*P[14][16] + PS72*P[3][16] + PS74*P[0][16] + P[5][16];
		nextP[6][16] = -PS165*P[14][16] + PS166*P[13][16] + PS60*P[0][16] + PS62*P[3][16] - PS70*P[15][16] - PS72*P[2][16] + PS74*P[1][16] + P[6][16];
		nextP[7][16] = P[4][16]*dt + P[7][16];
		nextP[8][16] = P[5][16]*dt + P[8][16];
		nextP[9][16] = P[6][16]*dt + P[9][16];
		nextP[10][16] = P[10][16];
		nextP[11][16] = P[11][16];
		nextP[12][16] = P[12][16];
		nextP[13][16] = P[13][16];
		nextP[14][16] = P[14][16];
		nextP[15][16] = P[15][16];
		nextP[16][16] = P[16][16];
		nextP[0][17] = -PS11*P[1][17] - PS12*P[2][17] - PS13*P[3][17] + PS6*P[10][17] + PS7*P[11][17] + PS9*P[12][17] + P[0][17];
		nextP[1][17] = PS11*P[0][17] - PS12*P[3][17] + PS13*P[2][17] - PS34*P[10][17] - PS7*P[12][17] + PS9*P[11][17] + P[1][17];
		nextP[2][17] = PS11*P[3][17] + PS12*P[0][17] - PS13*P[1][17] - PS34*P[11][17] + PS6*P[12][17] - PS9*P[10][17] + P[2][17];
		nextP[3][17] = -PS11*P[2][17] + PS12*P[1][17] + PS13*P[0][17] - PS34*P[12][17] - PS6*P[11][17] + PS7*P[10][17] + P[3][17];
		nextP[4][17] = -PS139*P[15][17] + PS140*P[14][17] - PS44*P[13][17] + PS60*P[2][17] + PS62*P[1][17] + PS72*P[0][17] - PS74*P[3][17] + P[4][17];
		nextP[5][17] = PS160*P[15][17] - PS162*P[13][17] - PS60*P[1][17] + PS62*P[2][17] - PS65*P[14][17] + PS72*P[3][17] + PS74*P[0][17] + P[5][17];
		nextP[6][17] = -PS165*P[14][17] + PS166*P[13][17] + PS60*P[0][17] + PS62*P[3][17] - PS70*P[15][17] - PS72*P[2][17] + PS74*P[1][17] + P[6][17];
		nextP[7][17] = P[4][17]*dt + P[7][17];
		nextP[8][17] = P[5][17]*dt + P[8][17];
		nextP[9][17] = P[6][17]*dt + P[9][17];
		nextP[10][17] = P[10][17];
		nextP[11][17] = P[11][17];
		nextP[12][17] = P[12][17];
		nextP[13][17] = P[13][17];
		nextP[14][17] = P[14][17];
		nextP[15][17] = P[15][17];
		nextP[16][17] = P[16][17];
		nextP[17][17] = P[17][17];
		nextP[0][18] = -PS11*P[1][18] - PS12*P[2][18] - PS13*P[3][18] + PS6*P[10][18] + PS7*P[11][18] + PS9*P[12][18] + P[0][18];
		nextP[1][18] = PS11*P[0][18] - PS12*P[3][18] + PS13*P[2][18] - PS34*P[10][18] - PS7*P[12][18] + PS9*P[11][18] + P[1][18];
		nextP[2][18] = PS11*P[3][18] + PS12*P[0][18] - PS13*P[1][18] - PS34*P[11][18] + PS6*P[12][18] - PS9*P[10][18] + P[2][18];
		nextP[3][18] = -PS11*P[2][18] + PS12*P[1][18] + PS13*P[0][18] - PS34*P[12][18] - PS6*P[11][18] + PS7*P[10][18] + P[3][18];
		nextP[4][18] = -PS139*P[15][18] + PS140*P[14][18] - PS44*P[13][18] + PS60*P[2][18] + PS62*P[1][18] + PS72*P[0][18] - PS74*P[3][18] + P[4][18];
		nextP[5][18] = PS160*P[15][18] - PS162*P[13][18] - PS60*P[1][18] + PS62*P[2][18] - PS65*P[14][18] + PS72*P[3][18] + PS74*P[0][18] + P[5][18];
		nextP[6][18] = -PS165*P[14][18] + PS166*P[13][18] + PS60*P[0][18] + PS62*P[3][18] - PS70*P[15][18] - PS72*P[2][18] + PS74*P[1][18] + P[6][18];
		nextP[7][18] = P[4][18]*dt + P[7][18];
		nextP[8][18] = P[5][18]*dt + P[8][18];
		nextP[9][18] = P[6][18]*dt + P[9][18];
		nextP[10][18] = P[10][18];
		nextP[11][18] = P[11][18];
		nextP[12][18] = P[12][18];
		nextP[13][18] = P[13][18];
		nextP[14][18] = P[14][18];
		nextP[15][18] = P[15][18];
		nextP[16][18] = P[16][18];
		nextP[17][18] = P[17][18];
		nextP[18][18] = P[18][18];
		nextP[0][19] = -PS11*P[1][19] - PS12*P[2][19] - PS13*P[3][19] + PS6*P[10][19] + PS7*P[11][19] + PS9*P[12][19] + P[0][19];
		nextP[1][19] = PS11*P[0][19] - PS12*P[3][19] + PS13*P[2][19] - PS34*P[10][19] - PS7*P[12][19] + PS9*P[11][19] + P[1][19];
		nextP[2][19] = PS11*P[3][19] + PS12*P[0][19] - PS13*P[1][19] - PS34*P[11][19] + PS6*P[12][19] - PS9*P[10][19] + P[2][19];
		nextP[3][19] = -PS11*P[2][19] + PS12*P[1][19] + PS13*P[0][19] - PS34*P[12][19] - PS6*P[11][19] + PS7*P[10][19] + P[3][19];
		nextP[4][19] = -PS139*P[15][19] + PS140*P[14][19] - PS44*P[13][19] + PS60*P[2][19] + PS62*P[1][19] + PS72*P[0][19] - PS74*P[3][19] + P[4][19];
		nextP[5][19] = PS160*P[15][19] - PS162*P[13][19] - PS60*P[1][19] + PS62*P[2][19] - PS65*P[14][19] + PS72*P[3][19] + PS74*P[0][19] + P[5][19];
		nextP[6][19] = -PS165*P[14][19] + PS166*P[13][19] + PS60*P[0][19] + PS62*P[3][19] - PS70*P[15][19] - PS72*P[2][19] + PS74*P[1][19] + P[6][19];
		nextP[7][19] = P[4][19]*dt + P[7][19];
		nextP[8][19] = P[5][19]*dt + P[8][19];
		nextP[9][19] = P[6][19]*dt + P[9][19];
		nextP[10][19] = P[10][19];
		nextP[11][19] = P[11][19];
		nextP[12][19] = P[12][19];
		nextP[13][19] = P[13][19];
		nextP[14][19] = P[14][19];
		nextP[15][19] = P[15][19];
		nextP[16][19] = P[16][19];
		nextP[17][19] = P[17][19];
		nextP[18][19] = P[18][19];
		nextP[19][19] = P[19][19];
		nextP[0][20] = -PS11*P[1][20] - PS12*P[2][20] - PS13*P[3][20] + PS6*P[10][20] + PS7*P[11][20] + PS9*P[12][20] + P[0][20];
		nextP[1][20] = PS11*P[0][20] - PS12*P[3][20] + PS13*P[2][20] - PS34*P[10][20] - PS7*P[12][20] + PS9*P[11][20] + P[1][20];
		nextP[2][20] = PS11*P[3][20] + PS12*P[0][20] - PS13*P[1][20] - PS34*P[11][20] + PS6*P[12][20] - PS9*P[10][20] + P[2][20];
		nextP[3][20] = -PS11*P[2][20] + PS12*P[1][20] + PS13*P[0][20] - PS34*P[12][20] - PS6*P[11][20] + PS7*P[10][20] + P[3][20];
		nextP[4][20] = -PS139*P[15][20] + PS140*P[14][20] - PS44*P[13][20] + PS60*P[2][20] + PS62*P[1][20] + PS72*P[0][20] - PS74*P[3][20] + P[4][20];
		nextP[5][20] = PS160*P[15][20] - PS162*P[13][20] - PS60*P[1][20] + PS62*P[2][20] - PS65*P[14][20] + PS72*P[3][20] + PS74*P[0][20] + P[5][20];
		nextP[6][20] = -PS165*P[14][20] + PS166*P[13][20] + PS60*P[0][20] + PS62*P[3][20] - PS70*P[15][20] - PS72*P[2][20] + PS74*P[1][20] + P[6][20];
		nextP[7][20] = P[4][20]*dt + P[7][20];
		nextP[8][20] = P[5][20]*dt + P[8][20];
		nextP[9][20] = P[6][20]*dt + P[9][20];
		nextP[10][20] = P[10][20];
		nextP[11][20] = P[11][20];
		nextP[12][20] = P[12][20];
		nextP[13][20] = P[13][20];
		nextP[14][20] = P[14][20];
		nextP[15][20] = P[15][20];
		nextP[16][20] = P[16][20];
		nextP[17][20] = P[17][20];
		nextP[18][20] = P[18][20];
		nextP[19][20] = P[19][20];
		nextP[20][20] = P[20][20];
		nextP[0][21] = -PS11*P[1][21] - PS12*P[2][21] - PS13*P[3][21] + PS6*P[10][21] + PS7*P[11][21] + PS9*P[12][21] + P[0][21];
		nextP[1][21] = PS11*P[0][21] - PS12*P[3][21] + PS13*P[2][21] - PS34*P[10][21] - PS7*P[12][21] + PS9*P[11][21] + P[1][21];
		nextP[2][21] = PS11*P[3][21] + PS12*P[0][21] - PS13*P[1][21] - PS34*P[11][21] + PS6*P[12][21] - PS9*P[10][21] + P[2][21];
		nextP[3][21] = -PS11*P[2][21] + PS12*P[1][21] + PS13*P[0][21] - PS34*P[12][21] - PS6*P[11][21] + PS7*P[10][21] + P[3][21];
		nextP[4][21] = -PS139*P[15][21] + PS140*P[14][21] - PS44*P[13][21] + PS60*P[2][21] + PS62*P[1][21] + PS72*P[0][21] - PS74*P[3][21] + P[4][21];
		nextP[5][21] = PS160*P[15][21] - PS162*P[13][21] - PS60*P[1][21] + PS62*P[2][21] - PS65*P[14][21] + PS72*P[3][21] + PS74*P[0][21] + P[5][21];
		nextP[6][21] = -PS165*P[14][21] + PS166*P[13][21] + PS60*P[0][21] + PS62*P[3][21] - PS70*P[15][21] - PS72*P[2][21] + PS74*P[1][21] + P[6][21];
		nextP[7][21] = P[4][21]*dt + P[7][21];
		nextP[8][21] = P[5][21]*dt + P[8][21];
		nextP[9][21] = P[6][21]*dt + P[9][21];
		nextP[10][21] = P[10][21];
		nextP[11][21] = P[11][21];
		nextP[12][21] = P[12][21];
		nextP[13][21] = P[13][21];
		nextP[14][21] = P[14][21];
		nextP[15][21] = P[15][21];
		nextP[16][21] = P[16][21];
		nextP[17][21] = P[17][21];
		nextP[18][21] = P[18][21];
		nextP[19][21] = P[19][21];
		nextP[20][21] = P[20][21];
		nextP[21][21] = P[21][21];

		// add process noise that is not from the IMU
		for (unsigned i = 16; i <= 21; i++) {
			nextP[i][i] += process_noise[i];
		}

	}

	// Don't do covariance prediction on wind states unless we are using them
	if (_control_status.flags.wind) {

		// calculate variances and upper diagonal covariances for wind states
		nextP[0][22] = -PS11*P[1][22] - PS12*P[2][22] - PS13*P[3][22] + PS6*P[10][22] + PS7*P[11][22] + PS9*P[12][22] + P[0][22];
		nextP[1][22] = PS11*P[0][22] - PS12*P[3][22] + PS13*P[2][22] - PS34*P[10][22] - PS7*P[12][22] + PS9*P[11][22] + P[1][22];
		nextP[2][22] = PS11*P[3][22] + PS12*P[0][22] - PS13*P[1][22] - PS34*P[11][22] + PS6*P[12][22] - PS9*P[10][22] + P[2][22];
		nextP[3][22] = -PS11*P[2][22] + PS12*P[1][22] + PS13*P[0][22] - PS34*P[12][22] - PS6*P[11][22] + PS7*P[10][22] + P[3][22];
		nextP[4][22] = -PS139*P[15][22] + PS140*P[14][22] - PS44*P[13][22] + PS60*P[2][22] + PS62*P[1][22] + PS72*P[0][22] - PS74*P[3][22] + P[4][22];
		nextP[5][22] = PS160*P[15][22] - PS162*P[13][22] - PS60*P[1][22] + PS62*P[2][22] - PS65*P[14][22] + PS72*P[3][22] + PS74*P[0][22] + P[5][22];
		nextP[6][22] = -PS165*P[14][22] + PS166*P[13][22] + PS60*P[0][22] + PS62*P[3][22] - PS70*P[15][22] - PS72*P[2][22] + PS74*P[1][22] + P[6][22];
		nextP[7][22] = P[4][22]*dt + P[7][22];
		nextP[8][22] = P[5][22]*dt + P[8][22];
		nextP[9][22] = P[6][22]*dt + P[9][22];
		nextP[10][22] = P[10][22];
		nextP[11][22] = P[11][22];
		nextP[12][22] = P[12][22];
		nextP[13][22] = P[13][22];
		nextP[14][22] = P[14][22];
		nextP[15][22] = P[15][22];
		nextP[16][22] = P[16][22];
		nextP[17][22] = P[17][22];
		nextP[18][22] = P[18][22];
		nextP[19][22] = P[19][22];
		nextP[20][22] = P[20][22];
		nextP[21][22] = P[21][22];
		nextP[22][22] = P[22][22];
		nextP[0][23] = -PS11*P[1][23] - PS12*P[2][23] - PS13*P[3][23] + PS6*P[10][23] + PS7*P[11][23] + PS9*P[12][23] + P[0][23];
		nextP[1][23] = PS11*P[0][23] - PS12*P[3][23] + PS13*P[2][23] - PS34*P[10][23] - PS7*P[12][23] + PS9*P[11][23] + P[1][23];
		nextP[2][23] = PS11*P[3][23] + PS12*P[0][23] - PS13*P[1][23] - PS34*P[11][23] + PS6*P[12][23] - PS9*P[10][23] + P[2][23];
		nextP[3][23] = -PS11*P[2][23] + PS12*P[1][23] + PS13*P[0][23] - PS34*P[12][23] - PS6*P[11][23] + PS7*P[10][23] + P[3][23];
		nextP[4][23] = -PS139*P[15][23] + PS140*P[14][23] - PS44*P[13][23] + PS60*P[2][23] + PS62*P[1][23] + PS72*P[0][23] - PS74*P[3][23] + P[4][23];
		nextP[5][23] = PS160*P[15][23] - PS162*P[13][23] - PS60*P[1][23] + PS62*P[2][23] - PS65*P[14][23] + PS72*P[3][23] + PS74*P[0][23] + P[5][23];
		nextP[6][23] = -PS165*P[14][23] + PS166*P[13][23] + PS60*P[0][23] + PS62*P[3][23] - PS70*P[15][23] - PS72*P[2][23] + PS74*P[1][23] + P[6][23];
		nextP[7][23] = P[4][23]*dt + P[7][23];
		nextP[8][23] = P[5][23]*dt + P[8][23];
		nextP[9][23] = P[6][23]*dt + P[9][23];
		nextP[10][23] = P[10][23];
		nextP[11][23] = P[11][23];
		nextP[12][23] = P[12][23];
		nextP[13][23] = P[13][23];
		nextP[14][23] = P[14][23];
		nextP[15][23] = P[15][23];
		nextP[16][23] = P[16][23];
		nextP[17][23] = P[17][23];
		nextP[18][23] = P[18][23];
		nextP[19][23] = P[19][23];
		nextP[20][23] = P[20][23];
		nextP[21][23] = P[21][23];
		nextP[22][23] = P[22][23];
		nextP[23][23] = P[23][23];

		// add process noise that is not from the IMU
		for (unsigned i = 22; i <= 23; i++) {
			nextP[i][i] += process_noise[i];
		}

	}

	// stop position covariance growth if our total position variance reaches 100m
	// this can happen if we lose gps for some time
	if ((P[7][7] + P[8][8]) > 1e4f) {
		for (uint8_t i = 7; i <= 8; i++) {
			for (uint8_t j = 0; j < _k_num_states; j++) {
				nextP[i][j] = P[i][j];
				nextP[j][i] = P[j][i];
			}
		}
	}

	// covariance matrix is symmetrical, so copy upper half to lower half
	for (unsigned row = 1; row < _k_num_states; row++) {
		for (unsigned column = 0 ; column < row; column++) {
			P[row][column] = P[column][row] = nextP[column][row];
		}
	}

	// copy variances (diagonals)
	for (unsigned i = 0; i < _k_num_states; i++) {
		P[i][i] = nextP[i][i];
	}

	// fix gross errors in the covariance matrix and ensure rows and
	// columns for un-used states are zero
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
	P_lim[0] = 1.0f;		// quaternion max var
	P_lim[1] = 1e6f;		// velocity max var
	P_lim[2] = 1e6f;		// positiion max var
	P_lim[3] = 1.0f;		// gyro bias max var
	P_lim[4] = 1.0f;		// delta velocity z bias max var
	P_lim[5] = 1.0f;		// earth mag field max var
	P_lim[6] = 1.0f;		// body mag field max var
	P_lim[7] = 1e6f;		// wind max var

	for (int i = 0; i <= 3; i++) {
		// quaternion states
		P[i][i] = math::constrain(P[i][i], 0.0f, P_lim[0]);
	}

	for (int i = 4; i <= 6; i++) {
		// NED velocity states
		P[i][i] = math::constrain(P[i][i], 1E-6f, P_lim[1]);
	}

	for (int i = 7; i <= 9; i++) {
		// NED position states
		P[i][i] = math::constrain(P[i][i], 1E-6f, P_lim[2]);
	}

	for (int i = 10; i <= 12; i++) {
		// gyro bias states
		P[i][i] = math::constrain(P[i][i], 0.0f, P_lim[3]);
	}

	// force symmetry on the quaternion, velocity and position state covariances
	makeSymmetrical(P, 0, 12);

	// the following states are optional and are deactivated when not required
	// by ensuring the corresponding covariance matrix values are kept at zero

	// accelerometer bias states
	if ((_params.fusion_mode & MASK_INHIBIT_ACC_BIAS) || _accel_bias_inhibit) {
		zeroRows(P, 13, 15);
		zeroCols(P, 13, 15);

	} else {
		// Find the maximum delta velocity bias state variance and request a covariance reset if any variance is below the safe minimum
		const float minSafeStateVar = 1e-9f;
		float maxStateVar = minSafeStateVar;
		bool resetRequired = false;

		for (uint8_t stateIndex = 13; stateIndex <= 15; stateIndex++) {
			if (P[stateIndex][stateIndex] > maxStateVar) {
				maxStateVar = P[stateIndex][stateIndex];

			} else if (P[stateIndex][stateIndex] < minSafeStateVar) {
				resetRequired = true;
			}
		}

		// To ensure stability of the covariance matrix operations, the ratio of a max and min variance must
		// not exceed 100 and the minimum variance must not fall below the target minimum
		// Also limit variance to a maximum equivalent to a 0.1g uncertainty
		const float minStateVarTarget = 5E-8f;
		float minAllowedStateVar = fmaxf(0.01f * maxStateVar, minStateVarTarget);

		for (uint8_t stateIndex = 13; stateIndex <= 15; stateIndex++) {
			P[stateIndex][stateIndex] = math::constrain(P[stateIndex][stateIndex], minAllowedStateVar, sq(0.1f * CONSTANTS_ONE_G * _dt_ekf_avg));
		}

		// If any one axis has fallen below the safe minimum, all delta velocity covariance terms must be reset to zero
		if (resetRequired) {
			float delVelBiasVar[3];

			// store all delta velocity bias variances
			for (uint8_t stateIndex = 13; stateIndex <= 15; stateIndex++) {
				delVelBiasVar[stateIndex - 13] = P[stateIndex][stateIndex];
			}

			// reset all delta velocity bias covariances
			zeroCols(P, 13, 15);

			// restore all delta velocity bias variances
			for (uint8_t stateIndex = 13; stateIndex <= 15; stateIndex++) {
				P[stateIndex][stateIndex] = delVelBiasVar[stateIndex - 13];
			}
		}

		// Run additional checks to see if the delta velocity bias has hit limits in a direction that is clearly wrong
		// calculate accel bias term aligned with the gravity vector
		float dVel_bias_lim = 0.9f * _params.acc_bias_lim * _dt_ekf_avg;
		float down_dvel_bias = 0.0f;

		for (uint8_t axis_index = 0; axis_index < 3; axis_index++) {
			down_dvel_bias += _state.accel_bias(axis_index) * _R_to_earth(2, axis_index);
		}

		// check that the vertical component of accel bias is consistent with both the vertical position and velocity innovation
		bool bad_acc_bias = (fabsf(down_dvel_bias) > dVel_bias_lim
				     && down_dvel_bias * _vel_pos_innov[2] < 0.0f
				     && down_dvel_bias * _vel_pos_innov[5] < 0.0f);

		// record the pass/fail
		if (!bad_acc_bias) {
			_fault_status.flags.bad_acc_bias = false;
			_time_acc_bias_check = _time_last_imu;

		} else {
			_fault_status.flags.bad_acc_bias = true;
		}

		// if we have failed for 7 seconds continuously, reset the accel bias covariances to fix bad conditioning of
		// the covariance matrix but preserve the variances (diagonals) to allow bias learning to continue
		if (_time_last_imu - _time_acc_bias_check > (uint64_t)7e6) {
			float varX = P[13][13];
			float varY = P[14][14];
			float varZ = P[15][15];
			zeroRows(P, 13, 15);
			zeroCols(P, 13, 15);
			P[13][13] = varX;
			P[14][14] = varY;
			P[15][15] = varZ;
			_time_acc_bias_check = _time_last_imu;
			_fault_status.flags.bad_acc_bias = false;
			ECL_WARN("EKF invalid accel bias - resetting covariance");

		} else {
			// ensure the covariance values are symmetrical
			makeSymmetrical(P, 13, 15);
		}

	}

	// magnetic field states
	if (!_control_status.flags.mag_3D) {
		zeroRows(P, 16, 21);
		zeroCols(P, 16, 21);

	} else {
		// constrain variances
		for (int i = 16; i <= 18; i++) {
			P[i][i] = math::constrain(P[i][i], 0.0f, P_lim[5]);
		}

		for (int i = 19; i <= 21; i++) {
			P[i][i] = math::constrain(P[i][i], 0.0f, P_lim[6]);
		}

		// force symmetry
		makeSymmetrical(P, 16, 21);
	}

	// wind velocity states
	if (!_control_status.flags.wind) {
		zeroRows(P, 22, 23);
		zeroCols(P, 22, 23);

	} else {
		// constrain variances
		for (int i = 22; i <= 23; i++) {
			P[i][i] = math::constrain(P[i][i], 0.0f, P_lim[7]);
		}

		// force symmetry
		makeSymmetrical(P, 22, 23);
	}
}

void Ekf::resetMagCovariance()
{
	// set the quaternion covariance terms to zero
	zeroRows(P, 0, 3);
	zeroCols(P, 0, 3);

	// set the magnetic field covariance terms to zero
	zeroRows(P, 16, 21);
	zeroCols(P, 16, 21);
	_mag_decl_cov_reset = false;

	// set the field state variance to the observation variance
	for (uint8_t rc_index = 16; rc_index <= 21; rc_index ++) {
		P[rc_index][rc_index] = sq(_params.mag_noise);
	}

	// save covariance data for re-use when auto-switching between heading and 3-axis fusion
	save_mag_cov_data();
}

void Ekf::resetWindCovariance()
{
	// set the wind  covariance terms to zero
	zeroRows(P, 22, 23);
	zeroCols(P, 22, 23);

	if (_tas_data_ready && (_imu_sample_delayed.time_us - _airspeed_sample_delayed.time_us < (uint64_t)5e5)) {
		// Use airspeed and zer sideslip assumption to set initial covariance values for wind states

		// calculate the wind speed and bearing
		float spd = sqrtf(sq(_state.wind_vel(0)) + sq(_state.wind_vel(1)));
		float yaw = atan2f(_state.wind_vel(1), _state.wind_vel(0));

		// calculate the uncertainty in wind speed and direction using the uncertainty in airspeed and sideslip angle
		// used to calculate the initial wind speed
		float R_spd = sq(math::constrain(_params.eas_noise, 0.5f, 5.0f) * math::constrain(_airspeed_sample_delayed.eas2tas, 0.9f, 10.0f));
		float R_yaw = sq(math::radians(10.0f));

		// calculate the variance and covariance terms for the wind states
		float cos_yaw = cosf(yaw);
		float sin_yaw = sinf(yaw);
		float cos_yaw_2 = sq(cos_yaw);
		float sin_yaw_2 = sq(sin_yaw);
		float sin_cos_yaw = sin_yaw * cos_yaw;
		float spd_2 = sq(spd);
		P[22][22] = R_yaw * spd_2 * sin_yaw_2 + R_spd * cos_yaw_2;
		P[22][23] = - R_yaw * sin_cos_yaw * spd_2 + R_spd * sin_cos_yaw;
		P[23][22] = P[22][23];
		P[23][23] = R_yaw * spd_2 * cos_yaw_2 + R_spd * sin_yaw_2;

		// Now add the variance due to uncertainty in vehicle velocity that was used to calculate the initial wind speed
		P[22][22] += P[4][4];
		P[23][23] += P[5][5];

	} else {
		// without airspeed, start with a small initial uncertainty to improve the initial estimate
		P[22][22] = sq(1.0f);
		P[23][23] = sq(1.0f);

	}
}
