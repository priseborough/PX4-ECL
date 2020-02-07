/****************************************************************************
 *
 *   Copyright (c) 2019 Paul Riseborough
 *
 ****************************************************************************/

/*
 * @file EKFGSF_yaw_estimator.cpp
 * Definition of functions used to provide a backup estimate of yaw angle
 * that doesn't use a magnetometer.
 * Uses a bank of 3-state EKF's organised as a Guassian sum filter
 * where states are velocity N,E (m/s) and yaw angle (rad)
 *
 * @author Paul Riseborough <p_riseborough@live.com.au>
 */

#include "ekf.h"
#include <ecl.h>
#include <mathlib/mathlib.h>
#include <cstdlib>

void Ekf::quatPredictEKFGSF(uint8_t model_index)
{
	// generate attitude solution using simple complementary filter for the selected model

	// Accelerometer correction
	// Project 'k' unit vector of earth frame to body frame
	// Vector3f k = quaterion.conjugate_inversed(Vector3f(0.0f, 0.0f, 1.0f));
	// Optimized version with dropped zeros
	Vector3f k(
		2.0f * (_ahrs_ekf_gsf[model_index].quat(1) * _ahrs_ekf_gsf[model_index].quat(3) - _ahrs_ekf_gsf[model_index].quat(0) * _ahrs_ekf_gsf[model_index].quat(2)),
		2.0f * (_ahrs_ekf_gsf[model_index].quat(2) * _ahrs_ekf_gsf[model_index].quat(3) + _ahrs_ekf_gsf[model_index].quat(0) * _ahrs_ekf_gsf[model_index].quat(1)),
		(_ahrs_ekf_gsf[model_index].quat(0) * _ahrs_ekf_gsf[model_index].quat(0) - _ahrs_ekf_gsf[model_index].quat(1) * _ahrs_ekf_gsf[model_index].quat(1) - _ahrs_ekf_gsf[model_index].quat(2) * _ahrs_ekf_gsf[model_index].quat(2) + _ahrs_ekf_gsf[model_index].quat(3) * _ahrs_ekf_gsf[model_index].quat(3))
	);

	// Perform angular rate correction using accel data and reduce correction as accel magnitude moves away from 1 g (reduces drift when vehicle picked up and moved).
	// During fixed wing flight, compensate for centripetal acceleration assuming coordinated turns and X axis forward
	Vector3f correction = {};
	Vector3f accel = _ahrs_accel;
	if (_ahrs_accel_norm > 0.5f * CONSTANTS_ONE_G && (_ahrs_accel_norm < 1.5f * CONSTANTS_ONE_G || _ahrs_turn_comp_enabled)) {
		if (_ahrs_turn_comp_enabled) {
			// turn rate is component of gyro rate about vertical (down) axis
			float turn_rate = _ahrs_ekf_gsf[model_index].R(2,0) * _ang_rate_delayed_raw(0)
					  + _ahrs_ekf_gsf[model_index].R(2,1) * _ang_rate_delayed_raw(1)
					  + _ahrs_ekf_gsf[model_index].R(2,2) * _ang_rate_delayed_raw(2);

			// use measured airspeed to calculate centripetal acceeration if available
			float centripetal_accel;
			if (_imu_sample_delayed.time_us - _airspeed_sample_delayed.time_us < 1000000) {
				centripetal_accel = _airspeed_sample_delayed.true_airspeed * turn_rate;
			} else {
				centripetal_accel = _params.EKFGSF_tas_default * turn_rate;
			}

			// project Y body axis onto horizontal and multiply by centripetal acceleration to give estimated
			// centrietal acceleratoin vector in earth frame due to coordinated turn
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
	float spinRate = _ang_rate_delayed_raw.length();
	if (spinRate < 0.175f) {
		_ahrs_ekf_gsf[model_index].gyro_bias -= correction * (_params.EKFGSF_gyro_bias_gain * _imu_sample_delayed.delta_ang_dt);

		for (int i = 0; i < 3; i++) {
			_ahrs_ekf_gsf[model_index].gyro_bias(i) = math::constrain(_ahrs_ekf_gsf[model_index].gyro_bias(i), -gyro_bias_limit, gyro_bias_limit);
		}
	}

	Vector3f rates = _ang_rate_delayed_raw - _ahrs_ekf_gsf[model_index].gyro_bias;

	// Feed forward gyro
	correction += rates;

	// Apply correction to state
	_ahrs_ekf_gsf[model_index].quat += _ahrs_ekf_gsf[model_index].quat.derivative1(correction) * _imu_sample_delayed.delta_ang_dt;

	// Normalize quaternion
	_ahrs_ekf_gsf[model_index].quat.normalize();

	// uodate body to earth frame rotation matrix
	_ahrs_ekf_gsf[model_index].R = quat_to_invrotmat(_ahrs_ekf_gsf[model_index].quat);

}

void Ekf::alignQuatEKFGSF()
{
	// The roll and pitch or tilt is obtained from the gravity vector and will be the same for all models
	// so only need to calculate this once

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
	for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++) {
		_ahrs_ekf_gsf[model_index].quat = R;
		_ahrs_ekf_gsf[model_index].quat.normalize();
		_ahrs_ekf_gsf[model_index].quat_initialised = true;
	}
}

void Ekf::alignQuatYawEKFGSF()
{
	// Align yaw angle for each model
	for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++) {
		if (fabsf(_ahrs_ekf_gsf[model_index].R(2, 0)) < fabsf(_ahrs_ekf_gsf[model_index].R(2, 1))) {
			// get the roll, pitch, yaw estimates from the rotation matrix using a  321 Tait-Bryan rotation sequence
			Eulerf euler_init(_ahrs_ekf_gsf[model_index].quat);

			// set the yaw angle
			euler_init(2) = wrap_pi(_ekf_gsf[model_index].X[2]);

			// update the quaternions and rotation matrix
			_ahrs_ekf_gsf[model_index].quat = Quatf(euler_init);
			_ahrs_ekf_gsf[model_index].R = quat_to_invrotmat(_ahrs_ekf_gsf[model_index].quat);

		} else {
			// Calculate the 312 Tait-Bryan sequence euler angles that rotate from earth to body frame
			// PX4 math library does not support this so are using equations from
			// http://www.atacolorado.com/eulersequences.doc
			Vector3f euler312;
			euler312(0) = atan2f(-_ahrs_ekf_gsf[model_index].R(0, 1), _ahrs_ekf_gsf[model_index].R(1, 1));  // first rotation (yaw)
			euler312(1) = asinf(_ahrs_ekf_gsf[model_index].R(2, 1)); // second rotation (roll)
			euler312(2) = atan2f(-_ahrs_ekf_gsf[model_index].R(2, 0), _ahrs_ekf_gsf[model_index].R(2, 2));  // third rotation (pitch)

			// set the yaw angle
			euler312(0) = wrap_pi(_ekf_gsf[model_index].X[2]);

			// Calculate the body to earth frame rotation matrix from the corrected euler angles
			float c2 = cosf(euler312(2));
			float s2 = sinf(euler312(2));
			float s1 = sinf(euler312(1));
			float c1 = cosf(euler312(1));
			float s0 = sinf(euler312(0));
			float c0 = cosf(euler312(0));

			_ahrs_ekf_gsf[model_index].R(0, 0) = c0 * c2 - s0 * s1 * s2;
			_ahrs_ekf_gsf[model_index].R(1, 1) = c0 * c1;
			_ahrs_ekf_gsf[model_index].R(2, 2) = c2 * c1;
			_ahrs_ekf_gsf[model_index].R(0, 1) = -c1 * s0;
			_ahrs_ekf_gsf[model_index].R(0, 2) = s2 * c0 + c2 * s1 * s0;
			_ahrs_ekf_gsf[model_index].R(1, 0) = c2 * s0 + s2 * s1 * c0;
			_ahrs_ekf_gsf[model_index].R(1, 2) = s0 * s2 - s1 * c0 * c2;
			_ahrs_ekf_gsf[model_index].R(2, 0) = -s2 * c1;
			_ahrs_ekf_gsf[model_index].R(2, 1) = s1;

			// update the quaternion
			_ahrs_ekf_gsf[model_index].quat = Quatf(_ahrs_ekf_gsf[model_index].R);
			_ahrs_ekf_gsf[model_index].quat.normalize();
		}
		_ahrs_ekf_gsf[model_index].quat_initialised = true;
	}
}

// predict states and covariance for specified model index
void Ekf::statePredictEKFGSF(uint8_t model_index)
{
	// generate an attitude reference using IMU data
	quatPredictEKFGSF(model_index);

	// we don't start running the EKF part of the algorithm until there are regular velocity observations
	if (!_ekf_gsf_vel_fuse_started) {
		return;
	}

	// Calculate the yaw state using a projection onto the horizontal that avoids gimbal lock
	if (fabsf(_ahrs_ekf_gsf[model_index].R(2, 0)) < fabsf(_ahrs_ekf_gsf[model_index].R(2, 1))) {
		// use 321 Tait-Bryan rotation to define yaw state
		_ekf_gsf[model_index].X[2] = atan2(_ahrs_ekf_gsf[model_index].R(1, 0), _ahrs_ekf_gsf[model_index].R(0, 0));
	} else {
		// use 312 Tait-Bryan rotation to define yaw state
		_ekf_gsf[model_index].X[2] = atan2(-_ahrs_ekf_gsf[model_index].R(0, 1), _ahrs_ekf_gsf[model_index].R(1, 1)); // first rotation (yaw)
	}

	// calculate delta velocity in a horizontal front-right frame
	Vector3f del_vel_NED = _ahrs_ekf_gsf[model_index].R * _imu_sample_delayed.delta_vel;
	float dvx =   del_vel_NED(0) * cosf(_ekf_gsf[model_index].X[2]) + del_vel_NED(1) * sinf(_ekf_gsf[model_index].X[2]);
	float dvy = - del_vel_NED(0) * sinf(_ekf_gsf[model_index].X[2]) + del_vel_NED(1) * cosf(_ekf_gsf[model_index].X[2]);

	// sum delta velocties in earth frame:
	_ekf_gsf[model_index].X[0] += del_vel_NED(0);
	_ekf_gsf[model_index].X[1] += del_vel_NED(1);

	// predict covariance - autocode from https://github.com/priseborough/3_state_filter/blob/flightLogReplay-wip/calcPupdate.txt

	// local copy required
	float P00 = _ekf_gsf[model_index].P[0][0];
	float P01 = _ekf_gsf[model_index].P[0][1];
	float P02 = _ekf_gsf[model_index].P[0][2];
	float P10 = _ekf_gsf[model_index].P[1][0];
	float P11 = _ekf_gsf[model_index].P[1][1];
	float P12 = _ekf_gsf[model_index].P[1][2];
	float P20 = _ekf_gsf[model_index].P[2][0];
	float P21 = _ekf_gsf[model_index].P[2][1];
	float P22 = _ekf_gsf[model_index].P[2][2];

	// Use fixed values for delta velocity and delta angle process noise variances
	const float dvxVar = sq(_params.EKFGSF_accel_noise * _imu_sample_delayed.delta_vel_dt); // variance of forward delta velocity - (m/s)^2
	const float dvyVar = dvxVar; // variance of right delta velocity - (m/s)^2
	const float dazVar = sq(_params.EKFGSF_gyro_noise * _imu_sample_delayed.delta_ang_dt); // variance of yaw delta angle - rad^2

	float t2 = sinf(_ekf_gsf[model_index].X[2]);
	float t3 = cosf(_ekf_gsf[model_index].X[2]);
	float t4 = dvy*t3;
	float t5 = dvx*t2;
	float t6 = t4+t5;
	float t8 = P22*t6;
	float t7 = P02-t8;
	float t9 = dvx*t3;
	float t11 = dvy*t2;
	float t10 = t9-t11;
	float t12 = dvxVar*t2*t3;
	float t13 = t2*t2;
	float t14 = t3*t3;
	float t15 = P22*t10;
	float t16 = P12+t15;

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
void Ekf::stateUpdateEKFGSF(uint8_t model_index)
{
	// fuse NE velocity observations
	float velObsVar = sq(fmaxf(_gps_sample_delayed.sacc, _params.gps_vel_noise));
	_ekf_gsf[model_index].innov[0] = _ekf_gsf[model_index].X[0] - _gps_sample_delayed.vel(0);
	_ekf_gsf[model_index].innov[1] = _ekf_gsf[model_index].X[1] - _gps_sample_delayed.vel(1);

	// copy covariance matrix to temporary variables
	float P00 = _ekf_gsf[model_index].P[0][0];
	float P01 = _ekf_gsf[model_index].P[0][1];
	float P02 = _ekf_gsf[model_index].P[0][2];
	float P10 = _ekf_gsf[model_index].P[1][0];
	float P11 = _ekf_gsf[model_index].P[1][1];
	float P12 = _ekf_gsf[model_index].P[1][2];
	float P20 = _ekf_gsf[model_index].P[2][0];
	float P21 = _ekf_gsf[model_index].P[2][1];
	float P22 = _ekf_gsf[model_index].P[2][2];

	// calculate Kalman gain K  and covariance matrix P
	// autocode from https://github.com/priseborough/3_state_filter/blob/flightLogReplay-wip/calcK.txt
	// and https://github.com/priseborough/3_state_filter/blob/flightLogReplay-wip/calcPmat.txt
	float t2 = P00*velObsVar;
 	float t3 = P11*velObsVar;
	float t4 = velObsVar*velObsVar;
	float t5 = P00*P11;
	float t9 = P01*P10;
	float t6 = t2+t3+t4+t5-t9;
	float t7;
	if (fabsf(t6) > 1e-6f) {
		t7 = 1.0f/t6;
	} else {
		t7 = 0.0f;
	}
	float t8 = P11+velObsVar;
	float t10 = P00+velObsVar;
	float K[3][2];

 	K[0][0] = -P01*P10*t7+P00*t7*t8;
	K[0][1] = -P00*P01*t7+P01*t7*t10;
	K[1][0] = -P10*P11*t7+P10*t7*t8;
	K[1][1] = -P01*P10*t7+P11*t7*t10;
	K[2][0] = -P10*P21*t7+P20*t7*t8;
	K[2][1] = -P01*P20*t7+P21*t7*t10;

	float t11 = P00*P01*t7;
	float t15 = P01*t7*t10;
	float t12 = t11-t15;
	float t13 = P01*P10*t7;
	float t16 = P00*t7*t8;
	float t14 = t13-t16;
	float t17 = t8*t12;
	float t18 = P01*t14;
	float t19 = t17+t18;
	float t20 = t10*t14;
	float t21 = P10*t12;
	float t22 = t20+t21;
	float t27 = P11*t7*t10;
	float t23 = t13-t27;
	float t24 = P10*P11*t7;
	float t26 = P10*t7*t8;
	float t25 = t24-t26;
	float t28 = t8*t23;
	float t29 = P01*t25;
	float t30 = t28+t29;
	float t31 = t10*t25;
	float t32 = P10*t23;
	float t33 = t31+t32;
	float t34 = P01*P20*t7;
	float t38 = P21*t7*t10;
	float t35 = t34-t38;
	float t36 = P10*P21*t7;
	float t39 = P20*t7*t8;
	float t37 = t36-t39;
	float t40 = t8*t35;
	float t41 = P01*t37;
	float t42 = t40+t41;
	float t43 = t10*t37;
	float t44 = P10*t35;
	float t45 = t43+t44;

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

	// calculatei innovation variance
	_ekf_gsf[model_index].S[0][0] = P00 + velObsVar;
	_ekf_gsf[model_index].S[1][1] = P11 + velObsVar;
	_ekf_gsf[model_index].S[0][1] = P01;
	_ekf_gsf[model_index].S[1][0] = P10;

	for (uint8_t obs_index = 0; obs_index < 2; obs_index++) {
		// apply the state corrections
		for (unsigned row = 0; row < 3; row++) {
			_ekf_gsf[model_index].X[row] -= K[row][obs_index] * _ekf_gsf[model_index].innov[obs_index];
		}
	}

	// Apply yaw correction to AHRS quaternion
	// TODO - This is an  expensive process due to the number of trig operations so a method of doing it more efficiently,
	// eg storing rotation matrix from the state prediction that doesn't include the yaw rotation should be investigated.
	// This matrix could then be multiplied with the yaw rotation to obtain the updated R matrix from which the updated
	// quaternion is calculated
	if (fabsf(_ahrs_ekf_gsf[model_index].R(2, 0)) < fabsf(_ahrs_ekf_gsf[model_index].R(2, 1))) {
		// using a 321 Tait-Bryan rotation to define yaw state
		// take roll pitch yaw from AHRS prediction
		Eulerf euler321(_ahrs_ekf_gsf[model_index].R);

		// replace the yaw angle using the EKF state estimate
		euler321(2) = _ekf_gsf[model_index].X[2];

		// update the quaternion used by the AHRS prediction algorithm
		_ahrs_ekf_gsf[model_index].quat = Quatf(euler321);

	} else {
		// Calculate the 312 sequence euler angles that rotate from earth to body frame
		// from the AHRS prediction and replace the yaw angle with the EKF state estimate
		// See http://www.atacolorado.com/eulersequences.doc
		Vector3f euler312;
		euler312(0) = _ekf_gsf[model_index].X[2]; // yaw rotation taken from 3-state EKF
		euler312(1) = asinf(_ahrs_ekf_gsf[model_index].R(2, 1)); // second rotation (roll) taken from AHRS
		euler312(2) = atan2f(-_ahrs_ekf_gsf[model_index].R(2, 0), _ahrs_ekf_gsf[model_index].R(2, 2));  // third rotation (pitch) taken from AHRS

		// Calculate the body to earth frame rotation matrix from the euler angles using a 312 rotation sequence
		float c2 = cosf(euler312(2));
		float s2 = sinf(euler312(2));
		float s1 = sinf(euler312(1));
		float c1 = cosf(euler312(1));
		float s0 = sinf(euler312(0));
		float c0 = cosf(euler312(0));

		Dcmf R_to_earth;
		R_to_earth(0, 0) = c0 * c2 - s0 * s1 * s2;
		R_to_earth(1, 1) = c0 * c1;
		R_to_earth(2, 2) = c2 * c1;
		R_to_earth(0, 1) = -c1 * s0;
		R_to_earth(0, 2) = s2 * c0 + c2 * s1 * s0;
		R_to_earth(1, 0) = c2 * s0 + s2 * s1 * c0;
		R_to_earth(1, 2) = s0 * s2 - s1 * c0 * c2;
		R_to_earth(2, 0) = -s2 * c1;
		R_to_earth(2, 1) = s1;

		// update the quaternion used by the AHRS prediction algorithm
		_ahrs_ekf_gsf[model_index].quat = Quatf(R_to_earth);

	}
}

void Ekf::initialiseEKFGSF()
{
	memset(&_ekf_gsf, 0, sizeof(_ekf_gsf));
	const float yaw_increment = 2.0f * M_PI_F / (float)N_MODELS_EKFGSF;
	for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++) {
		// evenly space initial yaw estimates in the region between +-Pi
		_ekf_gsf[model_index].X[2] = -M_PI_F + (0.5f * yaw_increment) + ((float)model_index * yaw_increment);

		// All filter modesl start with the same weight
		_ekf_gsf[model_index].W = 1.0f / (float)N_MODELS_EKFGSF;

		// Assume velocity within 0.5 m/s of zero at alignment
		if (_control_status.flags.gps) {
			_ekf_gsf[model_index].X[0] = _gps_sample_delayed.vel(0);
			_ekf_gsf[model_index].X[1] = _gps_sample_delayed.vel(1);
			_ekf_gsf[model_index].P[0][0] = sq(_gps_sample_delayed.sacc);
			_ekf_gsf[model_index].P[1][1] = _ekf_gsf[model_index].P[0][0];
		} else {
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
		// check for excessive acceleration.
		const float accel_norm_sq = _ahrs_accel.norm_squared();
		const float upper_accel_limit = CONSTANTS_ONE_G * 1.1f;
		const float lower_accel_limit = CONSTANTS_ONE_G * 0.9f;
		bool ok_to_align = ((accel_norm_sq > lower_accel_limit * lower_accel_limit &&
			  accel_norm_sq < upper_accel_limit * upper_accel_limit));
		if (ok_to_align) {
			initialiseEKFGSF();
			alignQuatEKFGSF();
			_ahrs_ekf_gsf_tilt_aligned = true;
		}
		return;
	}

	// calculate common values used by the AHRS prediction models
	_ahrs_accel_norm = _ahrs_accel.norm();
	_ahrs_turn_comp_enabled = _control_status.flags.fixed_wing && _params.EKFGSF_tas_default > FLT_EPSILON;
	if (_ahrs_accel_norm > CONSTANTS_ONE_G) {
		if (_ahrs_turn_comp_enabled && _ahrs_accel_norm <= 2.0f * CONSTANTS_ONE_G) {
			_ahrs_accel_fusion_gain = _params.EKFGSF_tilt_gain * sq(1.0f - (_ahrs_accel_norm - CONSTANTS_ONE_G)/CONSTANTS_ONE_G);
		} else if (_ahrs_accel_norm <= 1.5f * CONSTANTS_ONE_G) {
			_ahrs_accel_fusion_gain = _params.EKFGSF_tilt_gain * sq(1.0f - 2.0f * (_ahrs_accel_norm - CONSTANTS_ONE_G)/CONSTANTS_ONE_G);
		}
	} else if (_ahrs_accel_norm > 0.5f * CONSTANTS_ONE_G) {
		_ahrs_accel_fusion_gain = _params.EKFGSF_tilt_gain * sq(1.0f + 2.0f * (_ahrs_accel_norm - CONSTANTS_ONE_G)/CONSTANTS_ONE_G);
	}

	// AHRS prediction cycle for each model
	// This always runs
	for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
		statePredictEKFGSF(model_index);
	}

	// The 3-state EKF models only run when flying to avoid corrupted estimates due to operator handling and GPS interference
	if (_control_status.flags.gps && _gps_data_ready && _control_status.flags.in_air) {
		if (!_ekf_gsf_vel_fuse_started) {
			for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
				// use the first measurement to set the velocities and correspnding covariances
				_ekf_gsf[model_index].X[0] = _gps_sample_delayed.vel(0);
				_ekf_gsf[model_index].X[1] = _gps_sample_delayed.vel(1);
				_ekf_gsf[model_index].P[0][0] = sq(_gps_sample_delayed.sacc);
				_ekf_gsf[model_index].P[1][1] = _ekf_gsf[model_index].P[0][0];
				alignQuatYawEKFGSF();
			}
			_ekf_gsf_vel_fuse_started = true;
		} else {
			float total_w = 0.0f;
			float newWeight[N_MODELS_EKFGSF];
			for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
				// subsequent measuremnents are fused as direct state observations
				stateUpdateEKFGSF(model_index);

				// calculate weighting for each model assuming a normal distribution
				newWeight[model_index]= fmaxf(gaussianDensityEKFGSF(model_index) * _ekf_gsf[model_index].W, 0.0f);
				total_w += newWeight[model_index];
			}

			// normalise the weighting function
			if (_ekf_gsf_vel_fuse_started && total_w > 0.0f) {
				float total_w_inv = 1.0f / total_w;
				for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
					_ekf_gsf[model_index].W  = newWeight[model_index] * total_w_inv;
				}
			}

			// enforce a minimum weighting value
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
			float scale_factor = (unmodified_weights_sum - correction_sum - _params.EKFGSF_weight_min) / (unmodified_weights_sum - _params.EKFGSF_weight_min);
			for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
				if (!change_mask[model_index]) {
					_ekf_gsf[model_index].W = _params.EKFGSF_weight_min + scale_factor * (_ekf_gsf[model_index].W - _params.EKFGSF_weight_min);
				}
			}
		}
	} else if (_ekf_gsf_vel_fuse_started && !_control_status.flags.in_air) {
		// reset EKF states and wait to fly again
		initialiseEKFGSF();
		_ekf_gsf_vel_fuse_started = false;
	}

	// calculate a composite state vector as a weighted average of the states for each model
	memset(&X_GSF, 0, sizeof(X_GSF));
	for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index ++) {
		for (uint8_t state_index = 0; state_index < 3; state_index++) {
			X_GSF[state_index] += _ekf_gsf[model_index].X[state_index] * _ekf_gsf[model_index].W;
		}
	}
	/*
	// calculate a composite covariance matrix from a weighted average of the covariance for each model
	// models with larger innvoations are weighted less
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
}

float Ekf::gaussianDensityEKFGSF(uint8_t model_index)
{
	float t2 = _ekf_gsf[model_index].S[0][0] * _ekf_gsf[model_index].S[1][1];
	float t5 = _ekf_gsf[model_index].S[0][1] * _ekf_gsf[model_index].S[1][0];
	float t3 = t2 - t5; // determinant
	float t4; // determinant inverse
	if (fabsf(t3) > 1e-6f) {
		t4 = 1.0f/t3;
	} else {
		t4 = 1.0f/t2;
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

	normDist = expf(-0.5f * normDist);
	normDist *= sqrtf(t4)/ M_TWOPI_F;
	return normDist;
}

void Ekf::getDataEKFGSF(float *yaw_composite, float yaw[N_MODELS_EKFGSF], float innov_VN[N_MODELS_EKFGSF], float innov_VE[N_MODELS_EKFGSF], float weight[N_MODELS_EKFGSF])
{
	memcpy(yaw_composite, &X_GSF[2], sizeof(X_GSF[2]));
	for (uint8_t model_index = 0; model_index < N_MODELS_EKFGSF; model_index++) {
		yaw[model_index] = _ekf_gsf[model_index].X[2];
		innov_VN[model_index] = _ekf_gsf[model_index].innov[0];
		innov_VE[model_index] = _ekf_gsf[model_index].innov[1];
		weight[model_index] = _ekf_gsf[model_index].W;
	}
}

void Ekf::makeCovSymEKFGSF(uint8_t model_index)
{
	float P01 = 0.5f * (_ekf_gsf[model_index].P[0][1] + _ekf_gsf[model_index].P[1][0]);
	float P02 = 0.5f * (_ekf_gsf[model_index].P[0][2] + _ekf_gsf[model_index].P[2][0]);
	float P12 = 0.5f * (_ekf_gsf[model_index].P[1][2] + _ekf_gsf[model_index].P[2][1]);
	_ekf_gsf[model_index].P[0][1] = _ekf_gsf[model_index].P[1][0] = P01;
	_ekf_gsf[model_index].P[0][2] = _ekf_gsf[model_index].P[2][0] = P02;
	_ekf_gsf[model_index].P[1][2] = _ekf_gsf[model_index].P[2][1] = P12;
}

// Reset heading and magnetic field states
void Ekf::resetYawToEKFGSF()
{
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
		// use a 312 sequence

		// Calculate the 312 sequence euler angles that rotate from earth to body frame
		// See http://www.atacolorado.com/eulersequences.doc
		Vector3f euler312;
		euler312(0) = X_GSF[2];
		euler312(1) = asinf(_R_to_earth(2, 1)); // second rotation (roll)
		euler312(2) = atan2f(-_R_to_earth(2, 0), _R_to_earth(2, 2));  // third rotation (pitch)

		// Calculate the body to earth frame rotation matrix from the euler angles using a 312 rotation sequence
		float c2 = cosf(euler312(2));
		float s2 = sinf(euler312(2));
		float s1 = sinf(euler312(1));
		float c1 = cosf(euler312(1));
		float s0 = sinf(euler312(0));
		float c0 = cosf(euler312(0));

		Dcmf R_to_earth;
		R_to_earth(0, 0) = c0 * c2 - s0 * s1 * s2;
		R_to_earth(1, 1) = c0 * c1;
		R_to_earth(2, 2) = c2 * c1;
		R_to_earth(0, 1) = -c1 * s0;
		R_to_earth(0, 2) = s2 * c0 + c2 * s1 * s0;
		R_to_earth(1, 0) = c2 * s0 + s2 * s1 * c0;
		R_to_earth(1, 2) = s0 * s2 - s1 * c0 * c2;
		R_to_earth(2, 0) = -s2 * c1;
		R_to_earth(2, 1) = s1;

		// calculate initial quaternion states for the ekf
		// we don't change the output attitude to avoid jumps
		quat_after_reset = Quatf(R_to_earth);
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

	// update the yaw angle variance using half the nominal yaw separation ebtween models
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

	ECL_INFO_TIMESTAMPED("EKF emergency yaw reset");

	// stop using the magnetometer in the main EKF otherwise it's fusion could drag the yaw around
	// and cause another navigation failure
	_control_status.flags.mag_fault = true;
	ECL_WARN_TIMESTAMPED("EKF stopping magnetometer use");

}

// gets simple AHRS derived data which will be logged and used for algorithm development work
// returns false when no data available
bool Ekf::get_algo_test_data(float delAng[3],
		float *delAngDt,
		float delVel[3],
		float *delVelDt,
		float vel[3],
		float *velErr,
		bool *fuse_vel,
		float quat[4])
{
	bool vel_updated = _control_status.flags.gps && _gps_data_ready;

	delAng[0] = _imu_sample_delayed.delta_ang(0);
	delAng[1] = _imu_sample_delayed.delta_ang(1);
	delAng[2] = _imu_sample_delayed.delta_ang(2);

	memcpy(delAngDt, &_imu_sample_delayed.delta_vel_dt, sizeof(float));

	delVel[0] = _imu_sample_delayed.delta_vel(0);
	delVel[1] = _imu_sample_delayed.delta_vel(1);
	delVel[2] = _imu_sample_delayed.delta_vel(2);

	memcpy(delVelDt, &_imu_sample_delayed.delta_vel_dt, sizeof(float));

	vel[0] = _gps_sample_delayed.vel(0);
	vel[1] = _gps_sample_delayed.vel(1);
	vel[2] = _gps_sample_delayed.vel(2);

	memcpy(velErr, &_gps_sample_delayed.sacc, sizeof(float));

	memcpy(fuse_vel, &vel_updated, sizeof(bool));

	quat[0] =  _state.quat_nominal(0);
	quat[1] =  _state.quat_nominal(1);
	quat[2] =  _state.quat_nominal(2);
	quat[3] =  _state.quat_nominal(3);

	return _filter_initialised;
}

// request the EKF reset the yaw to the estimate from the internal EKF-GSF filter
// argment should be incremented only when a new reset is required
void Ekf::request_ekfgsf_yaw_reset(uint8_t counter)
{
	if (counter > _yaw_extreset_counter) {
		_yaw_extreset_counter = counter;
		_do_emergency_yaw_reset = true;
	}
}
