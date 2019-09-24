/****************************************************************************
 *
 *   Copyright (c) 2017 Estimation and Control Library (ECL). All rights reserved.
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
 * @file magcal.h
 *
 * Small 4 state EKF that is able to estimate sensor biases from a single axis rotation
 * Uses sequential fusion of magnetometer measurements for estimation of bias offsets
 * Uses known earth field
 * Requires 360 deg yaw rotation performed when on-ground
 *
 * @author Paul Riseborough
 */

#pragma once

#include <mathlib/mathlib.h>
#include <matrix/math.hpp>

class MagCal
{
public:
	MagCal() = default;
	~MagCal() = default;

	// no copy, assignment, move, move assignment
	MagCal(const MagCal &) = delete;
	MagCal &operator=(const MagCal &) = delete;
	MagCal(MagCal &&) = delete;
	MagCal &operator=(MagCal &&) = delete;

	void runCalibrator();
	void setGpsLoc(const gps_message &gps_loc);
	void setMagData(float (&data)[3]);

private:

	void sampleDataMagCal();
	void processDataMagCal();
	void calcQuatMagCal();

	struct {
		Vector3f mag_bias;		// XYZ body franme magnetomer bias (Ga)
		float yaw_offset;		// yaw angle offset (rad)
	} _mag_cal_states{};			///< states used by mag bias EKF

	float _mag_cov_mat[5][5] = {};		///< covariance matrix used by mag bias EKF
	bool _mag_cal_sampling_active = false;	///< true when the mag bias EKF is active
	uint64_t _mag_cal_sample_time_us{0};	///< last time a mag sample was fused (uSec)
	float _mag_bias_ekf_yaw_last{0.0f};	///< yaw angle when data last used (rad)

	struct mag_cal_obs {
		Vector3f mag_data; 		// XYZ body frame mag field data (Ga)
		Quatf quaternion;  		// quaternion describing the rotation from body to earth frame
		float time_step;		// time lapsed from the previous measurement (sec)
	};					///< data structure used to store observation data points used by the calibration EKF

	mag_cal_obs _mag_cal_fit_data[36];	///< array of calibration data points stored for processing by the calibration EKF
	uint8_t _mag_sample_index{0};		///< index that increments with each new data point is stored for processing
	uint8_t _mag_cal_iteration_index{0};	///< loop counter used to control how many many timnes the calibration EKF loops through the data set
	float _mag_cal_residual[3];		///< mag calibration residuals resulting from previous pass of calibration EKF through the stored data (Ga)
	bool _mag_cal_sampling_complete{false};	///< true when sampling is complete processing of stored data can commence.
	bool _mag_cal_complete{false};		///< true when processing of the current stored data is complete
	int8_t _mag_cal_direction{0};		///< used to remember the direction of yaw rotation 1=CW, -1=CCW
	Vector3f _mag_field_EF_rotated;		///< magnetic field in earth frame used by the calibrator after application of yaw rotation to match vehicle assumed zero yaw at start (Ga)
	uint8_t _retry_count{0};		///< loop counter used to control how many times we rotate the earth field and re-start the EKF processing if the previous attempt failed to converge
	Quatf _mag_cal_quat;			///< quaternion describing rotation from body to earth frame and calculated using only IMU data.
	Vector3f _mag_cal_gyro_bias;		///< gyro bias learned and used by the quaternion calculation
	bool _mag_cal_quat_initialised{false};	///< true when calibrator quaternion has been aligned
	const float _mag_cal_nsamples{36};	///< number of samples used for calibration.
	bool _restart_sampling{false};		///< true when the sampling and processing should be restarted
	float _mag_declination_gps{0.0f};       ///< magnetic declination returned by the geo library using the last valid GPS position (rad)
	float _mag_inclination_gps{0.0f};	///< magnetic inclination returned by the geo library using the last valid GPS position (rad)
	float _mag_strength_gps{0.0f};	        ///< magnetic strength returned by the geo library using the last valid GPS position (T)
	Vector3f _mag_vec_NED;			///< magnetic earth field vector in NED axes obtained calculated from the strength, declination ans inclination data (T)
	gps_message &gps_loc{};			///< Last supplied GPS location
	bool _mag_earth_field_set{false};	///< true when the GPS location has been set
	bool _mag_data_received{false};		///< true when magnetomer data has been received and is awaiting processing
	Vector3f _mag_data;			///< Vector containing latest mag vector from sensor in Gauss

};
