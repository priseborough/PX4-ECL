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
 * Function for fusing gps and baro measurements/
 *
 * @author Paul Riseborough <p_riseborough@live.com.au>
 * @author Siddharth Bharat Purohit <siddharthbharatpurohit@gmail.com>
 *
 */

#include "ekf.h"
#include "mathlib.h"

void Ekf::fuseOptFlow()
{
	// TODO
}

void Ekf::get_flow_innov(float flow_innov[2])
{
	memcpy(flow_innov, _flow_innov, sizeof(_flow_innov));
}


void Ekf::get_flow_innov_var(float flow_innov_var[2])
{
	memcpy(flow_innov_var, _flow_innov_var, sizeof(_flow_innov_var));
}

void Ekf::get_drag_innov(float drag_innov[2])
{
	memcpy(drag_innov, _drag_innov, sizeof(_drag_innov));
}


void Ekf::get_drag_innov_var(float drag_innov_var[2])
{
	memcpy(drag_innov_var, _drag_innov_var, sizeof(_drag_innov_var));
}

// calculate optical flow gyro bias errors
void Ekf::calcOptFlowBias()
{
	// reset the accumulators if the time interval is too large
	if (_delta_time_of > 1.0f) {
		_imu_del_ang_of.setZero();
		_delta_time_of = 0.0f;
		return;
	}

	// if accumulation time differences are not excessive and accumulation time is adequate
	// compare the optical flow and and navigation rate data and calculate a bias error
	if ((fabsf(_delta_time_of - _flow_sample_delayed.dt) < 0.1f) && (_delta_time_of > 0.01f)) {
		// calculate a reference angular rate
		Vector3f reference_body_rate;
		reference_body_rate = _imu_del_ang_of * (1.0f / _delta_time_of);

		// calculate the optical flow sensor measured body rate
		Vector3f of_body_rate;
		of_body_rate = _flow_sample_delayed.gyroXYZ * (1.0f / _flow_sample_delayed.dt);

		// calculate the bias estimate using  a combined LPF and spike filter
		_flow_gyro_bias(0) = 0.99f * _flow_gyro_bias(0) + 0.01f * math::constrain((of_body_rate(0) - reference_body_rate(0)),
				     -0.1f, 0.1f);
		_flow_gyro_bias(1) = 0.99f * _flow_gyro_bias(1) + 0.01f * math::constrain((of_body_rate(1) - reference_body_rate(1)),
				     -0.1f, 0.1f);
		_flow_gyro_bias(2) = 0.99f * _flow_gyro_bias(2) + 0.01f * math::constrain((of_body_rate(2) - reference_body_rate(2)),
				     -0.1f, 0.1f);
	}

	// reset the accumulators
	_imu_del_ang_of.setZero();
	_delta_time_of = 0.0f;
}

// calculate the measurement variance for the optical flow sensor (rad/sec)^2
float Ekf::calcOptFlowMeasVar()
{
	// calculate the observation noise variance - scaling noise linearly across flow quality range
	float R_LOS_best = fmaxf(_params.flow_noise, 0.05f);
	float R_LOS_worst = fmaxf(_params.flow_noise_qual_min, 0.05f);

	// calculate a weighting that varies between 1 when flow quality is best and 0 when flow quality is worst
	float weighting = (255.0f - (float)_params.flow_qual_min);

	if (weighting >= 1.0f) {
		weighting = math::constrain(((float)_flow_sample_delayed.quality - (float)_params.flow_qual_min) / weighting, 0.0f,
					    1.0f);

	} else {
		weighting = 0.0f;
	}

	// take the weighted average of the observation noie for the best and wort flow quality
	float R_LOS = sq(R_LOS_best * weighting + R_LOS_worst * (1.0f - weighting));

	return R_LOS;
}
