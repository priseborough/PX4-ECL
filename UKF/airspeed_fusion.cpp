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
 * @file airspeed_fusion.cpp
 * airspeed fusion methods.
 *
 * @author Carl Olsson <carlolsson.co@gmail.com>
 * @author Roman Bast <bapstroman@gmail.com>
 * @author Paul Riseborough <p_riseborough@live.com.au>
 *
 */
#include "../ecl.h"
#include "ukf.h"
#include "mathlib.h"

void Ukf::fuseAirspeed()
{
	// TODO
}

void Ukf::get_wind_velocity(float *wind)
{
	wind[0] = _ukf_states.data.wind_vel(0);
	wind[1] = _ukf_states.data.wind_vel(1);
}

void Ukf::get_wind_velocity_var(float *wind_var)
{
	wind_var[0] = P_UKF(21,21);
	wind_var[1] = P_UKF(22,22);
}

void Ukf::get_true_airspeed(float *tas)
{
	float tempvar = sqrtf(sq(_state.vel(0) - _state.wind_vel(0)) + sq(_state.vel(1) - _state.wind_vel(1)) + sq(_state.vel(2)));
	memcpy(tas, &tempvar, sizeof(float));
}

/*
 * Reset the wind states using the current airspeed measurement, ground relative nav velocity, yaw angle and assumption of zero sideslip
*/
void Ukf::resetWindStates()
{
	// get euler yaw angle
	Eulerf euler321(_ukf_states.data.quat);
	float euler_yaw = euler321(2);

	if (_tas_data_ready && (_imu_sample_delayed.time_us - _airspeed_sample_delayed.time_us < (uint64_t)5e5)) {
		// estimate wind using zero sideslip assumption and airspeed measurement if airspeed available
		_ukf_states.data.wind_vel(0) = _ukf_states.data.vel(0) - _airspeed_sample_delayed.true_airspeed * cosf(euler_yaw);
		_ukf_states.data.wind_vel(1) = _ukf_states.data.vel(1) - _airspeed_sample_delayed.true_airspeed * sinf(euler_yaw);

	} else {
		// If we don't have an airspeed measurement, then assume the wind is zero
		_ukf_states.data.wind_vel(0) = 0.0f;
		_ukf_states.data.wind_vel(1) = 0.0f;
	}
}
