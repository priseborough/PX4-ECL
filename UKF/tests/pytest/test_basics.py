###############################################################################
#
# Copyright (c) 2017 Estimation and Control Library (ECL). All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
# 3. Neither the name ECL nor the names of its contributors may be
#    used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
# AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
###############################################################################

"""Basic tests and benchmarks testing the ecl's UKF

@author: Peter Dürr <Peter.Duerr@sony.com>
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import

import pytest

from test_utils import ecl

# Pylint does not like pytest fixtures, disable the warning
# pylint: disable=redefined-outer-name
from test_utils import update_sensors
from test_utils import initialized_ukf  # pylint: disable=unused-import


def test_filter_initialized(initialized_ukf):
    """Make sure the UKF updates after a few IMU, Mag and Baro updates
    """
    ukf, _ = initialized_ukf
    assert ukf.update()


@pytest.mark.benchmark
def test_ukf_imu_update_benchmark(initialized_ukf, benchmark):
    """Benchmark the time for an ukf update using just IMU, mag and baro data
    """
    ukf, _ = initialized_ukf
    benchmark(ukf.update)


def test_status_on_imu_mag_baro(initialized_ukf):
    """Make sure the UKF with zero inputs has reasonable status
    """
    ukf, _ = initialized_ukf
    control_mode = ukf.get_control_mode()

    assert control_mode.yaw_align
    assert control_mode.gps is False
    assert control_mode.opt_flow is False
    assert control_mode.mag_hdg
    assert control_mode.mag_3D is False
    assert control_mode.mag_dec is False
    assert control_mode.in_air is False
    assert control_mode.wind is False
    assert control_mode.baro_hgt
    assert control_mode.rng_hgt is False
    assert control_mode.gps_hgt is False
    assert control_mode.ev_pos is False
    assert control_mode.ev_yaw is False
    assert control_mode.ev_hgt is False
    assert control_mode.fuse_beta is False
    assert control_mode.update_mag_states_only is False
    assert control_mode.fixed_wing is False


def test_converges_to_zero():
    """Make sure the UKF with zero inputs converges to / stays at zero
    """
    ukf = ecl.Ukf()

    time_usec = 1000
    dt_usec = 5000

    # Provide a few samples
    for _ in range(10000):
	update_sensors(ukf, time_usec, dt_usec)
	ukf.update()
        time_usec += dt_usec

    converged_pos = ukf.get_position()
    converged_vel = ukf.get_velocity()
    converged_accel_bias = ukf.get_accel_bias()
    converged_gyro_bias = ukf.get_gyro_bias()
    for i in range(3):
        assert converged_pos[i] == pytest.approx(0.0, abs=1e-6)
        assert converged_vel[i] == pytest.approx(0.0, abs=1e-6)
        assert converged_accel_bias[i] == pytest.approx(0.0, abs=1e-5)
        assert converged_gyro_bias[i] == pytest.approx(0.0, abs=1e-5)
