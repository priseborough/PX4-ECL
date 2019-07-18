% IMPORTANT - This script requires the Matlab symbolic toolbox
% Derivation of EKF equations for estimation of magnetomer offsets using
% knowledge of relative yaw angle, tilt and earth field strength, 
% declination and inclination
% Author:  Paul Riseborough
% Last Modified: 17 July 2019

% State vector:

% Body frame mag bias vector mx_b, my_b, mz_b (mGauss)
% Mag field scale factor m_scale - This is used to scale the earth field 
% Initial yaw angle yaw_init (rad) - This is the starting yaw angle

% Observations:

% Body frame mag field vector mx, my, mz (mGauss)

% Time varying parameters:

% Earth frame magnetic field vector mn, me, md (mGauss)
% Quaternions q0, q1, q2, q3


clear all;

%% define symbolic variables and constants
syms mx_bias my_bias mz_bias real % mag sensor bias : mGauss
syms mn me md real % earth field : mGausss
syms m_scale real % earth field scale factor 
syms yaw_init real % initial yaw angle : rad
syms mag_noise real % magnetometer process noise mGauss
syms q0 q1 q2 q3 real % quaternions representing rotation from earth to body frame
nStates = 5;

%% derive Jacobians for fusion of magnetomer  measurements
states=[mx_bias;my_bias;mz_bias;m_scale;yaw_init];

% rotate quaternions by initial yaw angle
quat = QuatMult([cos(yaw_init),0,0,sin(yaw_init)],[q0,q1,q2,q3]);

% derive the body to nav direction cosine matrix
Tbn = Quat2Tbn(quat);

% rotate earth field into body frame and add bias errors to give predicted
% observation
mag_obs = transpose(Tbn)*[m_scale*mn;m_scale*me;m_scale*md] + [mx_bias;my_bias;mz_bias];

% calculate the observation jacobians
H_MAG = jacobian(mag_obs,states);
ccode(H_MAG,'file','H_MAG.c');
fix_c_code('H_MAG.c');