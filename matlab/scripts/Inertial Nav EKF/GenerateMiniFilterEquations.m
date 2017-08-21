% IMPORTANT - This script requires the Matlab symbolic toolbox and takes ~3 hours to run

% Derivation of 2-D nav EKF using a local NE earth frame and 
% XY body fixed frame
% Sequential fusion of velocity measurements
% Fusion of true airspeed
% Sequential fusion of magnetic flux measurements
% 3 state architecture.
% IMU data is assumed to arrive at a constant rate with a time step of dt
% IMU delta angle and velocity data are used as control inputs,
% not observations

% Author:  Paul Riseborough

% State vector:
% yaw angle
% Velocity - m/sec (North, East)

% Observations:
% NE velocity - m/s

% Time varying parameters:
% Z delta angle measurements in body axes - rad
% XY delta velocity measurements in body axes - m/sec


%% define symbolic variables and constants
clear all;
reset(symengine);
syms daz 'real' % IMU delta angle measurements in body axes - rad
syms dvx dvy 'real' % IMU delta velocity measurements in body axes - m/sec
syms yaw 'real'
syms vn ve 'real' % NED velocity - m/sec
syms dt 'real' % IMU time step - sec
syms gravity 'real' % gravity  - m/sec^2
syms dazVar dvxVar dvyVar  'real'; % IMU delta angle and delta velocity measurement variances
syms R_VN R_VE 'real' % variances for NED velocity measurements - (m/sec)^2

%% define the state prediction equations

% derive the body to nav direction transformation matrix
Tbn = [cos(yaw) , -sin(yaw) ;...
       sin(yaw) ,  cos(yaw)];

% attitude update equation
yawNew = yaw + daz;

% velocity update equations
velNew = [vn;ve] + Tbn*[dvx;dvy];

% Define the state vector & number of states
stateVector = [yaw;vn;ve];
nStates=numel(stateVector);

% Define vector of process equations
newStateVector = [yawNew;velNew];

% derive the state transition matrix
F = jacobian(newStateVector, stateVector);
matlabFunction(F,'file','calcFmat.m');

%% derive the covariance prediction equations
% This reduces the number of floating point operations by a factor of 6 or
% more compared to using the standard matrix operations in code

% Error growth in the inertial solution is assumed to be driven by 'noise' in the delta angles and
% velocities, after bias effects have been removed. 

% derive the control(disturbance) influence matrix from IMU noise to state
% noise
G = jacobian(newStateVector, [daz;dvx;dvy]);

% derive the state error matrix
distMatrix = diag([dazVar dvxVar dvyVar]);
Q = G*distMatrix*transpose(G);
matlabFunction(Q,'file','calcQmat.m');

syms P11 P12 P13 P21 P22 P23 P31 P32 P33 'real'
P = [P11 P12 P13;...
     P21 P22 P23;...
     P31 P32 P33];

% Derive the predicted covariance matrix using the standard equation
P_update = F*P*transpose(F) + Q;
matlabFunction(P_update,'file','calcPupdate.m');