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

%% generate 'truth' trajectory data using 100Hz data
dt = 0.01;
time = [0:dt:30];
Ndata = length(time);

% IMU measurements
angRateZ_truth = zeros(1,Ndata);
accelX_truth = zeros(1,Ndata);
accelY_truth = zeros(1,Ndata);

% trajectory state history
velN_truth = zeros(1,Ndata);
velE_truth = zeros(1,Ndata);
yawAng_truth = zeros(1,Ndata);

% internal states
speed = 0.0;
yawAng = 0.0;
for index = 1:Ndata
    if (time(index) < 5)
        % stay stationary
        accelX_truth(index) = 0.0;
        angRateZ_truth(index) = 0.0;
    elseif (time(index) < 10)
        % accelerate at 2 m/s/s for 5 seconds in a straight line
        accelX_truth(index) = 5.0;
        angRateZ_truth(index) = 0.0;
    elseif (time(index) < 15)
        % turn right at 0.5 rad/sec
        accelX_truth(index) = 0.0;
        angRateZ_truth(index) = 0.5;
    elseif (time(index) < 20)
        % turn left at 0.5 rad/sec
        accelX_truth(index) = 0.0;
        angRateZ_truth(index) = -0.5;
    elseif (time(index) < 25)
        % decelerate at 5 m/s/s for 5 seconds in a straight line
        accelX_truth(index) = -5.0;
        angRateZ_truth(index) = 0.0;
    else
        % coast
        accelX_truth(index) = 0.0;
        angRateZ_truth(index) = 0.0;
    end
    speed = speed + dt*accelX_truth(index);
    yawAng = yawAng + dt*angRateZ_truth(index);
    yawAng_truth(index) = yawAng;
    velN_truth(index) = speed * cos(yawAng);
    velE_truth(index) = speed * sin(yawAng);
    accelY_truth(index) = speed * angRateZ_truth(index);
    
    angRateZ(index) = 0.003*randn + angRateZ_truth(index);
    accelX(index) = 0.35*randn + accelX_truth(index);
    accelY(index) = 0.35*randn + accelY_truth(index);
    
end

%% Implement UKF

% Step 1: Define UT Scaling parameters and weight vectors
L = 3; % Size of state vector [yawAng ; velN ; velE]
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;

% Step 2: define 1-Sigma noise assumptions
gyrNoise = 0.003;
accelNoise = 0.35;

% calculate state nise as a function of IMU noise
% this can also be calculated dynamically for additinal accuracy
Q = diag([gyrNoise*dt;accelNoise*dt;accelNoise*dt].^2);

% Assume 0.5 m/s velocity measurement noise
R = diag([0.5;0.5].^2);

% Step 3: init state and covariance
x = zeros(3,Ndata);
x(:,1) = [0.05;-0.1;0.1]; % set a small amount of initial state error
P0 = diag([0.05;0.0;0.0].^2);

% loop through filter
P = P0;
for k = 2:Ndata
    % generate the sigma points
    sP = chol(P,'lower');
    chi_p = [x(:,k-1), x(:,k-1)*ones(1,L)+sqrt(L+lambda)*sP, ...
        x(:,k-1)*ones(1,L)-sqrt(L+lambda)*sP];
    
    % Propagate each sigma point through the prediction
    for j = 1:7
        % propagate yaw angle using measured yaw rate
        chi_m(1,j) = chi_p(1,j) + 0.5*(angRateZ(k)+angRateZ(k-1))*dt;
        % propagate NE velocity using yaw angle and acceleration
        chi_m(2,j) = chi_p(2,j) + 0.5*dt*(accelX(k)*cos(chi_m(1,j)) +accelX(k-1)*cos(chi_p(1,j))) - 0.5*dt*(accelY(k)*sin(chi_m(1,j)) +accelY(k-1)*sin(chi_p(1,j)));
        chi_m(3,j) = chi_p(3,j) + 0.5*dt*(accelY(k)*cos(chi_m(1,j)) +accelY(k-1)*cos(chi_p(1,j))) + 0.5*dt*(accelX(k)*sin(chi_m(1,j)) +accelX(k-1)*sin(chi_p(1,j)));
    end
    
    % calculate mean of predicted state
    x_m = chi_m * wm;
    
    % calculate covariance of predicted state
    P_m = Q;
    for i=1:2*L+1
        P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
    end
    
    % perform Observation Transformation
    
    % propagate each sigma point through the transformation
    
    
end