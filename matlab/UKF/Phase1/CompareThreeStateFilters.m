% Implementation of 2-D nav UKF using a local NE earth frame and XY body fixed frame

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

% Time varying control parameters:
% Z delta angle measurements in body axes - rad
% XY delta velocity measurements in body axes - m/sec

%% derive EKF Jacobians
% run('GenerateMiniFilterEquations.m');

%% generate 'truth' trajectory data using 100Hz data
clear all;

dt = 0.01;
time = [0:dt:30];
Ndata = length(time);

% IMU truth measurements
angRateZ_truth = zeros(1,Ndata);
accelX_truth = zeros(1,Ndata);
accelY_truth = zeros(1,Ndata);

% IMU measurements with noise
angRateZ_meas = zeros(1,Ndata);
accelX_meas = zeros(1,Ndata);
accelY_meas = zeros(1,Ndata);

% trajectory state history
velN_truth = zeros(1,Ndata);
velE_truth = zeros(1,Ndata);
yawAng_truth = zeros(1,Ndata);

% IMU noise parameters
gyrNoise = 0.003; % rad/sec
accelNoise = 0.35; % m/sec^2

% internal states
speed = 0.0;
yawAng = 0.0;
oldSpeed = 0.0;
oldYawAng = 0.0;
lastGpsMeasTime = 0;
gpsIndex = 1;
for index = 2:Ndata
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
    
    % calculate truth trajectory using trapezoidal integration
    oldSpeed = speed;
    speed = speed + dt*0.5*(accelX_truth(index) + accelX_truth(index-1));
    oldYawAng = yawAng;
    yawAng = yawAng + dt*0.5*(angRateZ_truth(index) + angRateZ_truth(index-1));
    yawAng_truth(index) = yawAng;
    velN_truth(index) = 0.5*(speed * cos(yawAng) + oldSpeed * cos(oldYawAng));
    velE_truth(index) = 0.5*(speed * sin(yawAng) + oldSpeed * sin(oldYawAng));
    accelY_truth(index) = speed * angRateZ_truth(index);
    
    % generate IMU measurements with additive Gaussian noise;
    angRateZ_meas(index) = gyrNoise*randn + angRateZ_truth(index);
    accelX_meas(index) = accelNoise*randn + accelX_truth(index);
    accelY_meas(index) = accelNoise*randn + accelY_truth(index);
    
    %  collect GPS measurements at 10Hz
    if (time(index) - lastGpsMeasTime >= 0.1)
        lastGpsMeasTime = time(index);
        gpsTime(gpsIndex) = time(index);
        gpsVelN(gpsIndex) = velN_truth(index) + 0.5 * randn;
        gpsVelE(gpsIndex) = velE_truth(index) + 0.5 * randn;
        gpsIndex = gpsIndex + 1;
    end
    
end
Ngps = gpsIndex - 1;

%% Implement UKF

% Step 1: Define UT Scaling parameters and weight vectors
L = 6; % Size of augmented state vector [yawAng ; velN ; velE ; zGyrNoise ; xAccelNoise ; yAccelNoise]
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;

% Step 2: define 1-Sigma noise assumptions
% state process noise
gyrNoise = 0.003;
accelNoise = 0.35;

% Assume 0.5 m/s velocity measurement noise
R = diag([0.5;0.5].^2);

% define initial state estimation error
deg2rad = pi/180;
state_error_init = [45*deg2rad;-0.2;0.2];

% Initialise state and covariance
x_ukf = zeros(3,Ndata);
x_ukf(:,1) = [yawAng_truth(1);velN_truth(1);velE_truth(1)] + state_error_init;
P0 = diag(state_error_init.^2);
Q0 = diag([gyrNoise;accelNoise;accelNoise].^2);

% loop through filter
P = P0;
gpsIndex = 1;
nP = 3;
nQ = 3;
chi_m_a = zeros(3,(2*L+1));
sQ = chol(Q0,'lower');
for k = 2:Ndata
    % generate the sigma points for the augmented state vector
    sP = chol(P,'lower');
    sPA = [sP,zeros(nP,nQ);zeros(nQ,nP),sQ];
    % augmented state vector
    x_a_prev = [x_ukf(:,k-1);zeros(3,1)];
    chi_p_a = [x_a_prev, x_a_prev*ones(1,L)+sqrt(L+lambda)*sPA, ...
        x_a_prev*ones(1,L)-sqrt(L+lambda)*sPA];
    
    % Propagate each sigma point through the prediction for the augmented state vector
    for j = 1:(2*L+1)
        % propagate yaw angle using yaw rate plus noise
        yawRate = 0.5*(angRateZ_meas(k)+angRateZ_meas(k-1)) + chi_p_a(4,j);
        chi_m_a(1,j) = chi_p_a(1,j) + yawRate*dt;
        
        % calculate XY acceleration including noise at current time step
        accelNowX = accelX_meas(k) + chi_p_a(5,j);
        accelNowY = accelY_meas(k) + chi_p_a(6,j);
        
        % calculate XY acceleration including noise at previous time step
        accelPrevX = accelX_meas(k-1) + chi_p_a(5,j);
        accelPrevY = accelY_meas(k-1) + chi_p_a(6,j);
        
        % calculate NE acceleration including noise at current time step
        accelNowN = accelNowX*cos(chi_m_a(1,j)) - accelNowY*sin(chi_m_a(1,j));
        accelNowE = accelNowY*cos(chi_m_a(1,j)) + accelNowX*sin(chi_m_a(1,j));
        
        % calculate NE acceleration including noise at previous time step
        accelPrevN = accelPrevX*cos(chi_p_a(1,j)) - accelPrevY*sin(chi_p_a(1,j));
        accelPrevE = accelPrevY*cos(chi_p_a(1,j)) + accelPrevX*sin(chi_p_a(1,j));
        
        % propagate NE velocity using trapezoidanl integration of NE
        % acceleration
        chi_m_a(2,j) = chi_p_a(2,j) + 0.5*dt*(accelNowN + accelPrevN);
        chi_m_a(3,j) = chi_p_a(3,j) + 0.5*dt*(accelNowE + accelPrevE);
    end
    
    % calculate mean of predicted state
    x_m = chi_m_a(1:3,:) * wm;
    
    % calculate covariance of predicted state
    P = zeros(3,3);
    for i=1:2*L+1
        P = P + wc(i)*(chi_m_a(1:3,i) - x_m)*(chi_m_a(1:3,i) - x_m)';
    end
    
    % keep covariance matrix positive definite
    P = fix_covariance(P);
    
    % Check for GPS measurements
    if ((gpsIndex <= Ngps) && (gpsTime(gpsIndex) <= time(k)))
        % Propagate each sigma-point through observation transformation
        psi_m(1,:) = chi_m_a(2,:); % direct observation of North velocity
        psi_m(2,:) = chi_m_a(3,:); % direct observation of East velocity
        
        y_m = psi_m * wm; % Calculate mean of predicted output
        
        % Calculate covariance of predicted output
        % and cross-covariance between state and output
        Pyy = R;
        Pxy = zeros(3,2);
        for i = 1:2*L+1
            Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
            Pxy = Pxy + wc(i)*(chi_m_a(1:3,i) - x_m)*(psi_m(:,i) - y_m)';
        end
        
        % Measurement Update
        K = Pxy/Pyy; % Calculate Kalman gain
        x_m = x_m + K*([gpsVelN(gpsIndex);gpsVelE(gpsIndex)] - y_m); % Update state estimate
        P = P - K*Pyy*K'; % Update covariance estimate
        
        % keep covariance matrix positive definite
        P = fix_covariance(P);
        
        gpsIndex = gpsIndex + 1;
        
    end
    
    % log state estimate
    x_ukf(:,k) = x_m;
    
end

%% Implement EKF
% state process noise
gyrNoise = 0.003;
accelNoise = 0.35;

% Assume 0.5 m/s velocity measurement noise
R = diag([0.5;0.5].^2);

% Initialise state and covariance
x_ekf = zeros(3,Ndata);
x_ekf(:,1) = [yawAng_truth(1);velN_truth(1);velE_truth(1)] + state_error_init;

% loop through filter
gpsIndex = 1;
P = P0;
H = [0,1,0;0,0,1];
for k = 2:Ndata
    % predict state vector using trapezoidal integration of gyro and accel
    yawRate = 0.5*(angRateZ_meas(k)+angRateZ_meas(k-1));
    x_ekf(1,k) = x_ekf(1,k-1) + yawRate*dt;
    
    % calculate XY acceleration
    accelNowX = accelX_meas(k);
    accelNowY = accelY_meas(k);
    
    % calculate XY acceleration
    accelPrevX = accelX_meas(k-1);
    accelPrevY = accelY_meas(k-1);
    
    % calculate NE acceleration at current time step
    accelNowN = accelNowX*cos(x_ekf(1,k)) - accelNowY*sin(x_ekf(1,k));
    accelNowE = accelNowY*cos(x_ekf(1,k)) + accelNowX*sin(x_ekf(1,k));
    
    % calculate NE acceleration at previous time step
    accelPrevN = accelPrevX*cos(x_ekf(1,k-1)) - accelPrevY*sin(x_ekf(1,k-1));
    accelPrevE = accelPrevY*cos(x_ekf(1,k-1)) + accelPrevX*sin(x_ekf(1,k-1));
    
    % propagate NE velocity using trapezoidal integration of NE
    % acceleration
    x_ekf(2,k) =  x_ekf(2,k-1) + 0.5*dt*(accelNowN + accelPrevN);
    x_ekf(3,k) =  x_ekf(3,k-1) + 0.5*dt*(accelNowE + accelPrevE);
    
    % predict covariance
    Q = calcQmat((gyrNoise*dt)^2,(accelNoise*dt)^2,(accelNoise*dt)^2,x_ekf(1,k-1));
    F = calcFmat(0.5*dt*(accelNowX + accelPrevX),0.5*dt*(accelNowY + accelPrevY),x_ekf(1,k-1));
    P = F*P*transpose(F) + Q;
    
    % keep covariance matrix positive definite
    P = fix_covariance(P);
    
    % Check for GPS measurements and fuse if available
    if ((gpsIndex <= Ngps) && (gpsTime(gpsIndex) <= time(k)))
        % Use sequential fusion
        for axis=1:2
            % get innovation, observation variance and observation Jacobian
            if (axis == 1)
                innovation = x_ekf(2,k) - gpsVelN(gpsIndex);
                H = [0,1,0];
                obs_var = R(axis,axis);
            elseif (axis == 2)
                innovation = x_ekf(3,k) - gpsVelE(gpsIndex);
                H = [0,0,1];
                obs_var = R(axis,axis);
            end
            
            % calculate Kalman gain
            innovation_variance = H*P*transpose(H) + obs_var;
            K = (P * transpose(H)) / innovation_variance;
            
            % update states
            x_ekf(:,k) = x_ekf(:,k) - K * innovation;
            
            % update covariance
            P = P - K*innovation_variance*transpose(K);
            
            % keep covariance matrix positive definite
            P = fix_covariance(P);
            
        end
        
        gpsIndex = gpsIndex + 1;
        
    end
    
end

%% Plot data for UKF and EKF

close all;
figure;

subplot(3,1,1);
rad2deg = 180/pi;
plot(time,[yawAng_truth;x_ukf(1,:);x_ekf(1,:)]*rad2deg);
grid on;
title('UKF and EKF comparison');
xlabel('time (sec)');
ylabel('yaw angle (deg)');
legend('truth','UKF','EKF');

subplot(3,1,2);
plot(time,[velN_truth;x_ukf(2,:);x_ekf(2,:)]);
grid on;
xlabel('time (sec)');
ylabel('vel N (m/sec)');
legend('truth','UKF','EKF');

subplot(3,1,3);
plot(time,[velE_truth;x_ukf(3,:);x_ekf(3,:)]);
grid on;
xlabel('time (sec)');
ylabel('vel E (m/sec)');
legend('truth','UKF','EKF');
