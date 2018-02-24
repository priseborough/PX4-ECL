% Implementation of 6DoF navigation UKF using a local NED navigation frame and XYZ body fixed frame

% IMU data is assumed to arrive at a constant rate with a time step of dt
% IMU data is used as control inputs, not observations

% Author:  Paul Riseborough

% It uses a augmented state vector 12x1 = 6x1 vehicle states + 6x1 control input noise states:

% Vehicle States
% 3x1 Generalized Rodrigues Parameter (GRP) vector representing the attitude
% error with parameters chosen such that the vector is equivalent to the 
% delta angle for small angles - rad
% 3x1 Velocity - m/sec (North, East, Down)

% Control Input Noise States
% 3x1 gyro noise - rad/sec
% 3x1 accelerometer noise - m/s^2

% Observations 3x1:
% GPS NED velocity - m/s

% Time varying control parameters 6x1:
% XYZ angular rate measurements in body frame - rad/sec
% XYZ specific force measurements in body frame - m/sec^2

%% generate 'truth' trajectory data using 100Hz data for 2D motion on a horizontal plane
clear all;

deg2rad = pi/180;
gravity = 9.80665;
dt = 0.01;
time = [0:dt:30];
Ndata = length(time);

% IMU truth measurements
angRate_truth = zeros(3,Ndata);
accel_truth = zeros(3,Ndata);

% IMU measurements with noise
angRate_meas = zeros(3,Ndata);
accel_meas = zeros(3,Ndata);

% IMU noise parameters
gyrNoise = 0.003; % rad/sec
accelNoise = 0.35; % m/sec^2

% GPS velocity noise
gpsVelNoise = 0.5; % m/sec

% define initial speed and yaw angle
oldSpeed = 0.0;
oldYawAng = -45*deg2rad;
speed = oldSpeed;
yawAng = oldYawAng;

% trajectory state history
velNED_truth = zeros(3,Ndata);
yawAng_truth = zeros(1,Ndata);

yawAng_truth(1) = yawAng;
velN_truth(1) = speed * cos(yawAng);
velE_truth(1) = speed * sin(yawAng);

% use trapezoidal integration of acceleration and angular rate to calculate truth trajectory 
lastGpsMeasTime = 0;
gpsIndex = 1;
for index = 2:Ndata
    if (time(index) < 5)
        % stay stationary
        accel_truth(1,index) = 0.0;
        angRate_truth(3,index) = 0.0;
    elseif (time(index) < 10)
        % accelerate at 2 m/s/s for 5 seconds in a straight line
        accel_truth(1,index) = 5.0;
        angRate_truth(3,index) = 0.0;
    elseif (time(index) < 15)
        % turn right at 0.5 rad/sec
        accel_truth(1,index) = 0.0;
        angRate_truth(3,index) = 0.5;
    elseif (time(index) < 20)
        % turn left at 0.5 rad/sec
        accel_truth(1,index) = 0.0;
        angRate_truth(3,index) = -0.5;
    elseif (time(index) < 25)
        % decelerate at 5 m/s/s for 5 seconds in a straight line
        accel_truth(1,index) = -5.0;
        angRate_truth(3,index) = 0.0;
    else
        % coast
        accel_truth(1,index) = 0.0;
        angRate_truth(3,index) = 0.0;
    end
    
    % calculate truth trajectory using trapezoidal integration
    oldSpeed = speed;
    speed = speed + dt*0.5*(accel_truth(1,index) + accel_truth(1,index-1));
    oldYawAng = yawAng;
    yawAng = yawAng + dt*0.5*(angRate_truth(3,index) + angRate_truth(3,index-1));
    yawAng_truth(index) = yawAng;
    velNED_truth(1,index) = 0.5*(speed * cos(yawAng) + oldSpeed * cos(oldYawAng));
    velNED_truth(2,index) = 0.5*(speed * sin(yawAng) + oldSpeed * sin(oldYawAng));
    velNED_truth(3,index) = 0.0;
    accel_truth(2,index) = speed * angRate_truth(3,index);
    
    % assume constant 1g gravity and vehicle stays level
    accel_truth(3,index) = -gravity;
    
    % generate IMU measurements with additive Gaussian noise;
    angRate_meas(:,index) = gyrNoise*randn(3,1) + angRate_truth(:,index);
    accel_meas(:,index) = accelNoise*randn(3,1) + accel_truth(:,index);
    
    %  collect GPS measurements at 10Hz
    if (time(index) - lastGpsMeasTime >= 0.1)
        lastGpsMeasTime = time(index);
        gpsTime(gpsIndex) = time(index);
        gpsVelNED(:,gpsIndex) = velNED_truth(:,index) + gpsVelNoise * randn(3,1);
        gpsIndex = gpsIndex + 1;
    end
    
end
Ngps = gpsIndex - 1;

%% Implement UKF

% Define Unscented Transform scaling parameters and weight vectors
L = 12; % Size of augmented state vector [3x1 rotVec ; 3x1 velNED ; 3x1 gyrNoise ; 3x1 accelNoise]
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;

% Quaternons are replaced by three-dimensional generalized Rodrigues 
% parameter (GRP) vector. This vector represents the attitude error and is
% initialised to zero. Parameters are chosen so that the vector magnitude
% is equal to the magnitude of the error for small rotations
a = 1;
f = 2*(a+1);

% Assume constant GPS velocity measurement noise
R = diag([gpsVelNoise;gpsVelNoise;gpsVelNoise].^2);

% define initial state estimation error
tilt_error_std = 5*deg2rad;
yaw_error_std = 30*deg2rad;
vel_error_std = 0.2; 
state_error_init = [tilt_error_std;tilt_error_std;yaw_error_std;vel_error_std*[1;1;1]];

% Initial covariance matrices for vehicle states and control input noise
% states
P0 = diag(state_error_init.^2);
Q0 = diag([gyrNoise * [1;1;1] ; accelNoise * [1;1;1]].^2);

% initialise quaternions and velocity states
quat_state = [1;0;0;0];
vel_state = [0;0;0];

% initialise arrays used to store output
vel_ukf = zeros(3,Ndata);
quat_ukf = zeros(4,Ndata);
euler_ukf = zeros(3,Ndata);
vel_ukf(:,1) = vel_state;
quat_ukf(:,1) = quat_state;
euler_ukf(:,1) = QuatToEul(quat_state);

P = P0; % initialise the covariance matrix for the vehicle states

gpsIndex = 1; % counter used to access GPS data
nP = 6; % number of vehicle states
nQ = 6; % number of control input noise states
sigma_x_a = zeros(L,(2*L+1)); % array containing the sigma points for the augmented state vector

% Precompute the lower diagonal Cholesky decomposition for the input noise
% covariance matrix because it is time invariant
sQ = chol(Q0,'lower');

% Loop through IMU data
for k = 2:Ndata
    % Calculate the lower diagonal Cholesky decomposition for the vehicle
    % state covariance matrix. This requires nP^3 operations
    sP = chol(P,'lower');
    
    % Assemble the augmented covariance matrix
    sPA = [sP,zeros(nP,nQ);zeros(nQ,nP),sQ];
    
    % Expected value of augmented state vector from previous frame
    x_a_prev = [...
        zeros(3,1);... % attitude error vector starts from zero by definition
        vel_state;... % NED velocity
        zeros(nQ,1)... % IMU noise is zero mean
        ];
    
    % Generate sigma points for the augmented state vector
    % Note : For a real time implementation we should try to take advantage
    % of sparseness in sPA and the fact that sqrt(L+lambda) is constant
    sigma_x_a = [x_a_prev, x_a_prev*ones(1,L)+sqrt(L+lambda)*sPA, x_a_prev*ones(1,L)-sqrt(L+lambda)*sPA];
    % convert the attitude error vector sigma points to equivalent delta 
    % quaternions 
    dq(:,1) = [1;0;0;0];
    for s=2:(2*L+1)
        normsigmaX2 = sigma_x_a(1:3,s)'*sigma_x_a(1:3,s);
        dq(1,s) = (-a*normsigmaX2 + f*sqrt( f^2 + (1-a^2)*normsigmaX2 ) ) / ( f^2 + normsigmaX2 );
        dq(2:4,s) = (a + dq(1,s)) * sigma_x_a(1:3,s) / f;
    end
    
    % Apply the delta quaternions to the previous estimate to calculate the
    % quaternion sigma points. Although we could propagate the delta angles
    % through the vehicle state prediction, it is more accurate to use
    % quaternions and convert back to a set of GRP attitude error vectors
    % when the covariance information needs to be extracted.
    sigmaQuat(:,1) = quat_state;
    for s=2:(2*L+1)
        sigmaQuat(:,s) = QuatMult(sigmaQuat(:,1),dq(:,s));
    end
    
    % Propagate each sigma point through the vehicle state prediction
    % equations using trapezoidal integration. 
    for s = 1:(2*L+1)
        
        % Store the previous rotation from body to nav frame
        Tbn_prev = Quat2Tbn(sigmaQuat(:,s));
        
        % Calculate the delta angle including gyro process noise
        % use trapezoidal integration for IMU rates and forward euler
        % integration for IMU rate noise
        delAng = 0.5*(angRate_meas(:,k)+angRate_meas(:,k-1))*dt + sigma_x_a(7:9,s)*dt;
        
        % Convert to a quaternion
        delQuat = RotToQuat(delAng);
        
        % Rotate the sigma point quaternion forward from k-1 to k
        sigmaQuat(:,s) = QuatMult(sigmaQuat(:,s),delQuat);
        
        % Convert to a rotation matrix from body to nav frame
        Tbn_now = Quat2Tbn(sigmaQuat(:,s));
        
        % Calculate acceleration including noise at current and previous time step
        accelNow = accel_meas(:,k) + sigma_x_a(10:12,s);
        accelPrev = accel_meas(:,k-1) + sigma_x_a(10:12,s);
        
        % Rotate into earth frame and correct for gravity
        accelNowNED = Tbn_now * accelNow + [0;0;gravity];
        accelPrevNED = Tbn_prev * accelPrev + [0;0;gravity];
        
        % Propagate NED velocity states using trapezoidal integration of
        % acceleration
        sigma_x_a(4:6,s) = sigma_x_a(4:6,s) + 0.5*dt*(accelNowNED + accelPrevNED);
        
    end
    
    % Store the quaternion state estimate
    % The attitude error is zero mean by definition so we do not need to
    % calculate the mean from the sigma points
    quat_state = sigmaQuat(:,1);
    
    % Convert propagated quaternions to delta quaternions around the
    % expected value
    dq(:,1) = [1;0;0;0];
    for s= 2:(2*L+1)
        dq(:,s) = QuatMult(sigmaQuat(:,1).*[1;-1;-1;-1] , sigmaQuat(:,s));
    end
    
    % Convert error quaternions to attitude error vector
    sigma_x_a(1:3,1) = [0;0;0];
    for s=2:(2*L+1)
        sigma_x_a(1:3,s) = f * dq(2:4,s)/( a + dq(1,s) );
    end
    
    % Calculate mean of predicted states
    x_m = sigma_x_a(1:nP,:) * wm;
    
    % Store the velocity state estimate
    vel_state = x_m(4:6);
    
    % Calculate covariance of predicted states
    P = zeros(nP,nP);
    for i=1:2*L+1
        P = P + wc(i)*(sigma_x_a(1:nP,i) - x_m)*(sigma_x_a(1:nP,i) - x_m)';
    end
    
    % Keep the covariance matrix positive definite
    P = fix_covariance(P);
    
    % Check for GPS measurements
    if ((gpsIndex <= Ngps) && (gpsTime(gpsIndex) <= time(k)))
        % Propagate each sigma-point through observation transformation
        psi_m(1:3,:) = sigma_x_a(4:6,:); % direct observation of velocity
        
        % y_m = psi_m * wm; % Calculate mean of predicted output
        y_m = x_m(4:6); % same as psi_m * wm because we are directly observing states
        
        % Calculate covariance of predicted output
        % and cross-covariance between state and output
        Pyy = R;
        Pxy = zeros(nP,3);
        for i = 1:2*L+1
            Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
            Pxy = Pxy + wc(i)*(sigma_x_a(1:nP,i) - x_m)*(psi_m(:,i) - y_m)';
        end
        
        % Calculate Kalman gain
        K = Pxy/Pyy; % Calculate Kalman gain
        
        % Update the state estimate
        x_m = x_m + K*(gpsVelNED(:,gpsIndex) - y_m);
        
        % Update the covariance estimate
        P = P - K*Pyy*K'; % Update covariance estimate
        
        % Keep the covariance matrix positive definite
        P = fix_covariance(P);
        
        % Calculate the quaternion correction from the attitude error
        % vector
        normdp2 = x_m(1:3)'*x_m(1:3);
        dq(1,1) = ( -a * normdp2 + f * sqrt( f^2 + (1-a^2) * normdp2 ) ) / ( f^2 + normdp2 );
        dq(2:4,1) = x_m(1:3) * ( a + dq(1,1) ) / f;
        
        % apply the correction to the stored quaternion state and zero the
        % attitude error vector state
        quat_state = QuatMult(quat_state , dq(:,1));
        x_m(1:3) = [0;0;0];
        
        % store the updated velocity state estimate
        vel_state = x_m(4:6);
        
        gpsIndex = gpsIndex + 1;
        
    end
    
    % Note that if we needed to fuse in another observation on the same
    % prediction cycle, then we first need to regenerate the sigma points
    % which requires recalculating the Cholesky deconposition for the 
    % vehicle state covariance matrix P, assembling the augmented
    % covariance matrix, assembling the augmented state vector and
    % regenerating the points using the unscented transform.
    
    % log outputs
    euler_ukf(:,k) = QuatToEul(quat_state);
    quat_ukf(:,k) = quat_state;
    vel_ukf(:,k) = vel_state;
    
    
end

%% Plot data for UKF

close all;
figure;

subplot(3,1,1);
rad2deg = 180/pi;
plot(time,[yawAng_truth;euler_ukf(3,:)]*rad2deg);
grid on;
title('UKF test');
xlabel('time (sec)');
ylabel('yaw angle (deg)');
legend('truth','UKF');

subplot(3,1,2);
plot(time,[velNED_truth(1,:);vel_ukf(1,:)]);
grid on;
xlabel('time (sec)');
ylabel('vel N (m/sec)');
legend('truth','UKF');

subplot(3,1,3);
plot(time,[velNED_truth(2,:);vel_ukf(2,:)]);
grid on;
xlabel('time (sec)');
ylabel('vel E (m/sec)');
legend('truth','UKF');
