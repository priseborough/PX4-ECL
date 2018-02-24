function output = RunFilter(param,imu_data,mag_data,baro_data,gps_data,varargin)

% compulsory inputs

% param : parameters defined  by SetParameterDefaults.m
% imu_data : IMU delta angle and velocity data in body frame
% mag_data : corrected magnetometer field measurements in body frame
% baro_data : barometric height measurements
% gps_data : GPS NED pos vel measurements in local earth frame

% optional inputs

% rng _data : measurements for a Z body axis aligned range finder
% flow_data : XY axis optical flow angular rate measurements in body frame
% viso_data : ZED camera visula odometry measurements

nVarargs = length(varargin);
if nVarargs >= 2
    flowDataPresent = ~isempty(varargin{1}) && ~isempty(varargin{2});
    rng_data = varargin{1};
    flow_data = varargin{2};
    if flowDataPresent
        fprintf('Using optical Flow Data\n',nVarargs);
    end
else
    flowDataPresent = 0;
end

if nVarargs >= 3
    visoDataPresent = ~isempty(varargin{3});
    viso_data = varargin{3};
    if visoDataPresent
        fprintf('Using ZED camera odometry data\n',nVarargs);
    end
else
    visoDataPresent = 0;
end


%% Set initial conditions

% constants
deg2rad = pi/180;
gravity = 9.80665; % initial value of gravity - will be updated when WGS-84 position is known

% initialise the state vector and quaternion
[q_m, x_m, imu_start_index] = InitStates(param,imu_data,gps_data,mag_data,baro_data,param.control.minStartTime);

dt_imu_avg = 0.5 * (median(imu_data.gyro_dt) + median(imu_data.accel_dt));
indexStop = length(imu_data.time_us) - imu_start_index;

% indexStop = 15000; % Only use 60 seconds of data for initial testing

% create static structs for output data
output = struct('time_lapsed',[]',...
    'euler_angles',[],...
    'velocity_NED',[],...
    'position_NED',[],...
    'gyro_bias',[],...
    'accel_bias',[],...
    'mag_NED',[],...
    'mag_XYZ',[],...
    'wind_NE',[],...
    'dt',0,...
    'state_variance',[],...
    'innovations',[],...
    'magFuseMethod',[]);

% initialise the state covariance matrix
[P,Q0] = InitCovariance(param,dt_imu_avg,1,gps_data);

% Precompute the lower diagonal Cholesky decomposition for the input noise
% covariance matrix because it is time invariant
sQ = chol(Q0,'lower');

%% Main Loop

% control flags
gps_use_started = boolean(false);
gps_use_blocked = boolean(false);
viso_use_blocked = boolean(false);
flow_use_blocked = boolean(false);
sigmaPointsAreStale = boolean(false);

% array access index variables
imuIndex = imu_start_index; % incremetns with each IMU sample
logIndex = 0; % Increments each time otuptu data is logged
predIndex = 0; % Increments each minor prediction cycle and resets on each major
last_gps_index = 0;
gps_fuse_index = 0;
last_baro_index = 0;
baro_fuse_index = 0;
last_mag_index = 0;
mag_fuse_index = 0;
last_flow_index = 0;
flow_fuse_index = 0;
last_viso_index = 0;
viso_fuse_index = 0;

dtSum = 0; % accumulated time step of the filter (sec)
dtSumLog = 0; %accumulated time step used to thin log data (sec)
dtSumMinor = 0; % accumulated time step used prediction minor updates (sec)

output.magFuseMethod = param.fusion.magFuseMethod;
range = 0.1;

% variables used to control dead-reckoning timeout
last_drift_constrain_time = - param.control.velDriftTimeLim;
last_synthetic_velocity_fusion_time = 0;
last_valid_range_time = - param.fusion.rngTimeout;
local_time = imu_data.time_us(1)*1e-6;
last_known_pos = [0 0];
tempCount = 0;
for index = 1:indexStop
    
    % read IMU measurements
    delta_angle(:,1) = imu_data.del_ang(imuIndex,:);
    delta_velocity(:,1) = imu_data.del_vel(imuIndex,:);
    dtImu = 0.5 * (imu_data.accel_dt(imuIndex) + imu_data.gyro_dt(imuIndex));
    imuIndex = imuIndex+1;
    dtSum = dtSum + dtImu;
    dtSumLog = dtSumLog + dtImu;
    
    % check for time slips and increase dtImu value to match
    new_local_time = imu_data.time_us(imuIndex)*1e-6;
    if ((new_local_time - local_time) > 2.0 * dtImu)
        dtImu = min((new_local_time - local_time) , 1.0);
    end
    local_time = new_local_time;
    
    % initialise a new major update cycle by updating the covariance matrix
    % and recalculating the sigma points
    if (predIndex == 0)
        
        % Add additive state process noise variances for all stationary
        % process models accumulated since the last major update.
        % These are  IMU delta angle and delta velocity biases, earth and
        % body mag fields
        dAngBiasSigma = dtImu*dtSumMinor*param.prediction.dAngBiasPnoise;
        dVelBiasSigma = dtImu*dtSumMinor*param.prediction.dVelBiasPnoise;
        magSigmaNED = dtSumMinor*param.prediction.magPnoiseNED;
        magSigmaXYZ = dtSumMinor*param.prediction.magPnoiseXYZ;
        dtSumMinor = 0;
        processNoiseVariance = [dAngBiasSigma*[1 1 1], dVelBiasSigma*[1 1 1], magSigmaNED*[1 1 1], magSigmaXYZ*[1 1 1]].^2;
        for i=1:12
            stateIndex = i + 9;
            P(stateIndex,stateIndex) = P(stateIndex,stateIndex) + processNoiseVariance(i);
        end
        
        % calculate an array of sigma points for the augmented state vector
        sigma_x_a = CalcSigmaPoints(x_m, P , sQ , param);
        sigmaPointsAreStale = boolean(false);
        
        % convert the attitude error vector sigma points to equivalent delta
        % quaternions
        dq(:,1) = [1;0;0;0];
        for s=2:(2*param.ukf.L+1)
            normsigmaX2 = sigma_x_a(1:3,s)'*sigma_x_a(1:3,s);
            dq(1,s) = (-param.grp.a*normsigmaX2 + param.grp.f*sqrt( param.grp.f^2 + (1-param.grp.a^2)*normsigmaX2 ) ) / ( param.grp.f^2 + normsigmaX2 );
            dq(2:4,s) = (param.grp.a + dq(1,s)) * sigma_x_a(1:3,s) / param.grp.f;
        end
        
        
        
        % Apply the delta quaternions to the previous estimate to calculate the
        % quaternion sigma points. Although we could propagate the delta angles
        % through the vehicle state prediction, it is more accurate to use
        % quaternions and convert back to a set of GRP attitude error vectors
        % when the covariance information needs to be extracted.
        sigmaQuat(:,1) = q_m;
        for s=2:(2*param.ukf.L+1)
            sigmaQuat(:,s) = QuatMult(sigmaQuat(:,1),dq(:,s));
        end
        
        if (tempCount == 0)
            disp(sigmaQuat(:,2));
            tempCount = tempCount+1;
        end
    end
    
    % Minor step -  Propagate each sigma point through the vehicle state prediction
    % equations but do not update the covariance matrix and mean state
    % estimate
    for s = 1:(2*param.ukf.L+1)
        % augment the delta angles and velocitities with the IMU noise
        % states
        delAngAug = delta_angle + sigma_x_a(24:26,s);
        delVelAug = delta_velocity + sigma_x_a(27:29,s);
        
        % predict sigma points forward
        x_a = sigma_x_a(:,s);
        x_q = sigmaQuat(:,s);
        [x_q, x_a] = PredictStates(x_q, x_a, delAngAug, delVelAug, imu_data.accel_dt(imuIndex), gravity, gps_data.refLLH(1,1)*deg2rad);
        sigma_x_a(:,s) = x_a;
        sigmaQuat(:,s) = x_q;
    end
    
    % Record the number of minor steps and accumulated time since the last
    % major step
    dtSumMinor = dtSumMinor + dtImu;
    predIndex = predIndex + 1;
    
    %  Perform a major update step, update the covariance matrix and
    %  state vector and fuse any pending measurements.
    if (predIndex >= param.prediction.downSampleRatio)
        predIndex = 0;
        
        % Store the quaternion state estimate
        % The attitude error is zero mean by definition so we do not need to
        % calculate the mean from the sigma points
        q_m = sigmaQuat(:,1);
        
        % Convert propagated quaternions to delta quaternions around the
        % expected value
        dq(:,1) = [1;0;0;0];
        for s= 2:(2*param.ukf.L+1)
            dq(:,s) = QuatMult(sigmaQuat(:,1).*[1;-1;-1;-1] , sigmaQuat(:,s));
        end
        
        % Convert error quaternions to attitude error vector
        sigma_x_a(1:3,1) = [0;0;0];
        for s=2:(2*param.ukf.L+1)
            sigma_x_a(1:3,s) = param.grp.f * dq(2:4,s)/( param.grp.a + dq(1,s) );
        end
        
        % Calculate mean of predicted states from the sigma points
        x_m = sigma_x_a(1:param.ukf.nP,:) * param.ukf.wm;
        
        % Update the kinematic states (vel and pos)
        % The others use a stationary state  prediction model and will not
        % change mean value during the prediction step
        x_m(4:9) = x_m(4:9);
        
        % constrain states
        x_m  = ConstrainStates(x_m,dt_imu_avg);
        
        % Calculate covariance of predicted vehicle states from the sigma points using
        % the unscented transform
        P = zeros(param.ukf.nP,param.ukf.nP);
        for i=1:2*param.ukf.L+1
            P = P + param.ukf.wc(i)*(sigma_x_a(1:param.ukf.nP,i) - x_m)*(sigma_x_a(1:param.ukf.nP,i) - x_m)';
        end
        
        % Keep the covariance matrix positive definite
        P = fix_covariance(P);
        
        % output state data
        if (dtSumLog >= 0.1)
            logIndex = logIndex + 1;
            dtSumLog = 0;
            output.time_lapsed(logIndex) = local_time;
            output.euler_angles(logIndex,:) = QuatToEul(q_m)';
            output.velocity_NED(logIndex,:) = x_m(4:6)';
            output.position_NED(logIndex,:) = x_m(7:9)';
            output.gyro_bias(logIndex,:) = x_m(10:12)';
            output.accel_bias(logIndex,:) = x_m(13:15)';
            output.mag_NED(logIndex,:) = x_m(16:18);
            output.mag_XYZ(logIndex,:) = x_m(19:21);
            output.wind_NE(logIndex,:) = x_m(22:23);
            
            % output covariance data
            for i=1:23
                output.state_variances(logIndex,i) = P(i,i);
            end
            
            % output equivalent euler angle variances
            % TODO method to translate attitude error vector variances to Euler
            % angle variances
            for i=1:3
                output.euler_variances(logIndex,i) = P(i,i);
            end
            
        end
        
        % Get most recent GPS data that had fallen behind the fusion time horizon
        latest_gps_index = find((gps_data.time_us - 1e6 * param.fusion.gpsTimeDelay) < imu_data.time_us(imuIndex), 1, 'last' );
        
        % Check if GPS use is being blocked by the user
        if ((local_time < param.control.gpsOnTime) && (local_time > param.control.gpsOffTime))
            gps_use_started = false;
            gps_use_blocked = true;
        else
            gps_use_blocked = false;
        end
        
        % If we haven't started using GPS, check that the quality is sufficient before aligning the position and velocity states to GPS
        if (~isempty(latest_gps_index) && ~gps_use_started && ~gps_use_blocked)
            if ((gps_data.spd_error(latest_gps_index) < param.control.gpsSpdErrLim) && (gps_data.pos_error(latest_gps_index) < param.control.gpsPosErrLim))
                states(5:7) = gps_data.vel_ned(latest_gps_index,:);
                states(8:9) = gps_data.pos_ned(latest_gps_index,1:2);
                gps_use_started = true;
                last_drift_constrain_time = local_time;
            end
        end
        
        % Fuse GPS data when available if GPS use has started
        if (~isempty(latest_gps_index) && gps_use_started && ~gps_use_blocked && (latest_gps_index > last_gps_index))
            last_gps_index = latest_gps_index;
            gps_fuse_index = gps_fuse_index + 1;
            
            % fuse NED GPS velocity
            [q_m, x_m, P, velInnov,velInnovVar,sigmaPointsAreStale] = FuseVelocity(...
                sigma_x_a, ...
                q_m, ...
                x_m, ...
                P, ...
                gps_data.vel_ned(latest_gps_index,1:3), ...
                gps_data.spd_error(latest_gps_index), ...
                param, ...
                sigmaPointsAreStale);
            
            if (sigmaPointsAreStale)
                last_drift_constrain_time = local_time;
            end
            
            % constrain states
            x_m  = ConstrainStates(x_m,dt_imu_avg);
            
            % data logging
            output.innovations.vel_time_lapsed(gps_fuse_index) = local_time;
            output.innovations.vel_innov(gps_fuse_index,:) = velInnov';
            output.innovations.vel_innov_var(gps_fuse_index,:) = velInnovVar';
            
            % Update the sigma points using the current covariance
            if (sigmaPointsAreStale)
                sigma_x_a = CalcSigmaPoints(x_m, P , sQ , param);
                sigmaPointsAreStale = boolean(false);
            end
            
            % fuse NE GPS position
            [q_m, x_m, P, posInnov,posInnovVar,sigmaPointsAreStale] = FusePosition(...
                sigma_x_a, ...
                q_m, ...
                x_m, ...
                P, ...
                gps_data.pos_ned(latest_gps_index,1:2), ...
                gps_data.pos_error(latest_gps_index), ...
                param, ...
                sigmaPointsAreStale);
            
            if (sigmaPointsAreStale)
                last_drift_constrain_time = local_time;
            end
            
            % constrain states
            x_m  = ConstrainStates(x_m,dt_imu_avg);
            
            % data logging
            output.innovations.pos_time_lapsed(gps_fuse_index) = local_time;
            output.innovations.posInnov(gps_fuse_index,:) = posInnov';
            output.innovations.posInnovVar(gps_fuse_index,:) = posInnovVar';
        end
        
        % Check if drift is being corrected by some form of aiding and if not, fuse in a zero position measurement at 5Hz to prevent states diverging
        if ((local_time - last_drift_constrain_time) > param.control.velDriftTimeLim)
            if (gps_use_started)
                last_known_pos = [states(8),states(9)];
                states(5:7) = [0,0,0];
                gps_use_started = 'false';
            end
            
            if ((local_time - last_synthetic_velocity_fusion_time) > 0.2)
                
                % Update the sigma points using the current covariance
                if (sigmaPointsAreStale)
                    sigma_x_a = CalcSigmaPoints(x_m, P , sQ , param);
                    sigmaPointsAreStale = boolean(false);
                end
                
                % fuse static position at last known location
                % origin moves with us until we start using GPS
                [q_m, x_m, P, ~,~,sigmaPointsAreStale] = FusePosition(...
                    sigma_x_a, ...
                    q_m, ...
                    x_m, ...
                    P, ...
                    last_known_pos, ...
                    param.fusion.staticPosErr, ...
                    param, ...
                    sigmaPointsAreStale);
                
                % constrain states
                x_m  = ConstrainStates(x_m,dt_imu_avg);
                
                last_synthetic_velocity_fusion_time = local_time;
            end
        end
        
        % Fuse new Baro data that has fallen beind the fusion time horizon
        latest_baro_index = find((baro_data.time_us - 1e6 * param.fusion.baroTimeDelay) < imu_data.time_us(imuIndex), 1, 'last' );
        if (latest_baro_index > last_baro_index)
            last_baro_index = latest_baro_index;
            baro_fuse_index = baro_fuse_index + 1;
            
            % Update the sigma points using the current covariance
            if (sigmaPointsAreStale)
                sigma_x_a = CalcSigmaPoints(x_m, P , sQ , param);
                sigmaPointsAreStale = boolean(false);
            end
            
            % fuse baro height
            [q_m, x_m, P, hgtInnov, hgtInnovVar] = FuseBaroHeight(...
                sigma_x_a, ...
                q_m, ...
                x_m, ...
                P, ...
                baro_data.height(latest_baro_index), ...
                param, ...
                sigmaPointsAreStale);
            
            % constrain states
            x_m  = ConstrainStates(x_m,dt_imu_avg);
            
            % data logging
            output.innovations.hgt_time_lapsed(baro_fuse_index) = local_time;
            output.innovations.hgtInnov(baro_fuse_index) = hgtInnov;
            output.innovations.hgtInnovVar(baro_fuse_index) = hgtInnovVar;
        end
        
        % Fuse new mag data that has fallen behind the fusion time horizon
        latest_mag_index = find((mag_data.time_us - 1e6 * param.fusion.magTimeDelay) < imu_data.time_us(imuIndex), 1, 'last' );
        if (latest_mag_index > last_mag_index)
            last_mag_index = latest_mag_index;
            mag_fuse_index = mag_fuse_index + 1;
            
            % output magnetic field length to help with diagnostics
            output.innovations.magLength(mag_fuse_index) = sqrt(dot(mag_data.field_ga(latest_mag_index,:),mag_data.field_ga(latest_mag_index,:)));
            
            % Update the sigma points using the current covariance
            if (sigmaPointsAreStale)
                sigma_x_a = CalcSigmaPoints(x_m, P , sQ , param);
                sigmaPointsAreStale = boolean(false);
            end
            
            % fuse magnetometer data
            if (param.fusion.magFuseMethod == 0 || param.fusion.magFuseMethod == 1)
                [q_m, x_m, P, magInnov, magInnovVar, sigmaPointsAreStale] = FuseMagnetometer( ...
                    sigma_x_a, ...
                    q_m, ...
                    x_m, ...
                    P, ...
                    mag_data.field_ga(latest_mag_index,:), ...
                    param);
                
                % data logging
                output.innovations.mag_time_lapsed(mag_fuse_index) = local_time;
                output.innovations.magInnov(mag_fuse_index,:) = magInnov;
                output.innovations.magInnovVar(mag_fuse_index,:) = magInnovVar;
                
                if (param.fusion.magFuseMethod == 1)
                    % Update the sigma points using the current covariance
                    if (sigmaPointsAreStale)
                        sigma_x_a = CalcSigmaPoints(x_m, P , sQ , param);
                        sigmaPointsAreStale = boolean(false);
                    end
                    
                    % fuse in the local declination value
                    [q_m, x_m, P, sigmaPointsAreStale] = FuseMagDeclination( ...
                        sigma_x_a, ...
                        q_m, ...
                        x_m, ...
                        P, ...
                        param, ...
                        sigmaPointsAreStale);
                    
                end
                
            elseif (param.fusion.magFuseMethod == 2)
                % Update the sigma points using the current covariance
                if (sigmaPointsAreStale)
                    sigma_x_a = CalcSigmaPoints(x_m, P , sQ , param);
                    sigmaPointsAreStale = boolean(false);
                end
                
                % fuse magnetomer data as a single magnetic heading measurement
                [q_m, x_m, P, hdgInnov, hdgInnovVar, sigmaPointsAreStale] = FuseMagHeading( ...
                    sigma_x_a, ...
                    q_m, ...
                    x_m, ...
                    P, ...
                    mag_data.field_ga(latest_mag_index,:)', ...
                    param, ...
                    sigmaPointsAreStale);
                
                % log data
                output.innovations.mag_time_lapsed(mag_fuse_index) = local_time;
                output.innovations.hdgInnov(mag_fuse_index) = hdgInnov;
                output.innovations.hdgInnovVar(mag_fuse_index) = hdgInnovVar;
                
            end
            
        end
        
        % Check if optical flow use is being blocked by the user
        if ((local_time < param.control.flowOnTime) && (local_time > param.control.flowOffTime))
            flow_use_blocked = true;
        else
            flow_use_blocked = false;
        end
        
        % Attempt to use optical flow and range finder data if available and not blocked
        if (flowDataPresent && ~flow_use_blocked)
            
            % Get latest range finder data and gate to remove dropouts
            last_range_index = find((rng_data.time_us - 1e6 * param.fusion.rangeTimeDelay) < imu_data.time_us(imuIndex), 1, 'last' );
            if (rng_data.dist(last_range_index) < param.fusion.rngValidMax)
                range = max(rng_data.dist(last_range_index) , param.fusion.rngValidMin);
                last_valid_range_time = local_time;
            end
            
            % Fuse optical flow data that has fallen behind the fusion time horizon if we have a valid range measurement
            latest_flow_index = find((flow_data.time_us - 1e6 * param.fusion.flowTimeDelay) < imu_data.time_us(imuIndex), 1, 'last' );
            
            if (~isempty(latest_flow_index) && (latest_flow_index > last_flow_index) && ((local_time - last_valid_range_time) < param.fusion.rngTimeout))
                % Update the sigma points using the current covariance
                if (sigmaPointsAreStale)
                    sigma_x_a = CalcSigmaPoints(x_m, P , sQ , param);
                    sigmaPointsAreStale = boolean(false);
                end
                
                last_flow_index = latest_flow_index;
                flow_fuse_index = flow_fuse_index + 1;
                
                % fuse flow data
                flowRate = [flow_data.flowX(latest_flow_index);flow_data.flowY(latest_flow_index)];
                bodyRate = [flow_data.bodyX(latest_flow_index);flow_data.bodyY(latest_flow_index)];
                [q_m, x_m, P, flowInnov, flowInnovVar, sigmaPointsAreStale] = FuseOpticalFlow( ...
                    sigma_x_a, ...
                    q_m, ...
                    x_m, ...
                    P, ...
                    flowRate,...
                    bodyRate, ...
                    range, ...
                    param, ...
                    sigmaPointsAreStale);
                
                if (sigmaPointsAreStale)
                    last_drift_constrain_time = local_time;
                end
                
                % data logging
                output.innovations.flow_time_lapsed(flow_fuse_index) = local_time;
                output.innovations.flowInnov(flow_fuse_index,:) = flowInnov;
                output.innovations.flowInnovVar(flow_fuse_index,:) = flowInnovVar;
                
            end
            
        end
        
        % Check if optical flow use is being blocked by the user
        if ((local_time < param.control.visoOnTime) && (local_time > param.control.visoOffTime))
            viso_use_blocked = true;
        else
            viso_use_blocked = false;
        end
        
        % attempt to use ZED camera visual odmetry data if available and not blocked
        if (visoDataPresent && ~viso_use_blocked)
            
            % Fuse ZED camera body frame odmometry data that has fallen behind the fusion time horizon
            latest_viso_index = find((viso_data.time_us - 1e6 * param.fusion.bodyVelTimeDelay) < imu_data.time_us(imuIndex), 1, 'last' );
            if (latest_viso_index > last_viso_index)
                last_viso_index = latest_viso_index;
                viso_fuse_index = viso_fuse_index + 1;
                
                % convert delta positon measurements to velocity
                relVelBodyMea = [viso_data.dPosX(latest_viso_index);viso_data.dPosY(latest_viso_index);viso_data.dPosZ(latest_viso_index)]./viso_data.dt(latest_viso_index);
                
                % convert quality metric to equivalent observation error
                % (this is a guess)
                quality = viso_data.qual(latest_viso_index);
                bodyVelError = param.fusion.bodyVelErrorMin * quality + param.fusion.bodyVelErrorMax * (1 - quality);
                
                % fuse measurements
                [q_m, x_m, P, bodyVelInnov, bodyVelInnovVar, sigmaPointsAreStale] = FuseBodyVel( ...
                    sigma_x_a, ...
                    q_m, ...
                    x_m, ...
                    P, ...
                    relVelBodyMea, ...
                    bodyVelError, ...
                    param, ...
                    sigmaPointsAreStale);
                
                if (sigmaPointsAreStale)
                    last_drift_constrain_time = local_time;
                end
                
                % data logging
                output.innovations.bodyVel_time_lapsed(viso_fuse_index) = local_time;
                output.innovations.bodyVelInnov(viso_fuse_index,:) = bodyVelInnov;
                output.innovations.bodyVelInnovVar(viso_fuse_index,:) = bodyVelInnovVar;
                
            end
            
        end
        
    end
    
    % update average delta time estimate
    output.dt = dtSum / index;
    
end

end