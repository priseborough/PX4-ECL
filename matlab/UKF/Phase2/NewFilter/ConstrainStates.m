function [states]  = ConstrainStates(states,dt_imu_avg)

% attitude error state is zero mean by definition
states(1:3,1) = [0;0;0];
    
% constrain gyro bias states
limit = 5.0*pi/180*dt_imu_avg;
for i=10:12
    if (states(i) > limit)
        states(i) = limit;
    elseif (states(i) < -limit)
        states(i) = -limit;
    end
end

% constrain accel bias states
limit = 0.5*dt_imu_avg;
for i=13:15
    if (states(i) > limit)
        states(i) = limit;
    elseif (states(i) < -limit)
        states(i) = -limit;
    end
end

end