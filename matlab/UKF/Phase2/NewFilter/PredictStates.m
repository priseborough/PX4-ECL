function [quat_state, states]  = PredictStates( ...
    quat_state, ...
    states, ... % previous state vector (3x1 rotVec ; 3x1 velNED ; 3x1 posNED ; 3x1 dAngBias ; 3x1 dVelBias ; 3x1 magNED ; 3x1 magXYZ ; 2x1 velWindNE)
    delAng, ... % IMU delta angle measurements, 3x1 (rad)
    delVel, ... % IMU delta velocity measurements 3x1 (m/s)
    dt, ... % accumulation time of the IMU measurement (sec)
    gravity, ... % acceleration due to gravity (m/s/s)
    latitude) % WGS-84 latitude (rad)

% define persistent variables for previous delta angle and velocity which
% are required for sculling and coning error corrections
persistent prevDelAng;
if isempty(prevDelAng)
    prevDelAng = delAng;
end

persistent prevDelVel;
if isempty(prevDelVel)
    prevDelVel = delVel;
end

persistent Tbn_prev;
if isempty(Tbn_prev)
    Tbn_prev = Quat2Tbn(quat_state(1:4));
end

% Remove sensor bias errors
delAng = delAng - states(10:12);
delVel = delVel - states(13:15);

% Correct delta velocity for rotation and skulling
% Derived from Eqn 25 of:
% "Computational Elements For Strapdown Systems"
% Savage, P.G.
% Strapdown Associates
% 2015, WBN-14010
correctedDelVel= delVel + ...
    0.5*cross(prevDelAng + delAng , prevDelVel + delVel) + 1/6*cross(prevDelAng + delAng , cross(prevDelAng + delAng , prevDelVel + delVel)) +  1/12*(cross(prevDelAng , delVel) + cross(prevDelVel , delAng));

% Calculate earth delta angle spin vector
delAngEarth_NED(1,1) = 0.000072921 * cos(latitude) * dt;
delAngEarth_NED(2,1) = 0.0;
delAngEarth_NED(3,1) = -0.000072921 * sin(latitude) * dt;

% Apply corrections for coning errors and earth spin rate
% Coning correction from :
% "A new strapdown attitude algorithm", 
% R. B. MILLER, 
% Journal of Guidance, Control, and Dynamics
% July, Vol. 6, No. 4, pp. 287-291, Eqn 11 
correctedDelAng   = delAng - 1/12*cross(prevDelAng , delAng) - transpose(Tbn_prev)*delAngEarth_NED;

% Save current measurements
prevDelAng = delAng;
prevDelVel = delVel;

% Convert the rotation vector to its equivalent quaternion
deltaQuat = RotToQuat(correctedDelAng);

% Update the quaternions by rotating from the previous attitude through
% the delta angle rotation quaternion
quat_state = QuatMult(quat_state,deltaQuat);

% Normalise the quaternions
quat_state = NormQuat(quat_state);

% Calculate the body to nav cosine matrix
Tbn = Quat2Tbn(quat_state);
Tbn_prev = Tbn;

% transform body delta velocities to delta velocities in the nav frame
delVelNav = Tbn * correctedDelVel + [0;0;gravity]*dt;

% take a copy of the previous velocity
prevVel = states(4:6);

% Sum delta velocities to get the velocity
states(4:6) = states(4:6) + delVelNav(1:3);

% integrate the velocity vrctor to get the position using trapezoidal
% integration
states(7:9) = states(7:9) + 0.5 * dt * (prevVel + states(4:6));

end