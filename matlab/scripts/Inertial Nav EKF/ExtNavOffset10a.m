% This script requires the Matlab symbolic toolbox and takes
% ~30 minutes to run

% Derivation of state and covariance prediction equations for the
% estimation of the position scale factor, position offset and orientation
% offset between an external nav systems world frame and the autopilot
% system navigation frame where the world frame distance units have an arbitrary
% scale factor relative to SI units.

% 10 state architecture.

% Author: Paul Riseborough

% State vector:

% XYZ velocity in world frame (length/sec)
% XYZ position in world frame (length)
% XYZ position of the workd frame origin in navigation frame (m)
% Scale factor that converts from nav to word frame length units. 

% Observations:

% direct state observation of XYZ position in world frame (m)

% Time varying parameters (control input):

% quaternions decribing rotation from world frame to navigation frame
% XYZ accceleration in navigation frame (m/s^2)

clear all;
reset(symengine);

syms pos_x_w pos_y_w pos_z_w 'real' % world frame position (m) * scale
syms vel_x_w vel_y_w vel_z_w 'real' % world frame velocity (m/s) * scale
syms acc_x_w acc_y_w acc_z_w 'real' % world frame acceleration (m/s/s) * scale
syms pos_x_n pos_y_n pos_z_n 'real' % navigation frame position (m)
syms vel_x_n vel_y_n vel_z_n 'real' % navigation frame velocity (m/s)
syms acc_x_n acc_y_n acc_z_n 'real' % navigation frame acceleration (m/s/s)
syms scale 'real' % scale factor to convert from nav frame to world frame
syms ogn_x ogn_y ogn_z 'real' % position of world frame origin in nav frame
syms q0 q1 q2 q3 'real' % quaternions decribing rotation of NED navigation frame relative to world frame
syms dt 'real' % time step (sec)
syms acc_var 'real' % variance of nav frame accel measurements (m/s/s)^2

% define the quaternion rotation vector for the state estimate
quat = [q0;q1;q2;q3];

% derive the world to nav direction cosine matrix
Tnw = Quat2Tbn(quat);

% The nav frame accelerations are treated as control inputs
accel_nav = [acc_x_n;acc_y_n;acc_z_n];

% rotate nav accel into word frame and integrate to get velocity and
% position
vel_nav = [vel_x_n;vel_y_n;vel_z_n];
vel_nav_new = vel_nav + accel_nav * dt;
pos_nav = [pos_x_n;pos_y_n;pos_z_n];
origin = [ogn_x;ogn_y;ogn_z];
pos_nav_new = pos_nav + vel_nav * dt;

% static process models
scale_new = scale;
origin_new = origin;

state_vec = [vel_nav;pos_nav;origin;scale];
state_vec_new = [vel_nav_new;pos_nav_new;origin_new;scale_new];

nStates=numel(state_vec);

% derive the state transition matrix
F = jacobian(state_vec_new, state_vec);

% define a symbolic covariance matrix using strings to represent 
% '_l_' to represent '( '
% '_c_' to represent ,
% '_r_' to represent ')' 
% these can be substituted later to create executable code
for rowIndex = 1:nStates
    for colIndex = 1:nStates
        eval(['syms OP_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['P(',num2str(rowIndex),',',num2str(colIndex), ') = OP_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end

% calculate the control influence matrix from nav frame accel input errors
% to states
G = jacobian(state_vec_new, accel_nav);

% derive the state error matrix
distMatrix = diag([acc_var acc_var acc_var]);
Q = G*distMatrix*transpose(G);
[Q,SQ]=OptimiseAlgebra(Q,'SQ');

% Derive the predicted covariance matrix using the standard equation
PP = F*P*transpose(F) + Q;

% Collect common expressions to optimise processing
[PP,SPP]=OptimiseAlgebra(PP,'SPP');

%% observation model of position in world frame
% correct for offset and scale factor
pos_world = scale * (pos_nav - origin);

% correct for world frame orientation
pos_world = Tnw * pos_world;

% calculate jacobians and Kalman gains for world frame position observation
syms R_POS 'real' % position observation variance

H_X = jacobian(pos_world(1),state_vec); % measurement Jacobian
[H_X,SH_X]=OptimiseAlgebra(H_X,'SH_X'); % optimise processing
K_X = (P*transpose(H_X))/(H_X*P*transpose(H_X) + R_POS);
[K_X,SK_X]=OptimiseAlgebra(K_X,'SK_X'); % Kalman gain vector

H_Y = jacobian(pos_world(2),state_vec); % measurement Jacobian
[H_Y,SH_Y]=OptimiseAlgebra(H_Y,'SH_Y'); % optimise processing
K_Y = (P*transpose(H_Y))/(H_Y*P*transpose(H_Y) + R_POS);
[K_Y,SK_Y]=OptimiseAlgebra(K_Y,'SK_Y'); % Kalman gain vector

H_Z = jacobian(pos_world(3),state_vec); % measurement Jacobian
[H_Z,SH_Z]=OptimiseAlgebra(H_Z,'SH_Z'); % optimise processing
K_Z = (P*transpose(H_Z))/(H_Z*P*transpose(H_Z) + R_POS);
[K_Z,SK_Z]=OptimiseAlgebra(K_Z,'SK_Z'); % Kalman gain vector

%% convert to C code

fileName = strcat('SymbolicOutput',int2str(nStates),'.mat');
save(fileName);
SaveScriptCode(nStates);
ConvertToM(nStates);
ConvertToC(nStates);
