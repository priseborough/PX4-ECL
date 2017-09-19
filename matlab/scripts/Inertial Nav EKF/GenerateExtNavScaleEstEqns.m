% This script requires the Matlab symbolic toolbox and takes
% ~5 minutes to run

% Derivation of state and covariance prediction equations for the
% estimation of the position scale factor between an external nav systems 
% world frame and the autopilot system navigation frame where the world 
% frame distance units have an arbitrary scale factor relative to SI units.

% 7 state architecture.

% Author: Paul Riseborough

% State vector:

% XYZ velocity in world frame (length/sec)
% XYZ position in world frame (length)
% Natural log of scale factor that converts from nav to word frame length units. 

% Observations:

% direct state observation of XYZ position in world frame (m)

% Time varying parameters (control input):

% XYZ accceleration in world frame (m/s^2)

clear all;
reset(symengine);

syms pos_x pos_y pos_z 'real' % world frame position in relative units (m) * scale
syms vel_x vel_y vel_z 'real' % world frame velocity in reltive units (m/s) * scale
syms delVelWorld_x delVelWorld_y delVelWorld_z 'real' % world frame acceleration in SI units (m/s/s)
syms scaleLog 'real' % natural log of scale factor to convert from nav frame to world frame
syms delT 'real' % time step (sec)
syms delVelVar 'real' % variance of world frame delta velocity increments

% The world frame accelerations are treated as control inputs
delVel = [delVelWorld_x;delVelWorld_y;delVelWorld_z];

% rotate nav accel into world frame and integrate to get velocity and
% position
vel_world = [vel_x;vel_y;vel_z];
vel_world_new = vel_world + exp(scaleLog) * delVel;
pos_world = [pos_x;pos_y;pos_z];
pos_world_new = pos_world + vel_world * delT;

% static process models
scaleLog_new = scaleLog;

state_vec = [vel_world;pos_world;scaleLog];
state_vec_new = [vel_world_new;pos_world_new;scaleLog_new];

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
        eval(['syms extNavP_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['P(',num2str(rowIndex),',',num2str(colIndex), ') = extNavP_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end

% calculate the control influence matrix from nav frame accel input errors
% to states
G = jacobian(state_vec_new, delVel);

% derive the state error matrix
distMatrix = diag([delVelVar delVelVar delVelVar]);
Q = G*distMatrix*transpose(G);
[Q,SQ]=OptimiseAlgebra(Q,'SQ');

% Derive the predicted covariance matrix using the standard equation
Pnew = F*P*transpose(F) + Q;

% Collect common expressions to optimise processing
[Pnew,SPP]=OptimiseAlgebra(Pnew,'SPP');

%% convert to C code

fileName = strcat('SymbolicOutput',int2str(nStates),'.mat');
save(fileName);
SaveScriptCode(nStates);
ConvertToM(nStates);
ConvertToC(nStates);
