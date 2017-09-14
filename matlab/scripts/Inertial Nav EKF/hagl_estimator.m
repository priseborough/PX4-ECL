% script to derive EKF equations for height above ground estimator

clear all;
reset(symengine);

syms dt 'real'
syms vn ve vd 'real'; % NED velocity (m/s)
syms pd 'real'; % vertical position of vehicle along D axis (m)
syms ptd 'real'; % vertical position of terrain along D axis (m)
syms q0 q1 q2 q3 'real' % quaternions defining rotation from NED earth frame to XYZ body frame
syms dVelN dVelE dVelD 'real'; % NED accel integrated over dt (m/s)

% derive the truth body to nav direction cosine matrix
Tbn = Quat2Tbn(quat);

% process equations
velNew = [vn;ve;vd] + [dVelN;dVelE;dVelD];
pdNew = pd + vd * dt;
ptdNew = ptd;

% define prior and predicted state vectors
stateVector = [vn;ve;vd;pd;ptd];
stateVectorNew = [velNew;pdNew;ptdNew];

% define predicted LOS rate from optical flow sensor
% range is defined as distance from camera focal point to centre of sensor fov
range = (ptd - pd) / Tbn(3,3);

% calculate relative velocity in body frame
relVelBody = transpose(Tbn)*[vn;ve;vd];

% divide by range to get predicted angular LOS rates relative to X and Y
% axes. Note these are body angular rate motion compensated optical flow rates
losRateX = +relVelBody(2)/range;
losRateY = -relVelBody(1)/range;

% calculate the observation Jacobian for the X axis
H_LOSX = jacobian(losRateX,stateVector); % measurement Jacobian
H_LOSX = simplify(H_LOSX);
save('temp2.mat','H_LOSX');
ccode(H_LOSX,'file','H_LOSX.c');
fix_c_code('H_LOSX.c');

% calculate the observation Jacobian for the Y axis
H_LOSY = jacobian(losRateY,stateVector); % measurement Jacobian
H_LOSY = simplify(H_LOSY);
save('temp3.mat','H_LOSY');
ccode(H_LOSY,'file','H_LOSY.c');
fix_c_code('H_LOSY.c');

% calculate Kalman gain vector for the X axis
K_LOSX = (P*transpose(H_LOSX))/(H_LOSX*P*transpose(H_LOSX) + R_LOS); % Kalman gain vector
K_LOSX = simplify(K_LOSX);
ccode(K_LOSX,'file','K_LOSX.c');
fix_c_code('K_LOSX.c');

% calculate Kalman gain vector for the Y axis
K_LOSY = (P*transpose(H_LOSY))/(H_LOSY*P*transpose(H_LOSY) + R_LOS); % Kalman gain vector
K_LOSY = simplify(K_LOSY);
ccode(K_LOSY,'file','K_LOSY.c');
fix_c_code('K_LOSY.c');