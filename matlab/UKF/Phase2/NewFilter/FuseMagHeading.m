function [...
    q_m, ... quaternion state estimate
    x_m, ... % vehicle state estimate
    P, ... % vehicle state covariance matrix
    innovation,... % XYZ magnetomer innovations (Gauss)
    varInnov, ... % XYZ magnetomer innovation variances (Gauss^2)
    sigmaPointStale] ... % true when the covariance matrix has been modifed and the sigma points are stale
    = FuseMagHeading( ...
    sigma_x_a, ... % array of sigma points for the augmented state vector
    q_m, ... quaternion state estimate
    x_m, ... % vehicle state estimate
    P, ... % vehicle state covariance matrix
    measMagXYZ, ... % XYZ field measurements (Gauss)
    param, ...
    sigmaPointStale)

% convert the attitude error vector sigma points to equivalent delta
% quaternions
dq(:,1) = [1;0;0;0];
for s=2:(2*param.ukf.L+1)
    normsigmaX2 = sigma_x_a(1:3,s)'*sigma_x_a(1:3,s);
    dq(1,s) = (-param.grp.a*normsigmaX2 + param.grp.f*sqrt( param.grp.f^2 + (1-param.grp.a^2)*normsigmaX2 ) ) / ( param.grp.f^2 + normsigmaX2 );
    dq(2:4,s) = (param.grp.a + dq(1,s)) * sigma_x_a(1:3,s) / param.grp.f;
end

% Apply the delta quaternions to the previous estimate to calculate the
% quaternion sigma points.
sigmaQuat(:,1) = q_m;
for s=2:(2*param.ukf.L+1)
    sigmaQuat(:,s) = QuatMult(sigmaQuat(:,1),dq(:,s));
end

% calculate observation sigma points
psi_m = zeros(1,(2*param.ukf.L+1));
for s = 1:(2*param.ukf.L+1)
    % calculate the rotation matrix from earth to body frame
    Tnb = transpose(Quat2Tbn(sigmaQuat(:,s)));
    
    % rotate earth field into body frame and add bias to get predicted
    % NED field
    magNED = sigma_x_a(19:21,s) + Tnb * sigma_x_a(16:18,s);
    
    % convert predicted field to a predicted declination angle
    psi_m(1,s) = atan2(magNED(2) , magNED(1));
    
end

% Calculate mean of predicted output
y_m = psi_m * param.ukf.wm;

% rotate measured field into earth frame
magNED = Quat2Tbn(q_m) * measMagXYZ;
y_obs = atan2(magNED(2) , magNED(1));

% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = param.fusion.magHdgError^2;
Pxy = zeros(param.ukf.nP,1);
for i = 1:(2*param.ukf.L+1)
    Pyy = Pyy + param.ukf.wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
    Pxy = Pxy + param.ukf.wc(i)*(sigma_x_a(1:param.ukf.nP,i) - x_m)*(psi_m(:,i) - y_m)';
end

% Calculate the innovation and innovation variance
varInnov = Pyy;
innovation = y_m - y_obs;

% Perform a chi2 innovation consistency check
if (innovation^2 / Pyy > param.fusion.magHdgGate^2)
    return;
end

% Calculate Kalman gain
K = Pxy / Pyy; % Calculate Kalman gain

% Update the state estimate
x_m = x_m - K * innovation;

% Update the covariance estimate
P = P - K * Pyy * K'; % Update covariance estimate
sigmaPointStale = boolean(true);

% Keep the covariance matrix positive definite
P = fix_covariance(P);

% Calculate the quaternion correction from the attitude error
% vector
normdp2 = x_m(1:3)'*x_m(1:3);
dq(1,1) = ( -param.grp.a * normdp2 + param.grp.f * sqrt( param.grp.f^2 + (1-param.grp.a^2) * normdp2 ) ) / ( param.grp.f^2 + normdp2 );
dq(2:4,1) = x_m(1:3) * ( param.grp.a + dq(1,1) ) / param.grp.f;

% apply the correction to the stored quaternion state
q_m = QuatMult(q_m , dq(:,1));

% renormalise
q_m = NormQuat(q_m);

end