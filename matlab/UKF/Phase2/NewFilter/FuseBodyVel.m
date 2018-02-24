function [...
    q_m, ... quaternion state estimate
    x_m, ... % vehicle state estimate
    P, ... % vehicle state covariance matrix
    innovation,... % XYZ magnetomer innovations (Gauss)
    varInnov, ... % XYZ magnetomer innovation variances (Gauss^2)
    sigmaPointStale] ... % true when the covariance matrix has been modifed and the sigma points are stale
    = FuseBodyVel( ...
    sigma_x_a, ... % array of sigma points for the augmented state vector
    q_m, ... quaternion state estimate
    x_m, ... % vehicle state estimate
    P, ... % vehicle state covariance matrix
    relVelBodyMea, ... % XYZ velocity measured by the camera (m/sec)
    bodyVelError, ... % velocity observation error (m/sec)
    param, ...
    sigmaPointStale)

varInnov = zeros(1,3);

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
psi_m = zeros(3,(2*param.ukf.L+1));
for s = 1:(2*param.ukf.L+1)
    Tnb = transpose(Quat2Tbn(sigmaQuat(:,s)));
    psi_m(:,s) = Tnb * sigma_x_a(4:6,s);    
end

% Calculate mean of predicted output
y_m = psi_m * param.ukf.wm;

% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = diag([1;1;1]*bodyVelError^2);
Pxy = zeros(param.ukf.nP,3);
for i = 1:(2*param.ukf.L+1)
    Pyy = Pyy + param.ukf.wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
    Pxy = Pxy + param.ukf.wc(i)*(sigma_x_a(1:param.ukf.nP,i) - x_m)*(psi_m(:,i) - y_m)';
end

% Calculate the innovation and innovation variance
varInnov(1) = Pyy(1,1);
varInnov(2) = Pyy(2,2);
varInnov(3) = Pyy(3,3);
innovation = y_m - relVelBodyMea;

% calculate the normalised innovation
invPyy = inv(Pyy);
innovNorm2 = innovation' * invPyy * innovation;

% Perform a chi2 innovation consistency check
if (innovNorm2 > param.fusion.bodyVelGate^2)
    return;
end

% Calculate Kalman gain
K = Pxy * invPyy; % Calculate Kalman gain

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