function [...
    q_m, ... quaternion state estimate
    x_m, ... % vehicle state estimate
    P, ... % vehicle state covariance matrix
    sigmaPointStale] ... % true when the covariance matrix has been modifed and the sigma points are stale
    = FuseMagDeclination( ...
    sigma_x_a, ... % array of sigma points for the augmented state vector
    q_m, ... quaternion state estimate
    x_m, ... % vehicle state estimate
    P, ... % vehicle state covariance matrix
    param, ...
    sigmaPointStale)

% calculate observation sigma points
psi_m = zeros(1,(2*param.ukf.L+1));
for s = 1:(2*param.ukf.L+1)
    % declination angle is taken from horizontal projection of earth field
    psi_m(1,s) = atan2(sigma_x_a(17,s),sigma_x_a(16,s));
end

% Calculate mean of predicted output
y_m = psi_m * param.ukf.wm;

% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = 0.5^2;
Pxy = zeros(param.ukf.nP,1);
for i = 1:(2*param.ukf.L+1)
    Pyy = Pyy + param.ukf.wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
    Pxy = Pxy + param.ukf.wc(i)*(sigma_x_a(1:param.ukf.nP,i) - x_m)*(psi_m(:,i) - y_m)';
end

% Calculate the innovation and innovation variance
innovation = y_m - param.fusion.magDeclDeg*(pi/180);

% Perform a chi2 innovation consistency check
if (innovation^2 / Pyy > param.fusion.magFieldGate^2)
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