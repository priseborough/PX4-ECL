function [...
    q_m, ... quaternion state estimate
    x_m, ... % vehicle state estimate
    P, ... % vehicle state covariance matrix
    innovation,... % D position innovations (m)
    varInnov, ... % D position innovation variance (m^2)
    sigmaPointStale] ... % true when the covariance matrix has been modifed and the sigma points are stale
    = FuseBaroHeight( ...
    sigma_x_a, ... % array of sigma points for the augmented state vector
    q_m, ... quaternion state estimate
    x_m, ... % vehicle state estimate
    P, ... % vehicle state covariance matrix
    measHgt, ... % height measurement (m)
    param, ... % parameter struct
    sigmaPointStale)

% Propagate each sigma-point through observation transformation
psi_m(1,:) = sigma_x_a(9,:); % direct observation of vertical position

% y_m = psi_m * wm; % Calculate mean of predicted output
y_m = x_m(9); % same as psi_m * wm because we are directly observing state

% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = param.fusion.baroHgtNoise^2;
Pxy = zeros(param.ukf.nP,1);
for i = 1:(2*param.ukf.L+1)
    Pyy = Pyy + param.ukf.wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
    Pxy = Pxy + param.ukf.wc(i)*(sigma_x_a(1:param.ukf.nP,i) - x_m)*(psi_m(:,i) - y_m)';
end

% Calculate the innovation and innovation variance
varInnov = Pyy;
innovation = y_m + measHgt;

% Apply an innovation consistency check
if (innovation^2 / (param.fusion.baroHgtGate^2 * varInnov)) > 1.0
    return;
end

% Calculate Kalman gain
K = Pxy/Pyy; % Calculate Kalman gain

% Update the state estimate
x_m = x_m - K*innovation;

% Update the covariance estimate
P = P - K*Pyy*K'; % Update covariance estimate
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