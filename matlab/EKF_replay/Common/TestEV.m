clear all;

% define the orientation of the vehicle relative to the EKF NED frame
rpy_ekf = [0;0;-pi/4];
quat_ekf = EulToQuat(rpy_ekf);
quat_ekf = NormQuat(quat_ekf);

% define the orientation of the vehicle relative to the EV frame
rpy_ev = [0;0;-pi/2];
quat_ev = EulToQuat(rpy_ev);
quat_ev = NormQuat(quat_ev);

% calculate the delta quaternion
quat_ev_inv = [quat_ev(1);-quat_ev(2);-quat_ev(3);-quat_ev(4)];
quat_delta = QuatMult(quat_ekf,quat_ev_inv);
quat_delta = NormQuat(quat_delta);

% convert to a delta angle so it can be filtered
ang_delta_filt = QuatToDeltaAngle(quat_delta);

% convert back to quaternion
quat_delta_filt = RotToQuat(ang_delta_filt);
quat_delta_filt = NormQuat(quat_delta_filt);

% convert to a rotation matrix to rotate EV measurements to EKF frame
Tev_ekf = Quat2Tbn(quat_delta_filt);

