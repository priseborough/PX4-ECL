function [ sigma_x_a ] = CalcSigmaPoints(states, P , sQ , param)

    % Calculate the lower diagonal Cholesky decomposition for the vehicle
    % state covariance matrix. This requires nP^3 operations
    sP = chol(P,'lower');
    
    % Assemble the augmented covariance matrix
    sPA = [sP,zeros(param.ukf.nP,param.ukf.nQ);zeros(param.ukf.nQ,param.ukf.nP),sQ];
    
    % Expected value of augmented state vector from previous frame
    x_a_prev = [...
        states;... % vehicle states
        zeros(param.ukf.nQ,1)... % IMU noise is zero mean
        ];
    
    % Generate sigma points for the augmented state vector
    % Note : For a real time implementation we should try to take advantage
    % of sparseness in sPA and the fact that sqrt(L+lambda) is constant
    sigma_x_a = [x_a_prev, x_a_prev*ones(1,param.ukf.L)+sqrt(param.ukf.L+param.ukf.lambda)*sPA, x_a_prev*ones(1,param.ukf.L)-sqrt(param.ukf.L+param.ukf.lambda)*sPA];


end

