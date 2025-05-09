function [hgt, jac, obs_vec, hgt_init, flag] = ...
    conf3d_nls(arbitrary, control, hgt_init, gimbal_tol)
%% CONF3D_NLS (Another) Nonlineasr LS solution for 3d conformal trans.
%
%   This function performs a 3D conformal coordinate transformation via
%   nonlinear least squares as presented in Wolf et al. and other texts on
%   photogrammetry. The initial approximations for the transformation are
%   input as a homogeneous transformation matrix. Gimbal tolerance is an
%   interval +/- 1 around the 3,1 element of the rotation matrix inside of
%   which X-Z gimbal lock is flagged.
%
%   INPUT:
%       arbitrary   [3, n] matrix of column vectors (coordinates)
%         control   [3, n] matrix of column vectors (coordinates)
%        hgt_init   [4, 4] homogeneous transformation matrix
%      gimbal_tol   tolerance +/- 1 for flagging gimbal lock
%    
%   OUTPUT:
%             hgt   [4, 4] homogeneous transform matrix
%             jac   [n, 3] Jacobian matrix of final iteration
%            Kvec   [n, 1] observation vector of final iteration
%            flag   bool indicating detection of gimbal lock

    arguments
        arbitrary (3,:) {mustBeNumeric}
        control (3,:) {mustBeNumeric}
        hgt_init (4,4) {mustBeNumeric}
        gimbal_tol (1,1) {mustBeNumeric} = 1e-4
    end

    n = size(arbitrary, 2);

    % Step 0. Pull 7 parameters out of initial HGT
    rotm_init = hgt_init(1:3, 1:3)';
    scale = norm(rotm_init(:, 1));
    rotm = rotm_init / scale;
    translation = hgt_init(1:3, 4);

    % Step 1. Nonlinear LS solution of rotation
    % Step 1a. Initial apprx of omega, phi, kappa
    [opk, gimbal_flag] = opkFromRotationMatrix(rotm, gimbal_tol);

    % Step 1b. Set up iteration
    L1 = 1;  % arbitrary initial value of L1 norm of NLS solution
    tol = 1e-10;
    iter = 1;
    iter_max = 20;

    % initialize Jacobian matrix
    jac = zeros(3 * n, 7);

    while L1 > tol && iter <= iter_max
        % Step 1c. Jacobian matrix  
        for ii = 1:n
            % note: on first pass, this rotm is non-orthonormal
            A = conf3d_partialDerivatives(scale, rotm, arbitrary(:, ii));
            jac(3 * ii - 2:3 * ii, :) = A;
        end

        % Step 1d. "Observed minus estimated" vector
        hgt_aprx = [scale * rotm(1:3, 1:3)' translation; 0 0 0 1];
        control_est = hgt_aprx * [arbitrary; ones(1, n)];
        obs_vec = reshape(control - control_est(1:3, :), [3 * n, 1]);

        % Step 1e. Solve and update
        delta = jac \ obs_vec;
       
        scale = scale + delta(1);
        opk = opk + delta(2:4);
        translation = translation + delta(5:7);

        % rotation matrix from current approximations of o, p, k  
        rotm = makehgtform('xrotate', opk(1), 'yrotate', opk(2), ...
            'zrotate', opk(3));
        rotm = rotm(1:3, 1:3)';

        L1 = norm(delta, 1);
        iter = iter + 1;
    end

    % Step 3. Form transformation matrix
    hgt = [scale * rotm(1:3, 1:3)' translation; 0 0 0 1];

    % Gimbal lock flag
    flag = false;
    
    if gimbal_flag
        flag = true;
    end

end
