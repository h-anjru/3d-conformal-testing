%% Lassiter direct linear initial approximation method

function [hgt, jac, Kvec] = lasConf3D_2(arbitrary, control)
%% LAS3DCONF2 (Another) Lassiter 3D conformal transformation method
%
%   This function performs a 3D conformal coordinate transformation
%   method developed by Lassiter. Similar to the method proposed by Wolf et
%   al. in Elements of Photogrammetry, the algorithm uses an iterative,
%   nonlinear least squares adjustment, but only to solve for the rotation
%   between the two coordinate systems. Scale and translation are direct
%   solutions.
%   
%   This method performs a direct linear estimation of the elements of the
%   homogeneous transformation matrix, i.e., the elements of the rotation
%   matrix multipled by the scale, and the translations in X, Y, and Z.
%   From these estimations, an estimation of scale is found to be the norm
%   of the first column vector of the scaled rotation matrix, and the
%   estimations of omega, phi, and kappa are then found in their usual
%   manner. As this direct linear method is solving for 12 unknowns, a
%   minimum of four points must be used.
%
%   INPUT:
%       arbitrary   [3, n] matrix of column vectors (coordinates)
%         control   [3, n] matrix of column vectors (coordinates)
%    
%   OUTPUT:
%             hgt   [4, 4] homogeneous transform matrix
%             jac   [n, 3] Jacobian matrix of final iteration
%            Kvec   [n, 1] observation vector of final iteration

    arguments
        arbitrary (3,:) {mustBeNumeric}
        control (3,:) {mustBeNumeric}
    end

    n = size(arbitrary, 2);

    % Step 1. Direct linear estimation of all parameters (n >= 4)
    coeff = zeros(n * 3, 9);
    
    for ii = 1:n
        coeff(3 * ii - 2, 1:3) = arbitrary(:, ii)';
        coeff(3 * ii - 1, 4:6) = arbitrary(:, ii)';
        coeff(3 * ii - 0, 7:9) = arbitrary(:, ii)';
    end
    
    coeff = [coeff repmat(eye(3), [n 1])];
    
    Lvec = reshape(control, [n * 3, 1]);
    
    params_init = coeff \ Lvec;
    rotm_init = reshape(params_init(1:9), [3, 3]);
    
    % Step 1b. Scale estimation
    scale = norm(rotm_init(:, 1));
    
    rotm = rotm_init / scale;
    
    % Step 1c. Translation estimation
    translation = params_init(10:12);
    
    % Step 2. Nonlinear LS solution of rotation
    % Step 2a. Initial apprx of omega, phi, kappa
    opk = opkFromRotationMatrix(rotm);
    
    % Step 2b. Set up iteration
    L1 = 1;  % arbitrary initial value of L1 norm of NLS solution
    tol = 1e-4;
    iter = 1;
    iter_max = 20;

    % initialize Jacobian matrix
    jac = zeros(3 * n, 7);
    
    while L1 > tol && iter <= iter_max
        % Step 2c. Jacobian matrix  
        for ii = 1:n
            % note: on first pass, this rotm is non-orthonormal
            A = partialDerivativesConf3D( ...
                scale, rotm, arbitrary(:, ii) ...
            );
            jac(3 * ii - 2:3 * ii, :) = A;
        end
        
        % Step 2d. "Observed minus estimated" vector
        las2_hgt_aprx = [scale * rotm(1:3, 1:3)' translation; 0 0 0 1];
        con_est = las2_hgt_aprx * [arbitrary; ones(1, n)];
        Kvec = reshape(control - con_est(1:3, :), [3 * n, 1]);
        
        % Step 5e. Solve and update
        delta = jac \ Kvec;
    
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
    
    % Step 3: Form transformation matrix
    hgt = [scale * rotm(1:3, 1:3)' translation; 0 0 0 1];

end
