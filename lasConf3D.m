function [hgt, jac, Kvec] = lasConf3D(arbitrary, control)
%% LAS3DCONF Lassiter 3D conformal transformation method
%
%   This function performs the 3D conformal coordinate transformation
%   method developed by Lassiter. Similar to the method proposed by Wolf et
%   al. in Elements of Photogrammetry, the algorithm uses an iterative,
%   nonlinear least squares adjustment, but only to solve for the rotation
%   between the two coordinate systems. Scale and translation are direct
%   solutions.
%
%   This method is novel in that the nonlinear least squares adjustment is
%   performed only to find the rotation. First, both sets of coordinates
%   are translated to be centered about their respective origins. Scale is 
%   then estimated as the mean of the ratios of each control point vector 
%   norm and its corresponding arbitrary vector norm. A direct linear 
%   estimation of the elements of the rotation matrix is performed, and 
%   from this matrix estimations of omega, phi, and kappa are derived. 
%   After the nonlinear least squares adjustment to find optimal values of
%   omega, phi, and kappa, the translation is found as the difference 
%   between the centroid of the control coordinates and the scaled and
%   rotated centroid of the arbitrary coordinates.
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

    % Step 1. Find centroids and subtract from respective sets of points
    centroid_arb = mean(arbitrary, 2);
    centroid_con = mean(control, 2);

    arb = arbitrary - repmat(centroid_arb, [1, n]);
    con = control - repmat(centroid_con, [1, n]);

    % Step 2. Solve for scale
    scale = mean(vecnorm(con, 2) ./ vecnorm(arb, 2));

    % Step 3. Scale arbitrary coords
    arb_scaled = scale * arb;

    % Step 4. Direct linear estimation of rotation
    coeff = zeros(n * 3, 9);

    for ii = 1:n
        coeff(3 * ii - 2, 1:3) = arb_scaled(:, ii)';
        coeff(3 * ii - 1, 4:6) = arb_scaled(:, ii)';
        coeff(3 * ii - 0, 7:9) = arb_scaled(:, ii)';
    end

    Lvec = reshape(con, [n * 3, 1]);

    rotmvec_init = coeff \ Lvec;
    rotm_init = reshape(rotmvec_init, [3, 3]);

    % To keep the bottom left element <= 1, the initial rotation matrix is
    % divided by the norm of the first column vector. (The estimation of
    % the phi rotation angle is the arcsin of the bottom left element, and
    % must be <= 1.)
    rotm_norm = vecnorm(rotm_init);

    rotm = rotm_init / rotm_norm(1);

    % Step 5. Nonlinear LS solution of rotation
    % Step 5a. Initial apprx of omega, phi, kappa
    opk = opkFromRotationMatrix(rotm);

    % Step 5b. Set up iteration
    L1 = 1;  % arbitrary initla value of L1 norm of NLS solution
    tol = 1e-4;
    iter = 1;
    iter_max = 20;

    % initialize Jacobian matrix
    jac = zeros(3 * n, 3);

    while L1 > tol && iter <= iter_max
        % Step 5c. Jacobian matrix   
        for ii = 1:n
            % note: on first pass, this rotm is non-orthonormal
            A = partialDerivativesConf3D(1, rotm, arb_scaled(:, ii));
            jac(3 * ii - 2:3 * ii, :) = A(:, 2:4);
        end

        % Step 5d. "Observed minus estimated" vector
        con_est = rotm(1:3, 1:3)' * arb_scaled;
        Kvec = reshape(con - con_est, [n * 3, 1]);

        % Step 5e. Solve and update
        delta = jac \ Kvec;

        opk = opk + delta;

        % rotation matrix from current approximations of o, p, k  
        rotm = makehgtform('xrotate', opk(1), ...
            'yrotate', opk(2), 'zrotate', opk(3));
        rotm = rotm(1:3, 1:3)';

        L1 = norm(delta, 1);
        iter = iter + 1;
    end

    % Step 6. Find translation
    scale_times_rotm = scale * rotm(1:3, 1:3);
    translation = centroid_con - scale_times_rotm' * centroid_arb;

    % Step 7. Form transformation matrix
    hgt = [scale_times_rotm' translation; 0 0 0 1];

end
