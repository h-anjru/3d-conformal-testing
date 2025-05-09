function hgt_init = conf3d_dle(arbitrary, control)
%% CONF3D_DLE Direct linear estimation of 3D conformal trans. parameters
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
%        hgt_init   [4, 4] homogeneous transform matrix
%            flag   bool indicating if gimbal lock is detected

    arguments
        arbitrary (3,:) {mustBeNumeric}
        control (3,:) {mustBeNumeric}
    end

    n = size(arbitrary, 2);

    % Step 1. Linear least squares solution for 12 HGT elements
    coeff_matrix = zeros(3*n, 9);
    
    for ii = 1:n
        coeff_matrix(3 * ii - 2, 1:3) = arbitrary(:, ii)';
        coeff_matrix(3 * ii - 1, 4:6) = arbitrary(:, ii)';
        coeff_matrix(3 * ii - 0, 7:9) = arbitrary(:, ii)';
    end

    coeff_matrix = [coeff_matrix repmat(eye(3), [n 1])];

    obs_vec = reshape(control, [3*n, 1]);

    params_init = coeff_matrix \ obs_vec;
    rotm_init = reshape(params_init(1:9), [3, 3]);

    % Step 2. Scale estimation
    scale = norm(rotm_init(:, 1));

    rotm = rotm_init / scale;

    % Step 3. Translation estimation
    translation = params_init(10:12);

    % Step 4. Form initial homogeneous transform matrix
    hgt_init = [scale * rotm(1:3, 1:3)' translation; 0 0 0 1];

end
