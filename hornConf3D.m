function [hgt, M, N, V, D] = hornConf3D(arbitrary, control)
%% HORN3DCONF  Horn absolute orientation (quaternion method)
%
%   This function performs BKP Horn's algorithm for a closed-form solution 
%   to the absolute orientation problem (i.e., a 3D conformal 
%   transformation). The output is a homogeneous transformation matrix
%   which describes the transformation from the aribtary to the control
%   coordinate system.
%
%   INPUT:
%       arbitrary   [3, n] matrix of column vectors (coordinates)
%         control   [3, n] matrix of column vectors (coordinates)
% 
%   OUTPUT:
%             hgt   [4, 4] homogeneous transform matrix
%
%   REFERENCE: 
%   Berthold K. P. Horn, "Closed-form solution of absolute orientation 
%   using unit quaternions," J. Opt. Soc. Am. A 4, 629-642 (1987)
%   doi:10.1364/JOSAA.4.000629

    arguments
        arbitrary (3,:) {mustBeNumeric}
        control (3,:) {mustBeNumeric}
    end

    % Step 1. Find centroids and subtract from respective sets of points
    centroid_arb = mean(arbitrary, 2);
    centroid_con = mean(control, 2);

    n = size(arbitrary, 2);

    arbitary_translated = arbitrary - repmat(centroid_arb, [1, n]);
    control_translated = control - repmat(centroid_con, [1, n]);

    % Step 2. Find matrix M, the sums of squares of corresponding points
    M = arbitary_translated * control_translated';

    % Step 3. Find matrix N (subsection 4.A)
    t = num2cell(reshape(M, [9, 1]));
    [xx, yx, zx, xy, yy, zy, xz, yz, zz] = deal(t{:});

    N = [
        xx + yy + zz  yz - zy       zx - xz       xy - yx
        yz - zy       xx - yy - zz  xy + yx       zx + xz
        zx - xz       xy + yx      -xx + yy - zz  yz + zy
        xy - yx       zx + xz       yz + zy      -xx - yy + zz
    ];

    % Step 4. Find eigenvector whose eigenvalue has most positive value
    [V, D] = eig(N);

    [~, idx] = max(diag(D));  % index of most positive eigenvalue
    unit_quaternion = V(:, idx);

    horn_rotm = quat2rotm(unit_quaternion');

    % Step 5. Find scale
    horn_scale = sum(vecnorm(control_translated)) / ...
        sum(vecnorm(arbitary_translated));

    % Step 6. Find translation
    scale_times_rotm =  horn_scale * horn_rotm;
    translation = centroid_con - scale_times_rotm * centroid_arb;

    % Step 7. Form transformation matrix
    hgt = [scale_times_rotm translation; 0 0 0 1];

end
