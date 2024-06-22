function [opk, gimbal_flag] = opkFromRotationMatrix(rotm, gimbal_tol)
%% OPKFROMROTATIONMATRIX  Derive omega, phi, kappa from a rotation matrix
%
%   Find the omega, phi, and kappa rotation angles from a passive rotation 
%   matrix. The fucntion will also check for potential gimbal lock, which 
%   occurs when the absolute values of the phi rotation angle equals 90°. 
%   As the angle phi is determined by taking the arcsin of the bottom-left 
%   element of the passive rotation matrix, the gimbal lock test determines
%   whether that element is near 1 (within a default tolerance of 1e-4).
%
%   If potential gimbal lock is detected, the omega rotation angle is set
%   to zero, and all remaining rotation after the phi rotation is expressed
%   as kappa.
%
%   INPUT:
%             rotm  [3, 3] passive rotation matrix
%       gimbal_tol  +/- range about 1 that flags rotation as gimbal locked
%
%   OUTPUT:
%              opk  [3, 1] vector of omega, phi, kappa [radians]
%      gimbal_flag  logical indicating whether gimbal lock detected

    arguments
        rotm (3,3) {mustBeNumeric}
        gimbal_tol (1,1) {mustBeNumeric} = 1e-4
    end

    % test for gimbal lock
    gimbal_flag = abs(rotm(3, 1)) > 1 - gimbal_tol && ...
        abs(rotm(3, 1)) < 1 + gimbal_tol;

    if gimbal_flag
        o = 0;
        k = atan2(-rotm(1, 2), rotm(2, 2));
    else
        o = atan2(-rotm(3, 2), rotm(3, 3));
        k = atan2(-rotm(2, 1), rotm(1, 1));
    end
    
    p = asin(rotm(3, 1));

    opk = [o p k]';

end
