function jac = conf3d_partialDerivatives(scale, rotm, vec)
% CONF3D_PARTIALDERIVAIVES Generate Jacobian matrix for 3d conf. NLS
%   Reference: Wolf et al. Elements of Photogrammetry, 4th ed.

    arguments
        scale (1,1) double {mustBePositive}
        rotm (3,3) double {mustBeNumeric}
        vec (3,1) double {mustBeNumeric}
    end

    opk = opkFromRotationMatrix(rotm);

    j11 = rotm(:, 1)' * vec;
    j12 = 0;
    j13 = [-sin(opk(2)) * cos(opk(3)), ...
        sin(opk(2)) * sin(opk(3)), ...
        cos(opk(2))] * vec * scale;
    j14 = (rotm(2, 1) * vec(1) - rotm(1, 1) * vec(2)) * scale;

    j21 = rotm(:, 2)' * vec;
    j22 = -rotm(:, 3)' * vec * scale;
    j23 = sin(opk(1)) * [1 -1 1] .* rotm(:, 1)' * vec * scale;
    j24 = (rotm(2, 2) * vec(1) - rotm(1, 2) * vec(2)) * scale;

    j31 = rotm(:, 3)' * vec;
    j32 = (rotm(:,2)' * vec) * scale;
    j33 = -cos(opk(1)) * rotm(:, 1)' * vec * scale;
    j34 = (rotm(2, 3) * vec(1) - rotm(1, 3) * vec(2)) * scale;

    jac = [
        j11 j12 j13 j14 1 0 0
        j21 j22 j23 j24 0 1 0
        j31 j32 j33 j34 0 0 1
    ];

end
