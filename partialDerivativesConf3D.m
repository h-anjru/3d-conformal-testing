function A = partialDerivativesConf3D(scale, rotm, vec)
    arguments
        scale (1,1) double {mustBePositive}
        rotm (3,3) double {mustBeNumeric}
        vec (3,1) double {mustBeNumeric}
    end

    opk = opkFromRotationMatrix(rotm);

    a11 = rotm(:, 1)' * vec;
    a12 = 0;
    a13 = [-sin(opk(2)) * cos(opk(3)), ...
        sin(opk(2)) * sin(opk(3)), ...
        cos(opk(2))] * vec * scale;
    a14 = (rotm(2, 1) * vec(1) - rotm(1, 1) * vec(2)) * scale;

    a21 = rotm(:, 2)' * vec;
    a22 = -rotm(:, 3)' * vec * scale;
    a23 = sin(opk(1)) * [1 -1 1] .* rotm(:, 1)' * vec * scale;
    a24 = (rotm(2, 2) * vec(1) - rotm(1, 2) * vec(2)) * scale;

    a31 = rotm(:, 3)' * vec;
    a32 = (rotm(:,2)' * vec) * scale;
    a33 = -cos(opk(1)) * rotm(:, 1)' * vec * scale;
    a34 = (rotm(2, 3) * vec(1) - rotm(1, 3) * vec(2)) * scale;

    A = [
        a11 a12 a13 a14 1 0 0
        a21 a22 a23 a24 0 1 0
        a31 a32 a33 a34 0 0 1
    ];

end
