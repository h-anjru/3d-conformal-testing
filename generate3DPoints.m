function [arbitrary, control, hgt, noise] = generate3DPoints( ...
    n, noise_scalar, s_rng, o_rng, p_rng, k_rng, TX_rng, TY_rng, TZ_rng)
%% GENERATE3DPOINTS  Generate a set of 3D points in two coordinate systems.
%
%   This function generates a set of points in an arbitrary coordinate 
%   system and applies to those points a transformation of scale, rotation,
%   and translation to determine those points' coordinates in a control
%   coordinate system. In addition to the aribitrary coordinates, the
%   control coordinates, and the transformation matrix, the function can
%   also generate a "noise" matrix, which can be added to either of the
%   point datasets to simulate real-world conditions. It is recommended to
%   add this noise to the arbitrary dataset.
%
%   The only required arugment is number of points n. All other
%   transformation parameters have default values. 
%
%   The default noise scalar is set to 0.
%
%   INPUT:
%                  n    number of points to be generated
%              s_rng    [1, 2] range for random scale for transformation
%              o_rng    [1, 2] range for random omega rotation [deg]
%              p_rng    [1, 2] range for random phi rotation [deg]
%              k_rng    [1, 2] range for random kappa rotation [deg]
%             TX_rng    [1, 2] range for random X translation [control]
%             TY_rng    [1, 2] range for random Y translation [control]
%             TZ_rng    [1, 2] range for random Z translation [control]
%       noise_scalar    (optional) scalaer for random noise
%
%   OUTPUT:
%          arbitrary    [3, n] aribitrary point coordinates as vectors
%            control    [3, n] control point coordinates as vectors
%                hgt    [4, 4] homogeneous transformation vector
%              noise    [3, n] noise component

    arguments
        n (1,1) {mustBeInteger, mustBePositive}
        noise_scalar (1,1) {mustBePositive} = 0
        s_rng (1,2) {mustBePositive} = [0.8 1.2]
        o_rng (1,2) {mustBeNumeric} = [-180 180]
        p_rng (1,2) {mustBeNumeric} = [-180 180]
        k_rng (1,2) {mustBeNumeric} = [-180 180]
        TX_rng (1,2) {mustBeNumeric} = [-1000 1000]
        TY_rng (1,2) {mustBeNumeric} = [-1000 1000]
        TZ_rng (1,2) {mustBeNumeric} = [-1000 1000]
    end
    
    % generate points in an arbitrary system
    arbitrary = randRange([-100 100], [3 n]);
    
    % generate true transformation parameters
    scale = randRange(s_rng);
    omega = randRange(rad2deg(o_rng));
    phi = randRange(rad2deg(p_rng));
    kappa = randRange(rad2deg(k_rng));
    TX = randRange(TX_rng);
    TY = randRange(TY_rng);
    TZ = randRange(TZ_rng);
    
    % form the transform matrix
    % Note: makehgtform() generates an active rotation (4 x 4 homogeneous)
    hgt = makehgtform('xrotate', omega, 'yrotate', phi, 'zrotate', kappa);
    hgt(1:3, 1:3) = scale * hgt(1:3, 1:3);
    hgt(1:3, 4) = [TX, TY, TZ];
    
    % transform arbitrary to control
    arbitrary_hmg = [arbitrary; ones(1, n)];
    control_hmg = hgt * arbitrary_hmg;
    control = control_hmg(1:3, :);  % no need to divide by w (w = 1)
    
    % generate noise component (recommended to add to arbitrary points)
    noise = noise_scalar * randn(size(arbitrary));
    % arbitrary_noise = arbitrary + noise;
end

function r = randRange(rng, size)
    arguments
        rng (1,2) double {mustBeNumeric}
        size (1,2) double {mustBeNumeric} = [1 1]
    end

    r = rng(1) + (rng(2) - rng(1)) * rand(size);
end
