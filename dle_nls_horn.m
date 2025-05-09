% This script performs a Monte Carlo simulation which perfrorms the 3D 
% conformal coordinate transformation via (1) nonlinear least squares 
% with initial approximations performed by the "direct linear 
% estimation" method and (2) the Horn closed-form quaternoin method.
% 
% Each transformation is performed on the same set of coordinate pairs.


% user input
trials = 10000;
common = 10;
check = 20;
noisyness = 0.05;

total = common + check;

% initialize flags for special cases
gimbal_flag = false;
flags = zeros(trials, 1);

% initialize result arrays
meanNorms_nls = zeros(trials, 1);
meanNorms_horn = zeros(trials, 1);

all_arb_noised = zeros([3, total, trials]);
all_con = zeros([3, total, trials]);

all_hgt_true = zeros([4, 4, trials]);
all_hgt_init = zeros([4, 4, trials]);
all_hgt_nls = zeros([4, 4, trials]);
all_hgt_horn = zeros([4, 4, trials]);

for ii = 1:trials
    % generate points
    [arb, con, hgt_true, noise] = generate3DPoints(total, noisyness);
    arb_noised = arb + noise;

    all_arb_noised(:, :, ii) = arb_noised;
    all_con(:, :, ii) = con;
    
    all_hgt_true(:, :, ii) = hgt_true;
    
    %% DLT approximation method
    hgt_init = conf3d_dle(arb_noised(:, 1:common), con(:, 1:common));

    [hgt_nls, jac_nls, Kvec_nls, hgt_init, gimbal_flag] = ...
        conf3d_nls(arb_noised(:, 1:common), con(:, 1:common), hgt_init);

    all_hgt_init(:, :, ii) = hgt_init;
    all_hgt_nls(:, :, ii) = hgt_nls;

    % DLT transform checkpoints
    chk_nls = hgt_nls * [arb_noised(:, common+1:total); ones(1, check)];
    
    % DLT get norms
    norms_nls = vecnorm(chk_nls(1:3, :) - con(:, common+1:total));
    meanNorms_nls(ii) = mean(norms_nls);

    % flag gimbal lock and reset
    flags(ii) = gimbal_flag;    
    
    %% Horn method
    [hgt_horn, M, N, V, D] = ...
        conf3d_horn(arb_noised(:, 1:common), con(:, 1:common));

    all_hgt_horn(:, :, ii) = hgt_horn;
    
    % Horn transform checkpoints
    chk_horn = hgt_horn * [arb_noised(:, common+1:total); ones(1, check)];
    
    % Horn get norms
    norms_horn = vecnorm(chk_horn(1:3, :) - con(:, common+1:total));
    meanNorms_horn(ii) = mean(norms_horn);

    % % break it off to explore
    % if mean(norms_nls) > 0.5 || gimbal_flag
    %     break
    % end

end

%% comparison stats
% performance ratio
performance_ratio = mean(meanNorms_nls) / mean(meanNorms_horn);

% check DLT init. aprx. quality (used later for plots)
[all_rot_true, all_scale_true] = params_from_hgt(all_hgt_true);
[all_rot_init, all_scale_init] = params_from_hgt(all_hgt_init);
[all_rot_nls, all_scale_nls] = params_from_hgt(all_hgt_nls);
[all_rot_horn, all_scale_horn] = params_from_hgt(all_hgt_horn);

% compare DLT init. aprx. to true
scale_check = all_scale_init - all_scale_true;
rot_check = all_rot_init - all_rot_true;
trans_check = all_hgt_init(1:3, 4, :) - all_hgt_true(1:3, 4, :);

%% time to plot all this nonsense!
bin_width = 1 / trials;  % seems to work well for larger trials

fig1 = figure(1);
hist_nls = histogram(meanNorms_nls);
ax1 = gca;
title(ax1, 'NLS mean of norms');

fig2 = figure(2);
hist_nls_y = histogram(meanNorms_nls);
ax2 = gca;
ax2.YLim = [0 log10(trials)];
title(ax2, 'NLS mean of norms (truncated y-axis)');

fig3 = figure(3);
hist_nls_x = histogram(meanNorms_nls);
ax3 = gca;
title(ax3, 'NLS mean of norms (truncated x-axis)');

fig4 = figure(4);
hist_horn = histogram(meanNorms_horn);
ax4 = gca;
title(ax4, 'Horn mean of norms');

fig5 = figure(5);
hist_scale = histogram(scale_check);
ax5 = gca;
title(ax5, 'DLT init. approx. scale discrepencies');

fig6 = figure(6);
hist_rot = histogram(rot_check);
ax6 = gca;
title(ax6, 'DLT init. approx. rotation parameter discrepencies');

fig7 = figure(7); 
hist_trans = histogram(trans_check);
ax7 = gca;
title(ax7, 'DLT init. approx. translation parameter discrepencies');

axes = [ax1 ax2 ax3 ax4 ax5 ax6 ax7];
hists = [hist_nls hist_nls_y hist_nls_x hist_horn];  % excluding DLTs

% bin width, figures 1-4
for ii = 1:length(hists)
    hists(ii).BinWidth = bin_width;
end

% text size, all figures
for ii = 1:length(axes)
    axes(ii).FontSize = 14;
end

% allow visual comparison between NLS and Horn
ax3.XLim = ax4.XLim;
ax3.YLim = ax4.YLim;
ax1.YLim = ax4.YLim;

function [rot, scale] = params_from_hgt(hgt)
% PARAMS_FROM_HGT Isolate transformation parameters from a homogeneous
% transformation matrix. Can accept 3D arrays (i.e., pages).

    rot = hgt(1:3, 1:3, :);
    rot = pagetranspose(rot);
    scale = vecnorm(rot);
    scale = scale(:, 1, :);
    rot = rot ./ scale;
end
