function ph = fn_stimulusPhaseJPC(spike_conv, fs, dyn, t_range_s, stim_t, stim_names, figDir, recID)
% FN_STIMULUSPHASEJPC  Phase of each stimulus as a ROTATION ANGLE in the jPC plane.
%
%   ph = FN_STIMULUSPHASEJPC(spike_conv, fs, dyn, t_range_s, stim_t, stim_names, figDir, recID)
%
% fn_dynamicalCycleMetrics reports phase as the fraction of cycle TIME elapsed, and cycles are cut
% where the trajectory returns, so phase 0 is not a fixed point on the attractor. Here the phase is
% the angle of the population state in the dominant rotational (jPC) plane, so the same angle means
% the same point on the cycle, within one recording.
%
% METHOD (self-contained re-fit; it does not call the jPCA_Aplysia code)
%   1. Bin spike_conv exactly as fn_dynamicalEpochs did (same bin, same window), subtract the mean.
%   2. PCA to 6 PCs (same dimension as the rot6 control) -> scores Z.
%   3. jPCA model: dZ/dt = M Z with M skew-symmetric, fitted by least squares on all 15 free entries.
%   4. Real Schur form of M gives three orthonormal rotational planes; the plane capturing most
%      variance is kept ("first jPC plane"). Axes are signed so rotation is counter-clockwise.
%   5. theta(t) = unwrap(atan2(z2, z1)). Per cycle: rotation swept (degrees) between its start and
%      end bin. Per stimulus: absolute angle, angle-based fraction of the containing cycle, and the
%      time-based fraction for comparison.
%
% OUTPUT ph (struct): ok, plane_var_frac (variance in the plane / total variance), rot6_pipeline
%   (dyn.rot6, for comparison), t (s), theta (unwrapped, rad), z (nb x 2), cycle_sweep_deg (nCyc x 1),
%   stim(k): name, t, cycle, phase_abs_rad, phase_abs_deg, phase_cycle_frac_angle,
%            phase_cycle_frac_time, cycle_sweep_deg
%
% NOTES
%   * Absolute angle is only comparable between stimuli within one recording (the plane's axes are
%     arbitrary). The jPC plane is fitted on the whole analysis window, so it is dominated by the
%     epoch with the most variance.
%   * The jPCA code in the toolbox also normalises the data; this re-fit does not, so plane_var_frac
%     is close to, not identical with, dyn.rot6. A large mismatch means "check", not "wrong".
%   * Cycle sweeps far from 360 degrees mean that cycle is not a clean rotation in this plane.

ph = struct('ok', false);
if nargin < 7, figDir = ''; end
if nargin < 8, recID = ''; end
if ~isfield(dyn, 'ok') || ~dyn.ok
    fprintf('  Stimulus phase (jPC plane): skipped (no cycles).\n'); return
end

%% 1. Bin as in fn_dynamicalEpochs
[nN, nFrames] = size(spike_conv);
bf      = max(round(fs * dyn.params.bin_ms / 1000), 1);
bin_eff = 1000 * bf / fs;
f1 = max(round(t_range_s(1) * fs) + 1, 1);
f2 = min(round(t_range_s(2) * fs), nFrames);
nb = floor((f2 - f1 + 1) / bf);
t0 = (f1 - 1) / fs;
X = zeros(nb, nN);
for b0 = 1:500:nb
    b1  = min(b0 + 499, nb);
    seg = double(spike_conv(:, f1 + (b0-1)*bf : f1 + b1*bf - 1));
    X(b0:b1, :) = reshape(mean(reshape(seg, nN, bf, b1-b0+1), 2), nN, b1-b0+1)';
end
X  = X - mean(X, 1);
dt = bin_eff / 1000;
tvec = t0 + ((1:nb)' - 0.5) * dt;

%% 2-3. 6-PC scores and skew-symmetric linear fit
K = 6;
[~, Z] = pca(X, 'NumComponents', K);
dZ = diff(Z) / dt;
Zm = Z(1:end-1, :);
pairs = nchoosek(1:K, 2);
A = zeros(numel(dZ), size(pairs, 1));
for p = 1:size(pairs, 1)
    i = pairs(p, 1); j = pairs(p, 2);
    Q = zeros(size(Zm));
    Q(:, i) =  Zm(:, j);
    Q(:, j) = -Zm(:, i);
    A(:, p) = Q(:);
end
a = A \ dZ(:);
M = zeros(K);
for p = 1:size(pairs, 1)
    M(pairs(p,1), pairs(p,2)) =  a(p);
    M(pairs(p,2), pairs(p,1)) = -a(p);
end

%% 4. Choose the rotational plane with the most variance
[U, T] = schur(M, 'real');
planes = {};
c = 1;
while c < K
    if abs(T(c+1, c)) > 1e-12 * max(1, norm(T, 'fro'))
        planes{end+1} = [c, c+1]; %#ok<AGROW>
        c = c + 2;
    else
        c = c + 1;
    end
end
if isempty(planes)
    fprintf('  Stimulus phase (jPC plane): skipped (no rotational plane found).\n'); return
end
bestVar = -inf; bestZ = [];
for q = 1:numel(planes)
    zq = Z * U(:, planes{q});
    v  = sum(var(zq, 0, 1));
    if v > bestVar, bestVar = v; bestZ = zq; end
end
z = bestZ;
theta = unwrap(atan2(z(:,2), z(:,1)));
if mean(diff(theta)) < 0                % make rotation counter-clockwise
    z(:,2) = -z(:,2);
    theta  = unwrap(atan2(z(:,2), z(:,1)));
end

%% 5. Per-cycle sweep and per-stimulus phase
sb   = dyn.startbin(:);
nCyc = size(dyn.cycle_t, 1);
sweep = nan(nCyc, 1);
for k = 1:nCyc
    b1 = min(max(sb(k), 1), nb); b2 = min(max(sb(k+1), 1), nb);
    sweep(k) = rad2deg(theta(b2) - theta(b1));
end

stim = struct('name', {}, 't', {}, 'cycle', {}, 'phase_abs_rad', {}, 'phase_abs_deg', {}, ...
    'phase_cycle_frac_angle', {}, 'phase_cycle_frac_time', {}, 'cycle_sweep_deg', {});
for k = 1:numel(stim_t)
    t = stim_t(k);
    s.name = stim_names{k}; s.t = t; s.cycle = NaN; s.phase_abs_rad = NaN; s.phase_abs_deg = NaN;
    s.phase_cycle_frac_angle = NaN; s.phase_cycle_frac_time = NaN; s.cycle_sweep_deg = NaN;
    b  = min(max(floor((t - t0) / dt) + 1, 1), nb);
    cc = find(dyn.cycle_t(:,1) <= t & t < dyn.cycle_t(:,2), 1, 'first');
    s.phase_abs_rad = mod(theta(b), 2*pi);
    s.phase_abs_deg = rad2deg(s.phase_abs_rad);
    if ~isempty(cc)
        b1 = min(max(sb(cc), 1), nb); b2 = min(max(sb(cc+1), 1), nb);
        s.cycle = cc;
        s.cycle_sweep_deg = sweep(cc);
        if theta(b2) ~= theta(b1)
            s.phase_cycle_frac_angle = (theta(b) - theta(b1)) / (theta(b2) - theta(b1));
        end
        s.phase_cycle_frac_time = (t - dyn.cycle_t(cc,1)) / (dyn.cycle_t(cc,2) - dyn.cycle_t(cc,1));
    end
    stim(end+1) = s; %#ok<AGROW>
end

ph.ok = true;
ph.plane_var_frac = bestVar / sum(var(X, 0, 1));
ph.rot6_pipeline  = dyn.rot6;
ph.t = tvec; ph.theta = theta; ph.z = z;
ph.cycle_sweep_deg = sweep; ph.stim = stim;

fprintf('  Stimulus phase (jPC plane): plane variance %.1f%% (pipeline rot6 %.1f%%) | median cycle sweep %.0f deg (IQR %.0f-%.0f)\n', ...
    100*ph.plane_var_frac, 100*dyn.rot6, median(sweep, 'omitnan'), prctile(sweep, 25), prctile(sweep, 75));
for k = 1:numel(stim)
    s = stim(k);
    if isnan(s.cycle)
        fprintf('    %s (%.1f s): angle %.0f deg, not inside any detected cycle\n', s.name, s.t, s.phase_abs_deg);
    else
        fprintf('    %s (%.1f s): cycle %d | angle in plane %.0f deg | %.2f of the rotation (%.2f of the time) | cycle sweeps %.0f deg\n', ...
            s.name, s.t, s.cycle, s.phase_abs_deg, s.phase_cycle_frac_angle, s.phase_cycle_frac_time, s.cycle_sweep_deg);
    end
end

%% Figure
if ~isempty(figDir)
    try
        fig = figure('Color', 'w', 'Position', [100 100 1000 420], 'Visible', 'off');
        set(fig, 'DefaultAxesFontName', 'Arial', 'DefaultTextFontName', 'Arial');
        ax1 = subplot(1, 2, 1); hold(ax1, 'on');
        plot(ax1, z(:,1), z(:,2), '-', 'Color', [0.8 0.8 0.8]);
        cols = lines(max(numel(stim), 2));
        for k = 1:numel(stim)
            if isnan(stim(k).cycle), continue; end
            cc = stim(k).cycle;
            b1 = min(max(sb(cc), 1), nb); b2 = min(max(sb(cc+1), 1), nb);
            plot(ax1, z(b1:b2,1), z(b1:b2,2), '-', 'Color', cols(k,:), 'LineWidth', 1.8);
            bs = min(max(floor((stim(k).t - t0) / dt) + 1, 1), nb);
            plot(ax1, z(bs,1), z(bs,2), 'v', 'MarkerFaceColor', cols(k,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 9);
            text(ax1, z(bs,1), z(bs,2), ['  ' stim(k).name], 'FontSize', 10);
        end
        axis(ax1, 'equal'); xlabel(ax1, 'jPC1'); ylabel(ax1, 'jPC2');
        title(ax1, 'Cycle containing each stimulus');
        pax = polaraxes('Position', [0.55 0.12 0.38 0.76]); hold(pax, 'on');
        for k = 1:numel(stim)
            polarplot(pax, [stim(k).phase_abs_rad stim(k).phase_abs_rad], [0 1], '-', 'Color', cols(k,:), 'LineWidth', 2);
            polarscatter(pax, stim(k).phase_abs_rad, 1, 60, cols(k,:), 'filled');
            text(pax, stim(k).phase_abs_rad, 1.12, stim(k).name, 'HorizontalAlignment', 'center');
        end
        rlim(pax, [0 1.2]); pax.RTickLabel = {};
        title(pax, 'Stimulus angle in the jPC plane');
        sgtitle(sprintf('Stimulus phase (jPC plane) | %s', recID), 'FontName', 'Arial');
        fn = fullfile(figDir, '06e_stimulus_phase_jpc.png');
        try, exportgraphics(fig, fn, 'Resolution', 300); catch, print(fig, fn, '-dpng', '-r300'); end
        close(fig);
    catch MEp
        warning('fn_stimulusPhaseJPC:PlotFailed', 'Phase figure skipped: %s', MEp.message);
    end
end
end
