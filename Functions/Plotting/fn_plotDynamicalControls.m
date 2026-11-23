function fn_plotDynamicalControls(dyn, surr, recording_ID, figuresDir) %#ok<INUSD>
% FN_PLOTDYNAMICALCONTROLS  Two single-question figures, one per slide.
%   Figure 1 (06c1_stimulus_alignment.png): do the changes in the cycles line up with the stimuli?
%     A  Block contrast for every cut of the cycle sequence into an early and a late block:
%        mean alignment within the two blocks minus mean alignment between them. The dot is the
%        strongest cut. Red dashed lines = P9, start of Evoked, C2, start of Recovery (landmarks).
%        This panel does not use clusters.
%     B  Are the times at which the cluster changes fall closer to the landmarks than chance?
%        Grey: mean distance from the landmarks to the nearest cluster change for randomly placed
%        change points. Red: observed.
%   Figure 2 (06c2_surrogate_controls.png): could the cycles and clusters arise without coordinated
%     population activity? Data (diamond) vs surrogates (every neuron circularly shifted by an
%     independent random amount).
%     A  Rotation strength: % variance in the first jPC plane (6 PCs).
%     B  Within - between cluster alignment: same-cluster minus different-cluster cycle pairs.
%        Clusters are built from the same matrix, so only the comparison with surrogates, which go
%        through the same procedure, is informative.
%     C  PCs needed for 80% of the variance.
%
%   FN_PLOTDYNAMICALCONTROLS(dyn, surr, recording_ID, figuresDir)

if ~isfield(dyn, 'ok') || ~dyn.ok, return; end
blue = [0.1 0.3 0.7]; gry = [0.5 0.5 0.5]; red = [0.8 0.1 0.1];

%% Figure 1 - stimulus alignment
fig1 = figure('Name', 'Dynamical epochs: stimulus alignment', 'Position', [100 100 1100 430], 'Visible', 'on');
axA = subplot(1, 2, 1); hold on;
bd = dyn.boundary;
if bd.ok
    plot(bd.t, bd.contrast, '-', 'Color', blue, 'LineWidth', 1.5);
    plot(bd.best_t, bd.best_contrast, 'o', 'MarkerFaceColor', blue, 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
    yl = ylim; yl(2) = yl(2) + 0.12 * diff(yl); ylim(yl);
    lm = dyn.landmarks;
    lmShort = strrep(lm.short, 'Per', 'C2');          % the start of the Perturbed epoch is C2
    for k = 1:numel(lm.t)
        if lm.t(k) > dyn.t_range(1) && lm.t(k) < dyn.t_range(2)
            xline(lm.t(k), '--', 'Color', red, 'LineWidth', 1.0);
            text(lm.t(k), yl(2) - mod(k, 2) * 0.07 * diff(yl), lmShort{k}, 'Color', red, ...
                'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                'BackgroundColor', 'w', 'Margin', 0.5);
        end
    end
    xlim(dyn.t_range);
else
    text(0.5, 0.5, 'Too few cycles', 'Units', 'normalized', 'HorizontalAlignment', 'center');
end
xlabel('Cut time (s)'); ylabel('Block contrast (alignment within - between blocks)');

axB = subplot(1, 2, 2); hold on;
lt = dyn.landmark_test;
if lt.ok
    histogram(lt.null, 30, 'FaceColor', gry, 'EdgeColor', 'w', 'Normalization', 'probability');
    xline(lt.observed, '-', 'Color', red, 'LineWidth', 2);
    yl = ylim; ylim([0 yl(2) * 1.2]);
    text(0.97, 0.93, sprintf('p = %.3f', lt.p), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 9);
else
    text(0.5, 0.5, 'No change points in window', 'Units', 'normalized', 'HorizontalAlignment', 'center');
end
xlabel('Mean distance, landmark to nearest cluster change (s)'); ylabel('Proportion of permutations');
fn_formatAxes(fig1, 'Arial', 10);
fn_letterAxes([axA axB], 14, 'Arial');
exportgraphics(fig1, fullfile(figuresDir, '06c1_stimulus_alignment.png'), 'Resolution', 500);
drawnow;   % figure stays open on screen

%% Figure 2 - surrogate controls
fig2 = figure('Name', 'Dynamical epochs: surrogate controls', 'Position', [140 100 1100 430], 'Visible', 'on');
if isstruct(surr) && isfield(surr, 'real')
    axC = subplot(1, 3, 1);
    drawComparison(axC, 100 * surr.real.rot6, 100 * surr.rot6(:), ...
        '% variance in first jPC plane (6 PCs)', surr.p.rot6, surr.n, blue, gry);
    axD = subplot(1, 3, 2);
    drawComparison(axD, surr.real.contrast, surr.contrast(:), ...
        'Within - between cluster alignment', surr.p.contrast, surr.n, blue, gry);
    axE = subplot(1, 3, 3);
    drawComparison(axE, surr.real.ndim80, surr.ndim80(:), ...
        'PCs for 80% of variance', surr.p.ndim80, surr.n, blue, gry);
    axes_all = [axC axD axE];
else
    axC = subplot(1, 1, 1);
    text(0.5, 0.5, 'Surrogates not run', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis off
    axes_all = axC;
end
fn_formatAxes(fig2, 'Arial', 10);
fn_letterAxes(axes_all, 14, 'Arial');
exportgraphics(fig2, fullfile(figuresDir, '06c2_surrogate_controls.png'), 'Resolution', 500);
drawnow;   % figure stays open on screen
end

function drawComparison(ax, realVal, surrVals, ylab, pval, nS, blue, gry)
% Data (diamond) vs surrogates (grey dots) on their own axis, with units and the p-value.
axes(ax); hold on; %#ok<LAXES>
sv = surrVals(isfinite(surrVals));
rs = rng; rng(1);
plot(2 + 0.3 * (rand(numel(sv), 1) - 0.5), sv, 'o', 'MarkerFaceColor', gry, 'MarkerEdgeColor', 'w', 'MarkerSize', 7);
rng(rs);
plot(1, realVal, 'd', 'MarkerFaceColor', blue, 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
set(ax, 'XTick', [1 2], 'XTickLabel', {'Data', sprintf('Surrogates (n=%d)', nS)});
xlim([0.4 2.6]);
allv = [sv(:); realVal];
pad = 0.15 * max(max(allv) - min(allv), eps);
ylim([min(allv) - pad, max(allv) + 2.2 * pad]);
ylabel(ylab);
if isfinite(pval)
    text(0.5, 0.97, sprintf('p = %.2f', pval), 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8);
end
end
