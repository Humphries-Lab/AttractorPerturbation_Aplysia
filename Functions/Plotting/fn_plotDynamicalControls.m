function fn_plotDynamicalControls(dyn, surr, recording_ID, figuresDir) %#ok<INUSD>
% FN_PLOTDYNAMICALCONTROLS  Controls for the dynamical-epoch analysis.
%   A  Block contrast (within - between time blocks) for every possible split of the cycle
%      sequence. Red dashed lines = P9 and the manual epoch starts; the dot is the strongest split.
%      This panel does not use clusters.
%   B  Permutation null for the mean distance between the landmarks and the nearest cluster change.
%   C  Rotation strength (% variance in the first jPC plane, 6 PCs): data vs surrogates.
%   D  Within - between CLUSTER alignment (mean alignment of same-cluster minus different-cluster
%      cycle pairs): data vs surrogates. Clusters are built from the same matrix, so only the
%      comparison with surrogates, which go through the same procedure, is informative.
%   E  PCs needed for 80% of the variance: data vs surrogates.
%   Surrogates = every neuron circularly shifted by an independent random amount.
%
%   FN_PLOTDYNAMICALCONTROLS(dyn, surr, recording_ID, figuresDir)

if ~isfield(dyn, 'ok') || ~dyn.ok, return; end
blue = [0.1 0.3 0.7]; gry = [0.5 0.5 0.5]; red = [0.8 0.1 0.1];

fig = figure('Name', 'Dynamical epochs: controls', 'Position', [100 100 1400 820], 'Visible', 'on');

%% A - block contrast over all splits
axA = subplot(2, 6, 1:3); hold on;
bd = dyn.boundary;
if bd.ok
    plot(bd.t, bd.contrast, '-', 'Color', blue, 'LineWidth', 1.5);
    plot(bd.best_t, bd.best_contrast, 'o', 'MarkerFaceColor', blue, 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
    yl = ylim; yl(2) = yl(2) + 0.12 * diff(yl); ylim(yl);
    lm = dyn.landmarks;
    for k = 1:numel(lm.t)
        if lm.t(k) > dyn.t_range(1) && lm.t(k) < dyn.t_range(2)
            xline(lm.t(k), '--', 'Color', red, 'LineWidth', 1.0);
            text(lm.t(k), yl(2) - mod(k, 2) * 0.07 * diff(yl), lm.short{k}, 'Color', red, ...
                'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                'BackgroundColor', 'w', 'Margin', 0.5);
        end
    end
    xlim(dyn.t_range);
else
    text(0.5, 0.5, 'Too few cycles', 'Units', 'normalized', 'HorizontalAlignment', 'center');
end
xlabel('Split time (s)'); ylabel('Block contrast (within - between time blocks)');

%% B - landmark permutation
axB = subplot(2, 6, 4:6); hold on;
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

%% C, D, E - data vs surrogates
if isstruct(surr) && isfield(surr, 'real')
    axC = subplot(2, 6, 7:8);
    drawComparison(axC, 100 * surr.real.rot6, 100 * surr.rot6(:), 'Rotation strength', ...
        '% variance in first jPC plane (6 PCs)', surr.p.rot6, surr.n, blue, gry);
    axD = subplot(2, 6, 9:10);
    drawComparison(axD, surr.real.contrast, surr.contrast(:), 'Cluster separation', ...
        'Within - between cluster alignment', surr.p.contrast, surr.n, blue, gry);
    axE = subplot(2, 6, 11:12);
    drawComparison(axE, surr.real.ndim80, surr.ndim80(:), 'Dimensionality', ...
        'PCs for 80% of variance', surr.p.ndim80, surr.n, blue, gry);
    axes_all = [axA axB axC axD axE];
else
    axC = subplot(2, 6, 7:12);
    text(0.5, 0.5, 'Surrogates not run', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis off
    axes_all = [axA axB axC];
end

fn_formatAxes(fig, 'Arial', 10);
fn_letterAxes(axes_all, 14, 'Arial');
exportgraphics(fig, fullfile(figuresDir, '06c_dynamical_controls.png'), 'Resolution', 500);
drawnow;   % figure stays open on screen
end

function drawComparison(ax, realVal, surrVals, ttl, ylab, pval, nS, blue, gry)
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
