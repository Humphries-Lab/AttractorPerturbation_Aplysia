function fn_plotDynamicalEpochs(dyn, recording_ID, figuresDir)
% FN_PLOTDYNAMICALEPOCHS  Summary figure for fn_dynamicalEpochs.
%   (a) cycle x cycle subspace alignment   (b) cycle clusters on the time axis
%   (c) cycles per cluster per protocol epoch   (d) mean alignment between epochs
%
%   In (a) and (b) the manual protocol epoch starts (Evoked, Perturbed = C2, Recovery) are
%   marked only to label cycles. They are not estimates of when the perturbation ends.
%   In (b) the two stimuli (P9, C2) are drawn in red and the manual epoch starts in grey.
%
%   FN_PLOTDYNAMICALEPOCHS(dyn, recording_ID, figuresDir)

if ~isfield(dyn, 'ok') || ~dyn.ok
    fprintf('  Dynamical epochs: nothing to plot (%s)\n', dyn.message);
    return
end

nSeg = numel(dyn.segment_s); nG = dyn.nGroups; nE = numel(dyn.ep_labels);
cmap = lines(max(nG, 1));
bnd  = dyn.ep_bounds(2:end, 1);     % manual protocol epoch starts (Evoked, Perturbed, Recovery)
red  = [0.85 0.1 0.1]; gry = [0.45 0.45 0.45];

fig = figure('Name', 'Dynamical epochs', 'Position', [100 100 1400 800], 'Visible', 'off');

% (a) cycle-by-cycle alignment
subplot(2,2,1);
imagesc(dyn.similarity); clim([0 1]); colormap(gca, hot); colorbar; axis square; hold on;
for k = 1:numel(bnd)
    idx = find(dyn.cycle_mid >= bnd(k), 1, 'first');
    if ~isempty(idx) && idx > 1
        xline(idx - 0.5, 'w-', 'LineWidth', 1); yline(idx - 0.5, 'w-', 'LineWidth', 1);
    end
end
xlabel('Cycle #'); ylabel('Cycle #'); title('Cycle-by-cycle subspace alignment');

% (b) clusters on the time axis
subplot(2,2,2); hold on;
for c = 1:nSeg
    t1 = dyn.cycle_t(c,1); t2 = dyn.cycle_t(c,2);
    patch([t1 t2 t2 t1], [0 0 1 1], cmap(dyn.grps(c),:), 'EdgeColor', 'w');
end
xl = dyn.t_range; xlim(xl); ylim([0 1.24]);

% stimuli (real events): P9 and C2 (= start of the Perturbed epoch)
stimT = dyn.landmarks.t(1); stimN = {'P9'};
iPer = find(strcmp(dyn.ep_labels, 'Perturbed'), 1);
if ~isempty(iPer), stimT(end+1) = dyn.ep_bounds(iPer, 1); stimN{end+1} = 'C2'; end

% manual epoch starts that are not stimuli (Evoked start, Recovery start): grey dotted
for k = 1:numel(bnd)
    if bnd(k) > xl(1) && bnd(k) < xl(2) && ~any(abs(stimT - bnd(k)) < 1e-9)
        xline(bnd(k), ':', 'Color', gry, 'LineWidth', 1.5);
    end
end
for k = 1:numel(stimT)
    if stimT(k) > xl(1) && stimT(k) < xl(2)
        xline(stimT(k), '-', 'Color', red, 'LineWidth', 1.8);
        text(stimT(k), 1.18, stimN{k}, 'Color', red, 'HorizontalAlignment', 'center', ...
            'FontSize', 9, 'FontWeight', 'bold');
    end
end
for e = 1:nE
    a = max(dyn.ep_bounds(e,1), xl(1)); b = min(dyn.ep_bounds(e,2), xl(2));
    if b > a
        text((a+b)/2, 1.07, dyn.ep_labels{e}, 'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end
hL = gobjects(1, nG);
for g = 1:nG, hL(g) = patch(NaN, NaN, cmap(g,:)); end
legend(hL, arrayfun(@(g) sprintf('Cluster %d', g), 1:nG, 'UniformOutput', false), 'Location', 'eastoutside');
set(gca, 'YTick', []); xlabel('Time (s)'); box on;
title('Dynamical epochs (cycle clusters)');
sub = {'Red: stimuli. Grey dotted: manual epoch starts (labels only)'};
idxR = find(strcmp(dyn.landmarks.name, 'start of Recovery'), 1);
if ~isempty(idxR) && ~isnan(dyn.landmarks.lag_to_nearest_change_s(idxR))
    sub{end+1} = sprintf('Nearest cluster change to the Recovery start: %+.0f s', dyn.landmarks.lag_to_nearest_change_s(idxR));
end
subtitle(sub, 'FontSize', 8);

% (c) cluster x protocol epoch counts
subplot(2,2,3);
imagesc(dyn.cluster_by_epoch); colormap(gca, flipud(gray)); axis tight;
for g = 1:nG
    for e = 1:nE
        text(e, g, sprintf('%d', dyn.cluster_by_epoch(g,e)), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
end
set(gca, 'XTick', 1:nE, 'XTickLabel', dyn.ep_labels, 'YTick', 1:nG, 'XTickLabelRotation', 30);
xlabel('Protocol epoch'); ylabel('Cluster'); title('Cycles per cluster in each protocol epoch');

% (d) alignment between protocol epochs
subplot(2,2,4);
imagesc(dyn.epoch_sim); clim([0 1]); colormap(gca, hot); colorbar; axis square;
for a = 1:nE
    for b = 1:nE
        if ~isnan(dyn.epoch_sim(a,b))
            text(b, a, sprintf('%.2f', dyn.epoch_sim(a,b)), 'HorizontalAlignment', 'center', ...
                'Color', 'w', 'FontWeight', 'bold');
        end
    end
end
set(gca, 'XTick', 1:nE, 'XTickLabel', dyn.ep_labels, 'YTick', 1:nE, 'YTickLabel', dyn.ep_labels, ...
    'XTickLabelRotation', 30);
title('Mean alignment between epochs (diagonal = within)');

sgtitle(sprintf('Dynamical epochs (jPCA cycles) | %s', recording_ID), 'Interpreter', 'none');
annotation(fig, 'textbox', [0.05 0.002 0.9 0.03], 'String', ...
    ['White lines in (a): start of Evoked, Perturbed (C2) and Recovery. Grey dotted lines in (b): start of Evoked (P9 + buffer) and ' ...
     'Recovery (C2 + delay). These are manual protocol labels, not estimates of when the perturbation ends.'], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9);
exportgraphics(fig, fullfile(figuresDir, '06_dynamical_epochs.png'), 'Resolution', 500);
close(fig);
end
