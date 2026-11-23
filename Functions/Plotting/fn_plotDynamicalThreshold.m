function fn_plotDynamicalThreshold(dyn, recording_ID, figuresDir) %#ok<INUSD>
% FN_PLOTDYNAMICALTHRESHOLD  Two single-question figures, one per slide.
%   Figure 1 (06b1_cluster_separation.png)
%     Are cycles in the same cluster more alike than cycles in different clusters?
%     Distribution of pairwise cycle alignment, same- vs different-cluster pairs, with the
%     shuffle-chance level (95th percentile) and the threshold that best separates the two.
%     (Clusters are built from the same alignment matrix, so some separation is built in.)
%   Figure 2 (06b2_drift_check.png)
%     Is Evoked-Recovery alignment lower than the time between the cycles predicts?
%     Alignment vs time between cycles. The line is fitted to same-epoch pairs only: solid over
%     the fitted range, dashed where it is extrapolated.
%
%   FN_PLOTDYNAMICALTHRESHOLD(dyn, recording_ID, figuresDir)

if ~isfield(dyn, 'ok') || ~dyn.ok
    return
end
blue = [0.1 0.3 0.7]; gry = [0.5 0.5 0.5]; red = [0.8 0.1 0.1];
st = dyn.sim_thresh; pr = dyn.pairs;

%% Figure 1 - same vs different cluster
fig1 = figure('Name', 'Cluster separation', 'Position', [100 100 560 430], 'Visible', 'on');
hold on;
edges = 0:0.05:1;
w = pr.sim(pr.same_cluster); b = pr.sim(~pr.same_cluster);
hB = gobjects(0); hW = gobjects(0); top = 0;
if ~isempty(b)
    hB = histogram(b, edges, 'FaceColor', gry, 'EdgeColor', 'w', 'Normalization', 'probability');
    top = max(top, max(hB.Values));
end
if ~isempty(w)
    hW = histogram(w, edges, 'FaceColor', blue, 'EdgeColor', 'w', 'Normalization', 'probability', 'FaceAlpha', 0.6);
    top = max(top, max(hW.Values));
end
xlim([0 1]); ylim([0 max(top, eps) * 1.35]);
if ~isnan(st.chance95)
    xline(st.chance95, ':', 'Chance', 'Color', 'k', 'LineWidth', 1.5, 'FontSize', 9, ...
        'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left');
end
if ~isnan(st.separation)
    xline(st.separation, '--', 'Threshold', 'Color', red, 'LineWidth', 1.5, 'FontSize', 9, ...
        'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left');
end
hs = gobjects(0); ls = {};
if ~isempty(hW), hs(end+1) = hW; ls{end+1} = 'Same cluster'; end
if ~isempty(hB), hs(end+1) = hB; ls{end+1} = 'Different cluster'; end
if ~isempty(hs), legend(hs, ls, 'Location', 'northeast', 'Box', 'off', 'FontSize', 9); end
if ~isnan(st.auc)
    text(0.97, 0.70, sprintf('AUC = %.2f', st.auc), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'FontSize', 9);
end
xlabel('Subspace alignment between two cycles'); ylabel('Proportion of cycle pairs');
fn_formatAxes(fig1, 'Arial', 10);
exportgraphics(fig1, fullfile(figuresDir, '06b1_cluster_separation.png'), 'Resolution', 500);
drawnow;   % figure stays open on screen

%% Figure 2 - drift check
fig2 = figure('Name', 'Drift check', 'Position', [140 100 600 430], 'Visible', 'on');
hold on;
m0 = pr.type == 0; m1 = pr.type == 1; m2 = pr.type == 2;
hp = gobjects(0); lp = {};
if any(m0), hp(end+1) = plot(pr.lag_s(m0), pr.sim(m0), '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 10); lp{end+1} = 'Other pairs'; end
if any(m1), hp(end+1) = plot(pr.lag_s(m1), pr.sim(m1), '.', 'Color', blue, 'MarkerSize', 12); lp{end+1} = 'Same epoch'; end
if any(m2), hp(end+1) = plot(pr.lag_s(m2), pr.sim(m2), '.', 'Color', red, 'MarkerSize', 12); lp{end+1} = 'Evoked-Recovery'; end
hasFit = isfield(dyn.drift, 'ok') && dyn.drift.ok;
if hasFit
    fitFun = @(x) dyn.drift.intercept + dyn.drift.slope_per_min / 60 * x;
    lagFit = [0 max(pr.lag_s(m1))];                    % range the line was fitted on (same-epoch pairs)
    lagExt = [max(pr.lag_s(m1)) max(pr.lag_s)];        % beyond it: extrapolated
    hp(end+1) = plot(lagFit, fitFun(lagFit), 'k-', 'LineWidth', 1.5);
    lp{end+1} = 'Fit to same-epoch pairs';
    if lagExt(2) > lagExt(1)
        hp(end+1) = plot(lagExt, fitFun(lagExt), '--', 'Color', gry, 'LineWidth', 1.5);
        lp{end+1} = 'Extrapolated';
    end
end
if ~isnan(st.chance95)
    hp(end+1) = yline(st.chance95, ':', 'Color', 'k', 'LineWidth', 1.2); lp{end+1} = 'Chance';
end
if ~isempty(hp), legend(hp, lp, 'Location', 'northeast', 'Box', 'off', 'FontSize', 9); end
xlabel('Time between cycles (s)'); ylabel('Subspace alignment between two cycles'); ylim([0 1]);
if hasFit
    title(sprintf('Evoked-Recovery: %.2f observed, %.2f expected', dyn.drift.obs_ER, dyn.drift.pred_ER), ...
        'FontWeight', 'normal', 'FontSize', 10);
end
fn_formatAxes(fig2, 'Arial', 10);
exportgraphics(fig2, fullfile(figuresDir, '06b2_drift_check.png'), 'Resolution', 500);
drawnow;   % figure stays open on screen
end
