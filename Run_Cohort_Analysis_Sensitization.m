% Run_Cohort_Analysis_Sensitization.m
% Pools and plots cohort-level trends for Sensitization experiments.
% Reads native '*_results_sensitization.mat' files directly.

clear; clc; close all;

%% 1. Setup
% Point this to your separate sensitization results folder
RESULTS_DIR = uigetdir(pwd, 'Select the Sensitization Results Folder');
if RESULTS_DIR == 0; return; end

files = dir(fullfile(RESULTS_DIR, '*_results_sensitization.mat'));
nAnimals = length(files);

if nAnimals == 0
    error('No sensitization results files found in %s', RESULTS_DIR);
end

%% 2. Extract Data Across the Cohort
fprintf('Pooling %d sensitization recordings...\n', nAnimals);

for i = 1:nAnimals
    filePath = fullfile(RESULTS_DIR, files(i).name);
    S = load(filePath);
    
    nStim = length(S.evoked_idx); % Usually 5
    
    % Initialize arrays on first pass
    if i == 1
        cohort.recording_IDs = cell(nAnimals, 1);
        cohort.PR_norm_baseline = nan(nAnimals, nStim);
        cohort.PR_norm_evoked   = nan(nAnimals, nStim);
        cohort.RR_evoked        = nan(nAnimals, nStim);
        cohort.cross_RR         = nan(nAnimals, nStim);
        cohort.onset_t          = nan(nAnimals, nStim);
    end
    
    cohort.recording_IDs{i} = S.cfg.recording_ID;
    
    % Extract dimensionality (PR)
    cohort.PR_norm_baseline(i, :) = S.PR_norm(S.baseline_idx);
    cohort.PR_norm_evoked(i, :)   = S.PR_norm(S.evoked_idx);
    
    % Extract self-recurrence density
    cohort.RR_evoked(i, :) = S.ev_recur_density;
    
    % Extract cross-recurrence (Evoked 1 vs Evoked i)
    if isfield(S, 'cross_recur_density')
        cohort.cross_RR(i, :) = S.cross_recur_density;
    end
    
    % Extract onset times (excluding failed detections)
    valid_onset = S.t_attractor_onset_all;
    valid_onset(~S.onset_detected_all) = NaN;
    cohort.onset_t(i, :) = valid_onset;
end

%% 3. Plotting Cohort Trends
stims = 1:nStim;
x_labels = arrayfun(@(x) sprintf('S%d',x), stims, 'UniformOutput', false);

figure('Name', 'Sensitization Cohort Trends', 'Position', [100, 100, 1200, 800], 'Color', 'w');

% Helper function for standard error
sem = @(x) nanstd(x, 1) ./ sqrt(sum(~isnan(x), 1));

% --- A. Dimensionality (PR_norm) ---
subplot(2,2,1); hold on;
plot(stims, cohort.PR_norm_evoked, 'Color', [0.7 0.7 0.7]); % Individual animals
errorbar(stims, nanmean(cohort.PR_norm_evoked, 1), sem(cohort.PR_norm_evoked), ...
    'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k');
xticks(stims); xticklabels(x_labels);
ylabel('Evoked PR_{norm}'); title('State Space Dimensionality');
grid on; box on;

% --- B. Attractor Density (RR) ---
subplot(2,2,2); hold on;
plot(stims, cohort.RR_evoked, 'Color', [0.7 0.7 0.7]); 
errorbar(stims, nanmean(cohort.RR_evoked, 1), sem(cohort.RR_evoked), ...
    'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xticks(stims); xticklabels(x_labels);
ylabel('Recurrence Rate (RR)'); title('Attractor Density');
grid on; box on;

% --- C. Cross-Recurrence (State Space Shift) ---
subplot(2,2,3); hold on;
plot(stims, cohort.cross_RR, 'Color', [0.7 0.7 0.7]); 
errorbar(stims, nanmean(cohort.cross_RR, 1), sem(cohort.cross_RR), ...
    'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
xticks(stims); xticklabels(x_labels);
ylabel('Cross-RR (vs S1)'); title('State Space Alignment to S1');
grid on; box on;

% --- D. Onset Time ---
subplot(2,2,4); hold on;
plot(stims, cohort.onset_t, 'Color', [0.7 0.7 0.7]); 
errorbar(stims, nanmean(cohort.onset_t, 1), sem(cohort.onset_t), ...
    'g-o', 'LineWidth', 2, 'MarkerFaceColor', 'g');
xticks(stims); xticklabels(x_labels);
ylabel('Onset Time (s)'); title('Attractor Onset Dynamics');
grid on; box on;

sgtitle('Sensitization Protocol: Cohort Summary', 'FontSize', 16, 'FontWeight', 'bold');

% 4. Save Cohort Data
save(fullfile(RESULTS_DIR, 'Cohort_Summary_Sensitization.mat'), 'cohort');
fprintf('Cohort analysis complete. Saved to %s\n', fullfile(RESULTS_DIR, 'Cohort_Summary_Sensitization.mat'));

%% 5. Export Cohort Data to Excel (Wide Format)
% Converts arrays to a flat table layout optimized for SPSS/R repeated measures.

% 1. Start with the recording IDs
T_out = table(cohort.recording_IDs, 'VariableNames', {'Recording_ID'});

% 2. Define the metrics to flatten and their column prefixes
var_names = {'PR_norm_baseline', 'PR_norm_evoked', 'RR_evoked', 'cross_RR', 'onset_t'};
prefixes  = {'PR_base', 'PR_evoked', 'RR', 'Cross_RR', 'Onset_Time'};

% 3. Loop through each metric and each stimulus to generate columns dynamically
for m = 1:length(var_names)
    metric_data = cohort.(var_names{m});
    for s = 1:nStim
        % Create column name, e.g., 'PR_evoked_S1'
        col_name = sprintf('%s_S%d', prefixes{m}, s);
        % Assign the column data (all animals for this stimulus)
        T_out.(col_name) = metric_data(:, s);
    end
end

% 4. Save to Excel
excelPath = fullfile(RESULTS_DIR, 'Cohort_Summary_Sensitization.xlsx');
writetable(T_out, excelPath);
fprintf('Excel export complete. Saved wide-format table to %s\n', excelPath);
