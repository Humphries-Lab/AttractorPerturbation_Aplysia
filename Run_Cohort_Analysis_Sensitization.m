% Run_Cohort_Analysis_Sensitization.m
% 
% Master script to extract, pool, and visualize cohort-level trends 
% for repeated-stimuli Sensitization experiments.

clear; clc; close all;

%% 1. Select Folder and Run Extraction
RESULTS_DIR = uigetdir(pwd, 'Select the Sensitization Results Folder');
if RESULTS_DIR == 0; return; end

% Call the cohort extraction function
T_cohort = fn_cohortAnalysis_sensitization(RESULTS_DIR);

% Load the saved cohort struct to easily plot the 3D matrices
load(fullfile(RESULTS_DIR, 'Cohort_Summary_Sensitization.mat'), 'cohort');

%% 2. Plot Cohort Trends
nStim = size(cohort.PR_evoked, 2);
stims = 1:nStim;
x_labels = arrayfun(@(x) sprintf('S%d',x), stims, 'UniformOutput', false);

figure('Name', 'Sensitization Cohort Summary', 'Position', [100, 100, 1200, 800], 'Color', 'w');
sem = @(x) nanstd(x, 1) ./ sqrt(sum(~isnan(x), 1)); % Standard Error calculation

% --- A. Dimensionality (PR_norm) ---
subplot(2,3,1); hold on;
plot(stims, cohort.PR_evoked, 'Color', [0.7 0.7 0.7]);
errorbar(stims, nanmean(cohort.PR_evoked, 1), sem(cohort.PR_evoked), 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k');
xticks(stims); xticklabels(x_labels);
ylabel('Evoked PR_{norm}'); title('Dimensionality');
grid on; box on;

% --- B. Attractor Density (RR) ---
subplot(2,3,2); hold on;
plot(stims, cohort.RR, 'Color', [0.7 0.7 0.7]); 
errorbar(stims, nanmean(cohort.RR, 1), sem(cohort.RR), 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xticks(stims); xticklabels(x_labels);
ylabel('RR'); title('Attractor Density');
grid on; box on;

% --- C. Onset Time ---
subplot(2,3,3); hold on;
plot(stims, cohort.Onset_Time, 'Color', [0.7 0.7 0.7]); 
errorbar(stims, nanmean(cohort.Onset_Time, 1), sem(cohort.Onset_Time), 'g-o', 'LineWidth', 2, 'MarkerFaceColor', 'g');
xticks(stims); xticklabels(x_labels);
ylabel('Onset Time (s)'); title('Convergence Speed');
grid on; box on;

% --- D. Subspace Alignment to S1 ---
subplot(2,3,4); hold on;
plot(stims, cohort.Align_to_S1, 'Color', [0.7 0.7 0.7]); 
errorbar(stims, nanmean(cohort.Align_to_S1, 1), sem(cohort.Align_to_S1), 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
xticks(stims); xticklabels(x_labels);
ylabel('Alignment (Cos \theta)'); title('Subspace Alignment to S1');
grid on; box on;

% --- E. Cross-Recurrence to S1 ---
subplot(2,3,5); hold on;
plot(stims, cohort.Cross_RR, 'Color', [0.7 0.7 0.7]); 
errorbar(stims, nanmean(cohort.Cross_RR, 1), sem(cohort.Cross_RR), 'm-o', 'LineWidth', 2, 'MarkerFaceColor', 'm');
xticks(stims); xticklabels(x_labels);
ylabel('Cross-RR vs S1'); title('Non-linear Identity vs S1');
grid on; box on;

sgtitle('Motor Program Dynamics Across Repeated Sensitization', 'FontSize', 16, 'FontWeight', 'bold');
