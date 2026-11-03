% Run_Cohort_Analysis_Sensitization.m
clear; clc; close all;

% 1. Get Directory
RESULTS_DIR = uigetdir(pwd, 'Select the Sensitization Results Folder');
if RESULTS_DIR == 0; return; end

% 2. Call the dedicated function (Matches main pipeline architecture)
T_cohort = fn_cohortAnalysis_sensitization(RESULTS_DIR);

disp('Cohort analysis complete. Data exported to Excel.');
