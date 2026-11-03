function T_out = fn_cohortAnalysis_sensitization(resultsDir)
    % FN_COHORTANALYSIS_SENSITIZATION Pools Sensitization results into a wide-format table.
    
    files = dir(fullfile(resultsDir, '*_results_sensitization.mat'));
    nAnimals = length(files);
    if nAnimals == 0
        error('No sensitization results files found in %s', resultsDir);
    end

    % Loop and extract
    for i = 1:nAnimals
        filePath = fullfile(resultsDir, files(i).name);
        S = load(filePath);
        nStim = length(S.evoked_idx);
        
        if i == 1
            cohort.recording_IDs = cell(nAnimals, 1);
            cohort.PR_base = nan(nAnimals, nStim);
            cohort.PR_evoked = nan(nAnimals, nStim);
            cohort.RR = nan(nAnimals, nStim);
            cohort.Cross_RR = nan(nAnimals, nStim);
            cohort.Onset_Time = nan(nAnimals, nStim);
        end
        
        cohort.recording_IDs{i} = S.cfg.recording_ID;
        cohort.PR_base(i, :) = S.PR_norm(S.baseline_idx);
        cohort.PR_evoked(i, :) = S.PR_norm(S.evoked_idx);
        cohort.RR(i, :) = S.ev_recur_density;
        if isfield(S, 'cross_recur_density'), cohort.Cross_RR(i, :) = S.cross_recur_density; end
        
        valid_onset = S.t_attractor_onset_all;
        valid_onset(~S.onset_detected_all) = NaN;
        cohort.Onset_Time(i, :) = valid_onset;
    end

    % Convert to wide-format table
    T_out = table(cohort.recording_IDs, 'VariableNames', {'Recording_ID'});
    var_names = {'PR_base', 'PR_evoked', 'RR', 'Cross_RR', 'Onset_Time'};
    
    for m = 1:length(var_names)
        metric_data = cohort.(var_names{m});
        for s = 1:nStim
            col_name = sprintf('%s_S%d', var_names{m}, s);
            T_out.(col_name) = metric_data(:, s);
        end
    end

    % Save Excel
    excelPath = fullfile(resultsDir, 'Cohort_Summary_Sensitization.xlsx');
    writetable(T_out, excelPath);
    fprintf('Saved sensitization cohort table to %s\n', excelPath);
end
