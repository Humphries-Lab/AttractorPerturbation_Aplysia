function T_out = fn_cohortAnalysis_sensitization(resultsDir)
    % FN_COHORTANALYSIS_SENSITIZATION Pools Sensitization results into a wide-format table.
    % Extracts Dimensionality (PR), Attractor Density (RR), Cross-Recurrence, 
    % Onset Time, and State Space Alignment (Alignment to S1).
    
    files = dir(fullfile(resultsDir, '*_results_sensitization.mat'));
    nAnimals = length(files);
    
    if nAnimals == 0
        error('No sensitization results files found in %s', resultsDir);
    end

    fprintf('Processing %d sensitization recordings...\n', nAnimals);

    % Loop and extract
    for i = 1:nAnimals
        filePath = fullfile(resultsDir, files(i).name);
        S = load(filePath);
        
        % Safely determine number of stimuli (fallback to 5)
        if isfield(S, 'evoked_idx')
            nStim = length(S.evoked_idx);
        else
            nStim = 5; 
        end
        
        % Initialize arrays on first pass
        if i == 1
            cohort.recording_IDs = cell(nAnimals, 1);
            cohort.PR_base       = nan(nAnimals, nStim);
            cohort.PR_evoked     = nan(nAnimals, nStim);
            cohort.RR            = nan(nAnimals, nStim);
            cohort.Cross_RR      = nan(nAnimals, nStim);
            cohort.Onset_Time    = nan(nAnimals, nStim);
            cohort.Align_to_S1   = nan(nAnimals, nStim);
        end
        
        % 1. ID Extraction
        if isfield(S, 'cfg') && isfield(S.cfg, 'recording_ID')
            cohort.recording_IDs{i} = S.cfg.recording_ID;
        else
            cohort.recording_IDs{i} = sprintf('Unknown_ID_%d', i);
        end
        
        % 2. Dimensionality (PR)
        if isfield(S, 'PR_norm') && isfield(S, 'baseline_idx') && isfield(S, 'evoked_idx')
            nBase = min(nStim, length(S.baseline_idx));
            nEvok = min(nStim, length(S.evoked_idx));
            
            if length(S.PR_norm) >= max(S.baseline_idx)
                cohort.PR_base(i, 1:nBase) = S.PR_norm(S.baseline_idx(1:nBase));
            end
            if length(S.PR_norm) >= max(S.evoked_idx)
                cohort.PR_evoked(i, 1:nEvok) = S.PR_norm(S.evoked_idx(1:nEvok));
            end
        end
        
        % 3. Attractor Density (RR)
        if isfield(S, 'ev_recur_density')
            nVals = min(nStim, length(S.ev_recur_density));
            cohort.RR(i, 1:nVals) = S.ev_recur_density(1:nVals);
        end
        
        % 4. Cross-Recurrence (Non-linear identity)
        if isfield(S, 'cross_recur_density')
            nVals = min(nStim, length(S.cross_recur_density));
            cohort.Cross_RR(i, 1:nVals) = S.cross_recur_density(1:nVals);
        end
        
        % 5. Onset Time & Detection Mask
        if isfield(S, 't_attractor_onset_all')
            nVals = min(nStim, length(S.t_attractor_onset_all));
            valid_onset = S.t_attractor_onset_all(1:nVals);
            
            % Apply detection mask to filter out failed onsets
            if isfield(S, 'onset_detected_all')
                mask = S.onset_detected_all(1:nVals);
                valid_onset(~mask) = NaN;
            end
            cohort.Onset_Time(i, 1:nVals) = valid_onset;
        end
        
        % 6. Subspace Alignment to S1 (Linear identity)
        if isfield(S, 'align_mat') && isfield(S, 'evoked_idx')
            % Ensure the alignment matrix is large enough to contain our evoked indices
            if size(S.align_mat, 1) >= max(S.evoked_idx) && size(S.align_mat, 2) >= max(S.evoked_idx)
                idx_S1 = S.evoked_idx(1); % Reference is Stimulus 1
                for s = 1:min(nStim, length(S.evoked_idx))
                    idx_curr = S.evoked_idx(s);
                    cohort.Align_to_S1(i, s) = S.align_mat(idx_S1, idx_curr);
                end
            end
        end
    end

    % Convert to wide-format table for repeated measures ANOVA / SPSS / R
    T_out = table(cohort.recording_IDs, 'VariableNames', {'Recording_ID'});
    var_names = {'PR_base', 'PR_evoked', 'RR', 'Cross_RR', 'Align_to_S1', 'Onset_Time'};
    
    for m = 1:length(var_names)
        metric_data = cohort.(var_names{m});
        for s = 1:size(metric_data, 2)
            col_name = sprintf('%s_S%d', var_names{m}, s);
            T_out.(col_name) = metric_data(:, s);
        end
    end

    % Save Excel
    excelPath = fullfile(resultsDir, 'Cohort_Summary_Sensitization.xlsx');
    writetable(T_out, excelPath);
    
    % Also save the raw struct for MATLAB plotting ease
    matPath = fullfile(resultsDir, 'Cohort_Summary_Sensitization.mat');
    save(matPath, 'cohort', 'T_out');
    
    fprintf('Success! Saved cohort table to %s\n', excelPath);
end
