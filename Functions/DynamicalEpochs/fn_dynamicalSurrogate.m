function surr = fn_dynamicalSurrogate(spike_conv, fs, t_range_s, cfg, ep_bounds, ep_labels, dyn_real)
% FN_DYNAMICALSURROGATE  Null control for fn_dynamicalEpochs (Exact resolution, maximum parallel throughput).

n = cfg.DYN_NSURR;
nanv = nan(1, n);
surr = struct('n', n, 'rot6', nanv, 'ndim80', nanv, 'n_cycles', nanv, 'n_clusters', nanv, ...
    'auc', nanv, 'contrast', nanv, 'n_changes', nanv);
surr.real = struct('rot6', dyn_real.rot6, 'ndim80', dyn_real.ndim80, ...
    'n_cycles', numel(dyn_real.segment_s), 'n_clusters', dyn_real.nGroups, ...
    'auc', dyn_real.sim_thresh.auc, ...
    'contrast', dyn_real.sim_thresh.within_mean - dyn_real.sim_thresh.between_mean, ...
    'n_changes', numel(dyn_real.change_t));

%% 1. Bin once at EXACT cfg.DYN_BIN_MS resolution (zero compromise to data or statistics)
[nN, nFrames] = size(spike_conv);
bf      = max(round(fs * cfg.DYN_BIN_MS / 1000), 1);
bin_eff = 1000 * bf / fs;
f1 = max(round(t_range_s(1) * fs) + 1, 1);
f2 = min(round(t_range_s(2) * fs), nFrames);
nb = floor((f2 - f1 + 1) / bf);
Xb = zeros(nN, nb);
for b0 = 1:500:nb
    b1  = min(b0 + 499, nb);
    seg = double(spike_conv(:, f1 + (b0-1)*bf : f1 + b1*bf - 1));
    Xb(:, b0:b1) = reshape(mean(reshape(seg, nN, bf, b1-b0+1), 2), nN, b1-b0+1);
end
fs_b      = 1000 / bin_eff;
t_range_b = [0, nb * bin_eff / 1000];

%% 2. Load cache of earlier surrogates (same settings, same seeds)
key = struct('bin_ms', bin_eff, 'min_seg', cfg.DYN_MIN_SEG_S, 'rec_thresh', cfg.DYN_REC_THRESH, ...
    'seed', cfg.DYN_SEED, 'nN', nN, 'nb', nb, 'fs', fs, 't_range', t_range_s);
cacheFile = '';
if isfield(cfg, 'RESULTS_DIR') && isfield(cfg, 'recording_ID')
    cacheFile = fullfile(cfg.RESULTS_DIR, sprintf('%s_surrogates.mat', cfg.recording_ID));
end
R = cell(1, n);
Rcached = {};
if ~isempty(cacheFile) && isfile(cacheFile)
    try
        C = load(cacheFile, 'key', 'R');
        if isequal(C.key, key)
            Rcached = C.R;
            for i = 1:min(numel(Rcached), n), R{i} = Rcached{i}; end
        else
            fprintf('  Surrogate cache ignored (settings or data changed).\n');
        end
    catch
        Rcached = {};
    end
end
todo = find(cellfun(@isempty, R));

fprintf('\n  Surrogate control (%d runs, each neuron circularly shifted independently):\n', n);
if ~isempty(todo)
    nT = numel(todo);
    fprintf('  Computing %d new surrogate(s); %d taken from the cache.\n', nT, n - nT);

    %% 3. Maximize Parallel 'Processes' Pool (unlock all CPU cores + prevent thread thrashing)
    usePool = false; nWorkers = 1; whyNot = '';
    nReq = -1;
    if isfield(cfg, 'DYN_PAR_WORKERS'), nReq = cfg.DYN_PAR_WORKERS; end
    if nReq == 0
        whyNot = 'cfg.DYN_PAR_WORKERS = 0';
    elseif ~license('test', 'Distrib_Computing_Toolbox')
        whyNot = 'no Parallel Computing Toolbox licence';
    elseif exist('gcp', 'file') ~= 2
        whyNot = 'Parallel Computing Toolbox not installed';
    else
        try
            % Detect physical cores and logical threads (via Java runtime)
            nPhys = feature('numcores');
            try
                nLogical = java.lang.Runtime.getRuntime().availableProcessors();
            catch
                nLogical = nPhys;
            end
            % If user didn't hardcode DYN_PAR_WORKERS, use all physical cores (or up to nLogical - 1)
            if nReq < 0
                nReq = max(nPhys, floor(0.75 * nLogical));
            end
            nReq = min(nReq, nT); % Never spawn more workers than remaining surrogates

            pool = gcp('nocreate');
            % If a ThreadPool or a smaller ProcessPool is open, restart with full nReq Processes workers
            if ~isempty(pool) && (isa(pool, 'parallel.ThreadPool') || pool.NumWorkers < nReq)
                delete(pool);
                pool = [];
            end
            if isempty(pool)
                clust = parcluster('Processes');
                clust.NumWorkers = nReq; % Override default profile cap
                pool = parpool(clust, nReq);
            end

            % Sync path and lock each worker to 1 internal thread so workers don't fight over CPU cache
            f1_sync = parfevalOnAll(pool, @path, 0, path);
            wait(f1_sync);
            f2_thr  = parfevalOnAll(pool, @maxNumCompThreads, 1, 1);
            wait(f2_thr);

            usePool = ~isempty(pool);
            if usePool, nWorkers = pool.NumWorkers; end
        catch MEpool
            whyNot = MEpool.message;
        end
    end

    if usePool
        fprintf('  Running on %d parallel Processes workers (100%% exact bin resolution).\n', nWorkers);
    else
        fprintf('  Running serially (%s).\n', whyNot);
    end

    %% 4. Broadcast Xb once to all workers & run with large chunks to eliminate barrier wait time
    tStart = tic;
    if usePool
        Xb_const = parallel.pool.Constant(Xb);
        chunk = max(3 * nWorkers, 12);
    else
        chunk = 1;
    end

    for c0 = 1:chunk:nT
        idx = c0:min(c0 + chunk - 1, nT);
        Rchunk = cell(1, numel(idx));
        todo_chunk = todo(idx);
        if usePool
            parfor kk = 1:numel(idx)
                Rchunk{kk} = fn_runOneSurrogate(Xb_const.Value, fs_b, t_range_b, cfg, ep_bounds, ep_labels, todo_chunk(kk));
            end
        else
            Rchunk{1} = fn_runOneSurrogate(Xb, fs_b, t_range_b, cfg, ep_bounds, ep_labels, todo_chunk(1));
        end
        for kk = 1:numel(idx), R{todo_chunk(kk)} = Rchunk{kk}; end
        fn_saveSurrCache(cacheFile, key, R, Rcached, n);
        el = toc(tStart);
        fprintf('    ... %d of %d done, %.1f min elapsed, about %.1f min left\n', ...
            idx(end), nT, el / 60, el / idx(end) * (nT - idx(end)) / 60);
    end
    fprintf('  Surrogates computed in %.1f min.\n', toc(tStart) / 60);
end

%% 5. Collect, print, p-values
for i = 1:n
    r = R{i};
    if isempty(r) || r.failed
        msg = ''; if ~isempty(r), msg = r.message; end
        warning('fn_dynamicalSurrogate:Failed', 'Surrogate %d failed: %s', i, msg);
        continue
    end
    surr.rot6(i) = r.rot6; surr.ndim80(i) = r.ndim80;
    if ~r.ok
        fprintf('    Surrogate %d: %d PCs for 80%% variance, rotation %.1f%%, no cycles found (%s)\n', ...
            i, r.ndim80, 100*r.rot6, r.message);
        continue
    end
    surr.n_cycles(i)   = r.n_cycles;
    surr.n_clusters(i) = r.n_clusters;
    surr.auc(i)        = r.auc;
    surr.contrast(i)   = r.contrast;
    surr.n_changes(i)  = r.n_changes;
    fprintf('    Surrogate %d: %d PCs for 80%% variance, rotation %.1f%%, %d cycles, %d clusters, within-between %.2f, AUC %.2f\n', ...
        i, r.ndim80, 100*r.rot6, r.n_cycles, r.n_clusters, r.contrast, r.auc);
end

surr.n_failed = sum(isnan(surr.n_cycles));
surr.p = struct('rot6', fn_pval(surr.rot6, surr.real.rot6, 'ge'), ...
                'ndim80', fn_pval(surr.ndim80, surr.real.ndim80, 'le'), ...
                'contrast', fn_pval(surr.contrast, surr.real.contrast, 'ge'));
if surr.n_failed > 0, fprintf('  %d of %d surrogates yielded no cycles.\n', surr.n_failed, n); end
fprintf('  Data vs surrogate mean: rotation (6 PCs) %.1f%% vs %.1f%% (p=%.2f) | PCs for 80%% variance %d vs %.0f (p=%.2f) | within-between %.2f vs %.2f (p=%.2f)\n', ...
    100*surr.real.rot6, 100*mean(surr.rot6, 'omitnan'), surr.p.rot6, ...
    surr.real.ndim80, mean(surr.ndim80, 'omitnan'), surr.p.ndim80, ...
    surr.real.contrast, mean(surr.contrast, 'omitnan'), surr.p.contrast);
end

function res = fn_runOneSurrogate(Xb, fs_b, t_range_b, cfg, ep_bounds, ep_labels, i)
% One surrogate: independent circular shift of every neuron, then the full pipeline.
c = cfg;
c.DYN_SURR_SHIFT = true;
c.DYN_SEED       = cfg.DYN_SEED + i;
c.DYN_NPERM      = 0;       % Unused inside surrogates (0 instead of 20 eliminates wasted inner loops)
c.DYN_NNULL      = 0;       % Unused inside surrogates
c.DYN_SKIP_ISROT = true;    % Extra jPCA in is_rotational is not used here
res = struct('i', i, 'ok', false, 'failed', false, 'message', '', 'rot6', NaN, 'ndim80', NaN, ...
    'n_cycles', NaN, 'n_clusters', NaN, 'auc', NaN, 'contrast', NaN, 'n_changes', NaN);
d = struct('ok', false, 'message', '');
try
    evalc('d = fn_dynamicalEpochs(Xb, fs_b, t_range_b, c, ep_bounds, ep_labels);');
catch ME
    res.failed = true; res.message = ME.message; return
end
if isfield(d, 'rot6'),   res.rot6   = d.rot6;   end
if isfield(d, 'ndim80'), res.ndim80 = d.ndim80; end
res.message = d.message;
if ~d.ok, return; end
res.ok         = true;
res.n_cycles   = numel(d.segment_s);
res.n_clusters = d.nGroups;
res.auc        = d.sim_thresh.auc;
res.contrast   = d.sim_thresh.within_mean - d.sim_thresh.between_mean;
res.n_changes  = numel(d.change_t);
end

function p = fn_pval(v, r, dirn)
nv = sum(~isnan(v));
if nv == 0 || isnan(r)
    p = NaN;
elseif strcmp(dirn, 'ge')
    p = (1 + sum(v >= r)) / (1 + nv);
else
    p = (1 + sum(v <= r)) / (1 + nv);
end
end

function fn_saveSurrCache(cacheFile, key, R, Rcached, n)
if isempty(cacheFile), return; end
Rsave = R;
if numel(Rcached) > n, Rsave(n+1:numel(Rcached)) = Rcached(n+1:end); end
try
    S.key = key; S.R = Rsave; %#ok<STRNU>
    save(cacheFile, '-struct', 'S');
catch MEsave
    warning('fn_dynamicalSurrogate:CacheSave', 'Could not save the surrogate cache: %s', MEsave.message);
end
end
