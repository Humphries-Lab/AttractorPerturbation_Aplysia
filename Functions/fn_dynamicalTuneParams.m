function tune = fn_dynamicalTuneParams(spike_conv, fs, cfg, t_range, ep_bounds, ep_labels, t_post, t_period, target)
% FN_DYNAMICALTUNEPARAMS  Data-driven choice of DYN_BIN_MS, DYN_MIN_SEG_S and DYN_REC_THRESH
% so that fn_dynamicalEpochs finds >= `target` cycles after the perturbation WITHOUT
% manufacturing them by over-splitting.
%
%   tune = FN_DYNAMICALTUNEPARAMS(spike_conv, fs, cfg, t_range, ep_bounds, ep_labels, t_post, t_period, target)
%
%   t_range   [t0 t1] window that the final fn_dynamicalEpochs call will analyse (dyn_range)
%   t_post    [t_C2 t_end]                  window in which >= target cycles are required
%   t_period  [t_recovery_start t_end]      window used to estimate the oscillation period
%             (default = t_post; prefer the recovery epoch so the evoked transient does not bias it)
%   target    minimum number of cycles after perturbation (default 10)
%
% STEPS
%   1. Dominant period T from the top 6 PCs of t_period (windowed FFT peak, cross-checked by
%      the autocorrelation). Cycle count is physically capped at diff(t_post)/T.
%   2. DYN_MIN_SEG_S = T/2 (a "cycle" shorter than half a period is a fragment);
%      DYN_BIN_MS    = T/40 (about 40 bins per cycle, never below 2 frames).
%   3. Re-runs fn_dynamicalEpochs over a grid of recurrence thresholds with those settings and
%      scores each by: cycles after perturbation >= target, cycle durations clustered around T,
%      few fragments (<T/2) and few stretches (>2T).
%
% Nothing is changed in cfg: copy the printed values into Config_UserSettings.m.

if nargin < 9 || isempty(target), target = 10; end
if nargin < 8 || isempty(t_period), t_period = t_post; end
[nN, nFrames] = size(spike_conv);

%% 1. Dominant period
bf = max(round(fs * 0.05), 1);                 % ~50 ms bins for spectral estimation
dt = bf / fs;
f1 = max(round(t_period(1) * fs) + 1, 1);
f2 = min(round(t_period(2) * fs), nFrames);
nb = floor((f2 - f1 + 1) / bf);
if nb < 40, error('fn_dynamicalTuneParams:ShortWindow', 't_period is too short to estimate a period.'); end
seg = double(spike_conv(:, f1 : f1 + nb*bf - 1));
X = reshape(mean(reshape(seg, nN, bf, nb), 2), nN, nb)';
[~, score] = pca(X - mean(X, 1));
K = min(6, size(score, 2));
S = detrend(score(:, 1:K));

w    = 0.5 - 0.5*cos(2*pi*(0:nb-1)'/(nb-1));
nfft = 2^nextpow2(4*nb);
F    = abs(fft(S .* w, nfft)).^2;
P    = sum(F(1:nfft/2+1, :), 2);
fr   = (0:nfft/2)' / (nfft * dt);
band = fr >= 2/(nb*dt) & fr <= 1/(8*dt);       % at least 2 cycles in window, at least 8 bins/cycle
Pb = P; Pb(~band) = 0;
[Ppk, ipk] = max(Pb);
T_fft = 1 / fr(ipk);
lo = ipk; while lo > 1 && Pb(lo-1) >= Ppk/2 && band(lo-1), lo = lo - 1; end
hi = ipk; while hi < numel(Pb) && Pb(hi+1) >= Ppk/2 && band(hi+1), hi = hi + 1; end
T_range = [1/fr(hi), 1/fr(max(lo, 2))];
peak_ratio = Ppk / median(P(band));

% autocorrelation cross-check (first peak after the first zero crossing)
Sc = S - mean(S, 1);
ac = real(ifft(sum(abs(fft(Sc, 2*nb)).^2, 2)));
ac = ac(1:floor(nb/2)); ac = ac / ac(1);
T_acf = NaN; acf_strength = NaN;
z = find(ac < 0, 1);
if ~isempty(z) && z + 2 <= numel(ac)
    loc = find(ac(z+1:end-1) > ac(z:end-2) & ac(z+1:end-1) >= ac(z+2:end), 1) + z;
    if ~isempty(loc), T_acf = (loc - 1) * dt; acf_strength = ac(loc); end
end
agree = ~isnan(T_acf) && abs(T_fft - T_acf) / T_fft < 0.25;
T = T_fft;

fprintf('\nTune: dominant period T = %.2f s (FFT; half-power %.2f-%.2f s; peak/median %.1f) | autocorr %.2f s (r=%.2f) | agree: %d\n', ...
    T, T_range(1), T_range(2), peak_ratio, T_acf, acf_strength, agree);
max_post = diff(t_post) / T;
max_per  = diff(t_period) / T;
fprintf('Tune: window after perturbation holds at most %.1f cycles (%.1f in the period-estimation window).\n', max_post, max_per);
if max_post < target
    warning('fn_dynamicalTuneParams:WindowTooShort', ...
        ['Only ~%.1f periods fit after the perturbation, fewer than the %d requested. No threshold can fix this ' ...
         'without splitting real cycles; a longer recovery or a faster oscillation is needed.'], max_post, target);
end

%% 2. Recommended bin / minimum segment
frame_ms = 1000 / fs;
bin_ms  = frame_ms * max(2, round(1000 * T / 40 / frame_ms));
min_seg = max(round(0.5 * T, 1), 2 * bin_ms / 1000);
fprintf('Tune: recommended DYN_BIN_MS = %.1f | DYN_MIN_SEG_S = %.1f (current %g | %g)\n', ...
    bin_ms, min_seg, cfg.DYN_BIN_MS, cfg.DYN_MIN_SEG_S);

%% 3. Recurrence-threshold sweep at those settings
if isfield(cfg, 'DYN_REC_SWEEP') && ~isempty(cfg.DYN_REC_SWEEP)
    grid = cfg.DYN_REC_SWEEP(:)';
else
    grid = cfg.DYN_REC_THRESH * [0.5 0.75 1 1.5 2];
end
grid = unique([grid, cfg.DYN_REC_THRESH]);

c = cfg; c.DYN_BIN_MS = bin_ms; c.DYN_MIN_SEG_S = min_seg;
c.DYN_NPERM = 100; c.DYN_NNULL = 100;          % nulls are not needed for tuning
iR = find(strcmp(ep_labels, 'Recovery'), 1);
n = numel(grid); nanv = nan(1, n);
res = struct('thresh', grid, 'n_post', nanv, 'n_rec', nanv, 'med_dur', nanv, ...
    'frac_short', nanv, 'frac_long', nanv, 'n_clusters', nanv, 'pass', false(1, n));

fprintf('  thresh | post-C2 | recovery | median dur (s) | <T/2 | >2T | clusters | pass\n');
for k = 1:n
    c.DYN_REC_THRESH = grid(k);
    d = struct('ok', false);
    try
        evalc('d = fn_dynamicalEpochs(spike_conv, fs, t_range, c, ep_bounds, ep_labels);');
    catch ME
        fprintf('  %6g | failed: %s\n', grid(k), ME.message); continue
    end
    if ~d.ok, fprintf('  %6g | no cycles\n', grid(k)); continue; end
    dur  = d.cycle_t(:,2) - d.cycle_t(:,1);
    post = d.cycle_mid >= t_post(1);
    res.n_post(k) = sum(post);
    if ~isempty(iR), res.n_rec(k) = sum(d.epoch_of_cycle == iR); end
    res.n_clusters(k) = d.nGroups;
    if any(post)
        res.med_dur(k)    = median(dur(post));
        res.frac_short(k) = mean(dur(post) < 0.5*T);
        res.frac_long(k)  = mean(dur(post) > 2*T);
        res.pass(k) = res.n_post(k) >= target && res.frac_short(k) <= 0.1 && res.frac_long(k) <= 0.1 ...
            && res.med_dur(k) >= 0.7*T && res.med_dur(k) <= 1.4*T;
    end
    fprintf('  %6g | %7g | %8g | %14.2f | %4.2f | %3.2f | %8g | %d\n', grid(k), res.n_post(k), res.n_rec(k), ...
        res.med_dur(k), res.frac_short(k), res.frac_long(k), res.n_clusters(k), res.pass(k));
end

best = NaN;
if any(res.pass)
    cand = find(res.pass);
    [~, j] = min(abs(res.med_dur(cand) - T));
    best = grid(cand(j));
    fprintf('Tune: use DYN_REC_THRESH = %g (>= %d cycles after perturbation, durations consistent with T).\n', best, target);
else
    fprintf(['Tune: no threshold gave >= %d cycles with durations consistent with T. Do not lower the threshold further ' ...
             'to force the count - check the window, T and the jPCA plane first.\n'], target);
end

tune = struct('T_fft', T_fft, 'T_range', T_range, 'T_acf', T_acf, 'periods_agree', agree, ...
    'max_cycles_post', max_post, 'bin_ms', bin_ms, 'min_seg_s', min_seg, ...
    'rec_thresh', best, 'sweep', res, 'target', target);
end
