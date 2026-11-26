# Aplysia Attractor Perturbation Analysis Toolbox
 
MATLAB toolbox for analysing population-level attractor dynamics after a
motor perturbation in *Aplysia* recordings. Couple of sister pipelines have been added which traces the motor activity without the perturbation and motor activity across P9 repeated stimuli for sensitization.
 
A single user-edited settings script, a folder of small reusable
functions (one function per file, each documenting the inputs/outputs and
which window of the time series it operates on), and a top-level pipeline
script that only calls those functions in sequence.
 
## Top-level scripts
 
1. **`Config_UserSettings.m`** - the *only* place parameters should be
   edited (file paths, sampling rate, epoch timing, window sizes,
   thresholds). Every other script/function reads its parameters from the
   `cfg` struct this script defines; nothing is hard-coded elsewhere.
2. **`Run_Attractor_Analysis.m`** - the main pipeline. Loads `cfg`, then
   calls the analysis and plotting functions in sequence (load/clean/smooth
   -> PCA -> windowed PCA & participation ratio -> subspace alignment ->
   recurrence-density analysis -> RR-based attractor epoch detection ->
   dynamical epochs from jPCA cycles (Section 6), with per-cycle metrics
   (6b) and the stimulus phase as a rotation angle in the jPC plane (6c)),
   saving all figures to `cfg.FIGURES_DIR` and a results `.mat` file to
   `cfg.RESULTS_DIR`. See "Dynamical epochs (jPCA cycles)" below for
   Section 6.
3. **`Run_Evoked_Analysis.m`** - companion pipeline for any single-stimulus
   recording without perturbation (Baseline + Evoked only, no Recovery epoch).
   Uses `cfg.t_P9_evoked`/`cfg.t_end_evoked` instead of `fn_getEpochTiming`,
   and altered `_evoked`-suffixed plotting wrappers (bottom of the script) calling
   upon the original `Functions/` to bypass `Recovery`/`C2`.
   Saves `<recording_ID>_results.mat` in the same format `Run_Cohort_Analysis.m`
   expects, but without any `_Recovery`/`_rr` columns. 
4. **`Run_Sensitization_Analysis.m`** - companion pipeline for a repeated-stimulus
   (sensitization) recording: five 5-minute files concatenated into one
   25-minute recording, each with one P9 at 90 s (file-relative), giving
   **1 Baseline + 5 Evoked epochs**. Protocol constants (`nStim`, `FILELEN_S`,
   `T_P9_INFILE_S`) are set at the top of the script; everything else comes
   from `Config_UserSettings.m`, and the same `Functions/` library is used
   unmodified. Runs the same stages as `Run_Evoked_Analysis.m`, but
   per stimulus, plus across-stimulus readouts: PR vs stimulus, subspace
   alignment (to Evoked_1, to Baseline, and step-to-step Evoked_i -> Evoked_i+1,
   each with the chance level), self-recurrence per stimulus, the full
   Evoked_i x Evoked_j cross-recurrence matrix (plus step-to-step i -> i+1
   transitions, shown on a fixed colour scale `RR_CLIMS` shared across
   recordings), and attractor onset latency per stimulus (threshold crossing on
   the fixed-epsilon recurrence trace). All plotting functions are local
   to the script (bottom of file); figures are saved as 500 DPI PNGs to
   `cfg.FIGURES_DIR/<recording_ID>/` and results to
   `<recording_ID>_results_sensitization.mat`. Incomplete trailing stimuli
   (recording shorter than 5 x 300 s) are dropped with a warning. The jPCA
   dynamical-epoch analysis (Section 6 of `Run_Attractor_Analysis.m`) is not
   part of this pipeline.
5. **`Run_Cohort_Analysis.m`** - run this AFTER `Run_Attractor_Analysis.m` or `Run_Evoked_Analysis.m`
   has been run for every animal/recording. Pools all
   `<recording_ID>_results.mat` files into one table for cohort-level
   statistics. Sensitization results are **not** pooled by this script.
6. **`Run_Cohort_Analysis_Sensitization.m`** - the dedicated cohort pipeline for repeated-stimulus data. Run this AFTER `Run_Sensitization_Analysis.m`. It pools `*_results_sensitization.mat` files into a wide-format table (one row per recording, with columns for each stimulus) optimized for SPSS/R repeated-measures ANOVA. Generates cohort trend plots tracking motor program identity and complexity across stimuli.
 
## Folders
 
- **`Functions/`** - one function per file, grouped by analysis stage:
  - `DataIO/` - loading the peaks matrix, deriving epoch timing from the protocol
  - `Preprocessing/` - neuron quality filtering, kernel-width estimation, Gaussian smoothing
  - `Dimensionality/` - PCA, windowed PCA, participation ratio (fixed & sliding)
  - `SubspaceAlignment/` - the alignment metric, K selection, epoch/K-sweep/sliding alignment
  - `Recurrence/` - recurrence density, cross-recurrence density, epsilon calibration, attractor epoch detection
  - `DynamicalEpochs/` - jPCA cycles and everything built on them (Section 6): `fn_dynamicalEpochs` (jPCA, cycle detection, cycle-by-cycle alignment, clustering, drift check, landmark tests), `fn_dynamicalCycleMetrics` (per-cycle PR, alignment, time-fraction stimulus phase), `fn_stimulusPhaseJPC` (stimulus phase as an angle in the jPC plane), `fn_dynamicalSweep` (recurrence-threshold sweep), `fn_dynamicalSurrogate` (shifted-neuron surrogate control), and their plotting functions `fn_plotDynamicalEpochs`, `fn_plotDynamicalThreshold`, `fn_plotDynamicalCycleMetrics`, `fn_plotDynamicalControls`
  - `Plotting/` - one function per figure; each takes only the data it needs to draw and a save path
  - `Utils/` - small generic helpers (pass/fail check, ternary, run-length detection, synthetic subspace construction)
- **`Validation/`** - synthetic validation ("unit test") scripts for the two
  key metrics, checked against known theoretical values. These are
  development/QA scripts, not part of the analysis pipeline itself (they
  use no real recording data):
  - `Validate_SubspaceAlignment.m`
  - `Validate_RecurrenceDensity.m`
- **`CohortAnalysis/`** - the functions that power the top-level cohort scripts.
  Sit outside `Functions/` because nothing in the single-recording pipelines calls them directly.
  - `fn_cohortAnalysis.m` - Pools single-stimulus/attractor `.mat` files.
  - `fn_cohortAnalysis_sensitization.m` - Pools sensitization files, extracts 3D matrices, and flattens them to wide-format Excel tables with dynamic missing-field padding.
- **`Results/`**, **`Figures/`** - default output locations (see `Config_UserSettings.m`).
 
## How do I use this with my own data?
 
1. Make sure your recording is saved as a `.mat` file containing a
   variable `peaks` (nNeurons x nFrames, logical or numeric event matrix).
2. Open `Config_UserSettings.m` and set `DATA_FILE`, `RESULTS_DIR`,
   `FIGURES_DIR`, `recording_ID`, `protocol` (`'12min'` or `'20min'`), and
   `fs` to match your recording. Adjust epoch timing, window sizes, and
   thresholds if needed (all fields are documented in that script).
3. Run `Run_Attractor_Analysis.m` (or `Run_Evoked_Analysis.m` for a single P9, or
   `Run_Sensitization_Analysis.m` for repeated P9 stimuli).
4. To add a new protocol duration, add a `case` to `fn_getEpochTiming.m`.
5. To add a new figure, write a new function in `Functions/Plotting/` that
   takes the data it needs plus a save path, and call it from
   `Run_Attractor_Analysis.m` - keep plotting code out of the analysis
   functions, and analysis code out of the plotting functions.
 
## Dynamical epochs (jPCA cycles): Section 6, 6b, 6c of `Run_Attractor_Analysis.m`
 
Sections 1-5 define epochs from the *protocol* (Baseline / Evoked / Perturbed /
Recovery) and compare them with epoch-level statistics. Section 6 asks the
complementary, protocol-free question: does the population's own cycle-by-cycle
dynamics change, and if so does the change line up with the perturbation? It
works on single cycles of the population trajectory, which also makes it
possible to separate a lasting change from slow drift.
 
### Switches (`Config_UserSettings.m`)
 
| Setting | Meaning |
|---|---|
| `cfg.DYN_ENABLE` | `true` runs Section 6 (and 6b, 6c); `false` skips them. |
| `cfg.DYN_WINDOW` | Analysis window: `'full'` = 0 to `t_end`; `'post_evoked'` = `t_evoked_start` to `t_end`. |
| `cfg.DYN_REC_SWEEP` | Recurrence thresholds (%) for the threshold sweep, e.g. `[5 10 15]`. Empty skips the sweep. |
| `cfg.DYN_NSURR` | Number of shifted-neuron surrogates (e.g. 99). `0` skips the surrogate control. |
 
Other `cfg.DYN_*` fields are documented in `Config_UserSettings.m`. The
protocol epochs used to label cycles are built in the script from the epoch
timing (`ep_bounds` / `ep_labels`: Baseline, Evoked, Perturbed, Recovery).
 
### Section 6 - `fn_dynamicalEpochs`
 
1. **Preprocessing.** `spike_conv` is binned (about 50 ms), mean-subtracted and
   reduced to the 6 leading PCs (the fixed dimension used throughout Section 6,
   matching the `rot6` control). The number of PCs needed for 80% of the
   variance is reported alongside.
2. **Rotational plane (jPCA).** Fit `dx/dt = M x` with `M` skew-symmetric.
   *Rotation strength* = fraction of variance in the first jPC plane.
3. **Cycles from recurrence.** A cycle closes where the trajectory recurs in
   the jPC plane (recurrence threshold 10%, minimum cycle length 10 s).
4. **Alignment and clustering.** For every pair of cycles,
   `S_ij = Var(X_i in the subspace of cycle j) / Var(X_i)`, symmetrised. Cycles
   are clustered on `S`; the clusters are the *dynamical epochs* and are compared
   with the protocol epochs (cycles per cluster in each protocol epoch, mean
   alignment between protocol epochs, whether Evoked and Recovery share a
   dominant cluster).
5. **Readouts and tests.**
   - *Change points* along the cycle sequence (times at which the cluster
     label changes), the nearest change to each protocol landmark, and a
     permutation test of whether change points lie closer to the landmarks than
     chance.
   - *Drift check:* a line `S = a + b * dt` is fitted to same-epoch cycle
     pairs; the expected Evoked<->Recovery alignment is the mean of `a + b * dt`
     over Evoked-Recovery pairs, compared with the observed alignment. The
     percentage of Evoked-Recovery pairs whose time separation is longer than any
     same-epoch pair (extrapolated) is reported.
   - *Strongest single split* of the cycle sequence into two time blocks
     (contrast = within-block minus between-block alignment), and the percentile
     of the contrast at each protocol landmark among all splits.
   - *Cycles after the perturbation:* the number of cycles after C2 (and in the
     Recovery epoch only) against a target of 10; reported as attained, never
     forced.
6. **Controls.** The recurrence-threshold sweep (`cfg.DYN_REC_SWEEP`) repeats
   the analysis at other thresholds and prints one summary table
   (`fn_dynamicalSweep`). The surrogate control (`fn_dynamicalSurrogate`) circularly
   shifts every neuron independently (`x_i(t) -> x_i((t + tau_i) mod T)`),
   reruns the analysis, and compares rotation strength, PCs for 80% variance and
   within-minus-between cluster alignment with the data using
   `p = (1 + #{surrogates >= data}) / (1 + N)`. With `N = 99` the smallest possible
   p is 0.01.
 
Figures: `fn_plotDynamicalEpochs` (cycle-by-cycle alignment, cluster timeline,
cycles per cluster in each protocol epoch, mean alignment between epochs),
`fn_plotDynamicalThreshold` (the cycle-closing threshold) and
`fn_plotDynamicalControls` (sweep and surrogate panels).
 
### Section 6b - `fn_dynamicalCycleMetrics`
 
Per-cycle readouts, called with the stimulus times `[t_P9 t_C2]`, their labels,
and `t_evoked_start` (the reference cycle is the first cycle after the start of
Evoked):
 
- participation ratio per cycle (`PR/N`, dimensionality over time),
- alignment of each cycle with the previous cycle and with the reference cycle
  (movement of the attractor),
- the phase of each stimulus as the **fraction of cycle time** elapsed,
  `(t_stim - t_start) / (t_end - t_start)`, of the cycle that contains it.
 
Because cycles are cut where the trajectory returns, phase 0 is not a fixed point
on the attractor, so this time fraction is *not* comparable between cycles. 6c
gives an angle instead. Figure: `fn_plotDynamicalCycleMetrics`.
 
### Section 6c - `fn_stimulusPhaseJPC`
 
The phase of each stimulus as a **rotation angle in the jPC plane**, so the same
angle means the same point on the cycle within one recording. It is a
self-contained re-fit and does not call the jPCA code:
 
1. Bin `spike_conv` as in `fn_dynamicalEpochs` (same bin, same window) and
   subtract the mean.
2. PCA to 6 PCs -> scores `Z`.
3. Fit `dZ/dt = M Z` with skew-symmetric `M` by least squares (all 15 free
   entries).
4. The real Schur form of `M` gives three orthonormal rotational planes; the
   plane capturing most variance is kept (the "first jPC plane"). Axes are signed
   so the rotation is counter-clockwise.
5. `theta(t) = unwrap(atan2(z2, z1))`. Per cycle: rotation swept (degrees)
   between its first and last bin. Per stimulus: absolute angle, angle-based
   fraction of the containing cycle, and the time-based fraction for comparison.
 
Output `ph` (also saved in the results file): `ok`, `plane_var_frac`,
`rot6_pipeline`, `t`, `theta`, `z`, `cycle_sweep_deg`, and `stim(k)` with
`name`, `t`, `cycle`, `phase_abs_rad`, `phase_abs_deg`, `phase_cycle_frac_angle`,
`phase_cycle_frac_time`, `cycle_sweep_deg`. Figure: `06e_stimulus_phase_jpc.png`
(trajectory of the cycle containing each stimulus, and the stimulus angles on a
polar plot).
 
### Saved variables
 
`dyn` (Section 6), `cyc` (6b), `ph` (6c), `dyn_sweep` (threshold sweep) and
`dyn_surr` (surrogates) are saved in `<recording_ID>_results.mat`.
 
### How to read the command-window log for Section 6
 
| Log line / variable | Meaning |
|---|---|
| `PCs needed for 80% of variance` | Dimensionality of the binned activity. Compare with the surrogates (shifted neurons need many more). |
| `jPCA (fixed 6 PCs): first jPC plane captures X%` | **Rotation strength**: variance in the first jPC plane of the 6-PC subspace (also saved as `dyn.rot6`). |
| `Cycles detected: n (mean s)` | Number and mean length of recurrence-defined cycles. |
| `Clusters (dynamical epochs) found` | Number of clusters of the cycle-by-cycle alignment matrix. |
| `Evoked vs Recovery: n vs m cycles` | Cycles falling in each protocol epoch. |
| `shared cluster (E in R), (R in E)` | Overlap of cluster membership between Evoked and Recovery in each direction. |
| `alignment E<->R / within E / within R` | Mean cycle-pair alignment `S` between Evoked and Recovery, within Evoked, and within Recovery. Within >> between = the epochs are distinct. |
| `same dominant cluster` | 1 if Evoked and Recovery are mostly in the same cluster, 0 otherwise. |
| `Change points (s)` | Times at which the cluster label changes along the cycle sequence. |
| `nearest change to <landmark>` | Signed time from each protocol landmark to the closest change point (in seconds, with the mean cycle length for scale). |
| `Change points vs landmarks: mean |lag| ... permutation p` | Mean distance of change points to the landmarks against random expectation, with a permutation p. A large p means change points are not locked to the landmarks. |
| `Similarity: within-cluster / between-cluster / chance / threshold / AUC` | Alignment of same-cluster and different-cluster cycle pairs, the 95th percentile of the chance level, the separating threshold and the AUC. These are circular (the clusters are built from the same matrix); use the surrogate comparison for inference. |
| `Drift check: within-epoch alignment changes by ... per min; E<->R observed vs expected from drift alone` | Slope of alignment with time separation inside epochs, and the Evoked<->Recovery alignment expected from that slope alone, with the extrapolated share of pairs. An observed value well below the expected one is evidence for a change beyond drift. |
| `Strongest single split ... contrast ... percentile among all splits` | The time point that best divides the cycle sequence into two blocks, and the contrast at each landmark relative to all possible splits. |
| `Mean alignment between protocol epochs` | Epoch-by-epoch mean alignment; the diagonal is within-epoch. `NaN` on the diagonal means the epoch has a single cycle. |
| `Baseline vs Recovery` | The same comparison for Baseline and Recovery. |
| `Cycles after perturbation (t >= C2): n (Recovery epoch only: m)` | Cycle counts against the target of 10; "of which k start after the perturbation" excludes the cycle that straddles it. |
| `Recurrence-threshold sweep` table | Cycles, clusters, overlaps and alignments at each threshold in `cfg.DYN_REC_SWEEP`. |
| `Data vs surrogate mean: rotation / PCs for 80% / within-between` | Data against the shifted-neuron surrogates with p-values (floor 1/(1+N)). |
| `Cycle metrics: ... reference = cycle r` | 6b: number of cycles, ranges of `PR/N` and of step alignment, and the reference cycle. |
| `P9/C2: cycle c, phase p (deg) \| PR/N \| align to previous, next-to-this` | 6b: the stimulus time-fraction phase and the PR and alignment around that cycle. |
| `Stimulus phase (jPC plane): plane variance ... median cycle sweep ...` | 6c: variance in the plane re-fitted here (close to, not identical with, `rot6`; see below) and the typical rotation swept by a cycle. |
| `angle in plane ... of the rotation (... of the time) ... cycle sweeps ...` | 6c: stimulus angle, angle-based and time-based cycle fractions, and how many degrees the containing cycle sweeps. |
 
## Cohort-level analysis (across recordings)
 
### Standard & Evoked Cohorts
Once you have run `Run_Attractor_Analysis.m` or `Run_Evoked_Analysis.m` for every animal, run **`Run_Cohort_Analysis.m`** to pool the results. Edit `recording_IDs` near the top first to list exactly which recordings belong to be analysed in a cohort. This is deliberate rather than auto-pooling everything in that folder, thus it keeps a record of exactly what went into each cohort table. 
 
It calls `fn_cohortAnalysis(RESULTS_DIR, recording_IDs)` and pulls a fixed set of scalar metrics into one row-per-recording table, saving it as `CohortAnalysis_Summary.mat` and an Excel sheet which can be used further for statistics. 
 
### Sensitization Cohorts
To analyse repeated-measures sensitization data, run **`Run_Cohort_Analysis_Sensitization.m`**. It automatically pools all `*_results_sensitization.mat` files in your chosen directory and extracts parameters across all 5 stimuli. 
 
It generates a wide-format Excel table (`Cohort_Summary_Sensitization.xlsx`) specifically structured for repeated-measures ANOVA, where each row is one animal and each column is a metric/stimulus permutation (e.g., `PR_evoked_S1`, `PR_evoked_S2`).
 
### What each variable means
 
**Sensitization Variables (Motor Identity vs. Complexity)**
In sensitization experiments, it is necessary to separate the *strength* of the motor response from the *identity* of the motor program. The following metrics are extracted across all stimuli (S1-S5):
 
| Variable | Meaning |
|---|---|
| `PR_evoked_S#` | **Complexity:** Normalized participation ratio for Stimulus #. A drop indicates the network is collapsing into lower-dimensional space (fewer neurons firing in complex patterns). |
| `RR_S#` | **Attractor Density:** Self-recurrence density for Stimulus #. An increase indicates the network is locking onto a specific trajectory more tightly. |
| `Onset_Time_S#` | **Convergence Speed:** Time (s) it takes for the network to cross the recurrence density threshold for Stimulus #. |
| `Align_to_S1_S#` | **Linear Identity:** Subspace alignment of Stimulus # relative to the PCA subspace of Stimulus 1. If PR drops but alignment remains high, the animal is executing the exact same motor pattern, just less vigorously. A drop indicates recruitment of a functionally different motor program. |
| `Cross_RR_S#` | **Non-linear Identity:** Cross-recurrence density of Stimulus # against Stimulus 1. Similar to alignment, this tests whether the specific non-linear state space trajectory matches the initial response. |
 
**Standard Attractor Variables (Evoked vs Recovery)**
*(Note: Refer to fixed vs RR-defined epoch distinctions based on the pipeline run)*
 
| Variable | Meaning |
|---|---|
| `align_ER` | Raw alignment score (0-1) between the Evoked and Recovery population subspaces. Higher = the two epochs occupy more similar directions in neural state space. |
| `chance_lvl` | The alignment score expected by chance for two *randomly oriented* subspaces of the same dimensionality (`= K/N`). |
| `align_ER_corrected` | `align_ER` rescaled against `chance_lvl`: `(align_ER - chance_lvl) / (1 - chance_lvl)`. **Use this one, not `align_ER`, when comparing across animals**. |
| `PR_norm_...` | Normalized participation ratio in the specified epoch. Higher = higher-dimensional activity; lower = collapsed onto dominant directions. |
| `RQA_ev_RR` / `_re_RR` | Self-recurrence density within the specified epoch using an epsilon calibrated on that epoch's own trajectory. |
| `cross_recur_density`| How often Evoked-epoch states recur *within the Recovery epoch*, using one shared epsilon (calibrated on Evoked). |
| `t_attractor_onset` | Time (s) the population first detectably locked onto a recurring state after the stimulus. |
| `onset_detected` | `true` if lock-on was genuinely detected; `false` means detection failed and timing columns are fallbacks. |
 
**Dynamical-epoch variables (jPCA cycles; stored in `dyn`, `cyc`, `ph`)**
 
| Variable | Meaning |
|---|---|
| `dyn.rot6` | Rotation strength: fraction of variance in the first jPC plane of the 6-PC subspace. |
| `dyn.cycle_t`, `dyn.cycle_mid`, `dyn.startbin` | Start/end times (s), midpoints and start bins of the recurrence-defined cycles. |
| `dyn.epoch_of_cycle`, `dyn.ep_labels` | Protocol-epoch label of each cycle and the label names. |
| `cyc` | Per-cycle `PR/N`, alignment with the previous and reference cycle, and the time-fraction phase of each stimulus. |
| `ph.plane_var_frac` | Variance in the plane fitted by `fn_stimulusPhaseJPC` over the total variance (compare with `dyn.rot6`). |
| `ph.cycle_sweep_deg` | Rotation swept by each cycle in the jPC plane. A clean rotation sweeps close to 360 degrees. |
| `ph.stim(k).phase_abs_deg` | Absolute angle of stimulus `k` in the jPC plane (comparable between stimuli within one recording only). |
| `ph.stim(k).phase_cycle_frac_angle` / `_time` | Fraction of the containing cycle by swept angle / by elapsed time. |
| `dyn_sweep`, `dyn_surr` | Recurrence-threshold sweep and surrogate-control results. |
 
## Things to be aware of
 
- All functions assume `spike_conv` (Gaussian-smoothed activity) is
  `neurons x frames`, single precision, and that PCA bases (`V`, `mu`) are
  computed once and reused everywhere a shared reference subspace is
  needed.
- In `Run_Attractor_Analysis.m`, self-recurrence density for an epoch is
  always calibrated with that epoch's *own* epsilon (`fn_safeEpsilon`);
  cross-recurrence density between two epochs intentionally uses one shared
  epsilon (from the reference epoch).
- `fn_recurrenceDensity` and `fn_crossRecurrenceDensity` are both
  O(N^2) in the number of trajectory points; always downsample long
  trajectories first (see `fn_getEpochTrajectory`, and `cfg.MAX_WIN_PTS` /
  `cfg.MAX_FULL_PTS`).
- In `Run_Sensitization_Analysis.m`, the PCA basis is fit from the start of
  Evoked_1 to the end of the recording, and the alignment dimensionality/chance
  level are calibrated on Evoked_1, so alignment comparisons are relative to that
  first response (its self-alignment is 1 by construction). The recurrence epsilon
  is one fixed percentile of pairwise distances pooled across all P9-to-P9
  intervals, identical for every stimulus. Sliding-window metrics and the epsilon
  pooling intervals span file joins, so changes at file boundaries can be artifacts.
- **Section 6 (jPCA) caveats.**
  - *Rotation can be weak.* Rotation strength is a fraction of the variance of a
    6-PC subspace; check `ph.cycle_sweep_deg` (and `ph.plane_var_frac`) before
    interpreting angles. If cycles sweep far less than 360 degrees in the jPC
    plane, treat the angle as a coarse phase only.
  - *`phase_cycle_frac_angle` is only meaningful when the containing cycle sweeps
    roughly a full rotation.* For cycles that barely rotate it can be unstable or
    fall outside [0, 1]; ignore those values and use `phase_abs_deg` or the time
    fraction instead.
  - *Angles are comparable within a recording only.* The jPC axes are arbitrary
    and the plane is fitted on the whole analysis window (`cfg.DYN_WINDOW`), so it
    is dominated by the epoch with the most variance.
  - *`plane_var_frac` is close to, not identical with, `dyn.rot6`.* The jPCA code
    in the toolbox normalises the data and `fn_stimulusPhaseJPC` does not. A large
    mismatch means "check", not "wrong".
  - *The time-fraction phase (6b) is relative to where each cycle was cut* and is
    not comparable across cycles; use the angle (6c) for that.
  - *Clusters and drift.* Separation of Evoked and Recovery cycles can be produced
    by slow drift alone. Always read the drift-check line (observed vs expected
    alignment and the extrapolated share of pairs) together with the
    surrogate comparison, not the cluster contrast on its own. The same-cluster vs
    different-cluster contrast and the AUC are circular because the clusters are
    built from the same alignment matrix.
  - *Landmark locking.* Change points that are not closer to the protocol
    landmarks than chance (large permutation p) mean the cluster boundaries do
    not track the stimuli, even when Evoked and Recovery fall in different
    clusters.
  - *Number of cycles.* Fewer than 10 cycles after C2 is reported as attained, not
    forced; compare the threshold sweep before trusting borderline results.
  - *Surrogate p-values.* With 99 surrogates the smallest attainable p is 0.01,
    so `p = 0.01` means "no surrogate reached the data".
  - *Sensitization recordings* do not run Section 6.
 
## Problems?
 
Check the docstring at the top of the relevant function first - every
function in `Functions/` documents its inputs, outputs, and (where
relevant) which time window it expects. If something still looks wrong,
re-run the matching script in `Validation/` to confirm the underlying
metric is behaving as expected before suspecting the pipeline wiring.
 
