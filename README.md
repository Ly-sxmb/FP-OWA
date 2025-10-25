# FP-OWA
Pseudocode

A.Overview
Algorithm: FunctionalRankingPipeline
Input:
    n_samples     // number of samples (e.g., 100)
    t_length      // length of the time grid (e.g., 30)
    noise_type     // contamination type: gaussian/spike/amplitude/... (see Module B2)
    noise_params    // parameters for the chosen noise (e.g., sigma)
    n_runs       // number of simulation runs (e.g., 500)
    parallel_cores   // number of parallel cores

Output:
    Mean, standard error, and 95% CI of evaluation metrics comparing the three methods (FP-OWA, WLR, FPCA)

Steps:
    1) Generate functional data with trend and periodic components (Module B1), obtaining:
        – clean_fd: noise-free fd object
        – true_integrals: ground-truth integrals
        – true_rank: ground-truth ranks
    2) Add contamination to clean_fd according to the specified noise type (Module B2) to obtain polluted_fd.
    3) Run the three ranking methods (Module C):
        – FP-OWA(polluted_fd)
        – wlr(polluted_fd)
        – fpca(polluted_fd)
    4) Use the evaluation module (Module D) to compare the ranking performance of each method.
    5) Perform batch simulations (Module E): repeat steps 2–4 and record metrics and runtime.
    6) Aggregate results and compute 95% confidence intervals (Module F).

End

B Data Module
B1. Algorithm: GenerateFunctionalData
Input:
    n, t_length, noise_sd (default = 0, used only to generate “raw observations” noise)

Output:
    clean_fd     // fd object represented by B-spline basis (no additional smoothing)
    t_points     // equally spaced time points [0, 1]
    basis       // B-spline basis
    true_integrals  // integrals of the “noise-free true patterns”
    true_rank     // ranks of true_integrals (averaged ranks for ties)

Steps:

Construct an equally spaced time sequence t ← linspace(0, 1, t_length).

Define a B-spline basis (range [0,1], nbasis = 8, norder = 4).

For each i = 1..n:
    – Randomly generate a trend coefficient slope ∈ U[3.5, 4.5].
    – Randomly generate a phase shift phase ∈ U[0, 0.2].
    – Define the true curve: true_i(t) ← slope * t + sin(2π (t + phase)).
    – Stack all true_i(t) to form true_patterns.

Generate raw observations: raw_data ← true_patterns + N(0, noise_sd) (add pointwise Gaussian noise).

Convert true_patterns into an fd object true_patterns_fd.

Compute the true integrals:
    – Prepare a constant basis representing the function 1.
    – For each curve true_i, compute ∫ true_i(t) dt = inprod(true_i, 1).
    – Obtain true_integrals.

Compute true_rank based on true_integrals (using average ranking for ties).

Convert raw_data into an fd object clean_fd using Data2fd (without additional smoothing).

Return: (clean_fd, t_points, basis, true_integrals, true_rank)
End

B2. Algorithm: AddPollution
Input:
    fd_obj, type ∈ {gaussian, spike, amplitude, poisson, uniform, laplace, exponential}
    ... // parameters corresponding to the selected type (e.g., sigma for gaussian)

Output:
    polluted_fd

Steps (process the coefficient matrix coefs of fd_obj based on the specified type):
    – gaussian:   coefs ← coefs + N(0, sigma)
    – spike:     randomly select ratio × number of sample columns, add a fixed large value (e.g., +5) to each selected column
    – amplitude:  randomly select ratio × number of sample columns, multiply each selected column by a fixed amplification factor (e.g., ×2)
    – poisson:   coefs ← coefs + Pois(lambda)
    – uniform:   coefs ← coefs + U(lower, upper)
    – laplace:   coefs ← coefs + Laplace(mu, sigma)
    – exponential: coefs ← coefs + Exp(rate)

Return:
    An fd object reconstructed using the same basis.
End

C.Ranking Methods
C1. Algorithm: FPCA_
Input:
    fd_obj

Output:
    (scores, ranks, first_component, n_components_used, explained_var_head)

Steps:

Extract the sample × time matrix: X ← transpose(fd_obj.coefs).

Center the data by subtracting the mean across samples for each time point: Xc ← row-wise subtract mean over samples.

Compute the covariance matrix: Σ ← cov(Xc).

Perform eigen-decomposition Σ = V Λ Vᵀ and calculate the proportion of explained variance (explained_var).

Select the smallest number of principal components k such that the cumulative explained variance ≥ 80% (ensure k ≥ 1).

Compute principal component scores: S ← Xc · V[:, 1..k].

For each curve, compute its mean value over time: mᵢ ← mean of X[i] across columns.

Compute the sign of the correlation between the first principal component and the mean values: corr_sign ← sign(cor(S[:,1], m)).
  If corr_sign < 0, reverse the direction of S[:,1].

Define the ranking scores: scores ← S[:,1]; ranks ← rank(scores).

Return: scores, ranks, first principal component vector, number of components used, and the first five explained variance ratios.
End

C2.Algorithm: WLR
Input:
    fd_obj

Output:
    (scores, ranks, seg_points, weights, segment_ranks)

Steps:

Extract the sample × time matrix: X ← transpose(fd_obj.coefs).

Find intersection points (segment boundaries):
    For all pairs of curves (i, j):
        Compute diff = X[i,·] − X[j,·], and identify positions where the sign of diff changes.
    Collect all such positions → intersections.

Form segment boundaries: seg_points ← {1, intersections, T_len} (remove duplicates and sort).

Compute segment weights based on relative length:
    weights ← (segment_length / total_length) for each segment.

For each segment s:
    – Extract the corresponding submatrix of X for this segment.
    – Compute the mean value of each curve within the segment.
    – Determine within-segment ranks: segment_ranks[, s] ← rank(segment_means, ties = average).
    – Accumulate total scores: total_scores += weights[s] × segment_ranks[, s].

Compute final overall ranks: ranks ← rank(total_scores, ties = average).

Return: total_scores, ranks, seg_points, weights, segment_ranks
End

C3.Algorithm: FP-OWA_Modified
Input:
    fd_obj, t_length_for_eval, lambda_candidates (e.g., 10^seq(-6, 0))

Output:
    (scores, seg_points, chosen_lambda)

Steps:

// (1) Spline smoothing and λ selection via GCV

Construct an evaluation grid t_grid over the range of fd_obj.

Evaluate fd_obj on t_grid to obtain raw_mat.

For each λ in lambda_candidates:
    – Apply smooth.basis to raw_mat using fdPar(basis, Lfd=2, λ).
    – Record the corresponding GCV(λ) value (take the mean if it is a vector).
    Select the λ* that minimizes GCV; set smooth_fd ← smoothed result corresponding to λ*.

// (2) Automatic segmentation via DBSCAN (clustering on “time-point sample vectors”)
4) Construct timepoint_data ← transpose(raw_mat)  // each row represents one time point with n-dimensional sample values.
5) Standardize timepoint_data.
6) Determine eps using the k-distance plot:
    k_distances ← kth-neighbor-distance(timepoint_data, k = 4)
    eps ← elbow_point(k_distances)
    minPts ← 2 × dimension(timepoint_data)  // heuristic rule.
7) Apply DBSCAN(timepoint_data, eps, minPts) to obtain labels.
    Replace noise points (label = 0) with the nearest nonzero neighboring label (or default to 1 if none).

Generate segment boundaries where labels change:
    seg_points ← {t_grid[1], change points, t_grid[end]}.

// (3) Segment weights = normalized mean MBD (Modified Band Depth) within each segment
9) For each segment s in seg_points:
    – seg_fd ← WindowFD(smooth_fd, seg_start, seg_end).
    – mbd_vals ← MBD(seg_fd)  (handle errors by assigning zeros).
    – w_s ← mean(mbd_vals); if all zeros, assign equal weights.
    Normalize all {w_s} so that their sum equals 1, forming the weight vector mbd_weights.

// (4) Within-segment integral ranking and weighted aggregation
10) For each segment s:
    – seg_fd ← WindowFD(smooth_fd, seg_start, seg_end).
    – Compute the basis integrals over [a, b]:
        diag_integrals ← diag(inprod(basis, basis, rng = [a, b])).
    – For each curve j: integral_j ← Σ_k (coef_{k,j} × diag_integrals_k).
    – Compute within-segment ranks: segment_ranks[, s] ← rank(integral_j, ties = average).

Compute final ranking scores: scores ← segment_ranks · mbd_weights.

Return: scores, seg_points, λ*
End

D.Algorithm: EvaluateRanking
Input: true_ranks, pred_scores, true_integrals (optional)
Output: metrics = {kendall_tau, spearman_rho, mard, top_consistency, bottom_consistency,
                  integral_rmse (if true_integrals), integral_mae (if true_integrals)}
Steps:
    1) pred_ranks ← rank(pred_scores, ties=average)
    2) 
         - kendall_tau ← Kendall(true_ranks, pred_scores)
         - spearman_rho ← Spearman(true_ranks, pred_scores)
    3) 
         - mard ← mean(|true_ranks - pred_ranks|)

Return: metrics
End

E.Algorithm: RunSimulation
Input:
    n_runs, sample_size, noise_type, noise_params, t_length, parallel_cores

Output:
    results_table  // rows = runs × methods; columns = various evaluation metrics and runtime

Steps:

Initialize the parallel computing environment using the specified number of cores (parallel_cores).

For i in 1..n_runs (executed in parallel):
    a) Generate synthetic data: data ← GenerateFunctionalData(sample_size, t_length).
    b) Add noise contamination: polluted_fd ← AddPollution(data.clean_fd, noise_type, noise_params).
    c) Record runtime and execute the three methods:
        FP-OWA_res ← EFSRWM_Modified(polluted_fd, t_length)
        wlr_res   ← WLR(polluted_fd)
        fpca_res  ← FPCA_Ranking_Improved(polluted_fd)
    d) Compute evaluation metrics for each method:
        FP-OWA_metrics ← EvaluateRanking(data.true_rank, FP-OWA_res.scores, data.true_integrals)
        wlr_metrics  ← ...
        fpca_metrics  ← ...
    e) Assemble three result rows (method ∈ {FP-OWA, WLR, FPCA}), including their runtime.

Combine all result rows from every run into a single results_table.

Return: results_table
End

F. Algorithm: SummarizeWithCI
Input:
    results_table  // contains method, evaluation metrics, and runtime

Output:
    final_summary  // mean, standard error (SE), and 95% confidence interval (CI) for each method and metric

Steps:

Group results by method and calculate the mean of each metric (including runtime).

Group results by method and calculate the standard error (SE) for each metric: SE = standard deviation / sqrt(valid sample size).

Compute the 95% confidence interval using z = 1.96: CI = mean ± 1.96 × SE.

Generate a summary table that includes the mean, SE, and 95% CI for all metrics and methods.

Return: final_summary
End
