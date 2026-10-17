# Blood-demand forecasting reanalysis

Reproduction materials for *Test-set-guided model selection in blood-demand forecasting*.
Original study: [Twumasi and Twumasi (2022)](https://doi.org/10.1016/j.ijforecast.2021.10.008).
The input manifest pins the authors' source files and SHA-256 digests. Original inputs are downloaded separately.

## Reproduce

From this repository's root, using Python 3.9+ and R:

```sh
python scripts/fetch_inputs.py
python scripts/verify_results.py
python scripts/verify_extensions.py
```

To refit the models and regenerate the additional checks:

```sh
Rscript --vanilla scripts/run.R .
python scripts/fetch_grnn.py
Rscript --vanilla scripts/run_window_sensitivity.R .
Rscript --vanilla scripts/run_grnn.R .
Rscript --vanilla scripts/check_information_flow.R .
python scripts/verify_results.py
python scripts/verify_extensions.py
Rscript --vanilla scripts/plot_results.R .
```

The first verification recalculates stored results without fitting models. The R run refits all configurations and overwrites results; it can take several hours. Install the R packages `forecast`, `tsfknn`, `imputeTS`, `tsoutliers`, `jsonlite`, and `zoo` first. The recorded environment and dependency versions are in `results/session-info.txt`; matching versions are recommended for numerical comparison.

## Scope

The release contains 53 configurations, 901 origin-level errors and 9,010 forecast rows, covering 18 distinct test months. G/P denote global/prefix preprocessing; O/V denote test-guided/historical parameter selection. Historical horizons are matched, 6, 12 or 18 months; NNAR seeds are 1–5.

Version 1.1.0 retains those results and adds six window-count/spacing designs (408 metric rows), GRNN comparisons (153 rows), and an audit of all seven original forward implementations. The window comparisons share feasible candidate sets across designs. GRNN source version 0.1.0 reproduces the original mixed MIMO-selection/recursive-forecasting calls; consistent MIMO and recursive comparisons are reported separately. The historical package version is unrecorded. LSTM checks cover input availability only, without fitting a network.

The original NNAR comparison disables two unused interval simulations. Nine stored checks support point-forecast equivalence; optionally regenerate them with `Rscript --vanilla scripts/check_intervals.R .`. This is a point-forecast reanalysis, not a reproduction of prediction intervals. The index-51 preprocessing fallback is recorded in the diagnostics.

`reference/` contains the code audit, original NNAR/GRNN errors and reconstruction cross-checks. GRNN errors were extracted from the pinned `Main_Script_IJF.html`. `SHA256SUMS.json` identifies release files; reruns may change hashes across software environments. Third-party sources retain their original terms; the historical GRNN source is GPL-2.
