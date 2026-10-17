# Blood-demand forecasting reanalysis

Reproduction materials for *Test-set-guided model selection in blood-demand forecasting*.
Original study: [Twumasi and Twumasi (2022)](https://doi.org/10.1016/j.ijforecast.2021.10.008).
The input manifest pins the authors' source files and SHA-256 digests. Original inputs are downloaded separately.

## Reproduce

From this repository's root, using Python 3.9+ and R:

```sh
python scripts/fetch_inputs.py
python scripts/verify_results.py
Rscript --vanilla scripts/run.R .
python scripts/verify_results.py
Rscript --vanilla scripts/plot_results.R .
```

The first verification recalculates stored results without fitting models. The R run refits all configurations and overwrites results; it can take several hours. Install the R packages `forecast`, `tsfknn`, `imputeTS`, `tsoutliers`, `jsonlite`, and `zoo` first. The recorded environment and dependency versions are in `results/session-info.txt`; matching versions are recommended for numerical comparison.

## Scope

The release contains 53 configurations, 901 origin-level errors and 9,010 forecast rows, covering 18 distinct test months. G/P denote global/prefix preprocessing; O/V denote test-guided/historical parameter selection. Historical horizons are matched, 6, 12 or 18 months; NNAR seeds are 1–5.

The original NNAR comparison disables two unused interval simulations. Nine stored checks support point-forecast equivalence; optionally regenerate them with `Rscript --vanilla scripts/check_intervals.R .`. This is a point-forecast reanalysis, not a reproduction of prediction intervals. The index-51 preprocessing fallback is recorded in the diagnostics.

`reference/` contains published NNAR errors and earlier reconstruction outputs used for cross-checks. `SHA256SUMS.json` identifies the release files; rerunning models may change output hashes across software environments.
