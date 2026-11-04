# Explicit implementation decisions

The manuscript specification is authoritative. Additional numerical choices are listed here rather than silently presented as settings already reported by the article.

## Censoring

Manuscript-specified: six elements; ALR reference chosen among positive fully uncensored components by smallest log-IQR; all other ALR coordinates as predictors; RF=100 trees, maximum depth eight, seed 42; Mills-ratio upper truncation; relative Frobenius convergence <1e-5 with 100-iteration cap; only flagged values replaced.

Additional implementation decisions:

| Decision | Version 0.3 choice |
|---|---|
| Initial censored concentration | Its documented upper bound |
| Iteration order | Synchronous; every target uses previous ALR matrix |
| RF uncertainty for correction | Target-wise out-of-bag residual RMSE |
| Scale floor | 1e-8 in ALR units |
| Minimum observed training responses | 30 |
| RF sampling | Bootstrap, all predictors per split, leaf minimum one, split minimum two |
| Reference tie | Alphabetical elemental column name |
| Numerical extreme tail | Mills-ratio asymptotic expansion; no arbitrary concentration floor |

The OOB scale is a single target-wise predictive-error estimate, not calibrated local conditional uncertainty. The independent implementation is not represented as unchanged upstream ReX-CoDA code. These choices need author review and, where relevant, Methods clarification before claiming exact manuscript-code equivalence.

Nonconvergence produces a diagnostic file and no completed processed-concentration result. Unclassified ordinary missing values are rejected. No source-record meaning is inferred from negative or repeated small values.

Artificial validation uses complete confirmed-observed rows, five 80:20 splits (42–46), training-percentile bounds 10/20/30%, training-only reference/model fitting, identical masks for ALR-RF, DL/2 and DL/sqrt(2), log RMSE and signed error. Structural diagnostics use unscaled clr correlation changes and PCA configurations with up to five components, unit Frobenius normalization and orthogonal Procrustes. The dimension/normalization convention is an explicit choice not specified in the manuscript. Test truth is used only for scoring/reconstruction diagnostics, never fitting or reference choice. Natural censored truth and simultaneous-censoring accuracy are not established by these tests.

## Clustering and optimization

Both optimizers share inputs, method-specific domains/decoding, feasibility, initial populations per K/seed, fixed silhouette subset, budgets, final evaluation and selection. K=3–8 searches each receive 416 requests per optimizer/seed. This is not 416 total across all K. Fixed density baselines are reference configurations, not equal-budget searches.

K-means baseline uses 30 explicit k-means++ starts per K/seed and retains the least-inertia result. Explicit restarts allow storage of its selected initial centroids for stability checks. This is the manuscript's 30-start baseline design, not a claim of bitwise equivalence to another implementation's n_init RNG sequence. Selected initial centroids, not a silently substituted final centroid matrix, are exported.

Radius bounds are empirical linear-interpolation quantiles of neighbors including self. Discrete decoding uses NumPy round-to-even on exact halves. MOPSO is synchronous with 15% velocity caps, clipped boundaries and -0.5 velocity damping; NSGA-II uses bounded SBX and polynomial mutation. Exact crowding ties use stable input order in survival and seeded random tie-breaking in tournaments. These fine-grained RNG/tie choices are documented numerical decisions.

Feasible unique partitions are reevaluated on all assigned data. Shared normalization pools both optimizers/all seeds, plus baseline candidates for K-means. Final selection is per configuration's own nondominated set, then coverage, minimum size and decoded lexicographic ties. Complementary CHI/DBCV do not enter optimization. Noise is excluded from assigned-only indices and remains explicit in coverage.

Internal quality with different K or assigned coverage does not alone establish optimizer superiority or geological significance.

## Interpretation and area evaluation

Profiles use per-element whole-dataset median normalization, with cluster median and IQR. This normalization/plot convention is explicit. PCA selection is computed from maximum absolute Spearman association to the interpretation-only Cu/Au/Mo rank score; PC4 is never hardcoded. C-N uses continuous hinge regressions, 1–3 segments, minimum 20 distinct thresholds/segment, BIC counting breakpoints, six-unit acceptance improvement; sensitivity repeats 10/30. Numerical differential-evolution search does not certify a global optimum. An epsilon-relative RSS floor prevents roundoff-driven extra segments on exact fits.

C-N Cu classes and clustering targets are ranked as complete groups. Exact equal group scores are merged without secondary within-group splitting. Noise stays in the common domain and is not independently targeted. Cumulative unions prevent overlap double counting; eligible unique points are counted once. The largest complete prefix within each ceiling is selected; no skipping to lower-ranked groups or splitting polygons. Profiles and occurrence data are not optimization features/objectives.

DEM hydrology and final GIS styling are external operations. There is no undocumented fusion of the nine maps into one final prospectivity score.
