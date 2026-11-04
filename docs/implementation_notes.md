# Implementation decisions and limitations

## Shared comparison

Within each clustering method, MOPSO and NSGA-II receive matching processed inputs, search domains, decoding, constraint violation, random seeds, initial populations, request budgets and final selection. Initial K-means populations use repeated k-means++ draws; density starting populations use uniform coordinates in the declared log-radius domain. The latter distribution and detailed RNG sequencing are implementation choices, not recoverable original-run evidence.

Nearest neighbors explicitly query the training array and therefore include self. Quantiles use NumPy's default linear interpolation. Nearest-integer decoding uses NumPy round-to-even at exact halves; OPTICS index values map to 8,10,12,15,20. Label-invariant canonical hashes deduplicate partitions while retaining noise separately.

MOPSO updates synchronously, uses zero starting velocities, 15% velocity limits, clipped boundaries with -0.5 velocity damping, a 50-member feasible archive, crowding tournaments and sequential crowding recomputation on truncation. NSGA-II uses bounded SBX (pair probability .9; coordinate participation .5; eta 15), polynomial mutation probability 1/D and eta 20. Equal crowding values in survival use stable input order. These detailed choices must be checked against any original implementation if it is recovered.

The 416 budget is per optimizer, per seed, and per K for K-means. K-means therefore receives six K searches; do not present 416 as its total across K. Density baselines are fixed references, not equal-budget searches. Final selections with different K and coverage require qualified interpretation; internal indices alone do not demonstrate geological superiority.

## Evaluation

Noise is excluded from assigned-only metrics. Search SI is computed within the same fixed subset, intersected with assigned observations for each candidate; it is not an estimate using distances to the full dataset. Full SI and DBI are recalculated for each unique feasible partition. Normalization is pooled across both optimizers and all seeds; K-means additionally includes feasible baseline configurations. Per-optimizer selection considers its pooled nondominated candidates across K.

Infeasible candidates and runs are logged without relaxing constraints. JSON null denotes an invalid/nonfinite audit value; feasibility and status fields identify the reason. DBCV requires the optional hdbscan package, explicitly enabled in configuration; its installed version must be recorded before final analysis.

C-N uses a numerical differential-evolution breakpoint search with continuous hinge regressions and minimum distinct-threshold counts. BIC counts slope/intercept, hinge coefficients and breakpoint locations; improvement must be at least six to accept the next segment. Global optimum is not certified. Run the documented 10/20/30-threshold sensitivity before interpreting thresholds. This is not asserted identical to an unknown earlier segmentation implementation.

## Outstanding integration

ReX-CoDA reference: Emenaha et al. (2026), DOI 10.1007/s11004-026-10331-8; author repository linked by the article: https://github.com/ogonna-emenaha/ReX-CoDA . This development package does not substitute another imputer under that name. Residual uncertainty estimation, training/prediction protocol and back-transformation preserving measured concentrations need confirmation against the original code.

Additional pending work: artificial censoring validation and preprocessing scenarios; complete plotting commands; projected spatial I/O and eligibility audit; actual manuscript-output comparison; full-budget execution on authorized data. No claimed result is hardcoded from the manuscript tables.
