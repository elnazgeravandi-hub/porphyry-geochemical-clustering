# Verification scope

Automated verification uses artificial numerical fixtures, not the restricted manuscript dataset.

Covered checks:

- clr unit-conversion/centering equivalence and invalid input rejection.
- Label-permutation-invariant partitions, noise and ARI.
- Request budgets and feasible optimizer results at reduced test budgets.
- All nine configurations, DBCV, and regeneration of exported partitions from stored selected parameters/initial centroids.
- Shared normalization and tie rules.
- Oriented PCA, enrichment ties and C-N power-law model selection.
- Profile median/IQR normalization and computed plots.
- Spatial file loading, projected coordinates, sample-catchment links, complete tied group selection, union areas, occurrence deduplication and computed maps.
- Independent upper-truncated mean implementation, observed-value preservation, invariance to hidden censored responses, nonconvergence diagnostics, explicit primary limits and native-unit preservation.
- Training-only artificial-censoring comparisons and clr/PCA-Procrustes identity diagnostics.
- Twenty fixed-parameter 80% subsamples.

The final test count/result is recorded in the delivered version's validation_result.txt. Results do not establish original-data agreement, laboratory coding validity, truth of natural censored observations, superiority of imputation or optimizers, or full-budget manuscript reproduction. Fixed-configuration sensitivity orchestration requires scientifically defined local scenarios; it is not asserted executed on the original data.
