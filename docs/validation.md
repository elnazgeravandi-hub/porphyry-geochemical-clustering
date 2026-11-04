# Development validation

Automated checks use ephemeral artificial numerical fixtures, not manuscript data.

- clr/unit conversion equivalence and rejection of nonpositive values.
- Label-permutation-invariant partition hashing and separate noise handling.
- MOPSO and NSGA-II request accounting and feasible centroid partitions.
- End-to-end execution across all nine configurations with explicitly reduced test budgets.
- Shared-normalization selection and coverage tie-breaking.
- PCA orientation and enrichment tie handling.
- C-N rejection of spurious extra segments on an exact power-law fixture.
- Spatial whole-tie/union test exists but is skipped without the optional Shapely dependency.

Initial development test result: 8 passed, 1 skipped. Spatial integration, optional DBCV, full manuscript budgets and agreement with original numerical results are not validated by these tests.

A numerical RSS floor at machine epsilon times max(centered response sum of squares, 1) prevents BIC favoring extra segments based solely on floating-point differences near exact fits. This is an explicit numerical implementation detail.
