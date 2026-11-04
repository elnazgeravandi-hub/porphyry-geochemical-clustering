# Input contract

The CLI accepts a local CSV of **positive processed concentrations**, after verified censored-value treatment. It does not accept raw laboratory codes or infer censoring from small numerical values.

| Column | Requirement |
|---|---|
| sample_id | Unique nonmissing string; preserve order |
| Cu, Mo, Bi, Pb, Zn | Positive finite numeric concentrations, ppm |
| Au | Positive finite numeric concentration, ppb |

Columns are selected explicitly in the order Cu, Au, Mo, Bi, Pb, Zn. Other columns are not predictors. No catchment concentration averaging, API feature or geographic coordinate enters clustering.

Au is converted to ppm before clr. Fixed per-column unit conversion introduces a constant clr offset which cancels under column median centering; the equivalence is tested. Zero IQR is an explicit error, rather than an undocumented scale substitution.

For eventual imputation integration, retain separate source-derived censoring flags and element-specific bounds (Au 5 ppb, Bi 0.1 ppm as stated by the author). Do not reconstruct these flags from already imputed values. Verified unflagged measurements must remain unchanged.

Each sample may be linked by sample_id to an upstream catchment. Catchments may overlap. Spatial assessment requires a documented projected area CRS, study boundary and eligible unique occurrence inventory. These are not required by the clustering CLI.
