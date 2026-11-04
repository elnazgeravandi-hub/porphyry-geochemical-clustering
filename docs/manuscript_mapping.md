# Flowchart and manuscript mapping

| Manuscript/flowchart stage | Source | Output / qualification |
|---|---|---|
| Quality checks and documented censoring | preprocess.py | ID/units/value checks, explicit private mask/bounds; geological/source verification is external |
| ALR-RF/Mills-ratio treatment | censoring.py | New independent algorithm; declared numerical choices require review |
| clr, median/IQR | core.prepare | Positive processed inputs; one common fitted transformation per comparison |
| Three baselines | cli.py / core.py | K-means 30 starts; density m=8 and Q95(d8) |
| MOPSO / NSGA-II | optimizers.py | Shared comparison within each clustering algorithm |
| Budget decision | Evaluator.requests / run audits | 416 requests including initialization and cache hits |
| Full-data selection | cli.choose / metrics | Pooled normalization and own nondominated sets |
| Baseline validation / Table 6 | validation_summary.csv | SI, DBI, CHI, DBCV, K, counts, coverage/noise; coverage fraction x100 for percent |
| Elemental profiles / Figure 12 | reporting.elemental_profiles | Computed medians/IQR and plots; no perturbation |
| PCA, Cu C-N / Figure 13 | references.py / reporting.py | Computed component, breakpoints, scores/classes, optional catchment maps |
| Occurrence capture / Table 7, Figure 14 | evaluation.capture / reporting.py | Computed actual areas and capture; no table-value transcription |
| Sampling stability | stability.py | 20 subsamples, 80%, seeds 42–61, fixed preprocessing/configurations |
| Treatment sensitivity | sensitivity.py | User-documented processed scenarios; fixed-config comparison before optional separate reoptimization |
| Final GIS maps / Figures 6–11 | Exported sample labels + external GIS | Catchment linkage and cartography; no undocumented final fusion rule |

Table 4 is a separate source-reported concentration summary and is not reconstructed from imputed output. Source/laboratory checks, actual input hashes and numerical matching to manuscript results remain necessary before an exact-reproduction claim. Tests are artificial and use reduced optimizer budgets, not the manuscript experiment.
