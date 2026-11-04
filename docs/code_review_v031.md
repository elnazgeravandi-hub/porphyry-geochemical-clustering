# Code review and corrections: version 0.3.1

Scope: local source deliverable version 0.3.0, with corrections validated as version 0.3.1. This review does not claim a byte-for-byte inspection of the current remote GitHub repository.

## Findings corrected

| Finding | Correction | Verification |
|---|---|---|
| No feasible optimized partition could cause density method outcomes to disappear from the summary | Report the fixed baseline, including constraint violation and metric limitations, and explicit infeasible status for each optimizer | All-noise regression fixture checks six outcomes across two density methods |
| Fixed density baseline was discarded when optimizer feasibility limits were violated | Keep fixed reference baseline and report its coverage/noise; constraints remain unchanged for optimized candidates | Same all-noise regression, exported baseline labels |
| Reports/stability could use a changed input file with the same IDs | Verify SHA-256 against the model run manifest before downstream work | Changed row order is rejected |
| Spatial assessment could silently pair positional scores with differently ordered catchments | Require exact aligned ID order and label lengths | Misaligned spatial fixture rejected |
| Clipping could leave a boundary-only, zero-area catchment | Reject zero-area clipped linked geometry | Explicit check, alongside spatial I/O tests |
| CLI numeric parsing differed from downstream round-trip parsing | Use float_precision=round_trip consistently | Exported labels regenerated from selected parameters |
| Nonfinite DBCV could look like a valid reported result | Record missing value with explicit DBCV status | Ordinary assigned-subset DBCV covered in integration tests |
| Repeated partitions incurred repeated expensive full SI calculations | Cache full metrics by canonical partition within each clustering method | Nine-configuration integration tests |
| Final first-front selection used expensive repeated general nondominated sorting | Add O(n log n) first-front algorithm for finite feasible 2D objectives | Exact comparison against general sorting, including ties |

The faster first-front routine is used only in final selection/export. Evolutionary survival still uses general feasibility-first front sorting. Optimization objectives, primary constraints, domains, seeds and evaluation budgets are unchanged. Six new full-budget checks run population 16 with 25 updates for each optimizer/algorithm pair and assert exactly 416 evaluation requests.

## Final automated result

30 passed, 0 failed, 0 skipped (15.75 s), using artificial fixtures. Includes clustering, full request budgets, censoring invariants, computed profiles/reference plots, spatial file I/O/maps, whole tied-group capture, and 20-subsample stability. Does not run the restricted manuscript dataset or prove the published table values.

## Scientific limitations that code tests cannot remove

- The censoring code is an independently written ALR/RF/Mills-ratio reimplementation, not certified upstream ReX-CoDA. OOB residual scale and other unspecified numerical choices remain explicitly documented and require scientific alignment with Methods.
- Actual full-budget outcomes on the original sample matrix, laboratory censoring flags/bounds, occurrence eligibility and manuscript numerical agreement remain unverified.
- Different K and coverage require qualified comparisons. Equal optimizer conditions within a method are not equal computational budgets between fixed baselines and optimized searches.
- Original data are not included in this patch. Exact downstream checks require the original local processed-input artifact, not a newly rewritten CSV with the same apparent values.
