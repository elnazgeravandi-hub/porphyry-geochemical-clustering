# Porphyry Geochemical Clustering

Code for **Hybrid Clustering and Multi-Objective Evolutionary Optimization (MOPSO–NSGA-II) for Porphyry geochemical modeling and Mineral Prospectivity Enhancement**.

Version 0.3.0 implements the supplied manuscript Methods, parameter tables and Figure 5 workflow. The censoring algorithm was independently rewritten. This is newly written method code; numerical equivalence to the manuscript results has not yet been established using the restricted original inputs.

## Data restrictions

The authors do not have permission to publish or redistribute the original geochemical dataset and restricted spatial data. These are not included. No original notebook, stored notebook outputs, real sample rows or real occurrence coordinates are distributed. Local inputs and derived outputs must remain outside the public repository. See DATA_AVAILABILITY.md.

## Workflow

1. Input checks and explicit documented censoring flags/bounds.
2. ALR random-forest regression and Mills-ratio upper-truncation correction; six-part clr; common median/IQR scaling.
3. K-means, DBSCAN and OPTICS baselines plus MOPSO and NSGA-II optimization (nine configurations).
4. Fixed search subset; full-data Pareto reassessment and shared normalized distance-to-ideal selection.
5. SI, DBI, CHI, assigned-only DBCV, cluster sizes, coverage and noise reporting.
6. Elemental profiles, oriented PCA, continuous C-N regression and breakpoint sensitivity.
7. Verified sample-catchment links, geometric union, actual target area and unique occurrence capture at 5%, 10%, 20% ceilings.
8. Separate sampling stability and fixed-configuration preprocessing sensitivity, with optional separate reoptimization.

GIS catchment delineation from the DEM and final cartographic design remain external GIS tasks. No undocumented fusion or final-map selection rule is invented.

## Install

Use Python 3.12 and an isolated environment:

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install -e ".[report,spatial,dbcv,test]"
python -m pytest tests -q
```

`requirements.txt` pins core libraries. `pyproject.toml` declares optional reporting, spatial, DBCV and test dependencies. Run metadata records actual installed versions, including spatial dependencies. These are the new implementation's versions, not a claim about a historical environment.

## Local execution

The expected CSV columns are `sample_id,Cu,Au,Mo,Bi,Pb,Zn`; Au is ppb and the other elements are ppm. Use the private input formats in docs/input_schema.md.

```bash
python -m porphyry.preprocess --input private_data/raw_numeric.csv --censoring private_data/censoring.json --settings configs/censoring_numerics.json --output results/prepared_001 --validate-artificial
porphyry --input results/prepared_001/processed_concentrations.csv --config configs/manuscript.json --output results/models_001
python -m porphyry.reporting --input results/prepared_001/processed_concentrations.csv --labels-dir results/models_001 --output results/report_001 --spatial-config private_data/spatial.json
python -m porphyry.stability --input results/prepared_001/processed_concentrations.csv --models results/models_001 --output results/stability_001
```

Each output directory must be new. Use **processed concentrations**, never scaled-clr output, as input to clustering/reporting. Saved model preprocessing is held fixed for sampling stability.

## Scientific verification and additional choices

The manuscript does not fully specify residual-uncertainty estimation, censoring initialization, numerical floors or detailed random-number sequencing. The implementation makes these explicit in configs/censoring_numerics.json and docs/implementation_decisions.md. In particular, censoring correction uses OOB residual RMSE as an explicit additional implementation choice. This choice requires review or disclosure in Methods before this code is described as its exact computational realization.

Automated tests use artificial fixtures and reduced optimizer budgets; no manuscript numbers are fitted or hardcoded. Execute the full configuration on authorized original inputs and compare Tables 4/6/7 and Figures 6–14 before making an exact-reproduction claim. See docs/manuscript_mapping.md.
