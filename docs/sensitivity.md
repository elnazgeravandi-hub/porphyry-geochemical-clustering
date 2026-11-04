# Sampling and treatment sensitivity

Sampling stability:

```bash
python -m porphyry.stability --input results/prepared_001/processed_concentrations.csv --models results/models_001 --output results/stability_001
```

Twenty 80% subsamples without replacement (seeds 42–61). Preprocessing parameters are loaded from the full-data run and not refitted. Every selected configuration is refitted, including saved initial K-means centroids. ARI is reported with noise as a label and on common assigned observations, together with common counts/fractions, K, cluster sizes and coverage. Fewer than two represented clusters yields an uninformative ARI (empty value), not a fabricated score.

Treatment sensitivity requires a private JSON with a scenarios array. Each item has name (unique identifier), processed_input (private processed concentration CSV), scientific_basis (actual rationale) and scaling_mode (fixed_primary or refit_common). No undocumented alternative threshold is generated.

```bash
python -m porphyry.sensitivity --primary-input results/prepared_001/processed_concentrations.csv --models results/models_001 --scenarios private_data/scenarios.json --output results/sensitivity_001
```

Add --reoptimize-config configs/manuscript.json to perform separate full-budget optimization after fixed-configuration comparisons. Reoptimization uses newly fitted scenario preprocessing; this remains separate from the fixed stage's declared scaling mode. The choice of scaling mode is not fully specified by the manuscript and must be stated in reporting. Sample IDs must match exactly.

Alternative feasibility limits are not generated or searched automatically; scientifically specified constraint-sensitivity experiments require explicit configuration/code support before being claimed. Do not label an unrun sensitivity experiment completed.
