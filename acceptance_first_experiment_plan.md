# Acceptance-first experiment plan

Task: predict growth/non-growth at h=2/h=4 weeks for (state,lineage,time) units.

Must include:
- split arms: random, lineage-grouped, strict temporal-group (primary)
- baselines: heuristics, logistic regression, LightGBM/XGBoost
- null controls: label/feature permutation
- metrics: AUROC, AUPRC, balanced accuracy, Brier, ECE + 95% CI
- significance: paired bootstrap + BH correction
- ablations + failure analysis + reproducibility bundle
