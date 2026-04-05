---
name: open-phasea-ad-benchmark-repro
description: Reproduce v5 leakage-safe AD cross-cohort stress-test with consistent permutation-null reporting, AMP-AD feature ablations, and transductive ComBat sensitivity.
allowed-tools: Bash(python *), Bash(pip *), WebFetch
---

# Reproduction (v5)

## 0) Clone
```bash
git clone https://github.com/githubbermoon/bio-paper-track-open-phasea.git
cd bio-paper-track-open-phasea
git checkout main
```

## 1) Environment
```bash
python -m pip install --upgrade pip
python -m pip install numpy pandas scipy scikit-learn pycombat
```

## 2) Run benchmark
```bash
python src/train/run_open_phaseA_benchmark.py
```
Includes:
- feature modes: `var`, `de_ttest`, `agora_only`, `de_agora_intersection`
- primary arms: `target_only`, `source_only`, `source_plus_target_raw`
- sensitivity arm: `source_plus_target_combat_transductive`
- null outputs:
  - `null_label_permutation_mean_auroc` (primary reported null)
  - `null_label_permutation_avg100_prob` (probability-averaged sensitivity output)

## 3) Run stats
```bash
python src/eval/compute_open_phaseA_bootstrap.py
```
Produces bootstrap CIs and paired tests with BH correction.

## 4) Run local sensitivity checks
```bash
python src/eval/model_family_sensitivity.py
python src/eval/null_stability_check.py
```
Produces model-family robustness table and DE-1000 1000-permutation null stability summary.

## 5) AMP-AD open subset ingest
```bash
python src/ingest/fetch_ampad_open_subset.py
```

## 6) Expected artifacts
- outputs/metrics/open_phaseA_main_results.csv
- outputs/metrics/open_phaseA_predictions.csv
- outputs/stats/open_phaseA_null_distribution.csv
- outputs/stats/open_phaseA_auroc_ci.csv
- outputs/stats/open_phaseA_paired_tests.csv
- outputs/stats/open_phaseA_model_family_sensitivity.csv
- outputs/stats/open_phaseA_null_stability_de1000_perm1000.csv
- outputs/stats/open_phaseA_stats.json
- outputs/stats/open_phaseA_stats_manifest.json
- outputs/open_phaseA_data_manifest.json
- outputs/data/ampad_open_nominated_targets.csv
- outputs/tables/ampad_open_subset_summary.csv

## 7) Validation checks
```bash
python - <<'PY'
import json
from pathlib import Path
import pandas as pd

root = Path('.')
required = [
  'outputs/metrics/open_phaseA_main_results.csv',
  'outputs/metrics/open_phaseA_predictions.csv',
  'outputs/stats/open_phaseA_null_distribution.csv',
  'outputs/stats/open_phaseA_auroc_ci.csv',
  'outputs/stats/open_phaseA_paired_tests.csv',
  'outputs/stats/open_phaseA_model_family_sensitivity.csv',
  'outputs/stats/open_phaseA_null_stability_de1000_perm1000.csv',
  'outputs/stats/open_phaseA_stats.json',
  'outputs/open_phaseA_data_manifest.json',
  'outputs/data/ampad_open_nominated_targets.csv',
]
for f in required:
    assert (root / f).exists(), f'MISSING: {f}'

main = pd.read_csv(root/'outputs/metrics/open_phaseA_main_results.csv')
assert set(['var','de_ttest','agora_only','de_agora_intersection']).issubset(set(main['feature_mode'].unique()))
assert 'source_plus_target_combat_transductive' in set(main['arm'])
assert 'null_label_permutation_mean_auroc' in set(main['arm'])

stats = json.loads((root/'outputs/stats/open_phaseA_stats.json').read_text())
means = [v['null_perm_auroc_mean'] for v in stats.values()]
assert min(means) > 0.45 and max(means) < 0.55, means

paired = pd.read_csv(root/'outputs/stats/open_phaseA_paired_tests.csv')
assert 'target_only_vs_null_perm_mean' in set(paired['comparison'])
assert 'bh_adjusted_p' in paired.columns

manifest = json.loads((root/'outputs/open_phaseA_data_manifest.json').read_text())
assert 'ComBat transductive sensitivity' in manifest['batch_harmonization']

print('VALIDATION_OK')
PY
```
