---
name: open-phasea-ad-benchmark-repro
description: Reproduce the final leakage-safe AD cross-cohort benchmark (v7 packaging) with consistent permutation-null reporting, AMP-AD Agora feature ablations, and transductive ComBat sensitivity analysis.
allowed-tools: Bash(python *), Bash(pip *), WebFetch
---

# Reproduction (final v7 package)

## 0) Clone
```bash
git clone https://github.com/githubbermoon/bio-paper-track-open-phasea.git
cd bio-paper-track-open-phasea
git checkout main
```

## 1) Environment (frozen)
Use the exact frozen environment from this repo:
```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## 2) Fresh regeneration sequence (exact order)
Important: run AMP-AD open fetch before training, because train script uses the Agora CSV.
```bash
python src/ingest/fetch_ampad_open_subset.py
python src/train/run_open_phaseA_benchmark.py
python src/eval/compute_open_phaseA_bootstrap.py
python src/eval/model_family_sensitivity.py
python src/eval/null_stability_check.py
```

## 3) What the benchmark includes
- feature modes: `var`, `de_ttest`, `agora_only`, `de_agora_intersection`
- primary leakage-safe arms: `target_only`, `source_only`, `source_plus_target_raw`
- sensitivity-only arm: `source_plus_target_combat_transductive`
- primary null in tables: `null_label_permutation_mean_auroc`
- additional null sensitivity output: `null_label_permutation_avg100_prob`

## 4) Expected artifacts
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

## 5) Validation checks
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

main = pd.read_csv(root / 'outputs/metrics/open_phaseA_main_results.csv')
assert set(['var','de_ttest','agora_only','de_agora_intersection']).issubset(set(main['feature_mode'].unique()))
assert 'source_plus_target_combat_transductive' in set(main['arm'])
assert 'null_label_permutation_mean_auroc' in set(main['arm'])

stats = json.loads((root / 'outputs/stats/open_phaseA_stats.json').read_text())
means = [v['null_perm_auroc_mean'] for v in stats.values()]
assert min(means) > 0.45 and max(means) < 0.55, means

paired = pd.read_csv(root / 'outputs/stats/open_phaseA_paired_tests.csv')
assert 'target_only_vs_null_perm_mean' in set(paired['comparison'])
assert 'bh_adjusted_p' in paired.columns

manifest = json.loads((root / 'outputs/open_phaseA_data_manifest.json').read_text())
assert 'ComBat transductive sensitivity' in manifest['batch_harmonization']
assert 'generated_outputs' in manifest and len(manifest['generated_outputs']) >= 7

print('VALIDATION_OK')
PY
```
