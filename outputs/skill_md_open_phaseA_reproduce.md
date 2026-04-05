---
name: open-phasea-ad-benchmark-repro
description: Reproduce v4 leakage-safe AD cross-cohort transfer stress-test with DE feature selection, ComBat harmonization arm, bootstrap/BH statistics, and AMP-AD open subset summary.
allowed-tools: Bash(python *), Bash(pip *), WebFetch
---

# Reproduction (v4)

## 0) Clone repository
```bash
git clone https://github.com/githubbermoon/bio-paper-track-open-phasea.git
cd bio-paper-track-open-phasea
git checkout main
```

## 1) Environment
- Python 3.11+

Install dependencies:
```bash
python -m pip install --upgrade pip
python -m pip install numpy pandas scipy scikit-learn pycombat
```

## 2) Run benchmark
```bash
python src/train/run_open_phaseA_benchmark.py
```
This now includes:
- feature modes: `var`, `de_ttest`
- transfer arms: `target_only`, `source_only`, `source_plus_target_raw`, `source_plus_target_combat`
- null strategy: `null_label_permutation_avg100` (100 permuted-label models averaged)

## 3) Run statistics
```bash
python src/eval/compute_open_phaseA_bootstrap.py
```
Computes bootstrap CIs and paired tests with BH correction for:
- `source_plus_target_combat_vs_target_only`
- `target_only_vs_null_label_permutation_avg100`
- `source_only_vs_target_only`

## 4) Ingest AMP-AD open subset
```bash
python src/ingest/fetch_ampad_open_subset.py
```

## 5) Expected outputs
- outputs/metrics/open_phaseA_main_results.csv
- outputs/metrics/open_phaseA_predictions.csv
- outputs/stats/open_phaseA_auroc_ci.csv
- outputs/stats/open_phaseA_paired_tests.csv
- outputs/stats/open_phaseA_stats.json
- outputs/stats/open_phaseA_stats_manifest.json
- outputs/open_phaseA_data_manifest.json
- outputs/data/ampad_open_nominated_targets.csv
- outputs/tables/ampad_open_subset_summary.csv

## 6) Validation checks
```bash
python - <<'PY'
import json
from pathlib import Path
import pandas as pd

root = Path('.')
required = [
  'outputs/metrics/open_phaseA_main_results.csv',
  'outputs/metrics/open_phaseA_predictions.csv',
  'outputs/stats/open_phaseA_auroc_ci.csv',
  'outputs/stats/open_phaseA_paired_tests.csv',
  'outputs/stats/open_phaseA_stats.json',
  'outputs/open_phaseA_data_manifest.json',
  'outputs/data/ampad_open_nominated_targets.csv',
  'outputs/tables/ampad_open_subset_summary.csv',
]
for f in required:
    assert (root / f).exists(), f'MISSING: {f}'

main = pd.read_csv(root/'outputs/metrics/open_phaseA_main_results.csv')
assert 'feature_mode' in main.columns
assert set(['var','de_ttest']).issubset(set(main['feature_mode'].unique()))
assert 'source_plus_target_combat' in set(main['arm'])
assert 'null_label_permutation_avg100' in set(main['arm'])

paired = pd.read_csv(root/'outputs/stats/open_phaseA_paired_tests.csv')
assert 'bh_adjusted_p' in paired.columns
assert 'source_plus_target_combat_vs_target_only' in set(paired['comparison'])

stats = json.loads((root/'outputs/stats/open_phaseA_stats.json').read_text())
# Null AUROC should center near 0.5 across settings
means = [v['perm_auroc_mean'] for v in stats.values()]
assert min(means) > 0.45 and max(means) < 0.55, means

manifest = json.loads((root/'outputs/open_phaseA_data_manifest.json').read_text())
assert manifest['batch_harmonization'].startswith('ComBat')
assert 'GSE63060' in manifest['urls'] and 'GSE63061' in manifest['urls']

amp = pd.read_csv(root/'outputs/tables/ampad_open_subset_summary.csv')
kv = dict(zip(amp['metric'], amp['value']))
assert int(kv['input_data::RNA']) > 0
assert int(kv['input_data::Protein']) > 0
assert int(kv['input_data::Genetics']) > 0

print('VALIDATION_OK')
PY
```

Expected terminal output: `VALIDATION_OK`
