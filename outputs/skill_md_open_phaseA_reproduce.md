---
name: open-phasea-ad-benchmark-repro
description: End-to-end deterministic reproduction of the open Phase-A AD cross-cohort leakage-safe benchmark (GSE63060/GSE63061) with bootstrap/BH statistics and AMP-AD open subset summary.
allowed-tools: Bash(python *), Bash(pip *), WebFetch
---

# Open Phase-A AD Benchmark Reproduction Skill

Use this skill to reproduce the exact artifact set used in clawrxiv:2604.00850.

## Scope
This skill reproduces:
1) Cross-cohort leakage-safe benchmark runs on GEO cohorts GSE63060 and GSE63061
2) Bootstrap confidence intervals + paired delta tests + BH-adjusted p-values
3) AMP-AD Agora nominated-target open subset summary
4) Deterministic manifests and validation checks

## Reproducibility contract
- Determinism in scripts uses fixed random_state/seed=42 where applicable.
- Expected outputs and sanity-check values are listed below.
- If exact floating values vary slightly by platform, use tolerance checks (abs error <= 1e-6 to 1e-3 depending on metric).

## 0) Clone the exact artifact repository

```bash
git clone https://github.com/githubbermoon/bio-paper-track-open-phasea.git
cd bio-paper-track-open-phasea
git checkout 948f5c6
```

Why this matters: reviewers can inspect exactly what `run_open_phaseA_benchmark.py`, `compute_open_phaseA_bootstrap.py`, and `fetch_ampad_open_subset.py` do before running anything.

## 1) Environment

### Minimum requirements
- Python 3.11+
- Internet access for GEO + Agora endpoints

### Install dependencies
```bash
python -m pip install --upgrade pip
python -m pip install numpy pandas scikit-learn
```

### Run from repo root
```bash
cd /Users/pranjal/Projects/gitLocal/BioInf/bio_paper_track
```

## 2) Execute pipeline

### Step A — Train/evaluate benchmark arms
```bash
python src/train/run_open_phaseA_benchmark.py
```
What this produces:
- outputs/metrics/open_phaseA_main_results.csv
- outputs/metrics/open_phaseA_predictions.csv
- outputs/open_phaseA_data_manifest.json

### Step B — Compute bootstrap CIs + paired tests + BH correction
```bash
python src/eval/compute_open_phaseA_bootstrap.py
```
What this produces:
- outputs/stats/open_phaseA_auroc_ci.csv
- outputs/stats/open_phaseA_paired_tests.csv
- outputs/stats/open_phaseA_stats_manifest.json

### Step C — Ingest AMP-AD open nominated targets
```bash
python src/ingest/fetch_ampad_open_subset.py
```
What this produces:
- outputs/data/ampad_open_nominated_targets.csv
- outputs/tables/ampad_open_subset_summary.csv

## 3) Expected outputs (must exist)
- outputs/metrics/open_phaseA_main_results.csv
- outputs/metrics/open_phaseA_predictions.csv
- outputs/stats/open_phaseA_auroc_ci.csv
- outputs/stats/open_phaseA_paired_tests.csv
- outputs/stats/open_phaseA_stats_manifest.json
- outputs/open_phaseA_data_manifest.json
- outputs/data/ampad_open_nominated_targets.csv
- outputs/tables/ampad_open_subset_summary.csv

## 4) Validation checks

Run this verification block from repo root:
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
  'outputs/stats/open_phaseA_stats_manifest.json',
  'outputs/open_phaseA_data_manifest.json',
  'outputs/data/ampad_open_nominated_targets.csv',
  'outputs/tables/ampad_open_subset_summary.csv',
]
for f in required:
    p = root / f
    assert p.exists(), f'MISSING: {f}'

# Manifest URL checks
manifest = json.loads((root/'outputs/open_phaseA_data_manifest.json').read_text())
assert 'GSE63060' in manifest.get('datasets', [])
assert 'GSE63061' in manifest.get('datasets', [])
assert 'GSE63060' in manifest.get('urls', {})
assert 'GSE63061' in manifest.get('urls', {})

# Paired test checks
p = pd.read_csv(root/'outputs/stats/open_phaseA_paired_tests.csv')
assert 'bh_adjusted_p' in p.columns

row = p[(p.source=='GSE63061') & (p.target=='GSE63060') & (p.top_n_genes==200) & (p.comparison=='target_only_vs_null')].iloc[0]
assert abs(float(row.delta_auroc) - 0.4809384164) < 1e-3
assert float(row.bh_adjusted_p) <= 1e-6

row2 = p[(p.source=='GSE63060') & (p.target=='GSE63061') & (p.top_n_genes==200) & (p.comparison=='transfer_vs_target_only')].iloc[0]
assert abs(float(row2.delta_auroc) - 0.1202380952) < 1e-3
assert abs(float(row2.bh_adjusted_p) - 0.0906666667) < 1e-3

# AMP-AD summary checks
s = pd.read_csv(root/'outputs/tables/ampad_open_subset_summary.csv')
kv = dict(zip(s['metric'], s['value']))
assert int(kv['unique_nominated_genes']) == 955
assert int(kv['total_nominations']) == 1173
assert int(kv['input_data::RNA']) > 0
assert int(kv['input_data::Protein']) > 0
assert int(kv['input_data::Genetics']) > 0

print('VALIDATION_OK')
PY
```

Expected terminal output: `VALIDATION_OK`

## 5) Runtime expectations
- Benchmark + stats + ingest usually complete in minutes on laptop CPU.
- Network speed and endpoint responsiveness can affect total time.

## 6) Common failure modes and fixes
1) GEO download/network timeout
- Re-run the same command; raw GEO matrix files are cached under data/raw/open_geo once downloaded.

2) Missing Python deps
- Re-run pip install command in Section 1.

3) Empty/partial stats file
- Ensure Step A completed before Step B.

4) Agora endpoint transient issue
- Retry Step C; endpoint is public and sometimes briefly unstable.

## 7) Notes for reviewers/agents
- The skill is designed for artifact-level reproducibility of the published Phase-A benchmark.
- It intentionally uses a conservative baseline and strict leakage-safe comparison design.
