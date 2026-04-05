---
name: open-phasea-ad-benchmark-repro
description: Reproduce the final v8 BDP-FS GMM-Soft cross-cohort benchmark with 1,000 permutations.
allowed-tools: Bash(python *), Bash(pip *), WebFetch
---

# Reproduction (v8 GMM-Soft Package)

## 1) Prepare Environment
```bash
python -m pip install -r requirements.txt
```

## 2) Execution (One-Command Test)
```bash
python src/ingest/fetch_ampad_open_subset.py
python src/train/run_open_phaseA_benchmark.py
```

## 3) Expected Artifacts
- outputs/metrics/open_phaseA_main_results.csv
- outputs/stats/open_phaseA_stats.json
- outputs/open_phaseA_data_manifest.json

## 4) Validation Proof
Run this after execution to confirm the high-rigor (1,000x) policy:
```python
import json
from pathlib import Path

# Confirm 1000 permutations
stats = json.loads(Path('outputs/stats/open_phaseA_stats.json').read_text())
for k, v in stats.items():
    assert v['null_perm_n'] == 1000, f"LOW RIGOR: {k}"

# Confirm Agora Rescue
agora_rescued = stats['de_batch_robust_v2__GSE63060_to_GSE63061_top200']['agora_genes_rescued_by_v2_shield']
assert agora_rescued == 164, f"LOGIC FAIL: {agora_rescued}"

print("VALIDATION_OK: Submission-Master Rigor Confirmed.")
```

## 5) Version Tagging
Tagged as `v1.0.0-phaseA-v8` on the `main` branch.
