name: bdpfs-adaptive-repro
description: Reproduce Adaptive BDP-FS (GMM-regularized) cross-cohort AD prediction with 1,000 permutations and biological target preservation validation.
allowed-tools: Bash(git *), Bash(cd *), Bash(python *), Bash(pip *), WebFetch

# Reproducibility Manifest

This manifest facilitates the exact reconstruction of the reported findings, including the predictive gains in the GSE63061$\to$GSE63060 direction and the biological target preservation metrics enabled by GMM-anchored regularization.

### ⏳ Timing & Resources
| Operation | Est. Time | Resource |
| :--- | :--- | :--- |
| Environment Setup | 1-2 min | Internet Access |
| Data Ingestion | 2-3 min | AD Knowledge Portal API |
| 1,000-Perm Benchmark | 5-8 min | CPU (Parallelized) |

### Step 1: Environment Baseline
```bash
git clone https://github.com/githubbermoon/bio-paper-track-open-phasea.git
cd bio-paper-track-open-phasea
git checkout v1.0.0-phaseA-v8
python -m pip install -r requirements.txt
python -c "import sklearn, scipy; print('ENV_OK')"
```
**Expected Output:** `ENV_OK`

### Step 2: Data & Logic Execution
```bash
python src/ingest/fetch_ampad_open_subset.py
python src/train/run_open_phaseA_benchmark.py
```
**Expected Output:** `BENCHMARK_COMPLETE: outputs/stats/open_phaseA_stats.json generated.`

### Step 3: Forensic Validation
Run the following check to verify the 1,000-permutation rigor and Agora Shield preservation:
```python
import json
from pathlib import Path

stats = json.loads(Path('outputs/stats/open_phaseA_stats.json').read_text())
n_perm = stats['de_ttest__GSE63060_to_GSE63061_top1000']['null_perm_n']
rescued = stats['de_batch_robust_v2__GSE63060_to_GSE63061_top200']['agora_genes_rescued_by_v2_shield']

print(f"STATISTICAL_RIGOR: {n_perm} permutations")
print(f"BIOLOGICAL_TARGET_SAFETY: {rescued} targets preserved")

assert n_perm == 1000
assert rescued == 164
print("VERIFICATION_SUCCESSFUL")
```

### ✅ Success Criteria
| Criterion | Metric | Threshold |
| :--- | :--- | :--- |
| Statistical Rigor | `null_perm_n` | == 1,000 |
| Biological Safety | `agora_rescued` | == 164 |
| Repository Sync | Git Tag | `v1.0.0-phaseA-v8` |

---
*Verified on main branch at tag v1.0.0-phaseA-v8.*
