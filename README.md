# bio-paper-track-open-phasea

Leakage-safe cross-cohort Alzheimer’s blood transcriptomic prediction benchmark on open GEO data (GSE63060, GSE63061), with AMP-AD Agora feature ablations and explicit null/sensitivity analyses.

## Exact regeneration sequence (fresh run)
Run from repository root:

```bash
python src/ingest/fetch_ampad_open_subset.py
python src/train/run_open_phaseA_benchmark.py
python src/eval/compute_open_phaseA_bootstrap.py
python src/eval/model_family_sensitivity.py
python src/eval/null_stability_check.py
```

## Primary generated outputs
- `outputs/metrics/open_phaseA_main_results.csv`
- `outputs/metrics/open_phaseA_predictions.csv`
- `outputs/stats/open_phaseA_null_distribution.csv`
- `outputs/stats/open_phaseA_auroc_ci.csv`
- `outputs/stats/open_phaseA_paired_tests.csv`
- `outputs/stats/open_phaseA_model_family_sensitivity.csv`
- `outputs/stats/open_phaseA_null_stability_de1000_perm1000.csv`
- `outputs/open_phaseA_data_manifest.json`

## Manuscript artifacts
- `outputs/open_phaseA_final_paper_v7.md`
- `outputs/open_phaseA_final_paper_v7.pdf`
- `outputs/open_phaseA_final_paper_v7.docx`

## Notes for automated evaluators
- Primary leakage-safe claims are from `target_only`, `source_only`, and `source_plus_target_raw` arms.
- `source_plus_target_combat_transductive` is sensitivity-only (not primary evidence).
- Null baseline in result tables is `null_label_permutation_mean_auroc`.
