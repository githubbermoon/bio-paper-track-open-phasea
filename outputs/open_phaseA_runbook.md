# Open Phase-A Runbook

Workspace: /Users/pranjal/Projects/gitLocal/BioInf/bio_paper_track

## Command used
python src/train/run_open_phaseA_benchmark.py

## What it does
1) Downloads and caches GEO series_matrix files for GSE63060 and GSE63061
2) Parses AD/CTL labels from GEO status fields (excludes MCI/transition classes)
3) Runs leakage-safe cohort-direction benchmarks (target_only, source_only, source_plus_target, label-permutation null)
4) Runs ablation on top variable genes (top_n=200 and 1000)
5) Writes metrics + stats + manifest outputs

## Generated outputs
- outputs/metrics/open_phaseA_main_results.csv
- outputs/stats/open_phaseA_stats.json
- outputs/open_phaseA_data_manifest.json
- outputs/tables/open_phaseA_table_main.csv
- outputs/tables/open_phaseA_auroc_plot_data.csv
- outputs/open_phaseA_verified_summary.md
- outputs/clawrxiv_open_phaseA_claims.md
- outputs/data/ampad_open_nominated_targets.csv
- outputs/tables/ampad_open_subset_summary.csv
- outputs/ampad_open_subset_notes.md

## Current claim boundary
- Supported now: signal above null + leakage-safe protocol reproducibility
- Not yet universal: transfer gain (mixed by direction)
- India-specific individual-level claim: deferred to LASI-DAD phase
