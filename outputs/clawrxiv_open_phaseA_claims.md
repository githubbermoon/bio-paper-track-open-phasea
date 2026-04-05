# clawRxiv Open Phase-A Claim Table (current)

Data scope: public GEO cohorts (GSE63060, GSE63061)
Evidence files:
- outputs/metrics/open_phaseA_main_results.csv
- outputs/stats/open_phaseA_stats.json
- outputs/tables/open_phaseA_table_main.csv
- outputs/open_phaseA_verified_summary.md

| Claim ID | Claim | Status | Evidence |
|---|---|---|---|
| OP1 | Leakage-safe target-only models outperform label-permutation null controls | verified | open_phaseA_stats.json (delta target_only-null > 0 in both directions) |
| OP2 | Source+target transfer always outperforms target-only | rejected (not universal) | mixed sign deltas in open_phaseA_stats.json |
| OP3 | Cross-cohort source-only models retain above-chance discriminative signal | partial | open_phaseA_main_results.csv (direction dependent) |
| OP4 | Open-data pipeline is reproducible from cached GEO series_matrix inputs | verified | open_phaseA_data_manifest.json + script path |

Submission-safe wording now:
- Emphasize robust signal-above-null and strict split hygiene.
- Present transfer gain as cohort-direction dependent, not guaranteed.
