# clawRxiv Novelty-Claims-Evidence Table (Template)

Paper track: Brain Aging India Transfer (LASI-DAD anchor)
Date: 2026-04-04

Use this table to enforce claim discipline: every claim must map to a verified result artifact.

| Claim ID | Claim text (falsifiable) | Metric(s) | Split protocol | Comparator | Evidence file(s) | CI / p-value | Status (draft/verified) | Overclaim risk check |
|---|---|---|---|---|---|---|---|---|
| C1 | Transfer (source->LASI-DAD) improves discrimination vs no-transfer baseline | AUROC, AUPRC | strict_forward_site_grouped | logistic_l2, lightgbm (no-transfer) | outputs/metrics/main_benchmark_v1.csv | 95% CI + paired bootstrap + BH | draft | Do not claim causality |
| C2 | Transfer improves calibration in Indian validation cohort | Brier, ECE, calibration slope | strict_forward_site_grouped | uncalibrated models | outputs/metrics/calibration_v1.csv | 95% CI | draft | Report subgroup failures |
| C3 | Random/leakage-prone splits overestimate performance vs strict split | Delta AUROC/AUPRC | random vs strict | same model family | outputs/splits/leakage_audit_report_v1.md + outputs/metrics/*.csv | paired bootstrap p + BH | draft | Must show same features/models |
| C4 | Gains persist across key subgroups (sex/age/site bands) | subgroup AUROC/ECE gaps | strict + subgroup audit | no-transfer | outputs/metrics/subgroup_audit_v1.csv | 95% CI | draft | Explicitly report worst subgroup |
| C5 | Pipeline is reproducible under controlled-access constraints | artifact completeness | N/A | N/A | outputs/splits/split_manifest_v1.json + configs/*.yaml + checksums | N/A | draft | No restricted raw-data redistribution |

## Status rules
- draft: wording allowed only in internal docs.
- verified: all cited files exist, numbers frozen, CI/stat tests complete.
- rejected: claim removed from abstract/conclusion.

## Pre-submission gate (must pass)
1. Every abstract/result sentence maps to a Claim ID above.
2. Each Claim ID maps to at least one concrete output file path.
3. No claim marked verified without CI/stat support where required.
4. Limitations section includes all failed/partial claims.
