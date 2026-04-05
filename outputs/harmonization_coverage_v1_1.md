# Harmonization Coverage v1.1

Date: 2026-04-04
File: data/processed/harmonization_dictionary_v1_1.csv

Tier policy:
- Tier 1 (required for primary paper): demographics + time/site keys + dementia endpoint
- Tier 2 (strongly preferred): APOE4, education, imaging where available
- Tier 3 (optional multi-omics lift): PRS/RNA/protein modules

Cohort readiness snapshot (field-name placeholders pending data dictionary confirmation):
- LASI/LASI-DAD: strong for demographics/outcomes; omics likely sparse/partial
- ADNI: strong multimodal source (clinical + imaging + genomics)
- AMP-AD: strong omics source (RNA/protein)
- NIAGADS/ADSP: genomics-heavy supplemental source
- UKB: broad multimodal, but access/cost constraints

Non-negotiable leakage guards:
1) visit_date <= index_date for all features
2) no subject overlap across splits
3) site-aware grouped temporal splitting
4) scaler/imputer fit on train only

Blocking items before model runs:
- Replace placeholder field names with exact data-dictionary variables per cohort
- Confirm dementia endpoint definitions per cohort
- Confirm index date logic and incident vs prevalent label policy
