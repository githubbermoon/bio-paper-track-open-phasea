Open-Source Phase-A Dataset Choice (No waiting)
Date: 2026-04-04

Decision: Pivot to fully open / immediate-access data for a fast publishable computational-bioinformatics paper.

Best dataset stack for our objective:
1) AMP-AD / AD Knowledge Portal (open-access subsets first)
   Link: https://adknowledgeportal.synapse.org/
   Why: highest relevance to Alzheimer’s multi-omics; strongest scientific weight.

2) NCBI GEO Alzheimer’s transcriptomics cohorts (bulk + blood/brain expression)
   Link: https://www.ncbi.nlm.nih.gov/geo/
   Why: immediate public download, supports cross-cohort transfer and robustness tests.

3) OASIS-1 MRI aging/dementia cohort (open-access terms)
   Link: https://www.oasis-brains.org/
   Why: independent external validation modality (imaging) for transportability checks.

Scope change required:
- Remove India individual-level validation claim for Phase-A (cannot be done immediately with open data alone).
- Keep India as Phase-B extension once LASI-DAD access is approved.

Phase-A publishable claim:
- Leakage-safe cross-cohort transfer on open Alzheimer’s omics/imaging cohorts improves robustness/calibration vs no-transfer baselines, with strict temporal/group split auditing.

Immediate execution order:
A) Pull AMP-AD open subsets + 2 GEO cohorts
B) Harmonize common feature schema (demographics + diagnosis + omics modules)
C) Run no-transfer vs transfer benchmarks under strict split protocols
D) Null controls + ablations + calibration/failure audit
E) Write clawRxiv paper with verified numbers only
