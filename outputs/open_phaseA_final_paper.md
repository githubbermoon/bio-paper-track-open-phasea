# Leakage-Safe Cross-Cohort Alzheimer’s Blood Transcriptomic Prediction on Open Data: v5 with Consistent Permutation Nulls and AMP-AD Feature Ablations

**Pranjal**

## Abstract
Cross-cohort Alzheimer’s disease (AD) blood transcriptomic prediction is sensitive to cohort shift and vulnerable to over-interpretation when evaluation controls are inconsistent. We present a fully open v5 revision on GEO cohorts GSE63060 and GSE63061 that directly addresses reviewer blockers: (i) consistent null reporting based on permutation-AUROC distributions, (ii) explicit separation of leakage-safe primary analysis from transductive ComBat sensitivity analysis, and (iii) model-integrated AMP-AD feature ablations. We evaluate target_only, source_only, and pooled source+target raw training under leakage-safe target-holdout splits; ComBat on stacked train+test is retained as sensitivity only. Feature modes are variance, DE t-test, Agora-only, and DE-Agora intersection. Across all settings, mean permutation-null AUROC is near chance (0.4887-0.5132). In primary analyses, target_only outperforms permutation-null means in both directions and multiple feature settings after BH correction. Conservative conclusion: leakage-safe target-domain signal is reproducible, transfer gains are direction-dependent, and transductive harmonization outcomes should not be used as primary evidence.

## 1. Introduction
Public AD transcriptomic benchmarks can appear strong while still failing key validity checks: leakage-safe boundaries, null calibration consistency, and explicit handling of cross-study harmonization assumptions. This revision focuses on evaluation correctness and claim discipline rather than model novelty.

Primary goals:
1) verify leakage-safe target-domain signal against a consistent permutation null;
2) quantify directional transfer behavior under source_only and pooled source+target training;
3) test whether AMP-AD-informed feature restrictions change transfer patterns.

## 2. Data
### 2.1 Predictive cohorts
- GSE63060 (GPL6947) and GSE63061 (GPL10558) from NCBI GEO [1,2]
- AD vs CTL labels only; ambiguous statuses removed

### 2.2 AMP-AD open context and integration
We ingest open Agora nominated targets [3,4]. In v5, Agora is integrated into modeling through two explicit feature modes:
- agora_only: probe set restricted to probes mapping to Agora symbols
- de_agora_intersection: DE-ranked probes restricted to Agora-mapped probes

## 3. Methods
### 3.1 Leakage-safe primary arms
For each direction (A->B and B->A), using target holdout split (stratified, random_state=42):
- target_only: train on target-train, test on target-test
- source_only: train on source, test on target-test
- source_plus_target_raw: train on source + target-train, test on target-test

### 3.2 Transductive sensitivity arm (not primary)
- source_plus_target_combat_transductive: ComBat on stacked train+test features (no labels), then model fit/eval.
This arm is reported as sensitivity only and excluded from leakage-safe primary claims.

### 3.3 Feature modes
- var: top-N variance probes
- de_ttest: top-N absolute t-statistic probes (AD vs CTL on target-train only)
- agora_only: top-N probes from Agora-mapped subset
- de_agora_intersection: top-N DE probes within the overlap of DE-ranked and Agora-mapped probes

N in {200, 1000}.

### 3.4 Null policy (v5 consistency fix)
Primary null is the distribution of AUROC values from 100 label-permuted training runs per setting. Reported null in tables is the mean of that permutation-AUROC distribution (not AUROC of averaged probabilities).

### 3.5 Model and inference
Class-balanced logistic regression (liblinear) with median imputation and scaling. We report AUROC as primary metric; paired bootstrap deltas for arm-vs-arm comparisons; BH correction for multiplicity [5,6].

\newpage

## 4. Results
### 4.1 DE feature setting (top 200/1000)

| Direction | Top genes | target_only | source_only | source+target raw | null (perm mean AUROC) |
|---|---:|---:|---:|---:|---:|
| GSE63060->GSE63061 | 200  | 0.7208 | 0.7565 | 0.8089 | 0.4903 |
| GSE63060->GSE63061 | 1000 | 0.6958 | 0.7488 | 0.8179 | 0.4986 |
| GSE63061->GSE63060 | 200  | 0.8453 | 0.8365 | 0.9003 | 0.4887 |
| GSE63061->GSE63060 | 1000 | 0.8908 | 0.8636 | 0.9120 | 0.4961 |

### 4.2 AMP-AD integration settings

| Direction | Feature mode | Top genes | target_only | null (perm mean AUROC) | Delta (target-null) |
|---|---|---:|---:|---:|---:|
| GSE63060->GSE63061 | agora_only | 200  | 0.7292 | 0.4994 | +0.2297 |
| GSE63060->GSE63061 | de_agora_intersection | 1000 | 0.6643 | 0.5069 | +0.1574 |
| GSE63061->GSE63060 | agora_only | 200  | 0.8732 | 0.4927 | +0.3804 |
| GSE63061->GSE63060 | de_agora_intersection | 200  | 0.8952 | 0.4971 | +0.3981 |

### 4.3 Statistical highlights
- target_only_vs_null_perm_mean remains significant (BH<0.05) across multiple settings in both directions.
- In DE-1000, source_plus_target_raw_vs_target_only is positive and significant in GSE63060->GSE63061.
- Transductive ComBat sensitivity can improve some settings, but because it uses stacked train+test features, it is not used for leakage-safe primary claims.

### 4.4 Null calibration check
Across all analyzed settings, permutation-null AUROC means fall in 0.4887-0.5132, matching chance-level expectations and removing prior table/text inconsistency.

## 5. Discussion
v5 resolves two central validity issues from prior review cycles:
1) null reporting is now internally consistent (single primary null definition);
2) ComBat outcomes are clearly separated into transductive sensitivity analysis rather than primary leakage-safe evidence.

The core signal remains: target-domain AD prediction exceeds permutation null under strict split discipline. Transfer effects are directional and feature-mode dependent, so broad universal transfer claims remain unwarranted.

## 6. Limitations
This benchmark uses two cohorts and one baseline model family. Probe-to-symbol mapping for Agora integration depends on public platform annotations and may introduce mapping incompleteness. Non-transductive harmonization methods with guaranteed train-only parameterization should be evaluated in future versions.

## 7. Conclusion
In this open v5 benchmark, leakage-safe target-domain signal is reproducibly above a consistent permutation-null baseline. AMP-AD-integrated feature modes are now part of the predictive experiment space, and transfer gains remain context-specific rather than universal. Transductive ComBat results are reported as sensitivity only.

## 8. Reproducibility
Code and artifacts: https://github.com/githubbermoon/bio-paper-track-open-phasea

Run sequence:
1) `python src/train/run_open_phaseA_benchmark.py`
2) `python src/eval/compute_open_phaseA_bootstrap.py`
3) `python src/ingest/fetch_ampad_open_subset.py`

Core outputs:
- `outputs/metrics/open_phaseA_main_results.csv`
- `outputs/metrics/open_phaseA_predictions.csv`
- `outputs/stats/open_phaseA_null_distribution.csv`
- `outputs/stats/open_phaseA_auroc_ci.csv`
- `outputs/stats/open_phaseA_paired_tests.csv`
- `outputs/stats/open_phaseA_stats.json`
- `outputs/open_phaseA_data_manifest.json`

## References
[1] NCBI GEO, “GSE63060.” https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63060

[2] NCBI GEO, “GSE63061.” https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63061

[3] AD Knowledge Portal, “Agora.” https://agora.adknowledgeportal.org/

[4] AD Knowledge Portal API, “Nominated genes endpoint.” https://agora.adknowledgeportal.org/api/v1/genes/nominated

[5] B. Efron and R. J. Tibshirani, An Introduction to the Bootstrap. Chapman & Hall/CRC, 1993.

[6] Y. Benjamini and Y. Hochberg, “Controlling the false discovery rate: a practical and powerful approach to multiple testing,” JRSS-B, 57(1):289-300, 1995. doi:10.1111/j.2517-6161.1995.tb02031.x.

[7] W. E. Johnson, C. Li, and A. Rabinovic, “Adjusting batch effects in microarray expression data using empirical Bayes methods,” Biostatistics, 8(1):118-127, 2007. doi:10.1093/biostatistics/kxj037.

[8] G. K. Smyth, “Linear models and empirical bayes methods for assessing differential expression in microarray experiments,” Stat Appl Genet Mol Biol, 3:Article3, 2004. doi:10.2202/1544-6115.1027.
