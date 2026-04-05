# Leakage-Safe Cross-Cohort Alzheimer’s Blood Transcriptomic Prediction on Open Data: Consistent Permutation Nulls, AMP-AD Feature Ablations, and Sensitivity Analyses

**Pranjal**

## Abstract
Cross-cohort Alzheimer’s disease (AD) blood transcriptomic prediction is sensitive to cohort shift and can be misinterpreted without strict evaluation controls. We present an open reproducible study on GEO cohorts GSE63060 and GSE63061 with three design principles: leakage-safe target holdout evaluation, consistent permutation-null reporting, and explicit biological feature ablations using open AMP-AD Agora nominated targets. Primary arms are target_only, source_only, and pooled source+target raw training; a transductive ComBat arm is reported as sensitivity analysis only. Feature modes are variance, DE t-test, Agora-only, and DE-Agora intersection. Across settings, mean permutation-null AUROC remains near chance (0.4887-0.5132). In primary analyses, target_only exceeds permutation-null means in both directions with multiplicity-controlled significance in multiple settings. Additional local sensitivity checks show that conclusions are model-stable across logistic regression, linear SVM, and random forest for the DE-1000 setting, and that increasing null simulations to 1000 permutations preserves chance-centered null behavior. Conservative conclusion: robust target-domain signal is reproducible, while cross-cohort transfer gains are directional and should be interpreted with explicit batch-policy caveats.

## 1. Introduction
Public AD blood transcriptomic resources are valuable for reproducible benchmarking, but cross-cohort evaluation can overstate transferability when leakage controls, null calibration, and harmonization assumptions are not explicit. This work focuses on evaluation rigor and reproducibility, not on proposing a new classifier architecture.

We address three practical questions:
1) Does target-domain signal remain above a consistent permutation-null baseline under leakage-safe splits?
2) Are transfer effects (source_only or pooled source+target training) direction-stable across cohorts?
3) Do AMP-AD-informed feature restrictions materially change predictive behavior?

Recent work has reported increasingly complex machine learning approaches for AD blood transcriptomic prediction, including deep-learning/XAI feature-selection pipelines and digital-diagnosis signatures [11], [12]. Related multi-omics machine-learning frameworks also continue to improve apparent classification performance in AD-focused settings [13]. Against this landscape, our study is intentionally conservative: we prioritize leakage-safe split discipline, explicit null calibration, and transparent sensitivity labeling to reduce over-optimistic transport claims in small cross-cohort settings.

## 2. Data
### 2.1 GEO cohorts
- GSE63060 and GSE63061 from NCBI GEO [1], [2]
- AD vs CTL labels only; ambiguous status labels excluded
- Directional setup: GSE63060->GSE63061 and GSE63061->GSE63060

### 2.2 AMP-AD open biological context
We use open Agora nominated targets from AD Knowledge Portal [3], [4]. Agora signals are integrated directly into model feature-space ablations (Agora-only and DE-Agora intersection), rather than treated only as narrative context.

## 3. Methods
### 3.1 Primary leakage-safe protocol
For each direction, target cohort is split into target-train/target-test (stratified, random_state=42). Primary arms:
- target_only: train on target-train, evaluate on target-test
- source_only: train on source cohort, evaluate on target-test
- source_plus_target_raw: train on concatenated source + target-train, evaluate on target-test

No target-test labels are used in feature selection, model fitting, or null generation.

### 3.2 Sensitivity arm (not primary evidence)
- source_plus_target_combat_transductive: ComBat applied on stacked train+test features (no labels) prior to fitting [7].
This arm is retained to quantify harmonization sensitivity but is excluded from primary leakage-safe claims.

### 3.3 Feature modes
- var: top-N variance probes
- de_ttest: top-N probes by absolute AD-vs-CTL t-statistic on target-train only [8]
- agora_only: top-N probes mapped to Agora nominated symbols
- de_agora_intersection: top-N DE-ranked probes within the Agora-mapped subset

N in {200, 1000}.

### 3.4 Null definition and consistency policy
Primary null is the distribution of AUROC values from label-permuted target-train models (100 permutations per setting in main benchmark). The null reported in primary tables is the mean of this permutation-AUROC distribution. This avoids mixing incomparable null definitions across sections.

### 3.5 Statistical inference
- AUROC (primary), AUPRC, balanced accuracy, Brier
- Bootstrap CIs for AUROC [5]
- Paired bootstrap deltas for arm comparisons
- Benjamini-Hochberg multiplicity control [6]

### 3.6 Additional local sensitivity experiments
To address scope concerns, we add two local sensitivity checks:
1) Model-family sensitivity (DE-1000 target_only): logistic regression, linear SVM [9], random forest [10].
2) Null-stability sensitivity: 1000 permutations (DE-1000, both directions).

\newpage

## 4. Results
### 4.1 Primary DE setting (top 200/1000)

| Direction | Top genes | target_only | source_only | source+target raw | null (perm mean AUROC) |
|---|---:|---:|---:|---:|---:|
| GSE63060->GSE63061 | 200  | 0.7208 | 0.7565 | 0.8089 | 0.4903 |
| GSE63060->GSE63061 | 1000 | 0.6958 | 0.7488 | 0.8179 | 0.4986 |
| GSE63061->GSE63060 | 200  | 0.8453 | 0.8365 | 0.9003 | 0.4887 |
| GSE63061->GSE63060 | 1000 | 0.8908 | 0.8636 | 0.9120 | 0.4961 |

Target-domain signal remains clearly above chance-centered nulls. Transfer uplift exists but is directional.

### 4.2 AMP-AD feature ablations

| Direction | Feature mode | Top genes | target_only | null (perm mean AUROC) | Delta (target-null) |
|---|---|---:|---:|---:|---:|
| GSE63060->GSE63061 | agora_only | 200  | 0.7292 | 0.4994 | +0.2297 |
| GSE63060->GSE63061 | de_agora_intersection | 1000 | 0.6643 | 0.5069 | +0.1574 |
| GSE63061->GSE63060 | agora_only | 200  | 0.8732 | 0.4927 | +0.3804 |
| GSE63061->GSE63060 | de_agora_intersection | 200  | 0.8952 | 0.4971 | +0.3981 |

Agora-constrained feature spaces retain measurable signal, though performance varies by direction and feature policy.

### 4.3 Model-family sensitivity (local)
DE-1000 target_only AUROC:

| Direction | Logistic regression | Linear SVM | Random forest |
|---|---:|---:|---:|
| GSE63060->GSE63061 | 0.6958 | 0.6940 | 0.7277 |
| GSE63061->GSE63060 | 0.8908 | 0.8842 | 0.8761 |

The central conclusion (signal above null, directional transfer behavior) is stable across model families.

### 4.4 Null-stability sensitivity (1000 permutations, DE-1000)

| Direction | Null mean AUROC | Null SD | q05 | q95 |
|---|---:|---:|---:|---:|
| GSE63060->GSE63061 | 0.5006 | 0.0729 | 0.3821 | 0.6143 |
| GSE63061->GSE63060 | 0.4947 | 0.0888 | 0.3489 | 0.6452 |

Increasing null simulations to 1000 preserves chance-centered behavior and supports calibration robustness.

## 5. Discussion
Three findings are robust across this benchmark:
1) Leakage-safe target-domain signal is reproducibly above permutation-null baselines.
2) Cross-cohort transfer effects are directional and not universally positive.
3) Biological feature restrictions (Agora-only and DE-Agora intersection) remain informative but do not remove directional sensitivity.

The transductive ComBat arm can be informative for sensitivity analysis, but it is excluded from primary evidence because ComBat parameters are estimated on stacked train+test features, which violates a strict predictive boundary for target-holdout evaluation [7]. In other words, even without labels, test-distribution information enters the harmonization step and can inflate apparent transportability; we therefore treat this arm as diagnostic only.

In directional settings where source+target raw outperforms target_only, a practical explanation is that the gain from larger pooled training size can, in some cohort directions, outweigh uncorrected batch-noise penalties by better capturing shared disease-associated signal.

The relatively high target_only AUROC in GSE63061->GSE63060 (0.8908 for DE-1000) is interpreted cautiously as a property of this specific curated binary AD-vs-CTL setup and cohort composition, not as proof of universally easy blood-based AD classification.

A local predictive-harmonization audit was run to test train-only ComBat parameterization; the current pycombat implementation failed at train-fit/test-transform stage because test batches did not match fit-time category requirements. This supports retaining ComBat as a sensitivity-only arm until a strict train-only harmonization alternative is integrated.

## 6. Limitations
This study still uses two cohorts and moderate sample sizes. Main benchmark null estimation uses 100 permutations per setting; we therefore include explicit 1000-permutation stability checks for a representative DE-1000 setting. External prospective validation and additional harmonization strategies with strict train-only parameterization remain future work.

## 7. Conclusion
A reproducible leakage-safe evaluation pipeline on open AD blood transcriptomic cohorts shows stable target-domain signal above chance, with transfer gains that are conditional on direction and feature policy. Consistent null definitions and explicit sensitivity analyses improve interpretability and reduce claim inflation in cross-cohort settings.

## 8. Reproducibility
Code and artifacts: https://github.com/githubbermoon/bio-paper-track-open-phasea

Run sequence:
1) `python src/ingest/fetch_ampad_open_subset.py`
2) `python src/train/run_open_phaseA_benchmark.py`
3) `python src/eval/compute_open_phaseA_bootstrap.py`
4) `python src/eval/model_family_sensitivity.py`
5) `python src/eval/null_stability_check.py`

Core outputs:
- `outputs/metrics/open_phaseA_main_results.csv`
- `outputs/metrics/open_phaseA_predictions.csv`
- `outputs/stats/open_phaseA_null_distribution.csv`
- `outputs/stats/open_phaseA_auroc_ci.csv`
- `outputs/stats/open_phaseA_paired_tests.csv`
- `outputs/stats/open_phaseA_model_family_sensitivity.csv`
- `outputs/stats/open_phaseA_null_stability_de1000_perm1000.csv`
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

[9] C. Cortes and V. Vapnik, “Support-vector networks,” Machine Learning, 20:273-297, 1995. doi:10.1007/BF00994018.

[10] L. Breiman, “Random forests,” Machine Learning, 45:5-32, 2001. doi:10.1023/A:1010933404324.

[11] H. Lei et al., “Alzheimer's disease prediction using deep learning and XAI based interpretable feature selection from blood gene expression data,” Scientific Reports, 2026. PMID: 41667529.

[12] M. Altab et al., “A machine learning-enabled blood transcriptomic signature for digital diagnosis and subtyping of Alzheimer's disease,” npj Digital Medicine, 2026. PMID: 41491414.

[13] S. Kumar et al., “An integrative multiomics random forest framework for robust biomarker discovery,” GigaScience, 2026. PMID: 41363728.
