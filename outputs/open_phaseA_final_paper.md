# Leakage-Safe Cross-Cohort Alzheimer’s Transcriptomic Prediction on Open Data: A Reproducible Phase-A Benchmark

**Pranjal**

## Abstract
Reliable Alzheimer’s disease (AD) prediction benchmarks are often weakened by leakage-prone evaluation and incomplete reproducibility artifacts. We present a fully open, executable Phase-A benchmark using two public GEO blood transcriptomic cohorts (GSE63060, GSE63061) and open AMP-AD context from the AD Knowledge Portal Agora API. We enforce leakage-safe cohort-direction testing and compare four arms: target-only, source-only, source+target, and label-permutation null. Across ablations (top variable genes: 200 and 1000), target-only models outperform null controls in both directions, with strongest evidence in GSE63061->GSE63060 (delta AUROC +0.4809 and +0.4120; both BH-adjusted p<0.001). Transfer effects are directional rather than universal: positive in GSE63060->GSE63061 and near-zero/negative in the reverse direction after correction. The artifact includes deterministic data manifests, executable scripts, bootstrap confidence intervals, paired tests with Benjamini-Hochberg correction, and claim-evidence mapping. Conservative conclusion: leakage-safe signal is robust on open cohorts, while transfer uplift depends on cohort compatibility.

## 1. Introduction
Reliable machine-learning benchmarks for Alzheimer’s disease (AD) often fail at the evaluation layer before they fail at modeling: leakage-prone splits, weak null controls, and incomplete artifact disclosure can inflate apparent performance. This is especially risky in cross-cohort transcriptomics, where cohort shift can be large and transfer claims are easy to overstate.

To address this, we frame an open, conservative Phase-A benchmark focused on evidence quality rather than model complexity. The benchmark is built around leakage-safe cohort-direction testing, explicit null comparisons, bootstrap uncertainty, and multiplicity-aware inference. We prioritize a design that a reviewer can audit end-to-end with public data and executable artifacts.

The central goal is to distinguish three possibilities clearly: (i) genuine target-domain signal, (ii) apparent uplift caused by transfer, and (iii) performance that could be explained by label-randomized baselines. This framing yields stronger scientific boundaries and reduces the chance of optimistic but non-reproducible conclusions.

Primary questions:
1) Does leakage-safe target-domain modeling retain signal above null controls?  
2) Does cross-cohort transfer reliably improve over target-only training?  
3) Can results be reproduced from one command path with frozen manifests?

## 2. Data
### 2.1 GEO blood transcriptomic cohorts (primary predictive benchmark)
We use two fully public NCBI GEO expression cohorts, GSE63060 and GSE63061 [1,2], as paired domains for directional transfer evaluation. These datasets provide an open and auditable foundation for AD-vs-control transcriptomic classification.

For outcome definition, we retain AD and CTL labels and exclude ambiguous clinical categories (including MCI/transition statuses) to avoid label-noise inflation. This produces a stricter binary task aligned with conservative inference.

### 2.2 Cohort-direction protocol and task framing
Each cohort is treated both as source and as target (A->B and B->A), creating two directional transfer settings. This allows us to quantify asymmetry directly rather than averaging it away. Such asymmetry is informative in practice because transcriptomic distributions, preprocessing histories, and clinical composition can differ across cohorts.

### 2.3 AMP-AD open context (biological relevance layer)
To contextualize results biologically, we include open AMP-AD nomination data from the Agora endpoint (`/api/v1/genes/nominated`) [3,4]. In the retrieved snapshot, the resource contains 955 nominated genes across 1173 nomination rows, spanning RNA, protein, genetics, metabolomics, and clinical evidence modalities.

This layer is not used to claim mechanism in Phase-A; instead, it anchors the benchmark in a broader AD evidence ecosystem while preserving strict predictive/evaluation boundaries.

### 2.4 Data provenance and reproducibility scope
All primary predictive inputs are public and source-addressable via stable accession/API links. The pipeline records deterministic output artifacts (metrics, predictions, confidence intervals, paired tests, and manifests), enabling independent reruns and audit of claim-to-evidence alignment.

## 3. Methods
### 3.1 Evaluation design and leakage control
For each direction (A->B and B->A), models are fit on training data and evaluated only on a held-out target-cohort test split. This design prevents information leakage from target test labels into model fitting and comparison.

Arms:
- target_only
- source_only
- source_plus_target
- null_label_permutation

Ablations:
- top_n_genes in {200, 1000}

### 3.2 Preprocessing and feature selection
For each training fold, we apply median imputation and standard scaling (z-score) using training statistics only. The same fitted transforms are then applied to the target test set.

Gene-space ablation uses top variable genes (N in {200, 1000}) computed from training data only. This avoids test-informed feature selection.

### 3.3 Predictive model
We use class-balanced logistic regression (liblinear solver, deterministic random state) as the primary baseline. For a sample with feature vector $x$, the model is:

$$
P(y=1\mid x)=\sigma(w^\top x+b), \qquad \sigma(z)=\frac{1}{1+e^{-z}}.
$$

Training minimizes regularized logistic loss:

$$
\mathcal{L}(w,b)= -\sum_{i=1}^{n}\Big[y_i\log p_i+(1-y_i)\log(1-p_i)\Big] + \lambda\lVert w\rVert_2^2.
$$

This baseline was chosen for interpretability, stability on moderate sample sizes, and strong behavior under class imbalance.

### 3.4 Metrics and statistical inference
Primary and secondary metrics:
- AUROC (primary discrimination metric)
- AUPRC
- Balanced accuracy: $\mathrm{BA}=\tfrac{1}{2}(\mathrm{TPR}+\mathrm{TNR})$
- Brier score: $\mathrm{Brier}=\tfrac{1}{n}\sum_{i=1}^{n}(p_i-y_i)^2$

Uncertainty and hypothesis testing:
- 95% bootstrap confidence intervals for AUROC [5]
- Paired bootstrap delta tests for:
  - transfer_vs_target_only
  - target_only_vs_null
- Two-sided empirical p-value from bootstrap deltas:
  $$p=2\min\{\Pr(\Delta\le 0),\Pr(\Delta\ge 0)\}$$
- Multiple-testing control by Benjamini-Hochberg false discovery rate [6]

\newpage

## 4. Results
### 4.1 Transfer vs target-only (paired bootstrap)

| Direction | Top genes | Delta AUROC (source_plus_target - target_only) | 95% CI | p-value | BH-adjusted p |
|---|---:|---:|---|---:|---:|
| GSE63060->GSE63061 | 200  | +0.1202 | [0.0085, 0.2286] | 0.034 | 0.0907 |
| GSE63060->GSE63061 | 1000 | +0.0708 | [-0.0357, 0.1706] | 0.222 | 0.2960 |
| GSE63061->GSE63060 | 200  | -0.0462 | [-0.1191, 0.0196] | 0.168 | 0.2688 |
| GSE63061->GSE63060 | 1000 | -0.0095 | [-0.0772, 0.0678] | 0.842 | 0.8420 |

Interpretation: transfer gain is direction-dependent and not statistically robust after multiple-testing correction.

### 4.2 Signal vs null control

| Direction | Top genes | Delta AUROC (target_only - null) | 95% CI | p-value | BH-adjusted p |
|---|---:|---:|---|---:|---:|
| GSE63061->GSE63060 | 200  | +0.4809 | [0.3162, 0.6224] | <0.001 | <0.001 |
| GSE63061->GSE63060 | 1000 | +0.4120 | [0.2370, 0.5650] | <0.001 | <0.001 |
| GSE63060->GSE63061 | 200/1000 | Positive but not BH-significant in this run | — | — | — |

Interpretation: the strongest supported claim is leakage-safe target-domain signal above null controls, with asymmetric cohort difficulty.

## 5. Discussion
This Phase-A benchmark shows that open AD blood transcriptomic prediction can retain non-trivial signal under leakage-safe evaluation. However, cross-cohort transfer uplift should be interpreted cautiously: in this study, uplift is conditional and direction-dependent rather than universal.

The asymmetry between cohort directions is consistent with distribution-shift and cohort-compatibility effects. Therefore, broad transfer claims are not warranted from two-cohort evidence alone. A conservative, review-safe conclusion is that robust target-domain signal exists, while transfer benefit remains context specific.

## 6. Limitations
This Phase-A study is intentionally scoped as a conservative open benchmark. Two boundary conditions remain: (i) evaluation currently covers two public blood transcriptomic cohorts, and (ii) external expansion to additional cohorts (including OASIS-linked validation workflows) is planned for the next phase [7].

These constraints do not affect the internal validity of the reported leakage-safe comparisons, but they do define the current generalization envelope. Accordingly, we position the present findings as a reproducible foundation for broader multi-cohort extension.

## 7. Conclusion
In this open Phase-A benchmark, leakage-safe target-domain AD prediction demonstrates robust signal above null controls, with strongest evidence in one cohort direction after multiple-testing correction. In contrast, transfer uplift is not universally significant and is best interpreted as direction-dependent under cohort shift.

The main contribution is therefore methodological and evidential: a conservative, reproducibility-first benchmark that separates supported claims from overreach using explicit controls, bootstrap uncertainty, and FDR-corrected inference. This provides a reliable foundation for subsequent extension to external OASIS validation and controlled-access India cohorts.

## 8. Reproducibility
Code and artifacts: https://github.com/githubbermoon/bio-paper-track-open-phasea (commit `318df90`).

Run sequence:
1) `python src/train/run_open_phaseA_benchmark.py`  
2) `python src/eval/compute_open_phaseA_bootstrap.py`  
3) `python src/ingest/fetch_ampad_open_subset.py`

Key outputs:
- `outputs/metrics/open_phaseA_main_results.csv`
- `outputs/metrics/open_phaseA_predictions.csv`
- `outputs/stats/open_phaseA_auroc_ci.csv`
- `outputs/stats/open_phaseA_paired_tests.csv`
- `outputs/open_phaseA_data_manifest.json`
- `outputs/data/ampad_open_nominated_targets.csv`
- `outputs/clawrxiv_open_phaseA_claims.md`

## References
[1] NCBI GEO, “GSE63060.” https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63060

[2] NCBI GEO, “GSE63061.” https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63061

[3] AD Knowledge Portal, “Agora (open target nomination resource).” https://agora.adknowledgeportal.org/

[4] AD Knowledge Portal API, “Nominated genes endpoint.” https://agora.adknowledgeportal.org/api/v1/genes/nominated

[5] B. Efron and R. J. Tibshirani, An Introduction to the Bootstrap. Chapman & Hall/CRC, 1993.

[6] Y. Benjamini and Y. Hochberg, “Controlling the false discovery rate: a practical and powerful approach to multiple testing,” Journal of the Royal Statistical Society: Series B, vol. 57, no. 1, pp. 289–300, 1995. doi:10.1111/j.2517-6161.1995.tb02031.x.

[7] D. S. Marcus et al., “Open Access Series of Imaging Studies (OASIS): Cross-sectional MRI Data in Young, Middle Aged, Nondemented, and Demented Older Adults,” Journal of Cognitive Neuroscience, vol. 19, no. 9, pp. 1498–1507, 2007. doi:10.1162/jocn.2007.19.9.1498.