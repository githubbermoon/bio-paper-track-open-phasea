# Leakage-Safe Cross-Cohort Alzheimer’s Blood Transcriptomic Prediction: A Reproducible v4 Transfer Stress-Test with DE Features and ComBat

**Pranjal**

## Abstract
Cross-cohort Alzheimer’s disease (AD) transcriptomic prediction is highly sensitive to cohort shift, feature construction, and leakage-prone evaluation. We present a fully open, reproducible v4 stress-test on GEO blood cohorts GSE63060 and GSE63061 with strict target-holdout evaluation, directional transfer (A->B and B->A), and explicit null controls. Relative to earlier versions, v4 adds (i) differential-expression-based feature selection (DE t-test on target-train only), (ii) an explicit ComBat harmonization arm for pooled source+target training, and (iii) a 100x permutation-averaged null baseline to stabilize chance-level calibration. Across settings, null AUROC distributions center near 0.5 (mean range 0.4887-0.4986), addressing reviewer concerns about null behavior. In the DE-feature setting, target_only remains significantly above null in both directions (e.g., delta AUROC +0.4406 and +0.4765 in GSE63061->GSE63060; BH-adjusted p<0.001), while transfer effects remain direction-dependent. ComBat-pooled transfer shows a significant uplift in one direction (GSE63060->GSE63061, DE-1000: delta +0.1060, BH=0.016) but not uniformly across all settings. Conservative conclusion: robust leakage-safe target-domain signal is reproducible; transfer gains are conditional and should be reported with harmonization and strict uncertainty controls.

## 1. Introduction
Public AD blood transcriptomic cohorts are useful for reproducible benchmarking but vulnerable to overclaiming when transfer is evaluated without strict controls. In particular, three recurring issues affect interpretation: test leakage, underpowered null modeling, and unmodeled cross-study batch effects.

This work is positioned as an evaluation-and-reproducibility contribution rather than a new classifier architecture. We provide a transparent stress-test protocol that reports target-only, source-only, pooled transfer (raw and ComBat-harmonized), and permutation null arms under the same split.

Primary questions:
1) Does target-domain signal remain above chance under leakage-safe evaluation?
2) Does cross-cohort transfer improve target performance, and is that improvement direction-stable?
3) How do DE features and ComBat harmonization change transfer conclusions versus variance-only baselines?

## 2. Data
### 2.1 GEO predictive cohorts
- GSE63060 and GSE63061 (NCBI GEO series_matrix) [1,2]
- Binary labels retained: AD vs CTL
- Ambiguous statuses excluded (e.g., MCI/transition/non-definitive labels)

### 2.2 AMP-AD open context
We ingest Agora nominated targets [3,4] to retain AD biological context and reproducibility linkage. In the open snapshot used here, the resource includes 955 unique nominated genes across 1173 nomination rows, with strong RNA/Protein/Genetics representation.

### 2.3 Cohort-direction setup
Each cohort serves as source and target (GSE63060->GSE63061 and GSE63061->GSE63060). Target data are split into train/test (stratified, random_state=42); all model comparison is on target test only.

## 3. Methods
### 3.1 Leakage-safe evaluation design
For each direction and feature setting:
- target_only: train on target-train, evaluate on target-test
- source_only: train on source, evaluate on target-test
- source_plus_target_raw: train on concatenated source + target-train (no harmonization)
- source_plus_target_combat: ComBat-harmonized pooled training/evaluation features
- null_label_permutation_avg100: average prediction from 100 label-permuted target-only models

No target-test labels are used in fitting, feature ranking, or null permutation generation.

### 3.2 Feature selection
Feature modes:
- var: top-N variance genes
- de_ttest: top-N absolute t-statistic genes from AD vs CTL on target-train only (a lightweight DE proxy aligned with limma-style differential-expression practice [7,9])

N in {200, 1000}. Gene selection is always performed after split, using only target-train labels.

### 3.3 ComBat harmonization
For pooled transfer analysis, we apply ComBat [8] using study-of-origin batch labels (source vs target) on stacked expression matrices (feature-only adjustment, no label input to ComBat). We report raw pooled and ComBat-pooled outcomes separately.

### 3.4 Predictive model
Class-balanced logistic regression baseline (liblinear):

$$
P(y=1\mid x)=\sigma(w^\top x+b), \quad \sigma(z)=\frac{1}{1+e^{-z}}.
$$

$$
\mathcal{L}(w,b)= -\sum_{i=1}^{n}\left[y_i\log p_i+(1-y_i)\log(1-p_i)\right]+\lambda\|w\|_2^2.
$$

### 3.5 Inference
- AUROC (primary), AUPRC, balanced accuracy, Brier
- Bootstrap CI for AUROC [5]
- Paired bootstrap deltas for arm comparisons
- Benjamini-Hochberg correction [6]

\newpage

## 4. Results
### 4.1 Primary arm AUROC (DE features)

| Direction | Top genes | target_only | source_only | source+target (ComBat) | null (perm avg100) |
|---|---:|---:|---:|---:|---:|
| GSE63060->GSE63061 | 200  | 0.7208 | 0.7565 | 0.7149 | 0.3804 |
| GSE63060->GSE63061 | 1000 | 0.6958 | 0.7488 | 0.8018 | 0.4774 |
| GSE63061->GSE63060 | 200  | 0.8453 | 0.8365 | 0.8812 | 0.4047 |
| GSE63061->GSE63060 | 1000 | 0.8908 | 0.8636 | 0.9076 | 0.4142 |

### 4.2 Paired bootstrap tests (DE features)

| Direction | Top genes | Comparison | Delta AUROC | 95% CI | BH-adjusted p |
|---|---:|---|---:|---|---:|
| GSE63060->GSE63061 | 200  | target_only vs null_avg100 | +0.3405 | [0.1908, 0.4986] | <0.001 |
| GSE63060->GSE63061 | 1000 | source+target_combat vs target_only | +0.1060 | [0.0413, 0.1834] | 0.016 |
| GSE63061->GSE63060 | 200  | target_only vs null_avg100 | +0.4406 | [0.2547, 0.6132] | <0.001 |
| GSE63061->GSE63060 | 1000 | target_only vs null_avg100 | +0.4765 | [0.2957, 0.6332] | <0.001 |

Interpretation: target-only signal is robust above null in both directions. ComBat-pooled transfer improves one DE setting (A->B, 1000 genes), but transfer uplift is not universally significant.

### 4.3 Null baseline calibration check
Across all DE settings, 100x-permutation null AUROC means are 0.4887-0.4986 (q05/q95 ranges span around 0.5), consistent with expected chance-centered behavior.

### 4.4 AMP-AD context summary
Open Agora ingestion confirms non-trivial AD evidence density (955 genes; 1173 nominations; RNA/Protein/Genetics dominant modalities), supporting disease-context relevance for future mechanism-aware feature constraints.

## 5. Discussion
The v4 update directly addresses prior review criticisms by adding DE-based selection, explicit ComBat harmonization, and stabilized null estimation. This changes the manuscript from a generic transfer benchmark into a stricter diagnostic protocol for identifying where transfer claims hold and where they collapse.

Two stable findings remain:
1) leakage-safe target-domain signal is reproducibly above null;
2) transfer benefit is direction- and configuration-dependent.

Thus, cross-cohort uplift should be reported conditionally, with harmonization strategy and uncertainty bounds made explicit.

## 6. Limitations
This benchmark still uses two cohorts and a single baseline model family. ComBat is applied in a transductive feature-only harmonization regime; fully prospective external validation across additional cohorts remains future work. AMP-AD is used here as contextual evidence, not yet as a strict causal or mechanistic model component.

## 7. Conclusion
In this reproducible v4 open-data benchmark, AD blood transcriptomic prediction shows robust leakage-safe target-domain signal above permutation null baselines. DE feature selection plus ComBat harmonization improves specific transfer settings, but transfer gains remain conditional rather than universal. The practical contribution is a transparent, executable protocol that sharpens claim boundaries for cross-cohort AD modeling.

## 8. Reproducibility
Code and artifacts: https://github.com/githubbermoon/bio-paper-track-open-phasea

Run sequence:
1) `python src/train/run_open_phaseA_benchmark.py`
2) `python src/eval/compute_open_phaseA_bootstrap.py`
3) `python src/ingest/fetch_ampad_open_subset.py`

Core outputs:
- `outputs/metrics/open_phaseA_main_results.csv`
- `outputs/metrics/open_phaseA_predictions.csv`
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

[7] G. K. Smyth, “Linear models and empirical bayes methods for assessing differential expression in microarray experiments,” Stat Appl Genet Mol Biol, 3:Article3, 2004. doi:10.2202/1544-6115.1027.

[8] W. E. Johnson, C. Li, and A. Rabinovic, “Adjusting batch effects in microarray expression data using empirical Bayes methods,” Biostatistics, 8(1):118-127, 2007. doi:10.1093/biostatistics/kxj037.

[9] M. E. Ritchie et al., “limma powers differential expression analyses for RNA-sequencing and microarray studies,” Nucleic Acids Res, 43(7):e47, 2015. doi:10.1093/nar/gkv007.
