# Open Phase-A Statistical Summary (CI + paired bootstrap)

Evidence files:
- outputs/stats/open_phaseA_auroc_ci.csv
- outputs/stats/open_phaseA_paired_tests.csv

## Key results (AUROC deltas)
- GSE63060->GSE63061 (top200): transfer_vs_target_only delta = +0.1202, CI [0.0085, 0.2286], p=0.034, BH=0.0907
- GSE63060->GSE63061 (top1000): transfer_vs_target_only delta = +0.0708, CI [-0.0357, 0.1706], p=0.222, BH=0.2960
- GSE63061->GSE63060 (top200): transfer_vs_target_only delta = -0.0462, CI [-0.1191, 0.0196], p=0.168, BH=0.2688
- GSE63061->GSE63060 (top1000): transfer_vs_target_only delta = -0.0095, CI [-0.0772, 0.0678], p=0.842, BH=0.8420

## Null-control checks
- GSE63061->GSE63060 top200: target_only_vs_null delta = +0.4809, CI [0.3162, 0.6224], p<0.001, BH<0.001
- GSE63061->GSE63060 top1000: target_only_vs_null delta = +0.4120, CI [0.2370, 0.5650], p<0.001, BH<0.001
- GSE63060->GSE63061 null comparisons are positive but not BH-significant in this run.

## Submission-safe interpretation
- Strongest robust claim: leakage-safe predictive signal above null is supported (especially in GSE63061->GSE63060 direction).
- Transfer uplift is directional/conditional and not universally significant after multiple-testing correction.
