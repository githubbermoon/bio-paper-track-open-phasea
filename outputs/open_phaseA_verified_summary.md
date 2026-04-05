# Open Phase-A Results Summary (Verified from outputs)

- Data: GEO GSE63060 + GSE63061 (AD vs CTL only; MCI/transition excluded).
- Arms: target_only, source_only, source_plus_target, null_label_permutation.

## Directional findings
- GSE63060_to_GSE63061_top200: ΔAUROC(source_plus_target-target_only)=0.1214; ΔAUROC(target_only-null)=0.1036
- GSE63061_to_GSE63060_top200: ΔAUROC(source_plus_target-target_only)=-0.0462; ΔAUROC(target_only-null)=0.4824
- GSE63060_to_GSE63061_top1000: ΔAUROC(source_plus_target-target_only)=0.0690; ΔAUROC(target_only-null)=0.1363
- GSE63061_to_GSE63060_top1000: ΔAUROC(source_plus_target-target_only)=-0.0066; ΔAUROC(target_only-null)=0.4179

Claim status:
- Leakage-safe signal above null: supported in both directions (positive target_only-null deltas).
- Transfer lift: mixed (one direction positive, one negative/near-zero), so claim as conditional/not universal.
