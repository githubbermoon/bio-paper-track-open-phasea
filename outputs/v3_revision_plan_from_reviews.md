# V3 Revision Plan (Derived from v1 + v2 AI Reviews)

## Objective
Convert current manuscript from "rigorous but low-novelty/weakly positioned" to "technically sound, clearly scoped, and publication-credible" by directly addressing all repeated reviewer weaknesses.

## Consolidated reviewer weaknesses
1. Novelty framed as weak/trivial; contribution appears to be standard ML hygiene.
2. Non-standard wording/jargon (e.g., "Phase-A benchmark", "claim-evidence mapping") hurts field alignment.
3. Batch-effect concern: source+target merged arm interpreted as technically unsound without explicit harmonization strategy.
4. Results incompleteness in narrative: source_only arm underreported in main text/tables.
5. AMP-AD/Agora integration seen as superficial (descriptive but not model-influencing).
6. Literature context and external comparison are insufficient.

## Revision strategy (priority order)

### P0 — Must-fix before resubmission
A. Reframe contribution in standard language
- Remove "Phase-A" from title/body claims.
- Position paper as "leakage-safe cross-cohort transfer stress-test" and "reproducible baseline protocol".

B. Resolve batch-effect criticism explicitly
- State that direct cohort pooling is exploratory only.
- Downgrade source_plus_target from inferential claim status.
- Keep inferential focus on target_only vs null and source_only transfer diagnostics.
- Add explicit "no ComBat/harmonization in this version" boundary and future extension text.

C. Complete results reporting
- Add source_only results table in main Results section.
- Keep paired inferential table for transfer_vs_target_only and target_only_vs_null.

### P1 — High-value upgrades in this revision cycle
D. Strengthen related work / positioning
- Add concise related-work section: AD blood transcriptomics, cross-cohort transfer pitfalls, batch normalization standards.
- Add a comparison table structure (what this paper does vs common prior practices), even if direct metric comparability is limited.

E. Integrate AMP-AD with testable analysis
- Add an explicit post-hoc analysis block (current version: evidence-layer interpretation only).
- Define a concrete vNext integration path into feature strategy (e.g., constrained feature sets / enrichment tests).

### P2 — Next experimental extension (not blocking immediate revision)
F. Add harmonized transfer arm
- Implement ComBat/ComBat-like harmonization and evaluate alongside current arms.
- Report whether transfer effects remain after harmonization.

## Concrete manuscript edits to perform now
1. Title and abstract wording: remove "Phase-A" branding and overclaiming.
2. Methods: add dedicated subsection "Batch effects and cohort pooling policy".
3. Results: add "source_only" performance table in main body.
4. Discussion: convert from novelty claim to reliability/diagnostic contribution claim.
5. Limitations: keep concise, field-standard, non-defensive wording.
6. Reproducibility: keep repo + pinned commit + executable skill_md.

## Success criteria for v3
- No reviewer can claim missing source_only evidence.
- No reviewer can interpret source_plus_target as an unqualified causal/transfer claim.
- Core claim remains statistically supported: leakage-safe target_only > null in strongest direction.
- Manuscript language is field-standard and jargon-minimized.
