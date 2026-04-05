from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


def _bh_adjust(pvals: list[float]) -> list[float]:
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = np.array(pvals)[order]
    adj = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = min(prev, ranked[i] * n / rank)
        prev = val
        adj[i] = val
    out = np.empty(n)
    out[order] = adj
    return out.tolist()


def _bootstrap_ci(y: np.ndarray, p: np.ndarray, n_boot: int = 1000, seed: int = 42) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    vals = []
    n = len(y)
    idx_all = np.arange(n)
    for _ in range(n_boot):
        idx = rng.choice(idx_all, size=n, replace=True)
        yb = y[idx]
        pb = p[idx]
        if len(np.unique(yb)) < 2:
            continue
        vals.append(roc_auc_score(yb, pb))
    if not vals:
        base = roc_auc_score(y, p)
        return base, base, base
    arr = np.array(vals)
    return float(np.mean(arr)), float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975))


def run() -> None:
    root = Path(__file__).resolve().parents[2]
    pred_path = root / "outputs/metrics/open_phaseA_predictions.csv"
    out_ci = root / "outputs/stats/open_phaseA_auroc_ci.csv"
    out_p = root / "outputs/stats/open_phaseA_paired_tests.csv"

    df = pd.read_csv(pred_path)

    ci_rows = []
    p_rows = []

    for (src, tgt, top_n), g in df.groupby(["source", "target", "top_n_genes"]):
        arms = {a: x.copy() for a, x in g.groupby("arm")}

        # CIs per arm
        for arm, ga in arms.items():
            y = ga["y_true"].to_numpy(dtype=int)
            p = ga["y_prob"].to_numpy(dtype=float)
            base = roc_auc_score(y, p)
            mean_b, lo, hi = _bootstrap_ci(y, p)
            ci_rows.append(
                {
                    "source": src,
                    "target": tgt,
                    "top_n_genes": int(top_n),
                    "arm": arm,
                    "auroc": float(base),
                    "auroc_boot_mean": mean_b,
                    "auroc_ci95_lo": lo,
                    "auroc_ci95_hi": hi,
                }
            )

        # paired comparisons on same sample_ids
        comp_defs = [
            ("source_plus_target", "target_only", "transfer_vs_target_only"),
            ("target_only", "null_label_permutation", "target_only_vs_null"),
        ]

        for a1, a2, cname in comp_defs:
            g1 = arms[a1].set_index("sample_id").sort_index()
            g2 = arms[a2].set_index("sample_id").sort_index()
            common = g1.index.intersection(g2.index)
            g1 = g1.loc[common]
            g2 = g2.loc[common]
            y = g1["y_true"].to_numpy(dtype=int)
            p1 = g1["y_prob"].to_numpy(dtype=float)
            p2 = g2["y_prob"].to_numpy(dtype=float)

            delta_obs = roc_auc_score(y, p1) - roc_auc_score(y, p2)
            rng = np.random.default_rng(42)
            deltas = []
            n = len(y)
            idx_all = np.arange(n)
            for _ in range(1000):
                idx = rng.choice(idx_all, size=n, replace=True)
                yb = y[idx]
                if len(np.unique(yb)) < 2:
                    continue
                deltas.append(roc_auc_score(yb, p1[idx]) - roc_auc_score(yb, p2[idx]))
            arr = np.array(deltas) if deltas else np.array([delta_obs])
            p_two = float(2 * min((arr <= 0).mean(), (arr >= 0).mean()))
            p_rows.append(
                {
                    "source": src,
                    "target": tgt,
                    "top_n_genes": int(top_n),
                    "comparison": cname,
                    "delta_auroc": float(delta_obs),
                    "ci95_lo": float(np.quantile(arr, 0.025)),
                    "ci95_hi": float(np.quantile(arr, 0.975)),
                    "p_value": min(1.0, p_two),
                }
            )

    p_df = pd.DataFrame(p_rows)
    p_df = p_df.sort_values(["source", "target", "top_n_genes", "comparison"]).reset_index(drop=True)
    p_df["bh_adjusted_p"] = _bh_adjust(p_df["p_value"].tolist())

    pd.DataFrame(ci_rows).to_csv(out_ci, index=False)
    p_df.to_csv(out_p, index=False)

    summary = {
        "ci_file": str(out_ci),
        "paired_tests_file": str(out_p),
        "n_ci_rows": int(len(ci_rows)),
        "n_tests": int(len(p_df)),
    }
    (root / "outputs/stats/open_phaseA_stats_manifest.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    run()
