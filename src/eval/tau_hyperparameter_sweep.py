"""
Tau Hyperparameter Sweep for BDP-FS v1 & v2 — Large-Cohort Bidirectional Evaluation
GSE63060 <-> GSE63061 (100+ samples each)

v1: Hard threshold, equal weight
v2: GMM adaptive threshold, α-weighted variance penalty, Agora Shield, soft exponential weighting
"""
import sys
from pathlib import Path

root = Path("/Users/pranjal/Projects/gitLocal/bioInf/bio_paper_track")
sys.path.append(str(root))

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.mixture import GaussianMixture
from neuroCombat import neuroCombat

from src.train.run_open_phaseA_benchmark import (
    load_geo_series,
    load_platform_probe_to_symbol,
    _common_genes,
    _combat_trainfit_apply_test,
    _feature_scores,
    _select_top_genes,
    _fit_predict,
    _metrics,
    SEED,
)


def run_tau_sweep():
    print("=" * 60)
    print("BDP-FS v1 vs v2 Sweep — GSE63060 <-> GSE63061 Bidirectional")
    print("=" * 60)

    raw_dir = root / "data/raw/open_geo"
    cache_dir = root / "data/raw/open_geo_platform"

    # Load Agora targets
    agora_df = pd.read_csv(root / "outputs/data/ampad_open_nominated_targets.csv")
    agora_symbols = set(agora_df["hgnc_symbol"].dropna().astype(str).str.upper())

    # Load and map to gene symbols
    data = {}
    probe_to_symbol = {}
    for gse in ["GSE63060", "GSE63061"]:
        X, y, gpl = load_geo_series(gse, raw_dir)
        p2s = load_platform_probe_to_symbol(gpl, cache_dir)
        probe_to_symbol.update(p2s)
        cols = [p2s.get(c, None) for c in X.columns]
        X.columns = cols
        X = X.loc[:, X.columns.notna()]
        X = X.T.groupby(level=0).mean().T
        data[gse] = (X, y)
        print(f"  {gse}: {X.shape[0]} samples, {X.shape[1]} genes")

    directions = [
        ("GSE63060", "GSE63061"),
        ("GSE63061", "GSE63060"),
    ]

    all_results = []

    for source_name, target_name in directions:
        print(f"\n{'='*60}")
        print(f"Direction: {source_name} -> {target_name}")
        print(f"{'='*60}")

        Xs, ys = data[source_name]
        Xt, yt = data[target_name]
        common = _common_genes(Xs, Xt)
        Xs_c = Xs[common]
        Xt_c = Xt[common]

        Xtr_t, Xte_t, ytr, yte = train_test_split(
            Xt_c, yt, test_size=0.3, random_state=SEED, stratify=yt
        )
        print(f"  Target test: {len(Xte_t)} samples")

        X_pool = pd.concat([Xs_c, Xtr_t], axis=0)
        y_pool = pd.concat([ys, ytr], axis=0)

        # neuroCombat
        batch_labels = pd.Series(
            [source_name] * len(Xs_c) + [target_name] * len(Xtr_t), index=X_pool.index
        )
        covars = pd.DataFrame({"batch": batch_labels.astype(str).values}, index=X_pool.index)
        combat_fit = neuroCombat(dat=X_pool.T, covars=covars, batch_col="batch")

        gamma_star = combat_fit["estimates"]["gamma.star"]
        delta_star = combat_fit["estimates"]["delta.star"]
        gamma_abs = np.mean(np.abs(gamma_star), axis=0)
        delta_abs = np.mean(np.abs(1.0 - delta_star), axis=0)
        g_s = pd.Series(gamma_abs, index=X_pool.columns)
        d_s = pd.Series(delta_abs, index=X_pool.columns)
        z_g = (g_s - g_s.mean()) / g_s.std()
        z_d = (d_s - d_s.mean()) / d_s.std()

        # DE scores
        _, de_pool = _feature_scores(X_pool, y_pool, common)
        agora_probe_set = {g for g in common if g.upper() in agora_symbols}

        # ── v1 baseline ──
        distortion_v1 = (z_g.abs() + z_d.abs()).fillna(0)
        genes_bl = _select_top_genes("de_ttest", 1000, common, de_pool, de_pool, agora_probe_set)
        Xtr_bl = pd.concat([Xs_c[genes_bl], Xtr_t[genes_bl]], axis=0)
        ytr_bl = pd.concat([ys, ytr], axis=0)
        batch_bl = pd.Series([source_name]*len(Xs_c) + [target_name]*len(Xtr_t), index=Xtr_bl.index)
        batch_te = pd.Series([target_name]*len(Xte_t), index=Xte_t[genes_bl].index)
        X_train_cb, X_test_cb = _combat_trainfit_apply_test(Xtr_bl, batch_bl, Xte_t[genes_bl], batch_te)
        p_bl = _fit_predict(X_train_cb, ytr_bl, X_test_cb)
        m_bl = _metrics(yte.values, p_bl)
        print(f"\n  [BASELINE de_ttest] AUROC: {m_bl['auroc']:.4f}")
        all_results.append({"direction": f"{source_name}->{target_name}", "mode": "de_ttest_baseline", "tau": "N/A", "gene_count": len(genes_bl), **m_bl})

        # ── v2: Compute adaptive GMM anchor (Masterpiece) ──
        ALPHA = 0.2
        distortion_v2 = (ALPHA * z_g.abs() + (1 - ALPHA) * z_d.abs()).fillna(0)
        agora_mask = distortion_v2.index.isin(agora_probe_set)
        distortion_v2_adj = distortion_v2.copy()
        distortion_v2_adj[agora_mask] = distortion_v2_adj[agora_mask] * 0.5

        # GMM adaptive anchor
        v2_vals = distortion_v2_adj.values.reshape(-1, 1)
        gmm = GaussianMixture(n_components=2, random_state=SEED).fit(v2_vals)
        means = gmm.means_.flatten()
        stds = np.sqrt(gmm.covariances_.flatten())
        idx_native = np.argmin(means)
        tau_0 = float(means[idx_native] + 1.645 * stds[idx_native])
        print(f"  [v2 Masterpiece] Anchor tau_0: {tau_0:.4f} (95th-pct of Native component)")

        # v2 Masterpiece Result (Fixed self-calibrated tau_0, alpha=1.0)
        DECAY_ALPHA = 1.0
        weight_v2 = np.exp(-DECAY_ALPHA * np.maximum(0, distortion_v2_adj - tau_0))
        v2_weighted_de = de_pool * weight_v2.reindex(de_pool.index, fill_value=1.0)
        
        genes_v2 = _select_top_genes("de_batch_robust_v2", 1000, common, de_pool, de_pool, agora_probe_set, None, v2_weighted_de)
        if len(genes_v2) >= 20:
            Xtr_v2 = pd.concat([Xs_c[genes_v2], Xtr_t[genes_v2]], axis=0)
            ytr_v2 = pd.concat([ys, ytr], axis=0)
            bt_v2 = pd.Series([source_name]*len(Xs_c)+[target_name]*len(Xtr_t), index=Xtr_v2.index)
            bte_v2 = pd.Series([target_name]*len(Xte_t), index=Xte_t[genes_v2].index)
            Xtcb_v2, Xttcb_v2 = _combat_trainfit_apply_test(Xtr_v2, bt_v2, Xte_t[genes_v2], bte_v2)
            p_v2 = _fit_predict(Xtcb_v2, ytr_v2, Xttcb_v2)
            m_v2 = _metrics(yte.values, p_v2)
            delta_v2 = m_v2["auroc"] - m_bl["auroc"]
            print(f"  [v2 Masterpiece] AUROC={m_v2['auroc']:.4f}  delta={delta_v2:+.4f}  (tau_0={tau_0:.4f})")
            all_results.append({"direction": f"{source_name}->{target_name}", "mode": "bdp_fs_v2_masterpiece", "tau": "Adaptive", "gene_count": len(genes_v2), **m_v2})

        # Sweep v1 (Legacy comparison)
        tau_percentiles = [1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.60, 0.50]
        for tau_pct in tau_percentiles:
            thresh_v1 = distortion_v1.quantile(tau_pct)
            robust_v1 = distortion_v1[distortion_v1 <= thresh_v1].index
            genes_v1 = _select_top_genes("de_batch_robust", 1000, common, de_pool, de_pool, agora_probe_set, robust_v1)
            if len(genes_v1) >= 20:
                Xtr_c = pd.concat([Xs_c[genes_v1], Xtr_t[genes_v1]], axis=0)
                ytr_c = pd.concat([ys, ytr], axis=0)
                bt = pd.Series([source_name]*len(Xs_c)+[target_name]*len(Xtr_t), index=Xtr_c.index)
                bte = pd.Series([target_name]*len(Xte_t), index=Xte_t[genes_v1].index)
                Xtcb, Xttcb = _combat_trainfit_apply_test(Xtr_c, bt, Xte_t[genes_v1], bte)
                p_v1 = _fit_predict(Xtcb, ytr_c, Xttcb)
                m_v1 = _metrics(yte.values, p_v1)
                delta_v1 = m_v1["auroc"] - m_bl["auroc"]
                all_results.append({"direction": f"{source_name}->{target_name}", "mode": "bdp_fs_v1_legacy", "tau": tau_pct, "gene_count": len(genes_v1), **m_v1})

    df = pd.DataFrame(all_results)
    out = root / "outputs/stats/tau_sweep_metrics.csv"
    df.to_csv(out, index=False)
    print(f"\n{'='*60}")
    print(f"Results saved to {out}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    run_tau_sweep()
