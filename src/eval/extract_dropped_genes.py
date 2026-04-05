import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from scipy.stats import zscore, ttest_ind
from sklearn.model_selection import train_test_split

# Add root to sys.path
root = Path("/Users/pranjal/Projects/gitLocal/bioInf/bio_paper_track")
sys.path.append(str(root))

from src.train.run_open_phaseA_benchmark import (
    load_geo_series,
    load_platform_probe_to_symbol,
    _common_genes,
    SEED,
)
from neuroCombat import neuroCombat

def main():
    print("--- Identifying 'Dropped High-DE Genes' in GSE63060 -> GSE63061 ---")
    
    raw_dir = root / "data/raw/open_geo"
    cache_dir = root / "data/raw/open_geo_platform"

    # Load cohorts
    data = {}
    for gse in ["GSE63060", "GSE63061"]:
        X, y, gpl = load_geo_series(gse, raw_dir)
        p2s = load_platform_probe_to_symbol(gpl, cache_dir)
        
        # Map to gene symbols (averaging probes)
        X.columns = [p2s.get(c, c) for c in X.columns]
        X = X.T.groupby(level=0).mean().T
        data[gse] = (X, y)

    # 60 -> 61
    source_name = "GSE63060"
    target_name = "GSE63061"
    
    Xs_raw, ys = data[source_name]
    Xt_raw, yt = data[target_name]
    
    # Common genes
    common = sorted(list(set(Xs_raw.columns) & set(Xt_raw.columns)))
    Xs = Xs_raw[common]
    Xt = Xt_raw[common]
    
    # Split target-train (as in benchmark)
    Xtr_t, Xte_t, ytr, yte = train_test_split(
        Xt, yt, test_size=0.3, random_state=SEED, stratify=yt
    )
    
    # DE Scores (on target-train)
    def get_de_stats(X, y):
        ad_mask = (y == 1)
        ctl_mask = (y == 0)
        t_stats, _ = ttest_ind(X[ad_mask], X[ctl_mask], axis=0)
        return np.abs(t_stats)

    de_scores = get_de_stats(Xtr_t, ytr)
    de_series = pd.Series(de_scores, index=common)
    
    # Distortion Scores (on pool)
    X_pool = pd.concat([Xs, Xtr_t], axis=0)
    batch_combo = pd.Series([source_name]*len(Xs) + [target_name]*len(Xtr_t), index=X_pool.index)
    
    covars = pd.DataFrame({"batch": batch_combo.values}, index=X_pool.index)
    fit = neuroCombat(dat=X_pool.T, covars=covars, batch_col="batch")
    
    gamma_star = fit["estimates"]["gamma.star"]
    delta_star = fit["estimates"]["delta.star"]
    gamma_abs = np.mean(np.abs(gamma_star), axis=0)
    delta_abs = np.mean(np.abs(1.0 - delta_star), axis=0)
    
    z_gamma = zscore(gamma_abs)
    z_delta = zscore(delta_abs)
    
    # v2 Asymmetric Score
    alpha = 0.2
    d_v2 = alpha * np.abs(z_gamma) + (1 - alpha) * np.abs(z_delta)
    
    # GMM Threshold
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(d_v2.reshape(-1, 1))
    means = gmm.means_.flatten()
    tau = np.mean(means)
    
    # Identify "Dropped High-DE" genes
    dist_series = pd.Series(d_v2, index=common)
    
    # Top 500 DE genes by t-stat magnitude
    top_de_genes = de_series.sort_values(ascending=False).head(500).index
    dropped_candidates = [g for g in top_de_genes if dist_series[g] > tau]
    
    # Sort by DE significance
    result_df = pd.DataFrame({
        "gene": dropped_candidates,
        "de_score": de_series[dropped_candidates],
        "distortion_score": dist_series[dropped_candidates]
    }).sort_values(by="de_score", ascending=False)
    
    os.makedirs("outputs/stats", exist_ok=True)
    result_df.to_csv("outputs/stats/bdpfs_dropped_high_de_genes_60to61.csv", index=False)
    print(f"--- Identified {len(dropped_candidates)} dropped high-DE genes. Results saved to outputs/stats/bdpfs_dropped_high_de_genes_60to61.csv ---")

if __name__ == "__main__":
    main()
