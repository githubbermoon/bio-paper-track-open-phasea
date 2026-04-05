from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

from src.train.run_open_phaseA_benchmark import SEED, _common_genes, _fit_predict, load_geo_series


def _select_de_top(Xtr_t: pd.DataFrame, ytr: pd.Series, common: pd.Index, top_n: int = 1000) -> list[str]:
    Xt = Xtr_t[common]
    cases = Xt.loc[ytr == 1]
    ctrls = Xt.loc[ytr == 0]
    stat, _ = ttest_ind(cases.values, ctrls.values, axis=0, equal_var=False, nan_policy="omit")
    score = pd.Series(np.abs(np.nan_to_num(stat, nan=0.0)), index=Xt.columns)
    return score.sort_values(ascending=False).head(min(top_n, len(score))).index.tolist()


def _perm_auroc_distribution(Xtr: pd.DataFrame, ytr: pd.Series, Xte: pd.DataFrame, yte: pd.Series, n_perm: int = 1000) -> np.ndarray:
    rng = np.random.default_rng(SEED)
    vals = []
    for _ in range(n_perm):
        y_perm = ytr.iloc[rng.permutation(len(ytr))].copy()
        y_perm.index = ytr.index
        p = _fit_predict(Xtr, y_perm, Xte)
        vals.append(float(roc_auc_score(yte.values, p)))
    return np.array(vals, dtype=float)


def run() -> None:
    root = Path(__file__).resolve().parents[2]
    raw_dir = root / "data/raw/open_geo"
    out = root / "outputs/stats/open_phaseA_null_stability_de1000_perm1000.csv"

    data = {}
    for gse in ["GSE63060", "GSE63061"]:
        X, y, _ = load_geo_series(gse, raw_dir)
        data[gse] = (X, y)

    rows = []
    for source, target in [("GSE63060", "GSE63061"), ("GSE63061", "GSE63060")]:
        Xs, _ = data[source]
        Xt, yt = data[target]
        common = _common_genes(Xs, Xt)
        Xtr_t, Xte_t, ytr, yte = train_test_split(Xt[common], yt, test_size=0.3, random_state=SEED, stratify=yt)
        genes = _select_de_top(Xtr_t, ytr, common, top_n=1000)
        Xtr = Xtr_t[genes]
        Xte = Xte_t[genes]

        arr = _perm_auroc_distribution(Xtr, ytr, Xte, yte, n_perm=1000)
        rows.append(
            {
                "source": source,
                "target": target,
                "feature_mode": "de_ttest",
                "top_n_genes": 1000,
                "n_perm": 1000,
                "null_perm_auroc_mean": float(np.mean(arr)),
                "null_perm_auroc_std": float(np.std(arr)),
                "null_perm_auroc_q01": float(np.quantile(arr, 0.01)),
                "null_perm_auroc_q05": float(np.quantile(arr, 0.05)),
                "null_perm_auroc_q95": float(np.quantile(arr, 0.95)),
                "null_perm_auroc_q99": float(np.quantile(arr, 0.99)),
            }
        )

    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)


if __name__ == "__main__":
    run()
