from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC

from src.train.run_open_phaseA_benchmark import SEED, _common_genes, load_geo_series


def _select_de_top(Xtr_t: pd.DataFrame, ytr: pd.Series, common: pd.Index, top_n: int = 1000) -> list[str]:
    Xt = Xtr_t[common]
    cases = Xt.loc[ytr == 1]
    ctrls = Xt.loc[ytr == 0]
    stat, _ = ttest_ind(cases.values, ctrls.values, axis=0, equal_var=False, nan_policy="omit")
    score = pd.Series(np.abs(np.nan_to_num(stat, nan=0.0)), index=Xt.columns)
    return score.sort_values(ascending=False).head(min(top_n, len(score))).index.tolist()


def _fit_predict_scores(model_name: str, Xtr: pd.DataFrame, ytr: pd.Series, Xte: pd.DataFrame) -> np.ndarray:
    if model_name == "logistic_regression":
        model = Pipeline(
            [
                ("imp", SimpleImputer(strategy="median")),
                ("sc", StandardScaler(with_mean=False)),
                (
                    "clf",
                    LogisticRegression(
                        max_iter=2000,
                        class_weight="balanced",
                        solver="liblinear",
                        random_state=SEED,
                    ),
                ),
            ]
        )
        model.fit(Xtr, ytr)
        return model.predict_proba(Xte)[:, 1]

    if model_name == "linear_svm":
        model = Pipeline(
            [
                ("imp", SimpleImputer(strategy="median")),
                ("sc", StandardScaler(with_mean=False)),
                ("clf", LinearSVC(class_weight="balanced", random_state=SEED, max_iter=6000)),
            ]
        )
        model.fit(Xtr, ytr)
        return model.decision_function(Xte)

    if model_name == "random_forest":
        imp = SimpleImputer(strategy="median")
        Xtr_i = imp.fit_transform(Xtr)
        Xte_i = imp.transform(Xte)
        clf = RandomForestClassifier(
            n_estimators=500,
            random_state=SEED,
            class_weight="balanced",
            max_features="sqrt",
            n_jobs=-1,
        )
        clf.fit(Xtr_i, ytr)
        return clf.predict_proba(Xte_i)[:, 1]

    raise ValueError(model_name)


def run() -> None:
    root = Path(__file__).resolve().parents[2]
    raw_dir = root / "data/raw/open_geo"
    out = root / "outputs/stats/open_phaseA_model_family_sensitivity.csv"

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

        for m in ["logistic_regression", "linear_svm", "random_forest"]:
            s = _fit_predict_scores(m, Xtr, ytr, Xte)
            rows.append(
                {
                    "source": source,
                    "target": target,
                    "feature_mode": "de_ttest",
                    "top_n_genes": 1000,
                    "model": m,
                    "auroc": float(roc_auc_score(yte, s)),
                    "auprc": float(average_precision_score(yte, s)),
                }
            )

    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)


if __name__ == "__main__":
    run()
