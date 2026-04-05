from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, brier_score_loss, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


def _safe_metrics(y_true: np.ndarray, y_prob: np.ndarray) -> dict[str, float]:
    y_pred = (y_prob >= 0.5).astype(int)
    out = {
        "auroc": float("nan"),
        "auprc": float("nan"),
        "balanced_accuracy": float("nan"),
        "brier": float("nan"),
    }
    try:
        out["auroc"] = float(roc_auc_score(y_true, y_prob))
    except Exception:
        pass
    try:
        out["auprc"] = float(average_precision_score(y_true, y_prob))
    except Exception:
        pass
    try:
        out["balanced_accuracy"] = float(balanced_accuracy_score(y_true, y_pred))
    except Exception:
        pass
    try:
        out["brier"] = float(brier_score_loss(y_true, y_prob))
    except Exception:
        pass
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--target", default="dementia_dx")
    ap.add_argument("--time-col", default="collection_time")
    ap.add_argument("--out", default="outputs/metrics/main_benchmark_v1.csv")
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    df[args.time_col] = pd.to_datetime(df[args.time_col], errors="coerce")
    df = df.sort_values(args.time_col).reset_index(drop=True)

    split1 = int(len(df) * 0.7)
    split2 = int(len(df) * 0.85)
    train = df.iloc[:split1].copy()
    test = df.iloc[split2:].copy()

    y_train = train[args.target].astype(int).values
    y_test = test[args.target].astype(int).values

    feature_cols = [c for c in df.columns if c not in {args.target, args.time_col}]
    num_cols = [c for c in feature_cols if pd.api.types.is_numeric_dtype(df[c])]
    cat_cols = [c for c in feature_cols if c not in num_cols]

    pre = ColumnTransformer(
        [
            ("num", Pipeline([("imp", SimpleImputer(strategy="median")), ("sc", StandardScaler())]), num_cols),
            ("cat", Pipeline([("imp", SimpleImputer(strategy="most_frequent")), ("oh", OneHotEncoder(handle_unknown="ignore"))]), cat_cols),
        ],
        remainder="drop",
    )

    rows: list[dict[str, object]] = []

    lr = Pipeline([("pre", pre), ("model", LogisticRegression(max_iter=200, class_weight="balanced"))])
    lr.fit(train[feature_cols], y_train)
    lr_prob = lr.predict_proba(test[feature_cols])[:, 1]
    m = _safe_metrics(y_test, lr_prob)
    rows.append({"model": "logistic_l2", **m, "status": "ok"})

    try:
        from lightgbm import LGBMClassifier  # type: ignore

        lgbm = Pipeline(
            [
                ("pre", pre),
                (
                    "model",
                    LGBMClassifier(
                        n_estimators=200,
                        learning_rate=0.05,
                        num_leaves=31,
                        objective="binary",
                        random_state=42,
                    ),
                ),
            ]
        )
        lgbm.fit(train[feature_cols], y_train)
        lgbm_prob = lgbm.predict_proba(test[feature_cols])[:, 1]
        m2 = _safe_metrics(y_test, lgbm_prob)
        rows.append({"model": "lightgbm", **m2, "status": "ok"})
    except Exception:
        rows.append({"model": "lightgbm", "auroc": np.nan, "auprc": np.nan, "balanced_accuracy": np.nan, "brier": np.nan, "status": "not_available"})

    out = pd.DataFrame(rows)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)


if __name__ == "__main__":
    main()
