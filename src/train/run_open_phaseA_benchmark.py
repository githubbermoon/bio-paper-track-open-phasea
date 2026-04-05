from __future__ import annotations

import gzip
import io
import json
import re
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from pycombat import Combat
from scipy.stats import ttest_ind
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


SEED = 42


def _series_url(gse: str) -> str:
    prefix = f"{gse[:-3]}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/{gse}_series_matrix.txt.gz"


def _download_text(url: str, cache_path: Path) -> str:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if cache_path.exists():
        return gzip.decompress(cache_path.read_bytes()).decode("utf-8", errors="ignore")
    with urlopen(url, timeout=60) as resp:
        data = resp.read()
    cache_path.write_bytes(data)
    return gzip.decompress(data).decode("utf-8", errors="ignore")


def _extract_label(sample_chars: list[str]) -> int | None:
    s = " | ".join(sample_chars).lower()

    status_hits = re.findall(r"status\s*:\s*([a-z0-9 _-]+)", s)
    for status in status_hits:
        st = status.strip()
        if st == "ad" or "alzheimer" in st:
            return 1
        if st in {"ctl", "control", "cn", "normal"}:
            return 0
        if "ctl to ad" in st or "mci" in st or "borderline" in st or "other" in st:
            return None

    if any(k in s for k in ["alzheimer", "ad case", "dementia"]):
        return 1
    if any(k in s for k in ["control", "cognitively normal", "healthy", "status: ctl"]):
        return 0

    m = re.search(r"diagnosis[:=]\s*([a-z0-9 _-]+)", s)
    if m:
        d = m.group(1)
        if "alzheimer" in d or d.strip() == "ad":
            return 1
        if "control" in d or "normal" in d or d.strip() == "ctl":
            return 0
    return None


def load_geo_series(gse: str, raw_dir: Path) -> tuple[pd.DataFrame, pd.Series]:
    txt = _download_text(_series_url(gse), raw_dir / f"{gse}_series_matrix.txt.gz")
    lines = txt.splitlines()

    sample_geo = []
    sample_chars_rows: list[list[str]] = []
    table_start = None
    table_end = None

    for i, line in enumerate(lines):
        if line.startswith("!Sample_geo_accession"):
            sample_geo = [x.strip().strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            row = [x.strip().strip('"') for x in line.split("\t")[1:]]
            sample_chars_rows.append(row)
        elif line.strip() == "!series_matrix_table_begin":
            table_start = i + 1
        elif line.strip() == "!series_matrix_table_end":
            table_end = i
            break

    if table_start is None or table_end is None:
        raise RuntimeError(f"Could not parse table boundaries for {gse}")

    char_by_sample: list[list[str]] = [[] for _ in sample_geo]
    for row in sample_chars_rows:
        for j, val in enumerate(row):
            if j < len(char_by_sample):
                char_by_sample[j].append(val)

    labels = [_extract_label(ch) for ch in char_by_sample]

    tab = "\n".join(lines[table_start:table_end])
    expr = pd.read_csv(io.StringIO(tab), sep="\t")
    expr.rename(columns={expr.columns[0]: "ID_REF"}, inplace=True)

    value_cols = [c for c in expr.columns if c in sample_geo]
    expr = expr[["ID_REF", *value_cols]].copy()
    expr["ID_REF"] = expr["ID_REF"].astype(str)

    x = expr.set_index("ID_REF").T
    x.index.name = "gsm"
    x = x.apply(pd.to_numeric, errors="coerce")

    y = pd.Series(labels, index=sample_geo, name="label")
    keep = y.notna() & y.index.isin(x.index)
    x = x.loc[keep.index[keep]].copy()
    y = y.loc[keep.index[keep]].astype(int)

    return x, y


def _metrics(y_true: np.ndarray, prob: np.ndarray) -> dict[str, float]:
    pred = (prob >= 0.5).astype(int)
    return {
        "auroc": float(roc_auc_score(y_true, prob)),
        "auprc": float(average_precision_score(y_true, prob)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, pred)),
        "brier": float(brier_score_loss(y_true, prob)),
    }


def _fit_predict(X_train: pd.DataFrame, y_train: pd.Series, X_test: pd.DataFrame) -> np.ndarray:
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
    model.fit(X_train, y_train)
    return model.predict_proba(X_test)[:, 1]


def _common_genes(Xa: pd.DataFrame, Xb: pd.DataFrame) -> pd.Index:
    common = Xa.columns.intersection(Xb.columns)
    if len(common) == 0:
        raise RuntimeError("No common genes across cohorts")
    return common


def _select_top_genes(
    X_train_target: pd.DataFrame,
    y_train_target: pd.Series,
    common: pd.Index,
    top_n: int,
    feature_mode: str,
) -> list[str]:
    Xt = X_train_target[common]

    if feature_mode == "var":
        score = Xt.var(axis=0)
    elif feature_mode == "de_ttest":
        cases = Xt.loc[y_train_target == 1]
        ctrls = Xt.loc[y_train_target == 0]
        stat, _ = ttest_ind(cases.values, ctrls.values, axis=0, equal_var=False, nan_policy="omit")
        score = pd.Series(np.abs(np.nan_to_num(stat, nan=0.0)), index=Xt.columns)
    else:
        raise ValueError(f"Unknown feature_mode: {feature_mode}")

    score = score.replace([np.inf, -np.inf], 0.0).fillna(0.0)
    return score.sort_values(ascending=False).head(min(top_n, len(score))).index.tolist()


def _combat_train_test(
    X_train: pd.DataFrame,
    batch_train: pd.Series,
    X_test: pd.DataFrame,
    batch_test: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    ComBat adjustment over stacked train+test expression (feature-only, no labels).
    This is a transductive harmonization step to align batch distributions.
    """
    combat = Combat()
    X_all = pd.concat([X_train, X_test], axis=0)
    b_all = pd.concat([batch_train, batch_test], axis=0)
    X_all_adj = combat.fit_transform(X_all.values, b=b_all.values)
    X_all_df = pd.DataFrame(X_all_adj, index=X_all.index, columns=X_all.columns)
    Xtr_df = X_all_df.loc[X_train.index]
    Xte_df = X_all_df.loc[X_test.index]
    return Xtr_df, Xte_df


def _null_avg_predict(
    Xtr: pd.DataFrame,
    ytr: pd.Series,
    Xte: pd.DataFrame,
    yte: pd.Series,
    n_perm: int = 100,
) -> tuple[np.ndarray, dict]:
    rng = np.random.default_rng(SEED)
    prob_acc = np.zeros(len(Xte), dtype=float)
    aucs: list[float] = []

    for _ in range(n_perm):
        y_perm = ytr.iloc[rng.permutation(len(ytr))].copy()
        y_perm.index = ytr.index
        p = _fit_predict(Xtr, y_perm, Xte)
        prob_acc += p
        aucs.append(float(roc_auc_score(yte.values, p)))

    p_avg = prob_acc / n_perm
    auc_arr = np.array(aucs, dtype=float)
    return p_avg, {
        "n_perm": int(n_perm),
        "perm_auroc_mean": float(np.mean(auc_arr)),
        "perm_auroc_std": float(np.std(auc_arr)),
        "perm_auroc_q05": float(np.quantile(auc_arr, 0.05)),
        "perm_auroc_q95": float(np.quantile(auc_arr, 0.95)),
    }


def run() -> None:
    root = Path(__file__).resolve().parents[2]
    raw_dir = root / "data/raw/open_geo"
    out_metrics = root / "outputs/metrics/open_phaseA_main_results.csv"
    out_preds = root / "outputs/metrics/open_phaseA_predictions.csv"
    out_stats = root / "outputs/stats/open_phaseA_stats.json"
    out_manifest = root / "outputs/open_phaseA_data_manifest.json"

    gses = ["GSE63060", "GSE63061"]
    data = {}
    for gse in gses:
        X, y = load_geo_series(gse, raw_dir)
        data[gse] = (X, y)

    rows = []
    pred_rows = []
    stats: dict[str, dict[str, float | int | str]] = {}

    for feature_mode in ["var", "de_ttest"]:
        for top_n in [200, 1000]:
            for source, target in [("GSE63060", "GSE63061"), ("GSE63061", "GSE63060")]:
                Xs, ys = data[source]
                Xt, yt = data[target]

                common = _common_genes(Xs, Xt)
                Xtr_t, Xte_t, ytr, yte = train_test_split(
                    Xt[common], yt, test_size=0.3, random_state=SEED, stratify=yt
                )

                genes = _select_top_genes(Xtr_t, ytr, common, top_n=top_n, feature_mode=feature_mode)
                Xs_use = Xs[genes]
                Xtr = Xtr_t[genes]
                Xte = Xte_t[genes]

                p_target = _fit_predict(Xtr, ytr, Xte)
                m_target = _metrics(yte.values, p_target)
                rows.append(
                    {
                        "source": source,
                        "target": target,
                        "feature_mode": feature_mode,
                        "top_n_genes": top_n,
                        "arm": "target_only",
                        **m_target,
                    }
                )
                for sid, yt_i, pr_i in zip(Xte.index, yte.values, p_target):
                    pred_rows.append(
                        {
                            "source": source,
                            "target": target,
                            "feature_mode": feature_mode,
                            "top_n_genes": top_n,
                            "arm": "target_only",
                            "sample_id": sid,
                            "y_true": int(yt_i),
                            "y_prob": float(pr_i),
                        }
                    )

                p_source = _fit_predict(Xs_use, ys, Xte)
                m_source = _metrics(yte.values, p_source)
                rows.append(
                    {
                        "source": source,
                        "target": target,
                        "feature_mode": feature_mode,
                        "top_n_genes": top_n,
                        "arm": "source_only",
                        **m_source,
                    }
                )
                for sid, yt_i, pr_i in zip(Xte.index, yte.values, p_source):
                    pred_rows.append(
                        {
                            "source": source,
                            "target": target,
                            "feature_mode": feature_mode,
                            "top_n_genes": top_n,
                            "arm": "source_only",
                            "sample_id": sid,
                            "y_true": int(yt_i),
                            "y_prob": float(pr_i),
                        }
                    )

                X_combo_raw = pd.concat([Xs_use, Xtr], axis=0)
                y_combo_raw = pd.concat([ys, ytr], axis=0)
                p_combo_raw = _fit_predict(X_combo_raw, y_combo_raw, Xte)
                m_combo_raw = _metrics(yte.values, p_combo_raw)
                rows.append(
                    {
                        "source": source,
                        "target": target,
                        "feature_mode": feature_mode,
                        "top_n_genes": top_n,
                        "arm": "source_plus_target_raw",
                        **m_combo_raw,
                    }
                )
                for sid, yt_i, pr_i in zip(Xte.index, yte.values, p_combo_raw):
                    pred_rows.append(
                        {
                            "source": source,
                            "target": target,
                            "feature_mode": feature_mode,
                            "top_n_genes": top_n,
                            "arm": "source_plus_target_raw",
                            "sample_id": sid,
                            "y_true": int(yt_i),
                            "y_prob": float(pr_i),
                        }
                    )

                batch_combo = pd.Series([source] * len(Xs_use) + [target] * len(Xtr), index=X_combo_raw.index)
                batch_te = pd.Series([target] * len(Xte), index=Xte.index)
                X_combo_cb, Xte_cb = _combat_train_test(X_combo_raw, batch_combo, Xte, batch_te)
                p_combo_cb = _fit_predict(X_combo_cb, y_combo_raw, Xte_cb)
                m_combo_cb = _metrics(yte.values, p_combo_cb)
                rows.append(
                    {
                        "source": source,
                        "target": target,
                        "feature_mode": feature_mode,
                        "top_n_genes": top_n,
                        "arm": "source_plus_target_combat",
                        **m_combo_cb,
                    }
                )
                # Keep backward-compatible arm label used by downstream stats/paper
                rows.append(
                    {
                        "source": source,
                        "target": target,
                        "feature_mode": feature_mode,
                        "top_n_genes": top_n,
                        "arm": "source_plus_target",
                        **m_combo_cb,
                    }
                )
                for sid, yt_i, pr_i in zip(Xte.index, yte.values, p_combo_cb):
                    pred_rows.append(
                        {
                            "source": source,
                            "target": target,
                            "feature_mode": feature_mode,
                            "top_n_genes": top_n,
                            "arm": "source_plus_target_combat",
                            "sample_id": sid,
                            "y_true": int(yt_i),
                            "y_prob": float(pr_i),
                        }
                    )
                    pred_rows.append(
                        {
                            "source": source,
                            "target": target,
                            "feature_mode": feature_mode,
                            "top_n_genes": top_n,
                            "arm": "source_plus_target",
                            "sample_id": sid,
                            "y_true": int(yt_i),
                            "y_prob": float(pr_i),
                        }
                    )

                p_null_avg, null_meta = _null_avg_predict(Xtr, ytr, Xte, yte, n_perm=100)
                m_null = _metrics(yte.values, p_null_avg)
                rows.append(
                    {
                        "source": source,
                        "target": target,
                        "feature_mode": feature_mode,
                        "top_n_genes": top_n,
                        "arm": "null_label_permutation_avg100",
                        **m_null,
                    }
                )
                rows.append(
                    {
                        "source": source,
                        "target": target,
                        "feature_mode": feature_mode,
                        "top_n_genes": top_n,
                        "arm": "null_label_permutation",
                        **m_null,
                    }
                )
                for sid, yt_i, pr_i in zip(Xte.index, yte.values, p_null_avg):
                    pred_rows.append(
                        {
                            "source": source,
                            "target": target,
                            "feature_mode": feature_mode,
                            "top_n_genes": top_n,
                            "arm": "null_label_permutation_avg100",
                            "sample_id": sid,
                            "y_true": int(yt_i),
                            "y_prob": float(pr_i),
                        }
                    )
                    pred_rows.append(
                        {
                            "source": source,
                            "target": target,
                            "feature_mode": feature_mode,
                            "top_n_genes": top_n,
                            "arm": "null_label_permutation",
                            "sample_id": sid,
                            "y_true": int(yt_i),
                            "y_prob": float(pr_i),
                        }
                    )

                key = f"{feature_mode}__{source}_to_{target}_top{top_n}"
                stats[key] = {
                    "feature_mode": feature_mode,
                    "source": source,
                    "target": target,
                    "top_n_genes": int(top_n),
                    "delta_auroc_source_plus_target_combat_vs_target_only": float(
                        m_combo_cb["auroc"] - m_target["auroc"]
                    ),
                    "delta_auroc_source_plus_target_raw_vs_target_only": float(
                        m_combo_raw["auroc"] - m_target["auroc"]
                    ),
                    "delta_auroc_target_only_vs_null_avg100": float(m_target["auroc"] - m_null["auroc"]),
                    **null_meta,
                }

    out_metrics.parent.mkdir(parents=True, exist_ok=True)
    out_stats.parent.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(rows).to_csv(out_metrics, index=False)
    pd.DataFrame(pred_rows).to_csv(out_preds, index=False)
    out_stats.write_text(json.dumps(stats, indent=2), encoding="utf-8")

    manifest = {
        "datasets": gses,
        "source": "NCBI GEO series_matrix",
        "urls": {g: _series_url(g) for g in gses},
        "cached_files": [str(p) for p in sorted(raw_dir.glob("*.gz"))],
        "feature_modes": ["var", "de_ttest"],
        "null_policy": "100x label permutations averaged at prediction level",
        "batch_harmonization": "ComBat (pycombat) applied in source_plus_target_combat arm",
    }
    out_manifest.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    run()
