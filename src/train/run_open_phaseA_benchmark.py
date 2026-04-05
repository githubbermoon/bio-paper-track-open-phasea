from __future__ import annotations

import gzip
import io
import json
import re
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


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

    # Prefer explicit GEO status fields when present (e.g., status: ad / status: ctl)
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

    m = re.search(r"diagnosis[:=]\\s*([a-z0-9 _-]+)", s)
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

    labels = []
    for ch in char_by_sample:
        labels.append(_extract_label(ch))

    tab = "\n".join(lines[table_start:table_end])
    expr = pd.read_csv(io.StringIO(tab), sep="\t")
    expr.rename(columns={expr.columns[0]: "ID_REF"}, inplace=True)

    value_cols = [c for c in expr.columns if c in sample_geo]
    expr = expr[["ID_REF", *value_cols]].copy()
    expr["ID_REF"] = expr["ID_REF"].astype(str)

    # to sample x gene
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
    model = Pipeline([
        ("imp", SimpleImputer(strategy="median")),
        ("sc", StandardScaler(with_mean=False)),
        ("clf", LogisticRegression(max_iter=2000, class_weight="balanced", solver="liblinear")),
    ])
    model.fit(X_train, y_train)
    return model.predict_proba(X_test)[:, 1]


def _common_top_var(Xa: pd.DataFrame, Xb: pd.DataFrame, top_n: int = 1000) -> list[str]:
    common = Xa.columns.intersection(Xb.columns)
    if len(common) == 0:
        raise RuntimeError("No common genes across cohorts")
    v = pd.concat([Xa[common], Xb[common]], axis=0).var(axis=0).sort_values(ascending=False)
    return v.head(min(top_n, len(v))).index.tolist()


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
    stats: dict[str, dict[str, float]] = {}

    for top_n in [200, 1000]:
        for source, target in [("GSE63060", "GSE63061"), ("GSE63061", "GSE63060")]:
            Xs, ys = data[source]
            Xt, yt = data[target]

            genes = _common_top_var(Xs, Xt, top_n=top_n)
            Xs_use = Xs[genes]
            Xt_use = Xt[genes]

            Xtr, Xte, ytr, yte = train_test_split(
                Xt_use, yt, test_size=0.3, random_state=42, stratify=yt
            )

            # no-transfer baseline
            p_target = _fit_predict(Xtr, ytr, Xte)
            m_target = _metrics(yte.values, p_target)
            rows.append({"source": source, "target": target, "top_n_genes": top_n, "arm": "target_only", **m_target})
            for sid, yt_i, pr_i in zip(Xte.index, yte.values, p_target):
                pred_rows.append({"source": source, "target": target, "top_n_genes": top_n, "arm": "target_only", "sample_id": sid, "y_true": int(yt_i), "y_prob": float(pr_i)})

            # source only transfer
            p_source = _fit_predict(Xs_use, ys, Xte)
            m_source = _metrics(yte.values, p_source)
            rows.append({"source": source, "target": target, "top_n_genes": top_n, "arm": "source_only", **m_source})
            for sid, yt_i, pr_i in zip(Xte.index, yte.values, p_source):
                pred_rows.append({"source": source, "target": target, "top_n_genes": top_n, "arm": "source_only", "sample_id": sid, "y_true": int(yt_i), "y_prob": float(pr_i)})

            # source + target train transfer
            X_combo = pd.concat([Xs_use, Xtr], axis=0)
            y_combo = pd.concat([ys, ytr], axis=0)
            p_combo = _fit_predict(X_combo, y_combo, Xte)
            m_combo = _metrics(yte.values, p_combo)
            rows.append({"source": source, "target": target, "top_n_genes": top_n, "arm": "source_plus_target", **m_combo})
            for sid, yt_i, pr_i in zip(Xte.index, yte.values, p_combo):
                pred_rows.append({"source": source, "target": target, "top_n_genes": top_n, "arm": "source_plus_target", "sample_id": sid, "y_true": int(yt_i), "y_prob": float(pr_i)})

            # null control: permuted labels
            y_perm = ytr.sample(frac=1.0, random_state=42).reset_index(drop=True)
            y_perm.index = ytr.index
            p_null = _fit_predict(Xtr, y_perm, Xte)
            m_null = _metrics(yte.values, p_null)
            rows.append({"source": source, "target": target, "top_n_genes": top_n, "arm": "null_label_permutation", **m_null})
            for sid, yt_i, pr_i in zip(Xte.index, yte.values, p_null):
                pred_rows.append({"source": source, "target": target, "top_n_genes": top_n, "arm": "null_label_permutation", "sample_id": sid, "y_true": int(yt_i), "y_prob": float(pr_i)})

            stats[f"{source}_to_{target}_top{top_n}"] = {
                "delta_auroc_source_plus_target_vs_target_only": m_combo["auroc"] - m_target["auroc"],
                "delta_auroc_target_only_vs_null": m_target["auroc"] - m_null["auroc"],
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
    }
    out_manifest.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    run()
