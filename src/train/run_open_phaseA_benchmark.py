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


def _platform_annot_url(gpl: str) -> str:
    num = int(gpl.replace("GPL", ""))
    prefix = f"GPL{str(num)[:-3]}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{prefix}/{gpl}/annot/{gpl}.annot.gz"


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


def load_geo_series(gse: str, raw_dir: Path) -> tuple[pd.DataFrame, pd.Series, str]:
    txt = _download_text(_series_url(gse), raw_dir / f"{gse}_series_matrix.txt.gz")
    lines = txt.splitlines()

    sample_geo = []
    sample_chars_rows: list[list[str]] = []
    table_start = None
    table_end = None
    gpl = ""

    for i, line in enumerate(lines):
        if line.startswith("!Series_platform_id"):
            parts = line.split("\t")
            if len(parts) > 1:
                gpl = parts[1].strip().strip('"')
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

    return x, y, gpl


def load_platform_probe_to_symbol(gpl: str, cache_dir: Path) -> dict[str, str]:
    txt = _download_text(_platform_annot_url(gpl), cache_dir / f"{gpl}.annot.gz")
    lines = txt.splitlines()

    start = None
    end = None
    for i, line in enumerate(lines):
        if line.strip() == "!platform_table_begin":
            start = i + 1
        elif line.strip() == "!platform_table_end":
            end = i
            break
    if start is None or end is None:
        return {}

    df = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", dtype=str)
    cols = {c.lower(): c for c in df.columns}

    id_col = None
    for c in ["id", "id_ref", "probe_id", "ilmnid", "name"]:
        if c in cols:
            id_col = cols[c]
            break
    symbol_col = None
    for c in ["symbol", "gene symbol", "genesymbol", "gene_symbol", "symboli", "gene"]:
        if c in cols:
            symbol_col = cols[c]
            break
    if id_col is None:
        return {}

    mapping = {}
    if symbol_col is None:
        return mapping

    for _, r in df[[id_col, symbol_col]].dropna().iterrows():
        probe = str(r[id_col]).strip()
        sym = str(r[symbol_col]).strip()
        if not probe or not sym:
            continue
        sym = re.split(r"[;|,/ ]+", sym)[0].strip().upper()
        if sym and sym not in {"NA", "N/A", "---", "NULL", "NONE"}:
            mapping[probe] = sym
    return mapping


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


def _feature_scores(X_train_target: pd.DataFrame, y_train_target: pd.Series, common: pd.Index) -> tuple[pd.Series, pd.Series]:
    Xt = X_train_target[common]
    var_score = Xt.var(axis=0).replace([np.inf, -np.inf], 0.0).fillna(0.0)

    cases = Xt.loc[y_train_target == 1]
    ctrls = Xt.loc[y_train_target == 0]
    stat, _ = ttest_ind(cases.values, ctrls.values, axis=0, equal_var=False, nan_policy="omit")
    de_score = pd.Series(np.abs(np.nan_to_num(stat, nan=0.0)), index=Xt.columns).replace([np.inf, -np.inf], 0.0).fillna(0.0)
    return var_score, de_score


def _select_top_genes(
    mode: str,
    top_n: int,
    common: pd.Index,
    var_score: pd.Series,
    de_score: pd.Series,
    agora_probe_set: set[str],
) -> list[str]:
    if mode == "var":
        s = var_score
    elif mode == "de_ttest":
        s = de_score
    elif mode == "agora_only":
        s = var_score[var_score.index.isin(agora_probe_set)]
    elif mode == "de_agora_intersection":
        s = de_score[de_score.index.isin(agora_probe_set)]
    else:
        raise ValueError(f"Unknown feature_mode: {mode}")

    s = s.reindex(common.intersection(s.index)).dropna()
    if len(s) == 0:
        return []
    return s.sort_values(ascending=False).head(min(top_n, len(s))).index.tolist()


def _combat_transductive_stacked(
    X_train: pd.DataFrame,
    batch_train: pd.Series,
    X_test: pd.DataFrame,
    batch_test: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Sensitivity-only ComBat: fit on stacked train+test features (no labels).
    Kept as transductive analysis; NOT used for leakage-safe primary claims.
    """
    combat = Combat()
    X_all = pd.concat([X_train, X_test], axis=0)
    b_all = pd.concat([batch_train, batch_test], axis=0)
    X_all_adj = combat.fit_transform(X_all.values, b=b_all.values)
    X_all_df = pd.DataFrame(X_all_adj, index=X_all.index, columns=X_all.columns)
    return X_all_df.loc[X_train.index], X_all_df.loc[X_test.index]


def _null_distribution(
    Xtr: pd.DataFrame,
    ytr: pd.Series,
    Xte: pd.DataFrame,
    yte: pd.Series,
    n_perm: int = 100,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(SEED)
    aucs = []
    prob_acc = np.zeros(len(Xte), dtype=float)

    for _ in range(n_perm):
        y_perm = ytr.iloc[rng.permutation(len(ytr))].copy()
        y_perm.index = ytr.index
        p = _fit_predict(Xtr, y_perm, Xte)
        aucs.append(float(roc_auc_score(yte.values, p)))
        prob_acc += p

    return np.array(aucs, dtype=float), (prob_acc / n_perm)


def run() -> None:
    root = Path(__file__).resolve().parents[2]
    raw_dir = root / "data/raw/open_geo"
    annot_cache = root / "data/raw/open_geo_platform"

    out_metrics = root / "outputs/metrics/open_phaseA_main_results.csv"
    out_preds = root / "outputs/metrics/open_phaseA_predictions.csv"
    out_stats = root / "outputs/stats/open_phaseA_stats.json"
    out_null = root / "outputs/stats/open_phaseA_null_distribution.csv"
    out_manifest = root / "outputs/open_phaseA_data_manifest.json"

    gses = ["GSE63060", "GSE63061"]
    data = {}
    gpl_of = {}
    for gse in gses:
        X, y, gpl = load_geo_series(gse, raw_dir)
        data[gse] = (X, y)
        gpl_of[gse] = gpl

    probe_to_symbol = {}
    for gpl in sorted(set(gpl_of.values())):
        probe_to_symbol.update(load_platform_probe_to_symbol(gpl, annot_cache))

    agora_df = pd.read_csv(root / "outputs/data/ampad_open_nominated_targets.csv")
    agora_symbols = set(agora_df["hgnc_symbol"].dropna().astype(str).str.upper())

    rows = []
    pred_rows = []
    null_rows = []
    stats: dict[str, dict[str, float | int | str]] = {}

    feature_modes = ["var", "de_ttest", "agora_only", "de_agora_intersection"]

    for feature_mode in feature_modes:
        for top_n in [200, 1000]:
            for source, target in [("GSE63060", "GSE63061"), ("GSE63061", "GSE63060")]:
                Xs, ys = data[source]
                Xt, yt = data[target]

                common = _common_genes(Xs, Xt)
                Xtr_t, Xte_t, ytr, yte = train_test_split(
                    Xt[common], yt, test_size=0.3, random_state=SEED, stratify=yt
                )

                var_score, de_score = _feature_scores(Xtr_t, ytr, common)
                agora_probe_set = {
                    g for g in common if probe_to_symbol.get(g, "").upper() in agora_symbols
                }

                genes = _select_top_genes(feature_mode, top_n, common, var_score, de_score, agora_probe_set)
                if len(genes) < 20:
                    # skip empty/degenerate settings safely
                    continue

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
                        "selected_gene_count": len(genes),
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
                        "selected_gene_count": len(genes),
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

                # pooled raw
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
                        "selected_gene_count": len(genes),
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

                # pooled combat leakage-safe
                batch_combo = pd.Series([source] * len(Xs_use) + [target] * len(Xtr), index=X_combo_raw.index)
                batch_te = pd.Series([target] * len(Xte), index=Xte.index)
                X_combo_cb, Xte_cb = _combat_transductive_stacked(X_combo_raw, batch_combo, Xte, batch_te)
                p_combo_cb = _fit_predict(X_combo_cb, y_combo_raw, Xte_cb)
                m_combo_cb = _metrics(yte.values, p_combo_cb)
                rows.append(
                    {
                        "source": source,
                        "target": target,
                        "feature_mode": feature_mode,
                        "top_n_genes": top_n,
                        "selected_gene_count": len(genes),
                        "arm": "source_plus_target_combat_transductive",
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
                            "arm": "source_plus_target_combat_transductive",
                            "sample_id": sid,
                            "y_true": int(yt_i),
                            "y_prob": float(pr_i),
                        }
                    )

                # null distribution
                null_aucs, p_null_avg = _null_distribution(Xtr, ytr, Xte, yte, n_perm=100)
                null_mean = float(np.mean(null_aucs))
                null_q05 = float(np.quantile(null_aucs, 0.05))
                null_q95 = float(np.quantile(null_aucs, 0.95))

                # Row value now uses permutation-AUROC mean for consistency
                rows.append(
                    {
                        "source": source,
                        "target": target,
                        "feature_mode": feature_mode,
                        "top_n_genes": top_n,
                        "selected_gene_count": len(genes),
                        "arm": "null_label_permutation_mean_auroc",
                        "auroc": null_mean,
                        "auprc": np.nan,
                        "balanced_accuracy": np.nan,
                        "brier": np.nan,
                    }
                )

                # keep probability-averaged null predictions for optional paired sensitivity analysis
                rows.append(
                    {
                        "source": source,
                        "target": target,
                        "feature_mode": feature_mode,
                        "top_n_genes": top_n,
                        "selected_gene_count": len(genes),
                        "arm": "null_label_permutation_avg100_prob",
                        **_metrics(yte.values, p_null_avg),
                    }
                )
                for sid, yt_i, pr_i in zip(Xte.index, yte.values, p_null_avg):
                    pred_rows.append(
                        {
                            "source": source,
                            "target": target,
                            "feature_mode": feature_mode,
                            "top_n_genes": top_n,
                            "arm": "null_label_permutation_avg100_prob",
                            "sample_id": sid,
                            "y_true": int(yt_i),
                            "y_prob": float(pr_i),
                        }
                    )

                for i, auc in enumerate(null_aucs):
                    null_rows.append(
                        {
                            "source": source,
                            "target": target,
                            "feature_mode": feature_mode,
                            "top_n_genes": top_n,
                            "perm_index": i,
                            "null_perm_auroc": float(auc),
                        }
                    )

                key = f"{feature_mode}__{source}_to_{target}_top{top_n}"
                p_target_vs_null = float((null_aucs >= m_target["auroc"]).mean())
                stats[key] = {
                    "feature_mode": feature_mode,
                    "source": source,
                    "target": target,
                    "top_n_genes": int(top_n),
                    "selected_gene_count": int(len(genes)),
                    "target_auroc": float(m_target["auroc"]),
                    "delta_auroc_source_plus_target_raw_vs_target_only": float(m_combo_raw["auroc"] - m_target["auroc"]),
                    "delta_auroc_source_plus_target_combat_transductive_vs_target_only": float(
                        m_combo_cb["auroc"] - m_target["auroc"]
                    ),
                    "null_perm_n": int(len(null_aucs)),
                    "null_perm_auroc_mean": null_mean,
                    "null_perm_auroc_std": float(np.std(null_aucs)),
                    "null_perm_auroc_q05": null_q05,
                    "null_perm_auroc_q95": null_q95,
                    "p_target_gt_null_perm": p_target_vs_null,
                }

    out_metrics.parent.mkdir(parents=True, exist_ok=True)
    out_stats.parent.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(rows).to_csv(out_metrics, index=False)
    pd.DataFrame(pred_rows).to_csv(out_preds, index=False)
    pd.DataFrame(null_rows).to_csv(out_null, index=False)
    out_stats.write_text(json.dumps(stats, indent=2), encoding="utf-8")

    manifest = {
        "datasets": gses,
        "platforms": gpl_of,
        "source": "NCBI GEO series_matrix",
        "urls": {g: _series_url(g) for g in gses},
        "cached_files": [str(p) for p in sorted(raw_dir.glob("*.gz"))],
        "feature_modes": feature_modes,
        "agora_unique_symbols": int(len(agora_symbols)),
        "null_policy": "100x label permutations; primary null reported as mean permutation AUROC",
        "batch_harmonization": "ComBat transductive sensitivity arm on stacked train+test features (no labels); primary leakage-safe claims use non-ComBat arms",
    }
    out_manifest.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    run()
