from __future__ import annotations

import gzip
import io
import json
import re
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from neuroCombat import neuroCombat, neuroCombatFromTraining
from scipy.stats import ttest_ind, zscore
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, brier_score_loss, roc_auc_score
from sklearn.mixture import GaussianMixture
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
    try:
        txt = _download_text(_platform_annot_url(gpl), cache_dir / f"{gpl}.annot.gz")
    except Exception as e:
        print(f"Warning: Could not fetch annot for {gpl} ({e})")
        return {}
        
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


def load_gse97760_custom(raw_dir: Path) -> tuple[pd.DataFrame, pd.Series, str]:
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97760/suppl/GSE97760_loess.txt.gz"
    cache_path = raw_dir / "GSE97760_loess.txt.gz"
    txt = _download_text(url, cache_path)
    df = pd.read_csv(io.StringIO(txt), sep="\t")
    
    # Gene symbols extraction
    df = df[df["GeneSymbol"].notna() & (df["GeneSymbol"] != "")]
    df["GeneSymbol"] = df["GeneSymbol"].astype(str).apply(lambda x: re.split(r"[;|,/ ]+", x)[0].strip().upper())
    df = df[df["GeneSymbol"] != ""]
    
    # Extract sample columns and collapse
    sample_cols = [c for c in df.columns if c.startswith("Norm_")]
    expr = df[["GeneSymbol"] + sample_cols].groupby("GeneSymbol").mean()
    X = expr.T
    
    y_labels = []
    sample_ids = []
    for c in X.index:
        sample_ids.append(c)
        if "Control" in c:
            y_labels.append(0)
        elif "AP" in c:
            y_labels.append(1)
        else:
            y_labels.append(None)
            
    y = pd.Series(y_labels, index=sample_ids, name="label")
    keep = y.notna()
    X = X.loc[keep]
    y = y.loc[keep].astype(int)
    
    return X, y, "GPL16699"


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


def _feature_scores(X_train: pd.DataFrame, y_train: pd.Series, common: pd.Index) -> tuple[pd.Series, pd.Series]:
    Xc = X_train[common]
    var_score = Xc.var(axis=0).replace([np.inf, -np.inf], 0.0).fillna(0.0)

    cases = Xc.loc[y_train == 1]
    ctrls = Xc.loc[y_train == 0]
    stat, _ = ttest_ind(cases.values, ctrls.values, axis=0, equal_var=False, nan_policy="omit")
    de_score = pd.Series(np.abs(np.nan_to_num(stat, nan=0.0)), index=Xc.columns).replace([np.inf, -np.inf], 0.0).fillna(0.0)
    return var_score, de_score


def _select_top_genes(
    mode: str,
    top_n: int,
    common: pd.Index,
    var_score: pd.Series,
    de_score: pd.Series,
    agora_probe_set: set[str],
    robust_features: pd.Index = None,
    v2_weighted_de: pd.Series = None,
) -> list[str]:
    if mode == "var":
        s = var_score
    elif mode == "de_ttest":
        s = de_score
    elif mode == "agora_only":
        s = var_score[var_score.index.isin(agora_probe_set)]
    elif mode == "de_agora_intersection":
        s = de_score[de_score.index.isin(agora_probe_set)]
    elif mode == "de_batch_robust":
        s = de_score[de_score.index.isin(robust_features)]
    elif mode == "de_batch_robust_v2":
        if v2_weighted_de is not None:
            s = v2_weighted_de
        else:
            s = de_score
    else:
        raise ValueError(f"Unknown feature_mode: {mode}")

    s = s.reindex(common.intersection(s.index)).dropna()
    if len(s) == 0:
        return []
    return s.sort_values(ascending=False).head(min(top_n, len(s))).index.tolist()


def _neurocombat_from_training_fixed(dat: np.ndarray, batch: np.ndarray, estimates: dict) -> np.ndarray:
    """
    Compatibility wrapper for neuroCombatFromTraining.
    The upstream helper can fail on batch index casting; this mirrors the official
    transform math while using robust index extraction.
    """
    batch = np.array(batch, dtype="str")
    old_levels = np.array(estimates["batches"], dtype="str")

    missing_levels = np.setdiff1d(np.unique(batch), old_levels)
    if missing_levels.shape[0] != 0:
        raise ValueError(f"The batches {missing_levels} are not part of the training dataset")

    wh = [int(np.where(old_levels == x)[0][0]) for x in batch]

    var_pooled = estimates["var.pooled"]
    stand_mean = estimates["stand.mean"][:, 0]
    mod_mean = estimates["mod.mean"]
    gamma_star = estimates["gamma.star"]
    delta_star = estimates["delta.star"]
    n_array = dat.shape[1]

    stand_mean = stand_mean + mod_mean.mean(axis=1)
    stand_mean = np.transpose([stand_mean] * n_array)

    bayesdata = np.subtract(dat, stand_mean) / np.sqrt(var_pooled)
    gamma = np.transpose(gamma_star[wh, :])
    delta = np.transpose(delta_star[wh, :])
    bayesdata = np.subtract(bayesdata, gamma) / np.sqrt(delta)
    bayesdata = bayesdata * np.sqrt(var_pooled) + stand_mean
    return bayesdata


def _combat_trainfit_apply_test(
    X_train: pd.DataFrame,
    batch_train: pd.Series,
    X_test: pd.DataFrame,
    batch_test: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Leakage-safe ComBat using neuroCombat training estimates:
    - fit ComBat estimates on training data only
    - apply frozen estimates to test data
    """
    covars_train = pd.DataFrame({"batch": batch_train.astype(str).values}, index=X_train.index)

    train_fit = neuroCombat(
        dat=X_train.T,
        covars=covars_train,
        batch_col="batch",
    )
    X_train_adj = pd.DataFrame(train_fit["data"].T, index=X_train.index, columns=X_train.columns)

    # Try official helper first; fall back to fixed compatibility wrapper if needed.
    try:
        test_apply = neuroCombatFromTraining(
            dat=X_test.T,
            batch=batch_test.astype(str).values,
            estimates=train_fit["estimates"],
        )
        X_test_data = test_apply["data"]
    except Exception:
        X_test_data = _neurocombat_from_training_fixed(
            dat=X_test.T.values,
            batch=batch_test.astype(str).values,
            estimates=train_fit["estimates"],
        )

    X_test_adj = pd.DataFrame(X_test_data.T, index=X_test.index, columns=X_test.columns)
    return X_train_adj, X_test_adj


def _null_distribution(
    Xtr: pd.DataFrame,
    ytr: pd.Series,
    Xte: pd.DataFrame,
    yte: pd.Series,
    n_perm: int = 1000,
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

    # Map AddNeuroMed datasets to Gene Symbols and collapse duplicates
    for gse in gses:
        X, y = data[gse]
        cols = []
        for c in X.columns:
            sym = probe_to_symbol.get(c, "")
            cols.append(sym if sym else None)
        X.columns = cols
        X = X.loc[:, X.columns.notna()]
        X = X.T.groupby(level=0).mean().T
        data[gse] = (X, y)

    # Load Custom GSE97760 Target
    X_97, y_97, gpl_97 = load_gse97760_custom(raw_dir)
    data["GSE97760"] = (X_97, y_97)
    gses.append("GSE97760")
    gpl_of["GSE97760"] = gpl_97

    agora_df = pd.read_csv(root / "outputs/data/ampad_open_nominated_targets.csv")
    agora_symbols = set(agora_df["hgnc_symbol"].dropna().astype(str).str.upper())

    rows = []
    pred_rows = []
    null_rows = []
    stats: dict[str, dict[str, float | int | str]] = {}

    feature_modes = ["var", "de_ttest", "agora_only", "de_agora_intersection", "de_batch_robust", "de_batch_robust_v2"]

    agora_probe_set_global = set()

    # ── PRIMARY: Bidirectional GSE63060 ↔ GSE63061 ──
    directions = [
        ("GSE63060", "GSE63061"),
        ("GSE63061", "GSE63060"),
    ]

    for source_name, target_name in directions:
        Xs_raw, ys = data[source_name]
        Xt_raw, yt_raw = data[target_name]

        for feature_mode in feature_modes:
            for top_n in [200, 1000]:
                common = _common_genes(Xs_raw, Xt_raw)
                Xs = Xs_raw[common]
                Xt = Xt_raw[common]
                yt = yt_raw

                Xtr_t, Xte_t, ytr, yte = train_test_split(
                    Xt, yt, test_size=0.3, random_state=SEED, stratify=yt
                )

                X_pool_full = pd.concat([Xs, Xtr_t], axis=0)
                y_pool_full = pd.concat([ys, ytr], axis=0)

                # BDP-FS distortion scores (Z-score standardized)
                batch_combo_full = pd.Series(
                    [source_name] * len(Xs) + [target_name] * len(Xtr_t),
                    index=X_pool_full.index,
                )
                covars_full = pd.DataFrame({"batch": batch_combo_full.astype(str).values}, index=X_pool_full.index)
                train_fit_full = neuroCombat(dat=X_pool_full.T, covars=covars_full, batch_col="batch")

                gamma_star = train_fit_full["estimates"]["gamma.star"]
                delta_star = train_fit_full["estimates"]["delta.star"]
                gamma_abs = np.mean(np.abs(gamma_star), axis=0)
                delta_abs = np.mean(np.abs(1.0 - delta_star), axis=0)
                g_s = pd.Series(gamma_abs, index=X_pool_full.columns)
                d_s = pd.Series(delta_abs, index=X_pool_full.columns)
                z_g = (g_s - g_s.mean()) / g_s.std()
                z_d = (d_s - d_s.mean()) / d_s.std()

                # ── v1: Original BDP-FS (equal weight, hard 80th percentile) ──
                distortion_scores = (z_g.abs() + z_d.abs()).fillna(0)
                threshold = distortion_scores.quantile(0.80)
                robust_features = distortion_scores[distortion_scores <= threshold].index

                agora_probe_set = {g for g in common if g.upper() in agora_symbols}

                # ── v2: BDP-FS Adaptive Masterpiece (GMM-Anchored Soft Weighting) ──
                ALPHA = 0.2  # Asymmetric Composition: Penalize variance 4x more than mean
                distortion_v2 = (ALPHA * z_g.abs() + (1 - ALPHA) * z_d.abs()).fillna(0)

                # Agora Shield: discount penalty for nominated targets
                agora_mask = distortion_v2.index.isin(agora_probe_set)
                distortion_v2[agora_mask] = distortion_v2[agora_mask] * 0.5

                # GMM adaptive anchor selection
                v2_vals = distortion_v2.values.reshape(-1, 1)
                gmm = GaussianMixture(n_components=2, random_state=SEED).fit(v2_vals)
                means = gmm.means_.flatten()
                stds = np.sqrt(gmm.covariances_.flatten())
                
                # Anchor tau_0 = 95th percentile of the "Native" (lower mean) component
                idx_native = np.argmin(means)
                tau_0 = float(means[idx_native] + 1.645 * stds[idx_native])
                
                # Soft distortion weighting: Score = |t| * exp(-max(0, D - tau_0))
                DECAY_ALPHA = 1.0 # Gradient of the penalty
                weight_v2 = np.exp(-DECAY_ALPHA * np.maximum(0, distortion_v2 - tau_0))

                _, de_pool_raw = _feature_scores(X_pool_full, y_pool_full, common)
                v2_weighted_de = de_pool_raw * weight_v2.reindex(de_pool_raw.index, fill_value=1.0)

                # Arm-native feature scoring
                var_target, de_target = _feature_scores(Xtr_t, ytr, common)
                var_source, de_source = _feature_scores(Xs, ys, common)
                var_pool, de_pool = _feature_scores(X_pool_full, y_pool_full, common)

                genes_target = _select_top_genes(feature_mode, top_n, common, var_target, de_target, agora_probe_set, robust_features, v2_weighted_de)
                genes_source = _select_top_genes(feature_mode, top_n, common, var_source, de_source, agora_probe_set, robust_features, v2_weighted_de)
                genes_pool = _select_top_genes(feature_mode, top_n, common, var_pool, de_pool, agora_probe_set, robust_features, v2_weighted_de)

                if min(len(genes_target), len(genes_source), len(genes_pool)) < 20:
                    continue

                source_label = source_name
                target_label = target_name

                # target-only
                p_target = _fit_predict(Xtr_t[genes_target], ytr, Xte_t[genes_target])
                m_target = _metrics(yte.values, p_target)
                rows.append({"source": source_label, "target": target_label, "feature_mode": feature_mode, "top_n_genes": top_n, "selected_gene_count": len(genes_target), "arm": "target_only", **m_target})
                for sid, yt_i, pr_i in zip(Xte_t[genes_target].index, yte.values, p_target):
                    pred_rows.append({"source": source_label, "target": target_label, "feature_mode": feature_mode, "top_n_genes": top_n, "arm": "target_only", "sample_id": sid, "y_true": int(yt_i), "y_prob": float(pr_i)})

                # source-only
                p_source = _fit_predict(Xs[genes_source], ys, Xte_t[genes_source])
                m_source = _metrics(yte.values, p_source)
                rows.append({"source": source_label, "target": target_label, "feature_mode": feature_mode, "top_n_genes": top_n, "selected_gene_count": len(genes_source), "arm": "source_only", **m_source})
                for sid, yt_i, pr_i in zip(Xte_t[genes_source].index, yte.values, p_source):
                    pred_rows.append({"source": source_label, "target": target_label, "feature_mode": feature_mode, "top_n_genes": top_n, "arm": "source_only", "sample_id": sid, "y_true": int(yt_i), "y_prob": float(pr_i)})

                # pooled raw
                X_combo_raw = pd.concat([Xs[genes_pool], Xtr_t[genes_pool]], axis=0)
                y_combo_raw = pd.concat([ys, ytr], axis=0)
                p_combo_raw = _fit_predict(X_combo_raw, y_combo_raw, Xte_t[genes_pool])
                m_combo_raw = _metrics(yte.values, p_combo_raw)
                rows.append({"source": source_label, "target": target_label, "feature_mode": feature_mode, "top_n_genes": top_n, "selected_gene_count": len(genes_pool), "arm": "source_plus_target_raw", **m_combo_raw})
                for sid, yt_i, pr_i in zip(Xte_t[genes_pool].index, yte.values, p_combo_raw):
                    pred_rows.append({"source": source_label, "target": target_label, "feature_mode": feature_mode, "top_n_genes": top_n, "arm": "source_plus_target_raw", "sample_id": sid, "y_true": int(yt_i), "y_prob": float(pr_i)})

                # pooled ComBat train-fit/test-apply
                batch_combo = pd.Series([source_label] * len(Xs[genes_pool]) + [target_label] * len(Xtr_t[genes_pool]), index=X_combo_raw.index)
                batch_te = pd.Series([target_label] * len(Xte_t[genes_pool]), index=Xte_t[genes_pool].index)
                X_combo_cb, Xte_cb = _combat_trainfit_apply_test(X_combo_raw, batch_combo, Xte_t[genes_pool], batch_te)
                p_combo_cb = _fit_predict(X_combo_cb, y_combo_raw, Xte_cb)
                m_combo_cb = _metrics(yte.values, p_combo_cb)
                rows.append({"source": source_label, "target": target_label, "feature_mode": feature_mode, "top_n_genes": top_n, "selected_gene_count": len(genes_pool), "arm": "source_plus_target_combat_trainfit", **m_combo_cb})
                for sid, yt_i, pr_i in zip(Xte_cb.index, yte.values, p_combo_cb):
                    pred_rows.append({"source": source_label, "target": target_label, "feature_mode": feature_mode, "top_n_genes": top_n, "arm": "source_plus_target_combat_trainfit", "sample_id": sid, "y_true": int(yt_i), "y_prob": float(pr_i)})

                # null distribution (Submission Rigor: 1,000 permutations)
                null_aucs, p_null_avg = _null_distribution(Xtr_t[genes_target], ytr, Xte_t[genes_target], yte, n_perm=1000)
                null_mean = float(np.mean(null_aucs))
                null_q05 = float(np.quantile(null_aucs, 0.05))
                null_q95 = float(np.quantile(null_aucs, 0.95))

                rows.append({"source": source_label, "target": target_label, "feature_mode": feature_mode, "top_n_genes": top_n, "selected_gene_count": len(genes_target), "arm": "null_label_permutation_mean_auroc", "auroc": null_mean, "auprc": np.nan, "balanced_accuracy": np.nan, "brier": np.nan})
                rows.append({"source": source_label, "target": target_label, "feature_mode": feature_mode, "top_n_genes": top_n, "selected_gene_count": len(genes_target), "arm": "null_label_permutation_avg1000_prob", **_metrics(yte.values, p_null_avg)})
                for sid, yt_i, pr_i in zip(Xte_t[genes_target].index, yte.values, p_null_avg):
                    pred_rows.append({"source": source_label, "target": target_label, "feature_mode": feature_mode, "top_n_genes": top_n, "arm": "null_label_permutation_avg1000_prob", "sample_id": sid, "y_true": int(yt_i), "y_prob": float(pr_i)})

                for i, auc in enumerate(null_aucs):
                    null_rows.append({"source": source_label, "target": target_label, "feature_mode": feature_mode, "top_n_genes": top_n, "perm_index": i, "null_perm_auroc": float(auc)})

                # Agora Shield Proof: count of genes where (unshielded > v1_thresh) but (shielded <= v1_thresh)
                rescued_agora = int(((distortion_scores > threshold) & (distortion_v2 <= threshold) & distortion_v2.index.isin(agora_probe_set)).sum())

                key = f"{feature_mode}__{source_label}_to_{target_label}_top{top_n}"
                p_target_vs_null = float((null_aucs >= m_target["auroc"]).mean())
                stats[key] = {
                    "feature_mode": feature_mode, "source": source_label, "target": target_label,
                    "top_n_genes": int(top_n), "selected_gene_count": int(len(genes_target)),
                    "target_auroc": float(m_target["auroc"]),
                    "delta_auroc_source_plus_target_raw_vs_target_only": float(m_combo_raw["auroc"] - m_target["auroc"]),
                    "delta_auroc_source_plus_target_combat_trainfit_vs_target_only": float(m_combo_cb["auroc"] - m_target["auroc"]),
                    "null_perm_n": int(len(null_aucs)), "null_perm_auroc_mean": null_mean,
                    "null_perm_auroc_std": float(np.std(null_aucs)),
                    "null_perm_auroc_q05": null_q05, "null_perm_auroc_q95": null_q95,
                    "p_target_gt_null_perm": p_target_vs_null,
                    "agora_genes_rescued_by_v2_shield": rescued_agora,
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
        "null_policy": "1000x label permutations (Consensus Submission Rigor); primary null reported as mean permutation AUROC",
        "batch_harmonization": "ComBat neuroCombat train-fit/test-apply arm (frozen training estimates) used as leakage-safe primary arm",
    }
    out_manifest.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    run()
