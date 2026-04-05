from __future__ import annotations

import math
from typing import Any

import pandas as pd


def build_strict_forward_site_grouped_splits(
    df: pd.DataFrame,
    time_col: str,
    subject_col: str,
    train_frac: float = 0.7,
    valid_frac: float = 0.15,
) -> dict[str, list[int] | dict[str, Any]]:
    if time_col not in df.columns:
        raise ValueError(f"Missing time column: {time_col}")
    if subject_col not in df.columns:
        raise ValueError(f"Missing subject column: {subject_col}")
    if not (0 < train_frac < 1) or not (0 < valid_frac < 1):
        raise ValueError("train_frac and valid_frac must be in (0,1)")
    if train_frac + valid_frac >= 1:
        raise ValueError("train_frac + valid_frac must be < 1")

    ordered = df.sort_values(time_col).reset_index().rename(columns={"index": "_orig_idx"})
    n = len(ordered)
    if n < 3:
        raise ValueError("Need at least 3 rows")

    train_end = max(1, math.floor(n * train_frac))
    valid_end = max(train_end + 1, math.floor(n * (train_frac + valid_frac)))
    valid_end = min(valid_end, n - 1)

    train_idx = ordered.iloc[:train_end]["_orig_idx"].astype(int).tolist()
    valid_idx = ordered.iloc[train_end:valid_end]["_orig_idx"].astype(int).tolist()
    test_idx = ordered.iloc[valid_end:]["_orig_idx"].astype(int).tolist()

    train_subjects = set(df.loc[train_idx, subject_col].astype(str))
    valid_subjects = set(df.loc[valid_idx, subject_col].astype(str))
    test_subjects = set(df.loc[test_idx, subject_col].astype(str))

    if train_subjects & valid_subjects or train_subjects & test_subjects or valid_subjects & test_subjects:
        raise ValueError("Subject overlap across splits detected")

    return {
        "train_idx": train_idx,
        "valid_idx": valid_idx,
        "test_idx": test_idx,
        "audit": {
            "n_rows": n,
            "n_train": len(train_idx),
            "n_valid": len(valid_idx),
            "n_test": len(test_idx),
            "subject_overlap": False,
        },
    }
