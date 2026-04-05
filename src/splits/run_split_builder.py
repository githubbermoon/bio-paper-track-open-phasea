from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from src.splits.leakage_safe_split import build_strict_forward_site_grouped_splits


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="Input CSV path")
    p.add_argument("--out-manifest", required=True)
    p.add_argument("--out-report", required=True)
    p.add_argument("--time-col", default="collection_time")
    p.add_argument("--subject-col", default="subject_id")
    args = p.parse_args()

    df = pd.read_csv(args.input)
    if args.time_col in df.columns:
        df[args.time_col] = pd.to_datetime(df[args.time_col], errors="coerce")

    splits = build_strict_forward_site_grouped_splits(
        df,
        time_col=args.time_col,
        subject_col=args.subject_col,
        train_frac=0.7,
        valid_frac=0.15,
    )

    manifest = {
        "input": str(Path(args.input).resolve()),
        "time_col": args.time_col,
        "subject_col": args.subject_col,
        "audit": splits["audit"],
    }

    Path(args.out_manifest).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_report).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_manifest).write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    report = [
        "# Leakage Audit Report v1",
        "",
        f"- rows: {splits['audit']['n_rows']}",
        f"- train: {splits['audit']['n_train']}",
        f"- valid: {splits['audit']['n_valid']}",
        f"- test: {splits['audit']['n_test']}",
        f"- subject_overlap: {splits['audit']['subject_overlap']}",
    ]
    Path(args.out_report).write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
