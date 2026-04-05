from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from urllib.request import urlopen

import pandas as pd

API_URL = "https://agora.adknowledgeportal.org/api/v1/genes/nominated"


def _fetch_json(url: str) -> dict:
    with urlopen(url, timeout=60) as resp:
        return json.loads(resp.read().decode("utf-8"))


def run() -> None:
    root = Path(__file__).resolve().parents[2]
    out_data = root / "data/raw/open_ampad/agora_genes_nominated.json"
    out_csv = root / "outputs/data/ampad_open_nominated_targets.csv"
    out_summary = root / "outputs/tables/ampad_open_subset_summary.csv"
    out_md = root / "outputs/ampad_open_subset_notes.md"

    payload = _fetch_json(API_URL)
    items = payload.get("items", [])

    out_data.parent.mkdir(parents=True, exist_ok=True)
    out_data.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    rows = []
    modality_counter = Counter()
    source_counter = Counter()
    study_counter = Counter()

    for it in items:
        ensembl = it.get("ensembl_gene_id")
        hgnc = it.get("hgnc_symbol")
        for nom in it.get("target_nominations", []) or []:
            src = nom.get("source")
            study = nom.get("study")
            input_data = nom.get("input_data")
            if src:
                source_counter[src] += 1
            if study:
                study_counter[study] += 1
            if input_data:
                for tok in [x.strip() for x in str(input_data).split(",") if x.strip()]:
                    modality_counter[tok] += 1
            rows.append(
                {
                    "ensembl_gene_id": ensembl,
                    "hgnc_symbol": hgnc,
                    "source": src,
                    "team": nom.get("team"),
                    "study": study,
                    "input_data": input_data,
                    "details": nom.get("details"),
                }
            )

    df = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)

    summary_rows = []
    summary_rows.append({"metric": "unique_nominated_genes", "value": int(pd.Series([r['ensembl_gene_id'] for r in rows]).nunique())})
    summary_rows.append({"metric": "total_nominations", "value": int(len(rows))})
    for k, v in source_counter.most_common(10):
        summary_rows.append({"metric": f"source::{k}", "value": int(v)})
    for k, v in modality_counter.most_common(15):
        summary_rows.append({"metric": f"input_data::{k}", "value": int(v)})
    for k, v in study_counter.most_common(10):
        summary_rows.append({"metric": f"study::{k}", "value": int(v)})

    pd.DataFrame(summary_rows).to_csv(out_summary, index=False)

    md = [
        "# AMP-AD Open Subset Notes",
        "",
        f"API source: {API_URL}",
        f"Nominated genes retrieved: {len(items)}",
        f"Expanded nomination rows: {len(rows)}",
        "",
        "Top sources:",
    ]
    for k, v in source_counter.most_common(8):
        md.append(f"- {k}: {v}")
    md.append("")
    md.append("Top input_data modalities:")
    for k, v in modality_counter.most_common(10):
        md.append(f"- {k}: {v}")

    out_md.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    run()
