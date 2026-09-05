#!/usr/bin/env python
"""统计 PRI-General 最终 101,336 个复合物的蛋白频次和序列长度分布。

统计口径：
* RCSB_PDB：一个唯一 pdb_id 视为一个复合物；
* 其他来源：一个唯一 sequence_hash 视为一个复合物；
* 同一复合物中同一蛋白只计一次；PDB 复合物若含多个不同蛋白，则分别计入。
底层审计表可以有多条技术记录，但不会直接作为统计单位。
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import SOURCE_LABELS, SOURCE_ORDER, clean, complex_key, configure_plot_font, sequence_hash


ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "reports" / "final"
OUT = BUILD / "analysis"
AUDIT = ROOT / "data" / "final" / "audit" / "general_stage_audit.csv"

def percentile_summary(values: list[int]) -> dict[str, float | int]:
    if not values:
        return {
            "n": 0, "mean": None, "sd": None, "min": None,
            "p10": None, "p25": None, "p50": None, "p75": None,
            "p90": None, "p95": None, "max": None,
        }
    arr = np.asarray(values, dtype=float)
    return {
        "n": int(arr.size),
        "mean": float(arr.mean()),
        "sd": float(arr.std(ddof=1)) if arr.size > 1 else 0.0,
        "min": int(arr.min()),
        "p10": float(np.percentile(arr, 10)),
        "p25": float(np.percentile(arr, 25)),
        "p50": float(np.percentile(arr, 50)),
        "p75": float(np.percentile(arr, 75)),
        "p90": float(np.percentile(arr, 90)),
        "p95": float(np.percentile(arr, 95)),
        "max": int(arr.max()),
    }


def length_bins(values: list[int], bins: list[int], labels: list[str]) -> pd.DataFrame:
    counts = pd.cut(pd.Series(values), bins=bins, labels=labels, include_lowest=True, right=True).value_counts()
    return pd.DataFrame({"length_bin": labels, "count": [int(counts.get(label, 0)) for label in labels]})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    usecols = [
        "source_database", "pdb_id", "sequence_hash", "protein_sequence",
        "rna_sequence", "protein_name", "protein_accession", "stage_queue",
    ]

    complexes: dict[str, dict[str, object]] = {}
    protein_sequences: dict[str, str] = {}
    rna_sequences: dict[str, str] = {}
    protein_names: defaultdict[str, Counter[str]] = defaultdict(Counter)
    protein_accessions: defaultdict[str, Counter[str]] = defaultdict(Counter)

    for chunk in pd.read_csv(AUDIT, usecols=usecols, chunksize=50_000, low_memory=False):
        # pandas 可能按原始 CSV 的列顺序返回 usecols；这里显式重排，
        # 防止 itertuples 将蛋白、RNA、来源等字段错位。
        chunk = chunk.loc[:, usecols]
        queue = chunk["stage_queue"].fillna(False).astype(bool)
        for row_index, row in zip(chunk.index, chunk.itertuples(index=False)):
            source_raw, pdb_id_raw, pair_hash_raw, protein_raw, rna_raw, name_raw, accession_raw, _ = row
            source_raw = clean(source_raw)
            source = SOURCE_LABELS.get(source_raw, source_raw or "未知来源")
            pdb_id = clean(pdb_id_raw)
            pair_hash = clean(pair_hash_raw)
            protein = clean(protein_raw)
            rna = clean(rna_raw)
            if not protein or not rna:
                continue
            key = complex_key(source_raw, pdb_id, pair_hash)
            record = complexes.setdefault(key, {
                "source": source,
                "source_raw": source_raw,
                "proteins": set(),
                "rnas": set(),
                "all_proteins": set(),
                "all_rnas": set(),
            })
            protein_hash = sequence_hash(protein)
            rna_hash = sequence_hash(rna)
            record["all_proteins"].add(protein_hash)
            record["all_rnas"].add(rna_hash)
            if bool(queue.loc[row_index]):
                record["proteins"].add(protein_hash)
                record["rnas"].add(rna_hash)
            protein_sequences.setdefault(protein_hash, protein)
            rna_sequences.setdefault(rna_hash, rna)
            name = clean(name_raw)
            accession = clean(accession_raw)
            if name:
                protein_names[protein_hash][name] += 1
            if accession:
                protein_accessions[protein_hash][accession] += 1

    source_protein_complex_counts: Counter[str] = Counter()
    source_protein_source_counts: defaultdict[str, Counter[str]] = defaultdict(Counter)
    final_complex_keys = {key for key, record in complexes.items() if record["proteins"]}
    final_protein_counts: Counter[str] = Counter()
    final_protein_source_counts: defaultdict[str, Counter[str]] = defaultdict(Counter)
    for key, record in complexes.items():
        source = str(record["source"])
        for protein_hash in record["all_proteins"]:
            source_protein_complex_counts[protein_hash] += 1
            source_protein_source_counts[protein_hash][source] += 1
        if key in final_complex_keys:
            for protein_hash in record["proteins"]:
                final_protein_counts[protein_hash] += 1
                final_protein_source_counts[protein_hash][source] += 1

    protein_complex_counts = final_protein_counts
    protein_source_counts: defaultdict[str, Counter[str]] = defaultdict(Counter)
    protein_source_counts.update(final_protein_source_counts)
    rna_complex_counts: Counter[str] = Counter()
    rna_source_counts: defaultdict[str, Counter[str]] = defaultdict(Counter)
    protein_length_rows: list[dict[str, object]] = []
    rna_length_rows: list[dict[str, object]] = []
    source_complexes: defaultdict[str, set[str]] = defaultdict(set)

    for key in final_complex_keys:
        record = complexes[key]
        source = str(record["source"])
        source_complexes[source].add(key)
        for protein_hash in record["proteins"]:
            protein_length_rows.append({
                "complex_id": key,
                "source": source,
                "protein_hash": protein_hash,
                "length_aa": len(protein_sequences[protein_hash]),
            })
        for rna_hash in record["rnas"]:
            rna_complex_counts[rna_hash] += 1
            rna_source_counts[rna_hash][source] += 1
            rna_length_rows.append({
                "complex_id": key,
                "source": source,
                "rna_hash": rna_hash,
                "length_nt": len(rna_sequences[rna_hash]),
            })

    source_columns = SOURCE_ORDER
    protein_rows = []
    for protein_hash, count in protein_complex_counts.items():
        name = protein_names[protein_hash].most_common(1)[0][0] if protein_names[protein_hash] else "未标注"
        accession = protein_accessions[protein_hash].most_common(1)[0][0] if protein_accessions[protein_hash] else ""
        row = {
            "protein_hash": protein_hash,
            "protein_name": name,
            "protein_accession": accession,
            "protein_length_aa": len(protein_sequences[protein_hash]),
            "complex_count": int(count),
            "source_count": int(sum(1 for source in source_columns if protein_source_counts[protein_hash][source] > 0)),
        }
        for source in source_columns:
            row[source] = int(protein_source_counts[protein_hash][source])
        protein_rows.append(row)
    protein_frequency = pd.DataFrame(protein_rows).sort_values(
        ["complex_count", "protein_name"], ascending=[False, True]
    ).reset_index(drop=True)
    protein_frequency.insert(0, "rank", np.arange(1, len(protein_frequency) + 1))
    protein_frequency.to_csv(OUT / "protein_complex_frequency.csv", index=False, encoding="utf-8-sig")

    top20 = protein_frequency.head(20).copy()
    top20.to_csv(OUT / "protein_top20_source_breakdown.csv", index=False, encoding="utf-8-sig")

    source_rank = sorted(source_protein_complex_counts, key=lambda item: (-source_protein_complex_counts[item], item))
    source_rank_index = {protein_hash: rank for rank, protein_hash in enumerate(source_rank, start=1)}
    final_rank = sorted(final_protein_counts, key=lambda item: (-final_protein_counts[item], item))
    final_rank_index = {protein_hash: rank for rank, protein_hash in enumerate(final_rank, start=1)}
    before_top20_rows = []
    for rank, protein_hash in enumerate(final_rank[:20], start=1):
        name = protein_names[protein_hash].most_common(1)[0][0] if protein_names[protein_hash] else "未标注"
        accession = protein_accessions[protein_hash].most_common(1)[0][0] if protein_accessions[protein_hash] else ""
        before = int(source_protein_complex_counts[protein_hash])
        after = int(final_protein_counts[protein_hash])
        row = {
            "rank_before": source_rank_index.get(protein_hash, ""),
            "rank_after": rank,
            "protein_hash": protein_hash,
            "protein_name": name,
            "protein_accession": accession,
            "protein_length_aa": len(protein_sequences[protein_hash]),
            "source_complex_count": before,
            "final_complex_count": after,
            "retained_fraction": after / before if before else 0.0,
            "removed_complex_count": before - after,
        }
        for source in SOURCE_ORDER:
            row[f"source_{source}"] = int(source_protein_source_counts[protein_hash][source])
            row[f"final_{source}"] = int(final_protein_source_counts[protein_hash][source])
        before_top20_rows.append(row)
    top20_before_after = pd.DataFrame(before_top20_rows)
    top20_before_after.to_csv(OUT / "protein_top20_source_before_after.csv", index=False, encoding="utf-8-sig")

    configure_plot_font()
    if not top20_before_after.empty:
        fig, ax = plt.subplots(figsize=(11, 8))
        plot_data = top20_before_after.iloc[::-1].reset_index(drop=True)
        y = np.arange(len(plot_data))
        height = 0.36
        ax.barh(y - height / 2, plot_data["source_complex_count"], height=height,
                color="#b7c9dc", label="源数据（筛选前）")
        ax.barh(y + height / 2, plot_data["final_complex_count"], height=height,
                color="#e67e22", label="最终队列（筛选后）")
        ax.set_yticks(y)
        ax.set_yticklabels(plot_data["protein_name"].astype(str))
        ax.set_xlabel("包含该蛋白的复合物数")
        ax.set_title("Top 20 蛋白在源数据和最终队列中的变化")
        ax.grid(axis="x", color="#dddddd", linewidth=0.7)
        ax.legend(loc="lower right")
        fig.tight_layout()
        fig.savefig(OUT / "protein_top20_source_before_after.png", dpi=250, bbox_inches="tight", facecolor="white")
        plt.close(fig)

    protein_length = pd.DataFrame(protein_length_rows)
    rna_length = pd.DataFrame(rna_length_rows)
    length_rows = []
    for source in SOURCE_ORDER:
        pvals = protein_length.loc[protein_length["source"].eq(source), "length_aa"].astype(int).tolist()
        rvals = rna_length.loc[rna_length["source"].eq(source), "length_nt"].astype(int).tolist()
        for molecule, values in [("protein", pvals), ("RNA", rvals)]:
            summary = percentile_summary(values)
            summary.update({"source": source, "molecule": molecule})
            length_rows.append(summary)
    length_summary = pd.DataFrame(length_rows)[
        ["source", "molecule", "n", "mean", "sd", "min", "p10", "p25", "p50", "p75", "p90", "p95", "max"]
    ]
    length_summary.to_csv(OUT / "complex_component_length_summary_by_source.csv", index=False, encoding="utf-8-sig")
    protein_length.to_csv(OUT / "final_complex_protein_lengths.csv", index=False, encoding="utf-8-sig")
    rna_length.to_csv(OUT / "final_complex_rna_lengths.csv", index=False, encoding="utf-8-sig")

    p_bins = [39, 199, 399, 599, 799, 999, 1299, 1599, 2000]
    p_labels = ["40–199", "200–399", "400–599", "600–799", "800–999", "1000–1299", "1300–1599", "1600–2000"]
    r_bins = [9, 19, 49, 99, 199, 299, 399, 500]
    r_labels = ["10–19", "20–49", "50–99", "100–199", "200–299", "300–399", "400–500"]
    bin_rows = []
    for source in SOURCE_ORDER:
        pvals = protein_length.loc[protein_length["source"].eq(source), "length_aa"].astype(int).tolist()
        rvals = rna_length.loc[rna_length["source"].eq(source), "length_nt"].astype(int).tolist()
        for molecule, values, bins, labels in [("protein", pvals, p_bins, p_labels), ("RNA", rvals, r_bins, r_labels)]:
            table = length_bins(values, bins, labels)
            table.insert(0, "molecule", molecule)
            table.insert(0, "source", source)
            bin_rows.append(table)
    pd.concat(bin_rows, ignore_index=True).to_csv(OUT / "complex_component_length_bins_by_source.csv", index=False, encoding="utf-8-sig")

    counts = {source: len(source_complexes[source]) for source in SOURCE_ORDER}
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(AUDIT),
        "statistical_unit": {
            "RCSB_PDB": "unique pdb_id",
            "other_sources": "unique sequence_hash",
            "protein_frequency": "number of final complexes containing the protein once or more",
            "length_distribution": "unique protein/RNA component sequences within final complexes",
        },
        "final_complex_count": int(len(final_complex_keys)),
        "source_resolved_complex_count": int(len(complexes)),
        "complex_count_by_source": counts,
        "unique_proteins_in_complexes": int(len(protein_complex_counts)),
        "unique_rnas_in_complexes": int(len(rna_complex_counts)),
        "protein_top1_complex_fraction": float(protein_frequency.iloc[0]["complex_count"] / len(final_complex_keys)) if len(protein_frequency) else 0.0,
        "protein_top10_complex_fraction": float(protein_frequency.head(10)["complex_count"].sum() / len(final_complex_keys)) if len(protein_frequency) else 0.0,
    }
    (OUT / "final_complex_analysis_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    if not top20.empty:
        fig, ax = plt.subplots(figsize=(10, 7))
        labels = top20["protein_name"].astype(str).str.slice(0, 28)
        ax.barh(labels.iloc[::-1], top20["complex_count"].iloc[::-1], color="#7f8fa6")
        ax.set_xlabel("出现该蛋白的复合物数")
        ax.set_title("最终 PRI-General 中出现频次最高的 20 个蛋白")
        ax.grid(axis="x", color="#dddddd", linewidth=0.7)
        fig.tight_layout()
        fig.savefig(OUT / "top20_protein_complex_frequency.png", dpi=250, bbox_inches="tight", facecolor="white")
        plt.close(fig)

        plot_source = top20.set_index("protein_name")[source_columns]
        fig, ax = plt.subplots(figsize=(11, 7))
        plot_source.iloc[::-1].plot.barh(stacked=True, ax=ax, color=["#8da0cb", "#fc8d62", "#66c2a5", "#b3b3b3"])
        ax.set_xlabel("复合物数")
        ax.set_ylabel("")
        ax.set_title("Top 20 蛋白的来源组成")
        ax.grid(axis="x", color="#dddddd", linewidth=0.7)
        ax.legend(title="来源", loc="lower right", fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT / "top20_protein_source_breakdown.png", dpi=250, bbox_inches="tight", facecolor="white")
        plt.close(fig)

    colors = ["#4c78a8", "#f58518", "#54a24b", "#999999"]
    fig, axes = plt.subplots(2, 1, figsize=(11, 10))
    for ax, frame, column, xlabel, title, bins in [
        (axes[0], protein_length, "length_aa", "蛋白长度（aa）", "各来源蛋白长度分布", np.linspace(40, 2000, 50)),
        (axes[1], rna_length, "length_nt", "RNA长度（nt）", "各来源 RNA 长度分布", np.linspace(10, 500, 50)),
    ]:
        for source, color in zip(SOURCE_ORDER, colors):
            values = frame.loc[frame["source"].eq(source), column].astype(int)
            if len(values):
                ax.hist(values, bins=bins, density=True, alpha=0.38, label=source, color=color)
                ax.axvline(values.median(), color=color, linewidth=1.2, linestyle="--")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("密度")
        ax.set_title(title)
        ax.grid(axis="y", color="#dddddd", linewidth=0.7)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "protein_rna_length_distributions_by_source.png", dpi=250, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(OUT)


if __name__ == "__main__":
    main()
