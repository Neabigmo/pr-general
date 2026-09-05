"""Generate additional structure-level PRI-General analysis figures.

This script works from the current audit table and the current complex-level
waterfall.  It treats one non-PDB ``sequence_hash`` or one PDB ``pdb_id`` as
one complex, and never uses raw technical rows as the main counting unit.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (
    SOURCE_COLORS,
    SOURCE_LABELS,
    SOURCE_ORDER,
    clean,
    complex_key,
    configure_plot_font,
    sequence_hash,
)


ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "reports" / "final"
AUDIT = ROOT / "data" / "final" / "audit" / "general_stage_audit.csv"
WATERFALL = ROOT / "deliverables" / "PRI-General_complex_stage_waterfall_latest.csv"
OUT = BUILD / "analysis" / "extended_analysis"
DELIVERABLES = ROOT / "deliverables" / "analysis"

def flag(value: object) -> bool:
    return clean(value).lower() in {"true", "1", "yes", "y"}


def add_nonempty(counter: Counter[str], value: object) -> None:
    value = clean(value) or "未标注"
    counter[value] += 1


def load_complexes() -> tuple[dict[str, dict], Counter[str], Counter[str]]:
    usecols = [
        "source_database", "pdb_id", "sequence_hash", "protein_sequence", "rna_sequence",
        "protein_name", "organism", "taxonomy_id", "system_class", "system_subtype",
        "rna_class", "interaction_evidence", "experimental_structure_available",
        "stage_queue", "stage_d1", "stage_d2", "stage_d3", "stage_protein_capped",
        "stage_family_capped", "stage_exclusion_reason",
    ]
    complexes: dict[str, dict] = {}
    protein_names: defaultdict[str, Counter[str]] = defaultdict(Counter)
    rna_names: defaultdict[str, Counter[str]] = defaultdict(Counter)
    for chunk in pd.read_csv(AUDIT, usecols=usecols, chunksize=50_000, low_memory=False):
        chunk = chunk.loc[:, usecols]
        for row in chunk.itertuples(index=False, name=None):
            (
                source_raw, pdb_id_raw, pair_hash_raw, protein_raw, rna_raw,
                protein_name_raw, organism_raw, taxid_raw, system_raw, subtype_raw,
                rna_class_raw, evidence_raw, structure_raw, queue_raw,
                d1_raw, d2_raw, d3_raw, capped_raw, family_raw, exclusion_raw,
            ) = row
            source_raw = clean(source_raw)
            pdb_id = clean(pdb_id_raw)
            pair_hash = clean(pair_hash_raw)
            protein = clean(protein_raw)
            rna = clean(rna_raw)
            key = complex_key(source_raw, pdb_id, pair_hash)
            if not key or not protein or not rna:
                continue
            source = SOURCE_LABELS.get(source_raw, source_raw or "未知来源")
            rec = complexes.setdefault(
                key,
                {
                    "source": source,
                    "source_raw": source_raw,
                    "pdb_id": pdb_id,
                    "protein_hashes": set(),
                    "rna_hashes": set(),
                    "queue_protein_hashes": set(),
                    "queue_rna_hashes": set(),
                    "protein_lengths": {},
                    "rna_lengths": {},
                    "meta": {},
                },
            )
            protein_hash = sequence_hash(protein)
            rna_hash = sequence_hash(rna)
            rec["protein_hashes"].add(protein_hash)
            rec["rna_hashes"].add(rna_hash)
            rec["protein_lengths"].setdefault(protein_hash, len(protein))
            rec["rna_lengths"].setdefault(rna_hash, len(rna))
            if flag(queue_raw):
                rec["queue_protein_hashes"].add(protein_hash)
                rec["queue_rna_hashes"].add(rna_hash)
            for name, value in {
                "organism": organism_raw,
                "taxonomy_id": taxid_raw,
                "system_class": system_raw,
                "system_subtype": subtype_raw,
                "rna_class": rna_class_raw,
                "interaction_evidence": evidence_raw,
                "experimental_structure_available": structure_raw,
            }.items():
                if name not in rec["meta"] and clean(value):
                    rec["meta"][name] = clean(value)
            if clean(protein_name_raw):
                protein_names[protein_hash][clean(protein_name_raw)] += 1
            d1 = flag(d1_raw)
            d2 = flag(d2_raw)
            d3 = flag(d3_raw)
            capped = flag(capped_raw)
            family = flag(family_raw)
            if not d1:
                stage = "D1 序列/长度检查"
            elif not d2:
                stage = "D2 exact pair 去重"
            elif not d3:
                stage = "D3 RNA90 去冗余"
            elif not capped:
                stage = "单蛋白限额"
            elif not family:
                stage = "Protein30 家族限额"
            else:
                stage = "最终队列"
            rec.setdefault("stage_candidates", set()).add(stage)

    final_keys = {key for key, rec in complexes.items() if rec["queue_protein_hashes"]}
    final_protein_counts: Counter[str] = Counter()
    final_rna_counts: Counter[str] = Counter()
    for key in final_keys:
        rec = complexes[key]
        for value in rec["queue_protein_hashes"]:
            final_protein_counts[value] += 1
        for value in rec["queue_rna_hashes"]:
            final_rna_counts[value] += 1
    return complexes, final_protein_counts, final_rna_counts


def write_distribution(rows: list[dict], path: Path) -> None:
    pd.DataFrame(rows).to_csv(path, index=False, encoding="utf-8-sig")


def make_tables(complexes: dict[str, dict], protein_counts: Counter[str], rna_counts: Counter[str]) -> dict:
    final = {key: rec for key, rec in complexes.items() if rec["queue_protein_hashes"]}
    source_counts = Counter(rec["source"] for rec in final.values())
    total = sum(source_counts.values())
    source_rows = [
        {"source": source, "complex_count": source_counts[source], "fraction": source_counts[source] / total}
        for source in SOURCE_ORDER
    ]
    write_distribution(source_rows, OUT / "source_distribution.csv")

    categorical = {
        "system_class": "system_class_distribution.csv",
        "system_subtype": "system_subtype_distribution.csv",
        "rna_class": "rna_class_distribution.csv",
        "organism": "organism_distribution.csv",
        "interaction_evidence": "evidence_distribution.csv",
    }
    for field, filename in categorical.items():
        rows = []
        for source in SOURCE_ORDER:
            counts = Counter()
            for rec in final.values():
                if rec["source"] == source:
                    add_nonempty(counts, rec["meta"].get(field, ""))
            source_total = sum(counts.values())
            for value, count in counts.most_common():
                rows.append({"source": source, "category": value, "complex_count": count, "fraction_within_source": count / source_total if source_total else 0})
        write_distribution(rows, OUT / filename)

    issue_fields = {
        "organism": "organism",
        "taxonomy_id": "taxonomy_id",
        "system_class": "system_class",
        "rna_class": "rna_class",
        "interaction_evidence": "interaction_evidence",
        "experimental_structure_available": "experimental_structure_available",
    }
    quality_rows = []
    for field, label in issue_fields.items():
        missing = sum(1 for rec in final.values() if not clean(rec["meta"].get(field, "")))
        quality_rows.append({"field": label, "scope": "final complexes", "missing_complex_count": missing, "missing_fraction": missing / total if total else 0})
    quality_rows.extend([
        {"field": "protein_sequence", "scope": "final complexes", "missing_complex_count": 0, "missing_fraction": 0.0},
        {"field": "rna_sequence", "scope": "final complexes", "missing_complex_count": 0, "missing_fraction": 0.0},
    ])
    write_distribution(quality_rows, OUT / "quality_issues.csv")

    for molecule, counts, filename in [
        ("protein", protein_counts, "protein_complex_frequency.csv"),
        ("RNA", rna_counts, "rna_complex_frequency.csv"),
    ]:
        values = sorted(counts.items(), key=lambda x: (-x[1], x[0]))
        frame = pd.DataFrame([{"rank": i, "sequence_hash": key, "complex_count": count} for i, (key, count) in enumerate(values, 1)])
        frame.to_csv(OUT / filename, index=False, encoding="utf-8-sig")

    bins = [(1, 1, "1"), (2, 5, "2–5"), (6, 10, "6–10"), (11, 50, "11–50"), (51, 100, "51–100"), (101, 500, "101–500"), (501, 1000, "501–1000"), (1001, 10**9, ">1000")]
    for molecule, counts, filename in [("protein", protein_counts, "protein_frequency_bins.csv"), ("RNA", rna_counts, "rna_frequency_bins.csv")]:
        rows = []
        for low, high, label in bins:
            n_components = sum(low <= n <= high for n in counts.values())
            rows.append({"molecule": molecule, "frequency_bin": label, "component_count": n_components, "component_fraction": n_components / len(counts) if counts else 0, "complex_occurrences": sum(n for n in counts.values() if low <= n <= high)})
        pd.DataFrame(rows).to_csv(OUT / filename, index=False, encoding="utf-8-sig")

    waterfall = pd.read_csv(WATERFALL)
    waterfall["retention_from_input"] = waterfall["complexes"] / float(waterfall.iloc[0]["complexes"])
    waterfall["loss_from_previous"] = waterfall["complexes"].shift(1).fillna(waterfall.iloc[0]["complexes"]) - waterfall["complexes"]
    waterfall.to_csv(OUT / "stage_waterfall.csv", index=False, encoding="utf-8-sig")
    return {"final": final, "source_counts": source_counts, "total": total, "waterfall": waterfall}


def top_categories(path: Path, category_limit: int = 7) -> tuple[list[str], pd.DataFrame]:
    frame = pd.read_csv(path)
    totals = frame.groupby("category")["complex_count"].sum().sort_values(ascending=False)
    keep = list(totals.head(category_limit).index)
    frame = frame.copy()
    frame.loc[~frame["category"].isin(keep), "category"] = "其他"
    frame = frame.groupby(["source", "category"], as_index=False)["complex_count"].sum()
    return keep + (["其他"] if len(totals) > category_limit else []), frame


def plot_waterfall(ax, waterfall: pd.DataFrame) -> None:
    x = np.arange(len(waterfall))
    bars = ax.bar(x, waterfall["complexes"], color="#6B8FB3")
    ax.set_xticks(x, waterfall["stage"])
    ax.tick_params(axis="x", rotation=28, labelsize=8)
    ax.set_ylabel("复合物数")
    ax.set_title("筛选流程中的复合物数量变化")
    ax.grid(axis="y", alpha=0.22)
    for bar, value in zip(bars, waterfall["complexes"]):
        ax.text(bar.get_x() + bar.get_width() / 2, value, f"{int(value):,}", ha="center", va="bottom", fontsize=8)


def plot_frequency_bins(ax, filename: str, title: str) -> None:
    frame = pd.read_csv(OUT / filename)
    ax.bar(frame["frequency_bin"], frame["component_count"], color="#7A9E9F")
    ax.set_title(title)
    ax.set_ylabel("不同序列组件数")
    ax.set_xlabel("出现在多少个复合物中")
    ax.tick_params(axis="x", rotation=30)
    ax.grid(axis="y", alpha=0.22)
    for i, value in enumerate(frame["component_count"]):
        ax.text(i, value, f"{int(value):,}", ha="center", va="bottom", fontsize=8)


def plot_category_stack(ax, filename: str, title: str, ylabel: str, normalize: bool = True) -> None:
    categories, frame = top_categories(OUT / filename)
    pivot = frame.pivot(index="source", columns="category", values="complex_count").reindex(SOURCE_ORDER).fillna(0)
    if normalize:
        pivot = pivot.div(pivot.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    palette = ["#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2", "#B279A2", "#FF9DA6", "#BAB0AC"]
    left = np.zeros(len(pivot))
    for color, category in zip(palette, categories):
        values = pivot[category].to_numpy() if category in pivot else np.zeros(len(pivot))
        ax.barh(pivot.index, values, left=left, color=color, label=category)
        left += values
    ax.invert_yaxis()
    ax.set_title(title)
    ax.set_xlabel("本来源占比" if normalize else "复合物数")
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False, fontsize=7, ncol=2, loc="lower right")
    ax.grid(axis="x", alpha=0.22)
    if normalize:
        ax.set_xlim(0, 1)
        ax.set_xticks(np.linspace(0, 1, 6), [f"{x:.0%}" for x in np.linspace(0, 1, 6)])


def plot_organism(ax) -> None:
    frame = pd.read_csv(OUT / "organism_distribution.csv")
    totals = frame.groupby("category")["complex_count"].sum().sort_values(ascending=False).head(12)
    frame = frame[frame["category"].isin(totals.index)].copy()
    pivot = frame.pivot(index="category", columns="source", values="complex_count").fillna(0)
    pivot = pivot.reindex(totals.index)
    pivot.plot.barh(ax=ax, stacked=True, color=[SOURCE_COLORS[s] for s in SOURCE_ORDER])
    ax.invert_yaxis()
    ax.set_title("最终复合物的主要物种来源")
    ax.set_xlabel("复合物数")
    ax.set_ylabel("")
    ax.legend(frameon=False, fontsize=7, ncol=2, loc="lower right")
    ax.grid(axis="x", alpha=0.22)


def plot_evidence(ax) -> None:
    frame = pd.read_csv(OUT / "evidence_distribution.csv")
    totals = frame.groupby("category")["complex_count"].sum().sort_values(ascending=False)
    keep = list(totals.head(6).index)
    frame.loc[~frame["category"].isin(keep), "category"] = "其他"
    pivot = frame.groupby(["source", "category"], as_index=False)["complex_count"].sum().pivot(index="source", columns="category", values="complex_count").reindex(SOURCE_ORDER).fillna(0)
    pivot = pivot.div(pivot.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    pivot.plot.barh(ax=ax, stacked=True, color=["#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2", "#B279A2", "#BAB0AC"])
    ax.invert_yaxis()
    ax.set_title("interaction evidence 构成")
    ax.set_xlabel("本来源占比")
    ax.set_ylabel("")
    ax.legend(frameon=False, fontsize=7, ncol=2, loc="lower right")
    ax.grid(axis="x", alpha=0.22)
    ax.set_xlim(0, 1)
    ax.set_xticks(np.linspace(0, 1, 6), [f"{x:.0%}" for x in np.linspace(0, 1, 6)])


def plot_quality(ax) -> None:
    frame = pd.read_csv(OUT / "quality_issues.csv")
    frame = frame.sort_values("missing_complex_count", ascending=True)
    colors = ["#D55E00" if x else "#7A9E9F" for x in frame["missing_complex_count"]]
    ax.barh(frame["field"], frame["missing_complex_count"], color=colors)
    ax.set_title("最终复合物字段缺失情况")
    ax.set_xlabel("缺失的复合物数")
    ax.grid(axis="x", alpha=0.22)
    for i, value in enumerate(frame["missing_complex_count"]):
        ax.text(value + 20, i, f"{int(value):,}", va="center", fontsize=8)


def plot_length_bins(ax, molecule: str, title: str) -> None:
    source_path = BUILD / "analysis" / "complex_component_length_bins_by_source.csv"
    frame = pd.read_csv(source_path)
    frame = frame[frame["molecule"].eq(molecule)]
    pivot = frame.pivot(index="source", columns="length_bin", values="count").reindex(SOURCE_ORDER).fillna(0)
    pivot = pivot.div(pivot.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_yticks(np.arange(len(pivot.index)), pivot.index)
    ax.set_xticks(np.arange(len(pivot.columns)), pivot.columns, rotation=35, ha="right")
    ax.set_title(title)
    ax.set_xlabel("长度区间")
    ax.set_ylabel("来源")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            value = pivot.iloc[i, j]
            if value > 0.005:
                ax.text(j, i, f"{value:.0%}", ha="center", va="center", fontsize=7, color="white" if value > 0.48 else "#222222")
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("本来源占比")


def main() -> None:
    configure_plot_font()
    OUT.mkdir(parents=True, exist_ok=True)
    DELIVERABLES.mkdir(parents=True, exist_ok=True)
    complexes, protein_counts, rna_counts = load_complexes()
    stats = make_tables(complexes, protein_counts, rna_counts)

    fig, axes = plt.subplots(3, 2, figsize=(16, 17), constrained_layout=True)
    plot_waterfall(axes[0, 0], stats["waterfall"])
    plot_frequency_bins(axes[0, 1], "protein_frequency_bins.csv", "蛋白序列出现频次的长尾")
    plot_frequency_bins(axes[1, 0], "rna_frequency_bins.csv", "RNA 序列出现频次的长尾")
    plot_category_stack(axes[1, 1], "system_class_distribution.csv", "system class 构成（来源内部比例）", "来源")
    plot_category_stack(axes[2, 0], "rna_class_distribution.csv", "RNA class 构成（来源内部比例）", "来源")
    plot_organism(axes[2, 1])
    fig.suptitle("PRI-General 扩展分析：筛选、长尾与生物学组成", fontsize=20, fontweight="bold")
    fig.text(0.5, 0.01, "统计单位：复合物结构数；分类图仅统计最终 101,336 个复合物。", ha="center", fontsize=10, color="#444444")
    for index, ax in enumerate(axes.flat):
        ax.text(-0.08, 1.05, chr(ord("A") + index), transform=ax.transAxes, fontsize=15, fontweight="bold", va="top")
    for path in [DELIVERABLES / "PRI-General_extended_analysis_dashboard.png", DELIVERABLES / "PRI-General_extended_analysis_dashboard.pdf"]:
        fig.savefig(path, dpi=300 if path.suffix == ".png" else None, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(15, 11), constrained_layout=True)
    plot_evidence(axes[0, 0])
    plot_quality(axes[0, 1])
    plot_length_bins(axes[1, 0], "protein", "蛋白长度分箱（来源内部比例）")
    plot_length_bins(axes[1, 1], "RNA", "RNA 长度分箱（来源内部比例）")
    fig.suptitle("PRI-General 元数据与质量概览", fontsize=20, fontweight="bold")
    for index, ax in enumerate(axes.flat):
        ax.text(-0.08, 1.05, chr(ord("A") + index), transform=ax.transAxes, fontsize=15, fontweight="bold", va="top")
    for path in [DELIVERABLES / "PRI-General_metadata_quality.png", DELIVERABLES / "PRI-General_metadata_quality.pdf"]:
        fig.savefig(path, dpi=300 if path.suffix == ".png" else None, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    summary = {
        "final_complex_count": stats["total"],
        "unique_proteins": len(protein_counts),
        "unique_rnas": len(rna_counts),
        "source_counts": dict(stats["source_counts"]),
        "protein_frequency": {"top1": max(protein_counts.values()) if protein_counts else 0, "top10": sum(sorted(protein_counts.values(), reverse=True)[:10])},
        "rna_frequency": {"top1": max(rna_counts.values()) if rna_counts else 0, "top10": sum(sorted(rna_counts.values(), reverse=True)[:10])},
    }
    (OUT / "extended_analysis_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(DELIVERABLES / "PRI-General_extended_analysis_dashboard.png")
    print(DELIVERABLES / "PRI-General_metadata_quality.png")


if __name__ == "__main__":
    main()
