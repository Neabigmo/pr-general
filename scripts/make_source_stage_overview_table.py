#!/usr/bin/env python
"""生成各来源经过各轮筛选后的数字总览表。"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common import RAW_SOURCE_DATABASES, SOURCE_LABELS


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "deliverables"


def main() -> None:
    audit_cols = [
        "source_database",
        "sequence_hash",
        "pdb_id",
        "stage_d1",
        "stage_d2",
        "stage_d3",
        "stage_protein_capped",
        "stage_family_capped",
        "stage_queue",
    ]
    stage_names = [
        "input",
        "stage_d1",
        "stage_d2",
        "stage_d3",
        "stage_protein_capped",
        "stage_family_capped",
        "stage_queue",
    ]
    source_keys: dict[str, dict[str, set[str]]] = {}
    total_keys: dict[str, set[str]] = {stage: set() for stage in stage_names}

    def complex_key(frame: pd.DataFrame) -> pd.Series:
        """Return one key per complex, never one key per technical row.

        PDB complexes are identified by unique ``pdb_id``. Other sources use
        the sequence-pair hash because they do not have an experimental
        structure identifier.
        """
        source_name = frame["source_database"].fillna("").astype(str)
        pdb_id = frame["pdb_id"].fillna("").astype(str).str.strip()
        pair_hash = frame["sequence_hash"].fillna("").astype(str).str.strip()
        keys = pair_hash.where(pair_hash.ne(""), "MISSING")
        pdb_mask = source_name.eq("RCSB_PDB") & pdb_id.ne("")
        keys.loc[pdb_mask] = "PDB:" + pdb_id.loc[pdb_mask]
        return keys

    for chunk in pd.read_csv(ROOT / "data" / "final" / "audit" / "general_stage_audit.csv", usecols=audit_cols, chunksize=100_000, low_memory=False):
        keys = complex_key(chunk)
        for source_name, group_index in chunk.groupby("source_database", dropna=False).groups.items():
            source_name = str(source_name)
            source_keys.setdefault(source_name, {stage: set() for stage in stage_names})
            group_keys = keys.loc[group_index]
            source_keys[source_name]["input"].update(group_keys.tolist())
            total_keys["input"].update(group_keys.tolist())
            for stage in stage_names[1:]:
                mask = chunk.loc[group_index, stage].fillna(False).astype(bool)
                selected = group_keys.loc[mask]
                source_keys[source_name][stage].update(selected.tolist())
                total_keys[stage].update(selected.tolist())

    rows = []
    source_order = ["npinter_main", "npinter_bindingsite", "mirna_protein", "rcsb_pdb"]
    for source_id in source_order:
        source_name = RAW_SOURCE_DATABASES[source_id]
        stage = source_keys[source_name]
        rows.append([
            SOURCE_LABELS[source_name],
            int(len(stage["input"])),
            int(len(stage["stage_d1"])),
            int(len(stage["stage_d2"])),
            int(len(stage["stage_d3"])),
            int(len(stage["stage_protein_capped"])),
            int(len(stage["stage_family_capped"])),
            int(len(stage["stage_queue"])),
        ])
    total = ["合计"] + [len(total_keys[stage]) for stage in stage_names]

    headers = [
        "来源",
        "可解析\n复合物",
        "D1\n严格 QC",
        "D2\n完全去重",
        "D3\nRNA90 去冗余",
        "单蛋白\n限额",
        "Protein30 家族\n限额",
        "最终\n复合物",
    ]
    table_rows = [headers] + rows + [total]
    table_frame = pd.DataFrame(table_rows[1:], columns=headers)
    table_frame.to_csv(OUT / "PRI-General_source_stage_overview_latest.csv", index=False, encoding="utf-8-sig")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Microsoft YaHei", "SimHei", "Arial", "DejaVu Sans"],
        "axes.unicode_minus": False,
    })
    fig, ax = plt.subplots(figsize=(18, 5.2))
    ax.axis("off")
    ax.set_title("表：各来源经过各轮筛选后保留的复合物数", fontsize=20, pad=18, color="#222222")
    cell_text = [[str(value) if not isinstance(value, int) else f"{value:,}" for value in row] for row in table_rows]
    table = ax.table(
        cellText=cell_text,
        colWidths=[0.19, 0.11, 0.11, 0.11, 0.13, 0.11, 0.13, 0.10],
        cellLoc="center",
        loc="center",
        bbox=[0.01, 0.23, 0.98, 0.60],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10.5)
    for (row_index, col_index), cell in table.get_celld().items():
        cell.set_edgecolor("#777777")
        cell.set_linewidth(0.8)
        cell.set_facecolor("#FFFFFF")
        cell.PAD = 0.012
        cell.get_text().set_color("#222222")
        cell.get_text().set_va("center")
        if col_index == 0:
            cell.get_text().set_ha("left")
        if row_index == 0:
            cell.set_facecolor("#F1F1F1")
            cell.get_text().set_weight("bold")
        if row_index == len(table_rows) - 1:
            cell.set_facecolor("#F7F7F7")
            cell.get_text().set_weight("bold")
    ax.text(
        0.01,
        0.10,
        "说明：本表只统计复合物，不统计技术记录。PDB 按唯一 pdb_id 计数；其他来源按唯一蛋白–RNA序列组合计数。无法补齐两端序列的来源记录不进入筛选。",
        transform=ax.transAxes,
        fontsize=11.5,
        color="#555555",
        ha="left",
        va="center",
    )
    fig.savefig(OUT / "PRI-General_source_stage_overview_latest.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / "PRI-General_source_stage_overview_latest.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    pd.DataFrame({"stage": stage_names, "complexes": [len(total_keys[stage]) for stage in stage_names]}).to_csv(
        OUT / "PRI-General_complex_stage_waterfall_latest.csv", index=False, encoding="utf-8-sig"
    )
    print(OUT / "PRI-General_source_stage_overview_latest.png")


if __name__ == "__main__":
    main()
