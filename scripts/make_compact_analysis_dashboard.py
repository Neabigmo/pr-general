"""Generate a compact, structure-level PRI-General analysis dashboard.

The script consumes the current final-complex analysis CSV/JSON outputs.  It
does not modify the dataset or rerun filtering.  All counts are complex
counts; sequence lengths are summarized at the unique component level used by
the existing analysis.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from common import SOURCE_COLORS, SOURCE_ORDER, configure_plot_font


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = ROOT / "reports" / "final" / "analysis"
OUT_DIR = ROOT / "deliverables" / "analysis"


def load_inputs() -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    summary = json.loads(
        (ANALYSIS_DIR / "final_complex_analysis_summary.json").read_text(encoding="utf-8")
    )
    top20 = pd.read_csv(ANALYSIS_DIR / "protein_top20_source_before_after.csv")
    lengths = pd.read_csv(ANALYSIS_DIR / "complex_component_length_summary_by_source.csv")
    return summary, top20, lengths


def add_panel_label(ax, label: str) -> None:
    ax.text(
        -0.08,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=15,
        fontweight="bold",
        va="top",
    )


def plot_source_composition(ax, summary: dict) -> None:
    counts = summary["complex_count_by_source"]
    labels = [s for s in SOURCE_ORDER if s in counts]
    values = [counts[s] for s in labels]
    total = summary["final_complex_count"]
    y = np.arange(len(labels))
    bars = ax.barh(y, values, color=[SOURCE_COLORS[s] for s in labels], height=0.62)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("最终复合物数")
    ax.set_title("最终 PRI-General 的来源组成", pad=10)
    ax.grid(axis="x", alpha=0.22)
    ax.set_axisbelow(True)
    for bar, value in zip(bars, values):
        ax.text(
            value + total * 0.008,
            bar.get_y() + bar.get_height() / 2,
            f"{value:,}  ({value / total:.1%})",
            va="center",
            fontsize=10,
        )
    ax.set_xlim(0, max(values) * 1.23)
    add_panel_label(ax, "A")


def plot_top20(ax, top20: pd.DataFrame) -> None:
    data = top20.sort_values("final_complex_count", ascending=True).copy()
    y = np.arange(len(data))
    height = 0.36
    ax.barh(y + height / 2, data["source_complex_count"], height=height, color="#B7C7D9", label="筛选前来源池")
    ax.barh(y - height / 2, data["final_complex_count"], height=height, color="#E87817", label="最终复合物")
    ax.set_yticks(y, data["protein_name"])
    ax.set_xlabel("复合物数")
    ax.set_title("高频蛋白：来源池与最终保留量", pad=10)
    ax.grid(axis="x", alpha=0.22)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.tick_params(axis="y", labelsize=8.5)
    add_panel_label(ax, "B")


def plot_length_summary(ax, lengths: pd.DataFrame, molecule: str, title: str, xlabel: str) -> None:
    data = lengths[lengths["molecule"] == molecule].copy()
    data["source"] = pd.Categorical(data["source"], SOURCE_ORDER, ordered=True)
    data = data.sort_values("source")
    y = np.arange(len(data))
    colors = [SOURCE_COLORS[s] for s in data["source"].astype(str)]
    # p10-p90: light line; p25-p75: thick line; median: dot.
    for yi, (_, row), color in zip(y, data.iterrows(), colors):
        ax.plot([row["p10"], row["p90"]], [yi, yi], color=color, linewidth=2.2, alpha=0.55)
        ax.plot([row["p25"], row["p75"]], [yi, yi], color=color, linewidth=8, solid_capstyle="butt")
        ax.scatter(row["p50"], yi, s=50, color="#222222", zorder=3)
        ax.text(row["p50"], yi + 0.19, f"中位数 {row['p50']:.0f}", ha="center", va="bottom", fontsize=8)
    ax.set_yticks(y, data["source"].astype(str))
    ax.set_xlabel(xlabel)
    ax.set_title(title, pad=10)
    ax.grid(axis="x", alpha=0.22)
    ax.set_axisbelow(True)
    ax.legend(
        handles=[
            Patch(facecolor="#888888", alpha=0.55, label="P10–P90"),
            Patch(facecolor="#888888", label="P25–P75"),
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#222222", label="中位数", markersize=6),
        ],
        frameon=False,
        fontsize=8,
        loc="lower right",
    )
    add_panel_label(ax, "C" if molecule == "protein" else "D")


def main() -> None:
    configure_plot_font()
    summary, top20, lengths = load_inputs()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(16, 11), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, wspace=0.22, hspace=0.25)
    plot_source_composition(fig.add_subplot(grid[0, 0]), summary)
    plot_top20(fig.add_subplot(grid[0, 1]), top20)
    plot_length_summary(
        fig.add_subplot(grid[1, 0]),
        lengths,
        "protein",
        "各来源蛋白长度（最终复合物中的蛋白组件）",
        "蛋白长度（aa）",
    )
    plot_length_summary(
        fig.add_subplot(grid[1, 1]),
        lengths,
        "RNA",
        "各来源 RNA 长度（最终复合物中的 RNA 组件）",
        "RNA 长度（nt）",
    )
    fig.suptitle(
        f"PRI-General 数据分析总览（{summary['final_complex_count']:,} 个最终复合物）",
        fontsize=20,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.01,
        "统计单位：复合物结构数；长度图展示组件序列的分位数（点为中位数，粗线为 P25–P75，细线为 P10–P90）。",
        ha="center",
        fontsize=10,
        color="#444444",
    )
    png_path = OUT_DIR / "PRI-General_analysis_dashboard.png"
    pdf_path = OUT_DIR / "PRI-General_analysis_dashboard.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    findings = pd.DataFrame(
        [
            {
                "metric": "final_complex_count",
                "value": summary["final_complex_count"],
                "unit": "complexes",
                "note": "最终复合物总数",
            },
            {
                "metric": "unique_proteins_in_complexes",
                "value": summary["unique_proteins_in_complexes"],
                "unit": "proteins",
                "note": "最终复合物中出现的不同蛋白",
            },
            {
                "metric": "unique_rnas_in_complexes",
                "value": summary["unique_rnas_in_complexes"],
                "unit": "RNAs",
                "note": "最终复合物中出现的不同 RNA",
            },
            {
                "metric": "protein_top1_complex_fraction",
                "value": summary["protein_top1_complex_fraction"],
                "unit": "fraction",
                "note": "出现频次最高的蛋白占最终复合物比例",
            },
            {
                "metric": "protein_top10_complex_fraction",
                "value": summary["protein_top10_complex_fraction"],
                "unit": "fraction",
                "note": "前十个高频蛋白合计占最终复合物比例",
            },
        ]
    )
    for source, count in summary["complex_count_by_source"].items():
        findings.loc[len(findings)] = {
            "metric": f"source_complex_count::{source}",
            "value": count,
            "unit": "complexes",
            "note": f"{source} 最终复合物数",
        }
    findings.to_csv(OUT_DIR / "PRI-General_analysis_key_findings.csv", index=False, encoding="utf-8-sig")
    print(png_path)
    print(pdf_path)


if __name__ == "__main__":
    main()
