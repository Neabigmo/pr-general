"""Small shared helpers for the standalone PRI-General scripts."""

from __future__ import annotations

import hashlib
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd


SOURCE_LABELS = {
    "NPInter_v5": "NPInter 主表",
    "NPInter_v5_bindingsite": "NPInter binding-site",
    "NPInter_v5_miRNA": "NPInter miRNA",
    "RCSB_PDB": "RCSB PDB",
}
SOURCE_ORDER = ["NPInter 主表", "NPInter binding-site", "NPInter miRNA", "RCSB PDB"]
SOURCE_COLORS = {
    "NPInter 主表": "#4C78A8",
    "NPInter binding-site": "#F58518",
    "NPInter miRNA": "#54A24B",
    "RCSB PDB": "#B279A2",
}
RAW_SOURCE_DATABASES = {
    "npinter_main": "NPInter_v5",
    "npinter_bindingsite": "NPInter_v5_bindingsite",
    "mirna_protein": "NPInter_v5_miRNA",
    "rcsb_pdb": "RCSB_PDB",
}


def clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "null"} else text


def sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def pair_hash(protein: str, rna: str) -> str:
    return sequence_hash(f"{protein}|{rna}")


def complex_key(source_raw: str, pdb_id: str, pair_hash_value: str) -> str:
    if source_raw == "RCSB_PDB" and pdb_id:
        return f"PDB:{pdb_id}"
    return f"PAIR:{pair_hash_value}" if pair_hash_value else ""


@lru_cache(maxsize=8)
def load_pipeline_rules(root: Path) -> dict[str, Any]:
    """Read the existing pipeline policy once for this project root."""

    config = load_pipeline_config(root)
    protein_min, protein_max = config["protein_length"]
    rna_min, rna_max = config["rna_length"]
    return {
        "protein_length_min": int(protein_min),
        "protein_length_max": int(protein_max),
        "rna_length_min": int(rna_min),
        "rna_length_max": int(rna_max),
        "binding_site_window_max_nt": int(config["binding_site_window_max_nt"]),
        "protein_family_max_fraction": float(config["protein_family_max_fraction"]),
        "nominal_target_rows": int(config.get("nominal_target_rows", 115000)),
        "protein_caps": {
            key: int(value)
            for key, value in config.get("protein_caps", {
                "at_most_500": 500,
                "from_501_to_2000": 2000,
                "from_2001_to_10000": 8000,
                "above_10000": 10000,
            }).items()
        },
    }


@lru_cache(maxsize=8)
def load_pipeline_config(root: Path) -> dict[str, Any]:
    """Read the release-count settings used by the existing validator."""

    import yaml

    path = root.resolve() / "config" / "pipeline.yaml"
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def configure_plot_font() -> str:
    """Prefer an installed CJK font and fall back to matplotlib's default."""

    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    candidates = [
        Path(r"C:\Windows\Fonts\simhei.ttf"),
        Path(r"C:\Windows\Fonts\msyh.ttc"),
        Path(r"C:\Windows\Fonts\simsun.ttc"),
    ]
    family = "DejaVu Sans"
    for path in candidates:
        if path.exists():
            font_manager.fontManager.addfont(str(path))
            family = font_manager.FontProperties(fname=str(path)).get_name()
            break
    plt.rcParams.update({
        "font.family": family,
        "font.sans-serif": [family, "Microsoft YaHei", "SimHei", "Arial", "DejaVu Sans"],
        "axes.unicode_minus": False,
    })
    return family
