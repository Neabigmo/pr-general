from __future__ import annotations

from pathlib import Path
from typing import Any


class ConfigError(ValueError):
    pass


def _reject_generated_counts(value: Any, path: str = "config") -> None:
    if isinstance(value, dict):
        for key, child in value.items():
            if str(key).lower() in {"counts", "statistics", "generated_counts"}:
                raise ConfigError(f"generated count/statistics key is not allowed in {path}.{key}")
            _reject_generated_counts(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _reject_generated_counts(child, f"{path}[{index}]")


def load_config(path: str | Path) -> dict[str, Any]:
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover
        raise ConfigError("PyYAML is required to load the v2 policy") from exc
    config_path = Path(path)
    data = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ConfigError("config root must be a mapping")
    _reject_generated_counts(data)
    if data.get("schema_version") != "pri-general/v2":
        raise ConfigError("config schema_version must be pri-general/v2")
    for key in ("scope", "evidence", "structure", "splits", "caps", "paths"):
        if not isinstance(data.get(key), dict):
            raise ConfigError(f"missing config section: {key}")
    return data
