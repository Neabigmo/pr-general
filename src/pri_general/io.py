from __future__ import annotations

import csv
import json
import zipfile
from pathlib import Path
from typing import Iterable, Iterator, Mapping


class SourceReadError(RuntimeError):
    pass


def iter_table(path: str | Path) -> Iterator[dict]:
    source = Path(path)
    if not source.exists():
        raise SourceReadError(f"source does not exist: {source}")
    suffix = source.suffix.lower()
    if suffix == ".parquet":
        try:
            import pandas as pd
            frame = pd.read_parquet(source, engine="pyarrow")
        except Exception as exc:
            raise SourceReadError(f"cannot read parquet source {source}: {exc}") from exc
        yield from frame.to_dict(orient="records")
        return
    if suffix in {".jsonl", ".ndjson"}:
        with source.open(encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    yield json.loads(line)
        return
    if suffix == ".json":
        value = json.loads(source.read_text(encoding="utf-8"))
        if isinstance(value, list):
            yield from value
        elif isinstance(value, dict):
            yield value
        else:
            raise SourceReadError(f"JSON source must contain an object or list: {source}")
        return
    if suffix == ".zip":
        try:
            archive = zipfile.ZipFile(source)
            members = [item for item in archive.infolist() if not item.is_dir()]
            if len(members) != 1:
                raise SourceReadError(f"zip source must contain exactly one table: {source}")
            member = members[0]
            delimiter = "\t" if member.filename.lower().endswith((".tsv", ".txt")) else ","
            with archive, archive.open(member) as binary:
                handle = __import__("io").TextIOWrapper(binary, encoding="utf-8", errors="replace")
                yield from csv.DictReader(handle, delimiter=delimiter)
        except SourceReadError:
            raise
        except Exception as exc:
            raise SourceReadError(f"cannot read zipped table {source}: {exc}") from exc
        return
    if suffix in {".csv", ".tsv", ".txt"}:
        delimiter = "\t" if suffix in {".tsv", ".txt"} else ","
        with source.open(encoding="utf-8", newline="") as handle:
            yield from csv.DictReader(handle, delimiter=delimiter)
        return
    raise SourceReadError(f"unsupported table format: {source.suffix}")


def write_json(path: str | Path, value: Mapping) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_jsonl(path: str | Path, rows: Iterable[Mapping]) -> int:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with destination.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(dict(row), ensure_ascii=False, sort_keys=True) + "\n")
            count += 1
    return count
