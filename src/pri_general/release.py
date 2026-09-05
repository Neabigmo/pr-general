from __future__ import annotations

import hashlib
from pathlib import Path


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_checksums(paths: list[str | Path], destination: str | Path, base_dir: str | Path | None = None) -> None:
    lines = []
    base = Path(base_dir).resolve() if base_dir is not None else None
    for path in sorted((Path(item) for item in paths), key=lambda item: str(item)):
        label = path.resolve().relative_to(base).as_posix() if base is not None else path.as_posix()
        lines.append(f"{sha256_file(path)}  {label}")
    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
