"""Download released ENCODE eCLIP reproducible peak BED files.

The ENCODE search is restricted to released eCLIP ``bed narrowPeak`` files
whose biological replicate field is ``[1, 2]``.  Those are the reproducible
two-replicate peak assets; RBNS/RNAcompete matrices are not used here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import Request, urlopen


PORTAL = "https://www.encodeproject.org"
SEARCH_URL = (
    f"{PORTAL}/search/?type=File&assay_title=eCLIP&status=released&"
    "file_format=bed&format=json&limit=all&frame=object"
)


def request_json(url: str) -> dict:
    request = Request(url, headers={"Accept": "application/json", "User-Agent": "PRI-General eCLIP downloader"})
    with urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, path: Path) -> None:
    request = Request(url, headers={"User-Agent": "PRI-General eCLIP downloader"})
    with urlopen(request, timeout=120) as response, path.open("wb") as handle:
        while chunk := response.read(1024 * 1024):
            handle.write(chunk)


def fetch_record(item: dict, peaks_dir: Path) -> dict:
    accession = str(item["accession"])
    path = peaks_dir / f"{accession}.bed.gz"
    expected_md5 = item.get("md5sum") or ""
    status = "reused"
    if not path.exists() or (expected_md5 and md5(path) != expected_md5):
        download(urljoin(PORTAL, item["href"]), path)
        status = "downloaded"
    observed_md5 = md5(path)
    if expected_md5 and observed_md5 != expected_md5:
        raise RuntimeError(f"MD5 mismatch for {accession}: {observed_md5} != {expected_md5}")
    return {
        "accession": accession,
        "dataset": item.get("dataset"),
        "target": item.get("target"),
        "assembly": item.get("assembly"),
        "file_format": item.get("file_format"),
        "file_format_type": item.get("file_format_type"),
        "output_type": item.get("output_type"),
        "biological_replicates": item.get("biological_replicates"),
        "file_size": item.get("file_size"),
        "md5sum": expected_md5,
        "observed_md5": observed_md5,
        "href": item.get("href"),
        "local_path": f"peaks/{accession}.bed.gz",
        "download_status": status,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("data/downloads/encode_eclip"))
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    peaks_dir = output / "peaks"
    peaks_dir.mkdir(parents=True, exist_ok=True)

    graph = request_json(SEARCH_URL).get("@graph", [])
    selected = [
        item for item in graph
        if item.get("status") == "released"
        and item.get("file_format") == "bed"
        and item.get("file_format_type") == "narrowPeak"
        and sorted(item.get("biological_replicates") or []) == [1, 2]
        and not item.get("no_file_available", False)
    ]
    selected.sort(key=lambda item: item.get("accession", ""))
    manifest = {
        "source": "ENCODE",
        "assay": "eCLIP",
        "selection_query": SEARCH_URL,
        "selection_rule": "released + bed narrowPeak + biological_replicates=[1,2]",
        "files": [],
    }
    workers = max(1, min(int(args.workers), 16))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(fetch_record, item, peaks_dir) for item in selected]
        for future in as_completed(futures):
            manifest["files"].append(future.result())
        (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    manifest["files"].sort(key=lambda item: item["accession"])
    (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        "output_dir": str(output),
        "selected_files": len(selected),
        "downloaded_files": sum(item["download_status"] == "downloaded" for item in manifest["files"]),
        "reused_files": sum(item["download_status"] == "reused" for item in manifest["files"]),
        "compressed_bytes": sum(int(item.get("file_size") or 0) for item in manifest["files"]),
        "assemblies": {assembly: sum(item.get("assembly") == assembly for item in manifest["files"]) for assembly in sorted({item.get("assembly") for item in manifest["files"]})},
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
