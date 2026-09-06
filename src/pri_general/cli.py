from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import load_config
from .pipeline import build, load_source_manifest, report
from .validate import validate_source_manifest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="pri-general")
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name in ("build", "validate", "report"):
        command = subparsers.add_parser(name)
        command.add_argument("--config", default="config/pri_general_v2.yaml")
        command.add_argument("--release", default="v2.0.0")
        if name == "validate":
            command.add_argument("--strict", action="store_true")
            command.add_argument("--framework", action="store_true", help="validate policy and source shape without release data")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    root = Path.cwd()
    config = load_config(args.config)
    if args.command == "build":
        print(json.dumps(build(root, args.config, args.release), ensure_ascii=False, indent=2))
        return 0
    manifest = load_source_manifest(root / config["paths"]["source_manifest"])
    issues = validate_source_manifest(manifest, root, strict=getattr(args, "strict", False))
    if args.command == "validate":
        release_manifest = root / config["paths"]["release_root"] / args.release / "manifest.json"
        release_data = None
        if release_manifest.exists():
            release_data = json.loads(release_manifest.read_text(encoding="utf-8"))
        if not args.framework:
            if release_data is None:
                issues.append(f"release manifest missing: {release_manifest}")
            elif release_data.get("status") != "READY":
                issues.append(f"release status is {release_data.get('status')}")
            issues.extend(f"release blocker: {item}" for item in (release_data or {}).get("blockers", []))
            if release_data is not None:
                print(json.dumps(release_data, ensure_ascii=False, indent=2))
        if issues:
            print(json.dumps({"status": "INCOMPLETE", "issues": issues}, ensure_ascii=False, indent=2))
            return 1
        print(json.dumps({"status": "OK", "issues": []}, ensure_ascii=False, indent=2))
        return 0
    result = report(root, args.config, args.release)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0
