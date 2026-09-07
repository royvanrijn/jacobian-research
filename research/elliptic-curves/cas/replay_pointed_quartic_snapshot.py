#!/usr/bin/env python3
"""Run the historical MW16 checker against its self-contained source bundle.

The working tree now uses a different search API. Extract the hash-checked
executed sources into an isolated temporary root, then run the original
read-only checker there. Never weaken its source equality checks or rewrite
historical certificates to describe the new backend.
"""
import argparse
import base64
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from pointed_quartic_migration import ROOT, REGRESSION_REVISION


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle", type=Path, required=True)
    p.add_argument("--summary", type=Path, required=True)
    args = p.parse_args()
    bundle_path, summary_path = args.bundle.resolve(), args.summary.resolve()
    if sha256(bundle_path.read_bytes()).hexdigest() != json.loads(summary_path.read_text())["bundle_sha256"]:
        raise ArithmeticError("bundle differs from the historical summary")
    bundle = json.loads(gzip.decompress(bundle_path.read_bytes()))
    with tempfile.TemporaryDirectory(prefix="pointed-regression-") as directory:
        root = Path(directory)
        for name, record in bundle["files"].items():
            relative = Path(name)
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("unsafe snapshot path")
            data = base64.b64decode(record["base64"]) if "base64" in record else record["text"].encode()
            if sha256(data).hexdigest() != record["sha256"]:
                raise ArithmeticError("snapshot content checksum mismatch: "+name)
            destination = root/relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(data)
        checker = "elliptic-curves/cas/check_mw16_sensitivity_policy.py"
        # This additional policy checker postdates the calibration run, so
        # its original pinned version may not itself be in the source bundle.
        if not (root/checker).exists():
            (root/checker).write_bytes(subprocess.check_output(
                ["git", "show", REGRESSION_REVISION+":"+checker], cwd=ROOT))
        subprocess.run([sys.executable, str(root/checker), "--bundle", str(bundle_path),
                        "--summary", str(summary_path), "--check"], cwd=root, check=True)


if __name__ == "__main__":
    main()
