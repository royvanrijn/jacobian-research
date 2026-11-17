#!/usr/bin/env python3
"""Retain a fixed author-repository snapshot; no search or remote mutation."""

import hashlib
import json
from pathlib import Path
import urllib.request

HERE = Path(__file__).resolve().parent
COMMIT = "74fe52e57de6f91988f63ac105cbc157534ad966"
FILES = {
    "README.md": "f39c9a71b11eb2f66e4f59613d1c457eb3170053",
    "Wirt_Hm_Suite_Python/excel_reader.py": "ef07b921bbe1020b616fe1399a5b3351f5ee394d",
    "Wirt_Hm_Suite_Python/d_hm.py": "448ef01ee5e9a8c9b28e8adc2fcdaccad070ad58",
    "Wirt_Hm_Suite_Python/calc_wirt.py": "243be617f89bd2421342bbe7362f5049e86ecec9",
    "Wirt_Hm_Suite_Python/sym_hm.py": "f16b1b673d92c73da9ccc85046cb344812aa6b20",
    "all_data_A.xlsx": "3882c6da34c9057f8f13da31afb9d8973b7158cb",
}


def main():
    directory = HERE / "sources/wirt-hm"
    directory.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name, expected_blob in FILES.items():
        url = f"https://raw.githubusercontent.com/ThisSentenceIsALie/Wirt_Hm/{COMMIT}/{name}"
        path = directory / name.replace("/", "__")
        if not path.exists():
            with urllib.request.urlopen(url, timeout=30) as response:
                data = response.read()
            actual = hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()
            assert actual == expected_blob, f"Git blob mismatch: {name}"
            path.write_bytes(data)
        data = path.read_bytes()
        assert hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest() == expected_blob
        manifest.append({"path": path.name, "repository_path": name, "url": url,
                         "git_blob_sha1": expected_blob, "sha256": hashlib.sha256(data).hexdigest(),
                         "bytes": len(data)})
    receipt = {"commit": COMMIT, "files": manifest, "binary_fallback": "GitHub connector does not support binary downloads; read-only pinned raw URL used for XLSX"}
    (directory / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"commit": COMMIT, "files": len(manifest), "bytes": sum(x["bytes"] for x in manifest)}, indent=2))


if __name__ == "__main__":
    main()
