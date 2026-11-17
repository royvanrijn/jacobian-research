#!/usr/bin/env python3
"""Retain the pinned primary October 2026 release; verify Git object hashes."""

import hashlib
import json
from pathlib import Path
import urllib.request

HERE = Path(__file__).resolve().parent
COMMIT = "06cbbb7e2737e1a11085b64e0ead021d63c550c3"
FILES = [{"path":"README.md","blob":"cf773fa99860bb8df7bfcac169f647bccf26a809"},{"path":"data/certificates.json","blob":"543fac7e04c5807b2ccce2682aee124418c85739"},{"path":"data/still_open.json","blob":"01c7c332a0bf1565d7558aec4c769ff36cbd31a9"},{"path":"search/alexq.py","blob":"158f2a5af646155d779162c5e85d53e3719c33b7"},{"path":"search/alexq_exact.py","blob":"e292958ac5c979c6fb959cbe333a2b13097252ad"},{"path":"verify/verify.py","blob":"71db0e88b2df13b45fe24b345dce14b5dfd45fcc"}]


def main():
    target = HERE / "sources/current13"
    target.mkdir(exist_ok=True)
    manifest = {"repository": "anaxia-research/conjecture-to-choose",
                "commit": COMMIT, "commit_date": "2026-10-22T11:28:09Z", "files": []}
    for item in FILES:
        path = target / item["path"]
        url = f"https://raw.githubusercontent.com/{manifest['repository']}/{COMMIT}/{item['path']}"
        data = path.read_bytes() if path.exists() else urllib.request.urlopen(url, timeout=30).read()
        blob = hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()
        assert blob == item["blob"], item["path"]
        path.parent.mkdir(exist_ok=True, parents=True)
        if not path.exists():
            path.write_bytes(data)
        manifest["files"].append({**item, "sha256": hashlib.sha256(data).hexdigest(),
                                  "bytes": len(data), "url": url})
    output = target / "manifest.json"
    if output.exists():
        assert json.loads(output.read_text()) == manifest
    else:
        output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"status": "PINNED", "commit": COMMIT, "files": len(FILES)}))


if __name__ == "__main__":
    main()
