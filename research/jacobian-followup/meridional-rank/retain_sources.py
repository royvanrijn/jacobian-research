#!/usr/bin/env python3
"""Retain the main primary inputs with hashes; preserve existing bytes."""

import hashlib
import json
from pathlib import Path
import urllib.request

HERE = Path(__file__).resolve().parent
SOURCES = {
    "BBK21": ("bbk21.pdf", "https://arxiv.org/pdf/1907.02982"),
    "BKM24": ("bkm24.html", "https://arxiv.org/html/2208.09032"),
    "BKP26A": ("fishnets-v1.html", "https://arxiv.org/html/2411.10642"),
    "GGLW26": ("whitehead-v1.html", "https://arxiv.org/html/2611.15673v1"),
    "CESY14": ("quandle-colourings-v2.html", "https://arxiv.org/html/1312.3307v2"),
    "MUSICK12": ("musick-bridge-v3.pdf", "https://arxiv.org/pdf/1208.4233v3"),
    "KL04": ("rational-tangles-v2.html", "https://arxiv.org/html/math/0311499v2"),
    "BBKM23": ("many-twigs.pdf", "https://arxiv.org/pdf/2008.00763"),
    "MORIAH07": ("moriah-heegaard-2007.pdf", "https://msp.org/gtm/2007/12/gtm-2007-12-007s.pdf"),
    "STW24": ("core-groups-v2.html", "https://arxiv.org/html/2311.02048v2"),
    "BC13": ("burness-crestani-direct-powers.pdf", "https://seis.bristol.ac.uk/~tb13602/docs/BC_1.pdf"),
    "PS17": ("poudel-saveliev-2017.pdf", "https://msp.org/agt/2017/17-5/agt-v17-n5-p03-p.pdf"),
    "DUT22": ("dutra-whitehead-v1.pdf", "https://arxiv.org/pdf/2212.13081v1"),
    "COR14": ("cornwell-knot-groups-v5.pdf", "https://arxiv.org/pdf/1303.4943v5"),
    "COR17": ("cornwell-kch-v3.pdf", "https://arxiv.org/pdf/1310.7526v3"),
    "GGLW26PDF": ("whitehead-v1.pdf", "https://arxiv.org/pdf/2611.15673v1"),
    "ATLAS-A5": ("atlas-a5.html", "https://brauer.maths.qmul.ac.uk/Atlas/alt/A5/"),
    "ATLAS-L32": ("atlas-l32.html", "https://brauer.maths.qmul.ac.uk/Atlas/lin/L32/"),
}


def main():
    directory = HERE / "sources"
    directory.mkdir(exist_ok=True)
    manifest = []
    for key, (name, url) in SOURCES.items():
        target = directory / name
        row = {"id": key, "url": url, "retrieved_utc_date": "2026-11-16"}
        try:
            if not target.exists():
                request = urllib.request.Request(url, headers={"User-Agent": "Research source audit"})
                with urllib.request.urlopen(request, timeout=25) as response:
                    blob = response.read()
                if name.endswith(".pdf") and not blob.startswith(b"%PDF"):
                    raise ValueError("response is not a PDF")
                if name.endswith(".html") and b"<html" not in blob[:2000].lower():
                    raise ValueError("response is not HTML")
                target.write_bytes(blob)
            blob = target.read_bytes()
            row.update({"file": name, "bytes": len(blob), "sha256": hashlib.sha256(blob).hexdigest()})
        except Exception as error:
            row["unavailable"] = str(error)
        manifest.append(row)
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
