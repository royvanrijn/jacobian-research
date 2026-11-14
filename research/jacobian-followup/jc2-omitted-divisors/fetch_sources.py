"""Explicit source refresh; never run by verify.py."""
from pathlib import Path
import hashlib
import json
from urllib.request import Request, urlopen

HERE = Path(__file__).resolve().parent/'sources'
SOURCES = [
    ('fonollosa-2609.33491v1.html', 'https://arxiv.org/html/2609.33491v1',
     'Theorem1.3; Sections5.2-5.4, equations15-19: finite normalization, source open and every omitted prime'),
    ('stacks-normal-031T.html', 'https://stacks.math.columbia.edu/tag/031T',
     'Normal Noetherian domain is the intersection of its height-one localizations'),
]


def main():
    HERE.mkdir(exist_ok=True)
    rows = []
    for name, url, use in SOURCES:
        with urlopen(Request(url, headers={'User-Agent': 'jacobian-research source audit'}), timeout=25) as response:
            data = response.read()
        (HERE/name).write_bytes(data)
        rows.append(dict(file=name, url=url, use=use, bytes=len(data),
                         sha256=hashlib.sha256(data).hexdigest()))
    manifest = dict(as_of='2026-11-14', sources=rows,
                    prior_plane_import_manifest='../../jc2-boundary-attachment/sources/manifest.json',
                    theorem_imports_independently_reproved=False)
    (HERE/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(rows, indent=2))


if __name__ == '__main__':
    main()
