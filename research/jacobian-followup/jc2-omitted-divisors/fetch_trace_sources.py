"""Explicit source refresh; the mathematical checker never uses the network."""
from pathlib import Path
import hashlib
import json
import urllib.request

HERE = Path(__file__).resolve().parent / 'sources'
ROWS = [
    ('stacks-trace-0BVH.html', 'https://stacks.math.columbia.edu/tag/0BVH',
     'Trace pairing and discriminant of a finite locally free algebra'),
    ('stacks-dual-0BW0.html', 'https://stacks.math.columbia.edu/tag/0BW0',
     'Complementary module, trace dual and Dedekind different'),
    ('stacks-jacobian-0BWB.html', 'https://stacks.math.columbia.edu/tag/0BWB',
     'Different of a flat quasifinite local complete intersection'),
    ('wood-binary-1007.5508v1.html', 'https://arxiv.org/html/1007.5508v1',
     'Sections 2.1 and 3: binary-form ring basis, discriminant and geometric construction'),
]


def main():
    sources = []
    HERE.mkdir(exist_ok=True)
    for name, url, use in ROWS:
        request = urllib.request.Request(url, headers={'User-Agent': 'omission-research/1.0'})
        with urllib.request.urlopen(request, timeout=25) as response:
            data = response.read()
            resolved = response.url
        assert len(data) > 5000, (name, len(data))
        (HERE / name).write_bytes(data)
        sources.append(dict(file=name, url=url, resolved_url=resolved, use=use,
                            bytes=len(data), sha256=hashlib.sha256(data).hexdigest()))
        print(name, len(data), flush=True)
    manifest = dict(as_of='2026-11-14', sources=sources,
                    imports_independently_reproved=False,
                    scope='Trace foundations and binary-form order; local arguments are written in DISCRIMINANT_LEDGER.md')
    (HERE / 'trace_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
