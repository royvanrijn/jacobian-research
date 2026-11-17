"""Explicit refresh of the Orevkov primary source; no network in checker."""
from pathlib import Path
import hashlib
import json
import urllib.request

HERE = Path(__file__).resolve().parent/'sources'
URL = 'https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf'


def main():
    request = urllib.request.Request(URL, headers={'User-Agent': 'omitted-divisor-research/1.0'})
    with urllib.request.urlopen(request, timeout=25) as response:
        data = response.read()
    assert data.startswith(b'%PDF') and len(data) > 50000
    name = 'orevkov-three-sheeted-jc86.pdf'
    (HERE/name).write_bytes(data)
    manifest = dict(as_of='2026-11-14', sources=[dict(file=name, url=URL, bytes=len(data),
        sha256=hashlib.sha256(data).hexdigest(),
        use='Lemma4.2 global Euler identity; Lemma5.2 no-jump boundary germ has immersive image; Section4 contractions identify the affine Stein normalization')],
        imported_theorems_independently_reproved=False,
        source_visual_pages_zero_based=[6, 7, 8])
    (HERE/'residue_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(dict(file=name, bytes=len(data), sha256=manifest['sources'][0]['sha256'])))


if __name__ == '__main__':
    main()
