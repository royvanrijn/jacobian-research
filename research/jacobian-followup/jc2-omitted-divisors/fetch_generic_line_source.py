"""Explicit network refresh of the fixed Chau degree-ratio source; never run by replay."""
from pathlib import Path
import hashlib
import json
from urllib.request import Request, urlopen

HERE=Path(__file__).resolve().parent
URL='https://arxiv.org/html/0711.3894v1'


def main():
    data=urlopen(Request(URL,headers={'User-Agent':'jacobian-research source pin'}),timeout=30).read()
    assert b'0711.3894' in data and b'Theorem' in data
    source=HERE/'sources/chau-0711.3894v1.html'
    source.write_bytes(data)
    manifest=dict(as_of='2026-11-14',sources=[dict(path=source.name,url=URL,
        sha256=hashlib.sha256(data).hexdigest(),bytes=len(data),
        imported_result='Theorem3: every nonproper image component of a complex plane Keller map has a polynomial parametrization with coordinate-degree ratio degP/degQ',
        excluded_imports='No simple-polynomial hypothesis or inverse conclusion from Theorem2 is applied')])
    (HERE/'sources/generic_line_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'source_bytes':len(data),'source_sha256':manifest['sources'][0]['sha256']}))


if __name__=='__main__':
    main()
