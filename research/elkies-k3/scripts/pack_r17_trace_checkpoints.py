#!/usr/bin/env python3
"""Package a fresh67-trace constructor run for the portable proof replayers."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import zipfile


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-dir',type=Path,required=True);a=p.parse_args()
    records=sorted((a.input_dir/'traces').glob('*.json'))
    expected=[f'{i:03d}.json' for i in range(67)]
    if [f.name for f in records]!=expected:raise ValueError('The complete67-trace checkpoint roster is required.')
    archive=a.input_dir/'trace-records.zip';manifest=a.input_dir/'trace-records-manifest.json'
    if archive.exists() or manifest.exists():raise FileExistsError('Preserve the existing archive and manifest; use a fresh output directory.')
    rows=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for f in records:
            raw=f.read_bytes();name=f.relative_to(a.input_dir).as_posix()
            info=zipfile.ZipInfo(name,date_time=(2026,9,12,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED
            z.writestr(info,raw,compress_type=zipfile.ZIP_DEFLATED,compresslevel=9)
            rows.append({'path':name,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)})
    with manifest.open('x') as f:
        json.dump({'archive_sha256':sha256(archive.read_bytes()).hexdigest(),'members':rows},f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps({'archive_bytes':archive.stat().st_size,'traces':67}))


if __name__=='__main__':main()
