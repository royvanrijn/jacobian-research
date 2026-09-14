#!/usr/bin/env python3
"""Retain and byte-check scheduler evidence; never run a point search."""
import argparse
import gzip
import io
import json
from pathlib import Path
import tarfile

from cancellation_scheduler_prepare import OUT, RAW
from finite_cancellation_corpus import ROOT, digest, write


def verify():
    manifest=json.loads((OUT/'replay-manifest.json').read_text())
    archive=OUT/'replay-bundle.tar.gz'
    if digest(archive.read_bytes())!=manifest['archive_sha256']:
        raise ArithmeticError('archive hash differs')
    with tarfile.open(archive,'r:gz') as tar:
        members=tar.getmembers()
        if len(members)!=len(manifest['files_sha256']):raise ArithmeticError('wrong member count')
        seen=set()
        for member in members:
            name=member.name
            if name in seen or not member.isfile() or Path(name).is_absolute() or '..' in Path(name).parts:
                raise ArithmeticError('unsafe or duplicate archive member')
            seen.add(name)
            if digest(tar.extractfile(member).read())!=manifest['files_sha256'].get(name):
                raise ArithmeticError('member hash differs')
        if seen!=set(manifest['files_sha256']):raise ArithmeticError('missing archive member')
    print(json.dumps({'status':'PASS_BYTE_COMPLETE_BUNDLE','files':len(seen),'bytes':archive.stat().st_size}),flush=True)


def pack():
    archive=OUT/'replay-bundle.tar.gz'
    if archive.exists() or (OUT/'replay-manifest.json').exists():raise FileExistsError('preserve archive')
    plan=json.loads((OUT/'protocol.json').read_text())
    if json.loads((OUT/'supervision.json').read_text())['status']!='COMPLETE':
        raise ArithmeticError('primary execution not complete')
    files={p for p in RAW.rglob('*') if p.is_file() and p.name!='supervisor.lock'}
    files|={p for p in OUT.iterdir() if p.is_file()}
    files|={ROOT/name for name in plan['source_sha256']}
    supplemental={}
    if (OUT/'dependency-audit.json').exists():
        supplemental=json.loads((OUT/'dependency-audit.json').read_text())['extra_source_sha256']
        files|={ROOT/name for name in supplemental}
    files.add(Path(__file__).resolve())
    names={str(p.relative_to(ROOT)):p for p in files}
    for name,h in {**plan['source_sha256'],**supplemental}.items():
        if digest(names[name].read_bytes())!=h:raise ArithmeticError('frozen source changed')
    hashes={}
    with archive.open('wb') as f:
        with gzip.GzipFile(fileobj=f,mode='wb',mtime=0,filename='') as gz:
            with tarfile.open(fileobj=gz,mode='w') as tar:
                for name,path in sorted(names.items()):
                    raw=path.read_bytes();hashes[name]=digest(raw)
                    info=tarfile.TarInfo(name);info.size=len(raw);info.mode=0o644;info.mtime=0
                    tar.addfile(info,io.BytesIO(raw))
    write(OUT/'replay-manifest.json',{'status':'RETAINED_REPLAY_INPUTS','files_sha256':hashes,
        'archive_sha256':digest(archive.read_bytes()),'protocol_sha256':digest((OUT/'protocol.json').read_bytes()),
        'boundary':'Byte-complete raw primary and separately labelled development receipts, frozen inputs and implementation sources. Integrity checking is not an arithmetic proof; independent arithmetic receipts are inside the archive. Extract only into a separate replay workspace.'})
    verify()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify',action='store_true')
    verify() if parser.parse_args().verify else pack()
