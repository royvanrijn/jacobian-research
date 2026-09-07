#!/usr/bin/env python3
"""Portable exact j-incidence, generic section transport and precision diagnostics."""
import argparse
from hashlib import sha256
from pathlib import Path
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies
ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT/'artifacts/local/elliptic-curves'

def package(output):
    archive = output.with_suffix('.zip')
    if output.exists() or archive.exists():
        raise FileExistsError('preserve prior incidence evidence')
    directory = LOCAL/'compact-cross-family-incidence-v1'
    if cert.read(directory/'replay-result.json')['status'] != 'PASS' or [r['status'] for r in cert.read(directory/'transport-v3-results.json')['rows']] != ['PASS','PASS']:
        raise ArithmeticError('exact incidence or transport replay incomplete')
    names = ['compact_cross_family_j_incidence_v1.json','compact_cross_family_j_incidence_replay_v1.json',
             'compact_published_r17_generic_transport_v1.json','explicit_height_precision_audit_v1.json']
    paths = {Path(__file__).resolve()}
    for name in names:
        p = ART/name; paths.add(p)
        paths.update(ROOT/s for s in cert.read(p)['sources'])
    for folder in (directory,LOCAL/'explicit-height-precision-v1'):
        paths.update(p for p in folder.iterdir() if p.is_file())
    paths.update(ROOT/'elliptic-curves/cas'/name for name in ('audit_compact_published_r17_transport.sage','audit_compact_published_r17_transport_v2.sage','audit_compact_published_r17_transport_v3.sage'))
    paths = dependencies(paths)
    rows = []
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for path in sorted(paths):
            raw = path.read_bytes(); name = str(path.relative_to(ROOT))
            info = zipfile.ZipInfo(name,date_time=(2026,9,5,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16
            z.writestr(info,raw);rows.append({'path':name,'bytes':len(raw),'sha256':sha256(raw).hexdigest()})
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest()!=r['sha256'] for r in rows):
            raise ArithmeticError('archive member verification failed')
    cert.write(output,{'schema':'elliptic-curves.cross-family-incidence-evidence.v1',
        'builder_sha256':cert.hashed(Path(__file__).resolve()),'archive':str(archive.relative_to(ROOT)),
        'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'files':rows,
        'claim_boundary':'Exact384-pair projective j-incidence and unimodular generic transport proofs, with separate numerical height-precision comparison. Failed numerical proposal adapters and logs are retained. No new curve or rank gain.'})
    print('PACKAGED CROSS-FAMILY INCIDENCE',len(rows),'files;',archive.stat().st_size,'bytes',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    package(p.parse_args().output.resolve())
