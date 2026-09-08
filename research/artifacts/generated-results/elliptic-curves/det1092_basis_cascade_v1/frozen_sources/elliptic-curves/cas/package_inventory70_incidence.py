#!/usr/bin/env python3
"""Self-contained inputs for all840 incidence checks and the duplicate transport."""
from pathlib import Path
from hashlib import sha256
import zipfile
import certify_compact_r17_candidates as cert
from certify_inventory70_incidence import ROOT,ART,PROOFS,TRANSPORT,INDEX
from package_recorded_mod2_audit import dependencies
OUT=ART/'inventory70_incidence_evidence_v1.json'
CHECKERS=['replay_compact_cross_family_incidence.py','replay_latest7_cross_family_incidence.py','replay_latest8_cross_family_incidence.py','replay_latest23_cross_family_incidence.py','audit_compact_published_r17_transport_v3.sage','certify_inventory70_incidence.py','verify_inventory70_incidence_bundle.py']
def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve incidence evidence')
    paths={Path(__file__).resolve(),INDEX,TRANSPORT,*PROOFS,ART/'inventory70_cross_family_incidence_v1.json',ROOT/'elliptic-curves/notes/INVENTORY70_CROSS_FAMILY_INCIDENCE_2026-09-06.md'}
    paths.update(ROOT/'elliptic-curves/cas'/n for n in CHECKERS)
    for p in [*PROOFS,TRANSPORT]:paths.update(ROOT/n for n in cert.read(p)['sources'])
    paths.update(p for p in (ROOT/'artifacts/local/elliptic-curves/latest23-cross-family-v1').rglob('*') if p.is_file())
    paths=dependencies(paths);members=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(paths):
            raw=p.read_bytes();name=str(p.relative_to(ROOT));info=zipfile.ZipInfo(name,date_time=(2026,9,6,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16;z.writestr(info,raw);members.append({'path':name,'bytes':len(raw),'sha256':sha256(raw).hexdigest()})
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest()!=r['sha256'] for r in members):raise ArithmeticError('archive integrity differs')
    cert.write(OUT,{'schema':'elliptic-curves.inventory70-incidence-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'files':members,'claim_boundary':'Exact incidence and generic transport inputs for70 curves840 pairs; independent rank proofs are outside this small supplement.'});print('PACKAGED INVENTORY70 INCIDENCE',len(members),archive.stat().st_size,flush=True)
if __name__=='__main__':main()
