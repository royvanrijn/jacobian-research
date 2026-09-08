#!/usr/bin/env python3
"""Small self-contained bundle for the latest8 exact incidence checks."""
from pathlib import Path
import zipfile
from hashlib import sha256
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/latest8-cross-family-v1';OUT=ART/'latest8_cross_family_evidence_v1.json'

def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve incidence bundle')
    proof=ART/'latest8_cross_family_j_incidence_v1.json';replay=ART/'latest8_cross_family_j_incidence_replay_v1.json'
    if cert.read(replay)['status']!='PASS':raise ArithmeticError('exact incidence replay missing')
    paths={Path(__file__).resolve(),proof,replay,ROOT/'elliptic-curves/cas/replay_latest8_cross_family_incidence.py',ROOT/'elliptic-curves/cas/verify_latest8_incidence_bundle.py',ROOT/'elliptic-curves/notes/LATEST_EIGHT_CROSS_FAMILY_INCIDENCE_2026-09-06.md'};paths.update(ROOT/p for p in cert.read(proof)['sources']);paths.update(p for p in D.rglob('*') if p.is_file());paths=dependencies(paths);members=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(paths):
            raw=p.read_bytes();name=str(p.relative_to(ROOT));info=zipfile.ZipInfo(name,date_time=(2026,9,6,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16;z.writestr(info,raw);members.append({'path':name,'bytes':len(raw),'sha256':sha256(raw).hexdigest()})
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest()!=r['sha256'] for r in members):raise ArithmeticError('incidence archive integrity failure')
    cert.write(OUT,{'schema':'elliptic-curves.latest8-incidence-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'files':members,'claim_boundary':'Self-contained exact96-pair rational j-incidence replay. The inventory supplies the target equations; this bundle does not itself replay their rank proofs. Other families and additional points remain outside scope.'});print('PACKAGED LATEST8 INCIDENCE',len(members),archive.stat().st_size,flush=True)
if __name__=='__main__':main()
