#!/usr/bin/env python3
"""Standalone exact model and independent-point proofs for five rank26 curves."""
from hashlib import sha256
from pathlib import Path
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'compact192_five_rank26_evidence_v1.json'


def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve standalone five-curve evidence')
    proof=cert.read(ART/'compact192_rank26_models_v1.json')
    if proof['status']!='PASS' or len(proof['curves'])!=5:raise ArithmeticError('complete five-curve certificate required')
    paths={Path(__file__).resolve(),ROOT/'elliptic-curves/cas/verify_compact192_five_rank26.py',
           ROOT/'elliptic-curves/cas/export_compact192_rank26_models.py',
           ROOT/'elliptic-curves/notes/COMPACT192_FIVE_RANK26_CURVES_2026-09-06.md',
           ART/'compact192_rank26_models_v1.json',ART/'new_compact192_rank26_curves.sage'}
    for name,digest in proof['sources'].items():
        path=ROOT/name
        if cert.hashed(path)!=digest:raise ArithmeticError('five-curve proof source changed')
        paths.add(path)
    paths=dependencies(paths)
    members=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for path in sorted(paths):
            raw=path.read_bytes();name=str(path.relative_to(ROOT));info=zipfile.ZipInfo(name,date_time=(2026,9,6,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16
            z.writestr(info,raw);members.append({'path':name,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)})
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest()!=r['sha256'] for r in members):raise ArithmeticError('standalone archive integrity differs')
    cert.write(OUT,{'schema':'elliptic-curves.compact192-five-rank26-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'files':members,'required_base_archives':[],
                   'claim_boundary':'Standalone exact five-curve proof supplement requiring no base archive. Replays global minimality, exact transports and26-point independence for each curve, mutual rational nonisomorphism and absence from the embedded593-equation catalogue and549 prior equations, plus exact Sage-export agreement and execution. The full192 search history, score selection, trace tables, chronology and other curve proofs are not replayed by this small bundle. No new point search, exact rank, conductor, full saturation, record or universal novelty claim.'})
    print('STANDALONE FIVE26 BUNDLE',len(members),'files',archive.stat().st_size,'bytes',flush=True)


if __name__=='__main__':main()
