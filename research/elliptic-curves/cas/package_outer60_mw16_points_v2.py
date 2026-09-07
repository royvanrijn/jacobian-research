#!/usr/bin/env python3
"""Standalone exact histories, point clouds and equation comparisons for60 curves."""
from hashlib import sha256
from pathlib import Path
import zipfile
import certify_compact_r17_candidates as cert
import outer60_mw16_pari_batch as batch
from package_recorded_mod2_audit import dependencies
ROOT=batch.ROOT;CAS=batch.CAS;ART=batch.ART
OUT=ART/'outer60_mw16_point_evidence_v2.json'

def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve standalone outer60 point evidence')
    p=batch.protocol();done=cert.read(batch.D/'post-batch/ledger.json');summary=cert.read(ART/'outer60_mw16_experiment_v1.json')
    if done['status']!='PASS' or summary['status']!='PASS' or len(p['rows'])!=60:raise ArithmeticError('terminal60 proof cohort with no unresolved odd-prime upgrade required')
    for path in [ROOT/'artifacts/local/elliptic-curves/outer60-mw16-controller-v1/points.supervisor.json',ROOT/'artifacts/local/elliptic-curves/outer60-mw16-controller-v1/verify.supervisor.json']:
        s=cert.read(path)
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('all source jobs must be terminal before packaging')
    control=ROOT/'artifacts/local/elliptic-curves/outer60-mw16-controller-v1'
    paths={Path(__file__).resolve(),CAS/'verify_outer60_mw16_points_portable_v2.py',CAS/'certify_outer60_mw16_results.py',CAS/'replay_outer60_mw16_geometry.py',CAS/'audit_recorded_point_mod2_rank_v3.py',CAS/'audit_retained_cloud_modl.py',ART/'outer60_mw16_results_v1.json',ART/'outer60_mw16_experiment_v1.json',batch.extension.OUT,batch.extension.D/'protocol.json',batch.extension.D/'controller/ledger.json',batch.extension.D/'fresh-validation.json',control/'protocol.json',control/'ledger.json'}
    # Original search logs and all archive states are immutable after terminal supervision.
    paths.update(q for q in batch.D.rglob('*') if q.is_file() and '__pycache__' not in q.parts)
    proofs=[ART/'outer60_mw16_results_v1.json']
    for row in p['rows']:
        proofs.extend(ART/('outer60_mw16_'+row['id'].replace('-','_')+'_'+suffix+'_v1.json') for suffix in ('mod2','modl'))
        early=ART/('outer60_mw16_'+row['id'].replace('-','_')+'_early_modl_v1.json')
        if early.exists():proofs.append(early)
    for data in [p,*[cert.read(q) for q in proofs]]:
        for name,h in data['sources'].items():
            q=ROOT/name
            if cert.hashed(q)!=h:raise ArithmeticError('bound point-proof source changed: '+name)
            paths.add(q)
    paths.update(proofs);paths=dependencies(paths);members=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for q in sorted(paths):
            raw=q.read_bytes();name=str(q.relative_to(ROOT));info=zipfile.ZipInfo(name,date_time=(2026,9,6,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16
            z.writestr(info,raw);members.append({'path':name,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)})
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest()!=r['sha256'] for r in members):raise ArithmeticError('standalone point evidence integrity differs')
    cert.write(OUT,{'schema':'elliptic-curves.outer60-mw16-point-evidence.v2','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'files':members,'required_base_archives':[],'supersedes_failed_verifier':'elliptic-curves/cas/verify_outer60_mw16_points_portable.py','correction':'The v1 verifier addressed an obsolete outer60-mw16-r17-pari-v1 directory and failed before any logical proof stage; v2 uses the recorded outer60-mw16-pari-v1 protocol. Original searches, proofs and v1 evidence remain unchanged.','claim_boundary':'Standalone replay of all60 exact admission histories, complete point clouds modulo2,3,5, rational quartic maps, generic-section transports, and independent point proofs with embedded593-catalogue and917-prior-equation comparisons. Frozen parameter/score selection and validation are hash-bound context; the286812899-address annular scan,15360 scalar-score trial, canonical table construction and historical controls are not rerun by this point-only supplement. No new search, exact rank, full saturation, point absence or universal novelty.'})
    print('STANDALONE OUTER60 POINT BUNDLE',len(members),'FILES',archive.stat().st_size,'BYTES',flush=True)

if __name__=='__main__':main()
