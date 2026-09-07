#!/usr/bin/env python3
"""Standalone exact histories, point clouds and equation comparisons for64 curves."""
from hashlib import sha256
from pathlib import Path
import zipfile
import certify_compact_r17_candidates as cert
import full11952_late64_r17_pari_batch as batch
from package_recorded_mod2_audit import dependencies
ROOT=batch.ROOT;CAS=batch.CAS;ART=batch.ART
OUT=ART/'full11952_late64_point_evidence_v1.json'

def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve standalone full64 point evidence')
    p=batch.protocol();done=cert.read(batch.D/'post-batch/ledger.json');summary=cert.read(ART/'full11952_late64_experiment_v1.json')
    if done['status']!='PASS' or summary['status']!='PASS' or len(p['rows'])!=64:raise ArithmeticError('terminal64 proof cohort with no unresolved odd-prime upgrade required')
    for path in [ROOT/'artifacts/local/elliptic-curves/full11952-late64-controller-v1/points.supervisor.json',ROOT/'artifacts/local/elliptic-curves/full11952-late64-controller-v1/stream.supervisor.json']:
        s=cert.read(path)
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('all source jobs must be terminal before packaging')
    paths={Path(__file__).resolve(),CAS/'verify_full11952_late64_points_portable.py',CAS/'certify_full11952_late64_r17_results.py',CAS/'audit_recorded_point_mod2_rank_v3.py',CAS/'audit_retained_cloud_modl.py',ART/'full11952_late64_r17_results_v1.json',ART/'full11952_late64_experiment_v1.json',batch.extension.OUT,batch.extension.D/'protocol.json',batch.extension.D/'fresh-validation.json',batch.extension.D/'selection-replay.supervisor.json',batch.extension.D/'validation-check.supervisor.json',batch.CONTROL/'ledger.json',ROOT/'artifacts/local/elliptic-curves/native11952-height125-control-v1/125000/verification.json'}
    # Original search logs and all archive states are immutable after terminal supervision.
    paths.update(q for q in batch.D.rglob('*') if q.is_file() and '__pycache__' not in q.parts)
    proofs=[ART/'full11952_late64_r17_results_v1.json']
    for row in p['rows']:
        proofs.extend(ART/('full11952_late64_r17_'+row['id'].replace('-','_')+'_'+suffix+'_v1.json') for suffix in ('mod2','modl'))
        early=ART/('full11952_late64_r17_'+row['id'].replace('-','_')+'_early_modl_v1.json')
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
    cert.write(OUT,{'schema':'elliptic-curves.full11952-late64-point-evidence.v1','archive':str(archive.relative_to(ROOT)),'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'files':members,'required_base_archives':[],'claim_boundary':'Standalone replay of all64 exact admission histories, complete point clouds modulo2,3,5, rational quartic maps, generic-section transports, and independent point proofs with embedded593-catalogue and853-prior-equation comparisons. Frozen score selection and control gates are hash-bound context; the original20.89-billion scanner, million-row scoring,4096 scalar-score trial, cache construction and historical controls are not rerun by this point-only supplement. No new search, exact rank, full saturation, point absence or universal novelty.'})
    print('STANDALONE LATE64 POINT BUNDLE',len(members),'FILES',archive.stat().st_size,'BYTES',flush=True)

if __name__=='__main__':main()
