#!/usr/bin/env python3
"""Standalone point evidence for a censor-aware matched score comparison."""
from hashlib import sha256
from pathlib import Path
import zipfile
import certify_compact_r17_candidates as cert
import strata60_mw16_pari_batch as batch
from package_recorded_mod2_audit import dependencies
ROOT=batch.ROOT;CAS=batch.CAS;ART=batch.ART
OUT=ART/'strata60_mw16_point_evidence_v1.json'
CONTROL=batch.extension.D/'point-controller'
AUDIT=batch.extension.D/'accounting-controller'


def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve standalone comparison bundle')
    p=batch.protocol();report_path=ART/'retained_mw16_score_strata_experiment_v1.json'
    audit_path=ART/'retained_mw16_score_strata_accounting_replay_v1.json'
    report=cert.read(report_path);audit=cert.read(audit_path)
    if cert.read(CONTROL/'ledger.json')['status']!='COMPLETE_FIXED_COMPARISON_AND_ACCOUNTING' or cert.read(AUDIT/'ledger.json')['status']!='PASS' or audit['status']!='PASS':
        raise ArithmeticError('terminal comparison and independent accounting required')
    if len(p['rows'])!=60 or [r['id'] for r in p['rows']]!=[r['id'] for r in report['rows']]:raise ArithmeticError('all60 allocated curves required')
    paths={Path(__file__).resolve(),CAS/'verify_strata60_mw16_points_portable.py',report_path,audit_path,
           CONTROL/'protocol.json',CONTROL/'ledger.json',AUDIT/'protocol.json',AUDIT/'ledger.json',
           batch.extension.OUT,batch.extension.D/'protocol.json',batch.extension.D/'controller/ledger.json'}
    paths.update(q for q in batch.D.rglob('*') if q.is_file() and '__pycache__' not in q.parts)
    for data in (p,report,audit):
        for name,h in data['sources'].items():
            q=ROOT/name
            if cert.hashed(q)!=h:raise ArithmeticError('bound evidence changed: '+name)
            paths.add(q)
    cases=[]
    for i,row in enumerate(report['rows']):
        folder=batch.D/row['id'];case={k:row[k] for k in ('id','arm','family','parameter','rank_lower_bound','certified_gain','worker_status','verification_status','completed_boxes','attempted_boxes')}
        case['index']=i;case['replay_point_proof']=row['verification_status']=='PASS'
        if case['replay_point_proof']:
            verified=cert.read(folder/'verification.json')
            if verified['status']!='PASS' or verified['rank_lower_bound']!=row['rank_lower_bound']:raise ArithmeticError('certified report differs')
            for key,hashkey in [('cloud_path','cloud_sha256'),('odd_path','odd_sha256')]:
                path=ROOT/verified[key]
                if cert.hashed(path)!=verified[hashkey]:raise ArithmeticError('point certificate changed')
                paths.add(path)
                for name,h in cert.read(path)['sources'].items():
                    q=ROOT/name
                    if cert.hashed(q)!=h:raise ArithmeticError('point source changed')
                    paths.add(q)
            case.update(cloud_path=verified['cloud_path'],odd_path=verified['odd_path'])
        elif row['certified_gain'] is not None or row['rank_lower_bound'] is not None:
            raise ArithmeticError('unverified row claims rank')
        cases.append(case)
    paths=dependencies(paths);members=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for path in sorted(paths):
            raw=path.read_bytes();name=str(path.relative_to(ROOT));info=zipfile.ZipInfo(name,date_time=(2026,9,6,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16
            z.writestr(info,raw);members.append({'path':name,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)})
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest()!=r['sha256'] for r in members):raise ArithmeticError('archive integrity differs')
    cert.write(OUT,{'schema':'elliptic-curves.strata60-point-evidence.v1','status':'PASS',
        'builder_sha256':cert.hashed(Path(__file__).resolve()),'archive':str(archive.relative_to(ROOT)),
        'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'files':members,
        'required_base_archives':[],'allocated_curves':60,'cases':cases,
        'certified_rows':sum(r['replay_point_proof'] for r in cases),
        'unresolved_rows':sum(not r['replay_point_proof'] for r in cases),
        'claim_boundary':'Standalone admission-history, exact generic16 transport, rational geometry, returned-point provenance and finite-quotient rank replay for every certified row of the60-curve matched score comparison. All allocations, unresolved outcomes, terminal partial prefixes and immutable computation-accounting evidence remain included. An unresolved row is not a zero gain. The prior population scan and matching replay are bound context, not rerun. The local independent accounting audit is embedded; isolated point replay does not rerun that audit or measure search time again. No new search, point absence, exact rank, catalogue novelty or rank-density theorem.'})
    print('STANDALONE STRATA60',len(members),'FILES',archive.stat().st_size,'BYTES',flush=True)


if __name__=='__main__':main()
