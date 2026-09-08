#!/usr/bin/env python3
"""Bind the enlarged-sample trial and all196 own-curve witness charts."""
import argparse,json,sys
from pathlib import Path
import larger_specialized_parity_trial as trial
import audit_larger_specialized_parity as sample
import audit_larger_specialized_parity_trial as audit
import certify_compact_r17_candidates as cert
import replay_retention24_geometry as geometry
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=trial.ROOT;ART=trial.ART;D=trial.BATCH;CAS=trial.CAS
OUT=ART/'larger_specialized_parity_trial_v1.json'
HISTORIES=[ROOT/'artifacts/local/elliptic-curves/full11952-64-r17-pari-v1/11952-0962587/result.json',
    ROOT/'artifacts/local/elliptic-curves/full11952-specialized-followup-v1/new-20260906-186/result.json',
    ROOT/'artifacts/local/elliptic-curves/full11952-million-height-v1/new-20260906-186/result.json',
    D/'new-20260906-186/result.json']
UNION=D/'own-union.json';MOD2=ART/'larger_parity_own196_union_mod2_v1.json';MODL=ART/'larger_parity_own196_union_modl_v1.json'
def inputs():
    p=trial.protocol();ledger=cert.read(D/'ledger.json')
    if ledger['status']!='PASS' or len(ledger['rows'])!=2 or any(r['status']!='PASS' for r in ledger['rows']):raise ArithmeticError('both declared trials must have succeeded')
    if sample.expected()!=cert.read(sample.OUT):raise ArithmeticError('exact65536 sample replay differs')
    for i in range(2):audit.main(i,True)
    records=[cert.read(path) for path in HISTORIES];last=records[-1]
    if any(any(r[k]!=last[k] for k in ('family','parameter','curve')) or len(r['charts'])!=49 for r in records):raise ArithmeticError('four fixed same-curve49-chart histories required')
    for row,path in zip(ledger['rows'],[D/'native28/result.json',HISTORIES[-1]]):
        if row['result_sha256']!=cert.hashed(path):raise ArithmeticError('terminal trial binding differs')
    paths=[Path(__file__).resolve(),Path(audit.__file__).resolve(),Path(geometry.__file__).resolve(),sample.OUT,
        D/'protocol.json',D/'ledger.json',*HISTORIES,
        ART/'full11952_million_new_20260906_186_coverage_v1.json',
        ART/'larger_parity_native28_coverage_v1.json',ART/'larger_parity_new_20260906_186_coverage_v1.json']
    sources={str(q.relative_to(ROOT)):cert.hashed(q) for q in paths}
    return p,records,sources
def payload(records,sources):
    last=records[-1]
    return {'status':'POINT_ONLY_CONCATENATION_NOT_AN_ADMISSION_TRANSCRIPT','family':last['family'],'parameter':last['parameter'],
        'curve':last['curve'],'final_state':last['final_state'],
        'charts':[{'search':{'finite_curve_points':c['search']['finite_curve_points']}} for r in records for c in r['charts']],
        'sources':sources,'claim_boundary':'Witness union only; original exact histories are preserved separately. No synthetic admission history or new point search.'}
def build():
    if OUT.exists() or UNION.exists():raise FileExistsError('preserve enlarged-sample union proof')
    p,records,sources=inputs();checkpoint(UNION,payload(records,sources))
    for name,script,args in [('union-mod2','audit_recorded_point_mod2_rank_v3.py',['--input',str(UNION),'--input-sha256',cert.hashed(UNION),'--output',str(MOD2),'--prime-bound','997']),('union-mod2-check','audit_recorded_point_mod2_rank_v3.py',['--check',str(MOD2)]),('union-modl','audit_retained_cloud_modl.py',['--input',str(MOD2),'--output',str(MODL)]),('union-modl-check','audit_retained_cloud_modl.py',['--check',str(MODL)])]:
        s=run([sys.executable,str(CAS/script),*args],limits=Limits(300,1610612736),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('union proof failed/censored')
    checkpoint(OUT,expected())
def expected():
    p,records,sources=inputs();u=payload(records,sources)
    if cert.read(UNION)!=u:raise ArithmeticError('own196 union provenance differs')
    geometry.cloud_check(records[-1],u['charts'],MOD2,UNION);cloud=cert.read(MOD2);odd=cert.read(MODL)
    control=cert.read(D/'native28/result.json');old=cert.read(ROOT/'artifacts/local/elliptic-curves/native28-specialized-parity-control-v1/maps.json')
    target=control['charts'][-1]['centre']['parity'];old_masks={r['parity'] for r in old['sample']}
    native=cert.read(ART/'larger_parity_native28_coverage_v1.json');own=cert.read(ART/'larger_parity_new_20260906_186_coverage_v1.json')
    if native['mod2_lower_bound']!=29 or native['odd_modulus_lower_bounds']!={'3':29,'5':29}:raise ArithmeticError('control gate not certified')
    return {'schema':'elliptic-curves.larger-specialized-parity-trial.v1','status':'PASS','sources':sources,
        'certificates':{str(q.relative_to(ROOT)):cert.hashed(q) for q in (UNION,MOD2,MODL)},
        'control':{'completed_boxes':native['completed_boxes'],'lower_bound':29,'target_chart_parity':target,
            'target_class_absent_from_original2048':target not in old_masks},
        'own':{'completed_boxes':own['completed_boxes'],'trial_retained_points':own['retained_points'],
            'union_chart_inputs':196,'union_retained_points':len(cloud['points']),'mod2_lower_bound':cloud['rank_lower_bound'],
            'odd_modulus_lower_bounds':{str(r['modulus']):r['finite_column_rank'] for r in odd['audits']}},
        'claim_boundary':p['boundaries']+' The recovery chart class was checked retrospectively only after both searches froze and finished. All196 own-curve chart inputs have an exact point-union audit. No isolated portable supplement replay is yet claimed.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args()
    if a.check:
        if cert.read(OUT)!=expected():raise ArithmeticError('larger-sample report differs')
    else:build()
    d=cert.read(OUT);print('LARGER PARITY TRIAL',d['control'],d['own'],flush=True)
