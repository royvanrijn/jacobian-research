#!/usr/bin/env python3
"""Fixed six-parent, two-fibre calibration. No ranking or population sweep."""
import argparse
import sys
from pathlib import Path
import certify_compact_r17_candidates as cert
import retained_native19_trial_v3 as engine
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import run, Limits

ROOT=Path(__file__).resolve().parents[2]
CAS=ROOT/'elliptic-curves/cas'
ART=ROOT/'artifacts/generated-results/elliptic-curves'
BATCH=ROOT/'artifacts/local/elliptic-curves/mestre-parent-calibration-v1'
SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'

def sources():
    names=['mestre_parent_calibration.py','prepare_mestre_parent_calibration.sage',
           'mestre_parent_adapter.py','icarm_curve245_mestre.py',
           'probe_mestre_fermigier_two_section_local_continuation.py',
           'mestre_root_tuples.py','nagao_1994.py','audit_recorded_point_mod2_rank_v3.py']
    return {**engine.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names}}

def freeze():
    path=BATCH/'intake-protocol.json'
    if path.exists():raise FileExistsError('preserve frozen calibration')
    files=[ART/n for n in ('mestre_component_coherent_sections_v1.json',
           'mestre_component_label_audit_v1.json','mestre_parent_portfolio_intake_v1.json',
           'mestre_parent_and_label_independent_v1.json','parent_portfolio_coverage_v1.json')]
    rows=[{'id':f'u{u}-'+('unit' if t=='1' else 'outer'), 'outer_u':str(u),
           'fibre_T':t,'initial_rank':11,'family':f'mestre-parent-u{u}'}
          for u in (11,13,17,19,23,29) for t in ('1','1009/101')]
    checkpoint(path,{'schema':'elliptic-curves.mestre-parent-calibration-intake.v1',
      'sources':sources(),'inputs':{str(p.relative_to(ROOT)):cert.hashed(p) for p in files},
      'rows':rows,'intake_wall_seconds':120,'rss_bytes':1610612736,'maximum_workers':1,
      'gate':'Six symbolically verified Q-distinct K3 parents have certified11-section lower bounds. Test practical point visibility and completed exposure on T=1 and one higher-height rational fibre per parent. These are separate parent and fibre parameters, not new high-rank claims.',
      'endpoint':'Certified rank gains above the own11 seed and completed point-search exposure. No score selection, validation-prime scoring, target curve input or parent superiority inference.',
      'plan':{'sample_size':2047,'sample_policy':'all nonzero parity masks in the own11 seed',
        'charts':49,'height':125000,'seconds_per_chart':10,'rank_stop':False,
        'geometry_wall_seconds':120,'worker_wall_seconds':600,'replay_wall_seconds':300,
        'maximum_point_boxes':588,'target_rank':32,'following_campaign':None},
      'failure_policy':'Preserve each failed/censored stage. No replacement, automatic retry or next wave. Every intake and all maps must pass before any point search.',
      'selection':'Fixed outer u=11,13,17,19,23,29; fixed T=1,1009/101. No data-dependent parameter selection. The latter has rational parameter height1009, but height depends on the chosen base coordinate.',
      'centre_policy':'384-bit canonical heights rounded at10^6, unimodular LLL and numerical CVP with exact parity/norm checks; largest49 computed norms among2047 masks. Numerical representatives are not certified closest vectors. Integral Gauss charts and hyperellred only. All twelve map files precede all point searches.'})
    print('FROZEN12 FIBRES; MAX588 BOXES',flush=True)

def intake_protocol():
    p=cert.read(BATCH/'intake-protocol.json')
    if p['sources']!=sources() or any(cert.hashed(ROOT/n)!=h for n,h in p['inputs'].items()):
        raise ArithmeticError('frozen intake sources or inputs changed')
    return p

def intake_all():
    p=intake_protocol();out=BATCH/'intake-ledger.json'
    if out.exists():raise FileExistsError('preserve intake ledger')
    ledger={'status':'RUNNING','rows':[]};checkpoint(out,ledger)
    for index,row in enumerate(p['rows']):
        d=BATCH/row['id']
        s=run([SAGE,str(CAS/'prepare_mestre_parent_calibration.sage'),'intake','--index',str(index)],
          limits=Limits(p['intake_wall_seconds'],p['rss_bytes']),log_path=d/'intake.log',
          checkpoint_path=d/'intake.supervisor.json',cwd=ROOT)
        ok=s['outcome']=='completed' and s['returncode']==0
        ledger['rows'].append({'id':row['id'],'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s})
        checkpoint(out,ledger);print(row['id'],'intake',s['outcome'],s['returncode'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger)
            raise ArithmeticError('all12 exact seed intakes required')
    ledger['status']='PASS';checkpoint(out,ledger)

def prepare():
    p=intake_protocol();out=BATCH/'protocol.json'
    if out.exists():raise FileExistsError('preserve search protocol')
    ledger=cert.read(BATCH/'intake-ledger.json')
    if ledger['status']!='PASS' or len(ledger['rows'])!=12:raise ArithmeticError('complete12 intake ledger required')
    seeds={}
    for row in p['rows']:
        f=BATCH/row['id']/'seed.json';seed=cert.read(f)
        if seed['status']!='PASS' or len(seed['points'])!=11:raise ArithmeticError('own11 seed required')
        seeds[str(f.relative_to(ROOT))]=cert.hashed(f)
    checkpoint(out,{**p,**p['plan'],'schema':'elliptic-curves.mestre-parent-calibration.v1',
      'intake_protocol_sha256':cert.hashed(BATCH/'intake-protocol.json'),'seed_hashes':seeds,
      'gp_sha256':cert.hashed(Path('/usr/bin/gp'))})
    print('ALL12 INTAKES PASS; SEARCH RULES FROZEN',flush=True)

def protocol():
    p=cert.read(BATCH/'protocol.json');intake_protocol()
    if p['intake_protocol_sha256']!=cert.hashed(BATCH/'intake-protocol.json') or any(cert.hashed(ROOT/n)!=h for n,h in p['seed_hashes'].items()):raise ArithmeticError('frozen seed inputs changed')
    return p

def masks(p):return list(range(1,1<<11))

def configure(index):
    global ROW,D,SEED
    ROW=protocol()['rows'][index];D=BATCH/ROW['id'];SEED=D/'seed.json'
    engine.ROW,engine.D,engine.SEED=ROW,D,SEED
    engine.protocol,engine.masks=protocol,masks

def launch():
    p=protocol();out=BATCH/'ledger.json'
    if out.exists():raise FileExistsError('preserve search ledger')
    ledger={'status':'RUNNING_GEOMETRY','maps':[],'rows':[]};checkpoint(out,ledger)
    for index,row in enumerate(p['rows']):
        configure(index)
        s=run([SAGE,str(CAS/'prepare_mestre_parent_calibration.sage'),'maps','--index',str(index)],
          limits=Limits(p['geometry_wall_seconds'],p['rss_bytes']),log_path=D/'maps.log',
          checkpoint_path=D/'maps.supervisor.json',cwd=ROOT)
        ok=s['outcome']=='completed' and s['returncode']==0
        ledger['maps'].append({'id':row['id'],'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s})
        checkpoint(out,ledger);print(row['id'],'maps',s['outcome'],s['returncode'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger)
            raise ArithmeticError('all maps required before points')
    ledger['status']='RUNNING_POINTS';checkpoint(out,ledger)
    for index,row in enumerate(p['rows']):
        configure(index);entry={'id':row['id'],'status':'RUNNING','stages':[]};ledger['rows'].append(entry);checkpoint(out,ledger)
        for stage in ('worker','replay'):
            s=run([sys.executable,str(Path(__file__).resolve()),stage,'--index',str(index)],
              limits=Limits(p[stage+'_wall_seconds'],p['rss_bytes']),log_path=D/(stage+'.log'),
              checkpoint_path=D/(stage+'.supervisor.json'),cwd=ROOT)
            ok=s['outcome']=='completed' and s['returncode']==0
            entry['stages'].append({'name':stage,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s})
            checkpoint(out,ledger);print(row['id'],stage,s['outcome'],s['returncode'],flush=True)
            if not ok:
                entry['status']=ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger)
                raise ArithmeticError('point worker/replay failed or censored')
        result=cert.read(D/'result.json');entry.update(status='PASS',rank_lower_bound=result['rank_lower_bound'],result_sha256=cert.hashed(D/'result.json'));checkpoint(out,ledger)
    ledger['status']='PASS';checkpoint(out,ledger)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['freeze','intake_all','prepare','launch','worker','replay']);parser.add_argument('--index',type=int);a=parser.parse_args()
    if a.stage in ('worker','replay'):configure(a.index);getattr(engine,a.stage)()
    else:globals()[a.stage]()
