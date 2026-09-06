#!/usr/bin/env python3
"""Run one predeclared own-subgroup follow-up after compact192 inventory replay."""
import argparse
import sys
import time
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run, Limits

ROOT=Path(__file__).resolve().parents[2]
CAS=ROOT/'elliptic-curves/cas'
ART=ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL=ROOT/'artifacts/local/elliptic-curves'
D=LOCAL/'compact192-specialized-followup-controller-v1'


def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve one follow-up controller')
    names=['finish_compact192_specialized_followup.py','compact192_specialized_followup.py',
           'prepare_compact192_specialized_followup.sage','audit_compact192_specialized_followup.py',
           'replay_compact192_specialized_geometry.sage']
    paths=[*(CAS/n for n in names),ART/'native28_specialized_parity_adaptive_coverage_v1.json',
           LOCAL/'compact192-evidence-finalization-v1/protocol.json']
    checkpoint(D/'protocol.json',{
        'schema':'elliptic-curves.compact192-specialized-followup-controller.v1',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'wait_deadline_seconds':120000,'prepare_seconds':600,
        'maximum_roster':24,'trial_seconds_per_curve':1800,
        'audit_seconds_per_curve':960,'geometry_seconds':1200,
        'rss_bytes':4294967296,
        'scope':'After the fixed compact192 initial cohort has passed catalogue, full-cloud and complete V13 inventory replay, freeze all newly added V13 curves from that cohort with lower bound26or27. The entire eligible roster must have1through24 entries; no truncation or refill. This explicit single follow-up changes point exposure using their own additional independent directions. The previously certified known29 control is the visibility gate. Freeze all2048-mask numerical samples and49 maps per curve before any new point search, then execute at most1176 boxes at height125000 and ten seconds each with one point worker and target28. Perform exact histories, full-cloud proofs modulo2,3,5 and rational geometry. Any failure or censoring stops the controller without retry. No second wave, automatic height increase, inventory promotion, exact-rank or universal-novelty claim.'})


def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen follow-up controller sources differ')
    out=D/'ledger.json'
    if out.exists():raise FileExistsError('preserve one follow-up execution')
    ledger={'status':'WAITING_FOR_INVENTORY_REPLAY','rows':[]};checkpoint(out,ledger)
    deadline=time.monotonic()+p['wait_deadline_seconds']
    proof=ART/'new_high_rank_curve_index_v13_memory_replay_v1.json'
    while not proof.exists():
        upstream=cert.read(LOCAL/'compact192-evidence-finalization-v1/ledger.json')
        if upstream['status'] not in ('WAITING_FOR_FIXED_JOBS','RUNNING','PASS') or any(r['status']!='PASS' for r in upstream['rows']):raise ArithmeticError('upstream inventory evidence failure requires review')
        if time.monotonic()>deadline:raise TimeoutError('fixed inventory replay wait deadline')
        time.sleep(5)
    ledger['status']='RUNNING';checkpoint(out,ledger)
    def stage(name,executable,script,args,seconds):
        path=D/(name+'.supervisor.json')
        if path.exists():raise FileExistsError('preserve one follow-up stage')
        result=run([executable,str(CAS/script),*args],limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=path,cwd=ROOT)
        ok=result['outcome']=='completed' and result['returncode']==0
        ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':result});checkpoint(out,ledger)
        print('COMPACT192 SPECIALIZED FOLLOWUP',name,ledger['rows'][-1]['status'],flush=True)
        if not ok:raise ArithmeticError('follow-up failed or censored; no retry')
    stage('freeze',sys.executable,'compact192_specialized_followup.py',['prepare'],p['prepare_seconds'])
    frozen=cert.read(LOCAL/'compact192-specialized-followup-v1/protocol.json');count=len(frozen['rows'])
    if not 1<=count<=p['maximum_roster'] or frozen['maximum_point_boxes']!=49*count:raise ArithmeticError('fixed whole eligible roster differs')
    ledger.update(roster=frozen['rows'],point_protocol_sha256=cert.hashed(LOCAL/'compact192-specialized-followup-v1/protocol.json'),maximum_point_boxes=49*count);checkpoint(out,ledger)
    stage('trial',sys.executable,'compact192_specialized_followup.py',['launch'],p['trial_seconds_per_curve']*count)
    stage('clouds',sys.executable,'audit_compact192_specialized_followup.py',[],p['audit_seconds_per_curve']*count)
    stage('geometry','/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python','replay_compact192_specialized_geometry.sage',[],p['geometry_seconds'])
    audit=cert.read(LOCAL/'compact192-specialized-followup-v1/verification-ledger.json')
    if audit['status']!='PASS' or [r['id'] for r in audit['rows']]!=[r['id'] for r in frozen['rows']]:raise ArithmeticError('complete full-cloud roster required')
    ledger['stronger_odd_prime_bounds']=[{'id':r['id'],'mod2':r['rank_lower_bound'],'odd':r['odd_modulus_lower_bounds']} for r in audit['rows'] if max(r['odd_modulus_lower_bounds'].values())>r['rank_lower_bound']]
    ledger['results']=[{k:r[k] for k in ('id','rank_lower_bound','odd_modulus_lower_bounds','attempted_charts','completed_boxes','retained_points')} for r in audit['rows']]
    ledger['status']='PASS_STRONGER_ODD_BOUNDS_REQUIRE_REVIEW' if ledger['stronger_odd_prime_bounds'] else 'PASS';checkpoint(out,ledger)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','launch']);stage=parser.parse_args().stage
    try:globals()[stage]()
    except Exception as error:
        if stage=='launch' and (D/'ledger.json').exists():
            d=cert.read(D/'ledger.json');d.update(status='FAILED_OR_CENSORED',reason=str(error));checkpoint(D/'ledger.json',d)
        raise
