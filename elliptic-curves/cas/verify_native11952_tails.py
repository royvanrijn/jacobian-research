#!/usr/bin/env python3
"""Replay the tails and independently audit every retained point on this curve."""
from pathlib import Path
import sys
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas'
D=ROOT/'artifacts/local/elliptic-curves/native11952-metric49-tails-v1'
ART=ROOT/'artifacts/generated-results/elliptic-curves'


def verify():
    terminal=cert.read(D/'terminal.json');path=D/'candidate-00/result.json';data=cert.read(path)
    if cert.hashed(path)!=terminal['result_sha256']:raise ArithmeticError('terminal tail input changed')
    folder=D/'verification'
    if (folder/'protocol.json').exists():raise FileExistsError('preserve tail verification')
    coverage=ART/'native11952_tail_coverage_v1.json';cloud=ART/'native11952_all_retained_mod2_v1.json'
    checkpoint(folder/'protocol.json',{'tail_input_sha256':cert.hashed(path),
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),CAS/'replay_native11952_tails.py',CAS/'audit_recorded_point_mod2_rank_v2.py')},
        'chart_replay_wall_seconds':400,'cloud_build_wall_seconds':300,'cloud_check_wall_seconds':180,'rss_bytes':1610612736,
        'scope':'Replay all retained tail charts and contiguous coverage. Separately combine raw points from the original49 and retained tail charts for one exact finite-quotient audit. The combined point-only input is not a synthetic state-history replay.'})
    result={'status':'RUNNING','rows':[]};checkpoint(folder/'result.json',result)
    def stage(name,args,seconds):
        r=run([sys.executable,*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=ROOT)
        ok=r['outcome']=='completed' and r['returncode']==0
        result['rows'].append({'stage':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'result.json',result);print('TAIL VERIFICATION',name,result['rows'][-1]['status'],flush=True)
        return ok
    stage('charts',[str(CAS/'replay_native11952_tails.py'),'--input',str(path),'--output',str(coverage)],400)
    protocol=cert.read(D/'protocol.json');parent_path=ROOT/protocol['input_path'];parent=cert.read(parent_path)
    if cert.hashed(parent_path)!=protocol['input_sha256']:raise ArithmeticError('point-cloud parent changed')
    inputs=[(parent_path,parent),(path,data)]
    if any(r['curve']!=data['curve'] or r['generic_points']!=data['generic_points'] for _,r in inputs):raise ArithmeticError('point-cloud models differ')
    combined=D/'combined-point-cloud-only.json'
    checkpoint(combined,{'schema':'elliptic-curves.combined-retained-point-cloud.v1','status':'COMBINED_POINT_CLOUD_ONLY',
        'family':data['family'],'parameter':data['parameter'],'curve':data['curve'],'generic_points':data['generic_points'],
        'rank_lower_bound':data['rank_lower_bound'],'final_state':data['final_state'],
        'charts':[row for _,r in inputs for row in r['charts']],
        'source_inputs':[{'path':str(p.relative_to(ROOT)),'sha256':cert.hashed(p),'charts':len(r['charts'])} for p,r in inputs],
        'claim_boundary':'Concatenation solely for exact rational-point collection and independence. It is not a single chronological MWState transcript. Each original state history has a separate replay.'})
    if stage('cloud-build',[str(CAS/'audit_recorded_point_mod2_rank_v2.py'),'--input',str(combined),'--input-sha256',cert.hashed(combined),'--output',str(cloud)],300):
        stage('cloud-check',[str(CAS/'audit_recorded_point_mod2_rank_v2.py'),'--check',str(cloud)],180)
    result['status']='COMPLETE_DECLARED_REPLAY_ATTEMPTS';checkpoint(folder/'result.json',result)


if __name__=='__main__':verify()
