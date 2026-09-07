#!/usr/bin/env python3
"""Standalone complete-cloud proofs for each terminal fixed height control."""
import argparse,sys
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';LOCAL=ROOT/'artifacts/local/elliptic-curves';ART=ROOT/'artifacts/generated-results/elliptic-curves'
FOLDERS={100000:LOCAL/'native11952-height-pair-v1/100000',1000000:LOCAL/'native11952-height-pair-v1/1000000',125000:LOCAL/'native11952-height125-control-v1/125000'}
def verify(height):
    folder=FOLDERS[height];out=folder/'verification.json'
    if out.exists():raise FileExistsError('preserve height verification')
    path=folder/'result.json';d=cert.read(path)
    for name in ('worker','replay'):
        s=cert.read(folder/(name+'.supervisor.json'))
        if s['outcome']!='completed' or s['returncode']!=0 or cert.hashed(Path(s['log']))!=s['log_sha256']:raise ArithmeticError('terminal worker and exact history replay required')
    if d['height']!=height or d['status']!='COMPLETE_DECLARED_HEIGHT_ARM' or len(d['charts'])!=49:raise ArithmeticError('fixed49 attempts required')
    cloud=ART/f'native11952_height{height}_mod2_v1.json';rows=[]
    for label,args in [('cloud-build',['--input',str(path),'--input-sha256',cert.hashed(path),'--output',str(cloud),'--prime-bound','997']),('cloud-check',['--check',str(cloud)])]:
        s=run([sys.executable,str(CAS/'audit_recorded_point_mod2_rank_v3.py'),*args],limits=Limits(120,1610612736),log_path=folder/(label+'.log'),checkpoint_path=folder/(label+'.supervisor.json'),cwd=ROOT);rows.append({'name':label,'status':'PASS' if s['outcome']=='completed' and s['returncode']==0 else 'FAILED_OR_CENSORED','supervision':s});checkpoint(out,{'status':'RUNNING','rows':rows})
        if rows[-1]['status']!='PASS':raise ArithmeticError('height cloud stage failed')
    proof=cert.read(cloud);checkpoint(out,{'status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),CAS/'audit_recorded_point_mod2_rank_v3.py',path,cloud,folder/'replay.supervisor.json')},'height':height,'attempted_charts':49,'completed_boxes':sum(c['search']['status']=='bounded_search_complete' for c in d['charts']),'rank_lower_bound':proof['rank_lower_bound'],'point_count':len(proof['points']),'certificate':str(cloud.relative_to(ROOT)),'rows':rows,'claim_boundary':'Exactly certified points from a known-control attempt; censored boxes remain censored. No new curve, exact rank or prospective guarantee.'});print('VERIFIED HEIGHT CONTROL',height,proof['rank_lower_bound'],len(proof['points']),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--height',type=int,required=True,choices=sorted(FOLDERS));a=p.parse_args();verify(a.height)
