#!/usr/bin/env python3
"""Independent rational geometry and full-cloud proofs before control admission."""
import argparse,sys
from pathlib import Path
import larger_specialized_parity_trial as trial
import certify_compact_r17_candidates as cert
from replay_retention24_geometry import geometry,cloud_check,tuples
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=trial.ROOT;ART=trial.ART;CAS=trial.CAS
def main(index,check):
    trial.configure(index);p=trial.protocol();folder=trial.D;row=trial.ROW;path=folder/'result.json';raw=cert.read(path);maps=cert.read(folder/'maps.json')
    terminal=cert.read(folder/'replay.supervisor.json')
    if terminal['outcome']!='completed' or terminal['returncode']!=0:raise ArithmeticError('terminal exact history required')
    geometry(raw,maps,tuples(raw['initial_state']['state']['reductions']['points']),{**p,'maps_path':folder/'maps.json'})
    stem='larger_parity_'+row['id'].replace('-','_');mod2=ART/(stem+'_mod2_v1.json');odd=ART/(stem+'_modl_v1.json');out=ART/(stem+'_coverage_v1.json')
    if not check:
        if out.exists():raise FileExistsError('preserve larger parity cloud proof')
        for name,script,args in [('mod2-build','audit_recorded_point_mod2_rank_v3.py',['--input',str(path),'--input-sha256',cert.hashed(path),'--output',str(mod2),'--prime-bound','997']),('mod2-check','audit_recorded_point_mod2_rank_v3.py',['--check',str(mod2)]),('odd-build','audit_retained_cloud_modl.py',['--input',str(mod2),'--output',str(odd)]),('odd-check','audit_retained_cloud_modl.py',['--check',str(odd)])]:
            s=run([sys.executable,str(CAS/script),*args],limits=Limits(300,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=ROOT)
            if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('independent full-cloud proof failed')
    cloud_check(raw,raw['charts'],mod2,path);cloud=cert.read(mod2);other=cert.read(odd)
    result={'schema':'elliptic-curves.larger-parity-point-coverage.v1','status':'PASS','id':row['id'],
        'sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in (Path(__file__).resolve(),trial.BATCH/'protocol.json',path,folder/'maps.json',folder/'replay.supervisor.json',mod2,odd)},
        'attempted_charts':len(raw['charts']),'completed_boxes':sum(c['search']['status']=='bounded_search_complete' for c in raw['charts']),
        'retained_points':len(cloud['points']),'mod2_lower_bound':cloud['rank_lower_bound'],
        'odd_modulus_lower_bounds':{str(a['modulus']):a['finite_column_rank'] for a in other['audits']},
        'claim_boundary':p['boundaries']}
    if check:
        if cert.read(out)!=result:raise ArithmeticError('larger parity coverage differs')
    else:checkpoint(out,result)
    print('LARGER PARITY POINT AUDIT',row['id'],result['mod2_lower_bound'],result['odd_modulus_lower_bounds'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,required=True);p.add_argument('--check',action='store_true');a=p.parse_args();main(a.index,a.check)
