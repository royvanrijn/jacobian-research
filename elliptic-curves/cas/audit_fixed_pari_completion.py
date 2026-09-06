#!/usr/bin/env python3
"""Combine disjoint adaptive chart attempts and audit all retained curve points."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from audit_recorded_point_mod2_rank_v2 import build,check
from continue_fixed_pari_search import paths
ROOT=Path(__file__).resolve().parents[2];LOCAL=ROOT/'artifacts/local/elliptic-curves';ART=ROOT/'artifacts/generated-results/elliptic-curves'
BASE={'rank26':['new-rank26-fixed-tails-v1/combined-point-cloud-only.json','new-rank26-pari43-v1/candidate-00/result.json'],
      'small-conductor':['prospective-mw16-f5-next12-v1/a1-fibration-05/candidate-00/result.json','prospective-mw16-small-conductor-followup-v2/candidate-00/result.json']}

def main(case):
    d,parent,maps,rank,_,_=paths(case);p=cert.read(d/'protocol.json');seed_path=ROOT/p['seed_path'];path=d/'candidate-00/result.json';seed,data=cert.read(seed_path),cert.read(path)
    if cert.hashed(seed_path)!=p['seed_sha256'] or cert.hashed(path)!=cert.read(d/'terminal.json')['result_sha256']:raise ArithmeticError('terminal inputs changed')
    stages=cert.read(parent/'verification/result.json')['rows'];child=cert.read(d/'verification/charts-result.json')
    if len(stages)!=3 or any(r['status']!='PASS' for r in stages) or child['outcome']!='completed' or child['returncode']!=0:raise ArithmeticError('source histories have not replayed')
    indexed=[(i,r) for i,r in enumerate(seed['charts'])]+[(r['parent_chart'],r) for r in data['charts']]
    if [i for i,r in indexed]!=list(range(len(indexed))) or len(seed['charts'])!=p['start_chart']:raise ArithmeticError('overlapping or noncontiguous adaptive roster')
    geometry=cert.read(maps);coverage=[]
    for i,r in indexed:
        if r['centre']!=geometry['rows'][i]['centre'] or r['search']['height_bound']!=100000:raise ArithmeticError('fixed adaptive map roster differs')
        coverage.append({'chart':i,'status':r['search']['status'],'height':100000,'infinity_checked':r['search']['infinity_checked']})
    cov=ART/(case.replace('-','_')+'_fast301_completion_v1.json')
    if cov.exists():raise FileExistsError('preserve completion proof')
    checkpoint(cov,{'schema':'elliptic-curves.fixed-pari-adaptive-completion.v1','status':'PASS','case':case,'planned_adaptive_maps':301,'attempted_maps':len(indexed),'complete_finite_boxes':sum(r['status']=='bounded_search_complete' for r in coverage),'coverage':coverage,'sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in (Path(__file__).resolve(),seed_path,path,maps,parent/'verification/result.json',d/'verification/charts-result.json')},'claim_boundary':'Disjoint fixed adaptive map roster and replayed retained histories. Complete finite coverage trusts pinned PARI invocations; neither finite completion nor a miss supplies a rank upper bound.'})
    inputs=[(LOCAL/name,cert.read(LOCAL/name)) for name in BASE[case]]+[(seed_path,seed),(path,data)]
    if any(r['curve']!=data['curve'] or r['generic_points']!=data['generic_points'] for _,r in inputs):raise ArithmeticError('combined curve transports differ')
    combined=d/'all-retained-point-cloud-only.json'
    checkpoint(combined,{'schema':'elliptic-curves.combined-retained-point-cloud.v1','status':'COMBINED_POINT_CLOUD_ONLY','family':data['family'],'parameter':data['parameter'],'curve':data['curve'],'generic_points':data['generic_points'],'rank_lower_bound':data['rank_lower_bound'],'final_state':data['final_state'],'charts':[r for _,s in inputs for r in s['charts']],'source_inputs':[{'path':str(q.relative_to(ROOT)),'sha256':cert.hashed(q),'charts':len(r['charts'])} for q,r in inputs],'claim_boundary':'All listed historical and current raw point clouds, including points not kept by orbit compression. This concatenation is not a chronological state transcript; each source has its own exact replay.'})
    cloud=ART/(case.replace('-','_')+'_all_current_retained_mod2_v1.json');build(combined,cloud,997,cert.hashed(combined));check(cloud)
    print('COMPLETED ADAPTIVE ROSTER',case,len(indexed),'boxes;',len(cert.read(cloud)['points']),'retained points; rank >=',cert.read(cloud)['rank_lower_bound'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',choices=BASE,required=True);main(p.parse_args().case)
