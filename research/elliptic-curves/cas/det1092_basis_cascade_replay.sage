#!/usr/bin/env sage-python
"""Independent exact replay of a terminal four-basis cascade; no point search."""
import argparse,csv,hashlib,json,time
from pathlib import Path
from fractions import Fraction as F
from sage.all import ZZ,QQ,matrix,EllipticCurve,block_diagonal_matrix,identity_matrix
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
from research_runtime.search_state import raw_state
from pointed_quartic_search import PointedQuarticSearch
import pari_pointed_backend as backend
import audit_recorded_point_mod2_rank_v3 as mod2
import audit_retained_cloud_modl as modl
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'artifacts/local/elliptic-curves/det1092-basis-cascade-v1';PKG=ROOT/'artifacts/generated-results/elliptic-curves/det1092_basis_cascade_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):
 with p.open('x') as f:json.dump(d,f,indent=2,sort_keys=True);f.write('\n')
def replay(aid):
 started=time.monotonic();arm=D/aid;fixture=read(arm/'fixture.json');policy=read(arm/'protocol.json');roster=read(PKG/'roster.json');row=next(r for r in roster['arms'] if r['arm_id']==aid)
 assert sha(arm/'fixture.json')==row['fixture_sha256'] and sha(arm/'orbits.tsv')==row['orbits_sha256'] and sha(arm/'protocol.json')==row['protocol_sha256']
 assert all(sha(ROOT/k)==v for k,v in policy['sources'].items())
 supervisor=read(arm/'supervisor-result.json');assert supervisor['exit_code']==0 and supervisor['terminal_file_exists']
 terminal=read(arm/'run/terminal.json');U=matrix(ZZ,fixture['basis_transform']);inv=matrix(ZZ,fixture['inverse_transform']);assert abs(U.det())==1 and U*inv==identity_matrix(ZZ,17)
 model=tuple(map(F,fixture['curve']));curve=EllipticCurve(QQ,model);reference=[curve(*q) for q in fixture['reference_points']];basis=[tuple(map(F,q)) for q in fixture['points']]
 assert [sum((int(c)*p for c,p in zip(row,reference)),curve(0)) for row in U.rows()]==[curve(*q) for q in basis]
 original=ROOT/'artifacts/generated-results/elliptic-curves/curve302_parent_degree2_multisection_orbits_v1.tsv';orbits_checked=0
 with original.open() as left,(arm/'orbits.tsv').open() as right:
  a=csv.DictReader(left,delimiter='\t');b=csv.DictReader(right,delimiter='\t')
  for x,y in zip(a,b,strict=True):
   assert matrix(ZZ,1,17,list(map(int,y['parent_MW17_w'].split())))*U==matrix(ZZ,1,17,list(map(int,x['parent_MW17_w'].split())))
   assert all(x[k]==y[k] for k in x if k!='parent_MW17_w');orbits_checked+=1
 cache=QuotientOnlyReductionCache(MemoryFactStore());state=raw_state(model,tuple(basis),cache=cache,prime_bound=1000);assert state.rank==17
 chart_count=point_count=0;stages=[];tested=set();ledger=[]
 for stage in terminal['stages']:
  wd=arm/'run'/f"epoch-{stage['epoch']:02d}";selection=read(wd/'selection.json');transport=read(wd/'metric-transport.json');T=block_diagonal_matrix(U,identity_matrix(ZZ,state.rank-17))
  assert matrix(ZZ,transport['transport'])==T
  assert T*matrix(ZZ,transport['reference_gram'])*T.transpose()==matrix(ZZ,selection['rounded_gram'])
  assert selection['basis']==[list(map(str,q)) for q in state.basis]
  assert transport['reference_basis']==fixture['reference_points']+selection['basis'][17:]
  charts=[]
  for j in range(stage['charts']):
   path=wd/f'chart-{j:03d}.json';chart=read(path);assert chart['index']==j and chart['centre']==selection['centres'][j]
   key=tuple(chart['centre']['point']);assert key not in tested;tested.add(key)
   search=PointedQuarticSearch(state=state,centre={'coefficients':chart['centre']['representative']},coordinate_policy=chart['mapping']['coordinate_policy'])
   points=backend.replay(search,chart['mapping'],chart['search']);point_count+=len(points);chart_count+=1;charts.append(chart)
   ledger.append({'path':str(path.relative_to(ROOT)),'sha256':sha(path),'epoch':stage['epoch'],'index':j,'centre':list(key),'search_status':chart['search']['status'],'points':len(points)})
  snapshot=wd/f"cloud-{stage['charts']-1:03d}.json";snap=read(snapshot);auditpath=wd/stage['audit'];audit=read(auditpath)
  assert snap['charts']==charts and snap['final_state']==state.record()
  assert audit['input_sha256']==sha(snapshot) and sha(auditpath)==stage['audit_sha256']
  expected=[tuple(map(F,p)) for p in snap['final_state']['state']['reductions']['points']];seen={(x,abs(y)) for x,y in expected}
  for chart in charts:
   for raw in chart['search']['finite_curve_points']:
    q=(F(raw['x']),F(raw['y']));key=q[0],abs(q[1])
    if key not in seen:seen.add(key);expected.append(q)
  assert audit['points']==[list(map(str,p)) for p in expected]
  mod2.check(auditpath)
  if stage['after']>stage['before']:
   more=read(wd/'modl.json');assert more['input_sha256']==sha(auditpath) and more['points']==audit['points'];modl.check(wd/'modl.json')
   assert audit['independent_points'][:state.rank]==[list(map(str,p)) for p in state.basis]
   state=raw_state(model,tuple(tuple(map(F,q)) for q in audit['independent_points']),cache=cache,prime_bound=1000)
  assert state.rank==stage['after'];stages.append({'epoch':stage['epoch'],'before':stage['before'],'after':stage['after'],'charts':stage['charts'],'mod2_sha256':sha(auditpath),'modl_sha256':sha(wd/'modl.json') if (wd/'modl.json').exists() else None})
  print('REPLAY STAGE',aid,stage['epoch'],state.rank,flush=True)
 assert state.rank==terminal['final_rank_lower_bound'] and chart_count==terminal.get('charts',chart_count)
 result={'starting_commit':roster['starting_commit'],'arm_id':aid,'status':'PASS_EXACT_REPLAY','final_rank_lower_bound':state.rank,'orbits_transport_checked':orbits_checked,'charts_replayed':chart_count,'point_records_replayed':point_count,'stages':stages,'chart_bindings':ledger,'wall_seconds':time.monotonic()-started,'checker_sha256':sha(Path(__file__)),'terminal_sha256':sha(arm/'run/terminal.json'),'claim_boundary':'Exact rational map and point witnesses, finite group lower-bound certificates with torsion exclusion, and unimodular seed/orbit/metric transport. No upper bound, exhaustive CVP re-enumeration, or new point search in replay.'}
 write(PKG/f'replay-{aid}.json',result);print('PASS',aid,state.rank,chart_count,flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('arm');replay(ap.parse_args().arm)
