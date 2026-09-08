#!/usr/bin/env sage-python
"""Check portable finite-rank and exact birational certificates without raw runs."""
import gzip,hashlib,json,tempfile,time
from pathlib import Path
from fractions import Fraction as F
from sage.all import QQ,ZZ,matrix,identity_matrix,EllipticCurve,PolynomialRing
import audit_recorded_point_mod2_rank_v3 as mod2
import audit_retained_cloud_modl as modl
import pari_pointed_backend as backend
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
from research_runtime.search_state import raw_state
ROOT=Path(__file__).resolve().parents[2];PKG=ROOT/'artifacts/generated-results/elliptic-curves/det1092_basis_cascade_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
start=time.monotonic();summary=read(PKG/'summary.json');results=[];equations={}
for name in ['presentation_ledger','point_novelty_ledger','raw_archive']:
 assert sha(PKG/(name+'.json'))==summary[name+'_sha256']
# Universal inverse chord map, checked over a rational function field.
base=PolynomialRing(QQ,names=['s','x0','y0','A']);s,x0,y0,A=base.gens();R=PolynomialRing(base.fraction_field(),'w');w=R.gen()
B=y0**2-x0**3-A*x0;q=s**4-6*x0*s**2-8*y0*s-3*x0**2-4*A
X=(s**2-x0+w)/2;Y=-y0+s*(X-x0)
assert (Y**2-X**3-A*X-B).quo_rem(w**2-q)[1]==0
assert (Y+y0)/(X-x0)==s and 2*X-s*s+x0==w
S=PolynomialRing(base.fraction_field(),'z');z=S.gen()
assert (z**4-6*x0*z**2-8*y0*z-3*x0**2-4*A).discriminant()==256*(-16*(4*A**3+27*B**2))
scratch=ROOT/'artifacts/local/elliptic-curves/det1092-basis-cascade-v1-compact-replay';scratch.mkdir(exist_ok=True)
for row in summary['arms']:
 aid=row['arm_id'];p=PKG/f'evidence-{aid}.json.gz';assert sha(p)==row['evidence_sha256'];d=json.loads(gzip.decompress(p.read_bytes()))
 fixture=d['fixture'];model=tuple(map(F,fixture['curve']));E=EllipticCurve(QQ,model);assert E.discriminant()!=0
 U=matrix(ZZ,fixture['basis_transform']);assert abs(U.det())==1 and U*matrix(ZZ,fixture['inverse_transform'])==identity_matrix(ZZ,17)
 reference=[E(*p) for p in fixture['reference_points']];assert [sum((int(c)*p for c,p in zip(r,reference)),E(0)) for r in U.rows()]==[E(*p) for p in fixture['points']]
 cache=QuotientOnlyReductionCache(MemoryFactStore());state=raw_state(model,tuple(tuple(map(F,p)) for p in fixture['points']),cache=cache,prime_bound=1000);assert state.rank==17
 charts=points=0
 with tempfile.TemporaryDirectory(dir=scratch) as temp:
  temp=Path(temp)
  for epoch in d['stages']:
   stage=epoch['stage'];assert stage is not None
   assert epoch['final_state']==state.record()
   cloud=list(state.basis);seen={(x,abs(y)) for x,y in cloud}
   for chart in epoch['charts']:
    search=PointedQuarticSearch(state=state,centre={'coefficients':chart['centre']['representative']},coordinate_policy=chart['mapping']['coordinate_policy'])
    found=backend.replay(search,chart['mapping'],chart['search']);charts+=1;points+=len(found)
    coeff={k:chart['mapping'][k] for k in ['reduced_P','reduced_Q']}
    eqid=hashlib.sha256(json.dumps(coeff,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    equations[eqid]=coeff
    # Audit collection is transcript order, not set order.
    for raw in chart['search']['finite_curve_points']:
     p=(F(raw['x']),F(raw['y']));key=(p[0],abs(p[1]))
     if key not in seen:seen.add(key);cloud.append(p)
   audit=epoch['audit'];assert audit['points']==[list(map(str,p)) for p in cloud]
   ap=temp/'audit.json';ap.write_text(json.dumps(audit));mod2.check(ap)
   if stage['after']>stage['before']:
    ml=epoch['modl'];assert ml['points']==audit['points'];mp=temp/'modl.json';mp.write_text(json.dumps(ml));modl.check(mp)
    assert audit['independent_points'][:state.rank]==[list(map(str,p)) for p in state.basis]
    state=raw_state(model,tuple(tuple(map(F,p)) for p in audit['independent_points']),cache=cache,prime_bound=1000)
   assert stage['after']==state.rank
   print('PORTABLE STAGE',aid,stage['epoch'],state.rank,flush=True)
 assert state.rank==row['final_rank_lower_bound'];results.append({'arm_id':aid,'rank_lower_bound':state.rank,'charts_replayed':charts,'point_records_replayed':points,'status':'PASS'})
ledger=read(PKG/'presentation_ledger.json')
assert {r['equation_id']:r['coefficients_ascending'] for r in ledger['equations']}==equations
assert ledger['new_elliptic_curve_isomorphism_classes']==0
out={'starting_commit':summary['starting_commit'],'status':'PASS_PORTABLE_EXACT_REPLAY','arms':results,'distinct_tested_equations_bound_to_exact_maps':len(equations),'universal_chord_birational_identity':{'quartic':'w^2=s^4-6*x0*s^2-8*y0*s-3*x0^2-4*A','B':'y0^2-x0^3-A*x0','X':'(s^2-x0+w)/2','Y':'-y0+s*(X-x0)','inverse_s':'(Y+y0)/(X-x0)','inverse_w':'2*X-s^2+x0','polynomial_remainder_zero':True,'inverse_identities_verified':True,'quartic_discriminant_equals_256_times_curve_discriminant':True},'checker_sha256':sha(Path(__file__)),'summary_sha256':sha(PKG/'summary.json'),'wall_seconds':time.monotonic()-start}
with (PKG/'compact_replay.json').open('x') as f:json.dump(out,f,indent=2,sort_keys=True);f.write('\n')
print('PASS PORTABLE',results,flush=True)
