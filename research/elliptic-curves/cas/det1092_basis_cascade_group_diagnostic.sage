#!/usr/bin/env sage-python
"""Post-run exact subgroup comparison; no feedback into frozen cascades."""
import gzip,hashlib,json,time
from pathlib import Path
from fractions import Fraction as F
from sage.all import QQ,ZZ,RR,RealField,matrix,vector,EllipticCurve,pari
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
from research_runtime.search_state import raw_state
import memory_rank_certificate as memory
from mod2_reduction_independence import _primes_up_to
ROOT=Path(__file__).resolve().parents[2];PKG=ROOT/'artifacts/generated-results/elliptic-curves/det1092_basis_cascade_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
summary=read(PKG/'summary.json');assert all(r['status']=='TERMINAL' and r['replay_status']=='PASS_EXACT_REPLAY' for r in summary['arms'])
d={r['arm_id']:json.loads(gzip.decompress((PKG/f"evidence-{r['arm_id']}.json.gz").read_bytes())) for r in summary['arms']};base=d['original'];model=tuple(map(F,base['fixture']['curve']));E=EllipticCurve(QQ,model)
cache=QuotientOnlyReductionCache(MemoryFactStore());primes=[p for p in _primes_up_to(1000) if p!=2 and E.discriminant()%p!=0];reference=[tuple(map(F,q)) for q in base['fixture']['points']];unions=[];start=time.monotonic()
for rank in range(18,1+min(r['final_rank_lower_bound'] for r in summary['arms'])):
 points=list(reference);seen={(x,abs(y)) for x,y in points};arms={}
 for aid,b in d.items():
  ep=next(e for e in b['stages'] if e['stage']['after']==rank)
  extras=[tuple(map(F,q)) for q in ep['audit']['independent_points'][17:]];arms[aid]=[list(map(str,q)) for q in extras]
  for q in extras:
   key=q[0],abs(q[1])
   if key not in seen:seen.add(key);points.append(q)
 state=raw_state(model,tuple(points),cache=cache,prime_bound=1000);bound=state.rank
 torsion=state.record()['state']['no_two_torsion_prime'];proof=memory.checked_rank(model,state.basis,primes,torsion)
 # Independent re-evaluation of the compact exact signatures with their own primes.
 replay=memory.checked_rank(model,state.basis,[r['prime'] for r in proof['signatures']],proof['no_rational_2_torsion_prime']);assert proof==replay
 unions.append({'per_arm_rank':rank,'union_rank_lower_bound':bound,'distinct_input_points_modulo_sign':len(points),'independent_points':[list(map(str,p)) for p in state.basis],'rank_certificate':proof,'arm_extension_points':arms})
 print('POST RUN UNION',rank,bound,flush=True)
# Numerical heights propose rational words only; exact rational group identities decide.
basepoints=[E(*q) for q in base['stages'][-1]['audit']['independent_points']];n=len(basepoints);unknown=[]
for aid,b in d.items():
 for p in b['stages'][-1]['audit']['independent_points'][17:]:
  q=E(*p)
  if q not in basepoints and q not in unknown:unknown.append(q)
relations={};RF=RealField(384);pari.default('realprecision',120)
if unknown:
 points=basepoints+unknown;H=pari(E).ellheightmatrix([list(q.xy()) for q in points],precision=384);G=matrix(RF,n,n,[RF(str(H[i,j])) for i in range(n) for j in range(n)])
 for k,q in enumerate(unknown):
  v=G.solve_right(vector(RF,[RF(str(H[i,n+k])) for i in range(n)]));w=[QQ(str(F(str(x)).limit_denominator(1000000))) for x in v];key=str(q.xy())
  if any(abs(x-RF(y))>RF(2)**-160 for x,y in zip(v,w)):
   relations[key]={'status':'NO_EXACT_WORD_FOUND_IN_DECLARED_DENOMINATOR_BOUND'};continue
  denominator=ZZ(1)
  for x in w:denominator=denominator.lcm(x.denominator())
  check=sum((int(denominator*c)*p for c,p in zip(w,basepoints)),E(0))
  if check!=int(denominator)*q:
   relations[key]={'status':'PROPOSED_WORD_REJECTED_BY_EXACT_GROUP_LAW','proposed_word':list(map(str,w)),'clearing_denominator':int(denominator)};continue
  relations[key]={'status':'EXACT_RATIONAL_GROUP_IDENTITY','word':list(map(str,w)),'clearing_denominator':int(denominator)}
comparisons=[]
for aid,b in d.items():
 words=[];unknowns=[];U=matrix(ZZ,b['fixture']['basis_transform'])
 for r in U.rows():words.append(list(map(QQ,r))+[QQ(0)]*(n-17))
 for raw in b['stages'][-1]['audit']['independent_points'][17:]:
  q=E(*raw)
  if q in basepoints:
   w=[QQ(0)]*n;w[basepoints.index(q)]=QQ(1);words.append(w)
  elif relations[str(q.xy())]['status']=='EXACT_RATIONAL_GROUP_IDENTITY':words.append(list(map(QQ,relations[str(q.xy())]['word'])))
  else:unknowns.append(raw)
 if unknowns:
  comparisons.append({'arm_id':aid,'status':'UNRESOLVED_EXACT_RELATION','unresolved_points':unknowns});continue
 T=matrix(QQ,words);det=T.det() if T.nrows()==n else None
 status='EXACT_INCLUSION_IN_ORIGINAL_FINAL_RATIONAL_SPAN';index=None
 if det is not None and det:
  if T.denominator()==1:
   status='SAME_INTEGRAL_SUBGROUP' if abs(det)==1 else 'PROPER_INTEGRAL_SUBGROUP_OF_ORIGINAL';index=int(abs(det))
  elif T.inverse().denominator()==1:
   inverse=T.inverse();newpoints=[E(*q) for q in b['stages'][-1]['audit']['independent_points']]
   assert all(sum((int(c)*q for c,q in zip(r,newpoints)),E(0))==p for r,p in zip(inverse.rows(),basepoints))
   status='PROPER_INTEGRAL_OVERGROUP_OF_ORIGINAL';index=int(abs(1/det))
  else:status='COMMENSURABLE_FINAL_SUBGROUPS_WITH_EXACT_RATIONAL_WORDS'
 comparisons.append({'arm_id':aid,'status':status,'finite_index_when_included':index,'words_in_original_final_basis':list(map(lambda r:list(map(str,r)),T.rows())),'matrix_rank':int(T.rank()),'determinant':str(det) if det is not None else None})
 print('FINAL GROUP',aid,comparisons[-1]['status'],flush=True)
result={'starting_commit':summary['starting_commit'],'status':'COMPLETE_POST_RUN_EXACT_DIAGNOSTIC','intermediate_unions':unions,'final_subgroup_comparisons':comparisons,'proposed_word_relations':relations,'word_search_policy':'384-bit canonical height linear solve, denominator<=1000000, residual tolerance2^-160 for proposal only; exact multiple group identity required. No failed proposal is a nonmembership claim.','maximum_union_rank_lower_bound':max(r['union_rank_lower_bound'] for r in unions),'above31_independent_review_required':any(r['union_rank_lower_bound']>31 for r in unions),'wall_seconds':time.monotonic()-start,'checker_sha256':sha(Path(__file__)),'summary_sha256':sha(PKG/'summary.json'),'boundary':'Post-run unions and exact words never feed any primary worker. A union gain is not the cost of a prospectively selected hybrid search.'}
with (PKG/'group_diagnostic.json').open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
