#!/usr/bin/env python3
"""Independent replay: affine rational identities, pole intersections, modular RR rank.

Does not import the producer/evaluator or use Sage elliptic group operations.
All parent/evaluation reads occur after checking the17-arm terminal barrier.
"""
import hashlib,json,time,resource,sys
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector
ROOT=Path(__file__).resolve().parents[4];PKG=ROOT/'research/artifacts/generated-results/elliptic-curves/det1092_blind_mw16_v1'
def resolve(p):
 if p.exists():return p
 try:rel=p.relative_to(PKG/'runs')
 except ValueError:return p
 return ROOT/'research/artifacts/local/elliptic-curves/det1092_blind_mw16_v1/runs'/rel
def sha(p):return hashlib.sha256(resolve(p).read_bytes()).hexdigest()
COMPACT='--compact' in sys.argv
barrier=json.loads((PKG/'primary_complete.json').read_text());assert len(barrier['runs'])==17
if not COMPACT:
 for run in barrier['runs']:
  for p,h in run['files'].items():assert sha(PKG/p)==h
else:
 archive=json.loads((PKG/'raw_archive.json').read_text())
 assert archive['primary_barrier_sha256']==sha(PKG/'primary_complete.json')
 for run in archive['arms']:assert sha(PKG/run['compact_certificate'])==run['compact_sha256']
protocol=json.loads((PKG/'protocol.json').read_text());commit=protocol['starting_commit']
assert sha(PKG/'roster.json')==protocol['roster_sha256']
evaluation=json.loads((PKG/'evaluation.json').read_text());assert evaluation['primary_barrier_sha256']==sha(PKG/'primary_complete.json')
parent_rel='research/artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_v1.json'
pp=PKG/'frozen_sources'/parent_rel;assert sha(pp)==protocol['source_bindings'][parent_rel]
parent=json.loads(pp.read_text());R=PolynomialRing(QQ,'t');K=R.fraction_field()
def dec(x):return K(R(x['numerator']))/R(x['denominator'])
a1,a2,a3,a4,a6=list(map(dec,parent['a_invariants']))
def on(P):
 if P is None:return True
 x,y=P;return y*y+a1*x*y+a3*y==x**3+a2*x*x+a4*x+a6
def neg(P):return None if P is None else (P[0],-P[1]-a1*P[0]-a3)
def add(P,Q):
 if P is None:return Q
 if Q is None:return P
 x,y=P;u,v=Q
 if x==u:
  if Q==neg(P):return None
  assert P==Q
  m=(3*x*x+2*a2*x+a4-a1*y)/(2*y+a1*x+a3)
 else:m=(v-y)/(u-x)
 xx=m*m+a1*m-a2-x-u;yy=-(m+a1)*xx-(y-m*x)-a3
 return (xx,yy)
def mul(n,P):
 n=ZZ(n)
 if n<0:return mul(-n,neg(P))
 out=None
 while n:
  if n%2:out=add(out,P)
  n//=2
  if n:P=add(P,P)
 return out
def word(w,B):
 Q=None
 for n,P in zip(w,B):
  if n:Q=add(Q,mul(n,P))
 return Q
def h(P):
 if P is None:return QQ(0)
 x=P[0];finite=QQ(x.denominator().degree())/2
 infinity=max(QQ(0),QQ(x.numerator().degree()-x.denominator().degree())/2-2)
 return 4+2*(finite+infinity)
def pair(P,Q):
 if P is None or Q is None:return QQ(0)
 if P==Q:return h(P)
 x,y=P;u,v=neg(Q)
 if x==u:
  if (u,v)==neg(P):return (h(P)+h(Q))/2
  slope=(3*x*x+2*a2*x+a4-a1*y)/(2*y+a1*x+a3)
 else:slope=(v-y)/(u-x)
 xx=slope*slope+a1*slope-a2-x-u
 # The height only needs the abscissa: do not form an unused ordinate.
 return (h(P)+h(Q)-h((xx,K(0))))/2
B=[tuple(map(dec,p)) for p in parent['basis_weierstrass_coordinates']];assert all(map(on,B))
G=matrix(QQ,17,lambda i,j:pair(B[i],B[j]));assert G.det()==1092 and G==matrix(QQ,parent['generic_height_gram'])
b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6;b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
Delta=R(-b2*b2*b8-8*b4**3-27*b6**2+9*b2*b4*b6)
assert Delta.degree()==24 and Delta.gcd(Delta.derivative()).degree()==0
# A degree24 squarefree discriminant and coefficient weights give24 I1 and
# smooth infinity; chi=2, so no reducible-fibre corrections in the pole formula.
assert all(a.denominator().degree()==0 and a.numerator().degree()<=2*w for a,w in zip([a1,a2,a3,a4,a6],[1,2,3,4,6]))
records=[];total=time.monotonic()
for arm,ev in zip(json.loads((PKG/'roster.json').read_text())['arms'],evaluation['rows']):
 start=time.monotonic();aid=arm['arm_id'];assert ev['arm_id']==aid
 fixture=json.loads((PKG/arm['fixture']).read_text());assert sha(PKG/arm['fixture'])==arm['fixture_sha256']
 assert set(fixture)=={'starting_commit','arm_id','a_invariants','sections','gram','policy'}
 assert fixture['policy']==protocol['policy'] and fixture['starting_commit']==commit
 omitted=arm['omitted_basis_index_one_based']-1;keep=[j for j in range(17) if j!=omitted]
 points=[tuple(map(dec,p)) for p in fixture['sections']];assert points==[B[j] for j in keep]
 A=matrix(QQ,16,lambda i,j:pair(points[i],points[j]));assert A==matrix(QQ,fixture['gram']) and A.is_positive_definite() and A.det()!=0
 d=json.loads((PKG/'arms'/aid/'certificate.json').read_text()) if COMPACT else json.loads(resolve(PKG/'runs'/aid/'result.json').read_text())
 if COMPACT:assert d['original_result_sha256']==next(r for r in barrier['runs'] if r['arm_id']==aid)['files']['runs/'+aid+'/result.json']
 selection_path=(PKG/'arms'/aid/'selection.json') if COMPACT else resolve(PKG/'runs'/aid/'selection.json');selection=json.loads(selection_path.read_text()) if selection_path.exists() else None
 claims=[]
 for ci,c in enumerate(d.get('candidates',[])):
  Q=tuple(map(dec,c['coordinates']));assert on(Q)
  cross=vector(QQ,[pair(P,Q) for P in points]);schur=h(Q)-cross*A.solve_right(cross)
  rank=17 if schur>0 else 16
  assert schur>=0 and rank==c['rank'] and str(schur)==c['schur_complement']
  assert cross==vector(QQ,c['pairings']) and A.det()*schur==QQ(c['gram_determinant'])
  evaluated=ev['candidates'][ci];w=vector(ZZ,evaluated['full_basis_coordinates'])
  assert word(w,B)==Q and rank==(17 if w[omitted] else 16)
  claims.append({'candidate_id':ci,'rank':rank,'schur_complement':str(schur),'exact_affine_word_identity':True})
 centres={};members_checked=0;rr_rank_certificates=[];referenced=set()
 for stage in d.get('stages',[]):
  if not stage['status'].startswith('COMPLETE'):continue
  ci=stage['centre'];n=stage['n'];centre=selection['centres'][ci];w=vector(ZZ,centre['word'])
  assert w*A*w==centre['height'] and n in [(3*centre['height']+3)//4+1,(3*centre['height']+3)//4+2]
  if ci not in centres:centres[ci]=word(w,points)
  T=centres[ci];bounds=[n,n-4,n-6];flat=[]
  for member in stage['members']:
   fs=list(map(R,member['line_coefficients']));f0,f1,f2=fs
   assert all(f.degree()<=b for f,b in zip(fs,bounds))
   assert f0+f1*T[0]+f2*T[1]==0
   flat.append([f[j] for f,b in zip(fs,bounds) for j in range(b+1)])
   if not f2:
    assert member['status']=='VERTICAL_LINE_INHERITED' and not member['candidate_ids'];continue
   a=-f2*f2;b=f1*f1-a1*f1*f2-a2*f2*f2+T[0]*a
   c=2*f0*f1-a1*f0*f2-a3*f1*f2-a4*f2*f2+T[0]*b
   assert f0*f0-a3*f0*f2-a6*f2*f2+T[0]*c==0
   if COMPACT:
    def enc(v):
     v=K(v)
     return {'numerator':list(map(str,v.numerator().list())) or ['0'],'denominator':list(map(str,v.denominator().list()))}
    expanded={'residual_coefficients':list(map(enc,[c,b,a])),'discriminant':enc(b*b-4*a*c)}
    if 'square_root' in member['derived_fields']:expanded['square_root']=enc(K(b*b-4*a*c).sqrt())
    assert hashlib.sha256(json.dumps(expanded,sort_keys=True,separators=(',',':')).encode()).hexdigest()==member['derived_expansions_sha256']
    member.update(expanded)
   assert [K(c),K(b),K(a)]==list(map(dec,member['residual_coefficients']))
   disc=b*b-4*a*c;assert disc==dec(member['discriminant'])
   if member['status']=='NONSPLIT_OVER_Q_T':assert not disc.is_square() and not member['candidate_ids']
   else:
    assert member['status']=='SPLIT_OVER_Q_T';s=dec(member['square_root']);assert s*s==disc
    roots={}
    for sign in [1,-1]:
     x=(-b+sign*s)/(2*a);Q=(x,-(f0+f1*x)/f2);assert on(Q);roots[Q]=True
    got={tuple(map(dec,d['candidates'][j]['coordinates'])) for j in member['candidate_ids']}
    assert got==set(roots);referenced.update(member['candidate_ids'])
   members_checked+=1
  # Certify full RR kernel dimension by a finite-field lower bound on matrix
  # rank plus exact independent rational kernel rows, avoiding RREF replay.
  if stage['status']=='COMPLETE':
   assert len(flat)==stage['kernel_dimension']
   F=matrix(QQ,flat) if flat else matrix(QQ,0,sum(b+1 for b in bounds))
   assert F.rank()==len(flat)
   den=T[0].denominator().lcm(T[1].denominator())
   pol=[R(den*v)*R.gen()**j for v,b in zip([K(1),T[0],T[1]],bounds) for j in range(b+1)]
   M=matrix(QQ,max(p.degree() for p in pol)+1,len(pol),lambda i,j:pol[j][i])
   assert list(M.dimensions())==stage['matrix_shape'] and M*F.transpose()==0
   bound=M.ncols()-len(flat);prime=None;skipped=[]
   for p in [149,151,157,163,167,173,179,181,191,193]:
    try:r=matrix(GF(p),M).rank()
    except (ZeroDivisionError,TypeError,ValueError):skipped.append([p,'denominator']);continue
    if r==bound:prime=p;break
    skipped.append([p,'rank_drop',int(r)])
   assert prime is not None,('RR_RANK_REPLAY_UNRESOLVED',aid,ci,n)
   assert bound==stage['matrix_rank']
   rr_rank_certificates.append({'centre':ci,'n':n,'rank':bound,'prime':prime,'skips':skipped})
 assert referenced==set(range(len(d.get('candidates',[]))))
 assert bool([c for c in claims if c['rank']==17])==ev['success']
 record={'arm_id':aid,'status':'PASS','rank16_determinant':str(A.det()),'candidate_claims':claims,'RR_members_checked':members_checked,'RR_rank_certificates':rr_rank_certificates,'wall_seconds':time.monotonic()-start}
 records.append(record);print(aid,'PASS',len(claims),'candidate claims',round(record['wall_seconds'],3),flush=True)
out={'starting_commit':commit,'status':'PASS_INDEPENDENT_COMPACT_GENERIC_REPLAY' if COMPACT else 'PASS_INDEPENDENT_GENERIC_REPLAY','checker_sha256':sha(Path(__file__)),'protocol_sha256':sha(PKG/'protocol.json'),'primary_barrier_sha256':sha(PKG/'primary_complete.json'),'evaluation_sha256':sha(PKG/'evaluation.json'),'implementation':'Independent affine addition and tangent formulas; pole intersection heights; exact coefficient elimination; modular RR rank certificates plus rational kernel upper bounds. No producer/evaluator imports or elliptic group operations. Shared Sage exact rational arithmetic only.','arms':records,'successes':evaluation['successes'],'conclusion':evaluation['conclusion'],'resource':{'wall_seconds':time.monotonic()-total,'user_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime,'max_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}}
destination=PKG/('compact_replay.json' if COMPACT else 'replay.json')
if destination.exists():
 previous=json.loads(destination.read_text())
 for field in ['status','checker_sha256','protocol_sha256','primary_barrier_sha256','evaluation_sha256','successes','conclusion']:assert previous[field]==out[field]
 for before,after in zip(previous['arms'],out['arms']):
  assert {k:v for k,v in before.items() if k!='wall_seconds'}=={k:v for k,v in after.items() if k!='wall_seconds'}
else:
 with destination.open('x') as f:json.dump(out,f,indent=2,sort_keys=True);f.write('\n')
 destination.chmod(0o444)
print(out['status'],flush=True)
