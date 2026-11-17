"""NEW MK0 preregistration controls; never replay after registration."""
from pathlib import Path
import hashlib,json,time
import sympy as s
R=Path(__file__).resolve().parents[2];Hdir=Path(__file__).resolve().parent
receipt=Hdir/'receipt-pascal-energy-framing.json';cp=Hdir/'sources/pascal-energy-framing-context.json'
known={e['id']:e for e in json.loads((R/'MATH_STATUS.json').read_text())['entries']}
assert 'HC4MYMK0' not in known,'registered checker replay forbidden'
assert not receipt.exists(),'receipt exists; never replay'
c=json.loads(cp.read_text());start=time.perf_counter();checks=[];alg=0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,v,a=False):
 global alg
 assert bool(v),name
 checks.append(name);alg+=int(a)
def zero(M):return all(s.expand(q)==0 for q in M)
for p,d in c['retained_proof_and_source_sha256'].items():check('sealed '+p,sha(R/p)==d)
x,y,r,t=s.symbols('x y r t');z=s.Matrix([x,y]);w=s.Matrix([r,t]);coords=[x,y,r,t]
M,N=s.symbols('M N',integer=True,positive=True)
u=x+y**M;v=y+u**N
check('all-M-N target inverse recovers second',s.simplify(v-u**N)==y,True)
check('all-M-N target inverse recovers first',s.simplify(u-(v-u**N)**M)==x,True)
check('all-M-N target determinant',s.simplify(s.Matrix([u,v]).jacobian([x,y]).det())==1,True)
I=s.eye(2);O=s.zeros(2);J=O.row_join(I).col_join(I.row_join(O));U=I.row_join(s.I*I).col_join((I/2).row_join(-s.I*I/2))
check('Euclidean Hessian normalizing congruence',U.T*J*U==s.eye(4),True)
check('normalizing transpose identity',U.T==U.inv()*J,True)
# Symbolic rank-one product and degree recurrences, not an iterate campaign.
a1,a2,b1,b2,c1,c2=s.symbols('a1 a2 b1 b2 c1 c2');vv=s.Matrix([a1,a2]);rr=s.Matrix([[b1,b2]]);qq=s.Matrix([[c1,c2]])
check('noncommutative rank-one product identity',zero((vv*rr)*(vv*qq)-(rr*vv)[0]*vv*qq),True)
DD=s.symbols('D',integer=True,positive=True)
for n in [0,1,2,4]:
 check('all-D scalar energy degree sum '+str(n),s.expand(sum((DD-1)*DD**j for j in range(n+1))+DD**n-((DD+1)*DD**n-1))==0,True)
 check('all-D dual iterate degree sum '+str(n),s.expand(sum((DD-1)*DD**j for j in range(n))+1-DD**n)==0,True)
# Two small full energy compositions verify the order of matrix factors.
for mm,nn in [(2,2),(2,3)]:
 uu=x+y**mm;G=s.Matrix([uu,y+uu**nn]);A=G.jacobian([x,y]);F=s.Matrix(list(G)+list(A.T*w));energy=(w.T*A*z)[0];D=mm*nn
 check('first jet identity '+str((mm,nn)),G.subs({x:0,y:0})==s.zeros(2,1) and A.subs({x:0,y:0})==s.eye(2),True)
 check('energy initial degree '+str((mm,nn)),s.Poly(s.expand(energy),*coords).total_degree()==D+1,True)
 e1=s.expand(energy.subs(dict(zip(coords,F)),simultaneous=True))
 ordered=(w.T*A*A.subs({x:G[0],y:G[1]},simultaneous=True)*G)[0]
 check('full energy composition matrix order '+str((mm,nn)),s.expand(e1-ordered)==0,True)
 check('first energy iterate exact degree '+str((mm,nn)),s.Poly(e1,*coords).total_degree()==(D+1)*D,True)
 check('first energy difference exact degree '+str((mm,nn)),s.Poly(e1-energy,*coords).total_degree()==(D+1)*D,True)
 if (mm,nn)==(2,2):
  h=(w.T*G)[0];hc=s.Matrix([s.diff(h,vv0) for vv0 in coords]);check('normalized map equals swapped gradient',zero(F-J*hc),True)
  new=s.symbols('a b c d');XX=s.Matrix(new);sub=dict(zip(coords,U*XX));psi=s.expand(h.subs(sub,simultaneous=True));gp=s.Matrix([s.diff(psi,vv0) for vv0 in new]);conj=U.inv()*F.subs(sub,simultaneous=True)
  check('actual normalized gradient conjugacy',zero(gp-conj),True)
  check('actual normalized Hessian first jet',s.hessian(psi,new).subs(dict.fromkeys(new,0))==s.eye(4),True)
  epsi=sum(XX[j]*gp[j] for j in range(4))-psi
  check('energy transported under linear normalization',s.expand(epsi-energy.subs(sub,simultaneous=True))==0,True)
  # An actual polynomial inverse despite Pascal nontermination proved above.
  p0,q0=s.symbols('p0 q0');iy=q0-p0**nn;ix=p0-iy**mm
  check('actual base inverse composition',zero(G.subs({x:ix,y:iy},simultaneous=True)-s.Matrix([p0,q0])),True)
check('no new external import or primary cache',c['external_imports']==[] and c['new_primary_caches']==[])
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=alg,elapsed_seconds=time.perf_counter()-start,checks=checks,proof_sha256=sha(Hdir/'UNIVERSAL_KELLER_FRAMING_AND_PASCAL_ENERGY_NONTERMINATION.md'),checker_sha256=sha(Path(__file__)),context_sha256=sha(cp),registered_mathematical_checkers_replayed=[],independent_replay=False,formal_verification=False,external_review=False)
receipt.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
