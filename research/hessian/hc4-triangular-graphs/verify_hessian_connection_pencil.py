"""NEW MI0 preregistration controls; refuse registered replay."""
from pathlib import Path
import hashlib,json,time
import sympy as s
R=Path(__file__).resolve().parents[2];Hdir=Path(__file__).resolve().parent
receipt=Hdir/'receipt-hessian-connection-pencil.json';cp=Hdir/'sources/hessian-connection-pencil-context.json'
known={e['id']:e for e in json.loads((R/'MATH_STATUS.json').read_text())['entries']}
assert 'HC4MYMI0' not in known,'registered checker replay forbidden'
assert not receipt.exists(),'receipt exists; never replay'
c=json.loads(cp.read_text());start=time.perf_counter();checks=[];alg=0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,v,a=False):
 global alg
 assert bool(v),name
 checks.append(name);alg+=int(a)
def zero(M):return all(s.cancel(q)==0 for q in M)
for p,d in c['retained_proof_and_source_sha256'].items():check('sealed '+p,sha(R/p)==d)
x,y,r,tau,a,b,lam=s.symbols('x y r tau a b lam')
z=s.Matrix([x,y]);w=s.Matrix([r,tau]);coords=[x,y,r,tau]
u=x+a*y*y;v=y+b*u*u;K=s.Matrix([u,v]);A=K.jacobian(z);Ai=A.inv()
check('all-parameter base determinant',s.expand(A.det())==1,True)
U,V,Xi,Eta=s.symbols('U V Xi Eta');yi=V-b*U**2;xi=U-a*yi**2
check('base inverse first composition',zero(K.subs({x:xi,y:yi},simultaneous=True)-s.Matrix([U,V])),True)
check('base inverse second composition',zero(s.Matrix([xi,yi]).subs({U:u,V:v},simultaneous=True)-z),True)
f=s.Function('D')(x,y);h=f+r*u+tau*v;H=s.hessian(h,coords)
C=s.hessian(f,[x,y])+r*s.hessian(u,[x,y])+tau*s.hessian(v,[x,y])
Hb=C.row_join(A.T).col_join(A.row_join(s.zeros(2)))
Hi=s.zeros(2).row_join(Ai).col_join(Ai.T.row_join(-Ai.T*C*Ai))
check('arbitrary-D Hessian blocks',zero(H-Hb),True)
check('arbitrary-D inverse Hessian blocks',zero(H*Hi-s.eye(4)),True)
check('unit determinant by explicit block operations',s.cancel(H.det())==1,True)
check('gradient output block formula',zero(s.Matrix([s.diff(h,q) for q in coords])-s.Matrix(list(s.Matrix([s.diff(f,q) for q in z])+A.T*w)+list(K))),True)
px,py,pr,ps=s.symbols('px py pr ps');vel=s.Matrix([px,py,pr,ps]);zvel=vel[:2,0]
third=s.Matrix([sum(s.diff(h,coords[i],coords[j],coords[k])*vel[j]*vel[k] for j in range(4) for k in range(4)) for i in range(4)])
acc=-lam*Hi*third;k2=s.Matrix([(zvel.T*s.hessian(ki,[x,y])*zvel)[0] for ki in K])
check('full base acceleration independent of D and w',zero(acc[:2,0]+lam*Ai*k2),True)
check('full dual acceleration linear in dual position',all(s.cancel(s.diff(acc[i],p,q))==0 for i in [2,3] for p in [r,tau,pr,ps] for q in [r,tau,pr,ps]),True)
ua=s.cancel(acc[0]+2*a*py**2+2*a*y*acc[1]);up=px+2*a*y*py
check('all-parameter u geodesic equation',s.cancel(ua-2*a*(1-lam)*py**2)==0,True)
check('all-parameter y geodesic equation',s.cancel(acc[1]+2*b*lam*up**2)==0,True)
alpha,beta,cc,kappa=s.symbols('alpha beta cc kappa',nonzero=True)
T=s.symbols('T');up0=1/(kappa*(1-T));yp0=cc*up0
check('all-parameter logarithmic velocity u equation',s.cancel((s.diff(up0,T)-alpha*yp0**2).subs(kappa,alpha*cc**2))==0,True)
check('all-parameter logarithmic velocity y equation',s.cancel((s.diff(yp0,T)-beta*up0**2).subs({kappa:alpha*cc**2,beta:alpha*cc**3},simultaneous=True))==0,True)
Gu=s.Matrix([[0,0],[2*b*lam,0]]);Gy=s.Matrix([[0,-2*a*(1-lam)],[0,0]])
check('all-parameter curvature commutator',zero(Gu*Gy-Gy*Gu-4*a*b*lam*(1-lam)*s.diag(1,-1)),True)
# A full actual geodesic control, beyond the projected equations.
L=-s.log(1-T);curve=s.Matrix([L-L**2,-L,0,0]);vv=curve.diff(T)
h0=(h-f).subs({a:1,b:1});H0=s.hessian(h0,coords);Hi0=Hi.subs({a:1,b:1}).subs({s.diff(f,x,x):0,s.diff(f,x,y):0,s.diff(f,y,y):0})
third0=s.Matrix([sum(s.diff(h0,coords[i],coords[j],coords[k]).subs(dict(zip(coords,curve)),simultaneous=True)*vv[j]*vv[k] for j in range(4) for k in range(4)) for i in range(4)])
check('full Levi-Civita finite-time geodesic',zero(curve.diff(T,2)+s.Rational(1,2)*Hi0.subs(dict(zip(coords,curve)),simultaneous=True)*third0),True)
check('logarithmic base coordinate recovered',s.cancel(curve[0]+curve[1]**2-L)==0,True)
check('zero section preserved for D zero',zero(acc[2:,0].subs({r:0,tau:0,pr:0,ps:0}).subs({s.diff(f,x,x):0,s.diff(f,x,y):0,s.diff(f,y,y):0,s.diff(f,x,x,x):0,s.diff(f,x,x,y):0,s.diff(f,x,y,y):0,s.diff(f,y,y,y):0})),True)
for i,q in enumerate(coords):
 check('all-D preserved volume '+str(i),s.cancel(s.trace(Hi*H.diff(q)))==0,True)
check('no external import or primary cache',c['external_imports']==[] and c['new_primary_caches']==[])
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=alg,elapsed_seconds=time.perf_counter()-start,checks=checks,proof_sha256=sha(Hdir/'HESSIAN_CONNECTION_PENCIL_AND_INCOMPLETE_AUTOMORPHISMS.md'),checker_sha256=sha(Path(__file__)),context_sha256=sha(cp),registered_mathematical_checkers_replayed=[],independent_replay=False,formal_verification=False,external_review=False)
receipt.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
