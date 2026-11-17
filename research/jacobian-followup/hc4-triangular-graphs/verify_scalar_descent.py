"""LE0/1 exact identity/provenance controls,not a universal proof verifier."""
import hashlib,json,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(name,c):
 assert c,name
 checks.append(name)
def zero(name,e):check(name,s.cancel(s.expand(e))==0)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
x,y,a,b,d,w,N,U,V,z,t,q=s.symbols('x y a b d w N U V z t q')
zero('norm ratio derivative cancels in every degree',s.diff(w**d/N,w)*a*w+s.diff(w**d/N,N)*d*a*N)
for n in [2,3,4]:
 Z=s.symbols('z0:'+str(n));mean=sum(Z)/n
 def deriv(f):return sum((a*zz+b)*s.diff(f,zz) for zz in Z)
 zero(f'trace normalized mean solves same affine equation degree{n}',deriv(mean)-a*mean-b)
 for i,zz in enumerate(Z):zero(f'trace centering removes inhomogeneous term degree{n} sheet{i}',deriv(zz-mean)-a*(zz-mean))
zero('enlarged constants control',s.diff(z,V))
for e in [2,3,5]:
 zero(f'Kummer missing Galois hypothesis control e{e}',1/(e*z**(e-1))-z/(e*z**e))
 check(f'Kummer extension nonconstant geometric power e{e}',s.degree(z**e,z)==e)

def jac(f,g):return s.diff(f,x)*s.diff(g,y)-s.diff(f,y)*s.diff(g,x)
for k in [2,3,5]:
 P=x+(y+x**k)**2;Q=y+x**k
 zero(f'nonlinear Keller control Jacobian degree{k}',jac(P,Q)-1)
 A=x*y+x+1;C=y**2+x**2-y
 delta=lambda f:A*s.diff(f,x)+C*s.diff(f,y)
 H=delta(P);J=delta(Q)
 zero(f'one component divergence identity degree{k}',jac(H,Q)+jac(P,J)-s.diff(A,x)-s.diff(C,y))

# The constant-Jacobian hypothesis cannot be dropped from that identity.
P=x;Q=y**2+x*y;A=0;C=x*y;j=jac(P,Q)
delta=lambda f:A*s.diff(f,x)+C*s.diff(f,y)
zero('nonKeller identity retains logarithmic Jacobian correction',(jac(delta(P),Q)+jac(P,delta(Q)))/j-(s.diff(A,x)+s.diff(C,y))-delta(j)/j)
check('nonKeller correction genuinely nonzero',s.cancel(delta(j)/j)!=0)

h=s.Function('h');H=h(x*y);A=x*H;C=-y*H
delta=lambda f:A*s.diff(f,x)+C*s.diff(f,y)
zero('all Hxy complete replicas divergence zero',s.diff(A,x)+s.diff(C,y))
zero('all Hxy complete replicas invariant xy',delta(x*y))
for n in [1,2,3]:zero(f'nonLF iterate formula order{n}',delta(x*H**(n-1))-x*H**n)
ht=(x*y)**2+x*y+1;fx=x*s.exp(t*ht);fy=y*s.exp(-t*ht)
zero('entire complete replica flow invariant xy',fx*fy-x*y)
zero('entire complete replica flow unit Jacobian',s.simplify(jac(fx,fy))-1)
zero('entire complete replica flow x velocity',s.diff(fx,t)-fx*ht)
zero('entire complete replica flow y velocity',s.diff(fy,t)+fy*ht)
for m,n in [(1,1),(2,3),(5,2)]:
 P=x*h(x**m*y**n)
 zero(f'single component weighted identity weights{n},{-m}',n*x*s.diff(P,x)-m*y*s.diff(P,y)-n*P)
P=x+x*x*y
zero('nonsingular noncoordinate component gradient Bezout',(1-2*x*y)*s.diff(P,x)+4*y*y*s.diff(P,y)-1)
zero('same nonsingular component has one Euler eigenvalue',x*s.diff(P,x)-y*s.diff(P,y)-P)

P=x+y*y;Q=y;P1=2*P+3;Q1=-Q+P**2+1
sx=P1-Q1**2;sy=Q1
zero('one pencil compatible source automorphism P identity',sx+sy**2-P1)
zero('one pencil compatible source automorphism constant Jacobian',jac(sx,sy)+2)
zero('one pencil triangular mate difference',sy+Q-(P**2+1))
ALG=len(checks)
context=json.loads((HERE/'sources/scalar-descent-context.json').read_text())
for p,v in context['retained_proof_and_source_sha256'].items():check('retained bytes '+p,sha(ROOT/p)==v)
for p,v in context['new_primary_sha256'].items():check('new primary bytes '+p,sha(HERE/p)==v)
primary=(HERE/'sources/moskowicz-centralizer-1802.04685v2.html').read_text()
check('classical centralizer primary attribution and theorem retained','Theorem 2.3' in primary and 'Moskowicz' in primary and 'Cheng' in primary)
proof=HERE/'SCALAR_DIFFERENTIAL_RIGIDITY_AND_ONE_COMPONENT_DESCENT.md';flat=' '.join(proof.read_text().split())
for p in ['CURRENT LAST UNSOLVED IMPLICATION','L^D=C(P)','NO NONTRIVIAL INTERMEDIATE FIELD GALOIS OVER F','Dz=a*z+b','D^n z in F','div(delta) in F','JUST ONE','BOTH component derivatives','NOT locally finite','SIX-dimensional','NOT F-linear independence','FINITE ABELIAN','ENTIRE nonproperness curve','SAME GENUS','BOTH components are rational','UNPROVED','DEFERRED']:
 check('written scope '+p,p in flat)
check('no new external theorem import',context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[])
flags={k:v for k,v in context.items() if isinstance(v,bool)}
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=ALG,checks=checks,proof_sha256=sha(proof),checker_sha256=sha(Path(__file__)),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-scalar-descent.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(outcome=result['outcome'],checks=len(checks),algebraic_controls=ALG,elapsed_seconds=result['elapsed_seconds']),indent=2))
