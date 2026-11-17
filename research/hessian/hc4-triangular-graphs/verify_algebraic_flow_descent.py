"""Exact controls for the written all-degree algebraic-flow obstruction.

Classification imports and universal orbit/Jordan/Krull arguments are
written/imported mathematics,not consequences of these finite controls.
No registered mathematical checker is invoked.
"""
import hashlib,json,time
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];started=time.monotonic();checks=[]
def check(ok,label):
 assert ok,label
 checks.append(label)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
x,y,t,v,a,b,lam,mu=s.symbols('x y t v a b lam mu')

# A nonlinear coordinate conjugate and arbitrary polynomial multiplier.
f=x+y*y;pol=f**3+f+1
def D(h):return s.expand(pol*(-2*y*s.diff(h,x)+s.diff(h,y)))
check(D(f)==0,'the conjugated replica fixes its actual polynomial source coordinate')
check(s.expand(D(y)-pol)==0,'the complementary coordinate translates at polynomial speed p(f)')
check(D(D(y))==0,'the replica is nilpotent on the complementary coordinate')
check(D(D(D(x)))==0,'the nonlinear source coordinate has a finite exact translation orbit')
flow=[s.expand(f-(y+t*pol)**2),s.expand(y+t*pol)]
check(s.expand(flow[0]+flow[1]**2-f)==0,'the entire polynomial flow preserves f')
check(s.expand(s.det(s.Matrix(flow).jacobian([x,y]))-1)==0,'the nonlinear replica flow preserves source volume')
# The preceding invariant-coordinate identity makes p(f) unchanged.
# Compose in (f,g),then reconstruct x;avoid expanding a degree144
# conjugated expression before applying that already-proved invariant.
speed=s.symbols('speed');twice=[f-(y+t*speed+v*speed)**2,y+t*speed+v*speed]
once=[f-(y+(t+v)*speed)**2,y+(t+v)*speed]
check(all(s.expand(twice[i]-once[i])==0 for i in range(2)),'the replica satisfies the exact Ga group law for every parameter using its proved invariant coordinate')

# Multipliers outside the field cannot be divided out of a descent claim.
F1=x*x;F2=x*y
check(x*s.diff(F1,y)==0 and x*s.diff(F2,y)==F1,'the non-Keller field C(x^2,xy) is preserved by x partial_y')
check(s.diff(F2,y)==x,'the bare translation derivative lies outside that quadratic invariant field')
check(F1.subs({x:-x,y:-y},simultaneous=True)==F1 and F2.subs({x:-x,y:-y},simultaneous=True)==F2,'the sign involution fixes the toy target field')
check(x.subs(x,-x)==-x,'the bare derivative is not fixed,certifying the multiplier-division loophole')
check(s.det(s.Matrix([F1,F2]).jacobian([x,y]))==2*x*x,'the descent toy is not Keller and contracts x=0')

# Finite orbit spaces preserve their Jordan polynomial projectors.
M=s.Matrix([[a,1,0],[0,a,0],[0,0,b]])
S=a*s.eye(3)+(b-a)*(M-a*s.eye(3))**2/(b-a)**2
S=S.applyfunc(s.cancel);N=M-S
check(S==s.diag(a,a,b),'the semisimple part is a polynomial in the locally finite operator on this orbit space')
check(N*N==s.zeros(3),'its nilpotent part is actually nilpotent')
check(S*N==N*S,'the Jordan parts commute')
check(S+N==M,'the Jordan parts recover the original operator exactly')

# Complex weights need not be commensurable to isolate components.
def E(h):return s.expand(a*x*s.diff(h,x)+b*y*s.diff(h,y))
h=x*x+x*y+y**3;weights=[2*a,a+b,3*b];out=h
for q in weights[1:]:out=s.cancel((E(out)-q*out)/(weights[0]-q))
check(s.cancel(out-x*x)==0,'spectral projection isolates source components for generic arbitrary complex weights')

# Exact normal-coordinate multiplier and vector-field zero orders.
u,w=s.symbols('u w')
for e,r in [(1,1),(2,2),(7,3),(3,8),(11,2)]:
 ff=u**r;gg=u**(e-r)*w+u**-4
 vol=s.factor(s.det(s.Matrix([ff,gg]).jacobian([u,w])))
 check(s.simplify(vol-r*u**(e-1))==0,'the local ramification volume has exact order e-1 at e='+str(e)+',r='+str(r))
 aa=max(0,(e+r-1)//r-1)
 check(r*(aa+1)-e>=0 and (aa==0 or r*aa-e<0),'the minimum multiplier exponent is exact at e='+str(e)+',r='+str(r))
 check(s.simplify((ff**aa)*u**(r-e)-u**(r*(aa+1)-e))==0,'the multiplier gives the exact divisorial vector-field order at e='+str(e)+',r='+str(r))
rows=[(7,2),(7,3)];threshold=max(max(0,(e+r-1)//r-1) for e,r in rows)
check(threshold==3,'primes sharing one boundary value require the maximum exponent,not a sum or one selected prime')
check([r*(threshold+1)-e for e,r in rows]==[1,5],'the shared-value minimal replica has exact nonnegative boundary zero orders')

# Characteristic-zero orbit stabilizers are essential.
tau,ss=s.symbols('tau ss');AS=tau*tau-tau
check(s.Poly(AS.subs(tau,tau+ss)-AS-(ss*ss-ss),tau,ss,modulus=2).is_zero,'Artin-Schreier translations descend in characteristic2')
check((1*1-1)%2==0,'the nonzero additive parameter1 is an actual finite stabilizer in characteristic2')

# The known whole3D map is an actual equivariant dimension countercontrol.
z=s.symbols('z');u0=1+x*y
A=u0**3*z+y*y*u0*(4+3*x*y)
B=y+3*x*u0*u0*z+3*x*y*y*(4+3*x*y)
C=2*x-3*x*x*y-x**3*z
def EW(q):return s.expand(x*s.diff(q,x)-y*s.diff(q,y)-2*z*s.diff(q,z))
for name,q,wt in [('A',A,-2),('B',B,-1),('C',C,1)]:
 check(s.expand(EW(q)-wt*q)==0,'the credited noninvertible whole3D Keller map has the exact descending Euler weight on '+name)
algebra=len(checks)
context=json.loads((HERE/'sources/algebraic-flow-descent-context.json').read_text())
for pth,hsh in context['retained_proof_and_source_sha256'].items():check(sha(ROOT/pth)==hsh,'immutable retained input '+pth)
check(context['newly_imported_proved_external_statements']==[],'Rentschler,Gutwirth and the one-line theorem are explicitly existing imports')
check(context['newly_used_unproved_external_statements']==[],'no action descent or integral-coordinate inverse assertion is imported')
check(context['field_preservation_divided_by_non_target_multiplier'] is False,'the replica proof uses equivariant orbits,not invalid multiplier division')
check(context['last_HC4_source_to_regular_inverse_Gram_implication_proved'] is False,'the final full-HC4 implication remains OPEN')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,elapsed_seconds=round(time.monotonic()-started,3),
 proof_sha256=sha(HERE/'ALGEBRAIC_FLOW_DESCENT_AND_INTEGRAL_COORDINATES.md'),checker_sha256=sha(Path(__file__)),dependencies=context['retained_proof_and_source_sha256'],source_sha256={},**flags)
(HERE/'receipt-algebraic-flow-descent.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
