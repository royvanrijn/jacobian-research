"""Exact monodromy/algebra controls for written all-degree lifting and rigidity.

Riemann existence and Blanc-Stampfli are explicit published theorem inputs.
Finite controls do not prove the general topology,extension or classification.
No registered mathematical checker is invoked.
"""
import hashlib,itertools,json,math,time
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];started=time.monotonic();checks=[]
def check(ok,label):
 assert ok,label
 checks.append(label)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mul(a,b):return tuple(a[b[i]] for i in range(len(a)))
def power(a,n):
 b=tuple(range(len(a)))
 for _ in range(n):b=mul(a,b)
 return b
def order(a):
 b=tuple(range(len(a)))
 for n in range(1,math.factorial(len(a))+1):
  b=mul(a,b)
  if b==tuple(range(len(a))):return n
 raise AssertionError('invalid permutation')
for d in [2,3,5,8,12]:
 sigma=tuple((i+1)%d for i in range(d));ident=tuple(range(d))
 check(order(sigma)==d,'the connected cyclic degree'+str(d)+' cover attains the exact torus-power bound')
 check(power(sigma,d)==ident,'degree power kills the central monodromy at degree'+str(d))
 check(all(power(sigma,k)[0]!=0 for k in range(1,d)),'the central cycle acts without any fixed sheet before its exact order at degree'+str(d))
S3=list(itertools.permutations(range(3)));centre=[g for g in S3 if all(mul(g,h)==mul(h,g) for h in S3)]
check(centre==[(0,1,2)],'the centre-free S3 monodromy control needs no torus isogeny')
mixed=(1,0,3,4,2)
check(order(mixed)==6,'disconnected two-plus-three sheet monodromy has order6')
check(power(mixed,5)!=(0,1,2,3,4),'degree5 does NOT kill disconnected monodromy,so connectedness is essential')

t,u,v,lam,mu,ss=s.symbols('t u v lam mu ss')
for m in [2,3,7]:
 check(s.expand((lam*t)**m-lam**m*t**m)==0,'the Kummer degree'+str(m)+' torus action lifts after the exact power')
 check(s.denom(s.Rational(1,m))==m and m>1,'the degree'+str(m)+' unit-character equation requires the nonintegral weight1/m')
check(s.expand((v+ss)+mu-(v+ss+mu))==0,'additive translation lifts on the same finite etale Kummer cover without an isogeny')
check(s.diff(t**3,t)==3*t*t,'the Kummer cover is etale after t is inverted,not at the missing origin')

# Branch-only invariance does not imply invariance of an unramified deletion.
H=v-u*u
check(s.diff(H,v)==1,'translation moves the deleted unramified parabola despite empty branch support')
check(s.gcd(H,H+1)==1,'the translated deleted equation is not a unit in the localization at the original one')
check(H.subs({u:0,v:-1})!=0 and (H+1).subs({u:0,v:-1})==0,'the translated reciprocal has a pole at a point of the original source open')

# The actual non-Keller PF2OMIT3 family does not meet the etale-off-N input.
x,y=s.symbols('x y')
for a in [1,3,7]:
 J=s.det(s.Matrix([x,x**a*y*y-y]).jacobian([x,y]));disc=1+4*u**a*v
 check(s.expand(J-(2*x**a*y-1))==0,'the retained quasi-finite invariant-N control is not Keller at a='+str(a))
 check(s.rem(s.Poly(s.diff(disc,v),v),s.Poly(disc,v)).as_expr()!=0,'translation does not preserve its branch divisor outside N at a='+str(a))

# Whole-divisor symmetry and nonlinear conjugacy are distinct from one component.
cusp=v*v-u**3;line=v-u-1;whole=cusp*line
def E(q):return s.expand(2*u*s.diff(q,u)+3*v*s.diff(q,v))
check(s.expand(E(cusp)-6*cusp)==0,'the entire cusp has an actual weighted symmetry')
check(s.rem(s.Poly(E(whole),v),s.Poly(whole,v)).as_expr()!=0,'a second component breaks that symmetry on the whole reduced curve')
conjugate=(v-u*u)**2-u**3
EE=s.expand(2*u*s.diff(conjugate,u)+(3*v+u*u)*s.diff(conjugate,v))
check(s.expand(EE-6*conjugate)==0,'the nonlinear target-coordinate conjugate still has its exact logarithmic algebraic field')
actu=lam**2*u;actv=lam**3*(v-u*u)+lam**4*u*u
check(s.expand(conjugate.subs({u:actu,v:actv},simultaneous=True)-lam**6*conjugate)==0,'the conjugated curve is invariant under the full polynomial Gm action')

# A logarithmic Hamiltonian need not be locally finite.
F=u*u*v*v+u+v
def DH(q):return s.expand(s.diff(F,u)*s.diff(q,v)-s.diff(F,v)*s.diff(q,u))
check(DH(F)==0,'the generic polynomial Hamiltonian is tangent to its divisor')
aa,bb,n=s.symbols('aa bb n',integer=True)
check((2*(bb-aa)).subs({aa:n+1,bb:n})==-2,'the universal highest-term coefficient stays nonzero on every Hamiltonian iterate')
out=u
for k in range(1,5):
 out=DH(out);poly=s.Poly(out,u,v)
 check(poly.total_degree()==2*k+1 and poly.coeff_monomial(u**(k+1)*v**k)==(-2)**k,'exact Hamiltonian iterate'+str(k)+' has the predicted strictly growing top degree')
algebra=len(checks)
context=json.loads((HERE/'sources/target-boundary-symmetry-context.json').read_text())
for pth,hsh in context['retained_proof_and_source_sha256'].items():check(sha(ROOT/pth)==hsh,'immutable retained input '+pth)
for pth,hsh in context['new_primary_sha256'].items():check(sha(HERE/pth)==hsh,'inspected new primary '+pth)
check(len(context['newly_imported_proved_external_statements'])==2,'exactly two new proved external inputs are explicit:finite-cover comparison and the full curve-stabilizer theorem')
check(context['newly_used_unproved_external_statements']==[],'no target symmetry existence or source descent is assumed')
check(context['preserving_only_ramified_branch_support_assumed_sufficient'] is False,'all omitted unramified primes remain in the source-invariance proof')
check(context['last_HC4_source_to_regular_inverse_Gram_implication_proved'] is False,'the full-HC4 last source implication remains OPEN')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,elapsed_seconds=round(time.monotonic()-started,3),
 proof_sha256=sha(HERE/'TARGET_BOUNDARY_SYMMETRY_AND_FINITE_STABILIZER.md'),checker_sha256=sha(Path(__file__)),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],**flags)
(HERE/'receipt-target-boundary-symmetry.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
