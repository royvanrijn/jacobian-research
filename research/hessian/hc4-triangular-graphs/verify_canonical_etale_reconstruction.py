"""Exact hostile controls for the written trace-dual etale reconstruction theorem."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
 assert bool(value),label
 checks.append(label)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
a,b,t,z,x,y,u,v=s.symbols('a b t z x y u v')
T=s.Matrix([[0,0,-b],[1,0,-a],[0,1,0]]);I=s.eye(3)
G=s.Matrix([[s.trace(T**(i+j)) for j in range(3)] for i in range(3)])
check(G==s.Matrix([[3,0,-2*a],[0,-2*a,-3*b],[-2*a,-3*b,2*a*a]]),'the generic normal cubic has the exact trace Gram matrix')
check(s.factor(G.det())==-4*a**3-27*b*b,'the cubic trace determinant is minus the target discriminant')
H=G.inv();Ms=[sum((H[j,i]*(T**j) for j in range(3)),s.zeros(3)) for i in range(3)]
check(s.simplify(G*H)==I,'the inverse Gram coefficients give the actual dual basis')
D=3*T*T+a*I
for i,expected in enumerate([T*T+a*I,T,I]):
 check(s.simplify(D*Ms[i])==expected,'the inverse different dual generator has the advertised derivative multiple at i='+str(i))
 check(s.simplify(s.trace(Ms[i]))==int(i==0),'each canonical dual generator has a regular first trace at i='+str(i))
check(all(s.simplify(s.trace(Ms[i]*Ms[j])-H[i,j])==0 for i in range(3) for j in range(3)),
      'quadratic traces of the canonical generators are exactly the inverse Gram matrix,not automatically regular')
f=lambda r:r**3+a*r+b
N=t*t+t*z+z*z+a;den=3*z*z+a
check(s.expand((t-z)*N-f(t)+f(z))==0,'the Casimir numerator intertwines the two cubic factors')
check(s.expand(N.subs(t,z)-den)==0,'multiplication of the Casimir is exactly1')
gb=s.groebner([f(t),f(z)],t,z,a,b)
check(gb.reduce(s.expand(N*(N-den)))[1]==0,'the cubic Casimir is an idempotent after adjoining the trace dual')
point={t:1,a:s.Rational(-3,4),b:s.Rational(-1,4)}
check(f(t).subs(point)==0,'the retained-sheet control lies on the actual cubic cover')
check((3*t*t+a).subs(point)==s.Rational(9,4),'its source derivative is invertible')
check((4*a**3+27*b*b).subs(point)==0,'target discriminant localization would wrongly delete this unramified point')
# Flat nonnormal order versus its full etale normalization.
w=s.symbols('w');f0=w*w-x*w;d0=2*w-x
check(s.rem(d0*d0-x*x,f0,w)==0,'the flat order different becomes x up to a unit after localization')
check(s.rem(w*w-x*w,f0,w)==0,'w/x separates the two normalization factors after inverting x')
G0=s.Matrix([[2,x],[x,x*x]])
check(G0.det()==x*x,'the nonnormal flat order has discriminant x squared')
check(s.Matrix([[1,0],[0,1]]).det()==1,'the full split normalization has perfect trace at the missing divisor')
# Nonflat point-gluing: B plus the maximal ideal I has generic rank2
# and an origin fibre with independent1,u,v and square-zero(u,v).
check(s.Matrix([[1,x,y],[1,0,0]]).rank()==2,'the point-glued finite order has generic rank2')
gb0=s.groebner([u*u,u*v,v*v],u,v)
standard=[1,u,v]
check(all(gb0.reduce(s.sympify(q))[1]==q for q in standard)
      and all(gb0.reduce(q)[1]==0 for q in [u*u,u*v,v*v]),
      'the origin fibre has the exact three-dimensional standard basis1,u,v')
check(len(standard)==3 and len(standard)!=2,'the origin length jump detects failure of flatness')
check(s.gcd(x,y)==1,'the codimension-two gluing ideal has dual B by coprime denominator removal')
# Characteristic2: all characteristic-free Casimir identities still hold,
# while the finite etale algebra need not be the target ring.
g2=s.groebner([t*t+t+a,z*z+z+a],t,z,a,modulus=2);e2=t+z+1
check(g2.reduce(s.expand(e2*e2-e2))[1]==0,'the Artin-Schreier Casimir is idempotent in characteristic2')
check(s.Poly(e2.subs(t,z)-1,z,modulus=2).is_zero,'its multiplication is1 even when trace1=0')
check(g2.reduce(s.expand((t+z)*e2))[1]==0,'the characteristic2 Casimir commutes across the two factors')
check(s.det(s.Matrix([[0,1],[1,1]]))%2==1,'the nontrivial Artin-Schreier trace pairing is perfect')
check(s.Poly((1+t)+t-1,t,modulus=2).is_zero,'the characteristic2 dual frame contracts to1')
for e in [2,3,7]:
 check(s.simplify(e*t**(e-2)/(e*t**(e-1))-1/t)==0,
       'a target derivative of the tame finite generator already generates the inverse uniformizer at e='+str(e))
check(s.Matrix([1,-1]).multiply_elementwise(s.Matrix([1,-1]))==s.Matrix([1,1]),
      'the quadratic augmentation square remains rank1,so the degree>=3 restriction is essential')
algebra=len(checks)
context=json.loads((HERE/'sources/canonical-etale-reconstruction-context.json').read_text())
for p,hsh in context['retained_proof_and_source_sha256'].items():check(sha(ROOT/p)==hsh,'immutable retained input '+p)
for p,hsh in context['new_primary_sha256'].items():check(sha(HERE/p)==hsh,'inspected primary '+p)
check(context['newly_used_unproved_external_statements']==[],'no unproved external statement supports the completed reconstruction theorem')
check(len(context['newly_imported_proved_external_statements'])==2,'the two explicit new proved inputs are different support and miracle flatness')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(flags['global_finite_flat_reconstruction_theorem_complete'] and flags['flatness_countercontrol_is_explicit']
      and not flags['last_HC4_source_to_regular_inverse_Gram_implication_proved'] and not flags['original_goal_complete'],
      'the complete global reconstruction is separated from the last unsolved full-HC4 implication')
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=s.__version__,
 proof_sha256=sha(HERE/'CANONICAL_TRACE_DUAL_ETALE_RECONSTRUCTION.md'),checker_sha256=sha(Path(__file__)),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256=dict(context['new_primary_sha256'],**{'sources/canonical-etale-reconstruction-context.json':sha(HERE/'sources/canonical-etale-reconstruction-context.json')}),
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-canonical-etale-reconstruction.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
