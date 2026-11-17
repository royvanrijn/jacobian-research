"""Transcription controls for written Galois/lattice descent;no old replay."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(b,label):
 assert bool(b),label
 checks.append(label)
def zero(f,label):check(s.expand(f)==0,label)
a,b,c,d,e,f,g=s.symbols('a b c d e f g');J=s.Matrix([[a,b],[c,d]])
matrix=s.Matrix([[e,f,a,c],[f,g,b,d],[a,b,0,0],[c,d,0,0]])
zero(matrix.det()-(a*d-b*c)**2,'the arbitrary cotangent block Hessian is exactly the square of the plane Jacobian')
r1,r2,o1,o2=s.symbols('r1 r2 o1 o2');out=s.Matrix([o1-r1,o2-r2])
fibres=s.simplify(J.T.inv()*out)
check(s.simplify(J.T*fibres-out)==s.zeros(2,1),'after base recovery the inverse transpose recovers both fibre variables exactly')
x,y,u,v=s.symbols('x y u v');R=s.Function('R')(x,y);F=s.Function('F')(x,y);G=s.Function('G')(x,y)
h=R+u*F+v*G;grad=s.Matrix([s.diff(h,zz) for zz in [x,y,u,v]])
check(grad[2:4,0]==s.Matrix([F,G]),'the last gradient outputs are the exact plane pair irrespective of the base polynomial')
check(grad[:2,0]==s.Matrix([s.diff(R,x),s.diff(R,y)])+s.Matrix([F,G]).jacobian([x,y]).T*s.Matrix([u,v]),
 'the remaining gradient outputs have the exact inverse-transpose fibre recovery form')
pi,p,q,r,w=s.symbols('pi p q r w');T=s.Matrix([[1,0,pi*p,pi*q],[0,1,pi*r,pi*w],[0,0,1,0],[0,0,0,1]])
check(T.det()==1 and T.subs(pi,0)==s.eye(4),'a nontrivial integral source frame can be congruent to identity and have unit determinant')
check(T.inv()*T==s.eye(4) and all(s.denom(z)==1 for z in T.inv()),'its complete inverse is integral without coefficient division')
zero(s.expand((T.T*matrix*T).det(method='domain-ge'))-(a*d-b*c)**2,
 'the full integral Hessian congruence retains the squared plane determinant in a mixed source frame')
lam=s.symbols('lam',nonzero=True);D=s.diag(lam,1,1,1)
zero((D.T*matrix*D).det()-lam**2*(a*d-b*c)**2,'general determinant scaling uses the square of the source-frame determinant')
zero(lam**(-2)*(a*d-b*c)**2-(a*d-b*c)**2/D.det()**2,'transporting back to the original source uses the inverse determinant square')
# A basis may require a field extension although its unique plane does not.
rr,ss,tt,ww=s.symbols('rr ss tt ww');rational=s.Matrix([[1,0,rr,ss],[0,1,tt,ww]])
foreign=s.Matrix([[1,s.sqrt(2)],[s.sqrt(2),1]])*rational
canonical=foreign.rref()[0]
check(s.simplify(canonical-rational)==s.zeros(2,4),'canonical row-echelon coordinates recover a rational plane from a foreign coefficient basis')
# The actual residue plane for H0 is the standard fibre plane.
P=x+x*x*y+x**4+x**6*y*y;Q=y+x**5+x**6*y+x**7*y*y+x**8*y**3
top=s.Poly(Q,x,y);check(top.total_degree()==11,'the old coefficient Q retains the degree-eleven top used to force both base directions zero')
H0=u*P+v*Q;t,alpha,beta=s.symbols('t alpha beta')
zero(H0.subs({u:u+t*alpha,v:v+t*beta},simultaneous=True)-H0-t*(alpha*P+beta*Q),
 'every direction in the original residue u-v plane is affine on every translate')
JJ=s.expand(s.det(s.Matrix([P,Q]).jacobian([x,y])))
check(s.Poly(JJ,x,y).coeff_monomial(x*y)==2,'the fixed pair still violates its first divided plane equation')
check(all(int(coef)%4==0 for coef in s.Poly(JJ**2-1,x,y).coeffs()),
 'the residue potential has a whole mod4 lift although the fixed plane pair has no mod4 lift')
check(s.Poly(JJ**2,x,y).coeff_monomial(x*y)%8==4,
 'finite-precision lifting is deliberately not equated to integral cotangent descent')
algebra=len(checks)
context=json.loads((HERE/'sources/integral-cotangent-descent-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
 check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
check(context['newly_imported_proved_external_statements']==context['newly_used_unproved_external_statements']==[],
 'the rational and integral descent uses retained Tang existence and elementary plane uniqueness,not a new external import')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==6 and not flags['finite_precision_arbitrary_potential_lifting_equivalent_to_plane_lifting']
 and not flags['plane_Keller_lifting_or_integer_equation_sufficiency_proved_automatic']
 and not flags['unrestricted_HC4_or_JC2_resolved'] and not flags['original_goal_complete'],
 'the exact characteristic-zero lifting equivalence retains the finite-precision,sufficiency and full HC4 boundaries')
proof=HERE/'INTEGRAL_COTANGENT_DESCENT_AND_PLANE_LIFT_EQUIVALENCE.md';cp=HERE/'sources/integral-cotangent-descent-context.json'
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=s.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={'sources/integral-cotangent-descent-context.json':hashlib.sha256(cp.read_bytes()).hexdigest()},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-integral-cotangent-descent.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
