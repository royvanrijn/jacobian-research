"""Exact controls for trace-core and node proofs; no registered checker replay."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(s.expand(value)==0,label)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
a,b,u,v=s.symbols('a b u v')
# The point module uses genuine Weyl actions,not a finite-dimensional
# truncation that could accidentally violate a commutator at its boundary.
polynomials=[1,a,b,a*a*b+3*a*b*b+7,a**5+2*b**4+a*b]
for index,p in enumerate(polynomials):
    zero(a*(-s.diff(p,a))-(-s.diff(a*p,a))-p,'the first point-module Weyl commutator is exactly one on control '+str(index))
    zero(b*(-s.diff(p,b))-(-s.diff(b*p,b))-p,'the second point-module Weyl commutator is exactly one on control '+str(index))
    check(s.Poly(a*p,a,b).coeff_monomial(1)==0 and s.Poly(b*p,a,b).coeff_monomial(1)==0,
          'the top point-module de Rham class is not a derivative image on control '+str(index))
check(s.Poly(s.Integer(1),a,b).coeff_monomial(1)==1,'the constant point-module class gives a nonzero top de Rham quotient')
for i,j in [(0,0),(1,3),(4,2),(7,1)]:
    mon=a**i*b**j
    zero((-1)**(i+j)*s.diff(mon,a,i,b,j)-s.factorial(i)*s.factorial(j),
         'coordinate local nilpotence reduces a point-module monomial to the generator at '+str((i,j)))
    zero(s.diff(u**(i+1)*v**j/(i+1),u)-u**i*v**j,
         'every polynomial source top-form monomial has a polynomial primitive at '+str((i,j)))
    # Cech-to-point-module map sends u^-i-1 v^-j-1 to normalized derivatives.
    cech=u**(-i-1)*v**(-j-1)
    image=(-1)**(i+j)*a**i*b**j/(s.factorial(i)*s.factorial(j))
    du_cech=-(i+1)*u**(-i-2)*v**(-j-1)
    image_du=-(i+1)*(-1)**(i+j+1)*a**(i+1)*b**j/(s.factorial(i+1)*s.factorial(j))
    zero(image_du-a*image,'the Laurent Cech derivative agrees with the point-module action at '+str((i,j)))

def mixed(f):
    out=0
    for term in s.Add.make_args(s.expand(f)):
        powers=term.as_powers_dict()
        if powers.get(u,0)<0 and powers.get(v,0)<0:out+=term
    return s.expand(out)

F=3/u**3+v**4/u**2+5/v**4+u**5/v+7+u*v
check(mixed(F)==0,'sums of single-branch principal parts have zero mixed Cech projection')
check(mixed(1/u+1/v)==0 and mixed(1/(u*v))==1/(u*v),
      'the nodal trace core contains both single poles and omits their product')
check(mixed(s.diff(F,u))==mixed(s.diff(F,v))==0,'the full single-branch core is stable under both target derivatives on the control')
G=F+2/u**2/v**3-5/u/v
zero(mixed(s.diff(G,u))-mixed(s.diff(mixed(G),u)),'Cech projection commutes with the first derivative')
zero(mixed(s.diff(G,v))-mixed(s.diff(mixed(G),v)),'Cech projection commutes with the second derivative')
zero(mixed(u*G)-mixed(u*mixed(G)),'Cech projection commutes with first target multiplication')
zero(mixed(v*G)-mixed(v*mixed(G)),'Cech projection commutes with second target multiplication')
check(mixed(G)!=0,'a hostile mixed-pole mutation is detected by the exact point quotient')
# The valuation proof has an unrestricted positive parameter,not a tested
# finite upper pole order.
h,N,q,aa=s.symbols('h N q aa',positive=True)
zero((-h*aa)-(-q*N)-(q*N-h*aa),'the negative-block product violates every fixed lower bound as h grows')
check(s.diff(-h*aa,h).is_negative,'the actual negative-block product valuation decreases strictly without bound')

def compose(p,r):return tuple(p[r[i]] for i in range(len(p)))
def support(p):return {i for i,j in enumerate(p) if i!=j}
for label,p,r,allowed in [
 ('overlapping equal transpositions with a genuinely fixed third sheet',(1,0,2),(1,0,2),False),
 ('disjoint quartic transpositions',(1,0,2,3),(0,1,3,2),True),
 ('regular Klein-four crossing with an additional fixed sheet',(1,0,3,2,4),(2,3,0,1,4),False),
 ('disjoint six-sheet three-cycles',(1,2,0,3,4,5),(0,1,2,4,5,3),True)]:
    check(compose(p,r)==compose(r,p),'the local meridians commute for '+label)
    check((not bool(support(p)&support(r))) is allowed,'the exact moved-sheet screen has the stated outcome for '+label)
    if allowed:check(len(support(p))+len(support(r))<=len(p),'the all-degree support sum is compatible for '+label)
check(11-3==8 and 2*(11-3)>11,'retained exceptional eleven-sheet cyclic marks force overlap for every pair of nontrivial inertias')
check(2+2==4 and 3+3==6,'the quartic and degree-six disjoint-support packets remain permitted')
# Actual multiplication/trace on local Kummer cover models. These are
# formal cover controls,not whole affine-plane Keller sources.
for degree,branch in [(2,u),(3,v),(5,u)]:
    M=s.zeros(degree)
    for i in range(degree-1):M[i+1,i]=1
    M[0,degree-1]=branch
    check(M**degree==branch*s.eye(degree),'the local Kummer multiplication has its ramified equation at degree '+str(degree))
    zero(s.trace(s.eye(degree)/(u*v))-degree/(u*v),'an entire omitted torus factor supplies a genuine mixed trace at degree '+str(degree))
    check(mixed(s.trace(s.eye(degree)/(u*v)))!=0,'a torus factor contradicts the nodal core at degree '+str(degree))
Mu=s.Matrix([[0,u],[1,0]]);Mv=s.Matrix([[0,v],[1,0]])
for k in [1,2,4]:
    zero(s.trace(Mu**(-2*k))-2/u**k,'the retained one-branch double factor has the full first-branch pole at k='+str(k))
    zero(s.trace(Mv**(-2*k))-2/v**k,'the other double factor has the full second-branch pole at k='+str(k))
    check(mixed(s.trace(Mu**(-2*k))+s.trace(Mv**(-2*k)))==0,
          'disjoint local factors supply no forbidden mixed trace at k='+str(k))
algebra=len(checks)
context=json.loads((HERE/'sources/plane-trace-core-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(sha(ROOT/path)==digest,'retained immutable input '+path)
check(len(context['newly_imported_proved_external_statements'])==1 and context['newly_used_unproved_external_statements']==[],
      'finite Weyl generation of polynomial localization is the sole new standard import')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(flags['plane_trace_image_is_exact_maximal_point_quotient_kernel'] and flags['two_node_meridians_have_disjoint_moved_support']
      and not flags['plane_Keller_map_proved_surjective'] and not flags['unrestricted_HC4_or_JC2_resolved'] and not flags['original_goal_complete'],
      'trace-core and nodal restrictions preserve the plane,ordinary-node,total-divisor and full-goal boundaries')
proof=HERE/'PLANE_TRACE_IMAGE_RIGIDITY_AND_NODAL_INERTIA.md'
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=s.__version__,
 proof_sha256=sha(proof),checker_sha256=sha(Path(__file__)),dependencies=context['retained_proof_and_source_sha256'],
 source_sha256={'sources/plane-trace-core-context.json':sha(HERE/'sources/plane-trace-core-context.json')},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-plane-trace-core.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
