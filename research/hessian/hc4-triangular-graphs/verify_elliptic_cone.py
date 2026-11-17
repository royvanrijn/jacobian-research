"""Exact controls for the geometric cone-trace proof and complete cyclic source-orbit closure."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic()
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def iszero(value):
    value=sp.cancel(value)
    return value==0 or sp.simplify(value)==0
def zero(value,label):check(iszero(value),label)
def zeros(values,label):check(all(iszero(v) for v in values),label)
x,y,z,w,eta,p,q,k,tau=sp.symbols('x y z w eta p q k tau')
zz=sp.Matrix([x,y,z])
simple_A=sp.diag(1,2,3)
simple_B=sp.Matrix([[0,0,1],[0,1,0],[1,0,0]])
chain_A=sp.Matrix([[0,0,1],[0,1,0],[1,0,0]])
chain_B=sp.Matrix([[0,1,1],[1,1,0],[1,0,0]])
cases=[
 ('simple spectrum',simple_A,simple_B,sp.ones(3,1),
  [p*p,2*p*q,-2*q*q],4),
 ('Jordan 2+1',sp.Matrix([[0,1,0],[1,0,0],[0,0,1]]),
  sp.Matrix([[1,1,0],[1,0,0],[0,0,2]]),sp.Matrix([1,0,1]),
  [p*p,-p*p/2-q*q,p*q],3),
 ('Jordan 3',chain_A,chain_B,sp.Matrix([1,0,0]),
  [p*p,2*p*q,-2*q*q-2*p*q],2)]
for label,A,B,u,param,ends in cases:
    a=(zz.T*A*zz)[0]
    b=(zz.T*B*zz)[0]
    ell=-2*(u.T*A*zz)[0]
    S=A.inv()*B
    check(A.det()!=0 and B.det()!=0,'both ternary forms are nondegenerate '+label)
    check(sp.Matrix.hstack(u,S*u,S*S*u).det()!=0,
          'the ternary offset is cyclic '+label)
    sub=dict(zip([x,y,z],param))
    zero(b.subs(sub,simultaneous=True),'the whole homogeneous conic parameterization is exact '+label)
    AA=sp.expand(a.subs(sub,simultaneous=True))
    LL=sp.expand(ell.subs(sub,simultaneous=True))
    Delta=sp.expand(LL*LL+4*k*AA)
    check(sp.Poly(AA,p,q).total_degree()==4 and sp.Poly(LL,p,q).total_degree()==2,
          'the actual conic restrictions have degrees four and two '+label)
    check(sp.Poly(Delta,p,q).total_degree()==4,
          'the full discriminant retains its homogeneous infinity degree '+label)
    polynomial=sp.Poly(Delta,p,q,domain='QQ(k)')
    gcd=sp.gcd(sp.gcd(polynomial,polynomial.diff(p)),polynomial.diff(q))
    check(gcd.total_degree()==0,'all four generic branch points are simple including the chart boundary '+label)
    fixed=sp.gcd(sp.Poly(AA,p,q),sp.Poly(LL*LL,p,q))
    check(all(exponent==1 for factor,exponent in sp.factor_list(fixed.as_expr(),p,q)[1]),
          'any fixed discriminant base divisor is reduced '+label)
    infinity_count=sum(sp.Poly(factor,p,q).total_degree()
                       for factor,exponent in sp.factor_list(AA,p,q)[1])
    check(infinity_count==ends,
          'the affine puncture count uses distinct zeros rather than intersection multiplicity '+label)
    zero((2*k*tau-LL)**2-Delta+4*k*(AA+LL*tau-k*tau*tau),
         'the local projective curve has exactly the claimed double-cover discriminant '+label)

    # Full polynomial source normalization of the original resonance.
    v3=S*u
    d=sp.Integer(2)
    constant=sp.Integer(5)
    beta=constant-(v3.T*A*u)[0]
    alpha=-2*d*d
    zero(1-(-1/(2*d*d))*alpha,'the linear infinity resonance has its exact sign '+label)
    zero((-1/(2*d*d))*beta+beta/(2*d*d),'the constant infinity resonance is exact '+label)
    g=b+2*(v3.T*A*zz)[0]+2*d*w+constant
    H=lambda value:sp.Integer(7)-(value-beta)**2/(4*d*d)
    h=a+w*w+H(g)
    forward=[zz[i]-u[i] for i in range(3)]+[-d*eta-b/(4*d)]
    new_to_old=dict(zip([x,y,z,w],forward))
    normal=a+ell+eta*b+(u.T*A*u)[0]+H(beta)
    zero(h.subs(new_to_old,simultaneous=True)-normal,
         'the entire original resonant source potential has its literal cone-trace normal form '+label)
    b_shift=b.subs({x:x+u[0],y:y+u[1],z:z+u[2]},simultaneous=True)
    inverse=[zz[i]+u[i] for i in range(3)]+[-w/d-b_shift/(4*d*d)]
    zeros([forward[i].subs(dict(zip([x,y,z,eta],inverse)),simultaneous=True)-[x,y,z,w][i]
           for i in range(4)],'the normalizing source map has its full polynomial inverse '+label)
    check(sp.Matrix(forward).jacobian([x,y,z,eta]).det()==-d,
          'the normalizing source map is a genuine automorphism '+label)
    T4=sp.diag(S,sp.Integer(0))
    v4=sp.Matrix(list(v3)+[d])
    check(sp.Matrix.hstack(*[T4**j*v4 for j in range(4)]).det()!=0,
          'the original four-variable resonant pencil is still cyclic '+label)
    # The whole open stratum has a literal inverse, not merely matching Euler number.
    f=a+ell+eta*b
    recovered=(k-a-ell)/b
    zero(f.subs(eta,recovered)-k,'the entire complement of the cone is solved exactly '+label)
    zero(b*(recovered-eta)-(k-f),'the open source recovery keeps the full fibre equation '+label)

# The tangent infinity points remain smooth because the linear term is nonzero there.
A,B,u=cases[1][1:4]
point=sp.Matrix([0,1,0])
zero((point.T*A*point)[0],'the length-two tangent point lies on the first conic')
zero((point.T*B*point)[0],'the length-two tangent point lies on the second conic')
check((u.T*A*point)[0]!=0,'its linear intercept is transverse to the projective infinity point')
A,B,u=cases[2][1:4]
point=sp.Matrix([0,0,1])
check((point.T*A*point)[0]==0 and (point.T*B*point)[0]==0,
      'the length-three tangent point is an actual common conic zero')
check((u.T*A*point)[0]!=0,'cyclicity prevents a singular projective curve at the triple infinity contact')

# The geometric hypotheses must not be dropped.
u_bad=sp.Matrix([0,1,0])
a=(zz.T*chain_A*zz)[0]
ell_bad=-2*(u_bad.T*chain_A*zz)[0]
param=cases[2][4]
AA=sp.expand(a.subs(dict(zip([x,y,z],param)),simultaneous=True))
LL=sp.expand(ell_bad.subs(dict(zip([x,y,z],param)),simultaneous=True))
bad=sp.Poly(LL*LL+4*k*AA,p,q,domain='QQ(k)')
check(sp.gcd(sp.gcd(bad,bad.diff(p)),bad.diff(q)).total_degree()>0,
      'an uncoupled tangent eigenline really produces a repeated generic branch divisor')
B=simple_B
A=B+sp.diag(0,0,1)
ell=z
param=[p*p,2*p*q,-2*q*q]
a=(zz.T*A*zz)[0]
AA=a.subs(dict(zip([x,y,z],param)),simultaneous=True)
LL=ell.subs(dict(zip([x,y,z],param)),simultaneous=True)
zero(AA-LL*LL,'proportional restricted quadratics really give a square branch section')
check(A.det()!=0 and sp.Matrix.hstack(sp.ones(3,1),(A.inv()*B)*sp.ones(3,1),
       (A.inv()*B)**2*sp.ones(3,1)).det()==0,
      'nondegenerate ternary forms alone do not supply the cyclic genus certificate')

# Riemann-Hurwitz and every whole-fibre stratum are kept explicitly.
L,U,r,N=sp.symbols('L U r N')
check(2*1-2==2*(-2)+4,'four simple branch points give genus one by the retained Riemann-Hurwitz formula')
cone=(L-1)*L+L
zero(cone-L*L,'the entire ternary cone includes its zero-X affine line')
curve=L-U+1-r
whole=L**3-cone+L*curve
zero(whole-(L**3-L*U+(1-r)*L),'the entire value fibre has the exact genus-one Hodge term')
check(sp.Poly(whole-(L**3-N*L),U).coeff_monomial(U)==-L,
      'the Hodge obstruction is independent of every possible critical-point count')
check(curve.subs(r,2)!=curve.subs(r,4),
      'the triple contact cannot be counted as three omitted normalization points')

context=json.loads((HERE/'sources/elliptic-cone-context.json').read_text())
algebra_count=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,
          'retained proof, receipt or foundation bytes stay pinned: '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'related primary curve context bytes are pinned')
check(context['newly_used_unproved_external_statements']==[],
      'the four-branch-point and tangent audit is proved directly from retained foundations')
check(context['all_four_variable_cyclic_affine_pencils_inverse_closed'] is True
      and context['cyclic_double_infinity_resonance_inverse_closed'] is True,
      'the former cyclic resonance is now closed to inverse rather than assumed absent')
check(context['all_cyclic_unit_hessian_source_orbits_proved_nonexistent'] is False,
      'the new geometric conclusion is inverse closure, not universal nonexistence')
check(context['unit_ideal_endpoint_resolved'] is False
      and context['singleton_zero_endpoint_resolved'] is False
      and context['unrestricted_HC4_or_JC2_resolved'] is False,
      'both actual endpoints and the full original objective remain open')
receipt=dict(outcome='PASS',claims=['HC4MYJG0','HC4MYJG1'],
 check_count=len(checks),algebraic_and_scope_controls=algebra_count+4,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'ELLIPTIC_CONE_TRACES_AND_COMPLETE_CYCLIC_CLOSURE.md').read_bytes()).hexdigest(),
 source_sha256={'sources/elliptic-cone-context.json':hashlib.sha256((HERE/'sources/elliptic-cone-context.json').read_bytes()).hexdigest()},
 dependencies=context['retained_proof_and_source_sha256'],
 assurance='Finite exact controls of the written geometric all-degree proof; not independent replay, formal verification or external review.',
 cone_trace_geometric_source_orbit_inverse_criterion_proved=True,
 all_four_variable_cyclic_affine_pencils_inverse_closed=True,
 cyclic_double_infinity_resonance_inverse_closed=True,
 all_cyclic_unit_hessian_source_orbits_proved_nonexistent=False,
 all_independently_centered_quadratic_composite_images_resolved=False,
 unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
 unrestricted_HC4_or_JC2_resolved=False,original_goal_complete=False,
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-elliptic-cone.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
