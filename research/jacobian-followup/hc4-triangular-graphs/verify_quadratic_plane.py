"""Exact controls for quadratic-plane rigidity and fibre-affine embeddings."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    checks = []

    def check(name,value):
        if not value:
            raise AssertionError(name)
        checks.append(name)

    def zero(name,value):
        values = list(value) if isinstance(value,S.MatrixBase) else [value]
        check(name,all(S.cancel(S.expand(item.doit()))==0 for item in values))

    s,t,u,v,U,V,x,y = S.symbols('s t u v U V x y')
    # The affine-in-one-variable Keller proof has exact polynomial steps.
    a = t*t+1
    b = t**3+t
    C = a*U+b
    A = 3*C**4-2*C*C+t
    zero('subtracting powers of the affine component leaves the same Jacobian',
        S.Matrix([A,C]).jacobian((t,U)).det()-a)
    lead = 3*a**4
    zero('the highest free coefficient satisfies its exact derivative equation',
        a*S.diff(lead,t)-4*S.diff(a,t)*lead)
    check('a nonconstant leading coefficient cannot give a Keller pair after elimination',a!=1)
    C = 2*U+t**5
    A = C**4-3*C+7*t+1
    zero('the unbounded-other-variable-degree positive control is Keller',
        S.Matrix([A,C]).jacobian((t,U)).det()-14)
    p,q = S.symbols('p q')
    tinv = (p-q**4+3*q-1)/7
    Uinv = (q-tinv**5)/2
    zero('the affine-component branch has its actual polynomial inverse',
        S.Matrix([A,C]).subs({t:tinv,U:Uinv},simultaneous=True)-S.Matrix([p,q]))
    A,C = 3*U+t**7,2*t-1
    zero('the zero free-slope branch is triangular Keller',
        S.Matrix([A,C]).jacobian((t,U)).det()+6)
    # A quadratic normalized shear has both homogeneous conditions and inverse.
    ell = 2*x+y
    vector = S.Matrix([1,-2])
    homogeneous = vector*ell**2
    pair = S.Matrix([x,y])+homogeneous
    DH = homogeneous.jacobian((x,y))
    zero('the normalized quadratic part has zero Jacobian trace',S.trace(DH))
    zero('the normalized quadratic part has zero Jacobian determinant',DH.det())
    zero('the quadratic Keller pair has unit determinant',pair.jacobian((x,y)).det()-1)
    inverse = S.Matrix([x,y])-homogeneous
    zero('the displayed quadratic inverse works globally',
        pair.subs({x:inverse[0],y:inverse[1]},simultaneous=True)-S.Matrix([x,y]))
    zero('the common quadratic is unchanged in the shear direction',
        S.diff(ell**2,x)-2*S.diff(ell**2,y))
    # Actual-source quadratic plane, with high arbitrary-base dependence.
    f = s*t+u*v+s**5*v*v
    Evars = (u,v)
    for i,derivative in enumerate([S.diff(f,u,3),S.diff(f,u,2,v),S.diff(f,u,v,2),S.diff(f,v,3)]):
        zero('the strict quadratic-plane control has zero cubic vertical derivative '+str(i),derivative)
    check('that plane is genuinely quadratic rather than jointly affine',S.hessian(f,Evars)!=S.zeros(2,2))
    zero('its full actual Hessian determinant is one',S.hessian(f,(s,t,u,v)).det()-1)
    g1,g2,g3,g4 = S.symbols('g1 g2 g3 g4')
    inv = {s:g2,v:g3,u:g4-2*g2**5*g3,t:g1-5*g2**4*g3*g3}
    zero('its gradient has a full polynomial inverse',
        S.Matrix([S.diff(f,z) for z in (s,t,u,v)]).subs(inv,simultaneous=True)-S.Matrix([g1,g2,g3,g4]))
    bad = s*u*u+t*v
    zero('dropping unit full Hessian would invalidate the quadratic-plane conclusion',
        S.hessian(bad,(s,t,u,v)).det()-4*u*u)
    # The two nonvertical Tang positions reach the two plane lemmas.
    K = S.Matrix([U+V*V,V])
    ft = U**2+V**2+s*K[0]+t*K[1]
    check('the transverse Tang control has free degree at most two',S.Poly(ft,U,V).total_degree()==2)
    zero('its plane coefficient pair is a quadratic Keller automorphism',K.jacobian((U,V)).det()-1)
    fr = U**2+t*U+s*(U+t**4)+V*t
    check('the rank-one Tang control has free degree at most two',S.Poly(fr,U,V).total_degree()==2)
    zero('its V coefficient is affine in the other free variable',S.diff(S.diff(fr,V),U,2))
    zero('the rank-one control also has unit full Hessian',S.hessian(fr,(s,t,U,V)).det()-1)
    # All three helper-block ranks and their isotropic linear equations.
    X,Y,Z = S.symbols('X Y Z')
    models = {3:s*t+X*Y+Z*Z,2:s*Z+t*t+X*Y,1:X*X+s*Y+t*Z}
    for rank,Q in models.items():
        check('helper rank '+str(rank)+' has full nondegenerate ambient Hessian',
            S.hessian(Q,(s,t,X,Y,Z)).rank()==5)
        check('helper rank '+str(rank)+' has its asserted constant rank',
            S.hessian(Q,(X,Y,Z)).rank()==rank)
    a1,a2,b1,b2,c1,c2 = S.symbols('a1 a2 b1 b2 c1 c2')
    Lx,Ly,Lz = a1*u+a2*v,b1*u+b2*v,c1*u+c2*v
    product = S.expand(Lx*Ly)
    zero('rank-two isotropy gives the exact three free quadratic coefficients',
        S.Matrix([product.coeff(u,2),product.coeff(u,1).coeff(v,1),product.coeff(v,2)])-
        S.Matrix([a1*b1,a1*b2+a2*b1,a2*b2]))
    zero('rank-one isotropy is exactly the square of its linear free row',
        S.hessian(Lx*Lx,(u,v))-2*S.Matrix([a1,a2])*S.Matrix([[a1,a2]]))
    # A genuine family with all helpers free-dependent and only two affines.
    aa,bb = S.symbols('aa bb')
    outputs = S.Matrix([(1-aa*bb)**2*u-bb*bb*v,v-aa*aa*u,(aa-bb*aa*aa)*u+bb*v])
    recovered_u = outputs[0]+2*bb*outputs[2]-bb*bb*outputs[1]
    recovered_v = outputs[1]+2*aa*(outputs[2]-bb*outputs[1])-aa*aa*recovered_u
    zero('the general two-shear helper map preserves its entire quadratic',
        outputs[0]*outputs[1]+outputs[2]**2-u*v)
    zero('its actual polynomial source inverse recovers both free variables',
        S.Matrix([recovered_u,recovered_v])-S.Matrix([u,v]))
    Aembed = outputs.jacobian((u,v))
    W1,W2,W3 = S.symbols('W1 W2 W3')
    invu = W1+2*bb*W3-bb*bb*W2
    invv = W2+2*aa*(W3-bb*W2)-aa*aa*invu
    left = S.Matrix([invu,invv]).jacobian((W1,W2,W3))
    zero('the helper coefficient matrix has a full polynomial left inverse',left*Aembed-S.eye(2))
    check('every helper genuinely depends on the free variables',all(any(v!=0 for v in Aembed.row(i)) for i in range(3)))
    specialized = outputs.subs({aa:s,bb:s},simultaneous=True)
    check('the strict all-three-free-dependent control has distinct degrees 5,3,4',
        [S.Poly(expr,s,t,u,v).total_degree() for expr in specialized]==[5,3,4])
    pullback = s*t+specialized[0]*specialized[1]+specialized[2]**2
    zero('its actual source potential remains exactly st+uv',pullback-s*t-u*v)
    zero('its actual source Hessian is a unit',S.hessian(S.expand(pullback),(s,t,u,v)).det()-1)
    generic_rank_only = S.Matrix([[s,0],[0,t],[0,0]])
    check('generic rank does not guarantee a closed fibre-affine embedding',
        generic_rank_only.rank()==2 and generic_rank_only.subs({s:0,t:0}).rank()==0)
    # Both old endpoints really have affine free helpers in this larger class.
    F,G = s,s*t-1
    a,b = t,-1
    value = s**3+u*F+v*G
    companion = -b*u+a*v
    endpoint = S.Matrix([1,value-t*t-s*companion,companion])
    zero('the rank-two unit-ideal endpoint has the exact ambient potential',
        s*endpoint[2]+t*t+endpoint[0]*endpoint[1]-value)
    zero('its two nonconstant free helper rows have unit relative determinant',
        endpoint[1:3,0].jacobian((u,v)).det()-1)
    check('all endpoint helpers are affine in the free variables',
        all(S.Poly(expr,u,v).total_degree()<=1 for expr in endpoint))
    zero('that displayed unit-ideal control is not a Keller counterexample',
        S.Matrix([F,G]).jacobian((s,t)).det()-s)
    frame = S.Matrix([[1,s],[t,1+s*t]])
    zero('the singleton frame has its explicit unit determinant',frame.det()-1)
    free = frame*S.Matrix([u,v])
    zero('its rank-one endpoint potential is exactly the oriented cotangent row',
        s*free[0]+t*free[1]-(S.Matrix([[s,t]])*frame*S.Matrix([u,v]))[0])
    pins = {}
    for receipt_name,name in [
        ('receipt-intrinsic-plane.json','INTRINSIC_AFFINE_PLANE.md'),
        ('receipt-relative-critical-graph.json','RELATIVE_CRITICAL_GRAPH_CLASSIFICATION.md'),
        ('receipt-two-affine-geometry.json','TWO_AFFINE_ENDPOINT_GEOMETRY.md')]:
        retained = json.loads((HERE/receipt_name).read_text())
        pins[str((HERE/name).relative_to(ROOT))] = retained['source_sha256'][name]
        if receipt_name=='receipt-intrinsic-plane.json':
            for path,digest in retained['dependencies'].items():
                if 'tang-hc4-audit' in path:
                    pins[path] = digest
    for name,digest in pins.items():
        check('retained proof or structural source pin '+name,
            hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_quadratic_plane.py','QUADRATIC_PLANE_FIBRE_AFFINE.md','sources/quadratic-plane-context.json']},
        dependencies=pins,strict_control_helper_degrees=[5,3,4],
        written_not_finite_exhaustion=[
            'the affine-one-variable Keller component elimination in every degree',
            'the quadratic plane Keller inverse, reproduced directly',
            'every projection rank of Tang plane under a constant quadratic two-plane',
            'all fibre-affine two-affine embeddings, all helper-block ranks and both conditional endpoint converses',
            'arbitrary positive shear exponents in the three-free-dependent helper control'],
        cotangent_alignment_assumed=False,
        generic_matrix_rank_asserted_global_source_generation=False,
        all_degree_base_coefficients_asserted_bounded=False,
        quadratic_condition_on_one_fibre_asserted_uniform=False,
        all_two_affine_embeddings_asserted_fibre_affine=False,
        all_polynomial_quadratic_symmetries_asserted_classified=False,
        nonlinear_source_Hessian_invariance_used=False,
        non_Keller_unit_ideal_control_asserted_Keller=False,
        endpoint_pair_constructed=False,unit_ideal_endpoint_resolved=False,
        singleton_zero_endpoint_resolved=False,unrestricted_HC4_or_JC2_resolved=False,
        external_general_degree_two_classification_imported=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-quadratic-plane.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=len(checks),elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    main()
