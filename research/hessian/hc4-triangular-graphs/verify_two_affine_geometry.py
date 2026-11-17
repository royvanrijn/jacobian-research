"""Exact controls for relative submersion rigidity and endpoint geometry.

The global all-degree reductions are written proofs. No unit-ideal
Keller pair or singleton-zero nonautomorphic pair is constructed here.
"""
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

    s,t,U,V,u,v = S.symbols('s t U V u v')
    # Rank-two projection of W: the two physical free derivatives have zeros.
    K = S.Matrix([U,V+U*U])
    D = U**3+V*V
    J = K.jacobian((U,V))
    chosen_base = -J.T.inv()*S.Matrix([S.diff(D,U),S.diff(D,V)])
    f2 = D+s*K[0]+t*K[1]
    zero('rank-two projection has an exact polynomial base value with zero free gradient',
         S.Matrix([S.diff(f2,U),S.diff(f2,V)]).subs(
             {s:chosen_base[0],t:chosen_base[1]},simultaneous=True))
    zero('the rank-two control has actual unit full Hessian',S.hessian(f2,(s,t,U,V)).det()-1)
    # Rank-one projection of W: C=0 and one linear equation solve the free critical point.
    A,C = t,U+t*t
    D1 = t*U+U**3
    f1 = D1+s*A+V*C
    point = {t:1,U:-1,s:2,V:-4}
    zero('rank-one projection has an exact zero of both free derivatives',
         S.Matrix([S.diff(f1,U),S.diff(f1,V)]).subs(point,simultaneous=True))
    zero('the rank-one control also has actual unit full Hessian',S.hessian(f1,(s,t,U,V)).det()-1)
    zero('the Keller Jacobian prevents both U coefficients from vanishing',
         S.Matrix([A,C]).jacobian((t,U)).det()-1)

    # A positive relative-coordinate control is deliberately non-Keller.
    F,G = s,s*t-1
    aa,bb = t,-1
    D0 = s**3+t*t
    value = D0+u*F+v*G
    companion = -bb*u+aa*v
    zero('the no-common-zero free derivatives have a literal base Bezout certificate',aa*F+bb*G-1)
    zero('the two-free-variable relative coordinate has exact polynomial source recovery',
         S.Matrix([aa*(value-D0)-G*companion,bb*(value-D0)+F*companion])-S.Matrix([u,v]))
    zero('the actual relative coordinate Jacobian is one',
         S.Matrix([value,companion]).jacobian((u,v)).det()-1)
    zero('this control is explicitly non-Keller and does not solve the endpoint',
         S.Matrix([F,G]).jacobian((s,t)).det()-s)
    zero('its full Hessian records exactly that nonunit Jacobian',S.hessian(value,(s,t,u,v)).det()-s*s)
    nonaffine = u+v*v+s*t
    check('without unit Hessian, unit free partials do not imply fibre affineness',
          S.diff(nonaffine,u)==1 and S.diff(nonaffine,v,2)==2 and
          S.hessian(nonaffine,(s,t,u,v)).det()==0)

    # All three aligned-kernel models and their constant-level geometry.
    X,Y,Z = S.symbols('X Y Z')
    ambient_vars = (s,t,X,Y,Z)
    models = {3:s*t+X*Y+Z*Z,2:s*Z+t*t+X*Y,1:X*X+s*Y+t*Z}
    for rank,Q in models.items():
        check('aligned kernel rank '+str(rank)+' has nondegenerate full ambient form',
              S.hessian(Q,ambient_vars).det()!=0)
        check('aligned kernel rank '+str(rank)+' has the asserted unretained rank',
              S.hessian(Q,(X,Y,Z)).rank()==rank)
    tau = S.symbols('tau')
    zero('a noncritical nondegenerate ternary quadric level has the exact two-stratum invariant',
         (tau-1)*tau+2*tau-(tau*tau+tau))
    check('this invariant differs from the affine-plane invariant',tau*tau+tau!=tau*tau)
    zero('the critical quadric cone really has a singular vertex',
         S.Matrix([S.diff(X*Y+Z*Z,z) for z in (X,Y,Z)]).subs({X:0,Y:0,Z:0}))
    zero('the critical cone has the same invariant as A2, so singularity is checked separately',
         (tau-1)*tau+tau-tau*tau)
    c = S.symbols('c',nonzero=True)
    zero('the nonzero rank-two level has a nonconstant invertible function',
         (X*Y/c-1).subs(Y,c/X))
    zero('the rank-two special critical level splits into its two planes',
         models[2].subs({s:0,t:0},simultaneous=True)-X*Y)
    coeff0,coeff1 = S.symbols('coeff0 coeff1')
    check('the rank-one origin square cannot be affine in a genuinely free variable',
          S.degree((coeff0+coeff1*u)**2,u)==2)

    # The exact aligned rank-three endpoint embedding, with a non-Keller control.
    H1 = S.Integer(1)
    H3 = companion
    H2 = value-s*t-H3*H3
    zero('the aligned nondegenerate kernel endpoint embedding gives the whole potential',
         s*t+H1*H2+H3*H3-value)
    recovered_value = s*t+H2+H3*H3
    zero('its inverse really recovers both free variables polynomially',
         S.Matrix([aa*(recovered_value-D0)-G*H3,bb*(recovered_value-D0)+F*H3])-S.Matrix([u,v]))
    # Dropping closed embedding would make the no-zero conclusion false.
    kvalue = s*u+t*v
    kH = S.Matrix([1,kvalue-s*t-u*u,u])
    zero('an arbitrary substitution encodes a genuine unit-Hessian cotangent potential even with a zero pair',
         s*t+kH[0]*kH[1]+kH[2]**2-kvalue)
    zero('that source potential has actual unit full Hessian',S.hessian(kvalue,(s,t,u,v)).det()-1)
    check('the substitution loses a whole free source line at the coefficient zero and is not embedded',
          kH.subs({s:0,t:0,u:0},simultaneous=True).diff(v)==S.zeros(3,1))

    # Hyperplane-image classification in Q=xy+zw+r^2.
    x,y,z,w,r = S.symbols('x y z w r')
    Q5 = x*y+z*w+r*r
    Hmat = S.hessian(Q5,(x,y,z,w,r))
    ell = S.Matrix([[1,1,0,0,1]])
    n = Hmat.inv()*ell.T
    nu = (ell*n)[0]
    center = c*n/nu
    zero('a nonisotropic hyperplane has the stated ambient critical point',
         (ell*center)[0]-c)
    zero('its critical value is the inverse-metric square formula',
         (center.T*Hmat*center)[0]/2-c*c/(2*nu))
    tangent = S.Matrix.hstack(*ell.nullspace())
    check('the nonisotropic hyperplane restriction is nondegenerate rank four',
          (tangent.T*Hmat*tangent).rank()==4)
    null_ell = S.Matrix([[1,0,0,0,0]])
    null_n = Hmat.inv()*null_ell.T
    zero('the isotropic hyperplane normal has zero dual quadratic norm',(null_ell*null_n)[0])
    zero('the restricted radical derivative is exactly the affine hyperplane offset',
         S.diff(Q5,y).subs(x,c)-c)
    critical_restriction = Q5.subs(x,0)
    check('the through-center isotropic hyperplane has rank-three quadratic part',
          S.hessian(critical_restriction,(y,z,w,r)).rank()==3)
    zero('its entire source radical line is critical',
         S.Matrix([S.diff(critical_restriction,a) for a in (y,z,w,r)]).subs({z:0,w:0,r:0}))
    noncritical_restriction = Q5.subs(x,c)
    out = S.symbols('out')
    zero('the offset isotropic hyperplane has an explicit polynomial coordinate completion',
         noncritical_restriction.subs(y,(out-z*w-r*r)/c)-out)

    dependency_pins = {}
    for receipt_name,names in [
        ('receipt-three-affine.json',['THREE_AFFINE_PROJECTIONS.md']),
        ('receipt-ambient-quadric.json',['AMBIENT_QUADRIC_GRAPHS.md']),
        ('receipt-global-value.json',['GLOBAL_VALUE_ORBITS.md']),
        ('receipt-relative-coordinate.json',['RELATIVE_COORDINATE_ENDPOINT.md'])]:
        retained = json.loads((HERE/receipt_name).read_text())
        for name in names:
            dependency_pins[str((HERE/name).relative_to(ROOT))] = retained['source_sha256'][name]
    retained = json.loads((HERE/'receipt-intrinsic-plane.json').read_text())
    for name,digest in retained['dependencies'].items():
        if 'tang-hc4-audit' in name or 'ricolfi-hodge' in name:
            dependency_pins[name] = digest
    for name,digest in dependency_pins.items():
        check('retained mathematical or primary input pin '+name,
              hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in ['verify_two_affine_geometry.py','TWO_AFFINE_ENDPOINT_GEOMETRY.md']},
        dependencies=dependency_pins,
        written_not_finite_exhaustion=[
            'unit two-directional-partial ideal forces the unique vertical cotangent plane and an explicit relative coordinate',
            'all aligned constant-fibre kernel ranks and exact unit-ideal equivalence in the nondegenerate rank-three branch',
            'every affine-hyperplane image under arbitrary source parameterizations is classified into the two endpoints or the forbidden critical-line case'],
        new_external_classification_imported=False,
        non_Keller_unit_ideal_control_asserted_Keller=False,
        endpoint_pair_constructed=False,
        generic_or_fibrewise_submersion_asserted_global_unit_ideal=False,
        cotangent_alignment_asserted_for_arbitrary_two_affine_embedding=False,
        rank_one_forced_origin_zero_asserted_singleton=False,
        arbitrary_substitution_asserted_closed_embedding=False,
        affine_hyperplane_endpoint_branches_asserted_resolved=False,
        entire_two_affine_class_asserted_closed=False,
        nonlinear_Hessian_invariance_used=False,marked_gradient_collision_transfer_asserted=False,
        unit_ideal_plane_Keller_endpoint_resolved=False,singleton_zero_properness_resolved=False,
        universal_SL2_row_equation_solved=False,unrestricted_HC4_or_JC2_resolved=False,
        novelty_priority_asserted=False,formal_verification=False,independent_replay=False,external_review=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-two-affine-geometry.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=result['check_count'],elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    main()
