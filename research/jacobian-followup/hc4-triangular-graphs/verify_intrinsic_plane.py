"""Exact controls for intrinsic affine planes and pure-trace directions.

The all-degree uniqueness, Tang existence and Hodge recognition proofs
are written arguments, separate from these finite identities.
"""
import argparse
from collections import Counter
import hashlib
import itertools
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

    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    def zero(name, value):
        values = list(value) if isinstance(value,S.MatrixBase) else [value]
        check(name,all(S.cancel(item.doit())==0 for item in values))

    s,U,P,Q = S.symbols('s U P Q')
    variables = (s,U,P,Q)
    F = S.Function('F')(s)
    g1,g0,d1,d0 = [S.Function(name)(s) for name in ('g1','g0','d1','d0')]
    f = P*F+Q*(g1*U+g0)+d1*U+d0
    hh = S.hessian(f,variables)
    zero('intersecting affine planes have a line in common',
         S.Matrix([hh[2,2],hh[2,3],hh[3,3],hh[2,1],hh[1,1]]))
    zero('intersecting-plane full Hessian is the square of Fprime times g1',
         hh.det()-(S.diff(F,s)*g1)**2)
    zero('the residual triangular pair has precisely that Jacobian',
         S.Matrix([F,g1*U+g0]).jacobian((s,U)).det()-S.diff(F,s)*g1)

    alpha,beta,gamma = S.symbols('alpha beta gamma',nonzero=True)
    concrete = P*(alpha*s+beta)+Q*(gamma*U+s**5)+U*(s**4+s)+s**7
    outputs = S.Matrix([S.diff(concrete,z) for z in variables])
    ys,yu,yp,yq = S.symbols('ys yu yp yq')
    sb = (yp-beta)/alpha
    ub = (yq-sb**5)/gamma
    qb = (yu-(sb**4+sb))/gamma
    pb = (ys-5*sb**4*qb-ub*(4*sb**3+1)-7*sb**6)/alpha
    zero('intersecting-plane triangular pair gives an actual full gradient inverse',
         outputs.subs({s:sb,U:ub,P:pb,Q:qb},simultaneous=True)-S.Matrix([ys,yu,yp,yq]))

    a,b,c,d = S.symbols('a b c d')
    C = S.Matrix([[a,b],[c,d]])
    linear = S.Matrix(S.symbols('l0:4'))
    quad = (S.Matrix([[s,U]])*C*S.Matrix([P,Q]))[0]+(linear.T*S.Matrix(variables))[0]
    hquad = S.hessian(quad,variables)
    zero('complementary affine planes give a constant bilinear Hessian',
         hquad-S.BlockMatrix([[S.zeros(2),C],[C.T,S.zeros(2)]]).as_explicit())
    zero('complementary-plane Hessian determinant is the constant block square',
         hquad.det()-C.det()**2)
    targets = S.Matrix([ys,yu,yp,yq])
    zero('complementary-plane gradient has its exact affine inverse',
         hquad*(hquad.inv()*(targets-linear))+linear-targets)

    misleading = P*Q+s*U
    zero('two separate affine directions do not assert their mixed derivative vanishes',
         S.diff(misleading,P,2)+S.diff(misleading,Q,2))
    check('joint affineness requires the mixed condition as well',S.diff(misleading,P,Q)==1)

    T = S.Matrix([[1,1,0,1],[0,1,1,0],[0,0,1,1],[0,0,0,1]])
    check('oblique source coordinate change is invertible',T.det()==1)
    wbasis = T.inv()[:,2:4]
    normals = T[:2,:]
    zero('oblique jointly affine plane is annihilated by the quotient normals',normals*wbasis)
    check('the annihilator has dimension two',normals.rank()==2 and wbasis.rank()==2)
    check('two independent pure normals determine a two-dimensional common kernel',
          len(normals.nullspace())==2)
    check('a third independent pure normal cannot annihilate that same plane',T[:3,:].rank()==3 and (T[:3,:]*wbasis).rank()==1)
    dependent = normals.col_join(normals[0,:]+normals[1,:])
    check('three pairwise nonparallel normals can still span only two dimensions',
          dependent.rank()==2 and all(dependent[list(pair),:].rank()==2 for pair in itertools.combinations(range(3),2)))

    hu,hv = S.symbols('hu hv')
    t = hu*hv
    comp,genus,punctures,degree,m = S.symbols('comp genus punctures degree m')
    rank_one = t*t+(degree-comp)*t+genus*(hu+hv)+punctures-comp
    poly = S.Poly(rank_one,hu,hv)
    zero('a pure rank-one trace forces both curve obstructions to vanish',
         S.Matrix([poly.coeff_monomial(hu)-genus,poly.coeff_monomial(1)-(punctures-comp)]))
    zero('rank-zero traces have exactly the required nonpositive pure signature',
         (t-m)*t-(t*t-m*t))

    J = S.Function('J')(s,U)
    H = s*U*J
    unit = s*P+U*Q+H
    zero('an unrestricted mixed base polynomial keeps unit cotangent Hessian',
         S.hessian(unit,variables).det()-1)
    zero('the first pure trace has a free coordinate and a product',unit.subs(s,0)-U*Q)
    zero('the second independent pure trace has a free coordinate and a product',unit.subs(U,0)-s*P)
    zero('the third independent pure trace has an actual shear to a product',
         unit.subs(P,0).subs(Q,Q-s*J)-U*Q)
    zero('the fourth coordinate trace has its actual product shear too',
         unit.subs(Q,0).subs(P,P-U*J)-s*P)
    zero('all four zero-offset generic traces have negative pure coefficient',
         (t-1)*t-(t*t-t))
    ua,ub = S.symbols('ua ub')
    linear_addition = ua*P+ub*Q+3*s-2*U+7
    zero('all affine-linear additions preserve the full Hessian',
         S.hessian(unit+linear_addition,variables)-S.hessian(unit,variables))
    zero('the quotient row perturbations shift the plane pair exactly',
         S.Matrix([S.diff(unit-ua*P-ub*Q,P),S.diff(unit-ua*P-ub*Q,Q)])-S.Matrix([s-ua,U-ub]))

    # A genuine automorphic Keller pair is enough for an honest profile control.
    shifted = s-1
    zero('the profile control is a genuine unit-Jacobian plane pair',
         S.Matrix([shifted,U]).jacobian((s,U)).det()-1)
    check('a zero coefficient at offset zero does not imply an empty whole profile',
          S.gcd(shifted.subs(s,0),U)==1 and S.gcd(shifted.subs(s,1),U)==U)
    zero('a profile root is simple because the specialized row has a nonzero U derivative',S.diff(U,U)-1)
    zero('linear perturbations recover the exact translated zero point',
         S.Matrix([shifted-ua,U-ub]).subs({s:1+ua,U:ub},simultaneous=True))

    # Abstract finite sets only: no multi-point Keller origin fibre is asserted.
    point_set = {(0,0),(1,2),(1,3)}
    directions = ((1,0),(0,1),(1,1),(2,-1))
    profiles = []
    for ell in directions:
        profile = Counter(ell[0]*z[0]+ell[1]*z[1] for z in point_set)
        check('finite reduced-set profile preserves total count '+str(ell),sum(profile.values())==len(point_set))
        profiles.append(profile)
    check('one projection retains coincident-point counts',profiles[0]=={0:1,1:2})
    check('an oblique projection retains a different collision',profiles[-1]=={0:2,-1:1})
    candidates = set(itertools.product(profiles[0],profiles[1]))
    recovered = {z for z in candidates if all(profile.get(ell[0]*z[0]+ell[1]*z[1],0)>0
                  for ell,profile in zip(directions,profiles))}
    check('n plus one nonparallel complete profiles recover the finite set',recovered==point_set)
    check('two coordinate profiles alone retain spurious cross-pairings',candidates!=point_set)
    check('all four quotient directions are pairwise nonparallel',
          all(S.Matrix(pair).det()!=0 for pair in itertools.combinations(directions,2)))

    dependencies = {}
    for proof,receipt in [('TWO_TRACE_COMPATIBILITY.md','receipt-two-trace.json'),
                          ('HODGE_TRACE_CLOSURE.md','receipt-hodge-trace.json'),
                          ('CRITICAL_COMPONENTS.md','receipt-critical-components.json')]:
        retained = json.loads((HERE/receipt).read_text())
        digest = hashlib.sha256((HERE/proof).read_bytes()).hexdigest()
        check('retained proof pin '+proof,digest==retained['source_sha256'][proof])
        dependencies[str((HERE/proof).relative_to(ROOT))] = digest
    tang_pins = json.loads((ROOT/'jacobian-followup/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['jacobian-followup/tang-hc4-audit/STRUCTURAL_REVIEW.md',
                 'jacobian-followup/tang-hc4-audit/source-v1.pdf']:
        digest = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input '+name,digest==tang_pins[name])
        dependencies[name] = digest
    manifest = json.loads((HERE/'sources/hodge-trace-import.json').read_text())
    source = HERE/'sources'/manifest['retained_file']
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    check('retained primary Hodge source',digest==manifest['sha256'])
    dependencies[str(source.relative_to(ROOT))] = digest
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in [Path(__file__).name,'INTRINSIC_AFFINE_PLANE.md']},
        dependencies=dependencies,
        written_not_finite_exhaustion=['elementary two-plane uniqueness in both relative positions',
            'pure-direction space equals the unique-plane annihilator under nonautomorphy',
            'three linearly independent pure directions force gradient invertibility',
            'two directions determine the only possible plane and its three polynomial identities',
            'two-direction-only universal closure is equivalent to JC2',
            'all complete parallel profiles count the reduced residual plane fibre',
            'n plus one quotient profiles recover a reduced n-point set',
            'affine-linear additions preserve pure directions under nonautomorphy'],
        elementary_uniqueness_requires_Tang=False,
        plane_existence_under_nonautomorphy_requires_Tang=True,
        positive_Hodge_coefficient_required=False,
        pure_direction_space_dimension_two_asserted_without_nonautomorphy=False,
        pairwise_nonparallel_asserted_sufficient=False,
        two_pure_directions_alone_asserted_to_close_JC2=False,
        multiple_point_controls_asserted_realized_by_Keller_pair=False,
        one_zero_profile_coefficient_asserted_unit_ideal=False,
        full_singleton_profile_asserted_properness=False,
        nonlinear_Hessian_invariance_used=False,
        unit_ideal_plane_Keller_endpoint_resolved=False,singleton_zero_properness_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,
        formal_verification=False,independent_replay=False,external_review=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    rendered = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':
    main()
