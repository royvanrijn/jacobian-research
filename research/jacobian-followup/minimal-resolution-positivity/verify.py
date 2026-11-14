"""Exact regressions for the written positivity proof and model transport.

No broad search, missing-artifact reconstruction, or geometric proof by finite
sampling is performed. The symbolic identities below are all-parameter checks.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

import sympy as S

from positivity import Model, REQUIRED_GEOMETRY, evaluate

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
checks = []


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    checks.append(name)


def zero(name, expression):
    check(name, S.cancel(expression) == 0)


def ruled_pairing(a, b, c, d, n):
    return -n*a*c+a*d+b*c


def symbolic_checks():
    n, b, N, gamma, rho = S.symbols('n b N gamma rho', integer=True)
    k = (-2, -n-2)
    h = (1, b)
    r = (1, 3*b-n-2)
    p = lambda x, y: ruled_pairing(*x, *y, n)
    zero('ruled H square', p(h, h)-(2*b-n))
    zero('ruled determinant ramification class', r[1]-(k[1]+3*h[1]))
    zero('ruled ramification genus is zero', 1+(p(r, r)+p(k, r))/2)
    zero('nonample ramification fibre sum', r[1].subs(b, n)-(2*n-2))
    # A one-point ramification divisor would require e=2N-1>N for N>1.
    zero('one-point ramification excess', (2*N-2)-(N-1)-(N-1))
    check('one-point ramification impossible polynomial inequality',
          S.Poly((2*N-2)-(N-1), N).as_expr() == N-1)
    # R^2 and Q^2 from K^2=10-rho, H^2=N, K.H=2gamma-2-N.
    kk, kh = 10-rho, 2*gamma-2-N
    rr = S.expand(kk+6*kh+9*N)
    qq = S.expand(4*kk+12*kh+9*N)
    pa = S.expand(1+(2*kk+9*kh+9*N)/2)
    zero('general ramification arithmetic genus', pa-(2+9*gamma-rho))
    zero('general adjoint square', qq-(16+24*gamma-4*rho-3*N))
    zero('general ramification square', rr-(3*N+12*gamma-2-rho))
    # Restriction determinant identity with arbitrary symbolic tangential jets.
    u, z = S.symbols('u z')
    c = S.symbols('c0:15')
    v = sum(c[5*i+j]*u**i*z**j for i in range(3) for j in range(5))
    det = S.Matrix([[S.diff(u, u), S.diff(u, z)],
                    [S.diff(v, u), S.diff(v, z)]]).det()
    zero('smooth reduced inverse-line ramification restriction',
         det.subs(u, 0)-S.diff(v.subs(u, 0), z))
    # Genus change when a surviving curve meets the contracted exceptional v times.
    s, kc, v0 = S.symbols('s kc v0', integer=True)
    zero('non-SNC arithmetic genus transport',
         ((s+v0*v0)+(kc-v0))/2+1-((s+kc)/2+1)-v0*(v0-1)/2)


def numerical_controls():
    # F0 with the usual section+fibre A2 boundary. These are numerical
    # candidates, not actual Keller morphisms. Existing basic pole tests
    # accept them, whereas the new ruled/genus tests reject them.
    sys.path.insert(0, str(RESEARCH/'archive/non-elliptic/plane-jc/cas'))
    from intrinsic_a2_boundary import IntrinsicA2Boundary, audit_a2_boundary, audit_keller_pole_vector
    results = []
    for label, poles, required in (
        ('ruled_degree_four', (1, 2), 'Q_BOUNDARY_INTERSECTION:F:-1'),
        ('genus_one_degree_eight', (2, 2), 'GENERAL_LINE_GENUS_LT_TWO'),
    ):
        names, matrix = ('S', 'F'), [[0, 1], [1, 0]]
        old = audit_keller_pole_vector(
            audit_a2_boundary(IntrinsicA2Boundary(names, S.Matrix(matrix), (0, 0))),
            poles, require_nonproper=False)
        check(label+' passes the retained basic numerical pole test', old.passes)
        model = Model.from_boundary(names, matrix, (0, 0), poles)
        obs = model.positivity_obstructions()
        check(label+' new obstruction', required in obs)
        if label == 'genus_one_degree_eight':
            check(label+' passes positivity alone', not model.positivity_obstructions(False))
        results.append(dict(name=label, geometric_realization=False, old_basic_test_passes=old.passes,
                            invariants=model.invariants, new_obstructions=obs))

    # Blow up a point on the P2 boundary. Raw Q.E=-2 is universal and must
    # disappear before testing the relatively minimal model. Degree16
    # makes the numerical positivity tests on the P2 class pass;
    # this still does not realize an etale A2 source (the rank-one proof
    # in PROOF.md excludes it).
    blown = Model.from_boundary(('L_strict', 'E'), [[0, 1], [1, -1]], (0, 0), (4, 4))
    q = 2*blown.canonical_intersections+3*blown.hyperplane_intersections
    check('raw exceptional Q intersection is minus two', q[1] == -2)
    minimal, history = blown.minimalize()
    check('remove exactly the vertical exceptional', [x['curve'] for x in history] == ['E'])
    check('minimal P2 lattice recovered', minimal.intersection == S.Matrix([[1]]))
    check('minimal numerical positivity passes', not minimal.positivity_obstructions(False))
    check('rank-one whole-source corollary is additional',
          minimal.positivity_obstructions() == ['MINIMAL_PICARD_RANK_LT_TWO'])
    results.append(dict(name='extra_blowup_negative_control', raw_Q_dot_E=q[1],
                        blowdowns=history, minimal_invariants=minimal.invariants,
                        realization_certified=False))
    # A curve of intersection two acquires one arithmetic genus under
    # blowdown; treating it as a smooth rational curve would be wrong.
    # This lattice is a test of the local transformation identity only.
    zero('double intersection creates arithmetic genus one',
         (S.Integer(2)**2-S.Integer(2))/2-1)
    for bad in ('complete_whole_A2_boundary', 'actual_projective_morphism',
                'etale_on_whole_A2', 'target_line_missed'):
        geom = {key: True for key in REQUIRED_GEOMETRY}
        geom[bad] = False
        check('fail closed before matrix reading: '+bad,
              evaluate(dict(geometry=geom))['status'] == 'NOT_APPLICABLE')
    degree_one = Model.from_boundary(('L',), [[1]], (0,), (1,))
    check('degree-one automorphism outside scope',
          degree_one.positivity_obstructions() == ['OUTSIDE_SCOPE: geometric degree one'])
    check('degree-one evaluator does not report exclusion',
          evaluate(dict(geometry={key: True for key in REQUIRED_GEOMETRY},
                        names=['L'], matrix=[[1]], genera=[0], poles=[1]))['status'] == 'OUTSIDE_SCOPE')
    # Smooth blowdowns must not be replaced by contracting a -2 curve.
    try:
        Model.from_boundary(('S', 'F'), [[-2, 1], [1, 0]], (0, 0), (1, 2)).blow_down(0)
    except ValueError:
        check('minus-two contraction rejected', True)
    else:
        check('minus-two contraction rejected', False)
    return results


def r7_window():
    N, B = 6, 10
    gamma, s = S.symbols('gamma s', integer=True)
    zero('r7 Riemann-Hurwitz identity', (-2*N+B+N-s)-(4-s))
    candidates = [(g, B-N+2-2*g) for g in range(2, 1+(B-N+1)//2)
                  if B-N+2-2*g >= 1]
    check('r7 genus and residue sum uniquely forced', candidates == [(2, 2)])
    # These are all partitions of s=2 into positive f_i, with sum e_i*f_i=6.
    types = []
    for f in range(1, 3):
        for e in range(1, N+1):
            if f == 2 and e*f == N:
                types.append(((e, f),))
    for e1 in range(1, N+1):
        for e2 in range(e1, N+1):
            if e1+e2 == N:
                types.append(((e1, 1), (e2, 1)))
    check('four complete generic infinity decompositions',
          types == [((3, 2),), ((1, 1), (5, 1)), ((2, 1), (4, 1)), ((3, 1), (3, 1))])
    # Derive both rational bounds before enumerating their four integers.
    lower = S.ceiling(10-S.Rational(16, 6))
    upper = S.ceiling(S.Rational(46, 4))-1
    check('Hodge lower rank', lower == 8)
    check('strict positivity upper rank', upper == 11)
    rows = [dict(rho=r, canonical_square=10-r, adjoint_square=46-4*r,
                 ramification_square=40-r, ramification_arithmetic_genus=20-r)
            for r in range(int(lower), int(upper)+1)]
    check('r7 numerical rows', [tuple(row.values()) for row in rows] ==
          [(8, 2, 14, 32, 12), (9, 1, 10, 31, 11),
           (10, 0, 6, 30, 10), (11, -1, 2, 29, 9)])
    counts = [dict(rho=r, omitted_unramified=t, infinity_primes=len(infinity),
                   contracted=r-1-t-len(infinity), infinity=infinity)
              for r in range(8, 12) for t in (2, 3) for infinity in types]
    check('r7 contracted boundary lower bound', min(x['contracted'] for x in counts) == 2)
    check('r7 contracted boundary upper bound', max(x['contracted'] for x in counts) == 7)
    # Degree of a generic line on each retained parametrization.
    data = json.loads((HERE.parent/'jc2-normal-r7/classification.json').read_text())
    t, aa, bb, cc = S.symbols('t aa bb cc')
    check('exactly two classified nonimmersed inputs', len(data['cases']) == 2)
    for i, case in enumerate(data['cases'], 1):
        p, q = S.sympify(case['p'], locals={'t': t}), S.sympify(case['q'], locals={'t': t})
        check(f'r7 class{i} parametrization degrees', (S.degree(p, t), S.degree(q, t)) == (6, 10))
        check(f'r7 class{i} generic line intersection degree', S.degree(aa*p+bb*q+cc, t) == 10)
        # This is an exact normalization gcd, not reliance on the stored boolean.
        u = S.symbols('u')
        common = S.Poly(S.gcd(p.subs(t, u)-p, q.subs(t, u)-q), u, t).monic().as_expr()
        zero(f'r7 class{i} birational normalization gcd', common-(u-t))
    return dict(conditional_hypotheses='PF2R7BND1: degree6, sole ramified affine prime (e,f)=(2,1), one of two nonimmersed r7 images',
                affine_ramification_degree=B, general_line_genus=2,
                infinity_residue_sum=2, infinity_types=types, rows=rows,
                complete_count_options=counts, no_case_exclusion_claimed=True)


def application_status():
    # This list is a scope audit, not a fabricated completion of the packets.
    quartic = evaluate(dict(geometry={'etale_on_whole_A2': True, 'target_line_missed': True},
                            reason='Quartic Cox packets do not supply a complete projective boundary and H.'))
    f2 = evaluate(dict(geometry={'etale_on_whole_A2': True, 'target_line_missed': True},
                       reason='F2 (75,125) retains local/partial boundary data; global source attachment remains open.'))
    # Construct only the cheap, retained numerical package; no lower-band
    # elimination or missing generated artifact is requested.
    from log_boundary_compiler import frontier_72_108_case2_boundary_package_audit
    package = frontier_72_108_case2_boundary_package_audit()
    old = package.pole_audit
    boundary = old.boundary_audit.boundary
    diag = Model.from_boundary(boundary.names, boundary.intersection_matrix,
                               boundary.genera, old.pole_vector)
    reduced, history = diag.minimalize()
    diagnostic = dict(original_invariants=diag.invariants, minimal_invariants=reduced.invariants,
                      blowdowns=history, hypothetical_Keller_obstructions=reduced.positivity_obstructions(),
                      actual_Keller_application=False)
    ineligible = evaluate(dict(geometry={**{key: True for key in REQUIRED_GEOMETRY},
                                        'etale_on_whole_A2': False},
                              reason='Retained transformed (72,108) chart has Jacobian X^2.'))
    check('quartic data remain not applicable', quartic['status'] == 'NOT_APPLICABLE')
    check('live F2 data remain not applicable', f2['status'] == 'NOT_APPLICABLE')
    check('transformed Poisson chart rejected', ineligible['status'] == 'NOT_APPLICABLE')
    return dict(quartic=quartic, F2_75_125=f2,
                archived_72_108=dict(eligibility=ineligible, numerical_diagnostic=diagnostic))


def main():
    start = time.monotonic()
    manifest = json.loads((HERE/'sources/manifest.json').read_text())
    for rec in manifest['sources']:
        data = (HERE/'sources'/rec['file']).read_bytes()
        check('source bytes '+rec['file'], len(data) == rec['bytes'])
        check('source hash '+rec['file'], hashlib.sha256(data).hexdigest() == rec['sha256'])
    symbolic_checks()
    controls = numerical_controls()
    r7 = r7_window()
    applications = application_status()
    result = dict(status='PASS', checks=checks, check_count=len(checks),
                  seconds=round(time.monotonic()-start, 3), controls=controls,
                  applications=applications, geometric_proof_is_written=True,
                  r7_window=r7,
                  formal_verification=False, independent_external_review=False,
                  checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  implementation_sha256=hashlib.sha256((HERE/'positivity.py').read_bytes()).hexdigest(),
                  proof_sha256=hashlib.sha256((HERE/'PROOF.md').read_bytes()).hexdigest())
    result['r7_proof_sha256'] = hashlib.sha256((HERE/'R7_WINDOW.md').read_bytes()).hexdigest()
    (HERE/'receipt.json').write_text(json.dumps(result, indent=2, default=str)+'\n')
    print(json.dumps(dict(status=result['status'], checks=len(checks), seconds=result['seconds'],
                          applications={key: value.get('status', value.get('eligibility', {}).get('status'))
                                        for key, value in applications.items()}), indent=2))


if __name__ == '__main__':
    main()
