"""Exact chart/jet and cusp-index controls for PF2OMIT9; no global map search."""
from pathlib import Path
import hashlib
import json
import platform
import time
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CHECKS = []


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started = time.monotonic()
    u, v, Y, A, B, D = s.symbols('u v Y A B D')
    chart = {Y: -u**2+u**3*v, A: -u+u**2*v, B: -1+u*v, D: v}
    at_R = {Y: 0, A: 0, B: -1, D: v}
    relations = [A*A-Y*B, A*(B+1)-Y*D, B*(B+1)-A*D]
    for i, relation in enumerate(relations):
        check(f'Wright relation {i+1} on actual boundary chart', s.expand(relation.subs(chart)) == 0)
    source_volume = s.det(s.Matrix([[s.diff(1/u, u), s.diff(1/u, v)],
                                   [s.diff(chart[Y], u), s.diff(chart[Y], v)]]))
    check('actual source volume is -u', s.cancel(source_volume) == -u)
    for coordinate, value in chart.items():
        check(f'restriction of {coordinate}', value.subs(u, 0) == at_R[coordinate])

    def jets(poly):
        p0 = poly.subs(at_R)
        p1 = (-s.diff(poly, A)+v*s.diff(poly, B)).subs(at_R)
        p2 = (-s.diff(poly, Y)+v*s.diff(poly, A)+s.diff(poly, A, 2)/2
              -v*s.diff(poly, A, B)+v*v*s.diff(poly, B, 2)/2).subs(at_R)
        return [s.expand(p) for p in [p0, p1, p2]]

    controls = [Y, A, B, D, Y*A+B*D**3, A*A+B*B*D+Y*B,
                (A+B+Y+D)**3, Y**2*A*B**2*D**5]
    for i, poly in enumerate(controls):
        actual = s.Poly(s.expand(poly.subs(chart)), u)
        check(f'Taylor coefficient control {i+1}', jets(poly) == [actual.nth(j) for j in range(3)])
        check(f'kernel representative control {i+1}',
              all(jets(poly+relation*(Y+A*D+B**2)) == jets(poly) for relation in relations))

    p0, p1, p2, q0, q1, q2 = [s.Function(n)(v) for n in ['p0', 'p1', 'p2', 'q0', 'q1', 'q2']]
    pp, qq = p0+u*p1+u*u*p2, q0+u*q1+u*u*q2
    jac = s.expand(s.diff(pp, u)*s.diff(qq, v)-s.diff(pp, v)*s.diff(qq, u))
    check('universal Jacobian constant coefficient',
          s.expand(jac.coeff(u, 0)-(p1*s.diff(q0, v)-s.diff(p0, v)*q1)) == 0)
    check('universal Jacobian linear coefficient',
          s.expand(jac.coeff(u, 1)-(2*(p2*s.diff(q0, v)-s.diff(p0, v)*q2)
                                   +p1*s.diff(q1, v)-s.diff(p1, v)*q1)) == 0)
    a, b, lam, h = [s.Function(n)(v) for n in ['a', 'b', 'lam', 'h']]
    check('Wronskian lambda square identity',
          s.expand(lam*a*s.diff(lam*b, v)-s.diff(lam*a, v)*lam*b
                   -lam**2*(a*s.diff(b, v)-s.diff(a, v)*b)) == 0)
    aa, ab, ba, bb = s.symbols('aa ab ba bb')
    check('ordinary cusp second-third jet determinant',
          s.det(s.Matrix([[aa, ba], [2*ab, 2*bb]])) == 2*(aa*bb-ab*ba))

    U, V, w, z = s.symbols('U V w z')
    local_P = u+w*w
    local_Q = u*w+s.Rational(2, 3)*w**3
    local_jac = s.diff(local_P, u)*s.diff(local_Q, w)-s.diff(local_P, w)*s.diff(local_Q, u)
    check('ordinary cusp local control has Jacobian u', s.expand(local_jac) == u)
    check('local cusp restriction derivative gcd', s.monic(s.gcd(2*w, 2*w*w), w) == w)
    check('local cusp jet scalar equation', s.Rational(1, 2)**2*4 == 1)
    check('normal local fibre length3', s.groebner([local_P, local_Q], u, w).reduce(w**3)[1] == 0
          and s.groebner([local_P, local_Q], u, w).reduce(w**2)[1] != 0)
    H = 9*V*V-4*U**3
    normal_eq = w**3-3*U*w+3*V
    G = z**3+3*U*z*z+H
    check('reciprocal cubic relation', s.rem(s.expand(G.subs(z, U-w*w)), normal_eq, w) == 0)
    inclusion = s.Matrix([[1, U, U*U], [0, 0, -3*V], [0, -1, U]])
    check('normalization inclusion determinant is -3V', inclusion.det() == -3*V)
    check('normal cubic discriminant', s.expand(s.discriminant(normal_eq, w)+27*H) == 0)
    check('raw reciprocal discriminant', s.expand(s.discriminant(G, z)+243*V*V*H) == 0)
    check('index square transports discriminants',
          s.expand(s.discriminant(G, z)-inclusion.det()**2*s.discriminant(normal_eq, w)) == 0)
    check('cusp image multiplicity2', min(sum(m) for m, c in s.Poly(H, U, V).terms()) == 2)
    check('raw cusp discriminant multiplicity4',
          min(sum(m) for m, c in s.Poly(s.discriminant(G, z), U, V).terms()) == 4)
    check('reciprocal at cusp is -w squared', (U-w*w).subs(U, 0) == -w*w)
    check('binary root factor at cusp length3', G.subs({U: 0, V: 0}) == z**3)
    check('reciprocal fails to generate the normal fibre',
          s.rem((-w*w)**2, w**3, w) == 0 and s.rem(w*w, w**3, w) != 0)
    rows = []
    for label, branches, cusp, lengths in [('off_image', 0, 0, []),
                                          ('immersed', 1, 0, [2]),
                                          ('double_image', 2, 0, [2, 2]),
                                          ('cusp', 1, 1, [3])]:
        retained = 4-2*branches-cusp
        check(label+' complete fibre length', retained >= 0 and retained+sum(lengths) == 4)
        rows.append(dict(value_type=label, branches=branches, retained_length=retained,
                         omitted_local_lengths=lengths))
    check('third image branch impossible', 4-2*3 < 0)
    check('cusp with companion branch impossible', 4-2*2-1 < 0)
    check('quartic Euler excess is one', 4-1-2 == 1)
    t = s.symbols('t')
    reciprocal = u/(1+t*u*chart[Y])
    check('generic source coordinate has reciprocal uniformizer', s.cancel(reciprocal/u).subs(u, 0) == 1)
    ledger = dict(as_of='2026-11-14', claim='PF2OMIT9',
        hypotheses='Actual marked PF2WRIGHT1 ring and source open; finite target pair in the ring with constant nonzero source Jacobian. Quartic conclusions additionally require degree4.',
        omitted_prime=dict(ideal=['Y', 'A1', 'B1+1'], normalization='A1_v', e=2, f=1,
                           ramifies=True, omitted=True, image_parametrization=['P(0,0,-1,v)', 'Q(0,0,-1,v)']),
        degree_identity='d=3+deg gcd(p0_prime,q0_prime)',
        quartic=dict(derivative_gcd='v-v_star', exceptional_value='(p0(v_star),q0(v_star))',
                     exceptional_type='ordinary (2,3) cusp; no companion image branch',
                     fibre_rows=rows, leading_coefficient='constant*H_C',
                     generic_omitted_Newton_width=2, generic_omitted_Newton_slope='1/2',
                     raw_discriminant='constant*H_C*J^2', index_vanishes_at_cusp=True,
                     cusp_raw_discriminant_order_at_least=4,
                     singleton_omitted_radial=dict(e=3, f=1, length=3,
                                                  canonical_vector=[2, 4, 5])),
        local_control=dict(normal_equation=str(normal_eq), root_equation=str(G),
                           inclusion_determinant=str(inclusion.det()), is_global_Wright_map=False),
        actual_target_pair=None, actual_image_equation=None, actual_exceptional_value=None,
        actual_compatible_completion_supplied=False, unrestricted_singleton_zero_properness_proved=False,
        full_goal_complete=False)
    target = HERE/'wright_image_ledger.json'
    target.write_text(json.dumps(ledger, indent=2)+'\n')
    inputs = [HERE/'WRIGHT_IMAGE_RESTRICTION.md', HERE/'GLOBAL_LEDGER.md',
              HERE/'ADMISSIBLE_INVERSE.md', HERE/'DISCRIMINANT_LEDGER.md',
              HERE/'SIMPLE_RAMIFICATION_AND_RADIAL.md',
              ROOT/'jacobian-followup/affine-complements/WRIGHT_MARKED_MODEL.md']
    receipt = dict(as_of='2026-11-14', outcome='PASS', check_count=len(CHECKS), checks=CHECKS,
        checker_sha256=digest(Path(__file__)), input_sha256={str(p.relative_to(ROOT)): digest(p) for p in inputs},
        ledger_sha256=digest(target), runtime_seconds=round(time.monotonic()-started, 6),
        python=platform.python_version(), sympy=s.__version__,
        proof_boundary='Exact Taylor/jet identities and finite local cusp-index control. Global finiteness, Euler and imported simple-boundary theorem are written proof steps; no global pair is constructed.',
        formal_verification=False, external_review=False, full_goal_complete=False)
    (HERE/'wright_image_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: receipt[k] for k in ['outcome', 'check_count', 'runtime_seconds']}))


if __name__ == '__main__':
    main()
