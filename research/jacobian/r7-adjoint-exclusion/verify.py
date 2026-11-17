"""Exact algebraic endpoints for the written adjoint-fibration proof.

The geometric reduction and the cited Borisov assertions are not outputs
of this checker. All discovery certificates are small univariate Bezout
identities; Fraction replay does not call a polynomial-gcd implementation.
"""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import time

import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
checks = []


def check(name, value):
    if not value:
        raise AssertionError(name)
    checks.append(name)


def zero(name, value):
    check(name, S.cancel(value) == 0)


def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def add(a, b):
    result = [Fraction(0)]*max(len(a), len(b))
    for i, c in enumerate(a): result[i] += c
    for i, c in enumerate(b): result[i] += c
    return trim(result)


def mul(a, b):
    if not a or not b: return []
    result = [Fraction(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): result[i+j] += x*y
    return trim(result)


def frac(poly):
    return [Fraction(x) for x in poly]


def replay_unit(name, a, b, u, v):
    check(name, add(mul(frac(a), frac(u)), mul(frac(b), frac(v))) == [Fraction(1)])


def poly(coefficients, variable):
    return sum(S.Rational(c)*variable**i for i, c in enumerate(coefficients))


def source_checks():
    source = HERE.parent/'jc2-boundary-attachment/sources/borisov-stein-author.pdf'
    check('Borisov pinned source hash', hashlib.sha256(source.read_bytes()).hexdigest() ==
          'dc63edc476dcd8f2b9b35329bc3bbb415bfb15496fb8fa384ff22cfae24271f4')
    ledger = json.loads((ROOT/'MATH_STATUS.json').read_text())
    byid = {x['id']: x for x in ledger['entries']}
    for cid in ('JC2POS1', 'JC2POS2', 'JC2POS3', 'PF2R7CUSP1', 'PF2R7BND1'):
        check('retained dependency '+cid, byid[cid]['state'] == 'proved')


def geometry_regressions():
    rho, m, v = S.symbols('rho m v', integer=True)
    zero('adjoint square in genus two degree six', (10-rho)+2*(-4)+6-(8-rho))
    zero('canonical adjoint pairing', (10-rho)+(-4)-(6-rho))
    zero('movable Hodge upper bound', S.Rational(2)**2/6-S.Rational(2, 3))
    zero('degree-one fibre nef contradiction', 2*(-2)+3*1-(-1))
    check('Hodge and moving-pencil rank endpoint', [r for r in range(8, 12) if 6-r >= -2] == [8])
    zero('type1 raw adjoint intersection', 2*(-1)+3*1-1)
    zero('smooth blowdown adjoint transport', (2*(S.Symbol('kc')-v)+3*m)-(2*S.Symbol('kc')+3*m-2*v))
    check('first positive blowdown contact violates nefness', 1-2*1 < 0)
    # Filter the complete earlier list using the newly inspected f_i=1 theorem.
    old = [((3, 2),), ((1, 1), (5, 1)), ((2, 1), (4, 1)), ((3, 1), (3, 1))]
    surviving = [z for z in old if all(f == 1 for e, f in z)]
    check('three remaining infinity index pairs', surviving == old[1:])
    for epsilon in (0, 1):
        t = 3-epsilon
        c = 8-(1+t+2)
        check(f'cusp epsilon{epsilon} complete boundary count', c == 2+epsilon)
    # Minimal resolution of a*b=y^3: two crepant -2 curves.
    matrix = S.Matrix([[-2, 1], [1, -2]])
    check('affine A2 resolution determinant', matrix.det() == 3)
    check('affine A2 crepant canonical intersections', [-2-matrix[i, i] for i in range(2)] == [0, 0])
    check('epsilon0 needs more exceptional curves than available', 2+1 > 2)
    horizontal = []
    # Exhaust all positive integer (pole, fibre-degree) summands of 2.
    for pole in range(1, 3):
        for degree in range(1, 3):
            if pole*degree == 2: horizontal.append(((pole, degree),))
    horizontal.append(((1, 1), (1, 1)))
    injective = [z for z in horizontal if sum(degree for pole, degree in z) == 1]
    check('injective conic leaves one pole-two section', injective == [((2, 1),)])
    zero('horizontal exceptional section square', -2-1-(-3))
    check('double-line parameter degree', S.Rational(6, 2) == 3)
    # Homogeneous conic tangent to Z=0 at [0:1:0].
    X, Y, Z, a, b, c, z = S.symbols('X Y Z a b c z')
    conic = a*X**2+b*X*Z+Y*Z+c*Z**2
    zero('conic intersection at infinity is double', conic.subs(Z, 0)-a*X**2)
    zero('conic tangent at the fixed point', S.diff(conic, X).subs({X:0, Y:1, Z:0}))
    check('conic smooth at its infinity point', S.diff(conic, Z).subs({X:0, Y:1, Z:0}) == 1)
    # The tangency/Hermite relation eliminates all three parabola coefficients.
    P, V, Q, W, st, sv, aa, bb, cc = S.symbols('P V Q W st sv aa bb cc')
    right = (aa*P**2+bb*P+cc)-(aa*V**2+bb*V+cc)
    zero('quadratic Hermite contact identity', (P-V)*((2*aa*P+bb)+(2*aa*V+bb))-2*right)
    zero('nonaffine involution leading obstruction', S.Rational(10, 6)-2+S.Rational(1, 3))
    return dict(rank=8, remaining_infinity_pairs=surviving,
                remaining_cusp_epsilon=1, omitted_unramified_primes=2,
                exceptional_curves_at_A=3, horizontal_types=horizontal,
                injective_conic_horizontal_type=injective[0])


def algebraic_endpoints():
    t, h, x, A = S.symbols('t h x A')
    source = json.loads((HERE.parent/'jc2-normal-r7/classification.json').read_text())['cases']
    certs = json.loads((HERE/'certificates.json').read_text())
    check('both exact coefficient classes covered', len(certs) == len(source) == 2)
    output = []
    # Positive control: a genuine rational two-contact parabola family.
    # x=t^3,y=t^4 lies on y=2*x^2/(3*z)+z^2/3 with z=t^2.
    z = S.symbols('z')
    family = 2*x**2/(3*z)+z**2/3
    zero('positive two-contact family incidence', family.subs({x:t**3,z:t**2})-t**4)
    zero('positive two-contact family critical condition', S.diff(family,z).subs({x:t**3,z:t**2}))
    p0,q0=t**3,t**4
    e0=S.expand((p0-p0.subs(t,-t))*(S.diff(q0,t)*S.diff(p0,t).subs(t,-t)
                    +S.diff(q0,t).subs(t,-t)*S.diff(p0,t))
                    -2*(q0-q0.subs(t,-t))*S.diff(p0,t)*S.diff(p0,t).subs(t,-t))
    zero('positive involution satisfies Hermite identity',e0)
    check('positive family has the correct field degree three',
          max(S.degree(S.fraction(S.cancel(family))[0],z),
              S.degree(S.fraction(S.cancel(family))[1],z))==3)
    for i, (row, cert) in enumerate(zip(source, certs), 1):
        p, q = S.sympify(row['p']), S.sympify(row['q'])
        check(f'class{i} certificate parametrization matches', cert['p'] == row['p'] and cert['q'] == row['q'])
        dp, dq = S.diff(p, t), S.diff(q, t)
        sq = cert['squarefree']
        zero(f'class{i} squarefree input derivative', poly(sq['dp'], t)-dp)
        zero(f'class{i} squarefree input second derivative', poly(sq['ddp'], t)-S.diff(dp, t))
        replay_unit(f'class{i} independent squarefree Bezout replay', sq['dp'], sq['ddp'], sq['a'], sq['b'])
        sl = cert['slope']
        den, num = poly(sl['den'], t), poly(sl['num'], t)
        zero(f'class{i} reduced slope denominator', (t-1)*den-dp)
        zero(f'class{i} reduced slope numerator', (t-1)*num-dq)
        replay_unit(f'class{i} independent reduced-slope Bezout replay', sl['den'], sl['num'], sl['a'], sl['b'])
        numerator = S.Poly(S.expand(A*den*(x-p)**2+den*q+num*(x-p)), t)
        check(f'class{i} nonzero constant quadratic coefficient degree16', numerator.degree() == 16 and numerator.nth(16) == 6*A)
        check(f'class{i} zero quadratic coefficient degree14', S.degree(numerator.as_expr().subs(A, 0), t) == 14)
        # Modulo den the numerator is independent of A. Thus the lack of
        # cancellation holds for every A in C, not only generic A.
        zero(f'class{i} numerator reduction independent of A', S.rem(numerator.as_expr()-num*(x-p), den, t))
        res = S.resultant(den, x-p, t)
        check(f'class{i} generic target x prevents denominator cancellation', S.degree(res, x) == 4 and res != 0)
        v = 2*h-t
        eh = S.Poly(S.expand((p-p.subs(t, v))*(dq*dp.subs(t, v)+dq.subs(t, v)*dp)
                            -2*(q-q.subs(t, v))*dp*dp.subs(t, v)), t)
        re = cert['reflection']
        zero(f'class{i} t17 coefficient matches certificate', eh.nth(17)-poly(re['c17'], h))
        zero(f'class{i} t11 coefficient matches certificate', eh.nth(11)-poly(re['c11'], h))
        zero(f'class{i} universal reflection leading term', eh.nth(17)-3840*h**3)
        check(f'class{i} centered reflection nonzero', eh.nth(11).subs(h, 0) ==
              (-S.Rational(1280, 243) if i == 1 else -S.Rational(2048, 6075)))
        replay_unit(f'class{i} independent all-center reflection unit', re['c17'], re['c11'], re['bezout17'], re['bezout11'])
        # A tampered certificate must fail the independent engine.
        bad = list(re['bezout11']);bad[0] = str(Fraction(bad[0])+1)
        check(f'class{i} tampered Bezout negative control',
              add(mul(frac(re['c17']), frac(re['bezout17'])), mul(frac(re['c11']), frac(bad))) != [Fraction(1)])
        output.append(dict(case=i, constant_parabola_degrees=[14,16],
                           mobius_involutions='NONE', independent_fraction_replay='PASS',
                           reflection_at_zero=str(eh.nth(11).subs(h, 0))))
    return output


def main():
    start = time.monotonic()
    source_checks()
    geometric_regressions = geometry_regressions()
    endpoints = algebraic_endpoints()
    inputs = ['PROOF.md','SOURCE_AUDIT.md','REVIEW.md','README.md','OBJECTIVE.md',
              'certificates.json','verify.py']
    digests = {f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in inputs}
    result = dict(status='PASS', checks=checks, check_count=len(checks),
                  seconds=round(time.monotonic()-start,3),
                  geometric_regressions=geometric_regressions,
                  algebraic_endpoints=endpoints, input_sha256=digests,
                  geometry_is_written_proof=True,
                  borisov_finality_is_external_input=True,
                  checker_proves_external_finality=False,
                  formal_verification=False, external_review=False)
    (HERE/'receipt.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
    print(json.dumps(dict(status=result['status'],checks=len(checks),seconds=result['seconds'],
                          rank=8, endpoints=endpoints),indent=2))


if __name__ == '__main__':
    main()
