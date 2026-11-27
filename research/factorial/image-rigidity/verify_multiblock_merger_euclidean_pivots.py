#!/usr/bin/env python3
"""Exact, parameter-uniform first-four-pivot certificates; no SFC closure.

SymPy 1.14.0. The symbolic coefficient checks range over unrestricted
parameters by polynomial identities, rather than interpolation or sampling.
The separate rational checks test implementation and normalization.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from math import factorial
from pathlib import Path

import sympy as sp


def bcoef(m):
    return sp.prod(m-i for i in range(1, 4))/(2*m)


def ccoef(m):
    return sp.prod(m-i for i in range(1, 6))/(6*m**2)


def dcoef(m):
    return sp.prod(m-i for i in range(1, 8))/(24*m**3)


def inverse_gamma_cumulants():
    T = sp.Symbol('T')
    raw = [sp.Integer(1)]
    for k in range(1, 5):
        raw.append(1/sp.prod(T-j for j in range(k)))
    actual = [raw[1], raw[2]-raw[1]**2,
              raw[3]-3*raw[2]*raw[1]+2*raw[1]**3,
              raw[4]-4*raw[3]*raw[1]-3*raw[2]**2
              +12*raw[2]*raw[1]**2-6*raw[1]**4]
    expected = [1/T, 1/(T**2*(T-1)), 4/(T**3*(T-1)*(T-2)),
                6*(5*T-6)/(T**4*(T-1)**2*(T-2)*(T-3))]
    assert all(sp.cancel(a-b) == 0 for a, b in zip(actual, expected))


def positive_certificate(expr, variables, substitutions):
    numerator, denominator = sp.fraction(sp.factor(expr))
    polys = [sp.Poly(poly.subs(substitutions), *variables)
             for poly in (numerator, denominator)]
    out = {}
    for label, poly in zip(('numerator', 'denominator'), polys):
        assert all(c > 0 for c in poly.coeffs()), (label, poly)
        assert poly.TC() > 0
        terms = [[list(monomial), str(coefficient)]
                 for monomial, coefficient in poly.terms()]
        digest = hashlib.sha256(json.dumps(terms, separators=(',', ':')).encode()).hexdigest()
        out[label] = {'terms': len(terms), 'degrees': list(poly.degree_list()),
                      'minimum_coefficient': str(min(poly.coeffs())),
                      'constant': str(poly.TC()), 'terms_sha256': digest}
    return out


def third_certificates():
    n, b, c = sp.symbols('n b c')
    x, y, z = sp.symbols('x y z')
    N = n+1
    alpha = sp.factor(bcoef(N)-bcoef(n))
    beta = sp.factor(ccoef(N)-ccoef(n))

    def delta(t):
        return sp.factor(bcoef(N)/(N*t-1)-bcoef(n)/(n*t-1))

    def epsilon(t):
        return sp.factor(ccoef(N)/(N*t-1)-ccoef(n)/(n*t-1))

    def zeta(t):
        return sp.factor(ccoef(N)/((N*t-1)*(N*t-2))
                         - ccoef(n)/((n*t-1)*(n*t-2)))

    A = sp.factor((n-1)*alpha-alpha**2-bcoef(n)+beta)

    def B(t):
        return sp.factor((n-1-2*alpha)*delta(t)-bcoef(n)/(n*t-1)+3*epsilon(t))

    H, J, U, K, V = sp.symbols('H J U K V')
    ds = bcoef(N)*U-bcoef(n)*J
    qs = (n-1-2*alpha)*ds-bcoef(n)*J+3*(ccoef(N)*U-ccoef(n)*J)
    ts = ccoef(N)*V-ccoef(n)*K
    mu4 = n*H**2-alpha*H**2-ds
    mu6 = beta*H**3+3*H*(ccoef(N)*U-ccoef(n)*J)+4*ts \
          +n*H*mu4-bcoef(N)*(H**2+U)*H
    assert sp.cancel(H*mu6-mu4**2-(A*H**4+H**2*qs-ds**2+4*H*ts)) == 0

    expressions = {'1111': 24*A, '211': 12*A+2*B(b),
                   '22': 6*A+B(b)+B(c)-2*delta(b)*delta(c),
                   '31': 4*A+2*B(b)+4*zeta(b),
                   '4': A+B(b)-delta(b)**2+4*zeta(b)}
    answer = {key: positive_certificate(value, (x, y, z),
                                        {n: x+3, b: y+1, c: z+1})
              for key, value in expressions.items()}
    # The bound for the second pivot uses 0 <= delta(b) <= delta(1).
    for value in (delta(b), delta(1)-delta(b)):
        num, den = sp.fraction(sp.factor(value))
        for poly in (num, den):
            p = sp.Poly(poly.subs({n: x+2, b: y+1}), x, y)
            assert all(v >= 0 for v in p.coeffs())
        assert den.subs({n: 2, b: 1}) > 0
    return answer


def partitions(total, largest=None):
    if total == 0:
        yield ()
    else:
        for a in range(min(total, largest or total), 0, -1):
            for rest in partitions(total-a, a):
                yield (a,)+rest


def fourth_certificates():
    n = sp.Symbol('n')
    N = n+1
    H, J, K, L, U, V, W = sp.symbols('H J K L U V W')
    bs = sp.symbols('b0:6')
    xs = sp.symbols('x0:4')
    Bn, BN = bcoef(n)*(H**2+J), bcoef(N)*(H**2+U)
    Cn = ccoef(n)*(H**3+3*H*J+4*K)
    CN = ccoef(N)*(H**3+3*H*U+4*V)
    Dn = dcoef(n)*(H**4+6*H**2*J+3*J**2+16*H*K+6*L)
    DN = dcoef(N)*(H**4+6*H**2*U+3*U**2+16*H*V+6*W)
    mu4 = n*H**2-(BN-Bn)
    mu6 = CN-Cn+n*H*mu4-BN*H
    mu8 = -(DN-Dn)+n*H*mu6-BN*mu4+CN*H
    polynomial = sp.Poly(sp.expand((mu4-H**2)*mu8-mu6**2
                                   +2*H*mu4*mu6-mu4**3), H, J, K, L, U, V, W)

    def coefficient(part):
        r, zero = len(part), (0,)*len(part)
        polys = []
        for degree, m, kind in ((1, n, 0), (2, n, 1), (3, n, 2), (4, n, 3),
                                (2, N, 1), (3, N, 2), (4, N, 3)):
            poly = {}
            for i, a in enumerate(part):
                if a < degree:
                    continue
                bi = bs[i]
                value = 1 if kind == 0 else 1/(m*bi-1)
                if kind == 2:
                    value /= m*bi-2
                if kind == 3:
                    value = (5*m*bi-6)/((m*bi-1)**2*(m*bi-2)*(m*bi-3))
                key = list(zero)
                key[i] = degree
                poly[tuple(key)] = value
            polys.append(poly)

        def multiply(left, right):
            out = {}
            for ka, va in left.items():
                for kb, vb in right.items():
                    key = tuple(a+b for a, b in zip(ka, kb))
                    if all(a <= b for a, b in zip(key, part)):
                        out[key] = out.get(key, 0)+va*vb
            return out

        result = 0
        for powers, scalar in polynomial.terms():
            poly = {zero: sp.factor(scalar)}
            for base, power in zip(polys, powers):
                for _ in range(power):
                    poly = multiply(poly, base)
            result += poly.get(part, 0)
        return sp.factor(result)

    answer = {}
    for part in partitions(6):
        substitutions = {n: xs[0]+4}
        repeated = [i for i, a in enumerate(part) if a >= 2]
        substitutions.update({bs[i]: xs[j+1]+1 for j, i in enumerate(repeated)})
        answer[''.join(map(str, part))] = positive_certificate(
            coefficient(part), xs, substitutions)
    assert len(answer) == 11
    return answer


def exact_polynomial(m, blocks, slots=None):
    top = m//2
    values = [Q(1)]+[Q(0)]*top
    for i, (b, h) in enumerate(blocks):
        T = m*b if slots is None else slots[i]
        row = [Q(1)]
        for k in range(1, top+1):
            row.append(row[-1]*h/(k*(T-k+1)))
        values = [sum((values[k-j]*row[j] for j in range(k+1)), Q(0))
                  for k in range(top+1)]
    return sp.Poly.from_dict({(m-2*k,): (-1)**k*factorial(m)//factorial(m-2*k)*v
                             for k, v in enumerate(values)}, sp.Symbol('u'))


def rational_controls():
    controls = 0
    for blocks in ([(1, Q(1))], [(3, Q(7, 11)), (1, Q(5, 2))],
                   [(Q(7, 3), Q(3, 5)), (Q(13, 4), Q(7, 9)), (1, Q(0))],
                   [(1, Q(2)), (2, Q(3)), (4, Q(5)), (7, Q(11))]):
        H = sum((h/b for b, h in blocks), Q(0))
        for n in range(1, 17):
            lower, upper = exact_polynomial(n, blocks), exact_polynomial(n+1, blocks)
            fixed_slots = exact_polynomial(n, blocks, slots=[(n+1)*b for b, h in blocks])
            assert sp.expand(upper.diff().as_expr()/(n+1)-fixed_slots.as_expr()) == 0
            moments = [sp.Rational(1)]
            for k in range(1, 5):
                value = lower.nth(n-2*k) if n >= 2*k else 0
                for i in range(1, k+1):
                    if n+1 >= 2*i:
                        value -= upper.nth(n+1-2*i)*moments[k-i]
                moments.append(value)
            h, t, v, w = moments[1:]
            aa, bb = t-h**2, h*v-t**2
            cc = (t-h**2)*w-v**2+2*h*t*v-t**3
            expected = [h]
            if n >= 2:
                expected.append(aa/h)
            if n >= 3:
                expected.append(bb/(h*aa))
            if n >= 4:
                expected.append(h*cc/(aa*bb))
            left, right = upper, lower
            for step in range(1, min(4, n)+1):
                remainder = -left.rem(right)
                assert remainder.degree() == n-step
                assert remainder.LC() > 0
                assert remainder.LC() == expected[step-1]
                if step == 1:
                    assert remainder.LC() == H
                if step == 2:
                    assert remainder.LC() >= H
                left, right = right, remainder.monic()
            controls += 1
    # A real-rooted adjacent pair with the wrong order must be rejected.
    u = sp.Symbol('u')
    bad_upper, bad_lower = sp.Poly(u**3-u, u), sp.Poly(u**2-4, u)
    assert (-bad_upper.rem(bad_lower)).LC() < 0
    return controls


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    inverse_gamma_cumulants()
    result = {'scope': 'First four Euclidean pivots of multiblock Gaussian merger limits; '
                        'all starts, arbitrary block count, b_i>=1 and h_i>=0 with sum h_i>0. '
                        'This does not prove all-order interlacing or unrestricted SFC.',
              'third_pivot': third_certificates(),
              'fourth_pivot': fourth_certificates(),
              'rational_control_pairs': rational_controls(),
              'sympy_version': sp.__version__,
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'outcome': 'PASS', 'independent_replay': False,
              'formal_verification': False, 'external_review': False}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('PASS: 16 unrestricted symbolic coefficient certificates, '
          f'{result["rational_control_pairs"]} exact rational pair controls, negative control.')


if __name__ == '__main__':
    main()
