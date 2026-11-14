#!/usr/bin/env python3
"""Exact finite regressions for the written cusp-filtration proofs.

Standard library only. Universal rigidity, separability and dimension arguments
are written proofs in PROOF.md, not consequences of these finite controls.
"""
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from math import comb, gcd
from pathlib import Path
import json
import platform
import time

HERE = Path(__file__).resolve().parent
checks = []


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    checks.append(name)


def clean(poly):
    return {m: F(c) for m, c in poly.items() if c}


def add(*polys):
    out = defaultdict(F)
    for poly in polys:
        for monomial, coefficient in poly.items():
            out[monomial] += coefficient
    return clean(out)


def scale(poly, coefficient):
    return clean({m: c*coefficient for m, c in poly.items()})


def mul(left, right):
    out = defaultdict(F)
    for (i, j), c in left.items():
        for (r, s), d in right.items():
            out[i+r, j+s] += c*d
    return clean(out)


def power(poly, n):
    out = {(0, 0): F(1)}
    while n:
        if n & 1:
            out = mul(out, poly)
        poly = mul(poly, poly)
        n >>= 1
    return out


def derivative(poly, axis):
    out = {}
    for monomial, c in poly.items():
        if monomial[axis]:
            reduced = list(monomial)
            reduced[axis] -= 1
            out[tuple(reduced)] = c*monomial[axis]
    return clean(out)


def bracket(left, right):
    return add(mul(derivative(left, 0), derivative(right, 1)),
               scale(mul(derivative(left, 1), derivative(right, 0)), -1))


def normal(poly, a, b):
    """Exact remainder using x^a=y^b+w; keys are (x,y,w) exponents."""
    out = defaultdict(F)
    for (i, j), c in poly.items():
        quotient, remainder = divmod(i, a)
        for k in range(quotient+1):
            out[remainder, j+b*(quotient-k), k] += c*comb(quotient, k)
    return clean(out)


def from_normal(poly, a, b):
    out = defaultdict(F)
    for (i, j, k), c in poly.items():
        for r in range(k+1):
            out[i+a*r, j+b*(k-r)] += c*comb(k, r)*(-1)**(k-r)
    return clean(out)


def delta(poly, a, b):
    reduced = normal(poly, a, b)
    if not reduced:
        return None
    h = a*b-a-b
    return max(b*i+a*j+h*k for i, j, k in reduced)


def normal_basis(bound, a, b):
    h = a*b-a-b
    return [(i, j, k) for i in range(a)
            for j in range(bound//a+1)
            for k in range(bound//h+1)
            if b*i+a*j+h*k <= bound]


def rank(matrix, columns):
    rows = [[F(c) for c in row] for row in matrix]
    pivot_row = 0
    for column in range(columns):
        selected = next((r for r in range(pivot_row, len(rows))
                         if rows[r][column]), None)
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        pivot = rows[pivot_row][column]
        rows[pivot_row] = [c/pivot for c in rows[pivot_row]]
        for r in range(len(rows)):
            if r != pivot_row and rows[r][column]:
                factor = rows[r][column]
                rows[r] = [x-factor*y for x, y in zip(rows[r], rows[pivot_row])]
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


X = {(1, 0): F(1)}
Y = {(0, 1): F(1)}
ONE = {(0, 0): F(1)}


def derivation_kernel(a, b):
    """All possible Lx,Ly at filtered order <=0; impose the complete Lw bound."""
    h = a*b-a-b
    basis_x, basis_y = normal_basis(b, a, b), normal_basis(a, a, b)
    columns = []
    for monomial in basis_x:
        image = from_normal({monomial: F(1)}, a, b)
        columns.append(normal(scale(mul(power(X, a-1), image), a), a, b))
    for monomial in basis_y:
        image = from_normal({monomial: F(1)}, a, b)
        columns.append(normal(scale(mul(power(Y, b-1), image), -b), a, b))
    high = sorted({m for col in columns for m in col
                   if b*m[0]+a*m[1]+h*m[2] > h})
    matrix = [[col.get(m, 0) for col in columns] for m in high]
    n = len(columns)
    check(f'({a},{b}) complete preserving-derivation kernel has dimension1',
          n-rank(matrix, n) == 1)
    euler = [F(b) if m == (1, 0, 0) else F(0) for m in basis_x]
    euler += [F(a) if m == (0, 1, 0) else F(0) for m in basis_y]
    check(f'({a},{b}) actual Euler vector satisfies every high-degree equation',
          all(sum(c*v for c, v in zip(row, euler)) == 0 for row in matrix))
    return dict(source=[a, b], variables=n, constraints=len(matrix), kernel=1)


def normalized_shape(P, Q, a, b, p, q):
    if (p, q) == (a, b) and set(P) == {(1, 0)} and set(Q) == {(0, 1)}:
        return P[1, 0]**a == Q[0, 1]**b
    if (p, q) == (b, a) and set(P) == {(0, 1)} and set(Q) == {(1, 0)}:
        return P[0, 1]**b == Q[1, 0]**a
    return False


def diagnose(P, Q, a, b, p, q):
    if min(a, b, p, q) < 2 or gcd(a, b) != 1 or gcd(p, q) != 1:
        raise ValueError('Both cusp pairs must be coprime integers >=2')
    J = bracket(P, Q)
    if not J:
        return dict(algebraically_independent=False)
    A, B = delta(P, a, b), delta(Q, a, b)
    C = delta(add(power(P, p), scale(power(Q, q), -1)), a, b)
    j = delta(J, a, b)
    m = max(F(A, q), F(B, p))
    defect = C-(p*q-p-q)*m
    check(f'degree inequality {a,b,p,q,A,B,C,j}', defect >= j)
    if p*A != q*B:
        check(f'strict off-slope inequality {a,b,p,q,A,B}', defect > j)
    check(f'zero defect classification control {a,b,p,q,A,B,C,j}',
          (defect == 0) == normalized_shape(P, Q, a, b, p, q))
    return dict(algebraically_independent=True, source=[a, b], target=[p, q],
                component_degrees=[A, B], cusp_difference_degree=C,
                jacobian_degree=j, slope=str(m), defect=str(defect),
                jacobian_lower_bound_sharp=defect == j)


def graded_bracket(left, right, h):
    return scale(mul({(h+1, 0): F(1)}, bracket(left, right)), -1)


def graded_degree(poly, h):
    return max(i+h*j for i, j in poly)


def semigroup_contains(n, a, b):
    return any(n-a*i >= 0 and (n-a*i) % b == 0 for i in range(n//a+1))


def graded_symbol(a, b, p, q, e, U, V, H):
    h, s = a*b-a-b, p*q-p-q
    label = f'graded {a,b}->{p,q}, e={e}'
    check(label+' lies in actual source cusp ring',
          all(semigroup_contains(i, a, b) for poly in (U, V, H) for i, j in poly))
    check(label+' binomial relation', power(U, p) == power(V, q))
    check(label+' bracket U,V', graded_bracket(U, V, h) == {})
    check(label+' bracket V,H',
          graded_bracket(V, H, h) == scale(power(U, p-1), -p))
    check(label+' bracket U,H',
          graded_bracket(U, H, h) == scale(power(V, q-1), -q))
    check(label+' nonzero bracket certifies dominance in written argument',
          bool(graded_bracket(V, H, h)))
    check(label+' three exact graded degrees',
          [graded_degree(poly, h) for poly in (U, V, H)] == [e*q, e*p, e*s])
    return dict(source=[a, b], target=[p, q], scale=e,
                graded_degrees=[e*q, e*p, e*s],
                polynomial_Keller_lift_with_this_scaling='excluded by written theorem')


def main():
    start = time.monotonic()
    source_path = HERE/'sources/uniform-cusp-filtration.md'
    manifest = json.loads((HERE/'sources/manifest.json').read_text())
    source_bytes = source_path.read_bytes()
    check('pinned external source hash',
          sha256(source_bytes).hexdigest() == manifest['source_sha256'])
    check('external source retains provisional lifecycle',
          'Lifecycle: PROVISIONAL' in source_bytes.decode())
    # Check the author's own sealed mathematical body, independently of our manifest.
    marker = b'<!-- BODY-END -->'
    body_end = source_bytes.index(marker)+len(marker)
    if source_bytes[body_end:body_end+2] == b'\r\n':
        body_end += 2
    elif source_bytes[body_end:body_end+1] == b'\n':
        body_end += 1
    body = source_bytes[:body_end]
    check('external mathematical body byte count', len(body) == 9567)
    check('external mathematical body original seal',
          sha256(body).hexdigest() ==
          '2e8ff8d38a4452fa6bd4ec6f6405eda06fd2b80506a45a3d64e3cf97dcbba4ce')

    kernels = []
    for a in range(2, 9):
        for b in range(2, 9):
            if gcd(a, b) != 1:
                continue
            h = a*b-a-b
            w = add(power(X, a), scale(power(Y, b), -1))
            check(f'({a},{b}) actual cusp degree', delta(w, a, b) == h)
            check(f'({a},{b}) source binomial exact normal form',
                  normal(w, a, b) == {(0, 0, 1): F(1)})
            check(f'({a},{b}) y,cusp bracket',
                  bracket(Y, w) == scale(power(X, a-1), -a))
            check(f'({a},{b}) x,cusp bracket',
                  bracket(X, w) == scale(power(Y, b-1), -b))
            for i in range(2*a+1):
                poly = add({(i, b+1): F(i+1)}, {(i+1, 1): F(-2)})
                check(f'({a},{b}) normal form roundtrip {i}',
                      from_normal(normal(poly, a, b), a, b) == poly)
            kernels.append(derivation_kernel(a, b))

    controls = []
    for a, b in ((2, 3), (2, 5), (3, 4), (3, 5)):
        pairs = [
            ('identity', X, Y),
            ('swap', Y, X),
            ('translation', add(X, ONE), Y),
            ('x shear', add(X, power(Y, 2)), Y),
            ('y shear', X, add(Y, power(X, 2))),
            ('mixed automorphism', add(X, power(Y, 2)),
             add(Y, power(add(X, power(Y, 2)), 2))),
            ('matched scaling', scale(X, 2**b), scale(Y, 2**a)),
            ('unmatched scaling', scale(X, 2), Y),
            ('power n=3', power(X, 3), power(Y, 3)),
        ]
        for p, q in sorted({(a, b), (b, a), (2, 3), (2, 5), (5, 3)}):
            for label, P, Q in pairs:
                result = diagnose(P, Q, a, b, p, q)
                result['control'] = label
                controls.append(result)
        h = a*b-a-b
        for n in range(1, 9):
            result = diagnose(power(X, n), power(Y, n), a, b, a, b)
            check(f'power sharpness {a,b,n}',
                  result['cusp_difference_degree'] == n*h+(a+b)*(n-1)
                  and result['jacobian_degree'] == (a+b)*(n-1)
                  and result['jacobian_lower_bound_sharp'])
    check('dependent-pair guard',
          diagnose(power(X, 2), power(X, 3), 2, 5, 2, 3)
          == dict(algebraically_independent=False))
    translation = diagnose(add(X, ONE), Y, 2, 5, 2, 5)
    check('translation defeats component-degrees-only premise',
          translation['component_degrees'] == [5, 2]
          and translation['cusp_difference_degree'] == 5
          and translation['defect'] == '2')
    cubic = diagnose(power(X, 3), power(Y, 3), 2, 5, 2, 5)
    check('nonconstant Jacobian saturates sharp bound',
          cubic['cusp_difference_degree'] == 23
          and cubic['jacobian_degree'] == 14 and cubic['defect'] == '14')

    symbols = [
        graded_symbol(2, 5, 2, 5, 3, {(15, 0): F(1)}, {(6, 0): F(1)},
                      {(6, 1): F(1, 3)}),
        graded_symbol(2, 5, 2, 3, 5, {(15, 0): F(1)}, {(10, 0): F(1)},
                      {(2, 1): F(1, 5)}),
        graded_symbol(2, 5, 5, 3, 65, {(195, 0): F(1)}, {(325, 0): F(1)},
                      {(452, 1): F(1, 65)}),
    ]
    v = mul(X, power(Y, 5))
    # Two arbitrary nonzero leading-degree controls for the stated carrier class.
    for coefficients in ((1, 2, 3), (5, -1, 7)):
        R = add(scale(ONE, coefficients[0]), scale(v, coefficients[1]),
                scale(power(v, 2), coefficients[2]))
        carrier = mul(mul(X, power(add(v, scale(ONE, -1)), 2)), R)
        check(f'F2 carrier weight calibration {coefficients}',
              [delta(power(carrier, n), 2, 5) for n in (1, 3, 5)]
              == [65, 195, 325])

    f2_envelopes = []
    for cap, band, maximum, vertex, count in (
            (75, 15, 195, (15, 60), 706),
            (125, 25, 325, (25, 100), 1901)):
        support = [(i, j) for i in range(cap+1) for j in range(cap-i+1)
                   if 5*i-j <= band]
        check(f'F2 envelope lattice count {cap}', len(support) == count)
        check(f'F2 envelope unique auxiliary cusp maximum {cap}',
              max(5*i+2*j for i, j in support) == maximum
              and [m for m in support if 5*m[0]+2*m[1] == maximum] == [vertex])
        f2_envelopes.append(dict(total_degree_cap=cap, band_cap=band,
                                 lattice_points=count, maximum=maximum,
                                 unique_vertex=list(vertex)))
    check('F2 scalar normalization',
          F(-5, 9)**3 == F(-125, 729))
    check('F2 equal-slope exact threshold',
          F(195, 3) == F(325, 5) == 65
          and (5*3-5-3)*65 == 455)
    check('F2 source and target cusp types are distinct up to swap',
          {2, 5} != {5, 3})

    receipt = dict(
        as_of='2026-11-14', outcome='PASS', check_count=len(checks),
        checks=checks, preserving_derivation_kernels=kernels,
        polynomial_controls=controls, graded_controls=symbols,
        f2_support_envelopes=f2_envelopes,
        proof_boundary='Universal proofs are written; finite controls are not exhaustive Keller searches.',
        source_commit=manifest['commit'], source_sha256=manifest['source_sha256'],
        files={name: sha256((HERE/name).read_bytes()).hexdigest()
               for name in ('PROOF.md', 'verify.py', 'sources/manifest.json',
                            'sources/uniform-cusp-filtration.md')},
        support_source_files={name: sha256((HERE.parents[1]/name).read_bytes()).hexdigest()
            for name in ('archive/non-elliptic/plane-jc/F2_75_125_DERIVATION.md',
                         'archive/non-elliptic/plane-jc/F2_75_125_CARRIER_WRONSKIAN_CLASSIFIER.md',
                         'archive/non-elliptic/plane-jc/cas/classify_f2_75_125_layers.py')},
        python=platform.python_version(), elapsed_seconds=round(time.monotonic()-start, 4))
    (HERE/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: receipt[k] for k in ('outcome', 'check_count', 'elapsed_seconds')}))


if __name__ == '__main__':
    main()
