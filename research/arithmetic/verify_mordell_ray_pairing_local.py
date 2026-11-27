#!/usr/bin/env python3
"""Replay direct local Mordell-ray pairing certificates, using SageMath 10.9.

Run: sage -python research/arithmetic/verify_mordell_ray_pairing_local.py
The exact norm inputs are checked over Q(sqrt(6)); local operations retain
p-adic precision. No rank, global rational-point, or alternating-matrix reconstruction
routine supplies an entry. Global descent and local reciprocity remain written
mathematical inputs. --generate writes exploratory certificates, not proofs of
the global comparison or of a uniform run-time bound.
"""

import argparse
from itertools import product
import json
from pathlib import Path

from sage.all import (GF, Infinity, PolynomialRing, QQ, Qp, QuadraticField,
                      ZZ, inverse_mod, matrix, prod, vector)
from sage.rings.padics.precision_error import PrecisionError


DEFAULT_INPUT = Path(__file__).with_name('mordell_ray_pairing_local_exact.json')
HILBERT_THREE = [[0, 0, 0, 2], [0, 0, 1, 0],
                 [0, 2, 0, 0], [1, 0, 0, 0]]


def cubic_algebra(field, radicand):
    """Three-coordinate etale algebra; also works when the cubic splits.

    Explicit multiplication avoids generic quotient rings with inexact moduli.
    The inverse uses the exact cubic adjugate identity.
    """
    class Cubic:
        def __init__(self, value):
            if isinstance(value, Cubic):
                self.a = value.a
            elif isinstance(value, (tuple, list)):
                self.a = tuple(field(x) for x in (list(value) + [0] * 3)[:3])
            else:
                self.a = (field(value), field(0), field(0))

        def __add__(self, other):
            return Cubic(tuple(x + y for x, y in zip(self.a, Cubic(other).a)))

        __radd__ = __add__

        def __neg__(self):
            return Cubic(tuple(-x for x in self.a))

        def __sub__(self, other):
            return self + -Cubic(other)

        def __rsub__(self, other):
            return Cubic(other) + -self

        def __mul__(self, other):
            a, b, out = self.a, Cubic(other).a, [field(0)] * 3
            for i in range(3):
                for j in range(3):
                    out[(i + j) % 3] += a[i] * b[j] * (radicand if i + j >= 3 else 1)
            return Cubic(out)

        __rmul__ = __mul__

        def norm(self):
            a, b, c = self.a
            return a**3 + radicand*b**3 + radicand**2*c**3 - 3*radicand*a*b*c

        def inverse(self):
            a, b, c = self.a
            n = self.norm()
            return Cubic(((a*a-radicand*b*c)/n, (radicand*c*c-a*b)/n, (b*b-a*c)/n))

        def __truediv__(self, other):
            return self * Cubic(other).inverse()

        def __pow__(self, n):
            if n < 0:
                return self.inverse()**(-n)
            out, a = Cubic(1), self
            while n:
                if n % 2:
                    out = out * a
                a, n = a*a, n//2
            return out

    return Cubic, Cubic((0, 1, 0))


def symbolic_identities():
    """Universal polynomial checks on the curve, not a sampled specialization."""
    ring = PolynomialRing(QQ, names=('s', 'd', 'r', 'x', 'y'))
    s, d, r, x, y = ring.gens()
    field = ring.fraction_field()
    z, c = (d-1)/2, 3*r**3*s
    alpha, u = y-c, c*d
    sigma = lambda f: field(f)(s, d, z**2*r, x, y)
    sigma2 = lambda f: field(f)(s, d, z*r, x, y)
    g = field(alpha+s*z*r*x)/(x+6*z*r**2)
    tangent = y+u+s*d*r*x
    relations = ring.ideal([s**2-6, d**2+3, y**2-x**3-54*r**6]).groebner_basis()
    identities = (g*sigma(g)*sigma2(g)-alpha,
                  sigma(tangent)*g**3-alpha*tangent,
                  g**2*sigma(g)*(y-u)-alpha*tangent)
    for expression in identities:
        assert ring(field(expression).numerator()).reduce(relations) == 0
    assert ring(field(identities[2]+1).numerator()).reduce(relations) != 0


def local_setup(p, precision, negative_s):
    assert p != 2
    base = Qp(p, precision)
    ring = base['t']
    t = ring.gen()
    split = base(6).is_square()
    if p == 3:
        field = base.extension(t**2-6, names='s')
        s = field.gen()
        root = base(-2).sqrt()
        if root.residue() != 1:
            root = -root
        d = s/root
    elif split:
        s = base(6).sqrt()
        if int(s.residue()) > p//2:
            s = -s
        if base(-3).is_square():
            field, d = base, base(-3).sqrt()
        else:
            field = base.extension(t**2+3, names='d')
            d, s = field.gen(), field(s)
    else:
        field = base.extension(t**2-6, names='s')
        s, d = field.gen(), field(-3).sqrt()
    if negative_s:
        assert split
        s = -s
    degree = 2 if split and not base(-3).is_square() else 1
    return base, field, s, d, (d-1)/2, degree


def cube_root(a):
    field, v = a.parent(), a.valuation()
    if v % 3:
        return None
    pi = field.uniformizer()
    unit = a/pi**v
    if field.prime() == 3:
        sign = 1 if unit.residue() == 1 else -1
        if (unit/sign-1).valuation() < 4:
            return None
        root = sign*((unit/sign).log()/3).exp()
    else:
        ring = field.residue_field()['t']
        roots = (ring.gen()**3-unit.residue()).roots()
        if not roots:
            return None
        root = field.teichmuller(roots[0][0])
        # Newton iteration with a unit derivative; precision failure is an error.
        for _ in range(16):
            root = (2*root+unit/root**2)/3
    root *= pi**(v//3)
    assert (root**3/a-1).valuation() > 12
    return root


def point_from_seed(base, field, s, R, a, seed):
    x = base(QQ(seed['x']))
    square = x**3+54*R**2
    assert square and square.is_square()
    y = square.sqrt()
    unit_residue = lambda q: int((q/base(base.prime())**q.valuation()).residue())
    if unit_residue(y) != seed['y_unit_residue']:
        y = -y
    assert unit_residue(y) == seed['y_unit_residue']
    root = cube_root((field(y)-3*R*s)/a)
    assert root is not None, 'Local seed does not represent the specified Kummer class'
    return field(x), field(y), root


def find_point(base, field, s, R, a, search_bound):
    p = base.prime()
    for scale in range(3):
        for n in range(1, search_bound+1):
            exact_x = QQ(n)/p**(2*scale)
            x = base(exact_x)
            square = x**3+54*R**2
            if not square or not square.is_square():
                continue
            y0 = square.sqrt()
            for y in (y0, -y0):
                root = cube_root((field(y)-3*R*s)/a)
                if root is not None:
                    residue = int((y/base(p)**y.valuation()).residue())
                    return {'x': str(exact_x), 'y_unit_residue': residue}
    raise ValueError(f'Bounded point search exhausted at R={R}, p={p}; no impossibility claim')


def descent_scalar(field, s, d, z, R, coefficients, a, x, y, root):
    algebra, r = cubic_algebra(field, R)
    xi = algebra(coefficients)
    sigma = lambda f: algebra(tuple(c*z**(2*i) for i, c in enumerate(algebra(f).a)))
    sigma2 = lambda f: algebra(tuple(c*z**i for i, c in enumerate(algebra(f).a)))
    alpha, u = y-3*R*s, 3*R*s*d
    point_lift = (r*(s*z**2*x)+alpha)/(r**2*6+x)
    assert (point_lift.norm()/alpha-1).valuation() > 12
    theta = xi*root/point_lift
    assert (theta.norm()-1).valuation() > 12
    for h in (algebra(1), r, r**2):
        gamma = h+theta*sigma(h)+theta*sigma(theta)*sigma2(h)
        n = gamma.norm()
        if n and n.precision_absolute()-n.valuation() > 12:
            break
    else:
        raise PrecisionError('No resolved Hilbert-90 seed; increase precision')
    beta = a/((y-u)*n)
    assert beta.precision_absolute()-beta.valuation() > 12
    tangent = r*(s*d*x)+(y+u)
    quotient = sigma(xi)**2*sigma2(xi)/tangent
    error = quotient/(sigma(gamma)**3*beta)-1
    assert min(q.valuation() if q else Infinity for q in error.a) > 12
    return beta


def cube_test_three(a):
    assert a.precision_absolute()-a.valuation() > 4
    if a.valuation() % 3:
        return False
    unit = a/a.parent().uniformizer()**a.valuation()
    sign = 1 if unit.residue() == 1 else -1
    return (unit/sign-1).valuation() >= 4


def coordinates_three(a, z):
    lam = 1-z
    basis = [lam]+[1-lam**i for i in (1, 2, 3)]
    e0 = int(a.valuation() % 3)
    solutions = [(e0,)+e for e in product(range(3), repeat=3)
                 if cube_test_three(a/prod(b**j for b, j in zip(basis, (e0,)+e)))]
    assert len(solutions) == 1
    return solutions[0]


def hilbert(a, b, z):
    field = a.parent()
    if field.prime() == 3:
        aa = vector(GF(3), coordinates_three(a, z))
        bb = vector(GF(3), coordinates_three(b, z))
        return int(aa*matrix(GF(3), HILBERT_THREE)*bb)
    va, vb, pi = a.valuation(), b.valuation(), field.uniformizer()
    unit = (a/pi**va)**vb/(b/pi**vb)**va
    assert unit.precision_absolute() > 0
    residue = unit.residue()**((field.residue_field().order()-1)//3)
    return next(i for i in range(3) if z.residue()**i == residue)


def verify_global(case):
    R = ZZ(case['R'])
    assert R > 1 and R.is_squarefree()
    primes = R.prime_divisors()
    assert all(p % 3 == 2 and ZZ(-2).kronecker(p) == -1 for p in primes)
    L = QuadraticField(6, 's')
    s = L.gen()
    indices = [row['index'] for row in case['basis']]
    assert len(indices) == len(set(indices))
    assert set(indices) == set(primes) | ({0} if R % 9 in (1, 8) else set())
    for row in case['basis']:
        i = row['index']
        A, B = map(ZZ, row['representative'])
        a = A+B*s
        if i == 0:
            assert (A, B) == (5, 2)
        else:
            assert A % i == B % i == 0 and a.norm() == -i**3
            assert B % 3 == 0 and A % 9 in (1, 8)
        coeff = [L(QQ(c))+L(QQ(e))*s for c, e in row['norm_coefficients']]
        assert len(coeff) == 3
        algebra, _ = cubic_algebra(L, R)
        assert algebra(coeff).norm() == a, 'Invalid exact global norm certificate'
    return indices


def evaluate(case, generate=False, search_bound=500):
    indices = verify_global(case)
    R, precision = case['R'], case['precision']
    records = {(r['index'], r['p'], r['negative_s']): r for r in case.get('local', [])}
    expected_keys, local_records, result = set(), [], []
    for row in case['basis']:
        i, coefficients = row['index'], row['norm_coefficients']
        support = {3} | set(ZZ(R).prime_divisors())
        support |= {int(p) for pair in coefficients for q in pair
                    for p in QQ(q).denominator().prime_divisors() if p != 2}
        entries = [0]*len(indices)
        for p in sorted(support):
            for ns in ((False, True) if Qp(p)(6).is_square() else (False,)):
                key = (i, p, ns)
                expected_keys.add(key)
                if not generate:
                    assert key in records, f'Missing required local place {key}'
                base, field, s, d, z, degree = local_setup(p, precision, ns)
                a = field(row['representative'][0])+row['representative'][1]*s
                coeff = [field(QQ(A))+field(QQ(B))*s for A, B in coefficients]
                seed = (find_point(base, field, s, R, a, search_bound) if generate
                        else records[key]['point'])
                x, y, root = point_from_seed(base, field, s, R, a, seed)
                beta = descent_scalar(field, s, d, z, R, coeff, a, x, y, root)
                values = [hilbert(beta, field(other['representative'][0])+
                                   other['representative'][1]*s, z)
                          for other in case['basis']]
                weight = int(2*inverse_mod(degree, 3) % 3)
                record = {'index': i, 'p': int(p), 'negative_s': ns,
                          'point': seed, 'weight': weight, 'symbols': values}
                if p == 3:
                    record['scalar_coordinates'] = list(coordinates_three(beta, z))
                if not generate:
                    assert record == records[key], (key, record, records[key])
                local_records.append(record)
                entries = [(x+weight*y) % 3 for x, y in zip(entries, values)]
        result.append(entries)
    if not generate:
        assert set(records) == expected_keys, 'Missing or duplicate local place coverage'
        assert len(records) == len(case['local'])
        assert result == case['matrix']
    m = matrix(GF(3), result)
    assert m+m.transpose() == 0 and all(m[i, i] == 0 for i in range(m.nrows()))
    if generate:
        case['local'], case['matrix'] = local_records, result
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=DEFAULT_INPUT)
    parser.add_argument('--generate', action='store_true', help='Populate local certificates from exact norm inputs')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--search-bound', type=int, default=500)
    args = parser.parse_args()
    if args.generate:
        assert args.output is not None, '--generate requires an explicit output path'
    data = json.loads(args.input.read_text())
    symbolic_identities()
    for case in data['cases']:
        answer = evaluate(case, generate=args.generate, search_bound=args.search_bound)
        print(f'R={case["R"]}, basis={[r["index"] for r in case["basis"]]}: direct matrix {answer}', flush=True)
    if args.generate:
        args.output.write_text(json.dumps(data, indent=2)+'\n')
    else:
        # A bad norm must fail before any local computation is accepted.
        import copy
        bad = copy.deepcopy(data['cases'][0])
        bad['basis'][0]['norm_coefficients'][0][0] = str(QQ(bad['basis'][0]['norm_coefficients'][0][0])+1)
        try:
            verify_global(bad)
        except AssertionError:
            pass
        else:
            raise AssertionError('Corrupted norm certificate accepted')
        bad = copy.deepcopy(data['cases'][0])
        auxiliary = next(record for record in bad['local'] if (6*bad['R']) % record['p'])
        bad['local'].remove(auxiliary)
        try:
            evaluate(bad)
        except AssertionError as error:
            assert 'Missing required local place' in str(error)
        else:
            raise AssertionError('Missing auxiliary place accepted')
        print('PASS: universal descent identities, exact norms, every retained local place, '
              'wild coordinates, direct matrices, corrupted-norm and missing-place rejection. '
              'Global pairing comparison remains a written input.')


if __name__ == '__main__':
    main()
