#!/usr/bin/env python3
"""Exact algebra supporting MORDELL_RAY_CASSELS_TATE_PAIRING.md.

Python standard library only. This checks norm and Kummer identities, finite
field symbols, and matrix reconstruction; it does not prove the cohomological
comparison theorem or Chebotarev. No rank routine or floating point is used.
"""

from fractions import Fraction as F
import json
from math import isqrt
from pathlib import Path


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def neg(a):
    return -a[0], -a[1]


def mul(a, b):
    return a[0] * b[0] + 6 * a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def power(a, n):
    assert n >= 0
    out = (F(1), F(0))
    while n:
        if n & 1:
            out = mul(out, a)
        a = mul(a, a)
        n //= 2
    return out


def scale(a, c):
    return a[0] * c, a[1] * c


def norm(a):
    return a[0] ** 2 - 6 * a[1] ** 2


def verify_torsion_geometry():
    """Polynomial identities after factoring out powers of rho, rho^3=R.

    Coefficients lie in Q(s,d), s^2=6, d^2=-3. This checks the general
    S+T= sigma(T) relation and tangent slope, independently of any curve rank.
    """
    zero, one = (F(0), F(0)), (F(1), F(0))

    def badd(a, b):
        return add(a[0], b[0]), add(a[1], b[1])

    def bscale(a, c):
        return scale(a[0], c), scale(a[1], c)

    def bmul(a, b):
        return add(mul(a[0], b[0]), scale(mul(a[1], b[1]), -3)), add(
            mul(a[0], b[1]), mul(a[1], b[0]))

    s, d = ((F(0), F(1)), zero), (zero, one)
    z = bscale(badd((neg(one), zero), d), F(1, 2))
    slope = bscale(bmul(s, z), -1)  # S-to-T slope divided by rho.
    xsum = badd(bmul(slope, slope), (scale(one, 6), zero))
    assert xsum == bscale(z, -6)
    ysum = badd(bmul(slope, bscale(xsum, -1)), bscale(s, -3))
    assert ysum == bscale(bmul(s, d), 3)
    tangent = bscale(bmul(s, d), -1)
    assert bmul(tangent, bscale(bmul(s, d), 6)) == (scale(one, 108), zero)
    assert bmul(ysum, ysum) == (scale(one, -162), zero)


def quadratic(a):
    return tuple(map(F, a))


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def allowed(p):
    return is_prime(p) and p % 3 == 2 and pow(-2 % p, (p - 1) // 2, p) == p - 1


def chi9(n):
    """Additive cubic character modulo nine, with chi9(2)=1."""
    return {1: 0, 8: 0, 2: 1, 7: 1, 4: 2, 5: 2}[n % 9]


def pair_symbol(q, a, coefficients, norm_root, zeta_twist=0,
                negative_s=False, inverse_zeta=False):
    """Regular-prime residue formula for a symmetrized cubic norm solution.

    Inputs are exact pairs in Q(sqrt(6)); coefficients give xi=A+B*r+C*r^2.
    This finite-field implementation requires all inputs to be units at q.
    General exceptional-prime local normalization is a written obligation.
    """
    assert allowed(q)

    def plus(a, b):
        return (a[0] + b[0]) % q, (a[1] + b[1]) % q

    def times(a, b):
        return (a[0] * b[0] - a[1] * b[1]) % q, (
            a[0] * b[1] + a[1] * b[0] - a[1] * b[1]
        ) % q

    def fpow(a, n):
        out = (1, 0)
        while n:
            if n & 1:
                out = times(out, a)
            a = times(a, a)
            n //= 2
        return out

    def scalar(c):
        c = F(c)
        return c.numerator * pow(c.denominator, -1, q) % q, 0

    def smul(a, c):
        return times(a, scalar(c))

    s = next(t for t in range(q) if t * t % q == 6 % q)
    if negative_s:
        s = -s % q
    z = (0, 1)
    if inverse_zeta:
        z = fpow(z, 2)
    a = quadratic(a)
    a_local = plus(scalar(a[0]), scalar(a[1] * s))
    assert a_local != (0, 0)
    r = fpow(a_local, pow(3, -1, q - 1))
    rbar = smul(fpow(r, q * q - 2), norm_root)

    def xi(t, rr):
        out = (0, 0)
        for i, (A, B) in enumerate(map(quadratic, coefficients)):
            out = plus(out, smul(fpow(rr, i), A + B * t))
        return out

    def theta(t, rr):
        return times(xi(t, times(z, rr)), fpow(xi(t, times(fpow(z, 2), rr)), 2))

    eta = times(times(theta(s, r), theta(-s, rbar)), fpow(z, zeta_twist))
    assert eta != (0, 0)
    residue = fpow(eta, (q * q - 1) // 3)
    return next(j for j in range(3) if fpow(z, j) == residue)


def central_symbol(q, negative_s=False, inverse_zeta=False):
    """The original, uncorrected 5-23 symbol; retained public replay helper."""
    assert q not in (5, 23)
    return pair_symbol(q, (-35, 15),
                       [('-8/3', '14/3'), ('13/3', '1/3'), ('-64/15', '-7/5')],
                       -5, negative_s=negative_s, inverse_zeta=inverse_zeta)


def rank(rows):
    a = [[x % 3 for x in row] for row in rows]
    if not a:
        return 0
    r = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        a[r] = [(x * pow(a[r][c], -1, 3)) % 3 for x in a[r]]
        for i in range(len(a)):
            if i != r:
                u = a[i][c]
                a[i] = [(x - u * y) % 3 for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def reconstruct_matrix(n, radical_vectors, entry, extra_entries=()):
    """Exhaust all alternating forms; retain the specified block and radical."""
    from itertools import product

    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    solutions = []
    for entries in product(range(3), repeat=len(pairs)):
        m = [[0] * n for _ in range(n)]
        for (i, j), x in zip(pairs, entries):
            m[i][j], m[j][i] = x, -x % 3
        fixed = list(extra_entries)
        if entry is not None:
            fixed.append((0, 1, entry))
        if any(m[i][j] != value for i, j, value in fixed):
            continue
        if all(all(sum(x * y for x, y in zip(row, v)) % 3 == 0 for row in m)
               for v in radical_vectors):
            solutions.append(m)
    assert len(solutions) == 1, len(solutions)
    return solutions[0]


def main():
    verify_torsion_geometry()
    data = json.loads(Path(__file__).with_name('mordell_ray_pairing_exact.json').read_text())
    reps = {int(p): quadratic(a) for p, a in data['kummer_basis'].items()}
    for p, a in reps.items():
        if p == 0:
            assert a == (5, 2) and norm(a) == 1
            continue
        assert allowed(p)
        assert norm(a) == -p ** 3
        assert a[0].denominator == a[1].denominator == 1
        assert a[0] % 9 in (1, 8) and a[1] % 3 == 0

    a, b = reps[5], reps[23]
    A, B, C = map(quadratic, data['norm_solution_coefficients'])
    relative_norm = add(add(power(A, 3), mul(power(B, 3), a)),
                        add(mul(power(C, 3), power(a, 2)),
                            scale(mul(mul(mul(A, B), C), a), -3)))
    assert relative_norm == b
    assert add(relative_norm, (F(1), F(0))) != b  # Corrupt norm control.

    a = reps[0]
    A, B, C = map(quadratic, data['wild_norm_solution_coefficients'])
    wild_norm = add(add(power(A, 3), mul(power(B, 3), a)),
                    add(mul(power(C, 3), power(a, 2)),
                        scale(mul(mul(mul(A, B), C), a), -3)))
    assert wild_norm == reps[5]
    for m in (1, 2, 4, 5, 7, 8):
        for n in (1, 2, 4, 5, 7, 8):
            assert chi9(m * n) == (chi9(m) + chi9(n)) % 3

    symbols = {}
    for p, expected in data['central_symbols'].items():
        q = int(p)
        values = {central_symbol(q, ns, iz) for ns in (False, True) for iz in (False, True)}
        assert values == {expected}, (q, values)
        symbols[q] = expected
        corrected = {pair_symbol(q, reps[5], data['norm_solution_coefficients'], -5,
                                 zeta_twist=1, negative_s=ns, inverse_zeta=iz)
                     for ns in (False, True) for iz in (False, True)}
        assert corrected == {(expected + chi9(q)) % 3}
    assert symbols[503] == 0

    wild_symbols = {}
    for p, expected in data['wild_symbols'].items():
        q = int(p)
        values = {pair_symbol(q, reps[0], data['wild_norm_solution_coefficients'], 1,
                              negative_s=ns, inverse_zeta=iz)
                  for ns in (False, True) for iz in (False, True)}
        assert values == {expected}, (q, values)
        wild_symbols[q] = expected

    # Three rational classes at R=115*173 determine the previously unknown
    # ramified character correction uniquely. Other curves are new controls.
    correction = next(e for e in range(3) if (symbols[173] + e * chi9(173)) % 3 == 0)
    assert correction == data['tame_5_23_correction'] == 1

    def block(case, method):
        if method == 'tame_5_23':
            assert case['R'] % 115 == 0
            return (sum(symbols[p] + correction * chi9(p)
                        for p in case['primes'] if p not in (5, 23))
                    - symbols[503]) % 3
        assert method == 'wild_0_5'
        assert case['R'] % 5 == 0 and case['R'] % 9 in (1, 8)
        return (sum(wild_symbols[p] for p in case['primes'] if p != 5)
                - wild_symbols[47]) % 3

    witness_count = 0
    for case in data['curves']:
        R, primes = case['R'], case['primes']
        assert R == __import__('math').prod(primes)
        basis = case.get('basis', primes)
        epsilon = int(R % 9 in (1, 8))
        assert set(basis) == set(primes) | ({0} if epsilon else set())
        vectors = []
        for w in case['rational_kummer_witnesses']:
            x, y = map(F, w['point'])
            assert y * y == x ** 3 + 54 * R ** 2
            v = w['vector']
            assert len(v) == len(basis)
            rhs = power(quadratic(w['cube_root']), 3)
            for p, e in zip(basis, v):
                rhs = mul(rhs, power(reps[p], e))
            assert rhs == (y, F(-3 * R))
            vectors.append(v)
            witness_count += 1
        assert rank(vectors) == len(vectors)
        if case['known_block_from_symbol']:
            entry = block(case, case.get('block_formula', 'tame_5_23'))
        else:
            assert len(vectors) == len(basis)
            entry = None
        extra = [(basis.index(item['pair'][0]), basis.index(item['pair'][1]),
                  block(case, item['formula'])) for item in case.get('additional_blocks', [])]
        m = reconstruct_matrix(len(basis), vectors, entry, extra)
        assert m == case['normalized_pairing_matrix']
        assert rank(m) == case['pairing_rank']
        assert len(basis) - rank(m) == len(vectors)
        print(f"R={R}: exact normalized matrix {m}; rank {rank(m)}; "
              f"rational Kummer dimension {len(vectors)}")
    print(f'PASS: general torsion geometry, two exact norm identities, four embeddings per local symbol, '
          f'{witness_count} rational Kummer identities, and {len(data["curves"])} complete matrices. '
          'Global comparison and Chebotarev remain written proof inputs.')


if __name__ == '__main__':
    main()
