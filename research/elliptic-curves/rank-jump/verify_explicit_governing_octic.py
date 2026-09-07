#!/usr/bin/env python3
"""CAS-independent discriminant and Frobenius replay of the three frozen octics."""
import argparse
from fractions import Fraction
from pathlib import Path
import sys
import retrospective as r

INPUT = r.OUT/'rank_jump_explicit_governing_octic_v1.json'
OUTPUT = r.OUT/'rank_jump_explicit_governing_octic_verification_v1.json'


def determinant(a):
    a = [row[:] for row in a]; sign = 1; previous = 1; n = len(a)
    for k in range(n-1):
        if a[k][k] == 0:
            j = next(j for j in range(k+1, n) if a[j][k])
            a[j], a[k] = a[k], a[j]; sign = -sign
        pivot = a[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                value = a[i][j]*pivot-a[i][k]*a[k][j]
                assert value % previous == 0
                a[i][j] = value//previous
            a[i][k] = 0
        previous = pivot
    return sign*a[-1][-1]


def discriminant(h):
    n = len(h)-1; derivative = [i*h[i] for i in range(1, n+1)]
    size = 2*n-1; rows = []
    for coeff, count in ((h, n-1), (derivative, n)):
        for i in range(count):
            row = [0]*size
            row[i:i+len(coeff)] = reversed(coeff)
            rows.append(row)
    return (-1)**(n*(n-1)//2)*determinant(rows)//h[-1]


def finite_polynomials(p):
    def trim(a):
        a = [z % p for z in a]
        while a and not a[-1]: a.pop()
        return a
    def sub(a, b):
        return trim([(a[i] if i < len(a) else 0)-(b[i] if i < len(b) else 0)
                     for i in range(max(len(a), len(b)))])
    def rem(a, b):
        a = trim(a); b = trim(b); assert b
        inverse = pow(b[-1], -1, p)
        while len(a) >= len(b):
            c = a[-1]*inverse % p; offset = len(a)-len(b)
            for i, z in enumerate(b): a[i+offset] = (a[i+offset]-c*z) % p
            a = trim(a)
        return a
    def mul(a, b, modulus):
        out = [0]*(len(a)+len(b)-1)
        for i, v in enumerate(a):
            for j, w in enumerate(b): out[i+j] += v*w
        return rem(out, modulus)
    def power(a, n, modulus):
        ans = [1]
        while n:
            if n & 1: ans = mul(ans, a, modulus)
            a = mul(a, a, modulus); n //= 2
        return ans
    def gcd(a, b):
        while b: a, b = b, rem(a, b)
        return trim([z*pow(a[-1], -1, p) for z in a])
    return trim, sub, mul, power, gcd


def prime_replay(row, entry):
    p = entry['prime']; trim, sub, mul, power, gcd = finite_polynomials(p)
    cubic = trim(list(map(int, row['cubic_ascending'])))
    h = trim(list(map(int, row['integral_octic_ascending'])))
    assert len(gcd(sub(power([0, 1], p, cubic), [0, 1]), cubic)) == 1
    assert sub(power([0, 1], p**3, cubic), [0, 1]) == []
    assert len(gcd(h, trim([i*h[i] for i in range(1, len(h))]))) == 1
    xpk = [0, 1]; counts = {}; degrees = []
    for k in range(1, 9):
        xpk = power(xpk, p, h)
        total = len(gcd(sub(xpk, [0, 1]), h))-1
        old = sum(d*n for d, n in counts.items() if k % d == 0)
        assert (total-old) % k == 0
        counts[k] = (total-old)//k
        assert counts[k] >= 0
        degrees.extend([k]*counts[k])
    assert sum(degrees) == 8 and degrees == entry['octic_factor_degrees']
    # Tonelli--Shanks in F_p[z]/cubic, using a rational nonsquare since 3 is odd.
    size = p**3; odd = size-1; valuation = 0
    while odd % 2 == 0: odd //= 2; valuation += 1
    nonsquare = next(i for i in range(2, p) if pow(i, (p-1)//2, p) == p-1)
    def sqrt(a):
        assert power(a, (size-1)//2, cubic) == [1]
        z = power([nonsquare], odd, cubic)
        x = power(a, (odd+1)//2, cubic); t = power(a, odd, cubic); m = valuation
        while t != [1]:
            probe = t; i = 0
            while probe != [1]: probe = mul(probe, probe, cubic); i += 1
            assert i < m
            b = power(z, 2**(m-i-1), cubic)
            x = mul(x, b, cubic); z = mul(b, b, cubic); t = mul(t, z, cubic); m = i
        assert mul(x, x, cubic) == a
        return x
    def scalar(q):
        q = Fraction(q); return q.numerator*pow(q.denominator, -1, p) % p
    radicals = []
    for x, y in row['scaled_points']:
        a = sqrt(trim([scalar(x), -1]))
        if power(a, p*p+p+1, cubic) != [scalar(y)]: a = trim([-z for z in a])
        assert power(a, p*p+p+1, cubic) == [scalar(y)]
        radicals.append(a)
    value = power(sub(radicals[0], [-z for z in radicals[1]]), p*p+p+1, cubic)
    assert len(value) == 1 and value[0] == entry['radical_norm_mod_p']
    psi = int(pow(value[0], (p-1)//2, p) == p-1)
    assert psi == entry['psi'] == entry['independent_radical_psi']
    return {'prime': p, 'factor_degrees': degrees, 'psi': psi, 'radical_norm': value[0]}


def compute():
    # The bound is explicit: the largest frozen polynomial discriminant has 6118 digits.
    if hasattr(sys, 'set_int_max_str_digits'): sys.set_int_max_str_digits(20000)
    data = r.read(INPUT); assert data['status'] == 'PASS'
    for path, sha in data['bindings'].items(): assert r.digest((r.ROOT/path).read_bytes()) == sha
    rows = []
    for row in data['rows']:
        h = list(map(int, row['integral_octic_ascending']))
        assert h[-1] == 1
        disc = discriminant(h)
        assert disc != 0 and str(disc) == row['integral_octic_discriminant']
        primes = [prime_replay(row, e) for e in row['inert_prime_table']]
        rows.append({'id': row['id'], 'discriminant_verified': True, 'inert_prime_replays': primes})
    return {'schema': 'rank-jump.explicit-governing-octic-verification.v1', 'status': 'PASS',
            'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in (INPUT, Path(__file__), Path(r.__file__))},
            'methods': ['Bareiss determinant of the integer Sylvester matrix',
                        'distinct-degree factor counts by polynomial gcd and Frobenius powers',
                        'Tonelli--Shanks and prescribed norms in the cubic finite field'],
            'rows': rows, 'boundary': 'CAS-independent arithmetic replay; the governing interpretation uses the accompanying cochain proof.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('mode', choices=['build', 'check']); args = parser.parse_args()
    result = compute()
    if args.mode == 'build': r.write_new(OUTPUT, result)
    else: assert r.read(OUTPUT) == result
    print('PASS: three integer discriminants and 46 independent modular octic/radical replays')
