#!/usr/bin/env python3
"""Exact finite certificates for specific attempted-construction obstructions.

No arbitrary-N recurrence or nonexistence assertion is tested.
"""
if not __debug__:
    raise RuntimeError('Exact certificates require assertions; disable optimization.')

from math import gcd
import json

def trim(f):
    f = list(f)
    while len(f) > 1 and f[-1] == 0:
        f.pop()
    return f

def add(a, b):
    out = [0] * max(len(a), len(b))
    for i, x in enumerate(a): out[i] += x
    for i, x in enumerate(b): out[i] += x
    return trim(out)

def scale(a, c): return trim([c*x for x in a])

def mul(a, b):
    out = [0] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i+j] += x*y
    return trim(out)

def rem(a, b, p):
    a = trim([x % p for x in a])
    b = trim([x % p for x in b])
    while a != [0] and len(a) >= len(b):
        j = len(a)-len(b)
        q = a[-1]*pow(b[-1], -1, p) % p
        for i, x in enumerate(b): a[i+j] = (a[i+j]-q*x) % p
        a = trim(a)
    return a

K = [3,-2,1]
Z = [-3,6,-1]
d = [-3,0,1]
D = mul(d,d)
N = add(scale(mul(K,K),2),scale(D,-1))
assert N == [9,-24,26,-8,1]
assert add(mul(Z,Z),scale(mul(K,K),-3)) == scale(D,-2)
P = [9,0,2,-16,9]
Q = [81,-48,2,0,1]
assert add(mul(N,N),scale(mul(D,D),8)) == mul(P,Q)
assert gcd(*P) == gcd(*Q) == 1
assert P[-1] % 5 != 0 and Q[-1] % 5 != 0
assert (2*2+1) % 5 == 0  # 5 splits in Z[i].

factor_tests = 0
reductions = []
for f in (P,Q):
    h = [x*pow(f[-1], -1, 5) % 5 for x in f]
    reductions.append(h)
    for a in range(5):
        assert rem(h, [a,1], 5) != [0]
        factor_tests += 1
    for a in range(5):
        for b in range(5):
            assert rem(h, [a,b,1], 5) != [0]
            factor_tests += 1
assert reductions == [[1,0,3,1,1],[1,2,2,0,1]]

# Every stage can repair all previous bad primes yet leave two new ones.
Dold = 1
old_primes = []
repair_tests = 0
for p,q in [(3,7),(11,19),(23,31),(43,47)]:
    assert p % 4 == q % 4 == 3
    assert all(gcd(p*q,r) == 1 for r in old_primes)
    value = Dold*Dold*p*q
    assert (value-1) % 4 == 0
    u = (value-1)//4
    assert u > 0 and 1+4*u == value
    for r in old_primes:
        assert value % (r*r) == 0 and value % (r*r*r) != 0
        repair_tests += 1
    assert value % p == value % q == 0
    assert value % (p*p) != 0 and value % (q*q) != 0
    assert value//(Dold*Dold) == p*q
    repair_tests += 3
    old_primes.extend((p,q))
    Dold *= p*q

print(json.dumps({
    'status': 'PASS',
    'polynomial_identity_checks': 2,
    'exhaustive_F5_linear_quadratic_divisor_tests': factor_tests,
    'monic_F5_quartic_coefficients_low_to_high': reductions,
    'finite_square_repair_invariant_checks': repair_tests,
    'arithmetic': 'exact integers and F5 arithmetic',
    'arbitrary_N_recurrence_tested': False,
    'finite_checks_prove_all_E_impossibility': False,
}, indent=2))
