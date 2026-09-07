#!/usr/bin/env sage
"""Reproduce the retracted standard-P1 simple-pole P3 CRT experiment.

Those seeds are now known to give P2 the wrong residual-I2 component and to
be first-order Hensel-obstructed.  The ambient and P1/P2 are rational.  The historical missing generator was expected over
Q(sqrt(-6)); this script first ranks all choices of split-prime embeddings for
its pole.  Exact section reconstruction is added below the ranking gate.
"""

from itertools import product as itertools_product
from sage.all import *


K.<s> = QuadraticField(-6)
KT.<T> = PolynomialRing(K)
primes = (11, 29, 53)
modulus = prod(primes)


def short_values(residues, roots, count=12, coefficient_bound=12):
    root_crt = ZZ(CRT_list(list(roots), list(primes)))
    residue_crt = ZZ(CRT_list(list(residues), list(primes)))
    lattice = matrix(ZZ, [
        [modulus, 0, 0],
        [-root_crt, 1, 0],
        [residue_crt, 0, 1],
    ]).LLL()
    values = {}
    for coefficients in itertools_product(
        range(-coefficient_bound, coefficient_bound+1), repeat=3
    ):
        if coefficients == (0, 0, 0):
            continue
        a, b, c = map(ZZ, vector(ZZ, coefficients)*lattice)
        if c == 0 or any(c % prime == 0 for prime in primes):
            continue
        common = gcd((a, b, c))
        a, b, c = a//common, b//common, c//common
        if c < 0:
            a, b, c = -a, -b, -c
        value = K(a+b*s)/c
        score = a*a+6*b*b+c*c
        values[value] = min(score, values.get(value, infinity))
    return tuple(sorted(values.items(), key=lambda item: item[1])[:count])


def reduce_at(value, prime, root):
    rational, radical = K(value).list()
    field = GF(prime)
    return (
        field(rational.numerator())/field(rational.denominator())
        + field(root)*field(radical.numerator())/field(radical.denominator())
    )


def polynomial_square_root(value, prescribed_leading=None):
    if not value:
        return value.parent().zero()
    if value.degree() % 2:
        return None
    half = value.degree()//2
    if prescribed_leading is None:
        leading_roots = value[value.degree()].sqrt(all=True)
    else:
        leading_roots = (prescribed_leading,)
    for leading_root in leading_roots:
        if leading_root**2 != value[value.degree()]:
            continue
        coefficients = [value.base_ring().zero()]*(half+1)
        coefficients[half] = leading_root
        for index in range(half-1, -1, -1):
            target = half+index
            known = sum(
                coefficients[left]*coefficients[target-left]
                for left in range(index+1, half)
                if index < target-left < half
            )
            coefficients[index] = (value[target]-known)/(2*leading_root)
        root = value.parent()(coefficients)
        if root**2 == value:
            return root
    return None


ranked = []
for roots in itertools_product((4, 7), (9, 20), (10, 43)):
    for value, score in short_values((5, 26, 18), roots):
        ranked.append((score, roots, value))

for rank, (score, roots, value) in enumerate(sorted(ranked)[:32], 1):
    print(
        f"Q80P3CRT|stage=pole_rank|rank={rank}|score={score}|"
        f"roots={roots}|z={value}",
        flush=True,
    )

chosen_roots = (4, 9, 10)
x_residues = (
    (4, 23, 15),
    (7, 6, 14),
    (5, 14, 11),
    (10, 14, 44),
    (8, 5, 50),
    (0, 19, 2),
    (3, 25, 42),
)
for index, residues in enumerate(x_residues):
    values = short_values(residues, chosen_roots, count=8)
    print(
        f"Q80P3CRT|stage=x_rank|index={index}|values="
        + ";".join(f"{value}@{score}" for value, score in values),
        flush=True,
    )

z = -5-3*s
p = QQ(-105)/8
q = 18-2*p
e = (p-42)**2/36
rho = QQ(1)/49
node = QQ(51)/2401
A = T**2*(-3+p*T+q*T**2+(p-42)*T**3)
B = (
    2*T**3+(2*p+e-45)*T**4+(-9*p-4*e+186)*T**5
    +(12*p+6*e-299)*T**6+(-5*p-4*e+210)*T**7+e*T**8
)
Z = T-z
top_matrix = matrix(K, [[1, 1], [rho**5, rho**6]])
x_choices = tuple(
    tuple(value for value, _ in short_values(residues, chosen_roots, count=8))
    for residues in x_residues[:5]
)

filter_prime = 97
filter_roots = tuple(GF(filter_prime)(-6).sqrt(all=True))
filter_data = []
for filter_root in filter_roots:
    filter_ring = PolynomialRing(GF(filter_prime), "t")
    filter_T = filter_ring.gen()
    filter_A = filter_ring([
        reduce_at(coefficient, filter_prime, filter_root) for coefficient in A
    ])
    filter_B = filter_ring([
        reduce_at(coefficient, filter_prime, filter_root) for coefficient in B
    ])
    filter_z = reduce_at(z, filter_prime, filter_root)
    filter_data.append((filter_root, filter_ring, filter_A, filter_B, filter_T-filter_z))

y_residues = (
    (3, 4, 47), (1, 23, 34), (8, 3, 48), (10, 21, 43), (5, 21, 1),
    (3, 20, 7), (5, 16, 10), (5, 16, 34), (0, 12, 31), (4, 9, 10),
)
tested = 0
filter_passed = 0
hits = []
for low_coefficients in itertools_product(*x_choices):
    tested += 1
    known = sum(value*T**index for index, value in enumerate(low_coefficients))
    right = vector(K, [
        3*(1-z)**2-known(1),
        node*(rho-z)**2-known(rho),
    ])
    x5, x6 = top_matrix.solve_right(right)
    X = known+x5*T**5+x6*T**6
    survives = False
    for filter_root, filter_ring, filter_A, filter_B, filter_Z in filter_data:
        try:
            filter_X = filter_ring([
                reduce_at(coefficient, filter_prime, filter_root) for coefficient in X
            ])
        except ZeroDivisionError:
            continue
        filter_right = filter_X**3+filter_A*filter_X*filter_Z**4+filter_B*filter_Z**6
        if polynomial_square_root(filter_right) is not None:
            survives = True
            break
    if not survives:
        continue
    filter_passed += 1
    exact_right = X**3+A*X*Z**4+B*Z**6
    Y = polynomial_square_root(exact_right)
    if Y is None:
        continue
    signs = []
    marking_ok = True
    for prime, root, residues in zip(primes, chosen_roots, zip(*y_residues)):
        reduced = tuple(ZZ(reduce_at(Y[index], prime, root)) for index in range(10))
        target = tuple(residues)
        if reduced == target:
            signs.append(1)
        elif reduced == tuple(-value % prime for value in target):
            signs.append(-1)
        else:
            marking_ok = False
            break
    if marking_ok:
        hits.append((X, Y, tuple(signs)))

print(
    f"Q80P3CRT|stage=exact_search|tested={tested}|filter97={filter_passed}|"
    f"hits={len(hits)}",
    flush=True,
)
for X, Y, signs in hits:
    print(f"Q80P3CRT|HIT|z={z}|X={X}|Y={Y}|signs={signs}", flush=True)
print(
    "Q80P3CRT|SUMMARY|status="
    + ("EXACT_P3_OVER_QSQRT_MINUS6" if hits else "BOUNDED_CRT_BOX_NO_HIT"),
    flush=True,
)
