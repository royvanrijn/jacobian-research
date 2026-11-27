"""Exact fixtures for the shifted local action, mass and prime-unit realization."""

from collections import Counter
from itertools import product

from sage.all import CyclotomicField, QQ, ZZ, pari, power_mod

from verify_two_prime_mirror_intermediate_fields import automorphism_count, element_order_exponent, group_elements


def check_shifted_orbits():
    for exponents, shifts in (((2, 1), (0, 1)), ((3, 2), (0, 1)), ((2, 1), (1, 0)),
                              ((2, 1), (0, 0))):
        elements = set(group_elements(3, exponents))
        actual = Counter()
        while elements:
            start = next(iter(elements))
            point = start
            length = 0
            while point in elements:
                elements.remove(point)
                length += 1
                point = tuple((4*x+t) % (3**a) for x, t, a in zip(point, shifts, exponents))
            assert point == start
            actual[length] += 1
        phase = max((a for a, t in zip(exponents, shifts) if t), default=0)
        previous = 0
        expected = Counter()
        for s in range(max(exponents)+1):
            cumulative = 3**sum(min(a, s+1) for a in exponents) if s >= phase else 0
            orbits, remainder = divmod(cumulative-previous, 3**s)
            assert remainder == 0
            expected[3**s] = orbits
            previous = cumulative
        assert actual == expected

    for s in range(1, 4):
        zeta = CyclotomicField(3**(s+1)).gen()
        polynomial = (zeta+zeta**(-1)).minpoly()
        assert ZZ(pari(polynomial).nfdisc()).valuation(3) == ((2*s+1)*3**s-1)//2


def check_phase_masses():
    exponents = (2, 1, 1)
    shifts = (0, 1, 0)
    elements = group_elements(3, exponents)
    orders = {x: element_order_exponent(x, 3, exponents) for x in elements}

    def phase(x, level):
        return sum(t*(v//3**(a-level))
                   for v, a, t in zip(x, exponents, shifts) if a >= level) % 3

    sizes = [sum(order <= j for order in orders.values()) for j in range(3)]
    for a, b in ((1, 1), (2, 1)):
        for cutoff in (0, 1, 2):
            def allowed(x, j):
                return orders[x] <= j and (j <= cutoff or phase(x, j) == 0)

            injections = 0
            for x in elements:
                if orders[x] != a or not allowed(x, a):
                    continue
                line = {tuple(scalar*3**(a-1)*v % 3**e for v, e in zip(x, exponents))
                        for scalar in range(3)}
                for y in elements:
                    if orders[y] != b or not allowed(y, b):
                        continue
                    leading = tuple(3**(b-1)*v % 3**e for v, e in zip(y, exponents))
                    if leading not in line:
                        injections += 1
            h_a = sum(allowed(x, a) for x in elements)
            h_b = sum(allowed(x, b) for x in elements)
            assert injections == (h_a-sizes[a-1])*(h_b-3*sizes[b-1])
            assert injections % automorphism_count(3, a, b) == 0


def check_prime_and_unit():
    m, a, M = 9, 10, 9
    A = -a**m
    p0 = ZZ(40423)
    p, q = 3*p0, ZZ(4000121269)
    B = -(p+q)//2
    assert p0.is_prime(proof=True) and q.is_squarefree()
    assert p - q == 4*A and p.gcd(q) == 1
    assert p0.gcd(a) == q.gcd(a) == 1
    assert (p0-4*a**m) % 3**(M-1) == 0

    for dyadic, P in ((4, ZZ(17666533719169548857641)),
                      (20, ZZ(2644283537742303038521))):
        assert P.is_prime(proof=True) and P % m == 1
        gamma_frame = QQ(4)/p
        unit_frame = QQ(q)/p
        for frame, expected in ((gamma_frame, m), (unit_frame, 1)):
            residue = (frame.numerator()*power_mod(frame.denominator(), -1, P)) % P
            residue = power_mod(residue, (P-1)//m, P)
            assert power_mod(residue, expected, P) == 1
            if expected > 1:
                assert power_mod(residue, expected//3, P) != 1
        t = A*P**m
        F = t*t + B*P**m + 1
        assert F > 0 and F % 32 == dyadic
        assert (t+1).valuation(3) >= M
        assert (F+3*(t-1)**2).valuation(3) >= M
        assert F.valuation(3) == 1 and F//3 % 3 == 2
        assert all(F.valuation(prime) == 1 for prime in (p0, 347, 719, 16033))
        J = 2*A*t+B
        assert J % (p*q) == 0
        assert J*J-4*A*A*F == p*q


if __name__ == "__main__":
    check_shifted_orbits()
    check_phase_masses()
    check_prime_and_unit()
    print("RAM1–RAM3 affine orbit, filtered mass and genuine prime-unit fixtures verified.")
