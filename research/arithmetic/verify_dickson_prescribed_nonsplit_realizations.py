"""Exact finite fixtures for IG3a and IG3; the all-level proofs are written.

Run from the repository root with
    sage -python research/arithmetic/verify_dickson_prescribed_nonsplit_realizations.py

The field computations certify three root-field discriminants and one degree-81
normal-closure discriminant. They do not certify the asymptotic field count or
the assertions for arbitrary levels.
"""

from collections import Counter, deque
import json

from sage.all import NumberField, PolynomialRing, QQ, ZZ, pari, valuation
from sage.env import SAGE_VERSION


ring = PolynomialRing(QQ, "x")
x = ring.gen()


def dickson(degree):
    previous, current = ring(2), x
    for _ in range(1, degree):
        previous, current = current, x*current - previous
    return current


def coefficient(depth):
    return ZZ(2*3**(depth + 2) - 3)


def discriminant_data(level, depth):
    m = ZZ(3**level)
    r = level - depth
    assert level >= 3 and 1 <= depth <= level - 2
    local_conductors = {1: ZZ(0)}
    local_conductors.update({k: ZZ((3**(k - 1) + 1)//2) for k in range(2, r)})
    local_conductors[r] = ZZ((m + 1)//2)
    local_sum = sum(2*3**(k - 1)*a for k, a in local_conductors.items())
    C = ZZ(3**depth*local_sum)
    base = ZZ((2*level - 1)*(m//3) - 1)
    exponent = ZZ(m*base + 2*C)
    closed = (2*level + 1)*m**2/3 + (3**(2*level - depth - 2) - 21*3**depth)/QQ(4)
    assert closed in ZZ and exponent == closed
    return m, r, local_conductors, C, exponent


def verify_trace_depths():
    field = NumberField(x**2 + 3, "pi")
    pi = field.gen()
    zeta = (-1 + pi)/2
    eta = (1 + pi**3)/(1 - pi**3)
    assert eta.norm() == 1 and zeta**3 == 1 and zeta != 1
    # Independently evaluate both torsion orientations and several free units.
    for depth in (0, 1, 2, 3):
        for torsion in (1, 2):
            for unit in (1, 2, 4):
                gamma = zeta**torsion * eta**(3**depth*unit)
                T = gamma.trace()
                assert valuation(T - 2, 3) == 1
                assert valuation(T + 1, 3) == depth + 2
    assert x**3 - 3*x - 2 == (x - 2)*(x + 1)**2

    for level, depth in ((3, 1), (4, 1), (4, 2), (5, 3)):
        m = 3**level
        A = coefficient(depth)
        assert A % 2 == 1 and valuation(A, 3) == 1
        for parameter in (QQ(1), -QQ(1)/5):
            assert parameter % 3 == 1
            assert all(valuation(parameter, p) == 0 for p in (2, *A.prime_divisors()))
            T = 2 + A*parameter**m
            assert valuation(T - 2, 3) == 1
            assert valuation(T + 1, 3) == depth + 2
            assert (T**2 - 4) % 8 == 5
            assert ((T**2 - 4)/(-3)) % 3 == 1
            shifted = dickson(m)(x - 1) - T
            assert shifted.leading_coefficient() == 1
            assert valuation(shifted[0], 3) == 1
            assert all(valuation(c, 3) >= 1 for c in shifted.list()[:-1])


def compose(left, right, modulus):
    a, b = left
    c, d = right
    return ((a*c) % modulus, (a*d + b) % modulus)


def affine_power(element, exponent, modulus):
    result = (1, 0)
    while exponent:
        if exponent % 2:
            result = compose(result, element, modulus)
        element = compose(element, element, modulus)
        exponent //= 2
    return result


def verify_groups_and_character_restriction():
    summaries = []
    for level, depth in ((3, 1), (4, 1), (4, 2)):
        m, r, conductors, C, exponent = discriminant_data(level, depth)
        sigma, tau = (4, 1), (1, 3**depth)
        group, queue = {(1, 0)}, deque([(1, 0)])
        while queue:
            item = queue.popleft()
            for generator in (sigma, tau):
                product = compose(item, generator, m)
                if product not in group:
                    group.add(product)
                    queue.append(product)
        assert len(group) == 3**(level + r - 1)
        assert all(a % 3 == 1 and (b - (a - 1)//3) % 3**depth == 0 for a, b in group)
        translations = {item for item in group if item[0] == 1}
        assert translations == {(1, b) for b in range(0, m, 3**depth)}
        assert compose(sigma, tau, m) != compose(tau, sigma, m)
        for a, b in group:
            if a == 4:
                assert affine_power((a, b), m//3, m) != (1, 0)
                assert affine_power((a, b), m, m) == (1, 0)

        # Restrict all m global cyclic characters to the local subgroup.
        orders = Counter(ZZ(3**r)//ZZ(k).gcd(3**r) for k in range(m))
        assert orders[1] == 3**depth
        counted_C = sum(orders[3**k]*a for k, a in conductors.items())
        assert counted_C == C
        summaries.append({"n": level, "s": depth, "local_group_order": len(group),
                          "relative_exponent_per_base_prime": int(C),
                          "global_three_exponent": int(exponent)})
    return summaries


def verify_root_fields():
    summaries = []
    for level, depth, parameter in ((3, 1, QQ(1)), (3, 1, -QQ(1)/5), (4, 2, QQ(1))):
        m = 3**level
        A = coefficient(depth)
        T = 2 + A*parameter**m
        H = A*parameter*(4 + A*parameter**m)
        D = ZZ(H.numerator()*H.denominator()).squarefree_part()
        assert D % 4 == 1 and valuation(D, 3) == 1
        assert D % 8 == 5 and (-D//3) % 3 == 1
        denominator = parameter.denominator()
        polynomial = denominator**m*(dickson(m)(x/denominator) - T)
        assert polynomial.leading_coefficient() == 1
        assert all(c in ZZ for c in polynomial)
        delta = ((2*level + 1)*m - 1)//2
        actual = ZZ(pari(polynomial).nfdisc())
        expected = 3**delta*abs(D//3)**((m - 1)//2)
        assert abs(actual) == expected
        summaries.append({"n": level, "s": depth, "t": str(parameter), "A": int(A),
                          "D": str(D), "root_degree": int(m), "three_exponent": int(delta)})
    return summaries


def verify_local_normal_closure_fixture():
    # The raw trace Tr(zeta_3 eta^3)=683/343 gives a global degree-81
    # realization of the local normal closure at three. It is not an IG3
    # global-family parameter, since its quadratic field is Q(zeta_3).
    polynomial = 7**27*dickson(27)(x/7) - 683*7**24
    real_cyclotomic = pari.polsubcyclo(27, 9)
    if real_cyclotomic.type() == "t_VEC":
        real_cyclotomic = real_cyclotomic[0]
    composita = pari(polynomial).polcompositum(real_cyclotomic)
    assert [int(f.poldegree()) for f in composita] == [81, 81, 81]
    actual = abs(ZZ(composita[0].nfdisc()))
    assert actual == 3**282*7**72
    assert ZZ(pari(polynomial).nfdisc()).abs() == 3**94*7**24
    # At the two rationally split primes of the mirror base, the global
    # degree-486 normal closure has 486/81 completions of this degree.
    assert discriminant_data(3, 1)[4] == (486//81)*valuation(actual, 3)


if __name__ == "__main__":
    verify_trace_depths()
    groups = verify_groups_and_character_restriction()
    fields = verify_root_fields()
    verify_local_normal_closure_fixture()
    print(json.dumps({"sage": SAGE_VERSION, "pari": str(pari.version()),
                      "groups_and_conductors": groups, "root_fields": fields,
                      "local_normal_closure_degree": 81, "local_fixture_three_exponent": 282},
                     indent=2))
    print("IG3a–IG3 finite trace, affine, conductor and field-discriminant fixtures verified.")
