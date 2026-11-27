"""Exact examples for PA1–PA5 at exponents 5, 7, and 11."""
from sage.all import GF, QQ, ZZ, Matrix, PolynomialRing, pari

ring = PolynomialRing(QQ, "T")
T = ring.gen()


def fixture(exponent, p, h, residue_rows, expected_local_type):
    p, h = ZZ(p), ZZ(h)
    primes = (p-h, p, p+h)
    assert all(q.is_prime() for q in primes)
    assert h.valuation(exponent) >= 2
    A, B, C = primes[0]**exponent, 4*p**exponent, primes[2]**exponent
    E = (4*A*C - (A+C-B)**2) // 4
    assert E % (h*h) == 0
    form = E // (h*h)
    assert form > 0 and form.valuation(exponent) == 1
    assert (-form // exponent) % exponent == 2
    u1, u2 = (A+B-C)//2, (B+C-A)//2
    assert (u1-2*p**exponent).valuation(exponent) >= 3
    assert (u2-2*p**exponent).valuation(exponent) >= 3
    assert (u1*u1 + E) // 4 == (primes[0]*p)**exponent
    assert (u2*u2 + E) // 4 == (primes[2]*p)**exponent
    assert all(E % q != 0 for q in primes)

    rows = []
    for q, radical, root, exponents in residue_rows:
        field = GF(q)
        r, zeta = field(radical), field(root)
        assert r*r == -field(E)
        assert zeta != 1 and zeta**exponent == 1
        values = [(field(tr) + r)/2 for tr in (u1, u2)]
        assert tuple(x**((q-1)//exponent) for x in values) == tuple(zeta**j for j in exponents)
        rows.append(exponents)
    assert Matrix(GF(exponent), rows).rank() == 2

    n = primes[0]*p
    previous, current = ring(2), T
    for _ in range(2, exponent+1):
        previous, current = current, T*current - n*previous
    polynomial = current - u1
    assert polynomial.is_irreducible()
    discriminant = ZZ(pari(-form.squarefree_part()).quaddisc())
    e = exponent-1 if exponent % 4 == 1 else (exponent-1)//2
    expected_disc = exponent**(exponent-1-(exponent-1)//e) * (abs(discriminant)//exponent)**((exponent-1)//2)
    field = pari(polynomial).nfinit()
    assert ZZ(pari(polynomial).nfdisc()) == expected_disc
    assert sorted((int(q[2]), int(q[3])) for q in field.idealprimedec(exponent)) == sorted(expected_local_type)
    return primes, E, (u1, u2), discriminant, polynomial, field


_, _, _, D, quintic, _ = fixture(
    5, 419, 150, [(11, 4, 4, (1, 3)), (71, 6, 54, (2, 4))], [(1, 1), (4, 1)])
assert D == -16290203579467969335
assert quintic == T**5 - 563555*T**3 + 63518847605*T + 3288814754552
fixture(7, 1153, 294, [(113, 13, 49, (2, 1)), (127, 23, 4, (4, 4))], [(1, 1), (3, 2)])
fixture(11, 4297, 726, [(23, 10, 2, (4, 4)), (67, 29, 64, (0, 5))], [(1, 1), (5, 1), (5, 1)])

marks, E, (u1, u2), _, polynomial, field = fixture(
    5, 421, 150, [(31, 10, 16, (0, 3)), (61, 9, 9, (2, 3))], [(1, 1), (4, 1)])
assert polynomial == T**5 - 570455*T**3 + 65083781405*T + 3167699018548
first, middle, last = marks
expected = [(QQ(2)/(4-2**5), QQ(1)/2), (1, 1), (QQ(1)/2, QQ(2)/(4-2**5))]
normalized = [(QQ(middle**5)/u1, (u2-u1)//2), (QQ(first**5)/u1, (u2-u1)//2),
              ((u1-u2)//2, QQ(middle**5)/u2)]
for q, row, target in zip(marks, normalized, expected):
    local = GF(q)
    for value, comparison in zip(row, target):
        assert local(value)**((q-1)//5) == local(comparison)**((q-1)//5)
for q, splitting_type in zip(marks, [[(1, 5)], [(1, 1)]*5, [(1, 1)]*5]):
    assert sorted((int(p[2]), int(p[3])) for p in field.idealprimedec(q)) == splitting_type

print("Verified the affine reflection fixtures, exact discriminants and local prime decompositions.")
