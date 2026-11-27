"""Exact norm, discriminant, prime-type and compositum checks for PC1–PC4."""
from sage.all import GF, QQ, ZZ, Matrix, PolynomialRing, pari

ring = PolynomialRing(QQ, "T")
T = ring.gen()
rational = ring.fraction_field()
t = rational(T)
F = -3*t**6 + 24*t**4 + 24*t*t + 4
quadratic_ring = PolynomialRing(rational, "R")
R = quadratic_ring.gen()
mirror = rational.extension(R**2+3*F, "rho")
rho = mirror.gen()
Z = 24+4/(t*t)
L = 4*Z-28
xi = -L/3+4*rho/(9*t**3)
conjugate = -L/3-4*rho/(9*t**3)
assert xi*conjugate == (Z/3)**3
assert xi+conjugate == -2*L/3
assert 27*(4*(Z/3)**3-(2*L/3)**2) == 64*F/t**6

binary = PolynomialRing(QQ, names=("P", "H"))
P, H = binary.gens()
sextic = -3*P**6+24*P**4*H*H+24*P*P*H**4+4*H**6
coefficient = -36*(H*H+6*P*P)
constant = 288*P*H*H+1224*P**3
assert -4*coefficient**3-27*constant*constant == 3**6*64*sextic
assert (T**3-7*T+7).is_irreducible() and (T**3-7*T+7).discriminant() == 49

field7 = GF(7)
for p in field7:
    for h in field7:
        if not p*(p-h)*(p+h):
            continue
        value = -3*p**6+24*p**4*h**2+24*p*p*h**4+4*h**6
        assert value and value.is_square() == (h == 0)
        assert ((p-h)**2 == (p+h)**2) == (h == 0)


def fixture(p, h, expected_D, rows, expected_second_types):
    p, h = ZZ(p), ZZ(h)
    primes = p-h, p, p+h
    assert h % 18 == 0 and all(q.is_prime() and q > 7 for q in primes)
    f = -3*p**6+24*p**4*h*h+24*p*p*h**4+4*h**6
    D = f.squarefree_part()
    assert D == expected_D and D.valuation(3) == 1 and D % 8 == 5
    assert (D//3) % 3 == 2
    assert all(f % q for q in primes)
    u1, u2 = (p-h)**3+p**3-(p+h)**3, (p+h)**3+p**3-(p-h)**3
    chars = []
    for q, branch, zeta, expected in rows:
        local = GF(q)
        branch, zeta = local(branch), local(zeta)
        assert branch*branch == f and zeta != 1 and zeta**3 == 1
        elements = [(local(u)+branch)/2 for u in (u1, u2)]
        assert tuple(a**((q-1)//3) for a in elements) == tuple(zeta**i for i in expected)
        chars.append(expected)
    assert Matrix(GF(3), chars).rank() == 2

    first = pari((p+h)**3*(T**3+1)-((p-h)*T+p)**3).polredabs()
    second = pari(T**3-36*(h*h+6*p*p)*T+288*p*h*h+1224*p**3)
    assert first.nfdisc() == -D//3 and second.nfdisc() == 81*D
    index_square = second.poldisc() // second.nfdisc()
    assert index_square == 24**2 * (f//D)
    first_field, second_field = first.nfinit(), second.nfinit()
    assert sorted((int(q[2]),int(q[3])) for q in first_field.idealprimedec(3)) == [(1,1)]*3
    assert [(int(q[2]),int(q[3])) for q in second_field.idealprimedec(3)] == [(3,1)]
    for q, expected in zip(primes,expected_second_types):
        split = q == p or q % 7 in (1,6)
        assert expected == ([(1,1)]*3 if split else [(1,3)])
        assert sorted((int(v[2]),int(v[3])) for v in second_field.idealprimedec(q)) == expected

    compositum = first.polcompositum(second)[0]
    assert compositum.poldegree() == 9
    assert compositum.nfdisc() == -3**11*D**4


fixture(383,126,-367119482556963,[(109,46,63,(0,2)),(199,99,106,(2,0))],
        [[(1,3)],[(1,1)]*3,[(1,3)]])
fixture(137,126,675712099077,[(37,17,26,(2,0)),(73,22,8,(2,1))],
        [[(1,3)],[(1,1)]*3,[(1,3)]])
fixture(61,18,-37385037867,[(19,1,7,(2,0)),(307,118,289,(2,1))],
        [[(1,1)]*3,[(1,1)]*3,[(1,3)]])

print("Verified the mirror norm and the two companion fields, their prime types and degree-nine discriminants.")
