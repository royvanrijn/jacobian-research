"""Exact finite controls for AW1–AW4 (SageMath and PARI)."""

from math import lcm, prod

from sage.all import GF, Mod, QQ, ZZ, PolynomialRing, euler_phi, fundamental_discriminant, pari


ring = PolynomialRing(ZZ, "x")
x = ring.gen()


def dickson_polynomials(up_to):
    family = [ring(2), x]
    for n in range(2, up_to+1):
        family.append(x*family[n-1]-family[n-2])
    return family


def verify_field_discriminants():
    family = dickson_polynomials(15)
    for degree, trace, discriminant in [
        (3, 5, 21),
        (5, 12, 140),
        (5, 22, 120),
        (7, 23, 21),
        (9, 5, 21),
        (15, 17, 285),
    ]:
        radicand = ZZ(trace**2-4).squarefree_part()
        assert fundamental_discriminant(radicand) == discriminant
        wild = 1
        for prime, exponent in ZZ(degree).factor():
            assert ZZ(trace-2).valuation(prime) == 1
            assert ZZ(discriminant).valuation(prime) == 1
            wild *= prime**(exponent*degree)
            primary = prime**exponent
            other = degree//primary
            zero_orbit = exponent*primary+(primary-1)//2
            other_orbits = (other-1)//2*(2*exponent+1)*primary
            assert zero_orbit+other_orbits == exponent*degree+(degree-1)//2
        field_polynomial = family[degree]-trace
        assert field_polynomial.is_irreducible()
        field_discriminant = ZZ(pari(field_polynomial).nfdisc())
        assert field_discriminant == (-discriminant)**((degree-1)//2)*wild


def verify_conductor_arithmetic():
    for prime in (3, 5, 7, 11, 13, 19):
        delta = 1 if prime % 4 == 1 else 2
        for residue_degree in (1, 2) if delta == 2 else (1,):
            number_of_primes = delta//residue_degree
            conductors = []
            for j in range(1, 6):
                n = prime**j
                previous = n//prime
                cyclotomic = (j*(prime-1)-1)*previous
                base_exponent = residue_degree*(cyclotomic-delta+1)//delta
                conductor = ((prime+3)*previous+2*(delta-1))//(2*delta)
                primitive_permutation = (
                    j*n+(n-1)//2-(j-1)*previous-(previous-1)//2
                )
                assert primitive_permutation == number_of_primes*(
                    base_exponent+residue_degree*conductor
                )
                next_cyclotomic = ((j+1)*(prime-1)-1)*n
                next_base_exponent = residue_degree*(next_cyclotomic-delta+1)//delta
                tower_conductor = (n+delta-1)//delta
                assert next_base_exponent-prime*base_exponent == (
                    residue_degree*(prime-1)*tower_conductor
                )
                assert conductor <= tower_conductor
                assert (conductor == tower_conductor) == (prime == 3)
                conductors.append((n-previous)*conductor)
                assert sum(conductors) == (
                    (prime+3)*(n*n-1)//(2*delta*(prime+1))+(delta-1)*(n-1)//delta
                )


def verify_actual_moving_prime():
    q = ZZ(11690712001)
    B, m, c, auxiliary = ZZ(15), ZZ(45), ZZ(511), ZZ(479)
    assert q.is_prime(proof=True) and auxiliary.is_prime(proof=True)
    A = 2*(B*c-1)
    assert B*c-1 == 2**4*auxiliary
    assert q % m == 1 and pow(q, m, 8*A) == (c*c+2*A) % (8*A)
    n, remainder = divmod(q**m-c*c, A)
    assert remainder == 0 and n > 0 and n % 8 == 2
    d = n*(B*B*n+2)
    assert (B*n+c)**2-d == q**m
    assert d % q != 0
    for prime in (ZZ(3), ZZ(5)):
        assert d.valuation(prime) == 1 and (d//prime) % prime == 4 % prime
    fixed_class = QQ(c*(B*c-2))/(B*c-1)
    assert fixed_class.valuation(auxiliary) == -1
    residue = GF(q)
    assert 2*residue(B*n+c) == residue(fixed_class)
    phase = residue(fixed_class)**((q-1)//m)
    assert phase == 10994014720 and phase.multiplicative_order() == m
    assert (phase**(-2)).multiplicative_order() == m


def verify_composite_mirror_arithmetic():
    for p, exponent, m0, expected_residue, expected_primes in [
        (3, 2, 5, 2, 2),
        (5, 1, 9, 6, 1),
        (7, 2, 3, 2, 1),
        (7, 1, 5, 2, 2),
    ]:
        delta = 1 if p % 4 == 1 else 2
        d = int(Mod(p, m0).multiplicative_order())
        residue = lcm(delta, d)
        if delta == 2 and d % 4 == 0 and pow(p, d//2, m0) == m0-1:
            residue //= 2
        assert residue == expected_residue
        assert euler_phi(m0)//residue == expected_primes
        contributions = []
        for j in range(1, exponent+1):
            n = p**j
            base = ((p+3)*p**(j-1)+2*(delta-1))//(2*delta)
            composite = (p+3)*p**(j-1)//2
            assert delta*base-(delta-1) == composite
            contributions.append((n-n//p)*composite)
        assert sum(contributions) == (p+3)*(p**(2*exponent)-1)//(2*(p+1))


if __name__ == "__main__":
    verify_field_discriminants()
    verify_conductor_arithmetic()
    verify_actual_moving_prime()
    verify_composite_mirror_arithmetic()
    print("AW1–AW4 Dickson discriminants, ray conductor arithmetic, and genuine prime verified")
