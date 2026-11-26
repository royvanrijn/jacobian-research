"""Replay the finite-field certificates for the fourteen-shift CM test."""
from math import prod

from sage.all import GF, PolynomialRing, TermOrder

from verify_cm_newton_full_width_reconstruction_certificates import reduced_residual_rows


WEIGHTS = (2, 1, 2, 1, 2)
NAMES = ("r", "ell", "h", "v", "w")
PRIMES = (331, 337, 347, 349, 353, 359, 367, 401, 65521)


def weight(exponent):
    return sum(power * scale for power, scale in zip(exponent, WEIGHTS))


def replay(primes=PRIMES):
    rows = [row.dict() for row in reduced_residual_rows(5, 13, even_linear=True)]
    assert len(rows) == 7
    leading_forms = []
    for index, row in enumerate(rows[:4], start=6):
        assert max(map(weight, row)) == index
        # Terms involving v or w have smaller weight in r, ell, h.
        leading_forms.append({
            tuple(term[:3]): coefficient for term, coefficient in row.items()
            if not any(term[3:]) and weight(term) == index
        })

    for prime in primes:
        field = GF(prime)
        leading_ring = PolynomialRing(
            field, names=NAMES[:3], order=TermOrder("wdegrevlex", WEIGHTS[:3]),
        )
        # The first three leading forms alone have a curve at 359.
        prefix_length, leading_length = (4, 68) if prime == 359 else (3, 84)
        leading_ideal = leading_ring.ideal([
            leading_ring({term: field(coefficient) for term, coefficient in form.items()})
            for form in leading_forms[:prefix_length]
        ])
        # A zero-dimensional homogeneous ideal supported at the origin
        # has only the origin as a geometric zero.
        assert leading_ideal.dimension() == 0
        assert leading_ideal.vector_space_dimension() == leading_length
        assert all(polynomial.constant_coefficient() == 0 for polynomial in leading_ideal.gens())

        residual_ring = PolynomialRing(
            field, names=NAMES, order=TermOrder("wdegrevlex", WEIGHTS),
        )
        residual_ideal = residual_ring.ideal([
            residual_ring({tuple(term): field(coefficient) for term, coefficient in row.items()})
            for row in rows
        ])
        assert tuple(residual_ideal.groebner_basis()) == (residual_ring.one(),)
        print(f"{prime}: leading quotient has length {leading_length}, residual ideal contains 1", flush=True)

    print(f"Product of primes: {prod(primes)}")


if __name__ == "__main__":
    replay()
