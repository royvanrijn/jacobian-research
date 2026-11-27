"""Exact crossing and fractional successor tests for differential order twelve."""

from fractions import Fraction
from math import gcd, lcm


def check_first_coordinate_six():
    u = 6
    rows = set()
    for v in range(u + 1, u * (u - 1) + 1):
        factor = gcd(u, v)
        for numerator in range(1, factor):
            f1, f2 = numerator * u // factor, numerator * v // factor
            if f1 < 2:
                continue
            divisor = gcd(f2 - 1, f1 - 1)
            rho, sigma = (f2 - 1) // divisor, (1 - f1) // divisor
            for r in range(1, u):
                for t in range(r):
                    if rho * (u - r) + sigma * (v - t):
                        continue
                    delta = Fraction(v - t, rho) - t
                    x, y = r - delta * Fraction(sigma, rho), t + delta
                    if divisor == 1 and y >= 2 and x == y - Fraction(1, rho):
                        continue
                    rows.add((v, (r, t), (f1, f2), (rho, sigma)))
    assert rows == {
        (10, (1, 0), (3, 5), (2, -1)),
        (12, (2, 0), (2, 4), (3, -1)),
        (15, (1, 0), (4, 10), (3, -1)),
        (30, (1, 0), (5, 25), (6, -1)),
    }


def partitions(total, largest):
    if total == 0:
        yield ()
    for first in range(min(total, largest), 0, -1):
        for rest in partitions(total - first, first):
            yield (first, *rest)


def scalar_edges(state):
    ell, mu, beta, previous = state
    largest_total = (2 * mu * beta + 1) // mu
    for total in range(5, largest_total + 1, 2):
        scaled_slope = ell - Fraction(beta * total, total * mu - 1)
        index = scaled_slope.denominator
        mate = total - 2
        if index <= 2 * mu and (2 % index == 0 or mate % index == 0):
            slope = scaled_slope / ell
            if slope > previous:
                yield mate, slope


def multiple_root_shifts(state):
    ell, mu, beta, previous = state
    upper = 1 - Fraction(beta, ell * mu)
    for degree in range(2, beta * mu + 2):
        slope = 1 - Fraction(beta * degree, ell * mu * (degree - 1))
        orbit = (ell * slope).denominator
        if not previous < slope < upper or (degree - 1) % orbit:
            continue
        for zero_order in range(mu):
            if (mu - zero_order) % orbit:
                continue
            for parts in partitions((mu - zero_order) // orbit, mu):
                multiplicities = parts + ((zero_order,) if zero_order else ())
                count = orbit * len(parts) + bool(zero_order)
                if not 2 <= count <= degree:
                    continue
                if any(degree * m == mu for m in multiplicities):
                    continue
                if all(degree * m % mu == 0 for m in multiplicities):
                    continue
                # Proposition 4.3 permits any root of maximal multiplicity.
                # Choose zero when the maximum is attained there.
                root_order = max(multiplicities)
                new_ell = ell if zero_order == root_order else lcm(ell, slope.denominator)
                new_x = mu - Fraction(beta, ell) - slope * (mu - root_order)
                new_beta = new_ell * (root_order - new_x)
                assert new_beta > 0 and new_beta.denominator == 1
                yield degree, orbit, slope, (new_ell, root_order, int(new_beta), slope)


def check_six_thirty_successors():
    F = Fraction
    previous = F(1, 6)
    incoming = [(6, mu, 5 * mu - 6, previous) for mu in range(2, 6)]
    assert all(not list(scalar_edges(state)) for state in incoming)
    transitions = {
        (state[1], degree, orbit, new_state)
        for state in incoming
        for degree, orbit, _, new_state in multiple_root_shifts(state)
    }
    assert transitions == {
        (2, 3, 1, (6, 1, 1, F(1, 2))),
        (2, 5, 2, (12, 1, 3, F(7, 12))),
        (3, 4, 1, (6, 1, 1, F(1, 3))),
        (3, 4, 1, (6, 2, 5, F(1, 3))),
        (3, 7, 2, (6, 1, 2, F(5, 12))),
        (3, 10, 3, (18, 1, 7, F(4, 9))),
        (4, 22, 3, (6, 1, 3, F(7, 18))),
    }
    next_states = {state for _, _, _, state in transitions}
    assert all(not list(multiple_root_shifts(state)) for state in next_states)
    scalar = {edge for state in next_states for edge in scalar_edges(state)}
    assert scalar == {(3, F(7, 12)), (5, F(17, 24)), (5, F(5, 12)), (13, F(7, 12))}


def check_six_thirty_root_field():
    def double_root_cubic(c):
        return 64*c**3 - 216*c**2 + 324*c - 81

    def two_simple_roots_cubic(k):
        return 448*k**3 - 1323*k**2 + 3888*k - 3888

    # Compare two polynomials of degree at most three.
    for c in range(4):
        cleared = (
            78*c*(8*c+3)
            - (2*c-3)*(8*c+3)*(4*c+9)
            - 84*c*(4*c+9)
        )
        assert cleared == -double_root_cubic(c)

    # The resultant with the monic quadratic in h is the norm of
    # r*h+s modulo h**2 + coefficient*h + constant.
    # Both polynomials in k have degree at most four.
    for k in range(5):
        coefficient, constant = -4*k-9, 15*k
        r, s = 45*k+144, 16*k*k-261*k
        resultant = s*s - coefficient*r*s + constant*r*r
        assert resultant == 7*k*two_simple_roots_cubic(k)

    assert 64 % 11 and 448 % 11
    for residue in range(11):
        assert double_root_cubic(residue) % 11
        assert two_simple_roots_cubic(residue) % 11


def check_resultant_successor_denominators():
    exceptional = set()
    for multiplicity in range(2, 6):
        # The cancelled gcd divides 5*m - 6. If D <= 2*m then
        # N*m - 1 <= 2*m*(5*m - 6), so this loop covers every N.
        limit = 2*(5*multiplicity - 6)
        for total in range(5, limit + 1, 2):
            slope_on_lattice = Fraction(total*(6+multiplicity) - 6,
                                        total*multiplicity - 1)
            denominator = slope_on_lattice.denominator
            if denominator <= 2*multiplicity:
                mate = total - 2
                exceptional.add((multiplicity, mate, denominator))
                assert gcd(denominator, multiplicity) == 1
                assert 2 % denominator and mate % denominator
    assert exceptional == {(4, 7, 5), (5, 21, 6)}


if __name__ == "__main__":
    check_first_coordinate_six()
    check_six_thirty_successors()
    check_six_thirty_root_field()
    check_resultant_successor_denominators()
    print("Checked the crossings, scalar outcomes, cubic eliminants and resultant successor denominators.")
