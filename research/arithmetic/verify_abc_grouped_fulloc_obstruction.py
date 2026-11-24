"""Replay the symmetric 4/7 obstruction for grouped abc box bounds."""

from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import gcd, lcm


DENOMINATOR = 2100
LETTER_TOTAL = 700
EXPONENT = 1200
LIFT = {1: 168, 3: 300, 4: 200, 5: 25, 15: 5, 16: 2}
TOTAL = 3 * LETTER_TOTAL


def require(condition, message):
    if not condition:
        raise ValueError(message)


def geometry_minimum(indices):
    """Count how many letters select each degree in the symmetric geometry rows."""
    best = TOTAL
    for multiplicities in product(range(4), repeat=len(indices)):
        mass = sum(LIFT[i] * m for i, m in zip(indices, multiplicities, strict=True))
        excess = sum((i - 1) * LIFT[i] * m
                     for i, m in zip(indices, multiplicities, strict=True))
        best = min(best, max(TOTAL - mass, TOTAL - DENOMINATOR + excess))
    return Fraction(best, DENOMINATOR)


def collision_minimum(cutoff, indices):
    """Minimize the collision upper bound over both root groups and all locks."""
    masks = range(1 << len(indices))
    subsets = [tuple(i for bit, i in enumerate(indices) if mask & (1 << bit))
               for mask in masks]
    masses = [sum(LIFT[i] for i in group) for group in subsets]
    weights = [sum(i * LIFT[i] for i in group) for group in subsets]
    scale = lcm(*indices)
    root_heights = [scale * weight // gcd(*group) if group else 0
                    for weight, group in zip(weights, subsets, strict=True)]
    tail = LETTER_TOTAL - masses[-1]
    weighted_tail = (cutoff + 1) * tail
    require(tail >= 0 and weights[-1] + weighted_tail <= DENOMINATOR,
            ("height or tail row", cutoff))
    first_groups = [(left, right)
                    for left, right in combinations_with_replacement(masks, 2)
                    if not left & right]
    divisible = [sum(1 << bit for bit, i in enumerate(indices) if i % degree == 0)
                 for degree in range(2, cutoff + 1)]
    fourier = max(masses[mask] for mask in divisible)
    require(6 * EXPONENT + 3 * fourier <= 4 * TOTAL, ("Fourier row", cutoff))

    best = None
    for released, mass in enumerate(masses):
        frozen = LETTER_TOTAL - mass
        modulus = scale * (weights[-1] - weights[released] + weighted_tail)
        # Energy is doubled to clear the halves in its two height penalties.
        twice_energy, left, right = min(
            (2 * scale * (2 * masses[-1] - masses[left] - masses[right])
             + max(0, 2 * root_heights[left] - modulus)
             + max(0, 2 * root_heights[right] - modulus), left, right)
            for left, right in first_groups
        )
        difference_caps = [2 * (mass - masses[released & mask]) for mask in divisible]
        constant_cap = scale * min(mass, *difference_caps)
        subgroup_locks = [chosen for chosen in masks if chosen and chosen & released == chosen]
        diagonal = 4 * scale * (tail + frozen + LETTER_TOTAL)
        for locked, locked_mass in enumerate(masses):
            third_modulus = scale * weights[locked]
            value_lock = max(0, scale * weights[released] - third_modulus)
            subgroup_caps = [
                scale * (mass - masses[chosen]) + max(0, root_heights[chosen] - third_modulus)
                for chosen in subgroup_locks
            ]
            cap = min(constant_cap, value_lock, *subgroup_caps)
            off_diagonal = (4 * scale * (tail + frozen)
                            + 2 * scale * (LETTER_TOTAL + locked_mass)
                            + twice_energy + 2 * cap)
            bound = max(diagonal, off_diagonal)
            if best is None or bound < best[0]:
                best = bound, released, locked, left, right
    bound, released, locked, left, right = best
    return (Fraction(bound, 4 * scale * DENOMINATOR),
            tuple(subsets[i] for i in (released, locked, left, right)))


def main():
    require(TOTAL == DENOMINATOR and 0 <= EXPONENT <= 2 * LETTER_TOTAL,
            "budget and pair rows")
    require(sum(LIFT.values()) == LETTER_TOTAL
            and sum(i * weight for i, weight in LIFT.items()) == DENOMINATOR,
            "exact height of the lift")
    expected = {
        6: Fraction(4, 7),
        7: Fraction(4, 7),
        14: Fraction(4, 7),
        15: Fraction(2327, 4200),
        16: Fraction(2293, 4200),
    }
    for cutoff, wanted in expected.items():
        support = tuple(i for i in LIFT if i <= cutoff)
        geometry = geometry_minimum(support)
        collision, parameters = collision_minimum(cutoff, support)
        require(geometry == Fraction(4, 7), ("geometry minimum", cutoff, geometry))
        require(collision == wanted, ("collision minimum", cutoff, collision))
        print(f"cutoff={cutoff} geometry={geometry} collision={collision} "
              f"(release,lock,group1,group2)={parameters}")


if __name__ == "__main__":
    main()
