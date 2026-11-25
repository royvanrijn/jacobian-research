"""Check the symmetric obstructions for ideal residue energies and a partial modulus."""

from fractions import Fraction
from functools import cache
from itertools import combinations, product
from math import gcd


def check_profile(scale, numerators):
    support = tuple(numerators)
    masks = tuple(mask for size in range(len(support) + 1)
                  for selected in combinations(support, size)
                  for mask in (frozenset(selected),))
    nonempty = masks[1:]
    masses = {mask: sum(numerators[degree] for degree in mask) for mask in masks}
    heights = {mask: sum(degree * numerators[degree] for degree in mask) for mask in masks}
    full = frozenset(support)
    total = masses[full]
    if 3 * total != scale or heights[full] != scale:
        raise ValueError("expected unit budget and three height-one letters")

    def ideal_energy(modulus):
        return 2 * total - min(total, modulus)

    @cache
    def difference(released):
        return min(masses[released], *(
            2 * sum(numerators[i] for i in released if i % degree)
            for degree in range(2, max(support) + 1)
        ))

    @cache
    def target_cost(released):
        choices = [total, difference(released)]
        for locked in masks:
            fixed, modulus = masses[locked], heights[locked]
            choices.append(fixed + max(0, heights[released] - modulus))
            for group in nonempty:
                if group <= released:
                    root = Fraction(heights[group], gcd(*group))
                    choices.append(fixed + masses[released] - masses[group] + max(0, root - modulus))
        return min(choices)

    def collision(released):
        locked = full - released
        frozen, modulus = masses[locked], heights[locked]
        off_diagonal = frozen + Fraction(total + ideal_energy(modulus) + target_cost(released), 2)
        return max(frozen + total, off_diagonal)

    def shift(left, right):
        frozen_left, frozen_right = full - left, full - right
        diagonal = masses[frozen_left] + masses[frozen_right] + total
        shift_height = min(heights[left] - heights[frozen_right],
                           heights[right] - heights[frozen_left],
                           scale - heights[frozen_left] - heights[frozen_right])
        if shift_height < 0:
            return diagonal
        off_diagonal = diagonal - total + Fraction(
            total + shift_height + difference(left) + target_cost(right), 2
        )
        return max(diagonal, off_diagonal)

    def rational(left, right, locked):
        combined_height = Fraction(heights[left] + heights[right], gcd(*left, *right))
        frozen = 2 * total - masses[left] - masses[right] + masses[locked]
        return frozen + max(0, combined_height - heights[locked])

    fourier_mass = max(sum(numerators[i] for i in support if i % d == 0)
                       for d in range(2, max(support) + 1))
    bounds = {
        "pair": 2 * total,
        "Fourier": Fraction(4 * scale - 3 * fourier_mass, 6),
        "geometry": min(max(scale - sum(masses[mask] for mask in group),
                            sum(heights[mask] - masses[mask] for mask in group))
                        for group in product(masks, repeat=3)),
        "direct energy": min(masses[mask] + ideal_energy(heights[mask]) for mask in masks),
        "collision energy": min(map(collision, nonempty)),
        "paired shift": min(shift(a, b) for a, b in product(nonempty, repeat=2)),
        "partial-modulus rational": min(rational(a, b, lock)
                                       for a, b, lock in product(nonempty, nonempty, masks)),
    }
    return {family: Fraction(value, scale) for family, value in bounds.items()}


def main():
    first = check_profile(150, {1: 19, 4: 30, 11: 1})
    second = check_profile(12, {1: 1, 3: 1, 4: 2})
    old_families = first.keys() - {"partial-modulus rational"}
    if min(first[family] for family in old_families) != Fraction(8, 15):
        raise ValueError("incorrect ideal-energy obstruction")
    if first["partial-modulus rational"] != Fraction(7, 15):
        raise ValueError("incorrect rational cut on the ideal-energy obstruction")
    if min(second.values()) != Fraction(1, 2):
        raise ValueError("incorrect joint obstruction")
    for numbers, bounds in (("19,30,1 at degrees 1,4,11", first),
                            ("1,1,2 at degrees 1,3,4", second)):
        print(numbers)
        for family, exponent in bounds.items():
            print(f"  {family} {exponent}")


if __name__ == "__main__":
    main()
