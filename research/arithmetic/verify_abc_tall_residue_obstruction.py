"""Replay an asymmetric obstruction to the height-tagged abc bounds."""

from dataclasses import dataclass
from fractions import Fraction
from functools import cache
from itertools import combinations, permutations, product
from math import gcd


DENOMINATOR = 95920
# Clear all possible denominators in early grouped root heights.
ROOT_SCALE = 60
HEIGHT = ROOT_SCALE * DENOMINATOR
OBSTRUCTION = Fraction(2627, 4796)
MASSES = (
    (4508, 5169, 9221, 11550, 1266, 0, 0),
    (4504, 5380, 9552, 11335, 1332, 0, 0),
    (4504, 5380, 9552, 11335, 1332, 0, 0),
)
TAIL = 1 << 6
ALL = (TAIL << 1) - 1
MASKS = range(ALL + 1)
NONEMPTY = range(1, ALL + 1)
SUBSETS = tuple(tuple(i for i in range(7) if mask & (1 << i)) for mask in MASKS)
SUBMASKS = tuple(tuple(group for group in NONEMPTY if group & mask == group) for mask in MASKS)
GROUP_PAIRS = tuple((left, right) for left, right in combinations(NONEMPTY, 2)
                    if not left & right)
DIVISIBLE = tuple(sum(1 << i for i in range(6) if (i + 1) % degree == 0)
                  for degree in (2, 3, 4, 5, 6))


@dataclass(frozen=True)
class Letter:
    total: int
    mass: tuple[int, ...]
    high: tuple[int, ...]
    low: tuple[int, ...]
    root: tuple[int, ...]
    difference: tuple[int, ...]
    fourier: int


def require(condition, message):
    if not condition:
        raise ValueError(message)


def projection(numerators, tall):
    values = [ROOT_SCALE * value for value in numerators]
    mass = tuple(sum(values[i] for i in subset) for subset in SUBSETS)
    early = tuple(sum((i + 1) * values[i] for i in subset if i < 6)
                  for subset in SUBSETS)
    high = tuple(weight + (HEIGHT - early[ALL] if mask & TAIL else 0)
                 for mask, weight in enumerate(early))
    low = high if tall else tuple(weight + (7 * values[6] if mask & TAIL else 0)
                                  for mask, weight in enumerate(early))
    root = tuple(weight if mask == 0 or mask & TAIL
                 else weight // gcd(*(i + 1 for i in subset))
                 for mask, (weight, subset) in enumerate(zip(high, SUBSETS, strict=True)))
    difference = tuple(min(mass[mask], *(2 * (mass[mask] - mass[mask & multiples])
                                         for multiples in DIVISIBLE))
                       for mask in MASKS)
    require(early[ALL] + 7 * values[6] <= HEIGHT, "early height and tail mass must fit the height budget")
    return Letter(mass[ALL], mass, high, low, root, difference,
                  max(mass[multiples] for multiples in DIVISIBLE))


@cache
def twice_energy(letter, modulus):
    """Two times the smallest grouped-root energy exponent."""
    losses = [max(0, 2 * height - modulus) for height in letter.root]
    return min(4 * letter.total - 2 * (letter.mass[left] + letter.mass[right])
               + losses[left] + losses[right] for left, right in GROUP_PAIRS)


@cache
def lock_cost(free, target):
    """Minimize frozen target mass plus the conditional one-value cap."""
    root_locks = [min(mass + max(0, height - weight)
                      for mass, weight in zip(target.mass, target.low, strict=True))
                  for height in free.root]
    result = []
    for mask in MASKS:
        value_lock = min(mass + max(0, free.high[mask] - weight)
                         for mass, weight in zip(target.mass, target.low, strict=True))
        group_locks = (free.mass[mask] - free.mass[j] + root_locks[j]
                       for j in SUBMASKS[mask])
        result.append(min(free.difference[mask], target.total, value_lock, *group_locks))
    return tuple(result)


def geometry_minimum(letters, total):
    released = [tuple(zip(letter.mass, letter.high, strict=True)) for letter in letters]
    return min(max(total - (m_a + m_b + m_c),
                   total - HEIGHT + (w_a + w_b + w_c) - (m_a + m_b + m_c))
               for (m_a, w_a), (m_b, w_b), (m_c, w_c) in product(*released))


def collision_minimum(letters):
    """Return four times the minimum across ordered collision choices."""
    best = 12 * HEIGHT
    for a, b, c in permutations(letters):
        caps = lock_cost(b, c)
        for released in NONEMPTY:
            frozen = b.total - b.mass[released]
            modulus = b.low[ALL ^ released]
            diagonal = 4 * (frozen + c.total)
            off_diagonal = (4 * frozen + 2 * (c.total + caps[released])
                            + twice_energy(a, modulus))
            best = min(best, max(diagonal, off_diagonal))
    return best


def residue_minimum(letters):
    """Return four times the minimum of the direct residue-energy bounds."""
    return min(4 * c.mass[mask] + twice_energy(a, c.low[mask]) + twice_energy(b, c.low[mask])
               for a, b, c in permutations(letters) for mask in NONEMPTY)


def shift_minimum(letters):
    """Return four times the minimum of the paired shift bounds."""
    best = 12 * HEIGHT
    for a, b, c in permutations(letters):
        caps = lock_cost(b, c)
        for left, right in product(NONEMPTY, repeat=2):
            modulus_a = a.low[ALL ^ left]
            modulus_b = b.low[ALL ^ right]
            shift = min(a.high[left] - modulus_b, b.high[right] - modulus_a,
                        HEIGHT - modulus_a - modulus_b)
            frozen = a.total + b.total - a.mass[left] - b.mass[right]
            diagonal = 4 * (frozen + c.total)
            if shift < 0:
                bound = diagonal
            else:
                off_diagonal = 4 * frozen + 2 * (c.total + shift + a.difference[left] + caps[right])
                bound = max(diagonal, off_diagonal)
            best = min(best, bound)
    return best


def main():
    letters = tuple(projection(row, tall=index != 0) for index, row in enumerate(MASSES))
    total = sum(letter.total for letter in letters)
    require(total == HEIGHT, "radical budget")
    require(all(letter.low[ALL] == HEIGHT for letter in letters[1:]), "tall letter heights")
    require(all(a[3] + 2 * a[4] >= b[3] + 2 * b[4] for a, b in zip(MASSES, MASSES[1:])),
            "six-index symmetry ordering")
    bounds = {
        "pair": Fraction(total - max(letter.total for letter in letters), HEIGHT),
        "Fourier": Fraction(4 * total - sum(letter.fourier for letter in letters), 6 * HEIGHT),
        "geometry": Fraction(geometry_minimum(letters, total), HEIGHT),
        "grouped collisions": Fraction(collision_minimum(letters), 4 * HEIGHT),
        "direct residue": Fraction(residue_minimum(letters), 4 * HEIGHT),
        "paired shift": Fraction(shift_minimum(letters), 4 * HEIGHT),
    }
    require(min(bounds.values()) == OBSTRUCTION, "unexpected program minimum")
    require(all(bounds[name] == OBSTRUCTION
                for name in ("geometry", "grouped collisions", "direct residue")),
            "unexpected binding families")
    for name, bound in bounds.items():
        print(f"{name} minimum: {bound}")
    partial = ROOT_SCALE * 5775
    require(all(2 * letter.root[1 << index] <= 4 * partial
                for letter in letters[:2] for index in (2, 3)),
            "partial modulus must extinguish all root-height penalties")
    divisor_bound = Fraction(2 * sum(letter.total for letter in letters[:2]) + 2 * partial
                             - sum(letter.mass[0b1100] for letter in letters[:2]), 2 * HEIGHT)
    require(divisor_bound == Fraction(48763, DENOMINATOR), "partial-divisor bound")
    print(f"obstruction: {OBSTRUCTION}, partial-divisor bound: {divisor_bound}")


if __name__ == "__main__":
    main()
