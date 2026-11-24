"""Check the asymmetric FullOC witness and its sparse-root exclusion."""

from fractions import Fraction
from itertools import combinations, permutations, product


INDICES = tuple(range(1, 7))
LETTERS = "abc"
SUBSETS = tuple(
    tuple(i for i in INDICES if mask & (1 << (i - 1))) for mask in range(64)
)
MODULUS_SETS = tuple(subset for subset in SUBSETS if 1 not in subset)
LEGACY_CAPS = {
    (2, 3): (2, 3),
    (3, 4): (3, 4),
    (2, 5): (5,),
    (3, 5): (5,),
    (4, 5): (5,),
}
DENOMINATOR = 4 * 1464
ZERO = (0, 0)
ONE = (DENOMINATOR, 0)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def values(bases, slopes):
    return tuple((4 * base, 4 * slope) for base, slope in zip(bases, slopes))


def add(*terms):
    return sum(term[0] for term in terms), sum(term[1] for term in terms)


def multiply(factor, term):
    return factor * term[0], factor * term[1]


def subtract(left, right):
    return left[0] - right[0], left[1] - right[1]


def half(term):
    require(term[0] % 2 == term[1] % 2 == 0, "inexact half")
    return term[0] // 2, term[1] // 2


def dominates(left, right):
    """Check an affine inequality throughout 0 <= u <= 1 at its endpoints."""
    constant, slope = subtract(left, right)
    return constant >= 0 and constant + slope >= 0


def data():
    a = values((80, 62, 198, 142, 0, 0, 496), (-10, 4, 0, -10, 0, 0, -10))
    b = values((67, 68, 198, 130, 0, 0, 484), (3, -2, 0, 2, 0, 0, 2))
    coordinates = {s: dict(zip(INDICES, v[:6])) for s, v in zip(LETTERS, (a, b, b))}
    totals = {s: v[-1] for s, v in zip(LETTERS, (a, b, b))}
    fourier = dict(zip(LETTERS, values((336, 198, 198), (-12, 0, 0))))
    exponent, = values((854,), (-2,))
    return coordinates, totals, fourier, exponent


def check_root(coordinates, totals, fourier, exponent):
    total = add(*totals.values())
    require(total == (DENOMINATOR, -24), "budget must equal 1-u/244")
    for s in LETTERS:
        require(all(dominates(value, ZERO) for value in coordinates[s].values()),
                ("coordinate nonnegativity", s))
        require(dominates(fourier[s], ZERO), ("auxiliary nonnegativity", s))
        tail = subtract(totals[s], add(*coordinates[s].values()))
        require(dominates(tail, ZERO), ("tail nonnegativity", s))
        weighted = add(multiply(7, tail), *(
            multiply(i, coordinates[s][i]) for i in INDICES
        ))
        require(weighted == ONE, ("weighted height", s))
        for j in (2, 3, 5):
            grouped = add(*(coordinates[s][i] for i in INDICES if i % j == 0))
            require(dominates(fourier[s], grouped), ("Fourier", s, j))
    require(dominates(exponent, ZERO), "exponent nonnegativity")
    for a, b in combinations(LETTERS, 2):
        require(dominates(add(totals[a], totals[b]), exponent), ("pair", a, b))
    require(dominates(multiply(4, total), add(multiply(6, exponent), *fourier.values())),
            "summed Fourier")
    order = [add(coordinates[s][4], multiply(2, coordinates[s][5])) for s in LETTERS]
    require(all(dominates(left, right) for left, right in zip(order, order[1:])),
            "symmetry order")


def check_geometry(coordinates, totals, exponent):
    total = add(*totals.values())
    tabulated = [
        [
            (add(*(coordinates[s][i] for i in subset)),
             add(*(multiply(i - 1, coordinates[s][i]) for i in subset)))
            for subset in SUBSETS
        ]
        for s in LETTERS
    ]
    checked = 0
    for a, b, c in product(*tabulated):
        unweighted = add(a[0], b[0], c[0])
        weighted = add(a[1], b[1], c[1])
        # Either affine branch globally lower-bounds the maximum.
        rhs = max(subtract(total, unweighted), subtract(add(total, weighted), ONE))
        require(dominates(rhs, exponent), ("geometry", a, b, c))
        checked += 1
    require(checked == 64**3, "missing geometry instances")
    return checked


def check_collisions(coordinates, totals, exponent):
    tail = {s: subtract(totals[s], add(*coordinates[s].values())) for s in LETTERS}
    checked = 0
    for a, b in permutations(LETTERS, 2):
        c = next(s for s in LETTERS if s not in (a, b))
        modulus = [
            (add(*(coordinates[c][i] for i in subset)),
             add(*(multiply(i, coordinates[c][i]) for i in subset)), bool(subset))
            for subset in MODULUS_SETS
        ]
        for j, k in combinations(INDICES, 2):
            p, q = coordinates[a][j], coordinates[a][k]
            rest = subtract(subtract(totals[a], tail[a]), add(p, q))
            for subset in SUBSETS[1:]:
                mass = add(*(coordinates[b][i] for i in subset))
                frozen = subtract(totals[b], mass)
                released_height = add(*(multiply(i, coordinates[b][i]) for i in subset))
                frozen_height = subtract(ONE, released_height)
                # The lexicographic maximizers are affine lower bounds on
                # these positive parts, including if a sign changes inside.
                loss_p = max(ZERO, subtract(multiply(2, p), frozen_height))
                loss_q = max(ZERO, subtract(multiply(2, q), frozen_height))
                energy = add(multiply(2, rest), p, q, half(add(loss_p, loss_q)))
                differences = [mass]
                differences.extend(
                    multiply(2, add(*(coordinates[b][i] for i in subset if i % e)))
                    for e in range(2, 7)
                )
                differences.extend(multiply(2, coordinates[b][i])
                                   for i in LEGACY_CAPS.get(subset, ()))
                diagonal = add(tail[a], frozen, totals[c])
                for frozen_c, modulus_height, is_locked in modulus:
                    caps = differences
                    if is_locked:
                        caps = differences + [max(ZERO, subtract(released_height, modulus_height))]
                    common = add(tail[a], frozen, half(add(totals[c], frozen_c, energy)))
                    # Prove the collision inequality for every candidate cap
                    # and then take the minimum inside the outer maximum.
                    for cap in caps:
                        rhs = max(diagonal, add(common, half(cap)))
                        require(dominates(rhs, exponent), ("collision", a, b, j, k, subset))
                    checked += 1
    require(checked == 6 * 15 * 63 * 32, "missing collision instances")
    return checked


def check_sparse_cut(coordinates, totals, exponent):
    # A=b, B=c, C=a, with G1={2,4}, G2={3}, I_B={2,4}, and O_C empty.
    squares, cube = (2, 4), (3,)
    released_mass = add(*(coordinates["c"][i] for i in squares))
    released_height = add(*(multiply(i, coordinates["c"][i]) for i in squares))
    frozen_height = subtract(ONE, released_height)
    first_height = add(*(multiply(i // 2, coordinates["b"][i]) for i in squares))
    second_height = coordinates["b"][cube[0]]
    require(dominates(frozen_height, multiply(2, first_height)), "square penalty")
    require(dominates(frozen_height, multiply(2, second_height)), "cube penalty")
    energy = add(multiply(2, coordinates["b"][1]),
                 *(coordinates["b"][i] for i in squares + cube))
    frozen_a = subtract(totals["b"], add(*coordinates["b"].values()))
    frozen_b = subtract(totals["c"], released_mass)
    require(dominates(energy, totals["a"]), "off-diagonal dominance")
    # Every released B index is divisible by 2, so its difference cap is zero.
    bound = add(frozen_a, frozen_b, half(add(totals["a"], energy)))
    require(bound == (3280, -4), "unexpected sparse exponent")
    excess = subtract(exponent, bound)
    require(excess[0] > 0 and excess[0] + excess[1] > 0, "sparse row did not cut")
    return Fraction(excess[0] + excess[1], DENOMINATOR)


def main():
    coordinates, totals, fourier, exponent = data()
    check_root(coordinates, totals, fourier, exponent)
    geometry = check_geometry(coordinates, totals, exponent)
    collisions = check_collisions(coordinates, totals, exponent)
    excess = check_sparse_cut(coordinates, totals, exponent)
    print(f"verified {geometry} geometry and {collisions} FullOC collision instances")
    print(f"asymmetric witness valid on 243/244 <= lambda <= 1, sparse gap >= {excess}")


if __name__ == "__main__":
    main()
