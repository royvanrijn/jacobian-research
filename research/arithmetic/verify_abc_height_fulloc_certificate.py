"""Replay the abc exponent certificate with the two tall terms distinguished."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import gzip
from itertools import permutations, product
import json
import lzma
from math import gcd, lcm
from pathlib import Path
from types import ModuleType

from verify_abc_pair_deviation_certificate import (
    CHECKER_FILE, CHECKER_SHA256, require, source_bytes,
)


CERTIFICATE_FILE = Path(__file__).with_name("abc_height_fulloc_certificate.json.gz")
BUDGET_RETARGET_FILES = (
    (Fraction(91, 100), "abc_budget_910_delta.json.xz"),
    (Fraction(23, 25), "abc_budget_920_delta.json.xz"),
    (Fraction(37, 40), "abc_budget_925_delta.json.xz"),
    (Fraction(93, 100), "abc_budget_930_delta.json.xz"),
    (Fraction(187, 200), "abc_budget_935_delta.json.xz"),
    (Fraction(47, 50), "abc_budget_940_delta.json.xz"),
    (Fraction(189, 200), "abc_budget_945_delta.json.xz"),
    (Fraction(19, 20), "abc_budget_950_delta.json.xz"),
    (Fraction(191, 200), "abc_budget_955_delta.json.xz"),
    (Fraction(24, 25), "abc_budget_960_delta.json.xz"),
    (Fraction(193, 200), "abc_budget_965_delta.json.xz"),
    (Fraction(97, 100), "abc_budget_970_delta.json.xz"),
    (Fraction(39, 40), "abc_budget_975_delta.json.xz"),
    (Fraction(391, 400), "abc_budget_9775_delta.json.xz"),
    (Fraction(49, 50), "abc_budget_980_delta.json.xz"),
    (Fraction(1961, 2000), "abc_budget_9805_delta.json.xz"),
    (Fraction(981, 1000), "abc_budget_981_delta.json.xz"),
    (Fraction(393, 400), "abc_budget_9825_delta.json.xz"),
    (Fraction(983, 1000), "abc_budget_983_delta.json.xz"),
    (Fraction(123, 125), "abc_budget_984_delta.json.xz"),
    (Fraction(197, 200), "abc_budget_985_delta.json.xz"),
    (Fraction(493, 500), "abc_budget_986_delta.json.xz"),
)
UNIT_CERTIFICATE_PARTS = (
    "abc_height_fulloc_unit_certificate.part1.xz",
    "abc_height_fulloc_unit_certificate.part2.xz",
    "abc_height_fulloc_unit_certificate.part3.xz",
)
IMPROVED_UNIT_CERTIFICATE_PARTS = (
    "abc_height_fulloc_improved_unit_certificate.part1.xz",
    "abc_height_fulloc_improved_unit_certificate.part2.xz",
    "abc_height_fulloc_improved_unit_certificate.part3.xz",
    "abc_height_fulloc_improved_unit_certificate.part4.xz",
)
LETTERS = ("a", "b", "c")
EARLY = frozenset(range(1, 7))
TAIL = 7
ALL = EARLY | {TAIL}

# Each window specifies (multiplication, isogeny degree, height lower bound, moment).
MORDELL_HEIGHT_WINDOWS = (
    (1, 3, Fraction(7, 10), Fraction(1)),
    (2, 1, Fraction(3, 5), Fraction(13, 10)),
    (3, 1, Fraction(1, 4), Fraction(2)),
    (2, 3, Fraction(7, 40), Fraction(23, 10)),
    (4, 1, Fraction(2, 15), Fraction(13, 5)),
    (5, 1, Fraction(1, 10), Fraction(3)),
    (3, 3, Fraction(3, 40), Fraction(3)),
    (6, 1, Fraction(1, 16), Fraction(33, 10)),
    (7, 1, Fraction(1, 20), Fraction(18, 5)),
    (9, 1, Fraction(1, 40), Fraction(4)),
    (11, 1, Fraction(1, 60), Fraction(22, 5)),
    (13, 1, Fraction(1, 80), Fraction(47, 10)),
    (17, 1, Fraction(1, 140), Fraction(26, 5)),
)


# Each pair specifies a retained residual weight threshold and a 3-torsion moment.
MORDELL_LOCAL_WINDOWS = (
    (Fraction(1, 9), Fraction(1)),
    (Fraction(17, 200), Fraction(5, 4)),
    (Fraction(2, 25), Fraction(4, 3)),
    (Fraction(13, 200), Fraction(3, 2)),
    (Fraction(1, 27), Fraction(2)),
)


def mordell_moment_loss(k):
    if k <= 4:
        return (5 * k - 5) / 18
    return Fraction(5, 6) + Fraction(1597, 5000) * (k - 4)


def subset(values, allow_empty=False):
    require(isinstance(values, list) and (values or allow_empty)
            and all(type(value) is int and value in ALL for value in values)
            and values == sorted(set(values)), "invalid coordinate subset")
    return frozenset(values)


def load_checker():
    checker = ModuleType("zeng_fulloc_rows")
    exec(compile(source_bytes(CHECKER_FILE, CHECKER_SHA256),
                 str(CHECKER_FILE), "exec"), checker.__dict__)
    return checker


def initial_rows(checker, budget):
    # The short term is labelled a. Only the two tall terms are sorted.
    root = checker.base(budget)
    return root[:-2] + root[-1:]



def tail_fourier_rows(checker):
    """Use exact affine bounds on the unexposed Fourier mass."""
    rows = []
    primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
    first_hidden = {2: 8, 3: 9, 5: 25}
    for letter in LETTERS:
        for cutoff in (8, 9, 10, 12, 15, 18, 20, 24, 32):
            penalties = [(cutoff - first_hidden.get(q, q), q)
                         for q in primes if first_hidden.get(q, q) < cutoff]
            coefficients = checker.add(
                (cutoff, checker.var(letter)),
                *((i - cutoff, checker.var(f"{letter}{i}")) for i in EARLY),
                (-sum(weight for weight, _ in penalties), checker.var(f"F{letter}")),
                *((weight, checker.var(f"{letter}{i}"))
                  for weight, q in penalties for i in EARLY if i % q == 0),
            )
            rows.append(checker.row(coefficients, 1))
    return rows

def projected_signature_rows(checker, signatures):
    """Keep the provably removable mass in the six exact valuation classes."""
    rows = []
    for signature in signatures:
        terms = [(1, checker.var("D"))]
        for letter, divisor in zip(LETTERS, signature, strict=True):
            terms.append((-1, checker.var(letter)))
            if divisor is not None:
                terms.extend((1, checker.var(f"{letter}{degree}"))
                             for degree in range(divisor, TAIL, divisor))
        rows.append(checker.row(checker.add(*terms), 0))
    return rows


def residual_signature_rows(checker, euclidean=False):
    signatures = {(3, 3, 3)} if euclidean else {(4, 4, 4), (5, 5, 5)}
    if euclidean:
        signatures.update(permutations((2, 4, 4)))
    else:
        signatures.update(permutations((2, 4, 8)))
        signatures.update(permutations((3, 3, 9)))
    return projected_signature_rows(checker, sorted(signatures))


def square_cube_rows(checker):
    # The third letter is retained in full.
    return projected_signature_rows(checker, permutations((2, 3, None)))


def square_fifth_rows(checker):
    """Charge a quarter of the possible odd valuation mass in the quadratic kernel."""
    rows = []
    for square, fifth, frozen in permutations(LETTERS):
        terms = [(1, checker.var("D")), (1, checker.var(f"{fifth}5"))]
        terms.extend((-1, checker.var(letter)) for letter in LETTERS)
        terms.extend((-Fraction(1, 4), checker.var(letter)) for letter in (square, frozen))
        terms.extend((scale, checker.var(f"{letter}{degree}"))
                     for letter, scale in ((square, Fraction(5, 4)), (frozen, Fraction(1, 4)))
                     for degree in (2, 4, 6))
        rows.append(checker.row(checker.add(*terms), 0))
    return rows


def branch_rows(checker, node, signed=False, *, short="a"):
    """Return the disjunction of affine upper bounds for one arithmetic inequality."""
    require(short in LETTERS, "invalid short letter")
    zero = (checker.add(), Fraction(0))
    one = (zero[0], Fraction(1))

    def affine_sum(*terms):
        vector = checker.add(*[(scale, value[0]) for scale, value in terms])
        return vector, sum(scale * value[1] for scale, value in terms)

    def variable(name):
        return checker.var(name), Fraction(0)

    def row(affine):
        return checker.row(affine[0], -affine[1])

    def mass(letter, indices):
        if TAIL in indices:
            return affine_sum((1, variable(letter)), (-1, mass(letter, EARLY - indices)))
        return affine_sum(*[(1, variable(f"{letter}{i}")) for i in indices])

    def weighted_early(letter, indices, divisor=1):
        return affine_sum(*[(Fraction(i, divisor), variable(f"{letter}{i}"))
                            for i in indices if i in EARLY])

    def upper_weight(letter, indices):
        if TAIL in indices:
            return affine_sum((1, one), (-1, weighted_early(letter, EARLY - indices)))
        return weighted_early(letter, indices)

    def lower_weight(letter, indices):
        if TAIL not in indices:
            return weighted_early(letter, indices)
        if letter != short:
            return upper_weight(letter, indices)
        return affine_sum((1, weighted_early(letter, indices)), (7, mass(letter, {TAIL})))

    def root_height(letter, group):
        # TAIL is an aggregate, so a group containing it is used as a first power.
        if TAIL in group:
            return upper_weight(letter, group)
        return weighted_early(letter, group, gcd(*group))

    def max_rows(base, terms):
        branches = []
        for flags in product((0, 1), repeat=len(terms)):
            pairs = tuple(zip(flags, terms, strict=True))
            bound = row(affine_sum((1, base), *[
                (-scale, value) for selected, (scale, value) in pairs if selected
            ]))
            inequalities = [bound]
            if signed:
                inequalities.extend(row(affine_sum((-1 if selected else 1, value)))
                                    for selected, (_, value) in pairs)
            branches.append(inequalities)
        return branches

    def dominant_rows(base, alternatives):
        return [[row(affine_sum((1, base), (-1, bound)))] + [
            row(affine_sum((1, other), (-1, bound)))
            for index, other in enumerate(alternatives) if index != chosen
        ] for chosen, bound in enumerate(alternatives)]

    objective = variable("D")
    kind = node["type"]
    total = affine_sum(*[(1, variable(letter)) for letter in LETTERS])

    if kind in ("frozen-cofactor-plane", "frozen-cofactor-plane-strip"):
        require(signed and set(node) == {"type", "A", "B", "I", "J"},
                "unexpected frozen-cofactor plane metadata")
        a, b = node["A"], node["B"]
        require(a in LETTERS and b in LETTERS and a != b, "invalid plane letters")
        left, right = subset(node["I"]), subset(node["J"])
        require(left <= EARLY and gcd(*left) > 1, "invalid plane power group")
        released = affine_sum((1, mass(a, left)), (1, mass(b, right)))
        base = affine_sum((1, objective), (-1, total), (1, released))
        if kind == "frozen-cofactor-plane-strip":
            require(a == short, "the plane strip requires releasing the short term")
            gap = affine_sum((1, root_height(b, right)), (1, weighted_early(a, left)), (-1, one))
        else:
            gap = affine_sum((1, root_height(b, right)), (-1, lower_weight(a, ALL - left)))
        return max_rows(base, [(Fraction(1, gcd(*left)), gap)])

    if kind in ("cubic-two-quartics-strip", "cubic-two-quartics-strip-elementary"):
        require(signed and set(node) == {"type"}, "unexpected diagonal strip metadata")
        a = short
        b, c = (letter for letter in LETTERS if letter != short)
        x, y, z = variable(f"{a}3"), variable(f"{b}4"), variable(f"{c}4")
        frozen = affine_sum((1, total), (-1, x), (-1, y), (-1, z))
        base = affine_sum((1, objective), (-1, frozen))
        area = affine_sum((3, x), (1, y), (1, z), (-1, one))
        positive = row(affine_sum((-1, area)))
        lines = (
            affine_sum((Fraction(5, 8), area)),
            affine_sum((Fraction(7, 24), area), (Fraction(1, 3), x)),
        )
        guards = [positive]
        if kind == "cubic-two-quartics-strip-elementary":
            guards.append(row(affine_sum((1, base), (-1, area))))
        return [[row(area), row(base)]] + [
            branch + guards for branch in dominant_rows(base, lines)
        ]

    if kind in ("cubic-two-quartics", "cubic-two-quartics-lock"):
        require(signed and set(node) == {"type", "A", "B"},
                "unexpected cubic and two quartics metadata")
        a, b = node["A"], node["B"]
        require(a in LETTERS and b in LETTERS and a != b, "invalid diagonal letters")
        c = next(letter for letter in LETTERS if letter not in (a, b))
        x, y, z = variable(f"{a}3"), variable(f"{b}4"), variable(f"{c}4")
        # Lead with z^4, giving sqrt(x*y)/2, and use AM–GM with ratio 9/10.
        cover = affine_sum((Fraction(9, 40), x), (Fraction(5, 18), y))
        if kind == "cubic-two-quartics-lock":
            cover = affine_sum((Fraction(15, 16), cover))
        frozen = affine_sum((1, total), (-1, x), (-1, y), (-1, z))
        base = affine_sum((1, objective), (-1, frozen))
        alternatives = (affine_sum((2, cover)), *(
            affine_sum((1, cover), (Fraction(1, 3), side)) for side in (x, y, z)
        ))
        branches = dominant_rows(base, alternatives)
        if kind == "cubic-two-quartics-lock":
            weight = lower_weight(c, ALL - {4})
            weak = row(affine_sum((1, weight), (-Fraction(1, 8), one)))
            strong = row(affine_sum((-1, weight), (Fraction(1, 8), one)))
            return [[weak]] + [branch + [strong] for branch in branches]
        return branches

    if kind == "cubic-quartic-affine":
        require(signed and set(node) == {"type", "A", "B"},
                "unexpected cubic–quartic metadata")
        a, b = node["A"], node["B"]
        require(a in LETTERS and b in LETTERS and a != b, "invalid surface letters")
        c = next(letter for letter in LETTERS if letter not in (a, b))
        x, y = variable(f"{a}3"), variable(f"{b}4")
        u = affine_sum((1, one), (-3, x))
        # Quadratic tangent x(1-3x) at x=27/100, then AM–GM with y/4.
        cover = affine_sum((Fraction(2187, 10000), one),
                          (-Fraction(31, 50), x), (Fraction(1, 4), y))
        frozen = affine_sum((1, variable(c)), (1, variable(b)), (-1, y))
        alternatives = (affine_sum((2, cover)), *(
            affine_sum((1, cover), (Fraction(1, 3), side)) for side in (u, x, y)
        ))
        base = affine_sum((1, objective), (-1, frozen))
        return dominant_rows(base, alternatives)
    if kind in ("offdiag", "offdiag-strip", "offdiag-ray", "offdiag-ray-strip",
                "offdiag-lattice", "offdiag-lattice-strip", "offdiag-arch-lattice-strip",
                "offdiag-asymmetric", "offdiag-asymmetric-strip",
                "offdiag-mixed-powers", "offdiag-mixed-powers-strip",
                "offdiag-valuation-rays", "offdiag-valuation-rays-strip"):
        valuation_rays = kind in ("offdiag-valuation-rays", "offdiag-valuation-rays-strip")
        mixed_powers = kind in ("offdiag-mixed-powers", "offdiag-mixed-powers-strip") or valuation_rays
        asymmetric = kind in ("offdiag-asymmetric", "offdiag-asymmetric-strip") or mixed_powers
        archimedean_strip = kind in ("offdiag-arch-lattice-strip", "offdiag-asymmetric-strip",
                                    "offdiag-mixed-powers-strip", "offdiag-valuation-rays-strip")
        use_lattice = kind in ("offdiag-lattice", "offdiag-lattice-strip") or asymmetric or archimedean_strip
        use_ray = kind in ("offdiag-ray", "offdiag-ray-strip")
        groups = {"I1", "I2", "J1", "J2"} if asymmetric else {"I", "J"}
        fields = {"type", "C", "A", "O"} | groups | ({"R"} if use_ray or use_lattice else set())
        require(signed and set(node) == fields, "unexpected off-diagonal metadata")
        c, a = node["C"], node["A"]
        require(c in LETTERS and a in LETTERS and c != a, "invalid collision letters")
        b = next(letter for letter in LETTERS if letter not in (a, c))
        if asymmetric:
            left = (subset(node["I1"]), subset(node["I2"]))
            right = (subset(node["J1"]), subset(node["J2"]))
        else:
            left = (subset(node["I"]),) * 2
            right = (subset(node["J"]),) * 2
        require(mixed_powers or all((first | second) <= EARLY and gcd(*(first | second)) > 1
                                    for first, second in (left, right)),
                "invalid power groups in the off-diagonal estimate")
        locked = subset(node["O"])
        lock_weight = lower_weight(c, locked)
        gap = affine_sum((1, one), (-1, lock_weight))
        if kind in ("offdiag-strip", "offdiag-ray-strip", "offdiag-lattice-strip") or archimedean_strip:
            require(c == short, "the strip estimate requires locking the short addend")
            gap = upper_weight(c, ALL - locked)
        ray_letters = node["R"] if use_ray or use_lattice else []
        require(isinstance(ray_letters, list) and (not (use_ray or use_lattice) or ray_letters)
                and ray_letters in ([], [a], [b], sorted((a, b))), "invalid ray letters")
        weak, strong = [], []

        def off_diagonal(letter, groups):
            first, second = groups
            frozen = affine_sum((2, variable(letter)),
                                (-1, mass(letter, first)), (-1, mass(letter, second)))
            if letter not in ray_letters:
                return [affine_sum((1, frozen), (1, gap))]
            first_degree, second_degree = (1 if TAIL in group else gcd(*group)
                                           for group in (first, second))
            divisor = gcd(first_degree, second_degree)
            ray_degree = lcm(first_degree, second_degree) if mixed_powers else divisor
            if valuation_rays:
                ray_degree = min(first & second, default=None)
            root_sum = affine_sum((Fraction(1, divisor), upper_weight(letter, first)),
                                  (Fraction(1, divisor), upper_weight(letter, second)))
            test = affine_sum((1, lock_weight), (-1, root_sum))
            narrow = (frozen if ray_degree is None else
                      affine_sum((1, frozen), (Fraction(1, ray_degree), gap)))
            if use_lattice:
                area = (affine_sum((1, root_sum), (1, gap), (-1, one))
                        if archimedean_strip else affine_sum((-1, test)))
                return [narrow, affine_sum((1, frozen), (1, area))]
            weak.append([row(test)])
            strong.append(row(affine_sum((-1, test))))
            return [narrow]

        energy_a = off_diagonal(a, left)
        energy_b = off_diagonal(b, right)
        size_b = variable(b)
        alternatives = (
            size_b,
            *(affine_sum((Fraction(1, 2), energy), (Fraction(1, 2), size_b))
              for energy in energy_a),
            *(affine_sum((Fraction(1, 2), first), (Fraction(1, 2), second))
              for first, second in product(energy_a, energy_b)),
        )
        base = affine_sum((1, objective), (-1, mass(c, locked)))
        return weak + [branch + strong for branch in dominant_rows(base, alternatives)]

    if kind in ("mordell-236-height", "mordell-236-local"):
        local = kind == "mordell-236-local"
        fields = {"type", "S", "T", "window"}
        require(signed and set(node) == fields | (set() if local else {"form"}),
                "unexpected Mordell height metadata")
        square, cube = node["S"], node["T"]
        require(square in LETTERS and cube in LETTERS and square != cube,
                "invalid square and cube letters")
        sixth = next(letter for letter in LETTERS if letter not in (square, cube))
        square_even, sixth_even = (mass(letter, {2, 4, 6}) for letter in (square, sixth))
        residual = affine_sum((1, total), (-1, square_even), (-1, mass(cube, {3, 6})),
                              (-1, mass(sixth, {6})))
        kernel = affine_sum((1, variable(square)), (1, variable(sixth)),
                            (-1, square_even), (-1, sixth_even))
        if local:
            window = node["window"]
            require(type(window) is int and 0 <= window < len(MORDELL_LOCAL_WINDOWS),
                    "invalid local Mordell window")
            threshold, moment = MORDELL_LOCAL_WINDOWS[window]
            require(3**moment.numerator * (3 * threshold)**moment.denominator >= 1,
                    "insufficient local Mordell moment")
            residual_weight = affine_sum(*[
                (1, weighted_early(letter, {i for i in EARLY if i % degree}, degree))
                for letter, degree in ((square, 2), (cube, 3), (sixth, 6))
            ])
            margin = affine_sum((1, residual_weight), (-threshold, one))
            estimate = affine_sum((1, objective), (-1, residual),
                                  (-mordell_moment_loss(moment), kernel))
            return [[row(margin)], [row(affine_sum((-1, margin))), row(estimate)]]

        tau = (
            mass(sixth, {1}),
            affine_sum((1, mass(cube, {1, 4})), (1, mass(sixth, {2}))),
            affine_sum((1, mass(square, {1, 3, 5})), (1, mass(sixth, {3}))),
            affine_sum((1, mass(cube, {2, 5})), (1, mass(sixth, {4}))),
            mass(sixth, {5}),
        )
        height_forms = (
            (Fraction(1, 3),) * 5,
            (1, 0, 0, 0, 5),
            (Fraction(1, 3), Fraction(2, 3), 0, Fraction(4, 3), Fraction(5, 3)),
            (Fraction(1, 4), 0, Fraction(3, 4), 0, Fraction(5, 4)),
            (Fraction(1, 5), Fraction(2, 5), Fraction(2, 5), Fraction(4, 5), 1),
            tuple(Fraction(exponent, 6) for exponent in range(1, 6)),
        )
        form, window = node["form"], node["window"]
        require(type(form) is int and 0 <= form < len(height_forms)
                and type(window) is int and 0 <= window < len(MORDELL_HEIGHT_WINDOWS),
                "invalid Mordell height form or window")
        m, isogeny_degree, u, k = MORDELL_HEIGHT_WINDOWS[window]
        rank_moment = k - 1 if isogeny_degree == 3 else k
        require(isogeny_degree * m * m * u > 2
                and 3**rank_moment.numerator >= m**(2 * rank_moment.denominator),
                "insufficient Mordell height window or moment")
        moment_loss = mordell_moment_loss(k)
        height = affine_sum(*zip(height_forms[form], tau, strict=True))
        margin = affine_sum((1, height), (-u, one))
        estimate = affine_sum((1, objective), (-1, residual), (-moment_loss, kernel))
        return [[row(margin)], [row(affine_sum((-1, margin))), row(estimate)]]

    if kind == "holder-fixed-moments":
        require(signed and set(node) == {"type", "C", "O", "H"},
                "unexpected fixed Hölder metadata")
        c, moments = node["C"], node["H"]
        require(c in LETTERS and isinstance(moments, list) and len(moments) == 2,
                "invalid fixed Hölder letter or moments")
        lock = subset(node["O"])
        weight = lower_weight(c, lock)
        cap = affine_sum((1, mass(c, lock)))
        guards = []
        free_letters = (letter for letter in LETTERS if letter != c)
        for letter, assigned in zip(free_letters, moments, strict=True):
            if assigned is None:
                energy = variable(letter)
                guards.append(affine_sum((2, energy), (-1, weight)))
            else:
                require(isinstance(assigned, list) and all(
                    isinstance(pair, list) and len(pair) == 2
                    and type(pair[0]) is int and pair[0] in EARLY
                    and type(pair[1]) is int and 1 <= pair[1] <= 128
                    for pair in assigned
                ), "invalid fixed Hölder coordinate")
                indices = [index for index, _ in assigned]
                require(indices == sorted(set(indices))
                        and sum(Fraction(1, order) for _, order in assigned) <= 1,
                        "invalid fixed Hölder indices or reciprocal sum")
                for index, order in assigned:
                    guards.append(affine_sum((order, variable(f"{letter}{index}")), (-1, weight)))
                energy = affine_sum((2, variable(letter)), (-1, mass(letter, set(indices))))
            cap = affine_sum((1, cap), (Fraction(1, 2), energy))
        strong = [row(guard) for guard in guards]
        strong.append(row(affine_sum((1, objective), (-1, cap))))
        return [[row(affine_sum((-1, guard)))] for guard in guards] + [strong]

    if kind == "holder3":
        require(signed and set(node) == {"type", "C", "O"}, "unexpected Hölder metadata")
        c = node["C"]
        require(c in LETTERS, "invalid locked letter")
        a, b = [letter for letter in LETTERS if letter != c]
        locked = subset(node["O"])
        modulus = lower_weight(c, locked)

        def cases(letter):
            free = variable(letter)
            below = affine_sum((1, modulus), (-1, free))
            above = affine_sum((1, modulus), (-2, free))
            return (
                (affine_sum((2, free), (-Fraction(33, 40), modulus)), [row(below)]),
                (affine_sum((Fraction(27, 20), free), (-Fraction(7, 40), modulus)),
                 [row(affine_sum((-1, below))), row(above)]),
                (free, [row(affine_sum((-1, above)))]),
            )

        base = affine_sum((1, objective), (-1, mass(c, locked)))
        return [[row(affine_sum((1, base), (-Fraction(1, 2), left), (-Fraction(1, 2), right)))]
                + left_region + right_region
                for (left, left_region), (right, right_region) in product(cases(a), cases(b))]

    if kind in ("rlock", "rlock-strip", "rlock-slice", "rlock-strip-slice"):
        strip = kind in ("rlock-strip", "rlock-strip-slice")
        sliced = kind in ("rlock-slice", "rlock-strip-slice")
        fields = {"type", "C", "I", "J", "O", "g"}
        require(signed and set(node) == fields | ({"lead"} if sliced else set()),
                "unexpected rational lock metadata")
        c, degree = node["C"], node["g"]
        require(c in LETTERS and type(degree) is int and 1 <= degree <= 6,
                "invalid rational lock letter or degree")
        require(not strip or c == short, "the rational lock strip requires a short lock")
        a, b = [letter for letter in LETTERS if letter != c]
        left, right, locked = (subset(node[field]) for field in ("I", "J", "O"))
        require(degree == 1 or all(i in EARLY and i % degree == 0 for i in left | right),
                "invalid common power for rational lock")

        def height(letter, group, power):
            if power == 1:
                return upper_weight(letter, group)
            return weighted_early(letter, group, power)

        base = affine_sum((1, objective), (-1, variable(a)), (-1, variable(b)),
                          (1, mass(a, left)), (1, mass(b, right)), (-1, mass(c, locked)))
        clearance = (affine_sum((1, one), (-1, upper_weight(c, ALL - locked)))
                     if strip else lower_weight(c, locked))
        excess = affine_sum((1, height(a, left, degree)), (1, height(b, right, degree)),
                            (-1, clearance))
        if not sliced:
            return [[row(base), row(excess)],
                    [row(affine_sum((1, base), (-1, excess))), row(affine_sum((-1, excess)))]]

        lead = node["lead"]
        require(lead in (a, b), "invalid leading power for rational lock slices")
        p, q = (1 if TAIL in group else gcd(*group) for group in (left, right))
        root_a, root_b = height(a, left, p), height(b, right, q)
        curve = (affine_sum((Fraction(degree, p), root_b)) if lead == a
                 else affine_sum((Fraction(degree, q), root_a)))
        half_area = affine_sum((Fraction(1, 2), excess))
        return dominant_rows(base, (
            zero,
            affine_sum((1, half_area), (1, curve)),
            affine_sum((1, half_area), (1, root_a), (-1, clearance)),
            affine_sum((1, half_area), (1, root_b), (-1, clearance)),
        ))

    if kind == "tailgeo":
        require(set(node) == {"type", "I", "J", "K"}, "unexpected geometry metadata")
        released = [(letter, subset(node[field], True))
                    for letter, field in zip(LETTERS, ("I", "J", "K"), strict=True)]
        free_mass = affine_sum(*[(1, mass(letter, group)) for letter, group in released])
        excess_weight = affine_sum(*[(1, upper_weight(letter, group))
                                     for letter, group in released], (-1, free_mass))
        first = affine_sum((1, objective), (-1, total), (1, free_mass))
        second = affine_sum((1, objective), (-1, total), (1, one), (-1, excess_weight))
        return [[row(first)], [row(second)] + ([row(affine_sum((-1, first)))] if signed else [])]

    if kind == "residue":
        require(set(node) == {"type", "C", "O", "G", "H", "J", "K"},
                "unexpected residue-energy metadata")
        c = node["C"]
        require(c in LETTERS, "invalid locked letter")
        a, b = [letter for letter in LETTERS if letter != c]
        fixed_c = subset(node["O"])
        modulus = lower_weight(c, fixed_c)
        energy = affine_sum((2, variable(a)), (2, variable(b)))
        penalties = []
        for letter, first, second in ((a, "G", "H"), (b, "J", "K")):
            groups = (subset(node[first]), subset(node[second]))
            require(not groups[0].intersection(groups[1]), "overlapping residue groups")
            for group in groups:
                energy = affine_sum((1, energy), (-1, mass(letter, group)))
                penalties.append((Fraction(1, 4), affine_sum(
                    (2, root_height(letter, group)), (-1, modulus)
                )))
        base = affine_sum((1, objective), (-1, mass(c, fixed_c)), (-Fraction(1, 2), energy))
        return max_rows(base, penalties)

    require(kind == "tail" and set(node) == {"type", "A", "B", "G", "H", "I", "O", "cap"},
            "unexpected collision metadata")
    a, b = node["A"], node["B"]
    require(a in LETTERS and b in LETTERS and a != b, "invalid collision letters")
    c = next(letter for letter in LETTERS if letter not in (a, b))
    first, second, free_b, fixed_c = (
        subset(node[field], field == "O") for field in ("G", "H", "I", "O")
    )
    require(not first.intersection(second), "overlapping collision groups")
    free_mass_b = mass(b, free_b)
    frozen_b = affine_sum((1, variable(b)), (-1, free_mass_b))
    frozen_c = mass(c, fixed_c)
    modulus_b = lower_weight(b, ALL - free_b)
    modulus_c = lower_weight(c, fixed_c)
    energy = affine_sum((2, variable(a)), (-1, mass(a, first)), (-1, mass(a, second)))
    penalties = [(Fraction(1, 4), affine_sum((2, root_height(a, group)), (-1, modulus_b)))
                 for group in (first, second)]

    cap = node["cap"]
    if cap == ["target"]:
        difference = affine_sum((1, variable(c)), (-1, frozen_c))
    elif cap == ["baseline"]:
        difference = free_mass_b
    elif (type(cap) is list and len(cap) == 2 and cap[0] == "diff"
          and type(cap[1]) is int and 2 <= cap[1] <= 6):
        difference = affine_sum((2, mass(b, {i for i in free_b if i == TAIL or i % cap[1]})))
    elif cap == ["value"]:
        difference = zero
        penalties.append((Fraction(1, 2), affine_sum((1, upper_weight(b, free_b)), (-1, modulus_c))))
    elif type(cap) is list and len(cap) == 2 and cap[0] == "group":
        group = subset(cap[1])
        require(group.issubset(free_b), "residue group outside released coordinates")
        difference = mass(b, free_b - group)
        penalties.append((Fraction(1, 2), affine_sum((1, root_height(b, group)), (-1, modulus_c))))
    else:
        raise ValueError("invalid collision cap")

    diagonal = affine_sum((1, objective), (-1, frozen_b), (-1, variable(c)))
    off_diagonal = affine_sum(
        (1, objective), (-1, frozen_b), (-Fraction(1, 2), variable(c)),
        (-Fraction(1, 2), frozen_c), (-Fraction(1, 2), energy), (-Fraction(1, 2), difference)
    )
    off_branches = max_rows(off_diagonal, penalties)
    if signed:
        for branch in off_branches:
            branch.append(row(affine_sum((-1, diagonal))))
    return [[row(diagonal)]] + off_branches


def sparse_rows(rows):
    """Use one common denominator for each row and retain its nonzero coefficients."""
    result = []
    for row, bound in rows:
        denominator = lcm(bound.denominator, *(value.denominator for value in row))
        rhs = bound.numerator * (denominator // bound.denominator)
        coefficients = tuple((i, v.numerator * (denominator // v.denominator))
                             for i, v in enumerate(row) if v)
        result.append((denominator, rhs, coefficients))
    return tuple(result)


def check_leaf(checker, node, rows, limit):
    require(set(node) == {"type", "y", "bound"} and isinstance(node["y"], dict)
            and node["y"], "invalid rational leaf")
    terms = []
    for row_text, multiplier_text in node["y"].items():
        index, multiplier = int(row_text), Fraction(multiplier_text)
        require(str(index) == row_text and 0 <= index < len(rows) and multiplier > 0,
                "invalid leaf row or multiplier")
        denominator, rhs, coefficients = rows[index]
        terms.append((multiplier.numerator, multiplier.denominator * denominator, rhs, coefficients))

    common = lcm(*(denominator for _, denominator, _, _ in terms))
    vector = [0] * checker.N
    right = 0
    for numerator, denominator, rhs, coefficients in terms:
        scale = numerator * (common // denominator)
        right += scale * rhs
        for i, value in coefficients:
            vector[i] += scale * value
    require(vector[checker.D] == common and min(vector) >= 0, "invalid leaf dual coefficients")
    bound = Fraction(right, common)
    require(bound == Fraction(node["bound"]) <= limit, "invalid leaf dual bound")
    return bound


def certificate_root_rows(checker, budget, schema):
    rows = initial_rows(checker, budget)
    if schema >= 3:
        rows += tail_fourier_rows(checker)
    if schema >= 4:
        rows += residual_signature_rows(checker)
    if schema >= 5:
        rows += residual_signature_rows(checker, euclidean=True)
    if schema >= 6:
        rows += square_cube_rows(checker)
    if schema >= 7:
        rows += square_fifth_rows(checker)
    return rows


def verify(certificate, at_budget=None):
    require(set(certificate) == {"schema_version", "budget", "threshold", "tree"}
            and certificate["schema_version"] in (1, 2, 3, 4, 5, 6, 7), "unknown tall-term certificate schema")
    budget, threshold = Fraction(certificate["budget"]), Fraction(certificate["threshold"])
    require(0 < budget <= 1 and threshold > 0, "invalid certificate budget or threshold")
    budget_change = Fraction(0) if at_budget is None else at_budget - budget
    checker = load_checker()
    counts = Counter()
    signed = certificate["schema_version"] >= 2

    @lru_cache(maxsize=8192)
    def expand(metadata):
        return tuple(sparse_rows(rows) for rows in branch_rows(checker, json.loads(metadata), signed))

    def walk(node, rows):
        kind = node["type"]
        counts[kind] += 1
        if kind == "leaf":
            bound = check_leaf(checker, node, rows, threshold)
            if budget_change:
                bound += Fraction(node["y"].get("0", 0)) * budget_change
            return bound
        metadata = json.dumps({key: value for key, value in node.items() if key != "children"},
                              sort_keys=True, separators=(",", ":"))
        branches = expand(metadata)
        children = node.get("children")
        require(isinstance(children, list) and len(children) == len(branches),
                "missing affine branch")
        return max(walk(child, rows + added) for child, added in zip(children, branches, strict=True))

    upper = walk(certificate["tree"], sparse_rows(certificate_root_rows(
        checker, budget, certificate["schema_version"]
    )))
    result = {"budget": str(budget), "threshold": str(threshold), "upper": str(upper), **counts}
    if at_budget is not None:
        result["evaluated_budget"] = str(at_budget)
    return result


def retarget_certificate(archive, archive_hash, delta):
    """Reconstruct a pinned tree at a higher budget before the exact replay."""
    fields = {"schema_version", "baseline_archive_sha256", "baseline_tree_sha256",
              "source_budget", "source_bound", "target_budget", "target", "replacements"}
    require(set(delta) == fields and delta["schema_version"] == 1, "invalid budget delta header")
    old, new = Fraction(delta["source_budget"]), Fraction(delta["target_budget"])
    old_bound, limit = Fraction(delta["source_bound"]), Fraction(delta["target"])
    require(0 < old < new <= 1 and 0 < old_bound and 0 < limit,
            "invalid budget delta bounds")
    require(archive_hash == delta["baseline_archive_sha256"], "budget delta archive hash mismatch")
    matches = [certificate for certificate in archive
               if Fraction(certificate["budget"]) == old and Fraction(certificate["threshold"]) == old_bound]
    require(len(matches) == 1, "budget delta source tree is missing or ambiguous")
    source = matches[0]
    canonical = json.dumps(source, sort_keys=True, separators=(",", ":")).encode()
    require(sha256(canonical).hexdigest() == delta["baseline_tree_sha256"], "budget delta tree hash mismatch")

    checker = load_checker()
    schema = source["schema_version"]
    old_rows = certificate_root_rows(checker, old, schema)
    new_rows = certificate_root_rows(checker, new, schema)
    require(old_rows[1:] == new_rows[1:] and old_rows[0][0] == new_rows[0][0]
            and old_rows[0][1] == old and new_rows[0][1] == new,
            "the budget must affect only the first row's right-hand side")

    require(isinstance(delta["replacements"], list), "invalid budget delta replacements")
    changes = {}
    for entry in delta["replacements"]:
        require(set(entry) == {"path", "tree"} and isinstance(entry["path"], list)
                and all(type(step) is int and step >= 0 for step in entry["path"]),
                "invalid budget delta replacement path")
        path = tuple(entry["path"])
        require(path not in changes, "duplicate budget delta path")
        changes[path] = entry["tree"]

    def retarget(node, path=()):
        if node["type"] == "leaf":
            bound = Fraction(node["bound"]) + (new - old) * Fraction(node["y"].get("0", 0))
            if path in changes:
                require(bound > limit, "budget delta path does not target an over-limit leaf")
                return changes.pop(path)
            require(bound <= limit, "budget delta leaves an over-limit leaf unchanged")
            return {**node, "bound": str(bound)}
        return {**node, "children": [retarget(child, path + (i,))
                                    for i, child in enumerate(node["children"])]}

    tree = retarget(source["tree"])
    require(not changes, "budget delta does not name existing leaves")
    return {**source, "budget": str(new), "threshold": str(limit), "tree": tree}


def certificates(source_budget=None):
    compressed = CERTIFICATE_FILE.read_bytes()
    archive = json.loads(gzip.decompress(compressed))
    require(isinstance(archive, list) and archive, "empty certificate archive")
    for certificate in archive:
        if source_budget is None or Fraction(certificate["budget"]) == source_budget:
            yield certificate
    archive_hash = sha256(compressed).hexdigest()
    for budget, filename in BUDGET_RETARGET_FILES:
        if source_budget is None or source_budget == budget:
            delta = json.loads(lzma.decompress(CERTIFICATE_FILE.with_name(filename).read_bytes()))
            require(Fraction(delta["target_budget"]) == budget, "unexpected budget delta target")
            yield retarget_certificate(archive, archive_hash, delta)
    if source_budget in (None, Fraction(1)):
        for fragments in (UNIT_CERTIFICATE_PARTS, IMPROVED_UNIT_CERTIFICATE_PARTS):
            yield json.loads(b"".join(
                lzma.decompress(CERTIFICATE_FILE.with_name(name).read_bytes())
                for name in fragments
            ))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--at-budget", type=Fraction, help="reweight the selected trees at this radical budget")
    parser.add_argument("--source-budget", type=Fraction, help="replay only trees built at this radical budget")
    args = parser.parse_args()
    for option, value in (("--at-budget", args.at_budget), ("--source-budget", args.source_budget)):
        if value is not None and not 0 < value <= 1:
            parser.error(f"{option} must lie in (0, 1]")
    bounds = []
    for certificate in certificates(args.source_budget):
        result = verify(certificate, args.at_budget)
        print(json.dumps(result, sort_keys=True), flush=True)
        bounds.append(Fraction(result["upper"]))
        del certificate
    if not bounds:
        parser.error("no certificate has the requested source budget")
    if args.at_budget is not None:
        print(json.dumps({"budget": str(args.at_budget), "upper": str(min(bounds))}))


if __name__ == "__main__":
    main()
