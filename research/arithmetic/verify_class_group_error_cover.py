#!/usr/bin/env python3
# Finite tests for the maximal-subgroup cover in the norm class-error argument.
# This checks finite group combinatorics only, not the analytic sieve, class-prime
# theorem, principal-genus identification, or assertions about all class groups.
from argparse import ArgumentParser
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import gcd
from pathlib import Path
import json


class VerificationFailure(Exception):
    pass


def require(condition, message):
    if not condition:
        raise VerificationFailure(message)


class FiniteAbelianGroup:
    def __init__(self, moduli):
        self.moduli = tuple(moduli)
        require(all(n >= 2 for n in self.moduli), "Invalid cyclic modulus")
        self.elements = tuple(product(*(range(n) for n in self.moduli)))
        self.zero = tuple(0 for _ in self.moduli)
        self.size = len(self.elements)

    def add(self, a, b):
        return tuple((x + y) % n for x, y, n in zip(a, b, self.moduli))

    def scale(self, k, a):
        return tuple(k * x % n for x, n in zip(a, self.moduli))

    def order(self, a):
        s = self.zero
        for k in range(1, self.size + 1):
            s = self.add(s, a)
            if s == self.zero:
                return k
        raise VerificationFailure("Element order exceeds group size")

    def generated(self, generators):
        # Independent breadth-first closure, without the orientation formula.
        found = {self.zero}
        queue = [self.zero]
        for a in queue:
            for b in generators:
                c = self.add(a, b)
                if c not in found:
                    found.add(c)
                    queue.append(c)
        return frozenset(found)


def prime_divisors(n):
    answer = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            answer.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        answer.append(n)
    return answer


def doubled_subgroup(group):
    return frozenset(group.scale(2, a) for a in group.elements)


def maximal_subgroups(group):
    # Identify B=2G with the product of cyclic groups of orders n/gcd(n,2).
    # Enumerate every nonzero homomorphism B -> C_ell, then deduplicate kernels.
    bmod = tuple(n // gcd(n, 2) for n in group.moduli)
    coordinates = {}
    for c in product(*(range(n) for n in bmod)):
        b = tuple(2 * x % n for x, n in zip(c, group.moduli))
        require(b not in coordinates, "Noninjective coordinates for B")
        coordinates[b] = c
    B = doubled_subgroup(group)
    require(set(coordinates) == set(B), "Coordinates do not cover B")
    kernels = {}
    for ell in prime_divisors(len(B)):
        ranges = [range(ell) if n % ell == 0 else (0,) for n in bmod]
        for coefficients in product(*ranges):
            if not any(coefficients):
                continue
            M = frozenset(b for b, c in coordinates.items()
                          if sum(a * x for a, x in zip(coefficients, c)) % ell == 0)
            require(len(M) * ell == len(B), "Nonzero homomorphism not surjective")
            kernels[M] = ell
    return sorted(((ell, M) for M, ell in kernels.items()),
                  key=lambda row: (row[0], tuple(sorted(row[1]))))


def enumerate_subgroups(group, B):
    found = {frozenset((group.zero,))}
    queue = list(found)
    for H in queue:
        for b in sorted(B):
            K = group.generated(tuple(sorted(H)) + (b,))
            if K not in found:
                found.add(K)
                queue.append(K)
    return found


def inverse_families(group):
    families = []
    for g in group.elements:
        opposite = group.scale(-1, g)
        if group.scale(2, g) != group.zero and g < opposite:
            families.append((g, opposite))
    return families


def outside_families(group, families, M):
    return frozenset(i for i, (g, _) in enumerate(families)
                     if group.scale(2, g) not in M)


def validate_outside_families(group, families, M, selected):
    require(selected == outside_families(group, families, M),
            "Family selection must use 2g not in M")


def direct_orientation_sums(group, families, selected):
    # Explicitly choose j=0,...,ord(2g)-1 reserved-prime reversals per family.
    chosen = [families[i][0] for i in sorted(selected)]
    ranges = [range(group.order(group.scale(2, g))) for g in chosen]
    result = set()
    for reversals in product(*ranges):
        change = group.zero
        for j, g in zip(reversals, chosen):
            change = group.add(change, group.scale(-2 * j, g))
        result.add(change)
    return frozenset(result)


def check_group(name, moduli):
    group = FiniteAbelianGroup(moduli)
    B = doubled_subgroup(group)
    kernel = frozenset(g for g in group.elements if group.scale(2, g) == group.zero)
    require(len(B) * len(kernel) == group.size, "Doubling fibers have wrong size")
    families = inverse_families(group)
    require(2 * len(families) + len(kernel) == group.size,
            "Inverse families fail to partition the nonambiguous classes")
    maxima = maximal_subgroups(group)
    subgroups = enumerate_subgroups(group, B)
    direct_maxima = {H for H in subgroups if H != B
                     and not any(H < K < B for K in subgroups)}
    require({M for _, M in maxima} == direct_maxima,
            "Homomorphism kernels miss an explicitly enumerated maximal subgroup")
    require(bool(maxima) == (len(B) > 1), "Wrong trivial-B maximal-subgroup case")
    maximal_rows = []
    for ell, M in maxima:
        preimage = frozenset(g for g in group.elements if group.scale(2, g) in M)
        require(len(preimage) * ell == group.size, "Wrong maximal-subgroup preimage")
        require(kernel <= preimage, "Ambiguous class outside doubling preimage")
        outside = outside_families(group, families, M)
        validate_outside_families(group, families, M, outside)
        union = frozenset(g for i in outside for g in families[i])
        require(union == set(group.elements) - preimage, "Inverse-pair union mismatch")
        density = Fraction(len(outside), group.size)
        expected = Fraction(ell - 1, 2 * ell)
        require(density == expected, "Incorrect inverse-pair density")
        cap = sum(group.order(group.scale(2, families[i][0])) - 2 for i in outside)
        require(cap >= 0, "Negative deficiency cap")
        maximal_rows.append({"index_prime": ell, "subgroup": sorted(M),
                             "preimage_size": len(preimage),
                             "outside_inverse_pairs": len(outside),
                             "class_count_density": str(density),
                             "total_exponent_one_deficiency_cap": cap})
    saturated = 0
    covered = 0
    for flags in product((False, True), repeat=len(families)):
        rich = frozenset(i for i, flag in enumerate(flags) if flag)
        H = group.generated([group.scale(2, families[i][0]) for i in rich])
        if H == B:
            # This enumeration is independent of the breadth-first generation.
            require(direct_orientation_sums(group, families, rich) == B,
                    "Generating rich families fail to supply every orientation change")
            saturated += 1
        else:
            containing = [(ell, M) for ell, M in maxima if H <= M]
            require(containing, "Proper rich subgroup has no maximal cover")
            for ell, M in containing:
                require(not (rich & outside_families(group, families, M)),
                        "Outside family was rich in a containing maximal subgroup")
            covered += 1
    delta = min((Fraction(ell - 1, 2 * ell) for ell, _ in maxima), default=None)
    return {"name": name, "cyclic_moduli": list(moduli), "group_order": group.size,
            "double_subgroup_order": len(B), "doubling_kernel_order": len(kernel),
            "subgroups_enumerated_in_B": len(subgroups),
            "inverse_pair_families": len(families), "rich_subsets_tested": 2 ** len(families),
            "orientation_saturated_subsets": saturated,
            "proper_subgroup_covered_subsets": covered,
            "delta_star": None if delta is None else str(delta),
            "maximal_subgroups": maximal_rows}


def expect_rejected(name, operation):
    try:
        operation()
    except VerificationFailure:
        return {"name": name, "result": "rejected"}
    raise VerificationFailure("Negative control was accepted: " + name)


def negative_controls():
    group = FiniteAbelianGroup((8,))
    families = inverse_families(group)
    ell, M = maximal_subgroups(group)[0]
    wrong = frozenset(i for i, (g, _) in enumerate(families) if g not in M)
    controls = [expect_rejected("g_not_in_M_instead_of_2g_not_in_M",
        lambda: validate_outside_families(group, families, M, wrong))]
    expected = Fraction(ell - 1, 2 * ell)
    controls.append(expect_rejected("one_inverse_family_in_place_of_the_union_density",
        lambda: require(Fraction(1, group.size) == expected, "Union density lost")))
    group = FiniteAbelianGroup((3, 3))
    families = inverse_families(group)
    _, fixed_M = maximal_subgroups(group)[0]
    i = next(i for i, (g, _) in enumerate(families) if group.scale(2, g) not in fixed_M)
    H = group.generated([group.scale(2, families[i][0])])
    require(H != doubled_subgroup(group), "Fixed-maximal counterexample is not deficient")
    controls.append(expect_rejected("one_fixed_maximal_subgroup_in_place_of_all_maxima",
        lambda: require(H <= fixed_M, "A different maximal subgroup is needed")))
    return controls


def verify():
    cases = [("trivial", ()), ("C2", (2,)), ("C2xC2", (2, 2)),
             ("C3", (3,)), ("C4", (4,)), ("C6", (6,)), ("C8", (8,)),
             ("C2xC4", (2, 4)), ("C3xC3", (3, 3)), ("C4xC4", (4, 4)),
             ("C12", (12,)), ("C2xC8", (2, 8))]
    rows = [check_group(name, moduli) for name, moduli in cases]
    return {"scope": "Exact finite group enumeration for the listed cyclic products only. Not an analytic proof, an infinite class-group test, a class-prime distribution test, or a conductor-local ideal verification.",
            "arithmetic": "Python 3 standard library; exact integer operations and fractions; explicit checks remain active under optimization.",
            "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "reproduction_command": "python3 research/arithmetic/verify_class_group_error_cover.py --output research/arithmetic/class_group_error_cover_checks.json",
            "maximal_subgroup_method": "All nonzero homomorphisms from B=2G to C_ell for prime ell dividing |B|; kernels deduplicated and compared with independent exhaustive subgroup generation.",
            "orientation_method": "Breadth-first subgroup generation compared with explicit products of cyclic reversal ranges for every generating rich subset.",
            "groups_tested": len(rows), "rich_subsets_tested": sum(r["rich_subsets_tested"] for r in rows),
            "group_cases": rows, "negative_controls": negative_controls(), "result": "passed"}


def main():
    parser = ArgumentParser(description="Check the finite maximal-subgroup class-error cover.")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify()
    encoded = json.dumps(result, indent=2, sort_keys=True) + chr(10)
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.write_text(encoded)
        print(json.dumps({"output": str(args.output), "groups_tested": result["groups_tested"],
                          "rich_subsets_tested": result["rich_subsets_tested"], "result": result["result"]},
                         indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
