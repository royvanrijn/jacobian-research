"""Replay the abc exponent certificate with the two tall terms distinguished."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
from itertools import product
import json
from math import gcd
from pathlib import Path
from types import ModuleType

from verify_abc_pair_deviation_certificate import (
    CHECKER_FILE, CHECKER_SHA256, require, source_bytes,
)


CERTIFICATE_FILE = Path(__file__).with_name("abc_height_fulloc_certificate.json.gz")
LETTERS = ("a", "b", "c")
EARLY = frozenset(range(1, 7))
TAIL = 7
ALL = EARLY | {TAIL}


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

    objective = variable("D")
    kind = node["type"]
    total = affine_sum(*[(1, variable(letter)) for letter in LETTERS])


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

    if kind == "rlock":
        require(signed and set(node) == {"type", "C", "I", "J", "O", "g"},
                "unexpected rational lock metadata")
        c, degree = node["C"], node["g"]
        require(c in LETTERS and type(degree) is int and 1 <= degree <= 6,
                "invalid rational lock letter or degree")
        a, b = [letter for letter in LETTERS if letter != c]
        left, right, locked = (subset(node[field]) for field in ("I", "J", "O"))
        require(degree == 1 or all(i in EARLY and i % degree == 0 for i in left | right),
                "invalid common power for rational lock")

        def height(letter, group):
            if degree == 1:
                return upper_weight(letter, group)
            return weighted_early(letter, group, degree)

        base = affine_sum((1, objective), (-1, variable(a)), (-1, variable(b)),
                          (1, mass(a, left)), (1, mass(b, right)), (-1, mass(c, locked)))
        excess = affine_sum((1, height(a, left)), (1, height(b, right)),
                            (-1, lower_weight(c, locked)))
        return [[row(base), row(excess)],
                [row(affine_sum((1, base), (-1, excess))), row(affine_sum((-1, excess)))]]

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


def check_leaf(checker, node, rows, threshold):
    require(set(node) == {"type", "y", "bound"} and isinstance(node["y"], dict)
            and node["y"], "invalid rational leaf")
    vector = [Fraction(0)] * checker.N
    bound = Fraction(0)
    for name, text in node["y"].items():
        require(name.isdigit(), "invalid row index")
        index = int(name)
        multiplier = Fraction(text)
        require(str(index) == name and index < len(rows) and multiplier > 0,
                "invalid leaf multiplier")
        coefficients, right = rows[index]
        vector = [a + multiplier * b for a, b in zip(vector, coefficients, strict=True)]
        bound += multiplier * right
    require(vector[checker.D] == 1 and min(vector) >= 0, "invalid dual coefficients")
    require(bound == Fraction(node["bound"]) <= threshold, "invalid rational leaf bound")
    return bound


def verify(certificate, at_budget=None):
    require(set(certificate) == {"schema_version", "budget", "threshold", "tree"}
            and certificate["schema_version"] in (1, 2, 3), "unknown tall-term certificate schema")
    budget, threshold = Fraction(certificate["budget"]), Fraction(certificate["threshold"])
    require(0 < budget <= 1 and threshold > 0, "invalid certificate budget or threshold")
    budget_change = Fraction(0) if at_budget is None else at_budget - budget
    checker = load_checker()
    counts = Counter()
    signed = certificate["schema_version"] >= 2

    @lru_cache(maxsize=8192)
    def expand(metadata):
        return branch_rows(checker, json.loads(metadata), signed)

    def walk(node, rows):
        kind = node["type"]
        counts[kind] += 1
        if kind == "leaf":
            bound = check_leaf(checker, node, rows, threshold)
            return bound + Fraction(node["y"].get("0", 0)) * budget_change
        metadata = json.dumps({key: value for key, value in node.items() if key != "children"},
                              sort_keys=True, separators=(",", ":"))
        branches = expand(metadata)
        children = node.get("children")
        require(isinstance(children, list) and len(children) == len(branches),
                "missing affine branch")
        return max(walk(child, rows + added) for child, added in zip(children, branches, strict=True))

    root = initial_rows(checker, budget)
    if certificate["schema_version"] == 3:
        root += tail_fourier_rows(checker)
    upper = walk(certificate["tree"], root)
    result = {"budget": str(budget), "threshold": str(threshold), "upper": str(upper), **counts}
    if at_budget is not None:
        result["evaluated_budget"] = str(at_budget)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--at-budget", type=Fraction, help="reweight all trees at this radical budget")
    args = parser.parse_args()
    if args.at_budget is not None and not 0 < args.at_budget <= 1:
        parser.error("--at-budget must lie in (0, 1]")
    archive = json.loads(gzip.decompress(CERTIFICATE_FILE.read_bytes()))
    require(isinstance(archive, list) and archive, "empty certificate archive")
    bounds = []
    for certificate in archive:
        result = verify(certificate, args.at_budget)
        print(json.dumps(result, sort_keys=True), flush=True)
        bounds.append(Fraction(result["upper"]))
    if args.at_budget is not None:
        print(json.dumps({"budget": str(args.at_budget), "upper": str(min(bounds))}))


if __name__ == "__main__":
    main()
