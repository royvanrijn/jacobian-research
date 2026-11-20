#!/usr/bin/env python3
"""Exact rational tests of the constant-vector subtree-shuffling conjecture."""

from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from functools import cache
from hashlib import sha256
from itertools import combinations_with_replacement, permutations, product
from math import comb, factorial, prod
import json
from pathlib import Path


@cache
def partitions(total, length, least=0):
    """Nondecreasing tuples of nonnegative integers with the stated sum."""
    if length == 0:
        return ((),) if total == 0 else ()
    return tuple((value, *tail)
                 for value in range(least, total // length + 1)
                 for tail in partitions(total - value, length - 1, value))


def multinomial(items):
    return factorial(len(items)) // prod(factorial(n) for n in Counter(items).values())


def add_scaled(target, row, scale):
    for index, coefficient in row.items():
        value = target.get(index, 0) + coefficient * scale
        if value:
            target[index] = value
        else:
            target.pop(index, None)


def normal_form(coefficients, pivots):
    pending = dict(coefficients)
    remainder = {}
    while pending:
        index = max(pending)
        coefficient = pending[index]
        if index in pivots:
            add_scaled(pending, pivots[index], -coefficient)
        else:
            remainder[index] = pending.pop(index)
    return remainder


def insert_relation(row, pivots):
    reduced = normal_form(row, pivots)
    if reduced:
        pivot = max(reduced)
        leading = Fraction(reduced[pivot])
        pivots[pivot] = {j: Fraction(c) / leading for j, c in reduced.items()}


def dual_coordinate(free_index, pivots):
    """Recover a functional killing the relations and reading one free coefficient."""
    values = {free_index: Fraction(1)}
    for pivot, row in sorted(pivots.items()):
        value = -sum(c * values.get(j, 0) for j, c in row.items() if j != pivot)
        if value:
            values[pivot] = value
    return values


class ShuffleQuotient:
    """Graded one-generator symmetric d-ary algebra modulo the shuffle relations."""

    def __init__(self, arity, path_length):
        self.arity = arity
        self.path_length = path_length
        self.grade = [0]
        self.tree_ids = {(): 0}
        self.basis = [(0,)]
        self.pivots = [{}]
        self.inverse = [{0: 1}]
        self.unordered_counts = [1]

    def forests(self, total, length):
        """Multisets of quotient basis elements of the stated total grade."""
        for sizes in partitions(total, length):
            groups = (combinations_with_replacement(self.basis[size], amount)
                      for size, amount in Counter(sizes).items())
            for blocks in product(*groups):
                yield sum(blocks, ())

    def graft(self, children):
        return self.tree_ids[tuple(sorted(children))]

    def root_shuffles(self, side, bottom):
        """Sum over distinct assignments of off-path branches to ordered slots."""
        result = {}
        width = self.arity - 1
        counts = Counter(side)
        for group in combinations_with_replacement(sorted(counts), width):
            remaining = counts.copy()
            remaining.subtract(group)
            if any(n < 0 for n in remaining.values()):
                continue
            tail = tuple(remaining.elements())
            scale = multinomial(group)
            for branch, coefficient in self.lower_shuffle(tail, bottom).items():
                tree = self.graft((branch, *group))
                result[tree] = result.get(tree, 0) + scale * coefficient
                if not result[tree]:
                    del result[tree]
        return result

    @cache
    def lower_shuffle(self, side, bottom):
        if not side:
            return {bottom: 1}
        grade = self.grade[bottom] + sum(self.grade[i] for i in side) + len(side) // (self.arity - 1)
        return normal_form(self.root_shuffles(side, bottom), self.pivots[grade])

    def step(self, grade):
        assert grade == len(self.basis)
        arity, path = self.arity, self.path_length
        children_by_id = {}
        for children in self.forests(grade - 1, arity):
            identifier = len(self.grade)
            self.grade.append(grade)
            self.tree_ids[children] = identifier
            children_by_id[identifier] = children

        pivots = {}
        for bottom_grade in range(grade - path + 1):
            for side in self.forests(grade - path - bottom_grade, (arity - 1) * path):
                for bottom in self.basis[bottom_grade]:
                    insert_relation(self.root_shuffles(side, bottom), pivots)
        self.pivots.append(pivots)
        basis = tuple(i for i in children_by_id if i not in pivots)
        self.basis.append(basis)

        # Coefficient of t^grade in X+t*[G,...,G], with G the formal inverse.
        coefficient = {}
        for tree, children in children_by_id.items():
            n = multinomial(children) * prod(self.inverse[self.grade[c]].get(c, 0) for c in children)
            if n:
                coefficient[tree] = n
        self.inverse.append(normal_form(coefficient, pivots))

        count = sum(prod(comb(self.unordered_counts[size] + amount - 1, amount)
                         for size, amount in Counter(sizes).items())
                    for sizes in partitions(grade - 1, arity))
        self.unordered_counts.append(count)
        row = {"k": grade, "unordered": count, "quotient_dimension": len(basis),
               "reduced_coefficient_support": len(self.inverse[grade])}
        if self.inverse[grade]:
            free = min(self.inverse[grade])
            values = dual_coordinate(free, pivots)
            assert all(sum(c * values.get(j, 0) for j, c in p.items()) == 0 for p in pivots.values())
            evaluation = sum(c * values.get(j, 0) for j, c in coefficient.items())
            assert evaluation == self.inverse[grade][free]
            row["dual_evaluation"] = str(evaluation)
        return row


@cache
def planar_trees(arity, grade):
    if grade == 0:
        return ((),)
    return tuple(children
                 for sizes in partitions(grade - 1, arity)
                 for ordered in sorted(set(permutations(sizes)))
                 for children in product(*(planar_trees(arity, n) for n in ordered)))


def addresses(tree):
    yield ()
    for i, child in enumerate(tree):
        for address in addresses(child):
            yield (i, *address)


def subtree(tree, address):
    for index in address:
        tree = tree[index]
    return tree


@cache
def distinct_permutations(items):
    return set(permutations(items))


def replacements(tree, lookup, address=()):
    if address in lookup:
        return lookup[address]
    return tuple(replacements(child, lookup, (*address, j)) for j, child in enumerate(tree))


def unordered(tree):
    return tuple(sorted(unordered(child) for child in tree))


def literal_planar_projection(arity, path_length, grade):
    """Build the published shuffle classes by literal path and slot enumeration."""
    trees = planar_trees(arity, grade)
    shapes = sorted({unordered(tree) for tree in trees})
    ids = {shape: i for i, shape in enumerate(shapes)}
    columns = set()
    for tree in trees:
        for end in addresses(tree):
            if len(end) < path_length:
                continue
            slots = tuple((*end[:i], child)
                          for i in range(len(end) - path_length, len(end))
                          for child in range(arity) if child != end[i])
            side = tuple(sorted(subtree(tree, address) for address in slots))
            shuffle_class = frozenset(replacements(tree, dict(zip(slots, perm)))
                                      for perm in distinct_permutations(side))
            columns.add(shuffle_class)
    pivots = {}
    for column in columns:
        insert_relation(Counter(ids[unordered(tree)] for tree in column), pivots)
    target = Counter(ids[unordered(tree)] for tree in trees)
    return {"d": arity, "p": path_length, "k": grade, "planar": len(trees),
            "planar_classes": len(columns), "unordered": len(shapes),
            "quotient_dimension": len(shapes) - len(pivots),
            "constant_in_span": not normal_form(target, pivots)}


def compute(arity, path, limit):
    quotient = ShuffleQuotient(arity, path)
    rows = [quotient.step(k) for k in range(1, limit + 1)]
    return {"d": arity, "p": path, "rows": rows,
            "truncated_algebra_dimension": sum(map(len, quotient.basis))}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", type=Path)
    parser.add_argument("--scan", nargs=3, type=int, metavar=("D", "P", "K"))
    args = parser.parse_args()
    if args.scan:
        d, p, k = args.scan
        if d < 2 or p < 1 or k < 1:
            parser.error("require d >= 2, p >= 1 and k >= 1")
        print(json.dumps(compute(d, p, k), indent=2))
        return

    cases = [(2, 3, 11), (2, 4, 14), (3, 2, 5), (3, 3, 15)]
    results = [compute(*case) for case in cases]
    by_key = {(run["d"], run["p"], row["k"]): row for run in results for row in run["rows"]}
    cross_checks = []
    for d, p, k in [(2, 3, 5), (2, 3, 6), (2, 4, 6), (3, 2, 4), (3, 3, 4), (3, 3, 5)]:
        check = literal_planar_projection(d, p, k)
        row = by_key[d, p, k]
        assert check["unordered"] == row["unordered"]
        assert check["quotient_dimension"] == row["quotient_dimension"]
        assert check["constant_in_span"] == (not row["reduced_coefficient_support"])
        cross_checks.append(check)

    # A block of zero graded pieces propagates by grafting at the root.
    assert all(by_key[2, 3, k]["reduced_coefficient_support"] for k in range(1, 6))
    assert all(by_key[2, 3, k]["quotient_dimension"] == 0 for k in range(6, 12))
    assert all(by_key[3, 2, k]["quotient_dimension"] == 0 for k in range(2, 5))
    assert all(by_key[3, 3, k]["reduced_coefficient_support"] for k in range(1, 15))
    assert all(by_key[2, 4, k]["reduced_coefficient_support"] for k in range(1, 15))
    assert results[1]["truncated_algebra_dimension"] == 680
    assert by_key[3, 3, 14]["quotient_dimension"] == 35
    assert by_key[3, 3, 15]["quotient_dimension"] == 16
    assert not by_key[3, 3, 15]["reduced_coefficient_support"]
    assert results[3]["truncated_algebra_dimension"] == 211
    sources = [Path(__file__), Path(__file__).with_name("SUBTREE_SHUFFLE_AND_INVERSE_DEGREE.md")]
    receipt = {
        "schema_version": 1, "arithmetic": "rational", "runs": results,
        "literal_planar_controls": cross_checks,
        "input_hashes": {path.name: sha256(path.read_bytes()).hexdigest() for path in sources},
    }
    if args.write_receipt:
        args.write_receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print("PASS: exact rational shuffle quotients and independent literal planar controls")


if __name__ == "__main__":
    main()
