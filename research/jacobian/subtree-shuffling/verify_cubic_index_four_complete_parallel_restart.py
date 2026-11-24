#!/usr/bin/env python3
"""Exact complete cubic parallel-restart classification and inverse-support bound."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations_with_replacement
import json
from math import comb
from pathlib import Path
import platform
import time

if not __debug__:
    raise SystemExit("Refusing optimized Python: exact assertions are required.")

HERE = Path(__file__).resolve().parent
STEM = "cubic_index_four_complete_parallel_restart"
PROOF = "CUBIC_INDEX_FOUR_COMPLETE_PARALLEL_RESTART.md"
ASSURANCE = {"formal_verification": False, "external_review": False,
             "independent_replay": False, "worldwide_priority": False,
             "unrestricted_index_four_bound": False, "all_restart_rounds": False}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add(*polys):
    out = defaultdict(F)
    for poly in polys:
        for mon, coefficient in poly.items():
            out[mon] += coefficient
    return {mon: coefficient for mon, coefficient in out.items() if coefficient}


def scale(poly, scalar):
    return {mon: scalar * c for mon, c in poly.items() if scalar * c}


def multiply(left, right):
    out = defaultdict(F)
    for mon, c in left.items():
        for other, d in right.items():
            out[tuple(sorted(mon + other))] += c * d
    return {mon: c for mon, c in out.items() if c}


def diff(poly, variable):
    out = {}
    for mon, c in poly.items():
        multiplicity = mon.count(variable)
        if multiplicity:
            reduced = list(mon)
            reduced.remove(variable)
            out[tuple(reduced)] = c * multiplicity
    return out


def apply(poly, rows):
    return add(*(multiply(diff(poly, i), image) for i, image in rows.items()))


def monomial(*variables, coefficient=1):
    return {tuple(sorted(variables)): F(coefficient)}


def pack(poly):
    return [[list(mon), str(c)] for mon, c in sorted(poly.items())]


def column_kernel(columns):
    """Exact elimination with tracked combinations; every dependent column is retained."""
    pivots, null = {}, []
    for j, original in enumerate(columns):
        vector = dict(original)
        combination = {j: F(1)}
        while vector:
            key = max(vector)
            coefficient = vector[key]
            if key not in pivots:
                vector = scale(vector, 1 / coefficient)
                combination = scale(combination, 1 / coefficient)
                pivots[key] = (vector, combination)
                break
            pivot, source = pivots[key]
            vector = add(vector, scale(pivot, -coefficient))
            combination = add(combination, scale(source, -coefficient))
        else:
            assert combination[j] == 1 and all(i <= j for i in combination)
            null.append(combination)
    assert len(pivots) + len(null) == len(columns)
    return len(pivots), null


def kernel_space(number, degree, operators, weights):
    groups = defaultdict(list)
    for mon in combinations_with_replacement(range(number), degree):
        groups[sum(weights[i] for i in mon)].append(mon)
    assert sum(map(len, groups.values())) == comb(number + degree - 1, degree)
    basis, total_rank = [], 0
    for weight, monomials in sorted(groups.items()):
        columns = [{(tag, image_mon): c
                    for tag, rows in enumerate(operators)
                    for image_mon, c in apply({mon: F(1)}, rows).items()}
                   for mon in monomials]
        rank, null = column_kernel(columns)
        total_rank += rank
        for combination in null:
            poly = {monomials[i]: c for i, c in combination.items()}
            assert all(sum(weights[i] for i in mon) == weight for mon in poly)
            assert all(not apply(poly, rows) for rows in operators)
            basis.append(poly)
    assert total_rank + len(basis) == comb(number + degree - 1, degree)
    return total_rank, basis


class Substitution:
    def __init__(self, images):
        self.images = images
        self.term = lru_cache(maxsize=512)(self._term)

    def _term(self, mon):
        if not mon:
            return {(): F(1)}
        return multiply(self.images[mon[0]], self.term(mon[1:]))

    def __call__(self, poly):
        return add(*(scale(self.term(mon), c) for mon, c in poly.items()))


def matmul(left, right):
    rows = defaultdict(list)
    for (i, j), poly in right.items():
        rows[i].append((j, poly))
    out = {}
    for (i, k), poly in left.items():
        for j, other in rows[k]:
            key = (i, j)
            out[key] = add(out.get(key, {}), multiply(poly, other))
    return {key: poly for key, poly in out.items() if poly}


def read_literal(path):
    literal = json.loads(path.read_text())
    assert set(literal) == {"schema_version", "dimension", "weights", "H_terms"}
    assert literal["schema_version"] == 1 and literal["dimension"] == 58
    weights = literal["weights"]
    assert len(weights) == 58 and all(type(w) is int and w > 0 for w in weights)
    components = [{} for _ in range(58)]
    seen = set()
    for row in literal["H_terms"]:
        assert set(row) == {"output", "inputs", "coefficient"}
        i, mon = row["output"], tuple(row["inputs"])
        c = F(row["coefficient"])
        assert type(i) is int and 0 <= i < 58 and len(mon) == 3
        assert all(type(j) is int and 0 <= j < i for j in mon)
        assert mon == tuple(sorted(mon)) and c
        assert (i, mon) not in seen
        seen.add((i, mon))
        assert sum(weights[j] for j in mon) == weights[i]
        components[i][mon] = c
    assert len(seen) == 1248 and not components[0]
    assert weights[0] == 1 and min(weights[1:]) > 1
    assert max(weights[:48]) == weights[47] == 29
    mapping = {0: 0, 47: 1, **{i: i - 46 for i in range(48, 57)}}
    local = {}
    for i in range(48, 58):
        local[i - 46] = {tuple(sorted(mapping[j] for j in mon)): c
                         for mon, c in components[i].items()}
    return components, weights, local


def local_operators(components):
    delta = {i: diff(components[i], 1) for i in range(2, 11)}
    eta = {i: add(*(multiply(diff(components[i], j), diff(components[j], 0))
                    for j in range(2, i))) for i in range(2, 11)}
    fibres = {i: {k: add(*(multiply(diff(components[i], j),
                                   diff(components[j], k)) for j in range(2, i)))
                   for k in range(2, 11)} for i in range(2, 11)}
    epsilon = {}
    gamma = {}
    for i in range(2, 11):
        assert not apply(components[i], delta)
        assert all(not fibres[i][k] for k in range(5, 11))
        assert not add(*(multiply(fibres[i][k], diff(components[k], 1))
                         for k in range(2, 5)))
        epsilon[i] = add(*(multiply(fibres[i][k], diff(components[k], 0))
                           for k in range(2, 5)))
        quotient = {}
        for mon, c in fibres[i][2].items():
            assert mon.count(1) >= 2
            reduced = list(mon)
            reduced.remove(1)
            reduced.remove(1)
            quotient[tuple(reduced)] = c
        gamma[i] = quotient
    for i in range(2, 12):
        assert not apply(components[i], delta)
        assert not apply(components[i], epsilon)
    expected_gamma = {8: monomial(0, 1, coefficient=9),
                      9: monomial(1, 1, coefficient=-36),
                      10: add(monomial(1, 3, coefficient=-54),
                              monomial(0, 4, coefficient=36))}
    assert {i: q for i, q in gamma.items() if q} == expected_gamma
    return delta, epsilon, eta, gamma, fibres


def source_pins(repo, ledger_path):
    ledger = json.loads(ledger_path.read_text())
    assert set(ledger) == {"schema", "sources", "assurance"}
    assert ledger["schema"] == 1 and ledger["assurance"] == ASSURANCE
    roles = {row["role"] for row in ledger["sources"]}
    assert roles == {"proof", "checker", "literal_data", "literal_proof",
                     "pre_S_row_proof", "third_row_proof", "primary_context",
                     "cube_zero_inverse_corollary"}
    paths = {}
    for row in ledger["sources"]:
        assert set(row) == {"role", "path", "sha256", "scope"}
        path = (repo / row["path"]).resolve()
        assert path.is_relative_to(repo) and path.is_file()
        assert sha(path) == row["sha256"], "Source pin mismatch: " + str(path)
        paths[row["role"]] = path
    assert paths["checker"] == Path(__file__).resolve()
    assert paths["proof"] == HERE / PROOF
    return ledger, paths


def run(paths):
    components, literal_weights, local = read_literal(paths["literal_data"])
    delta, epsilon, eta, gamma, fibres = local_operators(local)
    old_weights = [literal_weights[i] for i in [0, 47, *range(48, 57)]]
    old_rank, helpers = kernel_space(11, 3, [delta, epsilon], old_weights)
    assert old_rank == 265 and len(helpers) == 21
    helper_weights = [sum(old_weights[i] for i in next(iter(q))) for q in helpers]
    positions = [i for i, weight in enumerate(helper_weights) if weight == 207]
    assert len(positions) == 1
    s_position = positions[0]
    h = helpers[s_position]
    first = next(iter(h))
    ratio = local[11][first] / h[first]
    assert scale(h, ratio) == local[11]
    helpers[s_position] = local[11]  # Canonical helper is S itself, not a scalar copy.
    s_index = 11 + s_position
    assert s_index == 31
    operators = [dict(delta), dict(epsilon), {}]
    omega = {2: monomial(1, 1), 3: monomial(0, 1, coefficient=-2),
             4: monomial(0, 0)}
    for i, helper in enumerate(helpers, 11):
        operators[0][i] = diff(helper, 1)
        operators[1][i] = apply(helper, eta)
        operators[2][i] = apply(helper, gamma)
        for k in range(2, 5):
            full_fibre = add(*(multiply(diff(helper, j), fibres[j][k])
                               for j in range(2, 11)))
            assert full_fibre == multiply(operators[2][i], omega[k])
    weights = old_weights + helper_weights
    for rows, shift in zip(operators, [29, 1, 89]):
        for i, poly in rows.items():
            assert all(sum(weights[j] for j in mon) == weights[i] - shift
                       for mon in poly)
    # These weight checks justify independent block elimination, not a degree bound.
    images = [monomial(0), monomial(1)]
    full_components = {**{i: local[i] for i in range(2, 11)},
                       **{i: q for i, q in enumerate(helpers, 11)}}
    for i in range(2, 32):
        images.append(add(monomial(i), Substitution(images)(full_components[i])))
    inverse = Substitution(images)
    forward_images = [add(monomial(i), scale(full_components.get(i, {}), -1))
                      for i in range(32)]
    forward = Substitution(forward_images)
    for i in range(32):
        assert inverse(forward_images[i]) == monomial(i)
        assert forward(images[i]) == monomial(i)
    expected = {1: (20, 12), 2: (439, 89), 3: (5462, 522)}
    profiles = {}
    rank_records = []
    all_bases = {}
    for degree in [1, 2, 3]:
        rank, basis = kernel_space(32, degree, operators, weights)
        assert (rank, len(basis)) == expected[degree]
        profile = {}
        for q in basis:
            assert not diff(q, s_index)
            substituted = inverse(q)
            assert substituted
            for mon in substituted:
                w_degree = mon.count(1)
                rest_degree = len(mon) - w_degree
                profile[w_degree] = max(profile.get(w_degree, -1), rest_degree)
        profiles[degree] = profile
        all_bases[degree] = basis
        rank_records.append({"degree": degree,
                             "domain_dimension": comb(31 + degree, degree),
                             "rank": rank, "kernel_dimension": len(basis),
                             "all_basis_vectors_ignore_S": True,
                             "inverse_W_degree": max(profile),
                             "inverse_total_degree": max(b + a for b, a in profile.items()),
                             "inverse_profile": [[b, a] for b, a in sorted(profile.items())],
                             "basis_sha256": hashlib.sha256(json.dumps(
                                 [pack(q) for q in basis], separators=(",", ":")).encode()).hexdigest()})
    assert profiles == {1: {0: 3, 1: 0},
                        2: {0: 6, 1: 5, 2: 4, 3: 3, 4: 4, 5: 3},
                        3: {0: 9, 1: 8, 2: 7, 3: 6, 4: 7, 5: 6, 6: 5, 7: 4}}
    # First-layer helpers obey the same final envelope as cubic restarts.
    assert all(mon.count(1) <= 7 and len(mon) <= 11 for image in images for mon in image)
    sharp = inverse(local[11])
    assert sharp[(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1)] == -2
    assert {mon: c for mon, c in sharp.items() if all(i < 2 for i in mon)} == \
        monomial(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, coefficient=-2)
    assert all(not apply(local[11], rows) for rows in operators)
    # Directly validate the retained sharp base, including covector independence.
    base_j = {(i, j): diff(components[i], j)
              for i in range(48) for j in range(48) if diff(components[i], j)}
    base_square = matmul(base_j, base_j)
    assert not matmul(base_square, base_j)
    r_column = base_square[(47, 1)]
    independence_coefficient = F(393016350032033013, 3404635167332)
    assert r_column[(0, 1, 1, 30)] == independence_coefficient
    coefficients = []
    for i, component in enumerate(components):
        value = F(i == 0)
        for mon, c in component.items():
            product = c
            for j in mon:
                product *= coefficients[j]
            value += product
        coefficients.append(value)
    c = F(8803690708443294969, 95329784685296)
    assert coefficients[47] == c and coefficients[57] == -2 * c**7
    # Negative controls reject an invalid helper and a direct attempted S reuse.
    corrupted_helper = add(local[11], monomial(2, 2, 2))
    assert apply(corrupted_helper, delta)
    invalid_restart = monomial(1, 1, s_index)
    assert apply(invalid_restart, operators[2])
    sharp_weighted_degree = max(len(mon) + 28 * mon.count(1) for mon in sharp)
    assert sharp_weighted_degree == 207 and sharp_weighted_degree > 206
    # A duplicate S permits its frozen difference, preventing an overstrong claim.
    duplicate_index = 32
    duplicated = [dict(rows) for rows in operators]
    for rows in duplicated:
        rows[duplicate_index] = rows[s_index]
    difference = add(monomial(0, 0, s_index),
                     monomial(0, 0, duplicate_index, coefficient=-1))
    assert diff(difference, s_index)
    assert all(not apply(difference, rows) for rows in duplicated)
    duplicate_inverse = Substitution(images + [add(monomial(duplicate_index), sharp)])
    assert duplicate_inverse(difference) == difference
    return {
        "coefficient_field": "QQ; scalar extension to every characteristic-zero field",
        "pre_S_cubic_domain_dimension": 286,
        "pre_S_cubic_rank": old_rank,
        "pre_S_cubic_kernel_dimension": len(helpers),
        "helper_basis": [pack(q) for q in helpers],
        "helper_weights_for_elimination_only": helper_weights,
        "S_helper_index": s_index,
        "exact_third_fibre_equalities": 63,
        "inverse_composition_coordinate_equalities": 64,
        "restart_kernels": rank_records,
        "constant_kernel_dimension": 1,
        "all_frozen_direction_extensions_covered": "Expand a total cubic in frozen variables; apply the degree0..3 kernel bounds coefficientwise.",
        "uniform_support_envelope": {"W_degree_at_most": 7, "ordinary_degree_at_most": 11},
        "resulting_inverse_degree_bound": "max(D,7*B+4), with B=deg G_w>=1",
        "sharp_literal_instance": {"base_dimension": 48, "base_inverse_degree": 29,
                                   "base_jacobian_cube_zero": True,
                                   "r_du_independence_coefficient": str(independence_coefficient),
                                   "W_inverse_line_coefficient": str(c),
                                   "S_inverse_line_factor": "-2*U^4*W^7",
                                   "inverse_degree": 207, "jacobian_index": 4},
        "negative_controls_rejected": ["Add p1^3 to a valid helper",
                                       "Append w^2*S as the restart",
                                       "Replace the sharp degree ceiling by206"],
        "duplicate_S_positive_control": "u^2*(S-S_copy) is allowed and its inverse value is a degree-three frozen difference",
        "scope": "All finite parallel pre-S cubic helper families, then one further parallel layer of arbitrary homogeneous cubics in the displayed old coordinates and helpers. No further-round or unrestricted-index-four conclusion."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=HERE.parents[2])
    parser.add_argument("--output", type=Path, default=HERE / (STEM + "_verification.json"))
    args = parser.parse_args()
    start = time.monotonic()
    ledger_path = HERE / (STEM + "_source_ledger.json")
    ledger, paths = source_pins(args.repo_root.resolve(), ledger_path)
    checks = run(paths)
    receipt = {"schema": 1, "status": "PASS", "python": platform.python_version(),
               "arithmetic": "Exact rational sparse polynomials and complete column elimination; standard library only.",
               "source_ledger_sha256": sha(ledger_path), "sources": ledger["sources"],
               "checks": checks, "assurance": ASSURANCE,
               "elapsed_seconds": round(time.monotonic() - start, 6)}
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print("PASS: complete21-helper space; restart kernels12/89/522; every cubic ignores S; inverse W-degree<=7 and totaldegree<=11; sharp207; duplicate-helper control.")
    print(args.output)


if __name__ == "__main__":
    main()
