#!/usr/bin/env python3
"""Exact algebraic inverse slice and compressed-jet controls, using only stdlib."""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import math
import time
from fractions import Fraction
from pathlib import Path


if not __debug__:
    raise SystemExit("Refusing optimized Python: exact assertions are required.")

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_sources(ledger_path):
    ledger = json.loads(ledger_path.read_text())
    assert ledger["schema_version"] == 1
    expected = {
        "proof", "checker", "seed_input", "arithmetic_module",
        "jet_proof_context", "shuffle_dictionary", "primary_seed_construction",
    }
    assert {entry["role"] for entry in ledger["sources"]} == expected
    assert len(ledger["sources"]) == len(expected)
    paths = {}
    for entry in ledger["sources"]:
        path = (HERE / entry["path"]).resolve()
        assert path.is_file() and sha(path) == entry["sha256"], entry["path"]
        paths[entry["role"]] = path
    assert paths["checker"] == Path(__file__).resolve()
    spec = importlib.util.spec_from_file_location(
        "pinned_index17_arithmetic", paths["arithmetic_module"]
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    seed = json.loads(paths["seed_input"].read_text())
    module.validate(seed)
    return module, seed, paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ledger", type=Path,
        default=HERE / "cubic_index_seventeen_inverse_series_source_ledger.json",
    )
    parser.add_argument(
        "--output", type=Path,
        default=HERE / "cubic_index_seventeen_inverse_series_verification.json",
    )
    args = parser.parse_args()
    started = time.monotonic()
    m, seed, paths = read_sources(args.ledger)
    terms = seed["H_terms"]
    add, mul = m.add, m.mul

    def term(*indices, coefficient=1):
        return {tuple(sorted(indices)): Fraction(coefficient)}

    def scale(poly, coefficient):
        return m.clean({monomial: value * coefficient
                        for monomial, value in poly.items()})

    def power(poly, exponent):
        result = term()
        for _ in range(exponent):
            result = mul(result, poly)
        return result

    def substitute(poly, images):
        result = {}
        for monomial, coefficient in poly.items():
            value = term(coefficient=coefficient)
            for variable in monomial:
                value = mul(value, images[variable])
            result = add(result, value)
        return result

    def evaluate_map(records, images):
        result = [{} for _ in range(18)]
        for record in records:
            monomial = tuple(record["inputs"])
            value = substitute(
                {monomial: Fraction(record["coefficient"])}, images
            )
            row = record["output"]
            result[row] = add(result[row], value)
        return result

    # Two independent formal variables: index0=t, index1=A.
    one, t, a = term(), term(0), term(1)
    t2 = power(t, 2)
    r = add(one, scale(mul(t2, a), 2))
    discriminant = add(one, scale(t2, -4))
    relation = add(
        mul(discriminant, mul(a, add(one, mul(t2, a)))),
        scale(one, -1),
    )

    def leading(poly):
        # Lexicographic A,t order; a singleton polynomial is a Groebner basis
        # of its principal ideal, so exact division certifies membership.
        return max(poly, key=lambda monomial:
                   (monomial.count(1), monomial.count(0)))

    def divide(poly, divisor):
        work, quotient, remainder = dict(poly), {}, {}
        lead = leading(divisor)
        lead_exponents = (lead.count(0), lead.count(1))
        while work:
            monomial = leading(work)
            coefficient = work[monomial]
            exponents = (monomial.count(0), monomial.count(1))
            if all(a0 >= b0 for a0, b0 in zip(exponents, lead_exponents)):
                difference = tuple(a0 - b0 for a0, b0 in
                                   zip(exponents, lead_exponents))
                multiplier = term(
                    *([0] * difference[0] + [1] * difference[1]),
                    coefficient=coefficient / divisor[lead],
                )
                quotient = add(quotient, multiplier)
                work = add(work, scale(mul(multiplier, divisor), -1))
            else:
                remainder[monomial] = coefficient
                del work[monomial]
        assert add(mul(quotient, divisor), remainder) == poly
        return quotient, remainder

    x = [{} for _ in range(12)]
    x[0], x[1] = a, one
    x[2] = scale(mul(t, mul(add(scale(r, 6), scale(one, 2)),
                           discriminant)), -1)
    x[3] = add(scale(mul(mul(t2, power(a, 2)), x[2]), -1),
               scale(mul(t, a), -3))
    x[4], x[5] = scale(t2, -1), scale(mul(t, power(a, 2)), -1)
    x[6], x[7] = scale(mul(mul(t, a), x[2]), -1), scale(t, -1)
    x[8], x[9] = mul(power(t, 3), a), mul(t2, a)
    x[10] = scale(mul(power(t, 4), power(a, 2)), -2)
    x[11] = add(scale(mul(t, a), -1),
                scale(mul(power(t, 3), power(a, 2)), -1))
    cubic_parts = [{} for _ in range(5)]
    for record in terms:
        if 12 <= record["output"] < 17:
            value = substitute(
                {tuple(record["inputs"]): Fraction(record["coefficient"])}, x
            )
            row = record["output"] - 12
            cubic_parts[row] = add(cubic_parts[row], value)
    inverse_slice = x + [scale(mul(t, c), -1) for c in cubic_parts] + [one]
    target = [one, one] + [{} for _ in range(15)] + [one]

    def slice_residual(records, images):
        evaluated = evaluate_map(records, images)
        return [add(images[i], mul(t, evaluated[i]), scale(target[i], -1))
                for i in range(18)]

    residual = slice_residual(terms, inverse_slice)
    divisions = [divide(poly, relation) for poly in residual]
    assert all(not remainder for _, remainder in divisions)
    assert [m.ev(poly, [Fraction(0), Fraction(1)]) for poly in inverse_slice] \
        == [m.ev(poly, [Fraction(0), Fraction(1)]) for poly in target]

    # The nilpotency identity is recomputed from every literal cubic.
    jacobian = m.jacobian(terms)
    matrix_power = jacobian
    sixteenth = None
    for exponent in range(2, 18):
        matrix_power = m.mmul(matrix_power, jacobian)
        if exponent == 16:
            sixteenth = matrix_power
    assert not matrix_power and sixteenth
    witness = seed["sixteenth_power_witness"]
    expected = {tuple(witness["inputs"]): Fraction(witness["coefficient"])}
    assert sixteenth[(3, 17)] == expected
    assert all(column == 17 for row, column in sixteenth)
    assert all(index in {0, 1, 2, 17} for monomial in expected
               for index in monomial)

    # A direct recursive inverse prefix is independent of the algebraic formula.
    limit = 32
    series = [[Fraction(0)] * (limit + 1) for _ in range(18)]
    for index in (0, 1, 17):
        series[index][0] = Fraction(1)
    for degree in range(1, limit + 1):
        for record in terms:
            i, j, k = record["inputs"]
            coefficient = Fraction(record["coefficient"])
            value = sum(
                (series[i][u] * series[j][v] * series[k][degree - 1 - u - v]
                 for u in range(degree)
                 for v in range(degree - u)),
                Fraction(0),
            )
            series[record["output"]][degree] -= coefficient * value
    central = [Fraction(math.comb(2*j + 2, j + 1), 2)
               for j in range(limit // 2 + 1)]
    for degree in range(limit + 1):
        if degree % 2 == 0:
            assert series[0][degree] == central[degree // 2] > 0
            assert series[5][degree] == 0
        else:
            j = (degree - 1) // 2
            convolution = sum((central[i] * central[j-i]
                               for i in range(j + 1)), Fraction(0))
            assert series[5][degree] == -convolution < 0
            assert series[0][degree] == 0

    jet_checks = []
    for order in range(5):
        # Full jet coordinates use index old+18*level. Retain only the
        # four indicated initial coordinates and the first17 in later levels.
        retained = [0, 1, 2, 17] + [
            index + 18 * level for level in range(1, order + 1)
            for index in range(17)
        ]
        renumber = {old: new for new, old in enumerate(retained)}
        assert len(retained) == 17 * order + 4
        reduced_terms = []
        for record in m.jet_terms(terms, 18, order):
            if record["output"] not in renumber:
                continue
            if not all(index in renumber for index in record["inputs"]):
                continue
            reduced_terms.append({
                "output": renumber[record["output"]],
                "inputs": sorted(renumber[index] for index in record["inputs"]),
                "coefficient": record["coefficient"],
            })
        reduced_jacobian = m.jacobian(reduced_terms)
        old_operator = m.jet_operator(jacobian, 18, order)
        reduced_operator = {}
        for (row, column), poly in old_operator.items():
            if row not in renumber or column not in renumber:
                continue
            replacement = {
                tuple(sorted(renumber[index] for index in monomial)): value
                for monomial, value in poly.items()
                if all(index in renumber for index in monomial)
            }
            if replacement:
                reduced_operator[(renumber[row], renumber[column])] = replacement
        assert reduced_jacobian == reduced_operator
        levels = [old // 18 for old in retained]
        for record in reduced_terms:
            out_level = levels[record["output"]]
            in_levels = [levels[index] for index in record["inputs"]]
            assert len(in_levels) == 3 and max(in_levels) < out_level
            assert sum(2 * level + 1 for level in in_levels) == 2 * out_level + 1
        jet_checks.append({
            "order": order, "dimension": len(retained),
            "cubic_monomials": len(reduced_terms),
            "jacobian_entries": len(reduced_jacobian),
            "coefficientwise_operator_identity": True,
        })

    rejected = []
    wrong_relation = add(
        mul(add(one, scale(t2, 4)), mul(a, add(one, mul(t2, a)))),
        scale(one, -1),
    )
    assert any(divide(poly, wrong_relation)[1] for poly in residual)
    rejected.append("wrong_discriminant_sign")
    missing_helper = list(inverse_slice)
    missing_helper[5] = {}
    assert any(divide(poly, relation)[1]
               for poly in slice_residual(terms, missing_helper))
    rejected.append("omitted_x6_helper")
    wrong_scale = list(inverse_slice)
    wrong_scale[10] = scale(wrong_scale[10], Fraction(1, 2))
    assert any(divide(poly, relation)[1]
               for poly in slice_residual(terms, wrong_scale))
    rejected.append("wrong_scale11")
    perturbed = copy.deepcopy(terms)
    perturbed[0]["coefficient"] = str(Fraction(perturbed[0]["coefficient"]) + 1)
    assert any(divide(poly, relation)[1]
               for poly in slice_residual(perturbed, inverse_slice))
    rejected.append("perturbed_literal_cubic")
    assert not {(i, j): poly for (i, j), poly in sixteenth.items() if j != 17}
    rejected.append("freezing_tau_removes_sixteenth_power_witness")

    receipt = {
        "status": "PASS",
        "arithmetic": "exact rational sparse polynomials; Python standard library",
        "source_ledger_sha256": sha(args.ledger),
        "checker_sha256": sha(Path(__file__)),
        "seed_input_sha256": sha(paths["seed_input"]),
        "literal_cubic_monomials": len(terms),
        "all_order_inverse_equations": len(divisions),
        "all_order_nonzero_identity_quotients": sum(bool(q) for q, _ in divisions),
        "all_order_identity_max_quotient_terms": max(len(q) for q, _ in divisions),
        "formal_branch_constant_verified": True,
        "seed_jacobian_seventeenth_power_zero": True,
        "seed_jacobian_sixteenth_power_entries": len(sixteenth),
        "restricted_sixteenth_power_witness": witness,
        "inverse_prefix_controls_through_grade": limit,
        "first_eight_even_coefficients": [str(value) for value in central[:8]],
        "first_eight_odd_coefficients": [str(series[5][2*j+1]) for j in range(8)],
        "compressed_jet_controls": jet_checks,
        "negative_controls_rejected": rejected,
        "universal_claim_support": (
            "Polynomial inverse identities plus the written binomial positivity, "
            "weight, coefficientwise operator and retained-column arguments; "
            "the finite prefix and small jets are controls only."
        ),
        "seconds": round(time.monotonic() - started, 6),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
