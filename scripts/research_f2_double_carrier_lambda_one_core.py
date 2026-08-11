#!/usr/bin/env python3
"""Extract the unresolved lambda=1 polynomial core of the F2 formal branch.

This is an exploratory finite-field reduction.  It reads a long formal jet,
certifies candidate linear recurrences on a withheld suffix, evaluates only
the recurrence coordinates whose denominator is nonzero at lambda=1, and
keeps every other active coordinate symbolic.  The complete 367-equation
arithmetic circuit is then specialized to those symbolic coordinates.

A zero-dimensional solution of the emitted core would be an exact point of
one finite-field reduction only.  Recurrence certification from finitely many
terms is not a proof that the formal branch specializes, and emptiness of one
chosen affine slice is not emptiness of the full characteristic-zero ideal.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
CAS = ROOT / "plane-jc" / "cas"
sys.path.insert(0, str(CAS))

from probe_f2_75_125_nonlinear_modular import (  # noqa: E402
    build_modular_presentation,
    initial_localized_point,
)
from sparse_circuit_modp import (  # noqa: E402
    evaluate,
    quadratic_key_embedding,
    solve_linearization,
)


def berlekamp_massey(sequence: list[int], prime: int) -> list[int]:
    """Return C with s_n+sum_i C_i*s_(n-i)=0 over GF(prime)."""

    connection = [1]
    previous = [1]
    length = 0
    shift = 1
    previous_discrepancy = 1
    for index in range(len(sequence)):
        discrepancy = sequence[index]
        discrepancy += sum(
            connection[offset] * sequence[index - offset]
            for offset in range(1, length + 1)
        )
        discrepancy %= prime
        if discrepancy == 0:
            shift += 1
            continue
        saved = connection[:]
        scalar = discrepancy * pow(previous_discrepancy, -1, prime) % prime
        required = len(previous) + shift
        if len(connection) < required:
            connection += [0] * (required - len(connection))
        for offset, value in enumerate(previous):
            connection[offset + shift] = (
                connection[offset + shift] - scalar * value
            ) % prime
        if 2 * length <= index:
            length = index + 1 - length
            previous = saved
            previous_discrepancy = discrepancy
            shift = 1
        else:
            shift += 1
    while len(connection) > 1 and connection[-1] % prime == 0:
        connection.pop()
    return connection


def recurrence_holds(
    connection: list[int], sequence: list[int], start: int, prime: int
) -> bool:
    length = len(connection) - 1
    return all(
        (
            sequence[index]
            + sum(
                connection[offset] * sequence[index - offset]
                for offset in range(1, length + 1)
            )
        )
        % prime
        == 0
        for index in range(max(start, length), len(sequence))
    )


def recurrence_value_at_one(
    connection: list[int], sequence: list[int], train: int, prime: int
) -> int:
    """Evaluate the candidate rational generating function at lambda=1."""

    denominator = sum(connection) % prime
    if denominator == 0:
        raise ZeroDivisionError("the candidate recurrence has a pole at lambda=1")
    numerator_sum = 0
    for degree in range(train):
        numerator_sum += sum(
            connection[offset] * sequence[degree - offset]
            for offset in range(min(degree, len(connection) - 1) + 1)
        )
    return numerator_sum % prime * pow(denominator, -1, prime) % prime


def polynomial_add(
    left: dict[int, int], right: dict[int, int], prime: int
) -> dict[int, int]:
    if not left:
        return right.copy()
    if not right:
        return left.copy()
    result = left.copy()
    for monomial, coefficient in right.items():
        value = (result.get(monomial, 0) + coefficient) % prime
        if value:
            result[monomial] = value
        else:
            result.pop(monomial, None)
    return result


def polynomial_scale(
    scalar: int, polynomial: dict[int, int], prime: int
) -> dict[int, int]:
    scalar %= prime
    if scalar == 0 or not polynomial:
        return {}
    if scalar == 1:
        return polynomial.copy()
    return {
        monomial: scalar * coefficient % prime
        for monomial, coefficient in polynomial.items()
    }


def polynomial_multiply(
    left: dict[int, int],
    right: dict[int, int],
    prime: int,
    maximum_terms: int,
) -> dict[int, int]:
    if not left or not right:
        return {}
    if len(left) == 1 and 0 in left:
        return polynomial_scale(left[0], right, prime)
    if len(right) == 1 and 0 in right:
        return polynomial_scale(right[0], left, prime)
    if len(left) * len(right) > maximum_terms * 8:
        raise RuntimeError(
            "symbolic multiplication forecast exceeds the configured term cap: "
            f"{len(left)}*{len(right)}"
        )
    result: dict[int, int] = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            # Exponents use disjoint three-bit fields.  Total circuit degree
            # is at most seven, so ordinary integer addition cannot carry
            # between fields.
            monomial = left_monomial + right_monomial
            value = (
                result.get(monomial, 0)
                + left_coefficient * right_coefficient
            ) % prime
            if value:
                result[monomial] = value
            else:
                result.pop(monomial, None)
    if len(result) > maximum_terms:
        raise RuntimeError(
            "a specialized circuit polynomial exceeded the configured term cap: "
            f"{len(result)}>{maximum_terms}"
        )
    return result


def monomial_degree(monomial: int, variable_count: int) -> int:
    return sum((monomial >> (3 * index)) & 7 for index in range(variable_count))


def monomial_support(monomial: int, variable_count: int) -> tuple[int, ...]:
    return tuple(
        index
        for index in range(variable_count)
        if (monomial >> (3 * index)) & 7
    )


def solve_linear_core(
    equations: list[dict[int, int]], variable_count: int, prime: int
) -> dict[str, object]:
    """Row-reduce a specialized affine-linear core over ``GF(prime)``."""

    augmented_column = variable_count
    rows: list[dict[int, int]] = []
    certificates: list[dict[int, int]] = []
    for equation_index, polynomial in enumerate(equations):
        row: dict[int, int] = {}
        for monomial, coefficient in polynomial.items():
            coefficient %= prime
            if not coefficient:
                continue
            if monomial == 0:
                row[augmented_column] = (-coefficient) % prime
                continue
            support = monomial_support(monomial, variable_count)
            if len(support) != 1 or monomial_degree(monomial, variable_count) != 1:
                raise ValueError("the specialized core is not affine-linear")
            index = support[0]
            row[index] = (row.get(index, 0) + coefficient) % prime
            if not row[index]:
                row.pop(index)
        rows.append(row)
        certificates.append({equation_index: 1})

    pivots: dict[int, tuple[dict[int, int], dict[int, int]]] = {}
    inconsistent_rows: list[int] = []
    inconsistency_certificates: list[dict[int, int]] = []
    for row_index, (row, certificate) in enumerate(zip(rows, certificates)):
        while True:
            columns = sorted(
                column
                for column, value in row.items()
                if column < variable_count and value
            )
            if not columns:
                if row.get(augmented_column, 0):
                    inconsistent_rows.append(row_index)
                    inverse = pow(row[augmented_column], -1, prime)
                    normalized_certificate = {
                        index: coefficient * inverse % prime
                        for index, coefficient in certificate.items()
                        if coefficient * inverse % prime
                    }
                    # The normalized row is 0=1, so the same combination of
                    # the original equations is the constant -1.
                    combined: dict[int, int] = {}
                    for index, coefficient in normalized_certificate.items():
                        for monomial, value in equations[index].items():
                            new_value = (
                                combined.get(monomial, 0)
                                + coefficient * value
                            ) % prime
                            if new_value:
                                combined[monomial] = new_value
                            else:
                                combined.pop(monomial, None)
                    if combined != {0: prime - 1}:
                        raise AssertionError(
                            "the affine-core inconsistency certificate failed"
                        )
                    inconsistency_certificates.append(normalized_certificate)
                break
            column = columns[0]
            pivot = pivots.get(column)
            if pivot is None:
                inverse = pow(row[column], -1, prime)
                normalized_row = {
                    key: value * inverse % prime
                    for key, value in row.items()
                    if value * inverse % prime
                }
                normalized_certificate = {
                    key: value * inverse % prime
                    for key, value in certificate.items()
                    if value * inverse % prime
                }
                pivots[column] = (normalized_row, normalized_certificate)
                break
            scalar = row[column]
            pivot_row, pivot_certificate = pivot
            for key, value in pivot_row.items():
                combined = (row.get(key, 0) - scalar * value) % prime
                if combined:
                    row[key] = combined
                else:
                    row.pop(key, None)
            for key, value in pivot_certificate.items():
                combined = (certificate.get(key, 0) - scalar * value) % prime
                if combined:
                    certificate[key] = combined
                else:
                    certificate.pop(key, None)

    solution = [0] * variable_count
    if not inconsistent_rows:
        for column in sorted(pivots, reverse=True):
            row, _certificate = pivots[column]
            solution[column] = (
                row.get(augmented_column, 0)
                - sum(
                    coefficient * solution[index]
                    for index, coefficient in row.items()
                    if index < variable_count and index != column
                )
            ) % prime
    rank = len(pivots)
    return {
        "coefficient_rank": rank,
        "augmented_rank": rank + bool(inconsistent_rows),
        "inconsistent_rows": inconsistent_rows,
        "inconsistency_certificates": inconsistency_certificates,
        "solution_dimension_if_consistent": (
            variable_count - rank if not inconsistent_rows else None
        ),
        "zero_free_parameter_solution": solution if not inconsistent_rows else None,
    }


def specialize_circuit(
    *,
    prime: int,
    rho: int,
    constants: dict[str, int],
    symbolic_names: tuple[str, ...],
    maximum_terms: int,
    presentation: tuple | None = None,
) -> tuple[list[dict[int, int]], list[str], dict[str, tuple[int, int]]]:
    dag, roots, groups, _a_node, labels = (
        build_modular_presentation() if presentation is None else presentation
    )
    embedding = quadratic_key_embedding(prime, rho)
    symbolic_index = {name: index for index, name in enumerate(symbolic_names)}

    reachable: set[int] = set(roots)
    stack = list(roots)
    while stack:
        index = stack.pop()
        node = dag.nodes[index]
        children: tuple[int, ...]
        if node[0] == "scale":
            children = (int(node[2]),)
        elif node[0] in ("add", "mul"):
            children = (int(node[1]), int(node[2]))
        else:
            children = ()
        for child in children:
            if child not in reachable:
                reachable.add(child)
                stack.append(child)

    references = Counter(roots)
    for index in reachable:
        node = dag.nodes[index]
        if node[0] == "scale":
            references[int(node[2])] += 1
        elif node[0] in ("add", "mul"):
            references[int(node[1])] += 1
            references[int(node[2])] += 1

    values: list[dict[int, int] | None] = [None] * len(dag.nodes)
    maximum_observed = 0
    for index, node in enumerate(dag.nodes):
        if index not in reachable:
            continue
        kind = node[0]
        children = ()
        if kind == "const":
            coefficient = embedding(node[1]) % prime
            polynomial = {0: coefficient} if coefficient else {}
        elif kind == "var":
            name = str(node[1])
            if name in symbolic_index:
                polynomial = {1 << (3 * symbolic_index[name]): 1}
            else:
                coefficient = constants.get(name, 0) % prime
                polynomial = {0: coefficient} if coefficient else {}
        elif kind == "scale":
            child = int(node[2])
            children = (child,)
            assert values[child] is not None
            polynomial = polynomial_scale(embedding(node[1]), values[child], prime)
        elif kind == "add":
            left, right = int(node[1]), int(node[2])
            children = (left, right)
            assert values[left] is not None and values[right] is not None
            polynomial = polynomial_add(values[left], values[right], prime)
        elif kind == "mul":
            left, right = int(node[1]), int(node[2])
            children = (left, right)
            assert values[left] is not None and values[right] is not None
            polynomial = polynomial_multiply(
                values[left], values[right], prime, maximum_terms
            )
        else:
            raise AssertionError(f"unknown circuit operation {kind!r}")
        values[index] = polynomial
        maximum_observed = max(maximum_observed, len(polynomial))
        for child in children:
            references[child] -= 1
            if references[child] == 0:
                values[child] = None

    equations = []
    for root in roots:
        assert values[root] is not None
        equations.append(values[root])
        references[root] -= 1
    print(f"SPECIALIZED_DAG_REACHABLE={len(reachable)}")
    print(f"SPECIALIZED_MAX_INTERMEDIATE_TERMS={maximum_observed}")
    return equations, labels, groups


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--jet",
        type=Path,
        default=ROOT / "tmp/f2_formal_homotopy_polefree3_order64.json",
    )
    parser.add_argument("--train", type=int, default=44)
    parser.add_argument("--maximum-terms", type=int, default=250_000)
    parser.add_argument(
        "--force-symbolic",
        action="append",
        default=[],
        help="keep a nominally stable recurrence coordinate symbolic (repeatable)",
    )
    parser.add_argument(
        "--exact-affine-slice",
        action="store_true",
        help=(
            "ignore recurrence values, fix the deterministic nonpivot/gauge "
            "coordinates to their affine lambda values, and retain every "
            "fixed-Jacobian pivot coordinate symbolically"
        ),
    )
    parser.add_argument("--artifact", type=Path)
    args = parser.parse_args()

    jet = json.loads(args.jet.read_text())
    prime = int(jet["field"]["prime"])
    rho = int(jet["field"]["rho"])
    series_by_name = {
        str(name): [int(value) % prime for value in values]
        for name, values in jet["variable_series"].items()
    }
    if not 1 < args.train < len(next(iter(series_by_name.values()))):
        raise ValueError("the train length must leave a nonempty withheld suffix")

    constants: dict[str, int] = {}
    unresolved = []
    recurrence_ledger = {}
    for name, sequence in sorted(series_by_name.items()):
        connection = berlekamp_massey(sequence[: args.train], prime)
        stable = recurrence_holds(connection, sequence, args.train, prime)
        denominator_at_one = sum(connection) % prime
        recurrence_ledger[name] = {
            "order": len(connection) - 1,
            "connection_low_to_high": connection,
            "withheld_term_count": len(sequence) - args.train,
            "withheld_terms_pass": stable,
            "denominator_at_lambda_one": denominator_at_one,
        }
        if stable and denominator_at_one and name not in args.force_symbolic:
            value = recurrence_value_at_one(
                connection, sequence, args.train, prime
            )
            constants[name] = value
            recurrence_ledger[name]["candidate_value_at_lambda_one"] = value
        else:
            unresolved.append(name)

    presentation = None
    slice_mode = "withheld_recurrence_candidate"
    if args.exact_affine_slice:
        slice_mode = "exact_fixed_jacobian_affine_slice"
        presentation = build_modular_presentation()
        dag, roots, _groups, _a_node, _labels = presentation
        seed = initial_localized_point(
            prime, rho, int(jet["field"]["y"])
        )
        evaluated = evaluate(
            dag,
            seed,
            prime,
            quadratic_key_embedding(prime, rho),
            with_jacobian=True,
        )
        gauge_variables = tuple(
            str(name)
            for name in jet["higher_order_gauge"][
                "variables_prescribed_zero_from_order_two"
            ]
        )
        pivot_solve = solve_linearization(
            evaluated,
            roots,
            prime,
            right_hand_side=[0] * len(roots),
            prescribed_values={name: 0 for name in gauge_variables},
        )
        if pivot_solve.inconsistent_rows:
            raise AssertionError("the exact affine-slice Jacobian chart changed")
        symbolic_names = tuple(sorted(pivot_solve.pivot_variables))
        symbolic_set = set(symbolic_names)
        constants = {
            name: sum(series_by_name.get(name, [0, 0])[:2]) % prime
            for name in evaluated.variable_names
            if name not in symbolic_set
        }
        unresolved = list(symbolic_names)
    else:
        symbolic_names = tuple(unresolved)
    equations, labels, groups = specialize_circuit(
        prime=prime,
        rho=rho,
        constants=constants,
        symbolic_names=symbolic_names,
        maximum_terms=args.maximum_terms,
        presentation=presentation,
    )

    equation_ledger = []
    variable_occurrences = Counter()
    for label, polynomial in zip(labels, equations):
        degree = max(
            (monomial_degree(monomial, len(symbolic_names)) for monomial in polynomial),
            default=-1,
        )
        support = set()
        for monomial in polynomial:
            support.update(monomial_support(monomial, len(symbolic_names)))
        for index in support:
            variable_occurrences[symbolic_names[index]] += 1
        equation_ledger.append(
            {
                "label": label,
                "term_count": len(polynomial),
                "total_degree": degree,
                "variable_count": len(support),
                "variables": [symbolic_names[index] for index in sorted(support)],
                "terms": [
                    [monomial, coefficient]
                    for monomial, coefficient in sorted(polynomial.items())
                ],
            }
        )

    maximum_equation_degree = max(
        (row["total_degree"] for row in equation_ledger), default=-1
    )
    linear_core = (
        solve_linear_core(equations, len(symbolic_names), prime)
        if maximum_equation_degree <= 1
        else None
    )
    if linear_core:
        linear_core["labelled_inconsistency_certificates"] = [
            [
                {
                    "equation_index": equation_index,
                    "label": labels[equation_index],
                    "coefficient": coefficient,
                }
                for equation_index, coefficient in sorted(certificate.items())
            ]
            for certificate in linear_core["inconsistency_certificates"]
        ]

    exact_candidate_replay = None
    if linear_core and not linear_core["inconsistent_rows"]:
        solution = linear_core["zero_free_parameter_solution"]
        assert isinstance(solution, list)
        candidate = constants.copy()
        candidate.update(
            {
                name: value
                for name, value in zip(symbolic_names, solution)
                if value
            }
        )
        dag, roots, _groups, _a_node, _labels = (
            build_modular_presentation() if presentation is None else presentation
        )
        evaluated_candidate = evaluate(
            dag,
            candidate,
            prime,
            quadratic_key_embedding(prime, rho),
            with_jacobian=False,
        )
        residuals = [evaluated_candidate.values[root] for root in roots]
        exact_candidate_replay = {
            "point_support_size": len(candidate),
            "nonzero_residual_count": sum(bool(value) for value in residuals),
            "is_exact_modular_point": not any(residuals),
        }

    payload = {
        "schema": "jc2-f2-double-carrier-lambda-one-core-v1",
        "status": "exploratory_finite_field_slice",
        "slice_mode": slice_mode,
        "source_jet": str(args.jet),
        "field": {"prime": prime, "rho": rho},
        "recurrence_train_terms": args.train,
        "recurrence_withheld_terms": len(next(iter(series_by_name.values())))
        - args.train,
        "candidate_constant_coordinate_count": len(constants),
        "symbolic_coordinate_count": len(symbolic_names),
        "symbolic_coordinates": list(symbolic_names),
        "symbolic_coordinate_equation_occurrences": dict(
            sorted(variable_occurrences.items())
        ),
        "recurrences": recurrence_ledger,
        "equation_groups": groups,
        "equation_count": len(equations),
        "nonzero_equation_count": sum(bool(polynomial) for polynomial in equations),
        "equation_degree_histogram": dict(
            sorted(
                Counter(row["total_degree"] for row in equation_ledger).items()
            )
        ),
        "equation_term_count_histogram": dict(
            sorted(Counter(row["term_count"] for row in equation_ledger).items())
        ),
        "linear_core": linear_core,
        "exact_candidate_replay": exact_candidate_replay,
        "equations": equation_ledger,
        "claim_boundary": (
            "finite recurrence prefixes and one finite-field affine slice do "
            "not prove specialization, characteristic-zero existence, or emptiness"
        ),
    }
    if args.artifact:
        args.artifact.parent.mkdir(parents=True, exist_ok=True)
        args.artifact.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"WROTE {args.artifact}")
    print(f"LAMBDA_ONE_STABLE_COORDINATES={len(constants)}")
    print(f"LAMBDA_ONE_SYMBOLIC_COORDINATES={len(symbolic_names)}")
    print(f"LAMBDA_ONE_NONZERO_EQUATIONS={payload['nonzero_equation_count']}")
    print(f"LAMBDA_ONE_DEGREE_HISTOGRAM={payload['equation_degree_histogram']}")
    if linear_core:
        print(
            "LAMBDA_ONE_LINEAR_RANKS="
            f"{linear_core['coefficient_rank']},"
            f"{linear_core['augmented_rank']}"
        )
        print(
            "LAMBDA_ONE_LINEAR_SOLUTION_DIMENSION="
            f"{linear_core['solution_dimension_if_consistent']}"
        )
    if exact_candidate_replay:
        print(
            "LAMBDA_ONE_EXACT_CANDIDATE_REPLAY="
            f"{exact_candidate_replay['is_exact_modular_point']}"
        )


if __name__ == "__main__":
    main()
