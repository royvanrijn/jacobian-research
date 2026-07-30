#!/usr/bin/env python3
"""Completed-invariant zero-fiber probes for bidegree (3,3).

This research checker joins two previously separate calculations:

* the Hilbert/Jacobian comparison between corrected moment parameters and
  systems augmented by the quadratic Casimirs q_2 and q_4;
* modular zero-fiber probes after q_2=0, on the normalized null-quadratic
  chart F_2=X^2 and on its boundary F_2=0.

The Hilbert coefficients and nonzero modular Jacobian ranks are exact
characteristic-zero certificates.  The Singular radical probes are only
finite-field evidence unless independently reconstructed over QQ.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
import json
from math import factorial, gcd
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from research_completed_moment_algebra import (
    casimir_projectors_mod,
    component_matrices_mod,
    deterministic_point,
    hilbert_numerator,
    invariant_values_mod,
    invariant_hilbert_coefficients,
    moments_mod,
    moment_jacobian_mod,
    monomial_value,
    rank_mod,
    weighted_exponents,
)
from explore_two_pair_sic_bidegree33_full_anchor import (
    PARAMETERS,
    Q_POLYNOMIALS,
    WEIGHTS,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_casimir_fiber.json"
)
DEFAULT_PRIME = 32003
QUOTIENT_DIMENSION = 13
AMBIENT_DIMENSION = 16
MOMENT_CUTOFF = 24
HILBERT_CUTOFF = 120

# The highest-weight quadratic F_2=X^2 has phase +1.  After removing the
# unmatched biform variable, its coefficient polynomial is (1+q)^2.
NULL_QUADRATIC_WEIGHT = 1
NULL_QUADRATIC_Q = (1, 2, 1)

# On F_2=X^2, synchronization with the same destabilizing root is the
# coordinate space s3=...=s6=t2=t3=t4=0.
NULL_QUADRATIC_FORBIDDEN = (3, 4, 5, 6, 9, 10, 11)


def convolve_mod(
    left: tuple[int, ...],
    right: tuple[int, ...],
    prime: int,
) -> tuple[int, ...]:
    answer = [0] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        if left_coefficient % prime == 0:
            continue
        for right_index, right_coefficient in enumerate(right):
            answer[left_index + right_index] = (
                answer[left_index + right_index]
                + left_coefficient * right_coefficient
            ) % prime
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def polynomial_powers_mod(
    polynomial: tuple[int, ...],
    maximum: int,
    prime: int,
) -> tuple[tuple[int, ...], ...]:
    powers = [(1,)]
    for _ in range(maximum):
        powers.append(convolve_mod(powers[-1], polynomial, prime))
    return tuple(powers)


def restricted_moment_terms_mod(
    order: int,
    quadratic_mode: str,
    prime: int,
) -> dict[tuple[int, ...], int]:
    """Return one restricted moment as a sparse polynomial modulo prime."""

    assert quadratic_mode in {"null", "zero"}
    factorials = [
        factorial(index) % prime for index in range(3 * order + 1)
    ]
    inverse_factorials = [
        pow(factorial(index) % prime, -1, prime)
        for index in range(order + 1)
    ]
    basis_powers = tuple(
        polynomial_powers_mod(polynomial, order, prime)
        for polynomial in Q_POLYNOMIALS
    )
    quadratic_powers = polynomial_powers_mod(
        NULL_QUADRATIC_Q, order, prime
    )
    parameter_order = (0, 6, 1, 5, 7, 11, 2, 4, 8, 10, 3, 9)
    exponents = [0] * len(PARAMETERS)
    answer: dict[tuple[int, ...], int] = defaultdict(int)
    order_factorial = factorials[order]

    @lru_cache(maxsize=None)
    def remaining_weight_bounds(
        position: int,
        degree_left: int,
    ) -> tuple[int, int]:
        remaining = [
            WEIGHTS[index] for index in parameter_order[position:]
        ]
        if quadratic_mode == "null":
            remaining.append(NULL_QUADRATIC_WEIGHT)
        if degree_left == 0:
            return 0, 0
        if not remaining:
            return 1, 0
        return (
            degree_left * min(remaining),
            degree_left * max(remaining),
        )

    def visit(
        position: int,
        used_degree: int,
        weight: int,
        shift: int,
        inverse_denominator: int,
        q_polynomial: tuple[int, ...],
    ) -> None:
        if position == len(parameter_order):
            quadratic_exponent = order - used_degree
            if quadratic_mode == "zero" and quadratic_exponent:
                return
            total_weight = weight
            product = q_polynomial
            total_shift = shift
            if quadratic_mode == "null":
                total_weight += (
                    NULL_QUADRATIC_WEIGHT * quadratic_exponent
                )
                total_shift += (
                    NULL_QUADRATIC_WEIGHT * quadratic_exponent
                )
                product = convolve_mod(
                    product,
                    quadratic_powers[quadratic_exponent],
                    prime,
                )
            if total_weight != 0:
                return

            scalar = (
                order_factorial
                * inverse_denominator
                * inverse_factorials[quadratic_exponent]
            ) % prime
            contraction = 0
            for q_degree, coefficient in enumerate(product):
                diagonal = total_shift + q_degree
                if 0 <= diagonal <= 3 * order:
                    contraction += (
                        coefficient
                        * factorials[3 * order - diagonal]
                        * factorials[diagonal]
                    )
            coefficient = scalar * contraction % prime
            if coefficient:
                exponent_tuple = tuple(exponents)
                answer[exponent_tuple] = (
                    answer[exponent_tuple] + coefficient
                ) % prime
            return

        parameter_index = parameter_order[position]
        parameter_weight = WEIGHTS[parameter_index]
        available = order - used_degree
        for exponent in range(available + 1):
            new_weight = weight + exponent * parameter_weight
            degree_left = available - exponent
            minimum, maximum = remaining_weight_bounds(
                position + 1, degree_left
            )
            target_weight = -new_weight
            if not minimum <= target_weight <= maximum:
                continue
            exponents[parameter_index] = exponent
            visit(
                position + 1,
                used_degree + exponent,
                new_weight,
                shift + max(parameter_weight, 0) * exponent,
                inverse_denominator
                * inverse_factorials[exponent]
                % prime,
                convolve_mod(
                    q_polynomial,
                    basis_powers[parameter_index][exponent],
                    prime,
                ),
            )
        exponents[parameter_index] = 0

    visit(0, 0, 0, 0, 1, (1,))
    return {
        exponent_tuple: coefficient % prime
        for exponent_tuple, coefficient in answer.items()
        if coefficient % prime
    }


def serialize_modular_polynomial(
    terms: dict[tuple[int, ...], int],
    prime: int,
) -> str:
    pieces = []
    for exponents, coefficient in sorted(terms.items()):
        signed = coefficient if coefficient <= prime // 2 else coefficient - prime
        factors = []
        for variable, exponent in zip(PARAMETERS, exponents, strict=True):
            if exponent == 1:
                factors.append(variable)
            elif exponent > 1:
                factors.append(f"{variable}^{exponent}")
        monomial = "*".join(factors)
        if monomial:
            if signed == 1:
                pieces.append(monomial)
            elif signed == -1:
                pieces.append(f"-{monomial}")
            else:
                pieces.append(f"{signed}*{monomial}")
        else:
            pieces.append(str(signed))
    expression = "+".join(pieces).replace("+-", "-")
    return expression or "0"


def primitive_modular_terms(
    terms: dict[tuple[int, ...], int],
    prime: int,
) -> dict[tuple[int, ...], int]:
    """Scale a nonzero modular polynomial to leading coefficient one."""

    if not terms:
        return terms
    leading = terms[min(terms)]
    inverse = pow(leading, -1, prime)
    return {
        exponents: coefficient * inverse % prime
        for exponents, coefficient in terms.items()
    }


def quadratic_gradient_row_mod(
    point: list[list[int]],
    component: int,
    prime: int,
) -> list[int]:
    """Differentiate q_(2*component)=tr(A_component^2) by evaluation."""

    projectors = casimir_projectors_mod(3, prime)
    size = 4
    values = []
    for row in range(size):
        for column in range(size):
            plus = [entries[:] for entries in point]
            minus = [entries[:] for entries in point]
            plus[row][column] = (plus[row][column] + 1) % prime
            minus[row][column] = (minus[row][column] - 1) % prime
            plus_component = component_matrices_mod(
                plus, 3, projectors, prime
            )[component]
            minus_component = component_matrices_mod(
                minus, 3, projectors, prime
            )[component]
            plus_value = sum(
                plus_component[i][j] * plus_component[j][i]
                for i in range(size)
                for j in range(size)
            ) % prime
            minus_value = sum(
                minus_component[i][j] * minus_component[j][i]
                for i in range(size)
                for j in range(size)
            ) % prime
            values.append((plus_value - minus_value) * pow(2, -1, prime) % prime)
    return values


def degree_system_record(
    label: str,
    moment_orders: tuple[int, ...],
    casimir_components: tuple[int, ...],
    hilbert: list[int],
    moment_rows: list[list[int]],
    quadratic_rows: dict[int, list[int]],
    prime: int,
) -> dict[str, object]:
    degrees = (
        tuple(moment_orders)
        + (2,) * len(casimir_components)
    )
    assert len(degrees) == QUOTIENT_DIMENSION
    predicted_top = sum(degrees) - AMBIENT_DIMENSION
    numerator = hilbert_numerator(hilbert, degrees)
    first_negative = next(
        (
            [index, coefficient]
            for index, coefficient in enumerate(numerator)
            if coefficient < 0
        ),
        None,
    )
    first_tail = next(
        (
            [index, coefficient]
            for index, coefficient in enumerate(
                numerator[predicted_top + 1 :],
                start=predicted_top + 1,
            )
            if coefficient
        ),
        None,
    )
    rows = [moment_rows[order - 1] for order in moment_orders]
    rows.extend(quadratic_rows[component] for component in casimir_components)
    return {
        "label": label,
        "moment_orders": list(moment_orders),
        "casimirs": [f"q_{2 * component}" for component in casimir_components],
        "degrees": list(degrees),
        "total_invariant_degree": sum(degrees),
        "degree_product": (
            product(degrees)
        ),
        "predicted_hilbert_numerator_top_degree": predicted_top,
        "first_negative_through_cutoff": first_negative,
        "first_nonzero_after_predicted_top_through_cutoff": first_tail,
        "hilbert_checked_through_degree": len(hilbert) - 1,
        "numerator_coefficient_sum_through_predicted_top": sum(
            numerator[: predicted_top + 1]
        ),
        "jacobian_rank_mod_prime": rank_mod(rows, prime),
    }


def product(values: tuple[int, ...]) -> int:
    answer = 1
    for value in values:
        answer *= value
    return answer


def degree_comparison(prime: int) -> dict[str, object]:
    hilbert = invariant_hilbert_coefficients(3, HILBERT_CUTOFF)
    point = deterministic_point(3, 0, prime)
    moment_rows = moment_jacobian_mod(
        point, 3, MOMENT_CUTOFF, prime
    )
    quadratic_rows = {
        component: quadratic_gradient_row_mod(point, component, prime)
        for component in (1, 2)
    }

    systems = [
        (
            "corrected_moments",
            tuple(range(1, 13)) + (14,),
            (),
        ),
        (
            "low_degree_q2",
            tuple(range(1, 13)),
            (1,),
        ),
    ]
    records = [
        degree_system_record(
            label,
            moment_orders,
            casimirs,
            hilbert,
            moment_rows,
            quadratic_rows,
            prime,
        )
        for label, moment_orders, casimirs in systems
    ]

    matched_searches = []
    for casimir_components in ((1,), (1, 2)):
        extra_count = QUOTIENT_DIMENSION - len(casimir_components) - 2
        extra_sum = (
            92 - 1 - 2 - 2 * len(casimir_components)
        )
        candidates = []
        degree_sum_matches = 0
        hilbert_matches = 0
        for extra_orders in combinations(
            range(3, MOMENT_CUTOFF + 1), extra_count
        ):
            if sum(extra_orders) != extra_sum:
                continue
            degree_sum_matches += 1
            moment_orders = (1, 2) + extra_orders
            record = degree_system_record(
                "candidate",
                moment_orders,
                casimir_components,
                hilbert,
                moment_rows,
                quadratic_rows,
                prime,
            )
            if (
                record["first_negative_through_cutoff"] is not None
                or record[
                    "first_nonzero_after_predicted_top_through_cutoff"
                ]
                is not None
            ):
                continue
            hilbert_matches += 1
            if record["jacobian_rank_mod_prime"] != QUOTIENT_DIMENSION:
                continue
            candidates.append(record)
        candidates.sort(
            key=lambda record: (
                max(record["moment_orders"]),
                record["degree_product"],
                record["moment_orders"],
            )
        )
        assert candidates
        selected = candidates[0]
        selected["label"] = (
            "matched_total_q2"
            if len(casimir_components) == 1
            else "matched_total_q2_q4"
        )
        records.append(selected)
        matched_searches.append({
            "casimirs": selected["casimirs"],
            "target_total_invariant_degree": 92,
            "moment_search_range": [3, MOMENT_CUTOFF],
            "degree_sum_matches": degree_sum_matches,
            "hilbert_compatible_matches": hilbert_matches,
            "hilbert_compatible_full_rank_matches": len(candidates),
            "selection_rule": (
                "minimize largest moment order, then degree product, "
                "then lexicographic moment-order tuple"
            ),
            "selected_label": selected["label"],
        })
    assert all(
        record["jacobian_rank_mod_prime"] == QUOTIENT_DIMENSION
        for record in records
    )
    return {
        "status": (
            "exact Hilbert coefficient arithmetic and characteristic-zero "
            "Jacobian rank certificates via one good prime"
        ),
        "prime": prime,
        "systems": records,
        "matched_total_degree_searches": matched_searches,
        "mu14_casimir_class_test": mu14_casimir_class_test(prime),
    }


def mu14_casimir_class_test(
    prime: int,
    extra_samples: int = 3,
) -> dict[str, object]:
    """Compare mu_14 and q_2^7 modulo the lower-moment monomial ideal.

    Every degree-14 polynomial in mu_1,...,mu_12,q_2 that vanishes after
    mu_1=...=mu_12=0 is included, not a bounded subset of monomials.
    Full column rank modulo one prime is therefore an exact nonrelation
    certificate over characteristic zero for this generated subalgebra.
    """

    weights = tuple(range(1, 13)) + (2,)
    exponents = weighted_exponents(weights, 14)
    pure_q2 = (0,) * 12 + (7,)
    assert pure_q2 in exponents
    lower_ideal_exponents = [
        exponent_tuple
        for exponent_tuple in exponents
        if exponent_tuple != pure_q2
    ]
    sample_count = len(lower_ideal_exponents) + 2 + extra_samples
    projectors = casimir_projectors_mod(3, prime)
    lower_matrix = []
    q2_column = []
    mu14_column = []
    for sample_index in range(1, sample_count + 1):
        point = deterministic_point(3, 10_000 + sample_index, prime)
        moments = moments_mod(point, 3, 14, prime)
        quadratics, _odd = invariant_values_mod(
            point, 3, projectors, prime
        )
        base_values = moments[:12] + [quadratics[1]]
        lower_matrix.append([
            monomial_value(base_values, exponent_tuple, prime)
            for exponent_tuple in lower_ideal_exponents
        ])
        q2_column.append(pow(quadratics[1], 7, prime))
        mu14_column.append(moments[13])

    rank_lower = rank_mod(lower_matrix, prime)
    rank_lower_q2 = rank_mod(
        [
            row + [q2]
            for row, q2 in zip(
                lower_matrix, q2_column, strict=True
            )
        ],
        prime,
    )
    rank_lower_mu14 = rank_mod(
        [
            row + [mu14]
            for row, mu14 in zip(
                lower_matrix, mu14_column, strict=True
            )
        ],
        prime,
    )
    rank_all = rank_mod(
        [
            row + [q2, mu14]
            for row, q2, mu14 in zip(
                lower_matrix,
                q2_column,
                mu14_column,
                strict=True,
            )
        ],
        prime,
    )
    return {
        "prime": prime,
        "weight": 14,
        "base_generators": [
            *[f"mu_{order}" for order in range(1, 13)],
            "q_2",
        ],
        "degree_14_monomials_total": len(exponents),
        "lower_moment_ideal_columns": len(lower_ideal_exponents),
        "samples": sample_count,
        "ranks": {
            "lower_moment_ideal": rank_lower,
            "plus_q2_power_7": rank_lower_q2,
            "plus_mu14": rank_lower_mu14,
            "plus_both": rank_all,
        },
        "classes_independent_mod_lower_moment_span": (
            rank_lower_q2 == rank_lower + 1
            and rank_lower_mu14 == rank_lower + 1
            and rank_all == rank_lower + 2
        ),
        "interpretation": (
            "If the two classes are independent, there is no identity "
            "mu_14=c*q_2^7 modulo the degree-14 part of the ideal "
            "(mu_1,...,mu_12) inside the generated algebra "
            "Q[mu_1,...,mu_12,q_2]. This does not decide equality of "
            "their zero divisors in the full invariant quotient."
        ),
    }


def run_singular_probe(
    singular: str,
    quadratic_mode: str,
    prime: int,
    maximum_order: int,
    power_bound: int,
    timeout: int,
) -> dict[str, object]:
    moment_orders = tuple(range(2, maximum_order + 1))
    moments = {
        order: primitive_modular_terms(
            restricted_moment_terms_mod(order, quadratic_mode, prime),
            prime,
        )
        for order in moment_orders
    }
    serialized = {
        order: serialize_modular_polynomial(terms, prime)
        for order, terms in moments.items()
    }
    ideal = ",".join(serialized.values())
    variables = ",".join(PARAMETERS)
    forbidden = (
        NULL_QUADRATIC_FORBIDDEN
        if quadratic_mode == "null"
        else ()
    )
    power_code = ""
    for index in forbidden:
        variable = PARAMETERS[index]
        suffix = str(index)
        power_code += f"""
poly h{suffix}={variable};
int exponent{suffix}=1;
while ((exponent{suffix}<={power_bound}) && (reduce(h{suffix},G)!=0))
{{
  h{suffix}=h{suffix}*{variable};
  exponent{suffix}=exponent{suffix}+1;
}}
if (reduce(h{suffix},G)==0)
{{
  print("POWER {variable} "+string(exponent{suffix}));
}}
else
{{
  print("POWER {variable} 0");
}}
"""

    if quadratic_mode == "null":
        synchronized = ",".join(PARAMETERS[index] for index in forbidden)
        containment_code = f"""
ideal J={synchronized};
ideal GJ=std(J);
print("SYNC_REMAINDER_SIZE "+string(size(reduce(I,GJ))));
"""
    else:
        containment_code = ""

    code = f"""
option(redSB);
option(prot);
ring R={prime},({variables}),dp;
ideal I={ideal};
ideal G=std(I);
print("META "+string(dim(G))+" "+string(mult(G))+" "+string(size(G)));
{containment_code}
{power_code}
exit;
"""
    with tempfile.TemporaryDirectory(
        prefix="bidegree33-casimir-"
    ) as temporary:
        path = Path(temporary) / f"{quadratic_mode}.sing"
        path.write_text(code, encoding="utf-8")
        try:
            completed = subprocess.run(
                [singular, "-q", str(path)],
                check=False,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as error:
            return {
                "quadratic_mode": quadratic_mode,
                "prime": prime,
                "moment_orders": list(moment_orders),
                "moment_term_counts": {
                    str(order): len(terms)
                    for order, terms in moments.items()
                },
                "status": "timeout",
                "timeout_seconds": timeout,
                "stdout_tail": (error.stdout or "")[-2000:],
                "stderr_tail": (error.stderr or "")[-2000:],
            }
    meta = re.search(r"(?m)^META (\d+) (\d+) (\d+)$", completed.stdout)
    powers = {
        match.group(1): int(match.group(2))
        for match in re.finditer(
            r"(?m)^POWER (\w+) (\d+)$", completed.stdout
        )
    }
    synchronization = re.search(
        r"(?m)^SYNC_REMAINDER_SIZE (\d+)$", completed.stdout
    )
    return {
        "quadratic_mode": quadratic_mode,
        "prime": prime,
        "moment_orders": list(moment_orders),
        "moment_term_counts": {
            str(order): len(terms)
            for order, terms in moments.items()
        },
        "status": (
            "completed" if completed.returncode == 0 and meta else "failed"
        ),
        "returncode": completed.returncode,
        "dimension": int(meta.group(1)) if meta else None,
        "multiplicity": int(meta.group(2)) if meta else None,
        "groebner_basis_size": int(meta.group(3)) if meta else None,
        "synchronized_linear_ideal_contains_moments": (
            synchronization is not None
            and int(synchronization.group(1)) == 0
        ),
        "forbidden_coordinate_power_memberships": powers,
        "power_search_bound": power_bound,
        "stdout_tail": completed.stdout[-4000:],
        "stderr_tail": completed.stderr[-2000:],
        "interpretation": (
            "finite-field evidence only; a completed radical equality on "
            "the null chart does not by itself reconstruct an identity "
            "over QQ"
        ),
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=DEFAULT_PRIME)
    parser.add_argument("--maximum-order", type=int, default=12)
    parser.add_argument("--power-bound", type=int, default=40)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument(
        "--run-singular-probes",
        action="store_true",
        help=(
            "also run the older direct standard-basis probes on F2=X^2 "
            "and F2=0; off by default"
        ),
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    payload: dict[str, object] = {
        "status": (
            "exact parameter comparison plus bounded finite-field "
            "zero-fiber probes"
        ),
        "degree_comparison": degree_comparison(arguments.prime),
    }
    if arguments.run_singular_probes:
        singular = shutil.which("Singular")
        if singular is None:
            raise RuntimeError(
                "Singular is required with --run-singular-probes"
            )
        payload["zero_fiber_probes"] = [
            run_singular_probe(
                singular,
                mode,
                arguments.prime,
                arguments.maximum_order,
                arguments.power_bound,
                arguments.timeout,
            )
            for mode in ("null", "zero")
        ]
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
