#!/usr/bin/env sage-python
"""Probe the NS0007 A1+A3+2A6 fibre ansatz over a finite field.

The selected MW1 source has semistable profile ``I2+I4+2I7+4I1`` and an
exact pole-zero generator.  Supports are normalized to ``0,1,lambda,infinity``.
For every lambda and normalized degree-eight A polynomial, the script solves
the twenty branch jets for the thirteen B coefficients and retains exact
squarefree residual discriminants.

This gate imposes the reducible fibres only.  It does not yet impose the
polynomial MW section, identify the NS0007 marking, or lift to characteristic
zero.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

from sage.all import GF, PolynomialRing, binomial, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_POLES = (
    ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-rank1-section-poles-v1.json"
)
CANDIDATES = {
    "ns0007": {
        "ns_id": "NS0007",
        "source_id": "NS0007-S025",
        "source_file": "elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-c-v1.json",
        "root_type": "A1+A3+2A6",
        "height": "11/4",
        "orders": (2, 4, 7, 7),
        "profile": "I2+I4+2I7+4I1",
        "corrections": ["1/2", "3/4", "0", "0"],
    },
    "ns0034": {
        "ns_id": "NS0034",
        "source_id": "NS0034-S008",
        "source_file": "elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-a-v1.json",
        "root_type": "A2+A3+A4+A7",
        "height": "19/8",
        "orders": (4, 8, 3, 5),
        "profile": "I4+I8+I3+I5+4I1",
        "corrections": ["3/4", "7/8", "0", "0"],
    },
}


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def local_square_root(unit_coefficients, root0):
    field = root0.parent()
    answer = [field.zero()] * len(unit_coefficients)
    answer[0] = root0
    for degree in range(1, len(answer)):
        known = sum(answer[left] * answer[degree - left] for left in range(1, degree))
        answer[degree] = (unit_coefficients[degree] - known) / (2 * root0)
    return answer


def truncated_product(left, right, precision):
    field = left[0].parent()
    answer = [field.zero()] * precision
    for i, left_value in enumerate(left[:precision]):
        for j, right_value in enumerate(right[: precision - i]):
            answer[i + j] += left_value * right_value
    return answer


def multiplicative_branch(a_series, sign=1):
    field = a_series[0].parent()
    unit = [-value / field(3) for value in a_series]
    if not unit[0] or not unit[0].is_square():
        return None
    root0 = unit[0].sqrt()
    if sign == -1:
        root0 = -root0
    root = local_square_root(unit, root0)
    return [
        2 * value
        for value in truncated_product(
            truncated_product(root, root, len(root)), root, len(root)
        )
    ]


def taylor_coefficients(coefficients, point, precision):
    field = point.parent()
    return [
        sum(
            coefficients[index] * field(binomial(index, jet)) * point ** (index - jet)
            for index in range(jet, len(coefficients))
        )
        for jet in range(precision)
    ]


def order_at(poly, point):
    if not poly:
        return None
    shifted = poly(poly.parent().gen() + point)
    return min(index for index, value in enumerate(shifted.list()) if value)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--candidate", choices=sorted(CANDIDATES), default="ns0007")
parser.add_argument("--source", type=Path)
parser.add_argument("--source-id")
parser.add_argument("--section-poles", type=Path, default=DEFAULT_POLES)
parser.add_argument("--output", type=Path)
parser.add_argument("--prime", type=int, default=5)
parser.add_argument("--examples", type=int, default=20)
parser.add_argument(
    "--fixed-lambda-value",
    type=int,
    help="scan only this lambda value (which must lie outside {0,1})",
)
parser.add_argument(
    "--fixed-A8",
    type=int,
    help="fix the nonzero leading A coefficient and scan the other seven",
)
parser.add_argument(
    "--support-configuration",
    choices=("split", "conjugate-i7"),
    default="split",
    help=(
        "split uses supports 0,1,lambda,infinity; conjugate-i7 uses the "
        "rational I2/I4 supports 0,infinity and one quadratic Frobenius orbit"
    ),
)
parser.add_argument(
    "--max-a-samples-per-lambda",
    type=int,
    default=0,
    help="truncate each lambda slice; zero exhausts all normalized A polynomials",
)
parser.add_argument(
    "--sample-stride",
    type=int,
    default=0,
    help="optional coprime affine stride through coefficient indices",
)
parser.add_argument("--sample-offset", type=int, default=0)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

candidate = CANDIDATES[args.candidate]
source_path = (
    args.source.resolve()
    if args.source
    else ROOT / "artifacts/generated-results" / candidate["source_file"]
)
pole_path = args.section_poles.resolve()
source_id = args.source_id or candidate["source_id"]
output_path = (
    args.output.resolve()
    if args.output
    else ROOT
    / "artifacts/generated-results"
    / (
        f"elkies-k3-lattice-foundry-{args.candidate}-source-ansatz-mod{args.prime}.json"
        if args.support_configuration == "split"
        else f"elkies-k3-lattice-foundry-{args.candidate}-source-ansatz-"
        f"conjugate-i7-mod{args.prime}.json"
    )
)
source_payload = json.loads(source_path.read_text())
source_entry = next(
    row
    for row in source_payload["sources"]
    if row["ns_id"] == candidate["ns_id"] and row["source_id"] == source_id
)
source = source_entry["source"]
assert source["root_type"] == candidate["root_type"]
assert source["root_lattice_primitive"] and source["torsion"] == 1
assert source["mw_height_gram"] == [[candidate["height"]]]
pole_payload = json.loads(pole_path.read_text())
pole_row = next(
    row
    for row in pole_payload["sources"]
    if row["source_artifact"] == relative(source_path)
    and row["source_id"] == source_id
)
assert pole_row["minimum_section_pole_order"] == 0

field = GF(args.prime)
if field.characteristic() in (2, 3):
    raise SystemExit("--prime must differ from 2 and 3")
fixed_lambda = (
    None if args.fixed_lambda_value is None else field(args.fixed_lambda_value)
)
if fixed_lambda is not None and fixed_lambda in (field.zero(), field.one()):
    raise SystemExit("--fixed-lambda-value must lie outside {0,1}")
fixed_a8 = None if args.fixed_A8 is None else field(args.fixed_A8)
if fixed_a8 is not None and not fixed_a8:
    raise SystemExit("--fixed-A8 must be nonzero")
ring = PolynomialRing(field, "t")
t = ring.gen()

examples = []
accounting = {
    "normalized_A_samples": 0,
    "branch_eligible_with_signs": 0,
    "hermite_compatible_with_signs": 0,
    "exact_prescribed_orders": 0,
    "squarefree_examples_with_signs": 0,
}
lambda_records = []
per_lambda_total = args.prime ** (7 if fixed_a8 is not None else 8)
if args.sample_stride and math.gcd(args.sample_stride, per_lambda_total) != 1:
    raise SystemExit("--sample-stride must be coprime to prime^8")


def coefficient_digits(index):
    digit_count = 7 if fixed_a8 is not None else 8
    digits = [0] * digit_count
    for position in range(digit_count - 1, -1, -1):
        digits[position] = index % args.prime
        index //= args.prime
    return tuple(digits)


def polynomial_factor_order(poly, factor):
    order = 0
    while poly and not poly.quo_rem(factor)[1]:
        poly = poly.quo_rem(factor)[0]
        order += 1
    return order


def conjugate_support_orbits():
    """Represent quadratic I7 support orbits modulo t -> c*t."""
    records = []
    for norm in field:
        support = t**2 - t + norm
        if support.is_irreducible():
            records.append(("trace-one", int(norm), support))
    for norm in field:
        support = t**2 + norm
        if support.is_irreducible():
            records.append(("trace-zero", int(norm), support))
            break
    return records


if args.support_configuration == "conjugate-i7":
    if args.candidate != "ns0007":
        raise SystemExit("conjugate-i7 is currently implemented only for NS0007")
    if fixed_lambda is not None:
        raise SystemExit("--fixed-lambda-value applies only to split supports")
    support_records = []
    order_zero, order_infinity, order_quadratic = 2, 4, 7
    for orbit_type, orbit_parameter, support in conjugate_support_orbits():
        prime = args.prime
        trace = int(-support[1]) % prime
        norm = int(support[0]) % prime

        def pair_add(left, right):
            return ((left[0] + right[0]) % prime, (left[1] + right[1]) % prime)

        def pair_neg(value):
            return ((-value[0]) % prime, (-value[1]) % prime)

        def pair_mul(left, right):
            return (
                (left[0] * right[0] - left[1] * right[1] * norm) % prime,
                (left[0] * right[1]
                + left[1] * right[0]
                + left[1] * right[1] * trace) % prime,
            )

        def pair_scale(scalar, value):
            return ((scalar * value[0]) % prime, (scalar * value[1]) % prime)

        def pair_inverse(value):
            denominator = (
                (value[0] ** 2
                + value[0] * value[1] * trace
                + value[1] ** 2 * norm) % prime
            )
            if not denominator:
                raise ZeroDivisionError("zero in the quadratic residue field")
            inverse_denominator = pow(denominator, -1, prime)
            return (
                (value[0] + value[1] * trace) * inverse_denominator % prime,
                -value[1] * inverse_denominator % prime,
            )

        def pair_power(value, exponent):
            answer = (1, 0)
            while exponent:
                if exponent & 1:
                    answer = pair_mul(answer, value)
                value = pair_mul(value, value)
                exponent //= 2
            return answer

        def pair_branch(a_series):
            inverse_three = pow(3, -1, prime)
            unit = [pair_scale(-inverse_three, value) for value in a_series]
            if unit[0] == (0, 0):
                return None
            root0 = extension_square_roots.get(unit[0])
            if root0 is None:
                return None
            answer = [(0, 0)] * len(unit)
            answer[0] = root0
            inverse_twice_root = pair_inverse(pair_scale(2, root0))
            for degree in range(1, len(answer)):
                known = (0, 0)
                for left in range(1, degree):
                    known = pair_add(
                        known, pair_mul(answer[left], answer[degree - left])
                    )
                answer[degree] = pair_mul(
                    pair_add(unit[degree], pair_neg(known)), inverse_twice_root
                )
            cube = []
            for degree in range(len(answer)):
                cube_value = (0, 0)
                for left in range(degree + 1):
                    for middle in range(degree - left + 1):
                        right = degree - left - middle
                        cube_value = pair_add(
                            cube_value,
                            pair_mul(
                                pair_mul(answer[left], answer[middle]), answer[right]
                            ),
                        )
                cube.append(pair_scale(2, cube_value))
            return cube

        def base_branch(a_series):
            inverse_three = pow(3, -1, prime)
            unit = [(-value * inverse_three) % prime for value in a_series]
            if not unit[0] or unit[0] not in base_square_roots:
                return None
            answer = [0] * len(unit)
            answer[0] = base_square_roots[unit[0]]
            inverse_twice_root = pow(2 * answer[0], -1, prime)
            for degree in range(1, len(answer)):
                known = sum(
                    answer[left] * answer[degree - left]
                    for left in range(1, degree)
                )
                answer[degree] = (unit[degree] - known) * inverse_twice_root % prime
            return [
                2
                * sum(
                    answer[left] * answer[middle] * answer[degree - left - middle]
                    for left in range(degree + 1)
                    for middle in range(degree - left + 1)
                )
                % prime
                for degree in range(len(answer))
            ]

        alpha = (0, 1)
        extension_elements = [
            (left, right) for left in range(prime) for right in range(prime)
        ]
        base_square_roots = {
            value * value % prime: value for value in range(prime)
        }
        extension_square_roots = {
            pair_mul(value, value): value for value in extension_elements
        }

        rows = []
        for jet in range(order_zero):
            rows.append(
                [
                    field(binomial(index, jet)) * field.zero() ** (index - jet)
                    if index >= jet
                    else field.zero()
                    for index in range(13)
                ]
            )
        for jet in range(order_quadratic):
            extension_row = [
                pair_scale(
                    int(binomial(index, jet)) % prime, pair_power(alpha, index - jet)
                )
                if index >= jet
                else (0, 0)
                for index in range(13)
            ]
            for coordinate in range(2):
                rows.append([value[coordinate] for value in extension_row])
        for jet in range(order_infinity):
            rows.append([field(index == 12 - jet) for index in range(13)])
        hermite = matrix(field, rows)
        if hermite.nrows() != 20 or hermite.ncols() != hermite.rank():
            raise ArithmeticError("unexpected conjugate-support Hermite rank")
        compatibility = hermite.left_kernel().basis_matrix()
        if compatibility.nrows() != 7:
            raise ArithmeticError("unexpected conjugate-support compatibility codimension")

        local = {key: 0 for key in accounting}
        samples = 0
        sample_limit = args.max_a_samples_per_lambda or per_lambda_total
        if args.sample_stride:
            coefficient_iterator = (
                (
                    (args.sample_offset + sample_index * args.sample_stride)
                    % per_lambda_total,
                    coefficient_digits(
                        (args.sample_offset + sample_index * args.sample_stride)
                        % per_lambda_total
                    ),
                )
                for sample_index in range(sample_limit)
            )
        else:
            coefficient_iterator = itertools.islice(
                enumerate(
                    itertools.product(
                        range(args.prime), repeat=7 if fixed_a8 is not None else 8
                    )
                ),
                sample_limit,
            )
        for coefficient_index, digits in coefficient_iterator:
            samples += 1
            local["normalized_A_samples"] += 1
            coefficients = [(-3) % prime] + list(digits)
            if fixed_a8 is not None:
                coefficients.append(int(fixed_a8))
            if not coefficients[8]:
                continue
            A = ring([field(value) for value in coefficients])
            zero_series = coefficients[:order_zero]
            alpha_series = []
            for jet in range(order_quadratic):
                value = (0, 0)
                for index in range(jet, len(coefficients)):
                    value = pair_add(
                        value,
                        pair_scale(
                            coefficients[index] * int(binomial(index, jet)) % prime,
                            pair_power(alpha, index - jet),
                        ),
                    )
                alpha_series.append(value)
            infinity_series = [
                coefficients[8 - jet] for jet in range(order_infinity)
            ]
            zero_branch = base_branch(zero_series)
            alpha_branch = pair_branch(alpha_series)
            infinity_branch = base_branch(infinity_series)
            if any(
                branch is None
                for branch in (zero_branch, alpha_branch, infinity_branch)
            ):
                continue
            for alpha_sign, infinity_sign in itertools.product((1, -1), repeat=2):
                signed_alpha = [
                    pair_scale(alpha_sign, value) for value in alpha_branch
                ]
                alpha_coordinates = []
                for value in signed_alpha:
                    alpha_coordinates.extend(value)
                target = vector(
                    field,
                    list(zero_branch)
                    + alpha_coordinates
                    + [infinity_sign * value % prime for value in infinity_branch],
                )
                local["branch_eligible_with_signs"] += 1
                if compatibility * target:
                    continue
                b_coefficients = list(hermite.solve_right(target))
                local["hermite_compatible_with_signs"] += 1
                B = ring(b_coefficients)
                discriminant_core = 4 * A**3 + 27 * B**2
                orders = (
                    order_at(discriminant_core, field.zero()),
                    polynomial_factor_order(discriminant_core, support),
                    24 - discriminant_core.degree(),
                )
                if orders != (order_zero, order_quadratic, order_infinity):
                    continue
                local["exact_prescribed_orders"] += 1
                divisor = t**order_zero * support**order_quadratic
                residual, remainder = discriminant_core.quo_rem(divisor)
                if remainder or residual.degree() != 4:
                    raise ArithmeticError("unexpected conjugate-support residual")
                if not residual(field.zero()) or residual.gcd(support).degree():
                    continue
                if residual.gcd(residual.derivative()).degree() != 0:
                    continue
                local["squarefree_examples_with_signs"] += 1
                if len(examples) < args.examples:
                    # Recover the node root directly; the branch constant is 2*r^3.
                    alpha_root = extension_square_roots[
                        pair_scale(-pow(3, -1, prime), alpha_series[0])
                    ]
                    if alpha_sign == -1:
                        alpha_root = pair_neg(alpha_root)
                    infinity_root = base_square_roots[
                        -infinity_series[0] * pow(3, -1, prime) % prime
                    ]
                    if infinity_sign == -1:
                        infinity_root = -infinity_root % prime
                    examples.append(
                        {
                            "support_orbit_type": orbit_type,
                            "support_orbit_parameter": orbit_parameter,
                            "quadratic_I7_support": str(support),
                            "coefficient_index_within_support": coefficient_index,
                            "branch_signs_at_quadratic_infinity": [
                                alpha_sign,
                                infinity_sign,
                            ],
                            "A_coefficients_low_to_high": [
                                int(value) for value in coefficients
                            ],
                            "B_coefficients_low_to_high": [
                                int(value) for value in b_coefficients
                            ],
                            "residual_discriminant_coefficients_low_to_high": [
                                int(value) for value in residual
                            ],
                            "multiplicative_tangent_characters": {
                                "I2_at_zero": (
                                    "split" if 3 % prime in base_square_roots else "nonsplit"
                                ),
                                "I4_at_infinity": (
                                    "split"
                                    if 3 * infinity_root % prime in base_square_roots
                                    else "nonsplit"
                                ),
                                "conjugate_I7_pair": (
                                    "split"
                                    if pair_scale(3, alpha_root)
                                    in extension_square_roots
                                    else "nonsplit"
                                ),
                            },
                            "geometric_fibre_profile": candidate["profile"],
                        }
                    )
        for key in accounting:
            accounting[key] += local[key]
        support_records.append(
            {
                "orbit_type": orbit_type,
                "orbit_parameter": orbit_parameter,
                "quadratic_support": str(support),
                "samples_consumed": samples,
                "exhausted": samples == per_lambda_total
                and not args.max_a_samples_per_lambda,
                "accounting": local,
                "sample_stride": args.sample_stride or 1,
                "sample_offset": args.sample_offset,
            }
        )

    exhausted = all(record["exhausted"] for record in support_records)
    output = {
        "schema": "elkies-k3.lattice-foundry-ns0007-source-ansatz-arithmetic-modp.v1",
        "status": (
            "PASS_EXACT_EXHAUSTIVE_ARITHMETIC_SOURCE_FIBRE_ANSATZ_WITH_EXAMPLES"
            if exhausted and examples
            else "PASS_EXACT_EXHAUSTIVE_ARITHMETIC_SOURCE_FIBRE_ANSATZ_EMPTY"
            if exhausted
            else "PASS_BOUNDED_ARITHMETIC_SOURCE_FIBRE_ANSATZ"
            if examples
            else "PASS_BOUNDED_NO_ARITHMETIC_SOURCE_FIBRE_ANSATZ"
        ),
        "prime": args.prime,
        "support_configuration": "conjugate-i7",
        "scan": {
            "quadratic_support_orbits": support_records,
            "exhausted": exhausted,
            "normalized_A_polynomials_per_orbit": per_lambda_total,
        },
        "accounting": accounting | {"stored_examples": len(examples)},
        "ansatz": {
            "short_weierstrass": "y^2=x^3+A(t)x+B(t)",
            "degree_bounds": {"A": 8, "B": 12},
            "normalization": (
                "A(0)=-3; rational I2/I4 supports at 0,infinity; the two I7 "
                "supports form one irreducible quadratic Frobenius orbit"
            ),
            "normalized_reducible_supports": [
                "0:I2",
                "irreducible quadratic:I7+I7",
                "infinity:I4",
            ],
            "hermite_conditions_over_GF_p": 20,
            "B_coefficient_rank": 13,
            "compatibility_equations_on_A": 7,
            "expected_fibre_stratum_dimension": 2,
            "minimum_section_pole_order_to_impose": 0,
            "section_component_corrections": candidate["corrections"],
        },
        "examples": examples,
        "source": {
            "artifact": relative(source_path),
            "artifact_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
            "source_id": source_id,
            "source_gram_sha256": source["gram_sha256"],
            "section_pole_artifact": relative(pole_path),
            "section_pole_artifact_sha256": hashlib.sha256(
                pole_path.read_bytes()
            ).hexdigest(),
        },
        "proof_boundary": {
            "proved": (
                "Every stored model has the exact geometric fibre profile over "
                "the algebraic closure, with the two I7 supports forming one "
                "quadratic Frobenius orbit. Tangent-cone characters are audited."
                if examples
                else "The scan exactly checks the declared bounded arithmetic charts."
            ),
            "not_proved": (
                "The pole-zero MW section, NS0007 marking, characteristic-zero "
                "lift, and neighbour route are not proved."
            ),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/probe_lattice_foundry_ns0007_source_ansatz_modp.sage "
            f"--candidate ns0007 --prime {args.prime} "
            "--support-configuration conjugate-i7"
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    if args.check:
        if output_path.read_text() != serialized:
            raise SystemExit("NS0007 arithmetic modular source-ansatz artifact is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "FOUNDRYNS0007ARITHMETICANSATZ|"
        f"p={args.prime}|orbits={len(support_records)}|"
        f"samples={accounting['normalized_A_samples']}|"
        f"compatible={accounting['hermite_compatible_with_signs']}|"
        f"squarefree={accounting['squarefree_examples_with_signs']}|"
        f"exhausted={int(exhausted)}|"
        f"status={'PASS' if exhausted or examples else 'BOUNDED_NEGATIVE'}",
        flush=True,
    )
    raise SystemExit(0)


lambda_values = [fixed_lambda] if fixed_lambda is not None else list(field)
for lambda_value in lambda_values:
    if lambda_value in (field.zero(), field.one()):
        continue
    rows = []
    order_zero, order_one, order_lambda, order_infinity = candidate["orders"]
    for point, precision in (
        (field.zero(), order_zero),
        (field.one(), order_one),
        (lambda_value, order_lambda),
    ):
        for jet in range(precision):
            rows.append(
                [
                    field(binomial(index, jet)) * point ** (index - jet)
                    if index >= jet
                    else field.zero()
                    for index in range(13)
                ]
            )
    for jet in range(order_infinity):
        rows.append([field(index == 12 - jet) for index in range(13)])
    hermite = matrix(field, rows)
    if hermite.nrows() != 20 or hermite.ncols() != hermite.rank():
        raise ArithmeticError("unexpected Hermite rank")
    compatibility = hermite.left_kernel().basis_matrix()
    if compatibility.nrows() != 7:
        raise ArithmeticError("unexpected compatibility codimension")

    local = {key: 0 for key in accounting}
    samples = 0
    sample_limit = args.max_a_samples_per_lambda or per_lambda_total
    if args.sample_stride:
        coefficient_iterator = (
            (
                (args.sample_offset + sample_index * args.sample_stride)
                % per_lambda_total,
                coefficient_digits(
                    (args.sample_offset + sample_index * args.sample_stride)
                    % per_lambda_total
                ),
            )
            for sample_index in range(sample_limit)
        )
    else:
        coefficient_iterator = itertools.islice(
            enumerate(
                itertools.product(
                    range(args.prime), repeat=7 if fixed_a8 is not None else 8
                )
            ),
            sample_limit,
        )
    for coefficient_index, digits in coefficient_iterator:
        samples += 1
        local["normalized_A_samples"] += 1
        coefficients = [field(-3)] + [field(value) for value in digits]
        if fixed_a8 is not None:
            coefficients.append(fixed_a8)
        if not coefficients[8]:
            continue
        A = ring(coefficients)
        series = (
            coefficients[:order_zero],
            taylor_coefficients(coefficients, field.one(), order_one),
            taylor_coefficients(coefficients, lambda_value, order_lambda),
            [coefficients[8 - jet] for jet in range(order_infinity)],
        )
        positive = tuple(multiplicative_branch(values, 1) for values in series)
        if any(branch is None for branch in positive):
            continue
        for signs in itertools.product((1, -1), repeat=3):
            branches = (
                positive[0],
                [signs[0] * value for value in positive[1]],
                [signs[1] * value for value in positive[2]],
                [signs[2] * value for value in positive[3]],
            )
            local["branch_eligible_with_signs"] += 1
            target = vector(field, sum((list(branch) for branch in branches), []))
            if compatibility * target:
                continue
            b_coefficients = list(hermite.solve_right(target))
            local["hermite_compatible_with_signs"] += 1
            B = ring(b_coefficients)
            discriminant_core = 4 * A**3 + 27 * B**2
            orders = (
                order_at(discriminant_core, field.zero()),
                order_at(discriminant_core, field.one()),
                order_at(discriminant_core, lambda_value),
                24 - discriminant_core.degree(),
            )
            if orders != candidate["orders"]:
                continue
            local["exact_prescribed_orders"] += 1
            divisor = (
                t**order_zero
                * (t - 1) ** order_one
                * (t - lambda_value) ** order_lambda
            )
            residual, remainder = discriminant_core.quo_rem(divisor)
            if remainder or residual.degree() != 4:
                raise ArithmeticError("unexpected residual discriminant")
            if any(residual(point) == 0 for point in (0, 1, lambda_value)):
                continue
            if residual.gcd(residual.derivative()).degree() != 0:
                continue
            local["squarefree_examples_with_signs"] += 1
            if len(examples) < args.examples:
                examples.append(
                    {
                        "lambda": int(lambda_value),
                        "sample_index_within_lambda": samples,
                        "coefficient_index_within_lambda": coefficient_index,
                        "branch_signs_at_one_lambda_infinity": list(signs),
                        "A_coefficients_low_to_high": [int(value) for value in coefficients],
                        "B_coefficients_low_to_high": [int(value) for value in b_coefficients],
                        "residual_discriminant_coefficients_low_to_high": [
                            int(value) for value in residual
                        ],
                        "geometric_fibre_profile": candidate["profile"],
                    }
                )
    for key in accounting:
        accounting[key] += local[key]
    lambda_records.append(
        {
            "lambda": int(lambda_value),
            "samples_consumed": samples,
            "exhausted": samples == per_lambda_total and not args.max_a_samples_per_lambda,
            "accounting": local,
            "sample_stride": args.sample_stride or 1,
            "sample_offset": args.sample_offset,
        }
    )

exhausted = all(record["exhausted"] for record in lambda_records)
output = {
    "schema": f"elkies-k3.lattice-foundry-{args.candidate}-source-ansatz-modp.v1",
    "status": (
        "PASS_EXACT_EXHAUSTIVE_MODULAR_SOURCE_FIBRE_ANSATZ_WITH_EXAMPLES"
        if exhausted and examples
        else "PASS_EXACT_EXHAUSTIVE_MODULAR_SOURCE_FIBRE_ANSATZ_EMPTY"
        if exhausted
        else "PASS_BOUNDED_MODULAR_SOURCE_FIBRE_ANSATZ"
        if examples
        else "PASS_BOUNDED_NO_MODULAR_SOURCE_FIBRE_ANSATZ"
    ),
    "prime": args.prime,
    "scan": {
        "lambda_slices": lambda_records,
        "exhausted": exhausted,
        "normalized_A_polynomials_per_lambda": per_lambda_total,
    },
    "accounting": accounting | {"stored_examples": len(examples)},
    "ansatz": {
        "short_weierstrass": "y^2=x^3+A(t)x+B(t)",
        "degree_bounds": {"A": 8, "B": 12},
        "normalization": "A(0)=-3; supports at 0,1,lambda,infinity",
        "normalized_reducible_supports": [
            f"0:I{order_zero}",
            f"1:I{order_one}",
            f"lambda:I{order_lambda}",
            f"infinity:I{order_infinity}",
        ],
        "hermite_conditions": 20,
        "B_coefficient_rank": 13,
        "compatibility_equations_on_A": 7,
        "expected_fibre_stratum_dimension": 2,
        "minimum_section_pole_order_to_impose": 0,
        "section_component_corrections": candidate["corrections"],
    },
    "examples": examples,
    "source": {
        "artifact": relative(source_path),
        "artifact_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "source_id": source_id,
        "source_gram_sha256": source["gram_sha256"],
        "section_pole_artifact": relative(pole_path),
        "section_pole_artifact_sha256": hashlib.sha256(pole_path.read_bytes()).hexdigest(),
    },
    "proof_boundary": {
        "proved": (
            "The scan exhausts the displayed normalized finite-field fibre chart "
            "and proves that it contains no squarefree model."
            if exhausted and not examples
            else "Every stored example has the exact displayed finite-field fibre profile."
            if examples
            else "The scan exactly checks the declared bounded coefficient slices."
        ),
        "not_proved": (
            f"The pole-zero MW section, {candidate['ns_id']} marking, rational parameterization, "
            "characteristic-zero lift, and neighbour route are not proved."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/probe_lattice_foundry_ns0007_source_ansatz_modp.sage "
        f"--candidate {args.candidate}"
        + (
            f" --fixed-lambda-value {int(fixed_lambda)}"
            if fixed_lambda is not None
            else ""
        )
        + (f" --fixed-A8 {int(fixed_a8)}" if fixed_a8 is not None else "")
    ),
}
if fixed_lambda is not None:
    output["ansatz"]["fixed_lambda"] = int(fixed_lambda)
if fixed_a8 is not None:
    output["ansatz"]["fixed_A8"] = int(fixed_a8)
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if output_path.read_text() != serialized:
        raise SystemExit(f"{candidate['ns_id']} modular source-ansatz artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    f"FOUNDRY{candidate['ns_id']}ANSATZ|"
    f"p={args.prime}|lambdas={len(lambda_records)}|"
    f"samples={accounting['normalized_A_samples']}|"
    f"compatible={accounting['hermite_compatible_with_signs']}|"
    f"squarefree={accounting['squarefree_examples_with_signs']}|"
    f"exhausted={int(exhausted)}|"
    f"status={'PASS' if exhausted or examples else 'BOUNDED_NEGATIVE'}",
    flush=True,
)
