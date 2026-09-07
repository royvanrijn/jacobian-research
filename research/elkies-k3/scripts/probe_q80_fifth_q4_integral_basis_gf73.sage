#!/usr/bin/env sage
"""Probe the cubic integral bases behind the fifth q=4 adjoint pencil.

This is not a Riemann--Roch or Jacobian conversion.  It treats the compact
fourth-child curve as a cubic extension of GF(73)(t)(v) and asks only for the
finite and infinite maximal-order bases.  Those bases encode the conductor
corrections missing from the naive line projection.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import FunctionField, GF, PolynomialRing, lcm, matrix


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT / "artifacts/generated-results/"
    "q80-fourth-q12-cm24-moving-cubic-gf73.json"
)
artifact_bytes = ARTIFACT.read_bytes()
assert hashlib.sha256(artifact_bytes).hexdigest() == (
    "c6560b3db2d1232866e9996fc727924090aa46293c2482885cf9f9dbf4c21c89"
)
data = json.loads(artifact_bytes)

parser = argparse.ArgumentParser()
parser.add_argument("--specialize", type=int)
parser.add_argument("--infinite", action="store_true")
parser.add_argument("--valuations", action="store_true")
parser.add_argument("--pole-space-bound", type=int)
parser.add_argument("--project-section", action="store_true")
parser.add_argument("--write-projection-sample", action="store_true")
arguments = parser.parse_args()

finite = GF(73, impl="modn")
if arguments.specialize is None:
    parameter_field = FunctionField(finite, "t")
    t = parameter_field.gen()
else:
    parameter_field = finite
    t = finite(arguments.specialize)

base = FunctionField(parameter_field, "v")
v = base.gen()
polynomial_ring = PolynomialRing(base, "x")
x = polynomial_ring.gen()
polynomial = polynomial_ring(sum(
    parameter_field(coefficient)*t**t_degree*v**v_degree*x**x_degree
    for t_degree, v_degree, x_degree, coefficient
    in data["moving_terms_T_v_x_coefficient"]
))
assert polynomial.degree() == 3
if arguments.specialize is not None:
    assert polynomial.is_irreducible()
extension = base.extension(polynomial, "a")
finite_order = extension.maximal_order()
finite_basis = tuple(finite_order.basis())
print(
    "Q80FIFTHINTEGRAL|"
    f"specialize={arguments.specialize}|finite_basis_size={len(finite_basis)}|"
    f"finite_basis={tuple(map(str, finite_basis))}|"
    "status=PASS_FINITE_MAXIMAL_ORDER",
    flush=True,
)
if arguments.infinite:
    infinite_order = extension.maximal_order_infinite()
    infinite_basis = tuple(infinite_order.basis())
    print(
        "Q80FIFTHINTEGRAL|"
        f"specialize={arguments.specialize}|"
        f"infinite_basis_size={len(infinite_basis)}|"
        f"infinite_basis={tuple(map(str, infinite_basis))}|"
        "status=PASS_INFINITE_MAXIMAL_ORDER",
        flush=True,
    )
if arguments.valuations:
    infinite_places = tuple(extension.places_infinite(None))
    value_rows = []
    for place in infinite_places:
        value_rows.append(
            (
                int(place.degree()),
                int(extension(v).valuation(place)),
                tuple(int(value.valuation(place)) for value in finite_basis),
            )
        )
    print(
        "Q80FIFTHINTEGRAL|"
        f"specialize={arguments.specialize}|infinite_valuation_rows={tuple(value_rows)}|"
        "status=PASS_INFINITY_VALUATIONS",
        flush=True,
    )
if arguments.pole_space_bound is not None:
    if not arguments.infinite:
        raise ValueError("--pole-space-bound requires --infinite")
    if arguments.specialize is None:
        raise ValueError("the bounded lattice intersection is specialization-only")

    def power_vector(value):
        coefficients = list(value.list())
        return coefficients + [base(0)]*(3-len(coefficients))

    finite_matrix = matrix(
        base, 3, 3,
        lambda row, column: power_vector(finite_basis[column])[row],
    )
    infinite_matrix = matrix(
        base, 3, 3,
        lambda row, column: power_vector(infinite_basis[column])[row],
    )
    # f has no finite poles iff f=B_f*p with p in k[v]^3.  It has poles no
    # worse than v at infinity iff f=v*B_inf*q with q integral at infinity.
    transition = (v*infinite_matrix).inverse()*finite_matrix
    candidate_columns = tuple(
        (basis_index, degree)
        for basis_index in range(3)
        for degree in range(arguments.pole_space_bound+1)
    )
    polynomial_parts = {}
    maximum_part_degree = 0
    for candidate_index, (basis_index, degree) in enumerate(candidate_columns):
        coordinates = transition.column(basis_index)*v**degree
        for row, coordinate in enumerate(coordinates):
            numerator = coordinate.numerator()
            denominator = coordinate.denominator()
            quotient, _ = numerator.quo_rem(denominator)
            polynomial_parts[candidate_index, row] = quotient
            maximum_part_degree = max(maximum_part_degree, quotient.degree())
    constraint_rows = []
    for row in range(3):
        for degree in range(1, maximum_part_degree+1):
            constraint_rows.append([
                polynomial_parts[candidate_index, row][degree]
                for candidate_index in range(len(candidate_columns))
            ])
    # Pure scalar elimination avoids a platform-specific FFLAS SIGILL seen
    # for one otherwise ordinary 27x21 specialization on Apple Silicon.
    def right_kernel_basis(rows, column_count):
        reduced = [list(map(finite, row)) for row in rows]
        pivot_columns = []
        pivot_row = 0
        for column in range(column_count):
            source = next(
                (row for row in range(pivot_row, len(reduced))
                 if reduced[row][column]),
                None,
            )
            if source is None:
                continue
            reduced[pivot_row], reduced[source] = (
                reduced[source], reduced[pivot_row]
            )
            inverse = reduced[pivot_row][column]**(-1)
            reduced[pivot_row] = [inverse*value for value in reduced[pivot_row]]
            for row in range(len(reduced)):
                if row == pivot_row or not reduced[row][column]:
                    continue
                scalar = reduced[row][column]
                reduced[row] = [
                    left-scalar*right
                    for left, right in zip(reduced[row], reduced[pivot_row])
                ]
            pivot_columns.append(column)
            pivot_row += 1
            if pivot_row == len(reduced):
                break
        free_columns = [
            column for column in range(column_count)
            if column not in pivot_columns
        ]
        basis = []
        for free in free_columns:
            vector = [finite(0)]*column_count
            vector[free] = finite(1)
            for row, pivot in enumerate(pivot_columns):
                vector[pivot] = -reduced[row][free]
            basis.append(vector)
        return basis, len(pivot_columns)

    kernel_rows, constraint_rank = right_kernel_basis(
        constraint_rows, len(candidate_columns)
    )
    pole_basis = []
    for kernel_row in kernel_rows:
        value = extension(0)
        for coefficient, (basis_index, degree) in zip(
            kernel_row, candidate_columns
        ):
            value += coefficient*v**degree*finite_basis[basis_index]
        pole_basis.append(value)
    # The raw kernel may contain polynomial representations of zero only if
    # the finite basis ceased to be a basis, which cannot happen.  Validate
    # the defining valuation inequalities directly.
    infinite_places = tuple(extension.places_infinite(None))
    assert all(
        value.valuation(place) >= extension(v).valuation(place)
        for value in pole_basis for place in infinite_places
    )
    print(
        "Q80FIFTHINTEGRAL|"
        f"specialize={arguments.specialize}|"
        f"pole_space_bound={arguments.pole_space_bound}|"
        f"constraint_shape={len(constraint_rows)}x{len(candidate_columns)}|"
        f"constraint_rank={constraint_rank}|nullity={len(kernel_rows)}|"
        f"pole_basis={tuple(map(str, pole_basis))}|"
        "status=PASS_MANUAL_ORDER_INTERSECTION",
        flush=True,
    )
    if arguments.project_section:
        if len(kernel_rows) != 3:
            raise ArithmeticError("expected the cubic pole space to have dimension 3")
        # The echelon basis is 1,v,w.  Re-embed the normalization as the
        # genuine plane cubic defined by these three functions, then project
        # from the selected section R.
        w_coordinate = pole_basis[2]
        minimal_ring = PolynomialRing(base, "wmin")
        minimal = minimal_ring(w_coordinate.minimal_polynomial())
        common_denominator = lcm([
            coefficient.denominator() for coefficient in minimal.list()
        ])
        plane_ring = PolynomialRing(finite, names=("vp", "wp", "up"))
        vp, wp, up = plane_ring.gens()
        plane_cubic = plane_ring(0)
        for degree, coefficient in enumerate(minimal.list()):
            polynomial = coefficient*common_denominator
            assert polynomial.denominator() == 1
            polynomial = polynomial.numerator()
            plane_cubic += sum(
                finite(scalar)*vp**v_degree*wp**degree
                for v_degree, scalar in enumerate(polynomial.list())
            )
        assert plane_cubic.total_degree() == 3

        parameter_value = finite(arguments.specialize)
        section_v = (34*parameter_value+30)/(parameter_value+4)
        x_coefficients = (12, 26, 13, 39, 49)
        section_x_function = sum(
            finite(coefficient)*v**degree
            for degree, coefficient in enumerate(x_coefficients)
        )

        def evaluate_rational(function, value):
            return function.numerator()(value)/function.denominator()(value)

        section_x = evaluate_rational(section_x_function, section_v)
        w_coefficients = power_vector(w_coordinate)
        # Combine in GF(73)(v) before specializing.  Individual integral-basis
        # coefficients can have removable poles at a section even when the
        # represented function is regular there; termwise evaluation would
        # incorrectly raise 0/0 at such samples.
        section_w_function = sum(
            coefficient*section_x_function**degree
            for degree, coefficient in enumerate(w_coefficients)
        )
        section_w = evaluate_rational(section_w_function, section_v)
        assert plane_cubic(vp=section_v, wp=section_w, up=0) == 0
        projected = plane_ring(
            plane_cubic(wp=section_w+up*(vp-section_v))
        )
        residual, projection_remainder = projected.quo_rem(vp-section_v)
        assert projection_remainder == 0
        assert residual.degree(vp) == 2
        print(
            "Q80FIFTHINTEGRAL|"
            f"specialize={arguments.specialize}|"
            f"plane_cubic={plane_cubic}|"
            f"section_v={section_v}|section_w={section_w}|"
            f"projection_residual={residual}|"
            "residual_degree_v=2|status=PASS_COMPENSATED_SECTION_PROJECTION",
            flush=True,
        )
        if arguments.write_projection_sample:
            def polynomial_terms(polynomial):
                return [
                    [*map(int, exponents), int(coefficient)]
                    for exponents, coefficient in sorted(polynomial.dict().items())
                ]

            all_section_x_coefficients = (
                (41, 12, 67, 72, 43),
                (24, 63, 53, 61, 41),
                (7, 47, 45, 69, 62),
                (12, 26, 13, 39, 49),
                (48, 60, 25, 12, 30),
            )
            section_v_numerators_denominators = (
                (8-33*parameter_value, parameter_value+16),
                (-27-13*parameter_value, parameter_value+36),
                (18-36*parameter_value, parameter_value-32),
                (30+34*parameter_value, parameter_value+4),
                (6+3*parameter_value, parameter_value+35),
            )
            raw_w_squared_linear = (
                plane_cubic.monomial_coefficient(vp*wp**2)*vp
                + plane_cubic.monomial_coefficient(wp**2)
            )
            depressed_translation = raw_w_squared_linear/finite(3)
            explicit_section_points = []
            for section_index, (
                coefficients, (v_numerator, v_denominator)
            ) in enumerate(zip(
                all_section_x_coefficients,
                section_v_numerators_denominators,
            )):
                if not v_denominator:
                    explicit_section_points.append({
                        "index": section_index,
                        "at_v_infinity": True,
                        "x_coefficients": list(coefficients),
                    })
                    continue
                point_v = v_numerator/v_denominator
                point_x_function = sum(
                    finite(coefficient)*v**degree
                    for degree, coefficient in enumerate(coefficients)
                )
                point_x = evaluate_rational(point_x_function, point_v)
                point_w_function = sum(
                    coefficient*point_x_function**degree
                    for degree, coefficient in enumerate(w_coefficients)
                )
                point_w_raw = evaluate_rational(point_w_function, point_v)
                point_w_depressed = point_w_raw+depressed_translation(vp=point_v)
                depressed_cubic = plane_ring(
                    plane_cubic.subs({wp: wp-depressed_translation})
                )
                assert depressed_cubic(vp=point_v, wp=point_w_depressed) == 0, (
                    section_index, point_v, point_x, point_w_depressed
                )
                explicit_section_points.append({
                    "index": section_index,
                    "at_v_infinity": False,
                    "x_coefficients": list(coefficients),
                    "v": int(point_v),
                    "x": int(point_x),
                    "w_raw": int(point_w_raw),
                    "w_depressed": int(point_w_depressed),
                })

            sample_artifact = {
                "schema": "q80-fifth-q4-compensated-projection-sample-gf73-v1",
                "prime": 73,
                "T": int(parameter_value),
                "pole_space_bound": int(arguments.pole_space_bound),
                "constraint_shape": [
                    len(constraint_rows), len(candidate_columns)
                ],
                "constraint_rank": int(constraint_rank),
                "plane_cubic_terms_v_w_u_coefficient": polynomial_terms(
                    plane_cubic
                ),
                "section_v": int(section_v),
                "section_w": int(section_w),
                "explicit_section_points": explicit_section_points,
                "projection_residual_terms_v_w_u_coefficient": polynomial_terms(
                    residual
                ),
                "residual_degree_v": int(residual.degree(vp)),
                "source_artifact": str(ARTIFACT.relative_to(ROOT)),
                "source_sha256": hashlib.sha256(artifact_bytes).hexdigest(),
                "reproduce": (
                    "sage elkies-k3/scripts/"
                    "probe_q80_fifth_q4_integral_basis_gf73.sage "
                    f"--specialize {int(parameter_value)} --infinite "
                    "--pole-space-bound 6 --project-section "
                    "--write-projection-sample"
                ),
            }
            sample_path = (
                ROOT / "artifacts/generated-results/"
                f"q80-fifth-q4-compensated-projection-gf73-T{int(parameter_value)}.json"
            )
            encoded = json.dumps(
                sample_artifact, indent=2, sort_keys=True, default=int
            )+"\n"
            sample_path.write_text(encoded)
            sample_hash = hashlib.sha256(encoded.encode()).hexdigest()
            print(
                "Q80FIFTHINTEGRAL|"
                f"sample_artifact={sample_path}|sha256={sample_hash}|"
                "status=PASS_PROJECTION_SAMPLE_WRITE",
                flush=True,
            )
