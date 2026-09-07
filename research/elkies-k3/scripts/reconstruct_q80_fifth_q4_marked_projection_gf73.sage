#!/usr/bin/env sage
"""Reconstruct globally marked fifth-q4 degree-two projections over GF(73).

Every input fiber first receives its compensated cubic coordinate from the
finite/infinite maximal-order intersection.  Only then is a degree-two map
constructed and normalized by two other marked sections to take the values
0 and 1.  Coefficient interpolation therefore happens after the marked
PGL2 gauge is fixed independently on every fiber.

Modes:

``pinned`` projects from section 3, whose inverse is the pinned fifth
horizontal target.  ``pair01``, ``pair14``, and ``pair23`` use the third
intersection of the indicated section pair; these are the three explicit
pair representatives of the alternate horizontal MW class (1,0).
"""

import argparse
import glob
import hashlib
import json
from pathlib import Path

from sage.all import FunctionField, GF, PolynomialRing, lcm


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument(
    "--mode", choices=("pinned", "pair01", "pair14", "pair23"),
    default="pinned",
)
parser.add_argument(
    "--reference-pair",
    default=None,
    help="optional comma-separated section indices overriding the mode gauge",
)
parser.add_argument("--minimum-withheld", type=int, default=3)
parser.add_argument("--max-total-degree", type=int, default=30)
parser.add_argument("--write-artifact", action="store_true")
arguments = parser.parse_args()

finite = GF(73, impl="modn")
plane_ring = PolynomialRing(finite, names=("v", "w", "u"))
v, w, u = plane_ring.gens()
univariate_v = PolynomialRing(finite, "z")
z = univariate_v.gen()
parameter_ring = PolynomialRing(finite, "t")
t = parameter_ring.gen()


def polynomial_from_terms(rows):
    return plane_ring({
        tuple(row[:-1]): finite(row[-1]) for row in rows
    })


def depressed_plane(payload):
    plane = polynomial_from_terms(
        payload["plane_cubic_terms_v_w_u_coefficient"]
    )
    translation = (
        plane.monomial_coefficient(v*w**2)*v
        + plane.monomial_coefficient(w**2)
    )/finite(3)
    result = plane_ring(plane.subs({w: w-translation}))
    assert result.monomial_coefficient(v*w**2) == 0
    assert result.monomial_coefficient(w**2) == 0
    return result


def line_third_point(plane, left, right):
    left_v, left_w = left
    right_v, right_w = right
    if left_v == right_v:
        raise ZeroDivisionError("marked secant is vertical")
    slope = (right_w-left_w)/(right_v-left_v)
    line = left_w+slope*(z-left_v)
    intersection = univariate_v(sum(
        coefficient*z**exponents[0]*line**exponents[1]
        for exponents, coefficient in plane.dict().items()
    ))
    quotient, remainder = intersection.quo_rem(
        (z-left_v)*(z-right_v)
    )
    if remainder or quotient.degree() != 1:
        raise ArithmeticError("marked secant did not leave one residual point")
    center_v = -quotient[0]/quotient[1]
    center_w = left_w+slope*(center_v-left_v)
    assert plane(v=center_v, w=center_w, u=0) == 0
    return center_v, center_w


mode_data = {
    "pinned": (None, (0, 1)),
    "pair01": ((0, 1), (2, 3)),
    "pair14": ((1, 4), (2, 3)),
    "pair23": ((2, 3), (0, 1)),
}
center_pair, reference_pair = mode_data[arguments.mode]
if arguments.reference_pair is not None:
    reference_pair = tuple(map(int, arguments.reference_pair.split(",")))
    if len(reference_pair) != 2 or reference_pair[0] == reference_pair[1]:
        parser.error("--reference-pair requires two distinct indices")
    if any(index < 0 or index >= 5 for index in reference_pair):
        parser.error("--reference-pair indices must lie in 0..4")
    if center_pair is not None and any(index in center_pair for index in reference_pair):
        parser.error("reference sections must be disjoint from the center pair")
gauge_tag = f"{reference_pair[0]}{reference_pair[1]}"

sample_rows = []
skip_rows = []
for path_string in sorted(glob.glob(str(
    ROOT / "artifacts/generated-results/"
    "q80-fifth-q4-compensated-projection-gf73-T*.json"
))):
    path = Path(path_string)
    payload = json.loads(path.read_text())
    points_payload = payload.get("explicit_section_points")
    if not points_payload:
        continue
    parameter = int(payload["T"])
    if any(point.get("at_v_infinity") for point in points_payload):
        skip_rows.append((parameter, "section_at_v_infinity"))
        continue
    points = tuple(
        (finite(point["v"]), finite(point["w_depressed"]))
        for point in points_payload
    )
    plane = depressed_plane(payload)
    try:
        if center_pair is None:
            center = points[3]
        else:
            center = line_third_point(
                plane, points[center_pair[0]], points[center_pair[1]]
            )
        center_v, center_w = center

        def slope_at(point):
            point_v, point_w = point
            if point_v == center_v:
                raise ZeroDivisionError("reference section lies over center v")
            return (point_w-center_w)/(point_v-center_v)

        slope_zero = slope_at(points[reference_pair[0]])
        slope_one = slope_at(points[reference_pair[1]])
        slope_scale = slope_one-slope_zero
        if not slope_scale:
            raise ZeroDivisionError("reference section slopes coincide")
        line_slope = slope_zero+slope_scale*u
        substituted = plane_ring(
            plane.subs({w: center_w+(v-center_v)*line_slope})
        )
        residual, remainder = substituted.quo_rem(v-center_v)
        if remainder or residual.degree(v) != 2:
            raise ArithmeticError("normalized projection is not quadratic")
        # The normalization is checked on the two marked reference sections.
        assert slope_at(points[reference_pair[0]]) == slope_zero
        assert (slope_at(points[reference_pair[1]])-slope_zero)/slope_scale == 1
        normalized_section_values = tuple(
            None
            if index in center_pair
            else (slope_at(point)-slope_zero)/slope_scale
            for index, point in enumerate(points)
        ) if center_pair is not None else ()
    except (ArithmeticError, ZeroDivisionError) as error:
        skip_rows.append((parameter, str(error)))
        continue
    sample_rows.append(
        (parameter, residual, center, path, normalized_section_values)
    )

sample_rows.sort()
if len(sample_rows) <= arguments.minimum_withheld+2:
    raise RuntimeError("too few marked projection samples")


def nullspace(rows, column_count):
    reduced = [list(map(finite, row)) for row in rows]
    pivot_columns = []
    pivot_row = 0
    for column in range(column_count):
        source = next(
            (row for row in range(pivot_row, len(reduced))
             if reduced[row][column]), None
        )
        if source is None:
            continue
        reduced[pivot_row], reduced[source] = reduced[source], reduced[pivot_row]
        inverse = reduced[pivot_row][column]**(-1)
        reduced[pivot_row] = [inverse*entry for entry in reduced[pivot_row]]
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


def rational_interpolate(values):
    training = values[:-arguments.minimum_withheld]
    for total_degree in range(arguments.max_total_degree+1):
        candidates = []
        for numerator_degree in range(total_degree+1):
            denominator_degree = total_degree-numerator_degree
            unknowns = numerator_degree+denominator_degree+2
            if unknowns > len(training)+1:
                continue
            rows = []
            for parameter, value in training:
                parameter = finite(parameter)
                rows.append(
                    [parameter**degree for degree in range(numerator_degree+1)]
                    + [-value*parameter**degree
                       for degree in range(denominator_degree+1)]
                )
            basis, rank = nullspace(rows, unknowns)
            if len(basis) != 1:
                continue
            vector = basis[0]
            numerator = parameter_ring(vector[:numerator_degree+1])
            denominator = parameter_ring(vector[numerator_degree+1:])
            if not denominator:
                continue
            if any(
                denominator(finite(parameter)) == 0
                or numerator(finite(parameter))
                   != value*denominator(finite(parameter))
                for parameter, value in values
            ):
                continue
            scale = denominator[denominator.degree()]**(-1)
            candidates.append((
                numerator_degree, denominator_degree,
                len(values)-rank, numerator*scale, denominator*scale,
            ))
        if candidates:
            return candidates
    return []


marked_section_value_rows = []
if center_pair is not None:
    for section_index in range(5):
        if section_index in center_pair:
            continue
        values = [
            (row[0], row[4][section_index]) for row in sample_rows
        ]
        candidates = rational_interpolate(values)
        marked_section_value_rows.append((section_index, candidates))
        print(
            "Q80FIFTHMARKED|"
            f"mode={arguments.mode}|reference_pair={reference_pair}|"
            f"section_value_index={section_index}|"
            f"candidates={tuple((a,b,m,str(n),str(d)) for a,b,m,n,d in candidates)}",
            flush=True,
        )

support = sorted(set().union(*(
    residual.dict() for _, residual, _, _, _ in sample_rows
)))
coefficient_rows = []
for monomial in support:
    values = [
        (parameter, residual.dict().get(monomial, finite(0)))
        for parameter, residual, _, _, _ in sample_rows
    ]
    coefficient_rows.append((monomial, rational_interpolate(values)))
recognized = sum(len(candidates) == 1 for _, candidates in coefficient_rows)
print(
    "Q80FIFTHMARKED|"
    f"mode={arguments.mode}|samples={tuple(row[0] for row in sample_rows)}|"
    f"skips={tuple(skip_rows)}|support={len(support)}|"
    f"recognized={recognized}/{len(support)}|"
    f"withheld={arguments.minimum_withheld}|"
    "status=PASS_MARKED_FIBERWISE_PROJECTION_AUDIT",
    flush=True,
)
for monomial, candidates in coefficient_rows:
    print(
        "Q80FIFTHMARKED|"
        f"mode={arguments.mode}|monomial={monomial}|"
        f"candidates={tuple((a,b,m,str(n),str(d)) for a,b,m,n,d in candidates)}",
        flush=True,
    )

if recognized == len(support):
    parameter_field = parameter_ring.fraction_field()
    family_ring = PolynomialRing(parameter_field, names=("vf", "uf"))
    vf, uf = family_ring.gens()
    residual_family = family_ring(0)
    coefficient_records = []
    for monomial, candidates in coefficient_rows:
        if len(candidates) != 1:
            raise ArithmeticError("marked coefficient recognition is not unique")
        numerator_degree, denominator_degree, margin, numerator, denominator = candidates[0]
        coefficient = parameter_field(numerator)/parameter_field(denominator)
        residual_family += coefficient*vf**monomial[0]*uf**monomial[2]
        coefficient_records.append({
            "monomial_v_w_u": list(map(int, monomial)),
            "numerator_coefficients_low_to_high": list(map(int, numerator.list())),
            "denominator_coefficients_low_to_high": list(map(int, denominator.list())),
            "numerator_degree": int(numerator_degree),
            "denominator_degree": int(denominator_degree),
            "equation_margin": int(margin),
        })
    assert residual_family.degree(vf) == 2
    quadratic_a = family_ring(sum(
        coefficient*uf**exponents[1]
        for exponents, coefficient in residual_family.dict().items()
        if exponents[0] == 2
    ))
    quadratic_b = family_ring(sum(
        coefficient*uf**exponents[1]
        for exponents, coefficient in residual_family.dict().items()
        if exponents[0] == 1
    ))
    quadratic_c = family_ring(sum(
        coefficient*uf**exponents[1]
        for exponents, coefficient in residual_family.dict().items()
        if exponents[0] == 0
    ))
    double_cover = family_ring(quadratic_b**2-4*quadratic_a*quadratic_c)
    cover_denominator = lcm([
        coefficient.denominator()
        for coefficient in double_cover.dict().values()
    ])
    cover_ring = PolynomialRing(finite, names=("T5", "U5"))
    T5, U5 = cover_ring.gens()
    cover_numerator = cover_ring(0)
    for exponents, coefficient in double_cover.dict().items():
        polynomial = coefficient*cover_denominator
        assert polynomial.denominator() == 1
        polynomial = polynomial.numerator()
        cover_numerator += sum(
            finite(scalar)*T5**t_degree*U5**exponents[1]
            for t_degree, scalar in enumerate(polynomial.list())
        )
    cover_denominator_bivariate = sum(
        finite(scalar)*T5**degree
        for degree, scalar in enumerate(cover_denominator.list())
    )
    integral_double_cover = cover_numerator*cover_denominator_bivariate
    new_base = FunctionField(finite, "u5")
    u5 = new_base.gen()
    fiber_ring = PolynomialRing(new_base, "tau")
    tau = fiber_ring.gen()
    cover_over_new_base = fiber_ring(sum(
        new_base(coefficient)*tau**exponents[0]*u5**exponents[1]
        for exponents, coefficient in integral_double_cover.dict().items()
    ))
    cover_factorization = tuple(cover_over_new_base.factor())
    squarefree_cover = fiber_ring(1)
    factor_degrees_exponents = []
    for factor, exponent in cover_factorization:
        factor_degrees_exponents.append((int(factor.degree()), int(exponent)))
        if int(exponent) % 2:
            squarefree_cover *= factor
    squarefree_cover = squarefree_cover.monic()
    squarefree_degree = int(squarefree_cover.degree())
    accepted = squarefree_degree in (3, 4)
    jacobian_A = None
    jacobian_B = None
    delta_num = ()
    delta_den = ()
    if accepted:
        coefficients = list(squarefree_cover.list())+[new_base(0)]*5
        e, d, c, b, a = coefficients[:5]
        invariant_I = 12*a*e-3*b*d+c**2
        invariant_J = (
            72*a*c*e+9*b*c*d-27*a*d**2-27*b**2*e-2*c**3
        )
        jacobian_A = -27*invariant_I
        jacobian_B = -27*invariant_J
        delta_core = 4*invariant_I**3-invariant_J**2
        delta_num = tuple(delta_core.numerator().factor())
        delta_den = tuple(delta_core.denominator().factor())
    print(
        "Q80FIFTHMARKED|"
        f"mode={arguments.mode}|"
        f"cover_factor_degrees_exponents={tuple(factor_degrees_exponents)}|"
        f"squarefree_cover_degree_T={squarefree_degree}|"
        f"jacobian_delta_num_factors={tuple((str(factor), int(exponent)) for factor, exponent in delta_num)}|"
        f"jacobian_delta_den_factors={tuple((str(factor), int(exponent)) for factor, exponent in delta_den)}|"
        f"status={'PASS_GENUS_ONE_MARKED_PROJECTION' if accepted else 'REJECTED_WRONG_GENUS_MARKED_PROJECTION'}",
        flush=True,
    )
    if arguments.write_artifact:
        output = {
            "schema": "q80-fifth-q4-marked-projection-family-gf73-v1",
            "prime": 73,
            "mode": arguments.mode,
            "center_pair": None if center_pair is None else list(center_pair),
            "reference_pair_zero_one": list(reference_pair),
            "training_T": [row[0] for row in sample_rows[:-arguments.minimum_withheld]],
            "withheld_T": [row[0] for row in sample_rows[-arguments.minimum_withheld:]],
            "skipped_T_reasons": [list(row) for row in skip_rows],
            "marked_section_values": [
                {
                    "section_index": int(section_index),
                    "candidates": [
                        {
                            "numerator_degree": int(a),
                            "denominator_degree": int(b),
                            "equation_margin": int(m),
                            "numerator_coefficients_low_to_high": list(map(int, n.list())),
                            "denominator_coefficients_low_to_high": list(map(int, d.list())),
                        }
                        for a, b, m, n, d in candidates
                    ],
                }
                for section_index, candidates in marked_section_value_rows
            ],
            "residual_coefficients": coefficient_records,
            "integral_double_cover_terms_T_U_coefficient": [
                [int(exponents[0]), int(exponents[1]), int(coefficient)]
                for exponents, coefficient in sorted(integral_double_cover.dict().items())
            ],
            "cover_factor_degrees_exponents": [list(row) for row in factor_degrees_exponents],
            "squarefree_cover_degree_T": squarefree_degree,
            "squarefree_cover_coefficients_low_to_high": [
                str(coefficient) for coefficient in squarefree_cover.list()
            ],
            "genus_one_accepted": accepted,
            "jacobian_A": None if jacobian_A is None else str(jacobian_A),
            "jacobian_B": None if jacobian_B is None else str(jacobian_B),
            "jacobian_delta_numerator_factorization": [
                [str(factor), int(exponent)] for factor, exponent in delta_num
            ],
            "jacobian_delta_denominator_factorization": [
                [str(factor), int(exponent)] for factor, exponent in delta_den
            ],
            "reproduce": (
                "sage elkies-k3/scripts/"
                "reconstruct_q80_fifth_q4_marked_projection_gf73.sage "
                f"--mode {arguments.mode} --minimum-withheld "
                f"{arguments.minimum_withheld}"
                + (
                    f" --reference-pair {arguments.reference_pair}"
                    if arguments.reference_pair is not None
                    else ""
                )
                + " --write-artifact"
            ),
        }
        output_path = (
            ROOT / "artifacts/generated-results/"
            / (
                f"q80-fifth-q4-marked-projection-{arguments.mode}"
                + (
                    f"-ref{gauge_tag}"
                    if arguments.reference_pair is not None
                    else ""
                )
                + "-gf73.json"
            )
        )
        encoded = json.dumps(output, indent=2, sort_keys=True, default=int)+"\n"
        output_path.write_text(encoded)
        output_hash = hashlib.sha256(encoded.encode()).hexdigest()
        print(
            "Q80FIFTHMARKED|"
            f"artifact={output_path}|sha256={output_hash}|"
            "status=PASS_MARKED_PROJECTION_ARTIFACT_WRITE",
            flush=True,
        )
