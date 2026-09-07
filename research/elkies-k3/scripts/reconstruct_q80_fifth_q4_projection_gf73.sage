#!/usr/bin/env sage
"""Interpolate the compensated fifth-q4 cubic/projection over GF(73)(T).

Input samples are produced independently by
``probe_q80_fifth_q4_integral_basis_gf73.sage``.  Each sample obtains the
three-dimensional pole space by a finite/infinite maximal-order lattice
intersection and projects the resulting genuine plane cubic from the pinned
degree-one section.  This script seeks only rational coefficient functions
that are overdetermined by the available samples; an exactly determined
Padé fit is never reported as a certificate.  The current independently
depressed fiberwise coordinate reconstructs but has squarefree branch degree
12, so the script reports a passing rejection diagnostic and never writes it
as a fifth-family artifact.
"""

import argparse
import glob
import hashlib
import json
from pathlib import Path

from sage.all import FunctionField, GF, PolynomialRing, lcm


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument("--minimum-withheld", type=int, default=2)
parser.add_argument("--max-total-degree", type=int, default=20)
parser.add_argument("--raw", action="store_true")
parser.add_argument("--write-artifact", action="store_true")
arguments = parser.parse_args()

paths = tuple(sorted(
    Path(path) for path in glob.glob(str(
        ROOT / "artifacts/generated-results/"
        "q80-fifth-q4-compensated-projection-gf73-T*.json"
    ))
))
samples = []
for path in paths:
    payload = json.loads(path.read_text())
    if payload.get("schema") != (
        "q80-fifth-q4-compensated-projection-sample-gf73-v1"
    ):
        continue
    samples.append((payload["T"], payload, path))
samples.sort()
if len(samples) < 4:
    raise RuntimeError("need at least four compensated projection samples")

finite = GF(73, impl="modn")
parameter_ring = PolynomialRing(finite, "t")
t = parameter_ring.gen()
coordinate_ring = PolynomialRing(finite, names=("v", "w", "u"))
v, w, u = coordinate_ring.gens()


def nullspace(rows, column_count):
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


def rational_interpolate(values):
    if len(values) <= arguments.minimum_withheld:
        return []
    training_values = values[:-arguments.minimum_withheld]
    candidates = []
    for total_degree in range(arguments.max_total_degree+1):
        for numerator_degree in range(total_degree+1):
            denominator_degree = total_degree-numerator_degree
            unknowns = numerator_degree+denominator_degree+2
            if unknowns > len(training_values)+1:
                continue
            rows = []
            for parameter, value in training_values:
                parameter = finite(parameter)
                value = finite(value)
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
                   != finite(value)*denominator(finite(parameter))
                for parameter, value in values
            ):
                continue
            scale = denominator[denominator.degree()]**(-1)
            numerator *= scale
            denominator *= scale
            candidates.append(
                (total_degree, numerator_degree, denominator_degree,
                 len(values)-rank, numerator, denominator)
            )
        if candidates:
            break
    return candidates


def raw_polynomial(payload, key):
    return coordinate_ring({
        tuple(row[:-1]): finite(row[-1])
        for row in payload[key]
    })


def canonical_polynomials(payload):
    plane = raw_polynomial(
        payload, "plane_cubic_terms_v_w_u_coefficient"
    )
    residual = raw_polynomial(
        payload, "projection_residual_terms_v_w_u_coefficient"
    )
    if arguments.raw:
        return plane, residual
    w_squared_linear = (
        plane.monomial_coefficient(v*w**2)*v
        + plane.monomial_coefficient(w**2)
    )
    assert w_squared_linear.degree(v) <= 1
    translation = w_squared_linear/finite(3)
    depressed = coordinate_ring(plane(w=w-translation))
    assert depressed.monomial_coefficient(v*w**2) == 0
    assert depressed.monomial_coefficient(w**2) == 0
    slope_shift = translation.monomial_coefficient(v)
    shifted_residual = coordinate_ring(residual(u=u-slope_shift))
    return depressed, shifted_residual


canonical_samples = [
    (parameter, *canonical_polynomials(payload), payload, path)
    for parameter, payload, path in samples
]

section_x_coefficients = (
    (41, 12, 67, 72, 43),
    (24, 63, 53, 61, 41),
    (7, 47, 45, 69, 62),
    (12, 26, 13, 39, 49),
    (48, 60, 25, 12, 30),
)


def section_v_values(parameter):
    parameter = finite(parameter)
    def ratio(numerator, denominator):
        return None if denominator == 0 else numerator/denominator
    return (
        ratio(8-33*parameter, parameter+16),
        ratio(-27-13*parameter, parameter+36),
        ratio(18-36*parameter, parameter-32),
        ratio(30+34*parameter, parameter+4),
        ratio(6+3*parameter, parameter+35),
    )


root_ring = PolynomialRing(finite, "root_w")
root_w = root_ring.gen()
section_root_rows = []
for parameter, plane, _, payload, _ in canonical_samples:
    roots_for_sample = []
    for section_v_value in section_v_values(parameter):
        if section_v_value is None:
            roots_for_sample.append(())
            continue
        polynomial = root_ring(sum(
            coefficient*section_v_value**exponents[0]*root_w**exponents[1]
            for exponents, coefficient in plane.dict().items()
        ))
        roots_for_sample.append(tuple(polynomial.roots(multiplicities=False)))
    section_root_rows.append((parameter, tuple(roots_for_sample)))
unique_section_roots = sum(
    len(roots) == 1
    for _, roots_for_sample in section_root_rows for roots in roots_for_sample
)
print(
    "Q80FIFTHINTERPOLATE|"
    f"section_root_unique={unique_section_roots}/{5*len(section_root_rows)}|"
    f"root_count_patterns={tuple(sorted(set(tuple(len(roots) for roots in row) for _, row in section_root_rows)))}|"
    "status=PASS_EXPLICIT_SECTION_ROOT_AUDIT",
    flush=True,
)


def term_dictionary(polynomial):
    return polynomial.dict()


plane_support = sorted(set().union(*(
    term_dictionary(plane) for _, plane, _, _, _ in canonical_samples
)))
residual_support = sorted(set().union(*(
    term_dictionary(residual) for _, _, residual, _, _ in canonical_samples
)))


def reconstruct_support(polynomial_index, support):
    rows = []
    for monomial in support:
        values = []
        for sample in canonical_samples:
            parameter = sample[0]
            terms = term_dictionary(sample[polynomial_index])
            values.append((parameter, terms.get(monomial, finite(0))))
        candidates = rational_interpolate(values)
        rows.append((monomial, candidates))
    return rows


plane_rows = reconstruct_support(
    1, plane_support
)
residual_rows = reconstruct_support(
    2, residual_support
)
all_rows = plane_rows+residual_rows
recognized = sum(bool(candidates) for _, candidates in all_rows)
print(
    "Q80FIFTHINTERPOLATE|"
    f"samples={tuple(parameter for parameter, _, _ in samples)}|"
    f"sample_count={len(samples)}|plane_support={len(plane_support)}|"
    f"residual_support={len(residual_support)}|"
    f"recognized={recognized}/{len(all_rows)}|"
    f"minimum_withheld={arguments.minimum_withheld}|"
    f"depressed_w={int(not arguments.raw)}|"
    "status=PASS_BOUNDED_RECOGNITION_AUDIT",
    flush=True,
)
for label, rows in (("plane", plane_rows), ("residual", residual_rows)):
    for monomial, candidates in rows:
        compact = tuple(
            (num_degree, den_degree, margin, str(numerator), str(denominator))
            for _, num_degree, den_degree, margin, numerator, denominator
            in candidates
        )
        print(
            "Q80FIFTHINTERPOLATE|"
            f"kind={label}|monomial={monomial}|candidates={compact}",
            flush=True,
        )

if recognized == len(all_rows) and not arguments.raw:
    parameter_field = parameter_ring.fraction_field()
    family_ring = PolynomialRing(
        parameter_field, names=("vf", "wf", "uf")
    )
    vf, wf, uf = family_ring.gens()

    def reconstructed_polynomial(rows):
        value = family_ring(0)
        for monomial, candidates in rows:
            if len(candidates) != 1:
                raise ArithmeticError("coefficient recognition is not unique")
            _, _, _, _, numerator, denominator = candidates[0]
            coefficient = parameter_field(numerator)/parameter_field(denominator)
            value += coefficient*vf**monomial[0]*wf**monomial[1]*uf**monomial[2]
        return value

    plane_family = reconstructed_polynomial(plane_rows)
    residual_family = reconstructed_polynomial(residual_rows)
    section_v = parameter_field((34*t+30)/(t+4))
    section_w = (
        residual_family.monomial_coefficient(vf*uf**2)/parameter_field(3)
    )
    projected_numerator = plane_family.subs(
        {wf: section_w+uf*(vf-section_v)}
    )
    assert projected_numerator == (vf-section_v)*residual_family
    assert residual_family.degree(vf) == 2
    quadratic_a = family_ring(sum(
        coefficient*uf**exponents[2]
        for exponents, coefficient in residual_family.dict().items()
        if exponents[0] == 2
    ))
    quadratic_b = family_ring(sum(
        coefficient*uf**exponents[2]
        for exponents, coefficient in residual_family.dict().items()
        if exponents[0] == 1
    ))
    quadratic_c = family_ring(sum(
        coefficient*uf**exponents[2]
        for exponents, coefficient in residual_family.dict().items()
        if exponents[0] == 0
    ))
    double_cover = family_ring(quadratic_b**2-4*quadratic_a*quadratic_c)
    double_cover_coefficients = {
        exponents[2]: coefficient
        for exponents, coefficient in double_cover.dict().items()
    }
    cover_denominator = lcm([
        coefficient.denominator()
        for coefficient in double_cover_coefficients.values()
    ])
    cover_ring = PolynomialRing(finite, names=("T5", "U5"))
    T5, U5 = cover_ring.gens()
    cover_numerator = cover_ring(0)
    for u_degree, coefficient in double_cover_coefficients.items():
        polynomial = coefficient*cover_denominator
        assert polynomial.denominator() == 1
        polynomial = polynomial.numerator()
        cover_numerator += sum(
            finite(scalar)*T5**t_degree*U5**u_degree
            for t_degree, scalar in enumerate(polynomial.list())
        )
    cover_denominator_bivariate = sum(
        finite(scalar)*T5**t_degree
        for t_degree, scalar in enumerate(cover_denominator.list())
    )
    integral_double_cover = cover_numerator*cover_denominator_bivariate

    new_base_field = FunctionField(finite, "u5")
    u5 = new_base_field.gen()
    fiber_ring = PolynomialRing(new_base_field, "tau")
    tau = fiber_ring.gen()
    integral_over_new_base = fiber_ring(sum(
        new_base_field(coefficient)*tau**exponents[0]*u5**exponents[1]
        for exponents, coefficient in integral_double_cover.dict().items()
    ))
    cover_factorization = tuple(integral_over_new_base.factor())
    squarefree_cover = fiber_ring(1)
    factor_degrees_exponents = []
    for factor, exponent in cover_factorization:
        factor_degrees_exponents.append((int(factor.degree()), int(exponent)))
        if int(exponent) % 2:
            squarefree_cover *= factor
    squarefree_cover = squarefree_cover.monic()
    if squarefree_cover.degree() not in (3, 4):
        print(
            "Q80FIFTHINTERPOLATE|"
            f"cover_factor_degrees_exponents={tuple(factor_degrees_exponents)}|"
            f"squarefree_cover_degree_T={squarefree_cover.degree()}|"
            "reason=fiberwise_W_gauge_not_coherent_global_pencil|"
            "status=REJECTED_WRONG_GENUS_FIBERWISE_GAUGE",
            flush=True,
        )
        if arguments.write_artifact:
            print(
                "Q80FIFTHINTERPOLATE|artifact_write=SKIPPED_REJECTED_MODEL",
                flush=True,
            )
        raise SystemExit(0)
    quartic_coefficients = list(squarefree_cover.list())+[new_base_field(0)]*5
    e, d, c, b, a = quartic_coefficients[:5]
    quartic_I = 12*a*e-3*b*d+c**2
    quartic_J = 72*a*c*e+9*b*c*d-27*a*d**2-27*b**2*e-2*c**3
    jacobian_A = -27*quartic_I
    jacobian_B = -27*quartic_J
    jacobian_delta_core = 4*quartic_I**3-quartic_J**2
    delta_numerator_factorization = tuple(
        jacobian_delta_core.numerator().factor()
    )
    delta_denominator_factorization = tuple(
        jacobian_delta_core.denominator().factor()
    )
    print(
        "Q80FIFTHINTERPOLATE|"
        f"section_v={section_v}|section_w={section_w}|"
        f"cover_factor_degrees_exponents={tuple(factor_degrees_exponents)}|"
        f"squarefree_cover_degree_T={squarefree_cover.degree()}|"
        f"jacobian_delta_num_factors={tuple((str(factor), int(exponent)) for factor, exponent in delta_numerator_factorization)}|"
        f"jacobian_delta_den_factors={tuple((str(factor), int(exponent)) for factor, exponent in delta_denominator_factorization)}|"
        "status=PASS_FIFTH_BINARY_QUARTIC_JACOBIAN",
        flush=True,
    )

    if arguments.write_artifact:
        def coefficient_record(monomial, candidate):
            _, numerator_degree, denominator_degree, margin, numerator, denominator = candidate
            return {
                "monomial": list(map(int, monomial)),
                "numerator_coefficients_low_to_high": list(map(int, numerator.list())),
                "denominator_coefficients_low_to_high": list(map(int, denominator.list())),
                "numerator_degree": int(numerator_degree),
                "denominator_degree": int(denominator_degree),
                "equation_margin": int(margin),
            }

        output = {
            "schema": "q80-fifth-q4-compensated-projection-family-gf73-v1",
            "prime": 73,
            "training_T": [int(parameter) for parameter, _, _ in samples[:-arguments.minimum_withheld]],
            "withheld_T": [int(parameter) for parameter, _, _ in samples[-arguments.minimum_withheld:]],
            "source_sample_paths": [str(path.relative_to(ROOT)) for _, _, path in samples],
            "coordinate_normalization": (
                "depress W by W_new=W_old+(coeff(V*W_old^2)*V+coeff(W_old^2))/3"
            ),
            "plane_cubic_coefficients": [
                coefficient_record(monomial, candidates[0])
                for monomial, candidates in plane_rows
            ],
            "projection_residual_coefficients": [
                coefficient_record(monomial, candidates[0])
                for monomial, candidates in residual_rows
            ],
            "section_v": str(section_v),
            "section_w": str(section_w),
            "projection_identity_verified": True,
            "integral_double_cover_terms_T_U_coefficient": [
                [int(exponents[0]), int(exponents[1]), int(coefficient)]
                for exponents, coefficient in sorted(integral_double_cover.dict().items())
            ],
            "cover_factor_degrees_exponents": [
                list(row) for row in factor_degrees_exponents
            ],
            "squarefree_cover_coefficients_low_to_high": [
                str(coefficient) for coefficient in squarefree_cover.list()
            ],
            "jacobian_A": str(jacobian_A),
            "jacobian_B": str(jacobian_B),
            "jacobian_delta_numerator_factorization": [
                [str(factor), int(exponent)]
                for factor, exponent in delta_numerator_factorization
            ],
            "jacobian_delta_denominator_factorization": [
                [str(factor), int(exponent)]
                for factor, exponent in delta_denominator_factorization
            ],
            "reproduce": (
                "sage elkies-k3/scripts/"
                "reconstruct_q80_fifth_q4_projection_gf73.sage "
                "--minimum-withheld 8 --write-artifact"
            ),
        }
        output_path = (
            ROOT / "artifacts/generated-results/"
            "q80-fifth-q4-compensated-projection-family-gf73.json"
        )
        encoded = json.dumps(output, indent=2, sort_keys=True, default=int)+"\n"
        output_path.write_text(encoded)
        output_hash = hashlib.sha256(encoded.encode()).hexdigest()
        print(
            "Q80FIFTHINTERPOLATE|"
            f"artifact={output_path}|sha256={output_hash}|"
            "status=PASS_FIFTH_FAMILY_ARTIFACT_WRITE",
            flush=True,
        )
