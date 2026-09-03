#!/usr/bin/env sage-python
"""Compile and screen all 43 norm-twelve shared-zero R17 pencils.

The exact genus-one bisection ledger contains 43 norm-twelve trace sections.
For each one this replay applies the same Brandhorst--Elkies degree-two
Riemann--Roch compiler used for ``norm12-orbit-11952``.  The expensive section
transport is deliberately deferred: the first pass exports the Weierstrass
equation, a primitive normalized j-map, an exact PGL2-invariant critical-value
fingerprint, and exact j-preimage equations for the selected ICARM curves.

If a rational j-preimage is found, this program fails closed.  Such a hit must
be followed by a chart-specific twist/isomorphism and saturated-section
transport certificate before the result is promoted.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import gzip
import hashlib
import json
from math import gcd
from pathlib import Path
import sys

from sage.all import GF, PolynomialRing, QQ, lcm
from sage.env import SAGE_VERSION


sys.set_int_max_str_digits(0)

ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SPLITTING = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
)
CLASSIFICATION = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-isotropic-frame-classification-v1.json"
)
CURVE273 = ROOT / "artifacts/generated-results/elliptic-curves/icarm_curve273_rank30_v1.json"
CURVE302 = ROOT / "artifacts/generated-results/elliptic-curves/icarm_curve302_rank31_v1.json.gz"
LINEAGE = ROOT / "artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_lineage_v1.json"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-record-lineage-atlas-v1.json"
)
SEPARATION_PRIMES = (
    131,
    137,
    151,
    157,
    167,
)
RATIONAL_ROOT_OBSTRUCTION_PRIMES = (
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29,
    31,
    37,
    41,
    43,
    47,
    53,
    59,
    61,
    67,
    71,
    73,
    79,
    83,
    89,
    97,
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def load_json(path: Path):
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as source:
            return json.load(source)
    return json.loads(path.read_text())


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def polynomial_text(poly) -> list[str]:
    if not poly:
        return ["0"]
    return [rational_text(poly[index]) for index in range(poly.degree() + 1)]


def primitive_integer_polynomials(polynomials):
    """Normalize a polynomial tuple modulo one common nonzero QQ scalar."""

    denominators = [
        coefficient.denominator()
        for polynomial in polynomials
        for coefficient in polynomial
    ]
    denominator_lcm = lcm(denominators) if denominators else 1
    integer_coefficients = [
        int(coefficient * denominator_lcm)
        for polynomial in polynomials
        for coefficient in polynomial
    ]
    content = 0
    for coefficient in integer_coefficients:
        content = gcd(content, abs(coefficient))
    if not content:
        raise ArithmeticError("cannot normalize an all-zero polynomial tuple")
    integer_coefficients = [coefficient // content for coefficient in integer_coefficients]
    first_nonzero = next(value for value in reversed(integer_coefficients) if value)
    if first_nonzero < 0:
        integer_coefficients = [-value for value in integer_coefficients]

    result = []
    offset = 0
    for polynomial in polynomials:
        width = polynomial.degree() + 1 if polynomial else 1
        result.append(integer_coefficients[offset : offset + width])
        offset += width
    return result


def coefficient_digest(polynomials) -> str:
    payload = json.dumps(polynomials, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def exact_square_root(poly):
    if not poly.is_square():
        raise ArithmeticError("expected an exact polynomial square")
    return poly.sqrt()


def weierstrass_invariants(ainvs):
    a1, a2, a3, a4, a6 = [QQ(value) for value in ainvs]
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    delta = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    if c4**3 - c6**2 != 1728 * delta or not delta:
        raise ArithmeticError("target Weierstrass invariant identity failed")
    return c4, c6, delta


def target_records():
    curve273 = load_json(CURVE273)
    curve302 = load_json(CURVE302)
    lineage = load_json(LINEAGE)
    records = [
        {
            "curve_id": 273,
            "rank_lower_bound": 30,
            "ainvs": curve273["curve"]["ainvs"],
            "source": relative(CURVE273),
        },
        {
            "curve_id": 302,
            "rank_lower_bound": 31,
            "ainvs": curve302["curve"]["ainvs"],
            "source": relative(CURVE302),
        },
    ]
    for record in lineage["rootless_k3_interpolation_input"]["fibres"]:
        records.append(
            {
                "curve_id": int(record["curve_id"]),
                "rank_lower_bound": int(record["rank_lower_bound"]),
                "ainvs": record["source_ainvs"],
                "source": relative(LINEAGE),
            }
        )
    expected = [273, 302, 351, 356, 376, 377, 385]
    if [record["curve_id"] for record in records] != expected:
        raise ArithmeticError("record-lineage target inventory changed")
    for record in records:
        c4, c6, delta = weierstrass_invariants(record["ainvs"])
        record.update(
            {
                "c4": rational_text(c4),
                "c6": rational_text(c6),
                "discriminant": rational_text(delta),
                "j": rational_text(c4**3 / delta),
            }
        )
    return records


def compile_chart(source_record, Aold, Bold, Rt, Ru):
    """Return the polynomial short Weierstrass model for one exact trace."""

    Kt = Rt.fraction_field()
    Ku = Ru.fraction_field()
    u = Ru.gen()
    trace_data = source_record["trace_section"]
    h = Rt([QQ(value) for value in trace_data["h_coefficients_low_to_high"]])
    Nx = Rt([QQ(value) for value in trace_data["Nx_coefficients_low_to_high"]])
    Ny = Rt([QQ(value) for value in trace_data["Ny_coefficients_low_to_high"]])
    M0 = Rt([QQ(value) for value in trace_data["M0_coefficients_low_to_high"]])
    xP, yP = Kt(Nx / h**2), Kt(Ny / h**3)
    if yP**2 != xP**3 + Aold * xP + Bold or (M0 * Nx + Ny) % h**2:
        raise ArithmeticError(f"{source_record['label']} lost its trace-section identities")

    columns = [(Rt.gen() ** degree * Nx) % h**2 for degree in range(8)]
    columns += [(-Rt.gen() ** degree * Ny) % h**2 for degree in range(2)]
    from sage.all import matrix, vector

    rr_matrix = matrix(QQ, 8, 10, lambda i, j: columns[j][i])
    rr_kernel = rr_matrix.right_kernel_matrix()
    if rr_matrix.rank() != 8 or rr_kernel.nrows() != 2:
        raise ArithmeticError(f"{source_record['label']} has the wrong RR-kernel dimension")
    ab = []
    for row in rr_kernel.rows():
        row = vector(QQ, row)
        a = Rt(list(row[:8]))
        b = Rt(list(row[8:]))
        if (a * Nx - b * Ny) % h**2:
            raise ArithmeticError("Riemann--Roch congruence failed")
        ab.append((a, b))
    (a0, b0), (a1, b1) = ab

    Stu = PolynomialRing(Ku, "t")
    Ftu = Stu.fraction_field()
    lift_t = lambda value: Stu([Ku(coefficient) for coefficient in Rt(value)])
    hh, NNx, NNy = map(lift_t, (h, Nx, Ny))
    AAold = lift_t(Aold)
    xxP, yyP = Ftu(NNx / hh**2), Ftu(NNy / hh**3)
    aa0, bb0, aa1, bb1 = map(lift_t, (a0, b0, a1, b1))
    numerator_m = aa1 - u * aa0
    denominator_m = u * bb0 - bb1
    slope_m = Ftu(numerator_m / (denominator_m * hh))
    radical = (
        slope_m**4
        - 6 * xxP * slope_m**2
        - 8 * yyP * slope_m
        - 3 * xxP**2
        - 4 * AAold
    )
    radical_numerator = Stu(radical.numerator())
    radical_denominator = Stu(radical.denominator())
    square_factor = radical_numerator.gcd(radical_numerator.derivative()).monic()
    quartic, remainder = radical_numerator.quo_rem(square_factor**2)
    if remainder or quartic.degree() != 4 or quartic.gcd(quartic.derivative()).degree():
        raise ArithmeticError(f"{source_record['label']} did not compile to a smooth quartic")
    denominator_sqrt = exact_square_root(radical_denominator)
    if radical != Ftu(quartic * (square_factor / denominator_sqrt) ** 2):
        raise ArithmeticError("quartic radical identity failed")

    e, d, c, b, a = [Ku(quartic[index]) for index in range(5)]
    invariant_I = 12 * a * e - 3 * b * d + c**2
    invariant_J = (
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d**2
        - 27 * b**2 * e
        - 2 * c**3
    )
    Araw = Ku(-27 * invariant_I)
    Braw = Ku(-27 * invariant_J)
    factors_A = list(Ru(Araw.denominator()).factor())
    factors_B = list(Ru(Braw.denominator()).factor())
    if (
        len(factors_A) != 1
        or factors_A[0][1] != 8
        or len(factors_B) != 1
        or factors_B[0][1] != 12
    ):
        raise ArithmeticError(f"{source_record['label']} has an unexpected Jacobian gauge")
    ell = factors_A[0][0].monic()
    if factors_B[0][0].monic() != ell:
        raise ArithmeticError("A and B have different Jacobian poles")
    gauge = ell**2
    Achild = Ru(Araw * gauge**4)
    Bchild = Ru(Braw * gauge**6)
    delta_core = Ru(4 * Achild**3 + 27 * Bchild**2)
    if (Achild.degree(), Bchild.degree(), delta_core.degree()) != (8, 12, 24):
        raise ArithmeticError(f"{source_record['label']} lost the K3 degree profile")
    if delta_core.gcd(delta_core.derivative()).degree() or Achild.gcd(delta_core).degree():
        raise ArithmeticError(f"{source_record['label']} has a non-nodal finite fibre")
    return {
        "rr_kernel": [[rational_text(value) for value in row] for row in rr_kernel.rows()],
        "A": Achild,
        "B": Bchild,
        "delta_core": delta_core,
    }


def reduce_polynomial(poly, ring, field):
    return ring([field(coefficient) for coefficient in poly])


def critical_value_fingerprints(A, B, Ru):
    """Return modular residual critical-value polynomials of j/1728.

    For z=4*A^3/(4*A^3+27*B^2), the non-forced critical divisor is
    C=3*A'*B-2*A*B'.  It has degree 18 here.  Its image polynomial is
    invariant under rational PGL2 changes of the source coordinate.  Reduction
    of its monic degree-18 eliminant at a good prime commutes with the exact
    resultant calculation.  Thus a difference at one common good prime is an
    exact certificate that two characteristic-zero maps are inequivalent,
    without constructing enormous rational resultants.
    """

    C = 3 * A.derivative() * B - 2 * A * B.derivative()
    if C.degree() != 18 or C.gcd(A * B).degree():
        raise ArithmeticError("residual critical divisor is not finite, degree 18, and disjoint")
    fingerprints = {}
    for prime in SEPARATION_PRIMES:
        field = GF(prime)
        Rz = PolynomialRing(field, "z")
        z = Rz.gen()
        S = PolynomialRing(Rz, "u")
        try:
            AA, BB, CC = [
                S([Rz(field(coefficient)) for coefficient in poly])
                for poly in (A, B, C)
            ]
        except (ZeroDivisionError, TypeError, ValueError):
            continue
        if (AA.degree(), BB.degree(), CC.degree()) != (8, 12, 18):
            continue
        rational_map_equation = 4 * AA**3 - z * (4 * AA**3 + 27 * BB**2)
        resultant = Rz(CC.resultant(rational_map_equation))
        if not resultant or resultant.degree() != 18:
            continue
        normalized = resultant.monic()
        fingerprints[str(prime)] = [
            int(normalized[index]) for index in range(normalized.degree() + 1)
        ]
    if not fingerprints:
        raise ArithmeticError("too few good critical-value separation primes")
    return C, fingerprints


def primitive_projective_polynomial(poly, degree):
    coefficients = [QQ(poly[index]) for index in range(degree + 1)]
    denominator_lcm = lcm([value.denominator() for value in coefficients])
    integers = [int(value * denominator_lcm) for value in coefficients]
    content = 0
    for value in integers:
        content = gcd(content, abs(value))
    if not content:
        raise ArithmeticError("cannot normalize the zero projective polynomial")
    integers = [value // content for value in integers]
    first_nonzero = next(value for value in reversed(integers) if value)
    if first_nonzero < 0:
        integers = [-value for value in integers]
    return integers


def modular_no_projective_root_prime(primitive):
    for prime in RATIONAL_ROOT_OBSTRUCTION_PRIMES:
        residues = [value % prime for value in primitive]
        if all(
            sum(
                coefficient * pow(value, index, prime)
                for index, coefficient in enumerate(residues)
            )
            % prime
            for value in range(prime)
        ) and residues[-1]:
            return prime
    return None


def rational_projective_roots(poly, degree):
    primitive = primitive_projective_polynomial(poly, degree)
    obstruction_prime = modular_no_projective_root_prime(primitive)
    if obstruction_prime is not None:
        return [], obstruction_prime
    roots = [
        (QQ(root).numerator(), QQ(root).denominator())
        for root, _multiplicity in poly.roots(QQ)
    ]
    if not primitive[-1]:
        roots.append((QQ(1), QQ(0)))
    return roots, None


def normalize_pgl2_matrix(entries):
    denominators = [QQ(value).denominator() for value in entries]
    scale = lcm(denominators)
    integers = [int(QQ(value) * scale) for value in entries]
    content = 0
    for value in integers:
        content = gcd(content, abs(value))
    integers = [value // content for value in integers]
    first_nonzero = next(value for value in integers if value)
    if first_nonzero < 0:
        integers = [-value for value in integers]
    return integers


def pgl2_from_three_images(image_zero, image_one, image_infinity):
    """Map 0,1,infinity to three projective rational points."""

    z0 = tuple(map(QQ, image_zero))
    z1 = tuple(map(QQ, image_one))
    zi = tuple(map(QQ, image_infinity))
    determinant = lambda left, right: left[0] * right[1] - left[1] * right[0]
    alpha = determinant(z0, z1)
    beta = -determinant(zi, z1)
    entries = [alpha * zi[0], beta * z0[0], alpha * zi[1], beta * z0[1]]
    normalized = normalize_pgl2_matrix(entries)
    if normalized[0] * normalized[3] == normalized[1] * normalized[2]:
        raise ArithmeticError("three distinct images produced a singular PGL2 matrix")
    return normalized


def homogeneous_substitution(poly, numerator, denominator, degree):
    result = numerator.parent()(0)
    for index in range(degree + 1):
        result += poly[index] * numerator**index * denominator ** (degree - index)
    return result


def verify_pgl2_equivalence(representative, member, transform, Ru):
    """Verify j_member(phi(x))=j_representative(x) exactly."""

    x = Ru.gen()
    a, b, c, d = map(QQ, transform)
    numerator = a * x + b
    denominator = c * x + d
    member_numerator = homogeneous_substitution(member["N"], numerator, denominator, 24)
    member_denominator = homogeneous_substitution(member["D"], numerator, denominator, 24)
    return (
        member_numerator * representative["D"]
        == representative["N"] * member_denominator
    )


def exact_pgl2_equivalence(representative, member, Ru):
    """Discover a rational base change from three exact rational fibres."""

    sample_points = [QQ(0), QQ(1), QQ(-1), QQ(2), QQ(-2), QQ(3), QQ(-3)]
    samples = []
    values = set()
    for point in sample_points:
        denominator = representative["D"](point)
        if not denominator:
            continue
        value = QQ(representative["N"](point) / denominator)
        if value in values:
            continue
        samples.append((point, value))
        values.add(value)
        if len(samples) == 3:
            break
    if len(samples) != 3:
        raise ArithmeticError("could not choose three rational landmark fibres")

    root_sets = []
    obstructions = []
    for _point, value in samples:
        equation = Ru(member["N"] - value * member["D"])
        roots, obstruction = rational_projective_roots(equation, 24)
        if not roots:
            return None, {
                "landmark_values": [rational_text(item[1]) for item in samples],
                "separating_prime": obstruction,
            }
        root_sets.append(roots)
        obstructions.append(obstruction)

    def mobius_source_to_standard(point):
        p0, p1, p2 = [item[0] for item in samples]
        # psi(p0)=0, psi(p1)=1, psi(p2)=infinity.
        scale = (p1 - p2) / (p1 - p0)
        return [scale, -scale * p0, QQ(1), -p2]

    source_standard = mobius_source_to_standard(samples[0][0])
    # The inverse of psi maps 0,1,infinity back to the selected source points.
    sa, sb, sc, sd = source_standard
    source_inverse = [sd, -sb, -sc, sa]

    for image_zero in root_sets[0]:
        for image_one in root_sets[1]:
            for image_infinity in root_sets[2]:
                if len({image_zero, image_one, image_infinity}) != 3:
                    continue
                standard_to_member = pgl2_from_three_images(
                    image_zero, image_one, image_infinity
                )
                a, b, c, d = map(QQ, standard_to_member)
                e, f, g, h = map(QQ, source_standard)
                transform = normalize_pgl2_matrix(
                    [a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h]
                )
                if verify_pgl2_equivalence(representative, member, transform, Ru):
                    return transform, {
                        "landmark_source_parameters": [
                            rational_text(item[0]) for item in samples
                        ],
                        "landmark_j_values": [
                            rational_text(item[1]) for item in samples
                        ],
                        "landmark_member_images": [
                            [rational_text(value) for value in image]
                            for image in (image_zero, image_one, image_infinity)
                        ],
                        "identity_verified": True,
                    }
    return None, {
        "landmark_values": [rational_text(item[1]) for item in samples],
        "separating_prime": None,
        "rational_landmark_preimages_found_but_no_PGL2_identity": True,
    }


def exact_target_preimage(A, B, target, Ru):
    child_c4 = -48 * A
    child_delta = -16 * (4 * A**3 + 27 * B**2)
    target_c4 = QQ(target["c4"])
    target_delta = QQ(target["discriminant"])
    preimage = Ru(child_c4**3 * target_delta - target_c4**3 * child_delta)
    if not preimage:
        raise ArithmeticError("a nonconstant atlas j-map became a constant target j-map")
    projective_degree = 24
    primitive = primitive_projective_polynomial(preimage, projective_degree)
    obstruction_prime = modular_no_projective_root_prime(primitive)
    rational_roots = []
    root_at_infinity = not primitive[-1]
    decision_method = "projective modular no-root obstruction"
    if obstruction_prime is None:
        rational_roots = sorted({QQ(root) for root, _multiplicity in preimage.roots(QQ)})
        decision_method = "exact QQ rational-root fallback"
    return {
        "primitive_projective_equation_coefficients_u_degree_0_through_24": [
            str(value) for value in primitive
        ],
        "finite_degree": int(preimage.degree()),
        "projective_degree": projective_degree,
        "decision_method": decision_method,
        "modular_no_projective_root_prime": obstruction_prime,
        "finite_rational_roots": [rational_text(root) for root in rational_roots],
        "root_at_infinity": root_at_infinity,
        "rational_match": bool(rational_roots or root_at_infinity),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    args.output = args.output.resolve()

    model = load_json(MODEL)
    splitting = load_json(SPLITTING)
    classification = load_json(CLASSIFICATION)
    targets = target_records()
    frame_by_label = {
        record["label"]: record["frame_class"]
        for record in classification["classification"]["records"]
    }
    source_records = [
        record
        for record in splitting["construction"]["records"]
        if int(record["trace_norm"]) == 12
    ]
    if len(source_records) != 43 or set(frame_by_label) != {
        record["label"] for record in source_records
    }:
        raise ArithmeticError("the certified 43-chart inventory changed")

    Rt = PolynomialRing(QQ, "t")
    Ru = PolynomialRing(QQ, "u")
    Aold = Rt([QQ(value) for value in model["A_coefficients_low_to_high"]])
    Bold = Rt([QQ(value) for value in model["B_coefficients_low_to_high"]])

    charts = []
    internal_models = {}
    rational_hits = []
    for index, source_record in enumerate(source_records, 1):
        label = source_record["label"]
        compiled = compile_chart(source_record, Aold, Bold, Rt, Ru)
        A, B, delta_core = compiled["A"], compiled["B"], compiled["delta_core"]
        j_numerator = Ru(6912 * A**3)
        j_denominator = delta_core
        if j_numerator.gcd(j_denominator).degree():
            raise ArithmeticError(f"{label} has a non-reduced j-map")
        normalized_pair = primitive_integer_polynomials([j_numerator, j_denominator])
        j_digest = coefficient_digest(normalized_pair)
        critical, branch_values = critical_value_fingerprints(A, B, Ru)
        branch_digest = coefficient_digest([branch_values])
        internal_models[label] = {"N": j_numerator, "D": j_denominator}

        matches = {}
        for target in targets:
            match = exact_target_preimage(A, B, target, Ru)
            matches[str(target["curve_id"])] = match
            if match["rational_match"]:
                rational_hits.append(
                    {
                        "chart": label,
                        "curve_id": target["curve_id"],
                        "finite_roots": match["finite_rational_roots"],
                        "root_at_infinity": match["root_at_infinity"],
                    }
                )

        charts.append(
            {
                "label": label,
                "frame_class": frame_by_label[label],
                "trace_vector": source_record["pinned_rank17_w"],
                "equation_complexity": source_record["equation_complexity"],
                "riemann_roch": {
                    "constraint_rank": 8,
                    "kernel_dimension": 2,
                    "kernel_rows_a0_through_a7_b0_b1": compiled["rr_kernel"],
                },
                "weierstrass_model": {
                    "equation": "Y^2=X^3+A(u)*X+B(u)",
                    "A_coefficients_low_to_high": polynomial_text(A),
                    "B_coefficients_low_to_high": polynomial_text(B),
                    "degrees_A_B_Delta": [8, 12, 24],
                    "finite_discriminant_squarefree": True,
                },
                "normalized_j_map": {
                    "formula": "j(u)=N(u)/D(u)",
                    "numerator_coefficients_low_to_high": [
                        str(value) for value in normalized_pair[0]
                    ],
                    "denominator_coefficients_low_to_high": [
                        str(value) for value in normalized_pair[1]
                    ],
                    "projective_degree": 24,
                    "sha256": j_digest,
                },
                "pgl2_invariant": {
                    "method": (
                        "monic critical-value resultant for the degree-18 residual "
                        "critical divisor of j/1728"
                    ),
                    "residual_critical_coefficients_low_to_high": polynomial_text(critical),
                    "good_prime_monic_critical_value_polynomials_low_to_high": branch_values,
                    "critical_value_polynomial_degree": 18,
                    "sha256": branch_digest,
                },
                "target_j_preimages": matches,
            }
        )
        print(
            f"R17RECORDLINEAGE|chart={index}/43|label={label}|frame={frame_by_label[label]}|"
            f"j={j_digest[:12]}|branch={branch_digest[:12]}|hits="
            f"{sum(int(row['rational_match']) for row in matches.values())}",
            flush=True,
        )

    equivalence_classes = []
    direct_rejections = []
    for chart in charts:
        label = chart["label"]
        member_values = chart["pgl2_invariant"][
            "good_prime_monic_critical_value_polynomials_low_to_high"
        ]
        matched = False
        for equivalence_class in equivalence_classes:
            representative_label = equivalence_class["representative"]
            representative_chart = next(
                row for row in charts if row["label"] == representative_label
            )
            representative_values = representative_chart["pgl2_invariant"][
                "good_prime_monic_critical_value_polynomials_low_to_high"
            ]
            common_primes = sorted(
                set(member_values) & set(representative_values), key=int
            )
            separating_prime = next(
                (
                    prime
                    for prime in common_primes
                    if member_values[prime] != representative_values[prime]
                ),
                None,
            )
            if separating_prime is not None:
                continue
            transform, solve = exact_pgl2_equivalence(
                internal_models[representative_label], internal_models[label], Ru
            )
            if transform is not None:
                equivalence_class["members"].append(
                    {
                        "label": label,
                        "representative_to_member_pgl2_matrix_a_b_c_d": transform,
                        "solve_certificate": solve,
                    }
                )
                matched = True
                break
            direct_rejections.append(
                {
                    "left": representative_label,
                    "right": label,
                    "critical_value_fingerprint_collision": True,
                    "direct_solve": solve,
                }
            )
        if not matched:
            equivalence_classes.append(
                {
                    "representative": label,
                    "members": [
                        {
                            "label": label,
                            "representative_to_member_pgl2_matrix_a_b_c_d": [1, 0, 0, 1],
                            "solve_certificate": {"identity_verified": True},
                        }
                    ],
                }
            )
    if sorted(
        member["label"]
        for equivalence_class in equivalence_classes
        for member in equivalence_class["members"]
    ) != sorted(internal_models):
        raise ArithmeticError("PGL2 quotient lost an atlas chart")

    counts = defaultdict(int)
    for chart in charts:
        counts[chart["frame_class"]] += 1
    if dict(counts) != {"published-R17": 33, "alternate-Q80": 10}:
        raise ArithmeticError("frame-class counts changed")

    payload = {
        "schema": "elkies-k3.r17-norm12-record-lineage-atlas.v1",
        "status": (
            "PASS_EXACT_43_CHART_COMPILATION_PGL2_SEPARATION_NO_RATIONAL_RECORD_MATCH"
            if not rational_hits
            else "RATIONAL_J_MATCH_REQUIRES_TWIST_AND_SECTION_FOLLOWUP"
        ),
        "atlas": {
            "chart_count": len(charts),
            "frame_class_counts": dict(sorted(counts.items())),
            "pgl2_equivalence_class_count": len(equivalence_classes),
            "pgl2_equivalence_classes": equivalence_classes,
            "critical_value_collision_direct_rejections": direct_rejections,
            "separation_theorem": (
                "Maps in one displayed class satisfy the stored exact rational-function "
                "identity after the displayed PGL2(Q) change. Distinct classes are "
                "separated either by a common-good-prime reduction of the exact residual "
                "critical-value eliminant or by the exact three-landmark PGL2 solve."
            ),
            "charts": charts,
        },
        "targets": [
            {
                key: value
                for key, value in target.items()
                if key not in {"c4", "c6", "discriminant"}
            }
            for target in targets
        ],
        "rational_j_matches": rational_hits,
        "conditional_followup": {
            "triggered": bool(rational_hits),
            "required_on_hit": [
                "separate the rational quadratic twist class",
                "construct and verify the Q-isomorphism when the twist is trivial",
                "transport and saturate the seventeen generic sections",
                "certify the exact exceptional quotient against the target point subgroup",
            ],
        },
        "method": {
            "compiler": "Brandhorst--Elkies degree-two shared-zero Riemann--Roch kernel",
            "j_normalization": (
                "primitive common-integer normalization of "
                "(6912*A(u)^3, 4*A(u)^3+27*B(u)^2), with positive last nonzero coefficient"
            ),
            "pgl2_quotient": (
                "good-prime reductions of the exact degree-18 residual critical-value "
                "eliminant followed by exact three-landmark PGL2 discovery and identity checks"
            ),
            "target_gate": (
                "form the exact primitive projective degree-24 equation "
                "c4_chart(u)^3*Delta_target-c4_target^3*Delta_chart(u)=0 over QQ; "
                "certify absence of rational points by a projective modular no-root prime, "
                "falling back to exact QQ rational-root calculation when needed"
            ),
        },
        "claim_boundary": {
            "proved": [
                "all 43 certified norm-twelve shared-zero divisors compile to polynomial (8,12) Weierstrass charts",
                "the displayed classes are the rational-PGL2 quotient of the 43 j-maps",
                "the displayed exact target-preimage equations and modular obstructions decide rational j-accessibility for the seven pinned ICARM curves",
            ],
            "not_proved": [
                "a J1 surface-automorphism classification beyond j-map inequivalence",
                "absence from rootless charts outside the certified 43-member shared-zero degree-two atlas",
                "any rank upper bound for an ICARM target",
            ],
        },
        "inputs": {
            relative(path): digest(path)
            for path in (MODEL, SPLITTING, CLASSIFICATION, CURVE273, CURVE302, LINEAGE)
        },
        "software_assumptions": {"sage_version": SAGE_VERSION},
        "reproducing_command": (
            "sage -python elkies-k3/scripts/compile_r17_norm12_record_lineage_atlas.sage"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != serialized:
            raise ArithmeticError("stored record-lineage atlas differs from exact replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        "R17RECORDLINEAGE|charts=43|pgl2_classes={}|targets=7|rational_hits={}|status={}|output={}".format(
            len(equivalence_classes),
            len(rational_hits),
            payload["status"],
            relative(args.output),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
