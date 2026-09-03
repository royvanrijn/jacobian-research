#!/usr/bin/env sage-python
"""Construct the 121 inherited bisections on the alternate-Q80 endpoint.

The direct norm12/orbit11952 fibration has fibre

    D = (3,2,w) in U + R17(-1).

An old published-R17 height-four section ``S_v`` has degree

    D.S_v = 5 - <w,v>.

Exactly 121 oriented height-four sections have degree two.  Restricting the
compiled pencil coordinate ``u=L1/L0`` to each such rational curve gives a
quadratic map P1_t -> P1_u.  This script records its exact quadratic relation,
canonical squareclass (including the rational constant squareclass), and the
lifted section on the alternate-Q80 equation over ``QQ(u,s)``, ``s^2=q(u)``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    EllipticCurve,
    PolynomialRing,
    QQ,
    ZZ,
    matrix,
    pari,
    vector,
)
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
TARGET = ROOT / "artifacts/generated-results/elkies-2026-published-r17-target.json"
SPLITTING = ROOT / "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
DIRECT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-inherited-bisection-covers-v1.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def load_matrix(path: Path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def polynomial_text(poly) -> list[str]:
    if not poly:
        return ["0"]
    return [rational_text(poly[index]) for index in range(poly.degree() + 1)]


def rational_function_record(value):
    value = value.parent()(value)
    return {
        "numerator_coefficients_low_to_high": polynomial_text(value.numerator()),
        "denominator_coefficients_low_to_high": polynomial_text(value.denominator()),
    }


def evaluate_polynomial(poly, value):
    result = value.parent()(0)
    for coefficient in reversed(list(poly)):
        result = result * value + value.parent()(coefficient)
    return result


def evaluate_rational(function, value):
    return evaluate_polynomial(function.numerator(), value) / evaluate_polynomial(
        function.denominator(), value
    )


def reconstruct_basis(ring, A, B, section_data):
    points = []
    for expected_index, record in enumerate(section_data["sections"]):
        if int(record["basis_index"]) != expected_index:
            raise ArithmeticError("published section basis order changed")
        x_coordinate = ring([QQ(value) for value in record["x_coefficients_low_to_high"]])
        if expected_index == 0:
            y_coordinate = ring([QQ(value) for value in record["y_coefficients_low_to_high"]])
        else:
            reference_x, reference_y = points[int(record["chord"]["reference_basis_index"])]
            slope = ring([QQ(value) for value in record["chord"]["slope_coefficients_low_to_high"]])
            y_coordinate = reference_y + slope * (x_coordinate - reference_x)
        if y_coordinate**2 != x_coordinate**3 + A * x_coordinate + B:
            raise ArithmeticError("published section failed its equation")
        points.append((x_coordinate, y_coordinate))
    return points


def find_record(payload, label):
    for record in payload["construction"]["records"]:
        if record["label"] == label:
            return record
    raise KeyError(label)


def exact_square_root(poly):
    if not poly.is_square():
        raise ArithmeticError("expected an exact polynomial square")
    return poly.sqrt()


def squarefree_rational_representative(value) -> ZZ:
    """Canonical representative of a nonzero element of QQ*/QQ*2.

    Denominator primes are moved to the numerator, since p^-1 and p differ
    by the rational square p^-2.  The result is a signed squarefree integer.
    """

    value = QQ(value)
    if not value:
        raise ArithmeticError("zero has no squareclass")
    result = ZZ(-1 if value < 0 else 1)
    for integer in (abs(ZZ(value.numerator())), ZZ(value.denominator())):
        for prime, exponent in integer.factor():
            if exponent % 2:
                result *= prime
    return result


def canonical_squareclass(poly, ring):
    """Return q_can and g with poly=g^2*q_can, preserving constants."""

    poly = ring(poly)
    if not poly:
        raise ArithmeticError("zero branch polynomial")
    factorization = poly.factor()
    constant = QQ(factorization.unit())
    constant_class = squarefree_rational_representative(constant)
    support = ring.one()
    for factor, exponent in factorization:
        if exponent % 2:
            support *= factor.monic()
    canonical = ring(constant_class * support)
    ratio = ring.fraction_field()(poly / canonical)
    numerator = ring(ratio.numerator())
    denominator = ring(ratio.denominator())
    multiplier = ring.fraction_field()(
        exact_square_root(numerator) / exact_square_root(denominator)
    )
    if multiplier**2 * canonical != poly:
        raise ArithmeticError("canonical squareclass multiplier identity failed")
    return canonical, multiplier, constant_class


def quadratic_element_record(value, base_field):
    coefficients = list(value)
    coefficients += [base_field(0)] * (2 - len(coefficients))
    if len(coefficients) != 2:
        raise ArithmeticError("quadratic extension element has unexpected degree")
    return {
        "constant": rational_function_record(base_field(coefficients[0])),
        "sqrt_coefficient": rational_function_record(base_field(coefficients[1])),
    }


def published_expression(entries) -> str:
    terms = []
    for index, coefficient in enumerate(entries, start=1):
        coefficient = int(coefficient)
        if not coefficient:
            continue
        magnitude = abs(coefficient)
        atom = f"P{index}" if magnitude == 1 else f"{magnitude}*P{index}"
        if not terms:
            terms.append(atom if coefficient > 0 else f"-{atom}")
        else:
            terms.append(("+" if coefficient > 0 else "-") + atom)
    return "".join(terms) if terms else "O"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    model = json.loads(MODEL.read_text())
    section_data = json.loads(SECTIONS.read_text())
    target = json.loads(TARGET.read_text())
    splitting = json.loads(SPLITTING.read_text())
    direct = json.loads(DIRECT.read_text())
    if direct["status"] != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
        raise ArithmeticError("canonical direct alternate-Q80 artifact is not certified")
    if direct["weierstrass_model"]["fibre_configuration"] != "24 I1":
        raise ArithmeticError("canonical alternate-Q80 model is not 24I1")
    if direct["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS":
        raise ArithmeticError("canonical alternate-Q80 section basis is not saturated")
    if direct["sections"]["rank"] != 17 or direct["sections"]["height_gram_determinant"] != 948:
        raise ArithmeticError("canonical alternate-Q80 rank-17 lattice changed")

    Rt = PolynomialRing(QQ, "t")
    Kt = Rt.fraction_field()
    t = Rt.gen()
    Aold = Rt([QQ(value) for value in model["A_coefficients_low_to_high"]])
    Bold = Rt([QQ(value) for value in model["B_coefficients_low_to_high"]])
    Delta_old = Rt(-16 * (4 * Aold**3 + 27 * Bold**2))
    basis_xy = reconstruct_basis(Rt, Aold, Bold, section_data)
    Eold = EllipticCurve(Kt, [Aold, Bold])
    published_basis = [Eold(Kt(x), Kt(y)) for x, y in basis_xy]

    source_record = find_record(splitting, "norm12-orbit-11952")
    trace_data = source_record["trace_section"]
    h = Rt([QQ(value) for value in trace_data["h_coefficients_low_to_high"]])
    Nx = Rt([QQ(value) for value in trace_data["Nx_coefficients_low_to_high"]])
    Ny = Rt([QQ(value) for value in trace_data["Ny_coefficients_low_to_high"]])
    xP, yP = Kt(Nx / h**2), Kt(Ny / h**3)
    if yP**2 != xP**3 + Aold * xP + Bold:
        raise ArithmeticError("norm12/orbit11952 trace section failed")

    stored_kernel = direct["riemann_roch"]["kernel_rows_a0_through_a7_b0_b1"]
    kernel_rows = [vector(QQ, [QQ(value) for value in row]) for row in stored_kernel]
    if len(kernel_rows) != 2:
        raise ArithmeticError("stored alternate-Q80 pencil is not two-dimensional")
    ab = []
    for row in kernel_rows:
        a = Rt(list(row[:8]))
        b = Rt(list(row[8:]))
        if (a * Nx - b * Ny) % h**2:
            raise ArithmeticError("stored alternate-Q80 pencil row failed regularity")
        ab.append((a, b))
    (a0, b0), (a1, b1) = ab

    Ru = PolynomialRing(QQ, "u")
    Ku = Ru.fraction_field()
    u = Ru.gen()
    Stu = PolynomialRing(Ku, "t")
    Ftu = Stu.fraction_field()
    lift_t = lambda value: Stu([Ku(coefficient) for coefficient in Rt(value)])
    hh, NNx, NNy = map(lift_t, (h, Nx, Ny))
    xxP, yyP = Ftu(NNx / hh**2), Ftu(NNy / hh**3)
    AAold = lift_t(Aold)
    aa0, bb0, aa1, bb1 = map(lift_t, (a0, b0, a1, b1))
    numerator_m = aa1 - u * aa0
    denominator_m = u * bb0 - bb1
    slope_m = Ftu(numerator_m / (denominator_m * hh))
    radical = slope_m**4 - 6 * xxP * slope_m**2 - 8 * yyP * slope_m - 3 * xxP**2 - 4 * AAold
    radical_numerator = Stu(radical.numerator())
    radical_denominator = Stu(radical.denominator())
    square_factor = radical_numerator.gcd(radical_numerator.derivative()).monic()
    quartic, remainder = radical_numerator.quo_rem(square_factor**2)
    denominator_sqrt = exact_square_root(radical_denominator)
    if remainder or quartic.degree() != 4:
        raise ArithmeticError("canonical alternate-Q80 quartic reconstruction failed")
    if radical != Ftu(quartic * (square_factor / denominator_sqrt) ** 2):
        raise ArithmeticError("canonical alternate-Q80 radical identity failed")

    t0 = Ku(-denominator_m[0] / denominator_m[1])
    normalization = Ftu(denominator_sqrt / denominator_m**2)
    normalization_u = Ku(Stu(normalization.numerator())[0] / Stu(normalization.denominator())[0])
    v0 = Ku(normalization_u * numerator_m(t0) ** 2 / (hh(t0) ** 2 * square_factor(t0)))
    if v0**2 != quartic(t0):
        raise ArithmeticError("canonical alternate-Q80 quartic point failed")

    Sz = PolynomialRing(Ku, "z")
    z = Sz.gen()
    shifted = Sz(quartic(t0 + z))
    ee, dd, cc, bbb, aaa = [Ku(shifted[index]) for index in range(5)]
    if ee != v0**2:
        raise ArithmeticError("shifted quartic constant is not the distinguished square")
    a1g = dd / v0
    a2g = cc - dd**2 / (4 * v0**2)
    a3g = 2 * v0 * bbb
    a4g = -4 * v0**2 * aaa
    b2g = a1g**2 + 4 * a2g
    gauge_record = direct["weierstrass_model"]["gauge"]
    gauge = Ku(
        Ru([QQ(value) for value in gauge_record["numerator_coefficients_low_to_high"]])
        / Ru([QQ(value) for value in gauge_record["denominator_coefficients_low_to_high"]])
    )
    Achild = Ru([QQ(value) for value in direct["weierstrass_model"]["A_coefficients_low_to_high"]])
    Bchild = Ru([QQ(value) for value in direct["weierstrass_model"]["B_coefficients_low_to_high"]])
    Delta_child = Ru([QQ(value) for value in direct["weierstrass_model"]["discriminant_coefficients_low_to_high"]])
    if Delta_child.gcd(Delta_child.derivative()).degree() != 0:
        raise ArithmeticError("alternate-Q80 discriminant lost squarefreeness")

    def point_on_child(t_section, x_section, y_section):
        slope_section = evaluate_rational(slope_m, t_section)
        xP_section = evaluate_rational(xxP, t_section)
        yP_section = evaluate_rational(yyP, t_section)
        if y_section + yP_section != slope_section * (x_section - xP_section):
            raise ArithmeticError("inherited curve does not lie on its compiled D-fibre")
        radical_root = 2 * x_section - (slope_section**2 - xP_section)
        W_section = radical_root * evaluate_polynomial(
            denominator_sqrt, t_section
        ) / evaluate_polynomial(square_factor, t_section)
        if W_section**2 != evaluate_polynomial(quartic, t_section):
            raise ArithmeticError("inherited curve failed the alternate quartic")
        z_section = t_section - t_section.parent()(t0)
        x_general = (
            2 * t_section.parent()(v0) * (W_section + t_section.parent()(v0))
            + t_section.parent()(dd) * z_section
        ) / z_section**2
        y_general = (
            4 * t_section.parent()(v0) ** 2 * (W_section + t_section.parent()(v0))
            + 2 * t_section.parent()(v0) * t_section.parent()(dd) * z_section
            + (
                2 * t_section.parent()(v0) * t_section.parent()(cc)
                - t_section.parent()(dd) ** 2 / (2 * t_section.parent()(v0))
            ) * z_section**2
        ) / z_section**3
        X_section = t_section.parent()(gauge**2) * 9 * (
            x_general + t_section.parent()(b2g) / 12
        )
        Y_section = t_section.parent()(gauge**3) * 27 * (
            y_general
            + (t_section.parent()(a1g) * x_general + t_section.parent()(a3g)) / 2
        )
        if Y_section**2 != X_section**3 + t_section.parent()(Achild) * X_section + t_section.parent()(Bchild):
            raise ArithmeticError("inherited lifted section failed alternate-Q80 equation")
        return X_section, Y_section

    pinned = load_matrix(PINNED)
    trace_w = vector(ZZ, direct["divisor"]["pinned_trace_vector_w"])
    transport = matrix(ZZ, direct["frame_certificate"]["transport_rows_D_D_plus_O_complement"])
    transport_inverse = transport.inverse()
    alternate_frame = matrix(ZZ, direct["frame_certificate"]["frame_gram"])
    basis_change = matrix(ZZ, target["pinned_identification"]["basis_change_matrix"])

    minim = matrix(ZZ, pari(pinned).qfminim(4)[2])
    candidates = []
    for column in minim.columns():
        for old_mw in (vector(ZZ, column), -vector(ZZ, column)):
            if old_mw * pinned * old_mw != 4:
                continue
            old_degree = 5 - trace_w * pinned * old_mw
            if old_degree != 2:
                continue
            published_mw = vector(ZZ, old_mw * basis_change.transpose())
            old_class = vector(ZZ, [1, 1] + list(old_mw))
            new_class = old_class * transport_inverse
            if new_class[1] != 2 or any(value not in ZZ for value in new_class):
                raise ArithmeticError("inherited curve did not have integral degree two")
            alternate_w = vector(ZZ, new_class[2:])
            if alternate_w * alternate_frame * alternate_w != 10:
                raise ArithmeticError("inherited curve is not a norm-ten alternate bisection")
            mask = sum((int(entry) % 2) << index for index, entry in enumerate(alternate_w))
            candidates.append((tuple(map(int, published_mw)), mask, old_mw, alternate_w))
    candidates.sort(key=lambda item: item[0])
    if len(candidates) != 121:
        raise ArithmeticError(f"expected 121 inherited bisections, obtained {len(candidates)}")
    if len({mask for _, mask, _, _ in candidates}) != 121:
        raise ArithmeticError("the inherited curves do not represent 121 distinct alternate classes")

    seed_vectors = {
        tuple([0, 0, -1, 0, 0, 0, 0, 1] + [0] * 9): "B1",
        tuple([-1, -1, 0, 0, 0, 0, -1] + [0] * 10): "B2",
    }
    seed_expected = {
        "B1": "-P3+P8",
        "B2": "-P1-P2-P7",
    }
    records = []
    seen_seed_labels = set()
    for position, (published_tuple, mask, old_mw, alternate_w) in enumerate(candidates):
        published_mw = vector(ZZ, published_tuple)
        expression = published_expression(published_tuple)
        seed_label = seed_vectors.get(published_tuple)
        if seed_label is not None:
            seen_seed_labels.add(seed_label)
            if expression != seed_expected[seed_label]:
                raise ArithmeticError(f"{seed_label} expression changed")
        label = seed_label or f"inherited-{mask:05x}"
        old_point = sum(
            (
                coefficient * point
                for coefficient, point in zip(published_mw, published_basis)
            ),
            Eold(0),
        )
        Xold, Yold = Kt(old_point[0]), Kt(old_point[1])
        if Xold.denominator() != 1 or Yold.denominator() != 1:
            raise ArithmeticError("old height-four point is not polynomial")
        L0 = a0 * (Xold * h**2 - Nx) + b0 * (Yold * h**3 + Ny)
        L1 = a1 * (Xold * h**2 - Nx) + b1 * (Yold * h**3 + Ny)
        old_base_map = Kt(L1 / L0)
        map_numerator = Rt(old_base_map.numerator())
        map_denominator = Rt(old_base_map.denominator())
        if max(map_numerator.degree(), map_denominator.degree()) != 2:
            raise ArithmeticError("inherited pencil restriction is not degree two")
        relation = Stu(lift_t(map_numerator) - u * lift_t(map_denominator))
        if relation.degree() != 2:
            raise ArithmeticError("inherited quadratic relation lost its leading term")
        leading, linear, constant = [Ku(relation[index]) for index in (2, 1, 0)]
        raw_discriminant = Ru(linear**2 - 4 * leading * constant)
        canonical_q, square_multiplier, constant_class = canonical_squareclass(
            raw_discriminant, Ru
        )
        if canonical_q.degree() != 2 or canonical_q.gcd(canonical_q.derivative()).degree() != 0:
            raise ArithmeticError("inherited cover is not a smooth quadratic cover")
        if canonical_q.gcd(Delta_child).degree() != 0:
            raise ArithmeticError("inherited cover branches over an alternate singular fibre")

        Rs = PolynomialRing(Ku, "S")
        S = Rs.gen()
        quadratic_field = Ku.extension(S**2 - Ku(canonical_q), names="s")
        s = quadratic_field.gen()
        t_section = (
            -quadratic_field(linear) + quadratic_field(square_multiplier) * s
        ) / (2 * quadratic_field(leading))
        if evaluate_rational(old_base_map, t_section) != quadratic_field(u):
            raise ArithmeticError("quadratic root does not invert the inherited base map")
        x_section = evaluate_rational(Xold, t_section)
        y_section = evaluate_rational(Yold, t_section)
        Xnew, Ynew = point_on_child(t_section, x_section, y_section)

        record = {
            "label": label,
            "seed_label": seed_label,
            "published_basis_expression": expression,
            "published_basis_w": list(map(int, published_mw)),
            "pinned_rank17_w": list(map(int, old_mw)),
            "alternate_rank17_w": list(map(int, alternate_w)),
            "lattice_orbit_mask": hex(mask),
            "old_height": 4,
            "alternate_fibre_degree": 2,
            "alternate_bisection_trace_norm": 10,
            "old_base_map_u_of_t": rational_function_record(old_base_map),
            "quadratic_cover": {
                "variable": "t",
                "relation": "a(u)*t^2+b(u)*t+c(u)=0",
                "leading_coefficients": polynomial_text(Ru(leading)),
                "linear_coefficients": polynomial_text(Ru(linear)),
                "constant_coefficients": polynomial_text(Ru(constant)),
                "discriminant_coefficients_low_to_high": polynomial_text(raw_discriminant),
            },
            "canonical_squareclass": {
                "equation": "s^2=q(u)",
                "q_coefficients_low_to_high": polynomial_text(canonical_q),
                "rational_constant_squareclass_signed_squarefree_integer": int(constant_class),
                "raw_discriminant_equals_multiplier_squared_times_q": True,
                "multiplier": rational_function_record(square_multiplier),
            },
            "lifted_section": {
                "coefficient_field": "QQ(u,s), s^2=q(u)",
                "old_parameter_t": quadratic_element_record(t_section, Ku),
                "X": quadratic_element_record(Xnew, Ku),
                "Y": quadratic_element_record(Ynew, Ku),
                "equation_verified": True,
                "deck_conjugation": "s -> -s",
                "anti_invariant_height": 12,
            },
        }
        records.append(record)
        print(
            f"ALT121|{position + 1}/121|{label}|mask={mask:#07x}|qdeg={canonical_q.degree()}",
            flush=True,
        )
    if seen_seed_labels != {"B1", "B2"}:
        raise ArithmeticError("one or both requested seed bisections were not recovered")

    result = {
        "schema": "elkies-k3.bisection-extension-input.v1",
        "artifact_schema": "elkies-k3.r17-norm12-11952-inherited-bisection-covers.v1",
        "status": "PASS_EXACT_121_INHERITED_ALTERNATE_Q80_BISECTION_COVERS",
        "base_parameter": "u",
        "invariant_mw_rank": 17,
        "construction": {
            "source_height_four_pairs": 1311,
            "source_oriented_height_four_sections": 2622,
            "degree_formula": "D.S_v=5-<w,v>",
            "degree_two_oriented_sections": 121,
            "distinct_alternate_translation_classes": 121,
            "alternate_bisection_norm": 10,
            "seed_bisections": seed_expected,
        },
        "bisections": records,
        "inputs": {
            relative(path): digest(path)
            for path in (MODEL, SECTIONS, PINNED, TARGET, SPLITTING, DIRECT)
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "pari_version": ".".join(map(str, pari.version())),
            "required_features": [
                "Sage exact rational function fields",
                "Sage relative quadratic function fields",
                "PARI qfminim",
            ],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "construct_r17_norm12_11952_inherited_bisections.sage"
        ),
        "proof_boundary": (
            "This exact replay enumerates all 121 old height-four curves of degree two "
            "over the alternate-Q80 base, proves that they occupy distinct norm-ten "
            "translation classes, constructs each quadratic cover with its rational "
            "constant squareclass retained, and verifies one lifted section on the "
            "canonical 24I1 equation. Collision and product-character conclusions are "
            "separate exact post-processing steps."
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != serialized:
            raise ArithmeticError("stored 121-cover artifact differs from exact replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        "R17ALTQ80INHERITED|height4=2622-oriented|degree2=121|classes=121|"
        f"output={relative(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
