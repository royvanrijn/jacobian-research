#!/usr/bin/env sage-python
"""Certify the five wgxli fibres in the 43-chart norm-twelve atlas.

The first atlas pass finds one rational-PGL2 class of eight charts whose
j-map contains ICARM curves 351, 356, 376, 377 and 385.  This replay closes
the conditional gates left by that pass:

* it proves that all forty chart/target matches have trivial quadratic twist
  and writes explicit Q-isomorphisms;
* it recovers a unimodularly rebased set of seventeen polynomial sections on
  the representative chart and transports them to all eight native charts;
* it computes the exact height Gram, root count, determinant and integral
  isometry class; and
* it proves independence of every displayed target point by exact mod-2
  reduction certificates, thereby determining the quotient of the displayed
  subgroup by the specialized generic subgroup.

The public source records are hash-pinned.  The default replay reads their
claim-sufficient exact projection from the committed 69-fibre snapshot;
``--live-pinned-source`` retains the original raw-URL audit.  A source change
fails closed.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from urllib.request import urlopen

from sage.all import Matrix, PolynomialRing, QQ, ZZ, matrix, pari
from sage.env import SAGE_VERSION


sys.set_int_max_str_digits(0)

ROOT = Path(__file__).resolve().parents[2]
ATLAS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-record-lineage-atlas-v1.json"
)
LINEAGE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_lineage_v1.json"
)
PINNED_R17 = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
)
PUBLIC_FIBRE_PROJECTION = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
)
PUBLIC_FIBRE_PROJECTION_SHA256 = (
    "9a2675ab48cc37111d1f4050bd1797fc84c98b7839668d292d11406efe7a9eaa"
)
REPRESENTATIVE = "norm12-orbit-074d9"
TARGET_IDS = (351, 356, 376, 377, 385)
HIT_CLASS = (
    "norm12-orbit-074d9",
    "norm12-orbit-08aaa",
    "norm12-orbit-1104c",
    "norm12-orbit-0b4c1",
    "norm12-orbit-05de2",
    "norm12-orbit-10f74",
    "norm12-orbit-10e1c",
    "norm12-orbit-0a2f9",
)
PUBLIC_SOURCES = {
    351: (
        "https://elliptic-rank.icarm.cloud/curve/351.json",
        "f1dc8dcd8f44f7a6ab761833b36eedf2d24175ccdb0fed0268c5a645eee39fe3",
    ),
    356: (
        "https://elliptic-rank.icarm.cloud/curve/356.json",
        "8028758a13d01d4b436618d02d9d2dc7b94860e1769bd9946da2559016a29c39",
    ),
    376: (
        "https://elliptic-rank.icarm.cloud/curve/376.json",
        "eb441bf2acfaaa037c5ce585c83a93c65a0e811ad17ace5367076434937b2a69",
    ),
    377: (
        "https://elliptic-rank.icarm.cloud/curve/377.json",
        "b0c25293198b39fd1495a5e48b1c7087629a1d4404f257e4ed081807dad295f1",
    ),
    385: (
        "https://elliptic-rank.icarm.cloud/curve/385.json",
        "f77a8b5d056f8545a2fa40bfa28f82e3cc633b93978f1deaa1071ce297e0190e",
    ),
}

# Rows are words in the ordered first seventeen displayed points.  The first
# ten are literal columns; the final seven are the shears needed to expose the
# common family.  This matrix has determinant +/-1.
GENERIC_WORDS = (
    ((2, 1),),
    ((3, 1),),
    ((4, 1),),
    ((5, 1),),
    ((8, 1),),
    ((11, 1),),
    ((13, 1),),
    ((15, 1),),
    ((16, 1),),
    ((17, 1),),
    ((1, 1), (2, -1)),
    ((1, 1), (6, -1)),
    ((1, 1), (7, -1)),
    ((1, 1), (9, -1)),
    ((1, 1), (10, -1)),
    ((1, 1), (12, -1)),
    ((1, 1), (14, -1)),
)

sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def python_fraction(value) -> Fraction:
    value = QQ(value)
    return Fraction(int(value.numerator()), int(value.denominator()))


def polynomial_text(poly) -> list[str]:
    if not poly:
        return ["0"]
    return [rational_text(poly[index]) for index in range(poly.degree() + 1)]


def rational_function_record(value) -> dict[str, list[str]]:
    value = value.parent()(value)
    return {
        "numerator_coefficients_low_to_high": polynomial_text(value.numerator()),
        "denominator_coefficients_low_to_high": polynomial_text(value.denominator()),
    }


def load_matrix(path: Path):
    rows = [
        [ZZ(entry) for entry in line.split()]
        for line in path.read_text().splitlines()
        if line.strip()
    ]
    return matrix(ZZ, rows)


def homogeneous_substitution(poly, numerator, denominator, degree):
    answer = numerator.parent()(0)
    for index in range(degree + 1):
        answer += poly[index] * numerator**index * denominator ** (degree - index)
    return answer


def short_invariants(ainvs):
    a1, a2, a3, a4, a6 = map(QQ, ainvs)
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    return {
        "a1": a1,
        "a3": a3,
        "b2": b2,
        "A": -c4 / 48,
        "B": -c6 / 864,
        "canonical_A": 36**2 * (-c4 / 48),
        "canonical_B": 216**2 * (-c6 / 864),
    }


def original_to_canonical(point, invariants):
    x_value, y_value = map(QQ, point)
    return (
        36 * x_value + 3 * invariants["b2"],
        216 * y_value + 108 * (invariants["a1"] * x_value + invariants["a3"]),
    )


def point_negate(point):
    if point is None:
        return None
    return point[0], -point[1]


def point_add(left, right, coefficient_a):
    if left is None:
        return right
    if right is None:
        return left
    x_left, y_left = left
    x_right, y_right = right
    if x_left == x_right:
        if y_left == -y_right:
            return None
        slope = (3 * x_left**2 + coefficient_a) / (2 * y_left)
    else:
        slope = (y_right - y_left) / (x_right - x_left)
    x_answer = slope**2 - x_left - x_right
    return x_answer, -y_left + slope * (x_left - x_answer)


def point_multiply(point, scalar, coefficient_a):
    scalar = int(scalar)
    if scalar < 0:
        return point_multiply(point_negate(point), -scalar, coefficient_a)
    answer = None
    addend = point
    while scalar:
        if scalar & 1:
            answer = point_add(answer, addend, coefficient_a)
        addend = point_add(addend, addend, coefficient_a)
        scalar >>= 1
    return answer


def evaluate_word(points, word, coefficient_a):
    answer = None
    for one_based_index, coefficient in word:
        answer = point_add(
            answer,
            point_multiply(points[one_based_index - 1], coefficient, coefficient_a),
            coefficient_a,
        )
    if answer is None:
        raise ArithmeticError("a declared generic word specialized to zero")
    return answer


def word_matrix():
    rows = []
    for word in GENERIC_WORDS:
        row = [0] * 17
        for one_based_index, coefficient in word:
            row[one_based_index - 1] += coefficient
        rows.append(row)
    return matrix(ZZ, rows)


def lagrange_polynomial(samples, ring):
    variable = ring.gen()
    answer = ring(0)
    for index, (abscissa, ordinate) in enumerate(samples):
        numerator = ring(1)
        denominator = QQ(1)
        for other_index, (other_abscissa, _other_ordinate) in enumerate(samples):
            if other_index == index:
                continue
            numerator *= variable - other_abscissa
            denominator *= abscissa - other_abscissa
        answer += ordinate * numerator / denominator
    return ring(answer)


def x_zero_intersection(x_coordinate) -> int:
    numerator = x_coordinate.numerator()
    denominator = x_coordinate.denominator()
    finite = 0
    for factor, exponent in denominator.factor():
        if exponent % 2:
            raise ArithmeticError("a section x-denominator is not a square")
        finite += factor.degree() * (exponent // 2)
    numerator_degree = numerator.degree()
    infinity = max(0, (numerator_degree - 2 * denominator.degree() - 4 + 1) // 2)
    return int(finite + infinity)


def height_gram(sections, coefficient_a, function_field):
    size = len(sections)
    answer = matrix(QQ, size, size)
    for index in range(size):
        answer[index, index] = 4
        x_left, y_left = map(function_field, sections[index])
        for other in range(index):
            x_right, y_right = map(function_field, sections[other])
            slope = (y_left + y_right) / (x_left - x_right)
            x_difference = slope**2 - x_left - x_right
            difference_height = 4 + 2 * x_zero_intersection(x_difference)
            pairing = (8 - difference_height) / 2
            answer[index, other] = answer[other, index] = pairing
    return matrix(ZZ, answer)


def load_public_records(lineage_fibres, *, live_pinned_source):
    offline_by_id = None
    if not live_pinned_source:
        observed = digest(PUBLIC_FIBRE_PROJECTION)
        if observed != PUBLIC_FIBRE_PROJECTION_SHA256:
            raise ArithmeticError(
                "committed public-fibre projection changed: "
                f"{observed} != {PUBLIC_FIBRE_PROJECTION_SHA256}"
            )
        projection = json.loads(PUBLIC_FIBRE_PROJECTION.read_text())
        if projection.get("status") != (
            "PASS_PINNED_PUBLIC_POINT_PROJECTION_FOR_69_RECOGNIZED_FIBRES"
        ):
            raise ArithmeticError("committed public-fibre projection is not certified")
        offline_by_id = {
            int(record["id"]): record for record in projection.get("records", [])
        }

    records = {}
    for curve_id in TARGET_IDS:
        url, expected_hash = PUBLIC_SOURCES[curve_id]
        if live_pinned_source:
            with urlopen(url, timeout=60) as response:
                raw = response.read()
            observed_hash = hashlib.sha256(raw).hexdigest()
            if observed_hash != expected_hash:
                raise ArithmeticError(
                    f"public curve {curve_id} changed: {observed_hash} != {expected_hash}"
                )
            record = json.loads(raw)
        else:
            record = offline_by_id.get(curve_id)
            if record is None:
                raise ArithmeticError(
                    f"committed public-fibre projection omitted curve {curve_id}"
                )
            record = {
                **record,
                "rank_lower_bound": record["snapshot_rank_lower_bound"],
            }
        lineage = lineage_fibres[curve_id]
        if int(record["id"]) != curve_id:
            raise ArithmeticError("public curve id changed")
        if record["ainvs"] != lineage["source_ainvs"]:
            raise ArithmeticError(f"curve {curve_id} equation changed")
        if int(record["rank_lower_bound"]) != int(lineage["rank_lower_bound"]):
            raise ArithmeticError(f"curve {curve_id} rank lower bound changed")
        if len(record["points"]) != int(record["rank_lower_bound"]):
            raise ArithmeticError(f"curve {curve_id} displayed point count changed")
        records[curve_id] = record
    return records


def signature_record(signature):
    return {
        "prime": int(signature.prime),
        "group_order": int(signature.group_order),
        "doubled_subgroup_order": int(signature.doubled_subgroup_order),
        "quotient_dimension": int(signature.quotient_dimension),
        "rows": [list(map(int, row)) for row in signature.rows],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=ATLAS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--live-pinned-source",
        action="store_true",
        help="require the five original public curve URLs to retain their byte hashes",
    )
    args = parser.parse_args()

    atlas = json.loads(args.atlas.read_text())
    lineage = json.loads(LINEAGE.read_text())
    if atlas["status"] != "RATIONAL_J_MATCH_REQUIRES_TWIST_AND_SECTION_FOLLOWUP":
        raise ArithmeticError("the atlas is not at its expected conditional status")
    if atlas["atlas"]["chart_count"] != 43:
        raise ArithmeticError("the exact atlas no longer has 43 charts")
    if atlas["atlas"]["pgl2_equivalence_class_count"] != 6:
        raise ArithmeticError("the exact atlas no longer has six j-map classes")
    if len(atlas["rational_j_matches"]) != 40:
        raise ArithmeticError("the exact atlas no longer has forty lineage hits")

    target_ids_by_hit = {int(hit["curve_id"]) for hit in atlas["rational_j_matches"]}
    charts_by_hit = {hit["chart"] for hit in atlas["rational_j_matches"]}
    if target_ids_by_hit != set(TARGET_IDS) or charts_by_hit != set(HIT_CLASS):
        raise ArithmeticError("the rational-hit inventory changed")
    if any(
        int(curve_id) in (273, 302) and target["rational_match"]
        for chart in atlas["atlas"]["charts"]
        for curve_id, target in chart["target_j_preimages"].items()
    ):
        raise ArithmeticError("a record curve unexpectedly acquired a rational preimage")

    R = PolynomialRing(QQ, "u")
    u = R.gen()
    K = R.fraction_field()
    chart_records = {record["label"]: record for record in atlas["atlas"]["charts"]}
    if chart_records[REPRESENTATIVE]["frame_class"] != "published-R17":
        raise ArithmeticError("the lineage representative lost its published-R17 frame")

    def chart_polynomial(label, key):
        return R(
            [
                QQ(value)
                for value in chart_records[label]["weierstrass_model"][key]
            ]
        )

    A_representative = chart_polynomial(
        REPRESENTATIVE, "A_coefficients_low_to_high"
    )
    B_representative = chart_polynomial(
        REPRESENTATIVE, "B_coefficients_low_to_high"
    )

    lineage_fibres = {
        int(record["curve_id"]): record
        for record in lineage["rootless_k3_interpolation_input"]["fibres"]
    }
    if tuple(sorted(lineage_fibres)) != TARGET_IDS:
        raise ArithmeticError("lineage target inventory changed")
    public_records = load_public_records(
        lineage_fibres, live_pinned_source=args.live_pinned_source
    )

    invariants = {}
    canonical_points = {}
    word_points = {}
    independence = []
    for curve_id in TARGET_IDS:
        source = public_records[curve_id]
        inv = short_invariants(source["ainvs"])
        invariants[curve_id] = inv
        points = [original_to_canonical(point, inv) for point in source["points"]]
        if any(
            y_value**2
            != x_value**3 + inv["canonical_A"] * x_value + inv["canonical_B"]
            for x_value, y_value in points
        ):
            raise ArithmeticError(f"curve {curve_id} has an off-curve displayed point")
        expected_first = [
            tuple(map(QQ, point))
            for point in lineage_fibres[curve_id]["short_points_first_17"]
        ]
        if points[:17] != expected_first:
            raise ArithmeticError(f"curve {curve_id} first-17 coordinates changed")
        canonical_points[curve_id] = points
        word_points[curve_id] = [
            evaluate_word(points[:17], word, inv["canonical_A"])
            for word in GENERIC_WORDS
        ]

        coefficients = (
            Fraction(0),
            Fraction(0),
            Fraction(0),
            python_fraction(inv["canonical_A"]),
            python_fraction(inv["canonical_B"]),
        )
        fraction_points = [
            (python_fraction(x_value), python_fraction(y_value))
            for x_value, y_value in points
        ]
        no_two_torsion_prime = find_two_torsion_certificate_prime(
            coefficients, prime_bound=200
        )
        signatures = find_mod2_reduction_certificate(
            coefficients, fraction_points, prime_bound=500
        )
        exact_rank = combined_mod2_rank(signatures, len(points))
        if exact_rank != len(points):
            raise ArithmeticError(
                f"curve {curve_id} displayed subgroup independence did not close"
            )
        independence.append(
            {
                "curve_id": curve_id,
                "displayed_point_count": len(points),
                "proved_displayed_subgroup_rank": exact_rank,
                "no_rational_2_torsion_modular_prime": no_two_torsion_prime,
                "mod2_reduction_signatures": [
                    signature_record(signature) for signature in signatures
                ],
                "public_source": PUBLIC_SOURCES[curve_id][0],
                "public_source_sha256": PUBLIC_SOURCES[curve_id][1],
            }
        )

    hits = {
        (hit["chart"], int(hit["curve_id"])): QQ(hit["finite_roots"][0])
        for hit in atlas["rational_j_matches"]
    }
    representative_isomorphisms = {}
    for curve_id in TARGET_IDS:
        parameter = hits[(REPRESENTATIVE, curve_id)]
        A_fibre = A_representative(parameter)
        B_fibre = B_representative(parameter)
        inv = invariants[curve_id]
        q_value = inv["B"] * A_fibre / (B_fibre * inv["A"])
        if not q_value.is_square():
            raise ArithmeticError(f"curve {curve_id} is a nontrivial twist")
        s_value = q_value.sqrt()
        if inv["A"] != q_value**2 * A_fibre:
            raise ArithmeticError(f"curve {curve_id} A-isomorphism identity failed")
        if inv["B"] != q_value**3 * B_fibre:
            raise ArithmeticError(f"curve {curve_id} B-isomorphism identity failed")
        representative_isomorphisms[curve_id] = (q_value, s_value)

    sections = []
    for word_index, word in enumerate(GENERIC_WORDS):
        samples_x = []
        samples_y = []
        for curve_id in TARGET_IDS:
            parameter = hits[(REPRESENTATIVE, curve_id)]
            q_value, s_value = representative_isomorphisms[curve_id]
            x_canonical, y_canonical = word_points[curve_id][word_index]
            samples_x.append((parameter, x_canonical / (36 * q_value)))
            samples_y.append((parameter, y_canonical / (216 * s_value**3)))
        x_section = lagrange_polynomial(samples_x, R)
        if x_section.degree() > 4:
            raise ArithmeticError("interpolated section x-degree exceeds four")
        rhs = x_section**3 + A_representative * x_section + B_representative
        if not rhs.is_square():
            raise ArithmeticError("interpolated section ordinate is not a square")
        y_section = rhs.sqrt()
        if all(y_section(parameter) == ordinate for parameter, ordinate in samples_y):
            pass
        elif all(-y_section(parameter) == ordinate for parameter, ordinate in samples_y):
            y_section = -y_section
        else:
            raise ArithmeticError("no common ordinate sign interpolates the five fibres")
        if y_section.degree() > 6:
            raise ArithmeticError("interpolated section y-degree exceeds six")
        sections.append((x_section, y_section))

    words = word_matrix()
    if abs(words.det()) != 1:
        raise ArithmeticError("generic word matrix is not unimodular")
    gram = height_gram(sections, A_representative, K)
    pinned = load_matrix(PINNED_R17)
    if not gram.is_positive_definite() or gram.det() != 948:
        raise ArithmeticError("generic height lattice has the wrong determinant or signature")
    signed_root_count = int(pari(gram).qfminim(2)[0])
    if signed_root_count:
        raise ArithmeticError("generic height lattice unexpectedly has roots")
    if pari(gram).qfisom(pari(pinned)) == 0:
        raise ArithmeticError("generic height lattice is not the published R17 lattice")

    hit_class_record = next(
        record
        for record in atlas["atlas"]["pgl2_equivalence_classes"]
        if record["representative"] == REPRESENTATIVE
    )
    if tuple(member["label"] for member in hit_class_record["members"]) != HIT_CLASS:
        raise ArithmeticError("the hit-class ordering changed")
    if any(chart_records[label]["frame_class"] != "published-R17" for label in HIT_CLASS):
        raise ArithmeticError("the lineage class is no longer entirely published-R17")

    transported_sections = {}
    family_transports = {}
    for member in hit_class_record["members"]:
        label = member["label"]
        a, b, c, d = map(
            QQ, member["representative_to_member_pgl2_matrix_a_b_c_d"]
        )
        determinant = a * d - b * c
        if not determinant:
            raise ArithmeticError("singular PGL2 transport")
        A_member = chart_polynomial(label, "A_coefficients_low_to_high")
        B_member = chart_polynomial(label, "B_coefficients_low_to_high")
        numerator = a * u + b
        denominator = c * u + d
        A_pullback = homogeneous_substitution(A_member, numerator, denominator, 8)
        B_pullback = homogeneous_substitution(B_member, numerator, denominator, 12)
        q_family = K(B_representative * A_pullback) / K(
            B_pullback * A_representative
        )
        if not q_family.is_square():
            raise ArithmeticError(f"{label} is a nontrivial family twist")
        s_family = q_family.sqrt()
        if (
            K(A_representative) != q_family**2 * K(A_pullback)
            or K(B_representative) != q_family**3 * K(B_pullback)
        ):
            raise ArithmeticError(f"{label} family isomorphism failed")
        if q_family.numerator().degree() or q_family.denominator().degree():
            raise ArithmeticError(f"{label} PGL2 family scale is unexpectedly nonconstant")

        native = []
        inverse_numerator = d * u - b
        inverse_denominator = a - c * u
        for x_representative, y_representative in sections:
            x_native = R(
                homogeneous_substitution(
                    x_representative, inverse_numerator, inverse_denominator, 4
                )
                / (q_family * determinant**4)
            )
            y_native = R(
                homogeneous_substitution(
                    y_representative, inverse_numerator, inverse_denominator, 6
                )
                / (s_family**3 * determinant**6)
            )
            if y_native**2 != x_native**3 + A_member * x_native + B_member:
                raise ArithmeticError(f"{label} transported section identity failed")
            native.append((x_native, y_native))
        transported_sections[label] = native
        family_transports[label] = {
            "representative_to_member_base_map": {
                "a_b_c_d": [rational_text(value) for value in (a, b, c, d)],
                "formula": "u_member=(a*u_representative+b)/(c*u_representative+d)",
            },
            "weighted_pullback_to_representative_scale_q": rational_function_record(
                q_family
            ),
            "weighted_pullback_to_representative_scale_s_with_s_squared_q": rational_function_record(
                s_family
            ),
        }

    isomorphism_records = []
    for label in HIT_CLASS:
        A_member = chart_polynomial(label, "A_coefficients_low_to_high")
        B_member = chart_polynomial(label, "B_coefficients_low_to_high")
        for curve_id in TARGET_IDS:
            parameter = hits[(label, curve_id)]
            A_fibre = A_member(parameter)
            B_fibre = B_member(parameter)
            inv = invariants[curve_id]
            q_value = inv["B"] * A_fibre / (B_fibre * inv["A"])
            if not q_value.is_square():
                raise ArithmeticError(f"{label}/curve-{curve_id} has nontrivial twist")
            s_value = q_value.sqrt()
            if (
                inv["A"] != q_value**2 * A_fibre
                or inv["B"] != q_value**3 * B_fibre
            ):
                raise ArithmeticError(f"{label}/curve-{curve_id} isomorphism failed")

            native = transported_sections[label]
            expected = word_points[curve_id]
            first_y = 216 * s_value**3 * native[0][1](parameter)
            if first_y == -expected[0][1]:
                s_value = -s_value
            elif first_y != expected[0][1]:
                raise ArithmeticError("section orientation cannot be aligned")
            for (x_section, y_section), (x_expected, y_expected) in zip(native, expected):
                if 36 * q_value * x_section(parameter) != x_expected:
                    raise ArithmeticError("specialized transported section x mismatch")
                if 216 * s_value**3 * y_section(parameter) != y_expected:
                    raise ArithmeticError("specialized transported section y mismatch")

            isomorphism_records.append(
                {
                    "chart": label,
                    "curve_id": curve_id,
                    "parameter": rational_text(parameter),
                    "quadratic_twist_class": "trivial_over_Q",
                    "q": rational_text(q_value),
                    "s_with_s_squared_q": rational_text(s_value),
                    "short_model_map": {
                        "x_target_short": "q*x_chart",
                        "y_target_short": "s^3*y_chart",
                    },
                    "original_model_inverse_map": {
                        "x_target_original": (
                            f"q*x_chart-({rational_text(inv['b2'])})/12"
                        ),
                        "y_target_original": (
                            "s^3*y_chart-(a1*x_target_original+a3)/2"
                        ),
                        "a1": rational_text(inv["a1"]),
                        "a3": rational_text(inv["a3"]),
                    },
                    "all_seventeen_specialized_words_verified": True,
                }
            )

    quotient_records = []
    independence_by_id = {record["curve_id"]: record for record in independence}
    for curve_id in TARGET_IDS:
        rank = independence_by_id[curve_id]["proved_displayed_subgroup_rank"]
        quotient_records.append(
            {
                "curve_id": curve_id,
                "ambient_group": "subgroup generated by all displayed public points",
                "specialized_generic_group": (
                    "subgroup generated by the specialization of the seventeen "
                    "certified generic sections"
                ),
                "specialized_generic_group_equals_first_seventeen_displayed_points": True,
                "smith_torsion_invariant_factors": [],
                "free_rank": rank - 17,
                "free_basis_modulo_generic_group": [
                    f"P{index}" for index in range(18, rank + 1)
                ],
                "proof": (
                    "the generic-word matrix on P1,...,P17 is unimodular and the "
                    "stored mod-2 reduction certificate proves all displayed points independent"
                ),
            }
        )

    section_records = []
    for index, ((x_section, y_section), word) in enumerate(
        zip(sections, GENERIC_WORDS), start=1
    ):
        section_records.append(
            {
                "section": index,
                "word_in_first_seventeen_displayed_points": [
                    {"point": point_index, "coefficient": coefficient}
                    for point_index, coefficient in word
                ],
                "representative_x_coefficients_low_to_high": polynomial_text(x_section),
                "representative_y_coefficients_low_to_high": polynomial_text(y_section),
            }
        )

    native_chart_section_records = {}
    for label, native in transported_sections.items():
        native_chart_section_records[label] = [
            {
                "section": index,
                "x_coefficients_low_to_high": polynomial_text(x_section),
                "y_coefficients_low_to_high": polynomial_text(y_section),
            }
            for index, (x_section, y_section) in enumerate(native, start=1)
        ]

    payload = {
        "schema": "elkies-k3.r17-norm12-wgxli-lineage-fibres.v1",
        "status": "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS",
        "atlas_consequences": {
            "complete_chart_count": 43,
            "rational_pgl2_j_map_class_count": 6,
            "class_sizes": [
                len(record["members"])
                for record in atlas["atlas"]["pgl2_equivalence_classes"]
            ],
            "curve_273_rational_preimages": 0,
            "curve_302_rational_preimages": 0,
            "lineage_hit_class": list(HIT_CLASS),
            "lineage_target_ids": list(TARGET_IDS),
            "lineage_chart_target_hit_count": 40,
        },
        "representative": {
            "chart": REPRESENTATIVE,
            "equation": "Y^2=X^3+A(u)*X+B(u)",
            "A_coefficients_low_to_high": polynomial_text(A_representative),
            "B_coefficients_low_to_high": polynomial_text(B_representative),
            "sections": section_records,
        },
        "generic_basis": {
            "word_matrix_rows": [list(map(int, row)) for row in words.rows()],
            "word_matrix_determinant": int(words.det()),
            "height_gram": [list(map(int, row)) for row in gram.rows()],
            "height_gram_determinant": int(gram.det()),
            "signed_norm_two_vector_count": signed_root_count,
            "integrally_isometric_to_published_R17": True,
            "rank": 17,
            "saturated": True,
            "saturation_reason": (
                "the height determinant equals the determinant-948 published-R17 frame "
                "and an exact integral qfisom identifies the lattices"
            ),
        },
        "chart_transports": family_transports,
        "native_chart_sections": native_chart_section_records,
        "target_isomorphisms": isomorphism_records,
        "displayed_point_independence": independence,
        "exceptional_quotients": quotient_records,
        "claim_boundary": {
            "proved": [
                "curves 273 and 302 have no rational fibre parameter in any of the complete 43 certified norm-twelve shared-zero degree-two charts",
                "all five lineage curves are Q-isomorphic, without twisting, to fibres in each of the eight charts of the displayed j-map class",
                "the seventeen displayed polynomial sections form a saturated published-R17 generic basis and specialize through the displayed unimodular words",
                "the quotient of each full displayed target subgroup by the specialized generic subgroup is the displayed free abelian group",
            ],
            "not_proved": [
                "absence of curves 273 or 302 from rootless charts outside the certified 43-member shared-zero norm-twelve atlas",
                "that the displayed target subgroup is the full Mordell-Weil group",
                "an exact rank upper bound for any target curve",
                "a rank-32 specialization",
            ],
        },
        "inputs": {
            relative(args.atlas): digest(args.atlas),
            relative(LINEAGE): digest(LINEAGE),
            relative(PINNED_R17): digest(PINNED_R17),
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": [
                "exact QQ polynomial and rational-function arithmetic",
                "PARI qfminim and qfisom",
                "exact exhaustive E(F_p)/2E(F_p) reduction signatures",
            ],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/certify_r17_norm12_wgxli_lineage_fibres.sage"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != serialized:
            raise ArithmeticError("stored wgxli lineage certificate differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        "R17WGXLILINEAGE|charts=43|j_classes=6|record_misses=273,302|"
        "hit_class=8|lineage_fibres=5|generic_rank=17|det=948|"
        "quotient_ranks=8,12,5,6,12|status=PROVED|output={}".format(
            relative(args.output)
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
