#!/usr/bin/env sage-python
"""Certify the integral MW/character-glue tuple for orbit 0x103b2.

The simultaneous split at t=1/25 has two logically separate outputs.

* On the genus-one double cover ``C: s^2=q(t)``, the pulled-back generic
  R17 lattice and the anti-invariant difference of the two bisection lifts
  have an exact integral C2-character saturation.
* On the single specialized elliptic curve over QQ, the seventeen generic
  points together with the new point generate a primitive rank-18 subgroup.

The first item is an integral Mordell--Weil height lattice and therefore has
an honest discriminant form.  The canonical height pairing on the second item
is real-valued, so it is recorded separately and is never presented as a
Nikulin discriminant form.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

from sage.all import (
    EllipticCurve,
    Genus,
    IntegralLattice,
    PolynomialRing,
    QQ,
    ZZ,
    block_diagonal_matrix,
    gcd,
    identity_matrix,
    lcm,
    matrix,
    pari,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
SPLITTING = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
)
RANK28 = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-rank28-genus-one-bisection-pilot-v1.json"
)
RANK28_FINGERPRINTS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/elkies_2026_rank28_construction_fingerprints_v1.json"
)
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-mw-glue-v1.json"
)

TARGET_LABEL = "norm12-orbit-103b2"
NEUTRAL_LABEL = "norm8-orbit-0f6b1"
FALSE_LOCAL_SURVIVOR_LABEL = "norm8-orbit-0c601"
SHELL_BOUNDS = (8, 10, 12)
SPECIALIZATION_SATURATION_PRIMES = (2, 3, 5, 7, 79813, 239999)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def display(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rows_as_strings(value) -> list[list[str]]:
    return [[str(entry) for entry in row] for row in value.rows()]


def discriminant_form_key(gram) -> dict:
    normal = Genus(gram).discriminant_form().normal_form()
    return {
        "invariants": list(map(int, normal.invariants())),
        "quadratic_gram": rows_as_strings(normal.gram_matrix_quadratic()),
        "value_module": str(normal.value_module_qf()),
    }


def smith_invariants(gram) -> list[int]:
    diagonal = gram.smith_form()[0].diagonal()
    return [abs(int(value)) for value in diagonal if abs(int(value)) > 1]


def load_gram():
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in GRAM.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def polynomial(coefficients, ring):
    return ring([QQ(value) for value in coefficients])


def primitive_polynomial(value):
    denominator = lcm([coefficient.denominator() for coefficient in value])
    coefficients = [ZZ(coefficient * denominator) for coefficient in value]
    content = gcd(coefficients)
    coefficients = [coefficient // content for coefficient in coefficients]
    if coefficients[-1] < 0:
        coefficients = [-coefficient for coefficient in coefficients]
    return value.parent()(coefficients)


def rational_square_root(value):
    value = QQ(value)
    if value < 0:
        return None
    numerator = math.isqrt(int(value.numerator()))
    denominator = math.isqrt(int(value.denominator()))
    if numerator * numerator != value.numerator():
        return None
    if denominator * denominator != value.denominator():
        return None
    return QQ(numerator) / denominator


def visible_lattice(base_gram, trace_coordinates):
    """Return the character sum, glue class, and saturated visible lattice."""

    trace = vector(ZZ, trace_coordinates)
    character = block_diagonal_matrix(2 * base_gram, matrix(ZZ, [[16]]))
    glue = vector(QQ, list(trace / 2) + [QQ(1) / 2])
    if glue * character * glue not in ZZ:
        raise ArithmeticError("the half-sum glue class is not even integral")
    if any(value.denominator() != 1 for value in glue * character):
        raise ArithmeticError("the half-sum glue class is not in the dual lattice")
    basis = identity_matrix(QQ, 18)
    basis[17] = glue
    saturated = (basis * character * basis.transpose()).change_ring(ZZ)
    deck = identity_matrix(ZZ, 18)
    deck[17] = vector(ZZ, list(trace) + [-1])
    if deck * saturated * deck.transpose() != saturated or deck**2 != 1:
        raise ArithmeticError("the deck involution failed on the saturated basis")
    anti_coordinates = vector(ZZ, list(-trace) + [2])
    if anti_coordinates * saturated * anti_coordinates != 16:
        raise ArithmeticError("the primitive anti-invariant line changed")
    if gcd(anti_coordinates) != 1:
        raise ArithmeticError("the anti-invariant generator is not primitive")
    eigensum_basis = identity_matrix(ZZ, 18)
    eigensum_basis[17] = anti_coordinates
    if abs(eigensum_basis.det()) != 2:
        raise ArithmeticError("the character-saturation quotient is not order two")
    return character, glue, saturated, deck, anti_coordinates


def shell_spectrum(gram) -> dict[str, int]:
    return {
        str(bound): int(pari(gram).qfminim(bound)[0]) for bound in SHELL_BOUNDS
    }


def survivor_indices(words: list[str], count: int) -> list[int]:
    values = [int(word, 16) for word in words]
    return [
        index
        for index in range(count)
        if (values[index // 64] >> (index % 64)) & 1
    ]


def specialization_height_record(model_payload, hit, target_cover):
    ring = PolynomialRing(QQ, "t")
    t = ring.gen()
    A = polynomial(model_payload["A_coefficients_low_to_high"], ring)
    B = polynomial(model_payload["B_coefficients_low_to_high"], ring)
    parameter = QQ(hit["t"])
    curve = EllipticCurve(QQ, [A(parameter), B(parameter)])
    generic = [curve(QQ(x), QQ(y)) for x, y in hit["generic_MW17_points"]]
    new = curve(*map(QQ, target_cover["point_on_R17_fibre"]))
    points = generic + [new]
    height = curve.height_pairing_matrix(points, precision=256)
    generic_gram = height[:17, :17]
    cross = matrix(height.base_ring(), 17, 1, [height[index, 17] for index in range(17)])
    projection = generic_gram.solve_right(cross)
    defect = height[17, 17] - (cross.transpose() * projection)[0, 0]

    # Reproduce the stable closest-vector comparison used for the rank-28
    # exceptional directions, now at t=1/25.
    corrections = []
    coset_values = []
    generic_vectors = []
    generic_values = []
    for scale in (10_000, 100_000, 1_000_000):
        rounded_generic = matrix(
            ZZ,
            [
                [ZZ((generic_gram[i, j] * scale).round()) for j in range(17)]
                for i in range(17)
            ],
        )
        generic_minimum = pari(rounded_generic).qfminim()
        generic_vector = [int(generic_minimum[2][row, 0]) for row in range(17)]
        generic_vectors.append(generic_vector)
        generic_values.append(int(generic_minimum[1]))

        rounded_full = matrix(
            ZZ,
            [
                [ZZ((height[i, j] * scale).round()) for j in range(18)]
                for i in range(18)
            ],
        )
        bound = ZZ((height[17, 17] * scale).round()) + 2
        enumeration = pari(rounded_full).qfminim(bound, None, 2)
        candidates = []
        for column in range(int(enumeration[2].ncols())):
            candidate = [int(enumeration[2][row, column]) for row in range(18)]
            if abs(candidate[-1]) != 1:
                continue
            if candidate[-1] < 0:
                candidate = [-value for value in candidate]
            value = int(vector(ZZ, candidate) * rounded_full * vector(ZZ, candidate))
            candidates.append((value, candidate[:-1]))
        value, correction = min(candidates)
        corrections.append(correction)
        coset_values.append(value)
    if len({tuple(value) for value in corrections}) != 1:
        raise ArithmeticError("closest specialization representative is unstable")
    if len({tuple(value) for value in generic_vectors}) != 1:
        raise ArithmeticError("generic first shell is unstable")
    correction = vector(ZZ, corrections[0] + [1])
    generic_vector = vector(ZZ, generic_vectors[0])

    return curve, points, {
        "canonical_height_status": "NUMERICAL_256_BIT_NOT_AN_INTEGRAL_LATTICE",
        "raw_new_point_height": str(height[17, 17]),
        "orthogonal_quotient_height_defect": str(defect),
        "orthogonal_projection_coefficients_in_generic_basis": [
            str(value) for value in projection.column(0)
        ],
        "closest_integral_coset_correction_in_generic_basis": corrections[0],
        "closest_integral_coset_height": str(correction * height * correction),
        "generic_first_shell_vector": generic_vectors[0],
        "generic_first_shell_height": str(
            generic_vector * generic_gram * generic_vector
        ),
        "stable_rounding_scales": [10_000, 100_000, 1_000_000],
        "rounded_coset_minima": [
            str(QQ(value) / scale)
            for value, scale in zip(coset_values, (10_000, 100_000, 1_000_000))
        ],
        "rounded_generic_minima": [
            str(QQ(value) / scale)
            for value, scale in zip(generic_values, (10_000, 100_000, 1_000_000))
        ],
    }


def load_specialization_curve_and_points():
    splitting = json.loads(SPLITTING.read_text())
    model_payload = json.loads(MODEL.read_text())
    hit = splitting["exact_stage"]["simultaneous_split_hits"][0]
    target_cover = next(
        cover for cover in hit["split_covers"] if cover["cover_label"] == TARGET_LABEL
    )
    ring = PolynomialRing(QQ, "t")
    t = ring.gen()
    A = polynomial(model_payload["A_coefficients_low_to_high"], ring)
    B = polynomial(model_payload["B_coefficients_low_to_high"], ring)
    parameter = QQ(hit["t"])
    curve = EllipticCurve(QQ, [A(parameter), B(parameter)])
    points = [curve(QQ(x), QQ(y)) for x, y in hit["generic_MW17_points"]]
    points.append(curve(*map(QQ, target_cover["point_on_R17_fibre"])))
    return curve, points


def isolated_specialization_saturation(prime: int) -> dict:
    """Run one eclib p-saturation in a disposable process.

    eclib retains large reduction tables between successive calls.  The two
    large possible saturation primes therefore run in isolated processes so
    the complete replay remains within the repository's memory envelope.
    """

    curve, points = load_specialization_curve_and_points()
    saturated_points, index, regulator = curve.saturation(
        points, min_prime=prime, max_prime=prime
    )
    if saturated_points != points and index == 1:
        raise ArithmeticError("eclib changed a basis while reporting index one")
    return {
        "prime": prime,
        "saturation_index": int(index),
        "regulator_approx": str(regulator),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--skip-specialization-saturation",
        action="store_true",
        help="skip the approximately two-minute full eclib saturation replay",
    )
    parser.add_argument(
        "--specialization-saturation-prime",
        type=int,
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args()

    if args.specialization_saturation_prime is not None:
        print(json.dumps(isolated_specialization_saturation(
            args.specialization_saturation_prime
        ), sort_keys=True))
        return

    splitting = json.loads(SPLITTING.read_text())
    rank28 = json.loads(RANK28.read_text())
    rank28_fingerprints = json.loads(RANK28_FINGERPRINTS.read_text())
    model = json.loads(MODEL.read_text())
    base_gram = load_gram()
    records = splitting["construction"]["records"]
    by_label = {record["label"]: record for record in records}
    target = by_label[TARGET_LABEL]
    hit = splitting["exact_stage"]["simultaneous_split_hits"][0]
    target_cover = next(
        cover for cover in hit["split_covers"] if cover["cover_label"] == TARGET_LABEL
    )
    neutral_cover = next(
        cover for cover in hit["split_covers"] if cover["cover_label"] == NEUTRAL_LABEL
    )
    if hit["t"] != "1/25" or target["lattice_orbit_mask"] != 0x103B2:
        raise ArithmeticError("the pinned simultaneous hit changed")

    trace = vector(ZZ, target["pinned_rank17_w"])
    if trace * base_gram * trace != 12:
        raise ArithmeticError("the target trace no longer has norm twelve")
    character, glue, saturated, deck, anti = visible_lattice(base_gram, trace)
    if saturated.det() != 497_025_024:
        raise ArithmeticError("the saturated visible determinant changed")
    if IntegralLattice(saturated).minimum() != 8:
        raise ArithmeticError("the visible lattice minimum changed")
    if int(pari(saturated).qfminim(2)[0]) != 0:
        raise ArithmeticError("the visible lattice acquired roots")

    character_smith, smith_left, unused = character.smith_form()
    glue_dual_numerator = vector(ZZ, glue * character)
    glue_smith = vector(ZZ, smith_left * glue_dual_numerator)
    glue_smith_coordinates = [
        {
            "index_zero_based": index,
            "coordinate": int(glue_smith[index] % abs(character_smith[index, index])),
            "modulus": abs(int(character_smith[index, index])),
        }
        for index in range(18)
        if abs(character_smith[index, index]) > 1
        and glue_smith[index] % abs(character_smith[index, index])
    ]

    target_spectrum = shell_spectrum(saturated)
    census_spectra = []
    for record in records:
        unused_character, unused_glue, visible, unused_deck, unused_anti = visible_lattice(
            base_gram, record["pinned_rank17_w"]
        )
        census_spectra.append(
            {
                "label": record["label"],
                "trace_norm": int(record["trace_norm"]),
                "signed_vector_counts_through_norm": shell_spectrum(visible),
            }
        )
    same_spectrum = [
        row["label"]
        for row in census_spectra
        if row["signed_vector_counts_through_norm"] == target_spectrum
    ]
    if same_spectrum != [TARGET_LABEL]:
        raise ArithmeticError("the norm-through-12 integral fingerprint is no longer unique")

    rank28_trace = rank28["traces"][0]
    unused_character, unused_glue, rank28_visible, unused_deck, unused_anti = visible_lattice(
        base_gram, rank28_trace["pinned_rank17_w"]
    )
    if pari(saturated).qfisom(pari(rank28_visible)) != 0:
        raise ArithmeticError("the target unexpectedly became isometric to the rank-28 trace lattice")

    ring = PolynomialRing(QQ, "t")
    t = ring.gen()
    q = polynomial(target["branch_polynomial_q_coefficients_low_to_high"], ring)
    h = polynomial(target["trace_section"]["h_coefficients_low_to_high"], ring)
    q_primitive = primitive_polynomial(q)
    h_primitive = primitive_polynomial(h)
    parameter = QQ(1) / 25
    q_value = q(parameter)
    q_root = rational_square_root(q_value)
    if q_root is None or str(q_root) != target_cover["cover_coordinate_s"]:
        raise ArithmeticError("the target quartic no longer splits at t=1/25")

    survivor_index_set = survivor_indices(
        hit["survivor_words_hex_low_to_high"], len(records)
    )
    survivor_labels = [records[index]["label"] for index in survivor_index_set]
    if survivor_labels != [FALSE_LOCAL_SURVIVOR_LABEL, NEUTRAL_LABEL, TARGET_LABEL]:
        raise ArithmeticError("the all-block local survivor set changed")
    exact_square_labels = []
    for record in records:
        branch = polynomial(record["branch_polynomial_q_coefficients_low_to_high"], ring)
        if rational_square_root(branch(parameter)) is not None:
            exact_square_labels.append(record["label"])
    if exact_square_labels != [NEUTRAL_LABEL, TARGET_LABEL]:
        raise ArithmeticError("the exact splitting set at t=1/25 changed")

    false_survivor = by_label[FALSE_LOCAL_SURVIVOR_LABEL]
    false_q = polynomial(
        false_survivor["branch_polynomial_q_coefficients_low_to_high"], ring
    )(parameter)
    if rational_square_root(false_q) is not None:
        raise ArithmeticError("the false local survivor became a rational square")

    curve, specialization_points, height_record = specialization_height_record(
        model, hit, target_cover
    )
    saturation_record = {
        "status": "SKIPPED_BY_EXPLICIT_FLAG",
        "rank_upper_bound_not_computed": True,
    }
    if not args.skip_specialization_saturation:
        prime_records = []
        for prime in SPECIALIZATION_SATURATION_PRIMES:
            completed = subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "--specialization-saturation-prime",
                    str(prime),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
            record = json.loads(completed.stdout.strip().splitlines()[-1])
            if record["saturation_index"] != 1:
                raise ArithmeticError(
                    f"the specialized subgroup is not {prime}-saturated"
                )
            prime_records.append(record)
        saturation_record = {
            "status": "PASS_FULL_ECLIB_SATURATION_OF_DISPLAYED_RANK18_SUBGROUP",
            "displayed_rank": 18,
            "saturation_index": 1,
            "eclib_index_bound": 574651205610,
            "index_bound_factorization": "2 * 3 * 5 * 79813 * 239999",
            "tamagawa_index_primes": [2, 3, 5, 7],
            "tested_saturation_primes": prime_records,
            "rank_upper_bound_not_computed": True,
            "meaning": (
                "The displayed rank-18 subgroup is primitive in E(Q). This does not prove "
                "that E(Q) has rank exactly 18."
            ),
        }
    elif (args.output if args.output.is_absolute() else ROOT / args.output).exists():
        pinned_output = args.output if args.output.is_absolute() else ROOT / args.output
        saturation_record = json.loads(pinned_output.read_text())[
            "specialized_MW_subgroup"
        ]
        if saturation_record.get("status") != (
            "PASS_FULL_ECLIB_SATURATION_OF_DISPLAYED_RANK18_SUBGROUP"
        ):
            raise ArithmeticError("the pinned specialization saturation is incomplete")

    rank28_raw = [
        Decimal(record["relative_height"]["raw_canonical_height"])
        for record in rank28_fingerprints["records"]
    ]
    rank28_defects = [
        Decimal(record["relative_height"]["quotient_height_defect"])
        for record in rank28_fingerprints["records"]
    ]

    quotient_profiles = {
        int(profile["relation_prime"]): profile
        for profile in hit["quotient_by_generic_MW17"]["relation_prime_profiles"]
        if profile["status"] == "PASS_EXACT_FINITE_QUOTIENT_PROFILE"
    }

    saturated_q = discriminant_form_key(saturated)
    result_status = "PASS_EXACT_VISIBLE_MW_GLUE_AND_SPECIALIZATION_RANK_AT_LEAST_18"
    result = {
        "schema": "elkies-k3.r17-norm12-103b2-mw-glue.v1",
        "status": result_status,
        "inputs": {
            display(path): digest(path)
            for path in (SPLITTING, RANK28, RANK28_FINGERPRINTS, MODEL, GRAM, Path(__file__))
        },
        "parameter_specialization": {
            "t": "1/25",
            "point": target_cover["point_on_R17_fibre"],
            "cover_coordinate_s": target_cover["cover_coordinate_s"],
            "exact_q_value": str(q_value),
            "exact_q_square_root": str(q_root),
            "neutral_split_cover": NEUTRAL_LABEL,
            "neutral_point_is_specialized_generic_section": "P16",
            "new_direction_cover": TARGET_LABEL,
        },
        "trace_class": {
            "label": TARGET_LABEL,
            "orbit_mask": int(target["lattice_orbit_mask"]),
            "orbit_hex": "0x103b2",
            "minimum_norm_in_R17_mod_2R17": 12,
            "pinned_R17_coordinates": target["pinned_rank17_w"],
            "published_R17_coordinates": target["published_basis_w"],
            "equation_complexity": target["equation_complexity"],
            "unique_regular_member": True,
        },
        "tuple_order": ["L", "G", "L_chi", "A_L", "q_L", "glue_subgroups"],
        "tuple": {
            "L": {
                "name": "orthogonal_character_sublattice_R17(2)_plus_<16>",
                "rank": 18,
                "determinant_absolute": int(character.det()),
                "gram": "R17(2) direct_sum <16>",
            },
            "G": {
                "name": "quadratic_base_change_deck",
                "structure": "C2",
                "order": 2,
                "action_on_saturated_basis_P1_through_P17_R": rows_as_strings(deck),
            },
            "L_chi": [
                {
                    "character": "+",
                    "basis": "pulled P1,...,P17",
                    "rank": 17,
                    "gram": "R17(2)",
                    "determinant_absolute": int((2 * base_gram).det()),
                },
                {
                    "character": "-",
                    "basis": ["T=R-sigma(R)=2R-tau"],
                    "coordinates_in_saturated_basis": list(map(int, anti)),
                    "rank": 1,
                    "gram": [[16]],
                    "primitive_in_declared_rational_span": True,
                },
            ],
            "A_L": smith_invariants(character),
            "q_L": discriminant_form_key(character),
            "glue_subgroups": [
                {
                    "context": "full_integral_saturation_in_the_declared_18_dimensional_character_span",
                    "subgroup": "<(tau/2,T/2)>",
                    "structure": "Z/2",
                    "order": 2,
                    "index": 2,
                    "trace_coset": "0x103b2 in R17/2R17",
                    "glue_vector_in_character_basis": [str(value) for value in glue],
                    "smith_coordinates_nonzero": glue_smith_coordinates,
                    "quadratic_values_mod_2Z": {
                        "invariant_projection": "6 = 0 mod 2Z",
                        "anti_invariant_projection": "4 = 0 mod 2Z",
                        "total": "10 = 0 mod 2Z",
                    },
                    "known_relation": "2R=tau+T",
                    "exhaustion_argument": (
                        "L+ and L- are primitive intersections with the two rational eigenspaces. "
                        "Any two nontrivial glue classes would differ in one primitive eigensublattice; "
                        "hence the exhibited order-two graph is the full saturation in this span."
                    ),
                }
            ],
        },
        "saturated_visible_MW_lattice": {
            "basis": ["pulled P1", "...", "pulled P17", "R"],
            "gram": rows_as_strings(saturated),
            "rank": 18,
            "determinant_absolute": int(saturated.det()),
            "smith_invariants": smith_invariants(saturated),
            "discriminant_form": saturated_q,
            "minimum": int(IntegralLattice(saturated).minimum()),
            "norm_two_root_count_signed": int(pari(saturated).qfminim(2)[0]),
            "signed_vector_counts_through_norm": target_spectrum,
            "lift_height": int(saturated[17, 17]),
            "lift_zero_intersection": 1,
            "anti_invariant_height": 16,
            "deck_character_saturation_index": 2,
            "scope": (
                "This is the complete integral saturation of the exhibited invariant rank-17 "
                "and anti-invariant rank-one rational span. Extra anti-invariant MW directions "
                "on the cover are not excluded."
            ),
        },
        "local_component_data": {
            "trace_pole_polynomial_degree": int(h.degree()),
            "trace_pole_polynomial_primitive_coefficients_low_to_high": list(
                map(int, h_primitive)
            ),
            "trace_pole_polynomial_irreducible_group": str(
                h_primitive.galois_group(pari_group=True)
            ),
            "branch_polynomial_degree": int(q.degree()),
            "branch_polynomial_primitive_coefficients_low_to_high": list(
                map(int, q_primitive)
            ),
            "branch_polynomial_discriminant": str(q_primitive.discriminant()),
            "branch_polynomial_discriminant_factorization": str(
                q_primitive.discriminant().factor()
            ),
            "branch_polynomial_galois_group": str(
                q_primitive.galois_group(pari_group=True)
            ),
            "branch_points": "four simple conjugate points on smooth source fibres",
            "branch_avoids_trace_poles": True,
            "branch_avoids_24I1_discriminant": True,
            "pulled_fibre_configuration": "48I1; no reducible-fibre component groups",
            "component_corrections": "zero (all fibres are irreducible I1)",
            "conjugate_lift_intersection": 4,
            "lift_self_intersection": -4,
            "height_identities": [
                "height(tau after degree-two pullback)=2*12=24",
                "height(T)=2*(R.sigma(R)-R^2)=2*(4-(-4))=16",
                "height(R)=(24+16)/4=10=8+2*(R.O), hence R.O=1",
            ],
        },
        "bridge_position": {
            "cover_level_integral_bridge": {
                "character_sum": "R17(2) direct_sum <16>",
                "glue_graph_order": 2,
                "glue_graph_generator": "(0x103b2/2,T/2)",
                "saturated_lift_R_height": 10,
            },
            "specialization_quotient_bridge": {
                "point_sequence": (
                    "specialized generic P1,...,P17; norm8-orbit-0f6b1; "
                    "norm12-orbit-103b2"
                ),
                "new_pivot_column_zero_based_mod_2": [
                    value
                    for value in quotient_profiles[2][
                        "finite_quotient_certificate"
                    ]["pivot_columns_zero_based"]
                    if value >= 17
                ],
                "new_pivot_column_zero_based_mod_3": [
                    value
                    for value in quotient_profiles[3][
                        "finite_quotient_certificate"
                    ]["pivot_columns_zero_based"]
                    if value >= 17
                ],
                "norm8_column_zero_based": 17,
                "norm12_column_zero_based": 18,
                "displayed_rank18_subgroup_saturation_index": saturation_record.get(
                    "saturation_index"
                ),
                "rank_upper_bound_not_computed": True,
            },
            "numerical_height_position": height_record,
        },
        "comparison_with_other_142_selected_covers_at_t_1_over_25": {
            "selected_cover_count": len(records),
            "other_cover_count": len(records) - 1,
            "other_covers_not_producing_a_new_quotient_direction": len(records) - 1,
            "nonsplit_other_covers": len(records) - len(exact_square_labels),
            "split_but_generic_dependent_other_covers": 1,
            "exact_square_labels": exact_square_labels,
            "all_twelve_sieve_prime_survivors": survivor_labels,
            "false_all_local_blocks_survivor": {
                "label": FALSE_LOCAL_SURVIVOR_LABEL,
                "q_value": str(false_q),
                "numerator_factorization": str(false_q.numerator().factor()),
                "denominator_factorization": str(false_q.denominator().factor()),
                "global_square": False,
            },
            "integral_fingerprint_unique_through_norm_12": True,
            "target_signed_shell_counts": target_spectrum,
            "distinct_shell_fingerprints_among_143": len(
                {
                    tuple(row["signed_vector_counts_through_norm"].items())
                    for row in census_spectra
                }
            ),
            "interpretation": (
                "The character dimensions, determinant and finite quadratic-form genus do not "
                "cause the arithmetic hit: every compiled genus-one bisection has the same basic "
                "C2 half-sum architecture. The decisive event at this fibre is q(1/25) being a "
                "rational square. The trace coset does, however, select a distinct integral "
                "visible lattice, uniquely fingerprinted in this 143-cover batch through norm 12."
            ),
        },
        "comparison_with_rank28_exceptional_directions": {
            "rank28_parameter": rank28["parameter"],
            "direction_count": int(rank28["target_count"]),
            "shared_trace_for_all_eleven": {
                "orbit_hex": rank28_trace["orbit_hex"],
                "minimum_norm": 8,
                "published_R17_coordinates": rank28_trace["published_basis_w"],
                "pinned_R17_coordinates": rank28_trace["pinned_rank17_w"],
            },
            "mask_xor_with_0x103b2": hex(
                int(rank28_trace["orbit_mask"]) ^ int(target["lattice_orbit_mask"])
            ),
            "mask_hamming_distance": (
                int(rank28_trace["orbit_mask"]) ^ int(target["lattice_orbit_mask"])
            ).bit_count(),
            "same_character_sublattice": "R17(2) direct_sum <16>",
            "same_character_glue_structure": "one isotropic Z/2 graph",
            "same_saturated_determinant_and_discriminant_form": True,
            "integrally_isometric_visible_lattices": False,
            "rank28_trace_signed_shell_counts": shell_spectrum(rank28_visible),
            "rank28_lift_height": 8,
            "target_lift_height": 10,
            "reason_for_lift_height_difference": (
                "The anti line has height 16 in both cases, but the pulled trace has height "
                "16 for norm eight and 24 for norm twelve. Thus R=(tau+T)/2 has height 8 "
                "for the rank-28 template and height 10 for 0x103b2."
            ),
            "rank28_raw_height_range": [str(min(rank28_raw)), str(max(rank28_raw))],
            "rank28_quotient_defect_range": [
                str(min(rank28_defects)),
                str(max(rank28_defects)),
            ],
            "target_specialization_height_position": height_record,
            "rank28_local_component_data": (
                "The eleven covers also have 48I1 and trivial component corrections; "
                "their norm-eight lifts have R.O=0 instead of R.O=1."
            ),
            "interpretation": (
                "The eleven rank-28 directions occupy eleven arithmetic quotient directions but "
                "all use one norm-eight trace coset and eleven fitted pencil members. Orbit "
                "0x103b2 instead uses a rigid norm-twelve member. It sits much closer to the "
                "specialized generic real span, and its shortest integral coset representative "
                "requires a nonzero generic correction."
            ),
        },
        "specialized_MW_subgroup": saturation_record,
        "claim_boundary": {
            "exact": (
                "The displayed cover-level rank-18 character span, its full integral saturation "
                "inside that rational span, Gram, Smith group, finite quadratic form, local fibre "
                "data, and t=1/25 square identity are exact. The displayed specialized rank-18 "
                "subgroup is primitive when the full saturation replay is enabled."
            ),
            "not_claimed": (
                "No rank upper bound is computed for the specialized elliptic curve, and no "
                "upper bound is computed for additional anti-invariant MW directions on the "
                "genus-one cover. The canonical-height bridge comparison is numerical."
            ),
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/certify_r17_norm12_103b2_mw_glue.sage"
        ),
    }

    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = args.output if args.output.is_absolute() else ROOT / args.output
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise SystemExit(f"stale or missing artifact: {output}")
        print("PASS R17 norm12 orbit 103b2 MW/glue certificate")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(serialized)
    print(output.relative_to(ROOT))


if __name__ == "__main__":
    main()
