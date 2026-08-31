#!/usr/bin/env sage-python
"""Search the promoted rank-nine paired base for split bisections.

status: ACTIVE_SEARCH
claim: bounded exact split-bisection search with finite-quotient rank lower bounds
inputs: promoted rank-nine base, published R17 model/sections, complete bisection batch
outputs: artifacts/generated-results/elkies-2026-rank9-paired-base-search.json
supersedes: none; extends verify_elkies_2026_rank19_rank9_base.sage

The search order is deliberately arithmetic:

1. replay PARI's rank interval and bounded small-prime saturation;
2. compute a high-precision canonical-height Gram and LLL-reduce it;
3. enumerate the resulting Mordell--Weil lattice by canonical height, but
   retain candidates primarily by the projective height of the rational
   function ``t(P)``;
4. specialize the seventeen published sections and the two sections defining
   the paired base;
5. test every other bisection quadratic for a rational square, materialize
   every resulting point, and certify the displayed rank by exact finite
   quotients;
6. only after the split count and rank certificate are known, attach a
   three-block Nagao score.

All height bounds and finalist counts are explicit bounded-search parameters.
Failure to find an extra split is not a rank upper bound.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from time import perf_counter

from sage.all import EllipticCurve, PolynomialRing, QQ, RealField, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
BISECTIONS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
BASE = ROOT / "artifacts/generated-results/elkies-2026-rank19-rank9-paired-base.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-rank9-paired-base-search.json"
FINITE_HELPER = ROOT / "elliptic-curves/cas/elliptic_candidate_record.py"
NAGAO_HELPER = ROOT / "elkies-k3/scripts/search_h92_q12o5867_rootless_nagao.py"
TARGET_MASKS = (42110, 43109)
THIRD_MODEL = (
    1,
    0,
    0,
    -70087047578007713577216,
    3865770423647395544516350651140096,
)
PAIRED_MODEL = (
    1,
    0,
    0,
    -60729194722297004073216,
    5758259762216167074332597509226496,
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def point_text(point) -> list[str]:
    if point.is_zero():
        return ["0"]
    return [rational_text(point[0]), rational_text(point[1])]


def evaluate(coefficients, value):
    answer = QQ(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + QQ(coefficient)
    return answer


def reconstruct_sections(section_document, t_value, curve):
    points = []
    for expected, record in enumerate(section_document["sections"]):
        if int(record["basis_index"]) != expected:
            raise ArithmeticError("the published section order changed")
        x_value = evaluate(record["x_coefficients_low_to_high"], t_value)
        if expected == 0:
            y_value = evaluate(record["y_coefficients_low_to_high"], t_value)
        else:
            chord = record["chord"]
            reference = points[int(chord["reference_basis_index"])]
            slope = evaluate(chord["slope_coefficients_low_to_high"], t_value)
            y_value = reference[1] + slope * (x_value - reference[0])
        point = curve(x_value, y_value)
        points.append(point)
    if len(points) != 17:
        raise ArithmeticError("expected seventeen published sections")
    return points


def lifted_point(record, t_value, square_root, curve):
    lifted = record["lifted_section"]
    x_value = evaluate(lifted["x0_coefficients"], t_value) + evaluate(
        lifted["x1_coefficients"], t_value
    ) * square_root
    y_value = evaluate(lifted["y0_coefficients"], t_value) + evaluate(
        lifted["y1_coefficients"], t_value
    ) * square_root
    return curve(x_value, y_value)


def linear_combination(curve, basis, coefficients):
    return sum((ZZ(coefficient) * point for coefficient, point in zip(coefficients, basis)), curve(0))


def height_matrix(curve, points, digits=100):
    gp = shutil.which("gp")
    if gp is None:
        raise FileNotFoundError("PARI/GP is required for the canonical-height Gram")
    coefficients = ",".join(rational_text(value) for value in curve.a_invariants())
    point_vector = ",".join(
        f"[{rational_text(point[0])},{rational_text(point[1])}]" for point in points
    )
    program = "\n".join(
        (
            f"default(realprecision,{digits});",
            f"E=ellinit([{coefficients}]);P=[{point_vector}];H=ellheightmatrix(E,P);",
            f'for(i=1,{len(points)},for(j=1,{len(points)},print("HEIGHT|",i,"|",j,"|",H[i,j])));',
            "quit",
        )
    ) + "\n"
    completed = subprocess.run(
        [gp, "-q", "-f"], input=program, text=True, capture_output=True, check=True
    )
    field = RealField(int((digits + 20) * 3.5))
    entries = {}
    for line in completed.stdout.splitlines():
        fields = line.split("|", 3)
        if len(fields) == 4 and fields[0] == "HEIGHT":
            entries[int(fields[1]) - 1, int(fields[2]) - 1] = field(fields[3])
    if len(entries) != len(points) ** 2:
        raise ArithmeticError("PARI emitted an incomplete canonical-height Gram")
    return matrix(field, len(points), len(points), entries)


def matrix_text(values, digits=90):
    return [[f"{value:.{digits}}" for value in row] for row in values.rows()]


def coordinate_bits(point) -> int:
    if point.is_zero():
        return 0
    return max(
        abs(ZZ(value.numerator())).nbits() + ZZ(value.denominator()).nbits()
        for value in (point[0], point[1])
    )


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


def fraction_point(point):
    return to_fraction(point[0]), to_fraction(point[1])


def to_fraction(value):
    value = QQ(value)
    return Fraction(int(value.numerator()), int(value.denominator()))


def best_finite_certificate(helper, model, points, prime_bound):
    fraction_points = tuple(fraction_point(point) for point in points)
    candidates = []
    for relation_prime in (2, 3):
        certificate = helper.build_finite_quotient_certificate(
            model,
            fraction_points,
            relation_prime=relation_prime,
            prime_bound=prime_bound,
        )
        candidates.append(certificate)
        if certificate["certified_independent"]:
            break
    certificate = max(
        candidates,
        key=lambda item: (
            int(item["combined_rank_over_relation_field"]),
            bool(item["torsion_witness"]),
            -int(item["relation_prime"]),
        ),
    )
    rank = int(certificate["combined_rank_over_relation_field"])
    pivots = tuple(int(value) for value in certificate["pivot_columns_zero_based"])
    if rank and certificate["torsion_witness"] is not None:
        independent_points = tuple(fraction_points[index] for index in pivots)
        replay = helper.build_finite_quotient_certificate(
            model,
            independent_points,
            relation_prime=int(certificate["relation_prime"]),
            prime_bound=prime_bound,
        )
        helper.verify_finite_quotient_certificate(model, independent_points, replay)
        if not replay["certified_independent"]:
            raise ArithmeticError("the pivot subset did not replay independently")
    else:
        replay = None
    return {
        "attempts": candidates,
        "best_relation_prime": int(certificate["relation_prime"]),
        "certified_rank_lower_bound": rank if replay is not None else 0,
        "independent_point_indices_zero_based": list(pivots) if replay is not None else [],
        "independent_subset_certificate": replay,
    }


def paired_cover_maps(third_curve, paired_curve):
    isogenies = [
        morphism
        for morphism in third_curve.isogenies_prime_degree(2)
        if morphism.codomain().is_isomorphic(paired_curve)
    ]
    if len(isogenies) != 1:
        raise ArithmeticError("the expected degree-two paired-base isogeny is not unique")
    isogeny = isogenies[0]
    codomain_isomorphism = (
        None if isogeny.codomain() == paired_curve else isogeny.codomain().isomorphism_to(paired_curve)
    )

    parameter_ring = PolynomialRing(QQ, "s")
    s = parameter_ring.gen()
    reciprocal_coefficients = (
        QQ(1346816601),
        QQ(-962433973020),
        QQ(-309051947898044),
        QQ(168863136988245440),
        QQ(46344697121074403584),
    )
    q, d, c, b, a = reciprocal_coefficients
    if not q.is_square() or q.sqrt() != 36699:
        raise ArithmeticError("the reciprocal pointed quartic changed")
    q = q.sqrt()
    pointed = EllipticCurve(
        QQ,
        [d / q, c - d**2 / (4 * q**2), 2 * q * b, -4 * q**2 * a, a * (d**2 - 4 * q**2 * c)],
    )
    paired_to_pointed = paired_curve.isomorphism_to(pointed)

    def image(point):
        answer = isogeny(point)
        return codomain_isomorphism(answer) if codomain_isomorphism is not None else answer

    def to_cover(point):
        if point.is_zero():
            return None
        X, Y = paired_to_pointed(point)[:2]
        if Y == 0:
            return None
        s_value = (4 * q**2 * (X + c) - d**2) / (2 * q * Y)
        if s_value == 0:
            return None
        reciprocal_y = (X * s_value**2 - d * s_value) / (2 * q) - q
        if reciprocal_y**2 != evaluate(reciprocal_coefficients, s_value):
            raise ArithmeticError("the inverse pointed-quartic identity failed")
        r_value = 1 / s_value
        denominator = 130 * r_value - 38636
        if denominator == 0:
            return None
        t_value = (289444 - r_value**2) / denominator
        u_value = 65 * t_value + r_value
        v_value = reciprocal_y * r_value**2 / denominator
        if u_value**2 != 4225 * t_value**2 + 38636 * t_value + 289444:
            raise ArithmeticError("the first paired conic identity failed")
        if v_value**2 != 1346816601 * t_value**2 + 7403338254 * t_value + 10921221529:
            raise ArithmeticError("the second paired conic identity failed")
        return t_value, u_value, v_value

    return image, to_cover


def enumerate_base_points(curve, basis, gram, integral_gram, height_bound, scale, to_paired, to_cover):
    pari.set_real_precision(300)
    # Include a 1e-6 height guard around the rounded integral form, then
    # enforce the declared cutoff again against the unrounded Gram.
    output = pari.qfminim(integral_gram, ZZ(height_bound) * scale + 100)
    representatives = output[2]
    candidates = {}
    evaluated = 0
    exceptional = 0
    for column in range(representatives.ncols()):
        base_vector = vector(ZZ, [ZZ(representatives[row, column]) for row in range(len(basis))])
        for sign in (1, -1):
            coefficients = sign * base_vector
            height = coefficients * gram * coefficients
            if height > height_bound + gram.base_ring()("1e-40"):
                continue
            point = linear_combination(curve, basis, coefficients)
            evaluated += 1
            cover = to_cover(to_paired(point))
            if cover is None:
                exceptional += 1
                continue
            t_value, u_value, v_value = cover
            numerator_bits = abs(ZZ(t_value.numerator())).nbits()
            denominator_bits = ZZ(t_value.denominator()).nbits()
            record = {
                "lll_coefficients": [int(value) for value in coefficients],
                "canonical_height": f"{height:.60}",
                "third_quotient_point": point_text(point),
                "t": rational_text(t_value),
                "u": rational_text(u_value),
                "v": rational_text(v_value),
                "t_numerator_bits": numerator_bits,
                "t_denominator_bits": denominator_bits,
                "t_projective_bits": max(numerator_bits, denominator_bits),
                "t_total_bits": numerator_bits + denominator_bits,
            }
            key = rational_text(t_value)
            priority = (
                record["t_projective_bits"],
                record["t_total_bits"],
                height,
                record["lll_coefficients"],
            )
            if key not in candidates or priority < candidates[key][0]:
                candidates[key] = (priority, record)
    records = [value[1] for value in candidates.values()]
    records.sort(
        key=lambda item: (
            item["t_projective_bits"],
            item["t_total_bits"],
            float(item["canonical_height"]),
            item["lll_coefficients"],
        )
    )
    best_bits = None
    pareto = []
    for record in sorted(records, key=lambda item: (float(item["canonical_height"]), item["t_projective_bits"])):
        if best_bits is None or record["t_projective_bits"] < best_bits:
            pareto.append(record)
            best_bits = record["t_projective_bits"]
    return records, pareto, {
        "qfminim_signed_vector_count": int(output[0]),
        "evaluated_signed_points": evaluated,
        "distinct_finite_t_values": len(records),
        "exceptional_or_infinite_map_points": exceptional,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=MODEL)
    parser.add_argument("--sections", type=Path, default=SECTIONS)
    parser.add_argument("--bisections", type=Path, default=BISECTIONS)
    parser.add_argument("--base", type=Path, default=BASE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--saturation-bound", type=int, default=1000)
    parser.add_argument("--height-bound", type=int, default=60)
    parser.add_argument("--t-finalists", type=int, default=100)
    parser.add_argument("--certificate-prime-bound", type=int, default=300)
    parser.add_argument(
        "--nagao-primes",
        default="19,41,43,61,71,73,79,83;89,107,113,127,131,137,139,151;157,163,167,173,179,181,191,193,197",
    )
    args = parser.parse_args()
    if args.saturation_bound < 3 or args.height_bound < 1 or args.t_finalists < 1:
        parser.error("all bounds must be positive, with saturation bound at least three")

    started = perf_counter()
    base_document = json.loads(args.base.read_text())
    model_document = json.loads(args.model.read_text())
    section_document = json.loads(args.sections.read_text())
    bisection_document = json.loads(args.bisections.read_text())
    if base_document.get("status") != "PASS_EXACT_RANK19_SURFACE_OVER_RANK_AT_LEAST_9_BASE":
        raise ArithmeticError("the promoted paired-base certificate is missing")
    if bisection_document.get("schema") != "elkies-k3.bisection-extension-input.v1":
        raise ArithmeticError("the complete equation-level bisection input changed")
    records = bisection_document["bisections"]
    if len(records) != 39120:
        raise ArithmeticError("the complete bisection count changed")
    by_mask = {int(record["lattice_orbit_mask"]): record for record in records}
    if any(mask not in by_mask for mask in TARGET_MASKS):
        raise ArithmeticError("a promoted paired-base cover is missing")

    third_curve = EllipticCurve(QQ, list(THIRD_MODEL))
    paired_curve = EllipticCurve(QQ, list(PAIRED_MODEL))
    old_points = [
        third_curve(QQ(point[0]), QQ(point[1]))
        for point in base_document["paired_base"]["points"]
    ]
    rank_result = pari(third_curve).ellrank(0, pari(old_points))
    rank_interval = [int(rank_result[0]), int(rank_result[1])]
    saturated_raw = pari(third_curve).ellsaturation(pari(old_points), args.saturation_bound)
    saturated = [third_curve(QQ(point[0]), QQ(point[1])) for point in saturated_raw]
    if len(saturated) != 9:
        raise ArithmeticError("PARI saturation changed the displayed rank")

    old_gram = height_matrix(third_curve, old_points)
    saturated_gram = height_matrix(third_curve, saturated)
    combined_gram = height_matrix(third_curve, saturated + old_points)
    pairings = combined_gram[:9, 9:]
    real_relation = saturated_gram.solve_right(pairings)
    old_in_saturated = matrix(
        ZZ,
        [[ZZ(real_relation[row, column].round()) for column in range(9)] for row in range(9)],
    )
    for column, point in enumerate(old_points):
        if linear_combination(third_curve, saturated, old_in_saturated.column(column)) != point:
            raise ArithmeticError("the saturation relation failed exact group-law replay")
    relation_index = abs(int(old_in_saturated.det()))

    # PARI qfminim becomes unstable on this nine-dimensional form once the
    # integer scale reaches 1e10. Eight retained decimals, together with the
    # explicit enumeration guard above, are ample for the integral cutoffs.
    scale = ZZ(10) ** 8
    integral_gram = matrix(
        ZZ,
        [[ZZ((scale * value).round()) for value in row] for row in saturated_gram.rows()],
    )
    lll_transform = matrix(ZZ, pari(integral_gram).qflllgram())
    if abs(lll_transform.det()) != 1:
        raise ArithmeticError("the height LLL transform is not unimodular")
    lll_basis = [
        linear_combination(third_curve, saturated, lll_transform.row(row))
        for row in range(9)
    ]
    lll_gram = lll_transform * saturated_gram * lll_transform.transpose()
    lll_integral_gram = lll_transform * integral_gram * lll_transform.transpose()

    to_paired, to_cover = paired_cover_maps(third_curve, paired_curve)
    enumerated, pareto, enumeration_summary = enumerate_base_points(
        third_curve,
        lll_basis,
        lll_gram,
        lll_integral_gram,
        args.height_bound,
        scale,
        to_paired,
        to_cover,
    )
    finalist_map = {record["t"]: record for record in enumerated[: args.t_finalists]}
    for record in pareto:
        finalist_map.setdefault(record["t"], record)
    finalists = list(finalist_map.values())
    finalists.sort(
        key=lambda item: (
            item["t_projective_bits"],
            item["t_total_bits"],
            float(item["canonical_height"]),
        )
    )

    finite_helper = load_module("elkies_pair_finite_helper", FINITE_HELPER)
    R = PolynomialRing(QQ, "t")
    A_coefficients = tuple(QQ(value) for value in model_document["A_coefficients_low_to_high"])
    B_coefficients = tuple(QQ(value) for value in model_document["B_coefficients_low_to_high"])
    specialization_results = []
    for candidate in finalists:
        t_value = QQ(candidate["t"])
        coefficient_a = evaluate(A_coefficients, t_value)
        coefficient_b = evaluate(B_coefficients, t_value)
        curve = EllipticCurve(QQ, [coefficient_a, coefficient_b])
        generic_points = reconstruct_sections(section_document, t_value, curve)
        # The paired-base equations use primitive squareclass representatives,
        # whereas each lifted-section record retains its original rational
        # square multiple. Recompute the latter root exactly here.
        selected_roots = {}
        for mask in TARGET_MASKS:
            q_value = evaluate(by_mask[mask]["residual_chord"]["q_coefficients"], t_value)
            if not q_value.is_square():
                raise ArithmeticError("a defining paired cover stopped splitting")
            selected_roots[mask] = q_value.sqrt()
        selected_points = [
            lifted_point(by_mask[mask], t_value, selected_roots[mask], curve)
            for mask in TARGET_MASKS
        ]
        known_points = generic_points + selected_points
        baseline = best_finite_certificate(
            finite_helper,
            tuple(to_fraction(value) for value in (0, 0, 0, coefficient_a, coefficient_b)),
            known_points,
            args.certificate_prime_bound,
        )

        split_records = []
        extra_points = []
        seen_points = {(point[0], point[1]) for point in known_points}
        for bisection in records:
            mask = int(bisection["lattice_orbit_mask"])
            if mask in TARGET_MASKS:
                continue
            q_value = evaluate(bisection["residual_chord"]["q_coefficients"], t_value)
            if not q_value.is_square():
                continue
            square_root = q_value.sqrt()
            point = lifted_point(bisection, t_value, square_root, curve)
            if (point[0], point[1]) not in seen_points:
                extra_points.append(point)
                seen_points.add((point[0], point[1]))
            split_records.append(
                {
                    "lattice_orbit_mask": mask,
                    "orbit_hex": f"0x{mask:05x}",
                    "q_value": rational_text(q_value),
                    "square_root": rational_text(square_root),
                    "point": point_text(point),
                }
            )
        all_points = known_points + extra_points
        total = best_finite_certificate(
            finite_helper,
            tuple(to_fraction(value) for value in (0, 0, 0, coefficient_a, coefficient_b)),
            all_points,
            args.certificate_prime_bound,
        )
        specialization_results.append(
            {
                **candidate,
                "specialized_model_a1_a2_a3_a4_a6": [
                    rational_text(value) for value in (0, 0, 0, coefficient_a, coefficient_b)
                ],
                "known_section_count": 19,
                "known_points": [point_text(point) for point in known_points],
                "known_rank_certificate": baseline,
                "other_quadratics_tested": len(records) - 2,
                "other_split_bisection_count": len(split_records),
                "materialized_distinct_extra_point_count": len(extra_points),
                "split_bisections": split_records,
                "total_displayed_point_count": len(all_points),
                "total_rank_certificate": total,
            }
        )

    nagao = load_module("elkies_pair_nagao_helper", NAGAO_HELPER)
    prime_blocks = nagao.parse_prime_blocks(args.nagao_primes)
    if len(prime_blocks) < 3:
        raise ArithmeticError("Nagao ranking requires at least three disjoint prime blocks")
    family = nagao.load_family_model(args.model)
    table_blocks, rejected_primes = nagao.build_residue_tables(family, prime_blocks)
    for result in specialization_results:
        numerator = int(QQ(result["t"]).numerator())
        denominator = int(QQ(result["t"]).denominator())
        score = nagao.Candidate(numerator, denominator, max(abs(numerator), denominator))
        for block in table_blocks:
            score = nagao.score_block(score, block, {})
        result["post_split_nagao"] = nagao.candidate_record(score)
    specialization_results.sort(
        key=lambda item: (
            -int(item["other_split_bisection_count"]),
            -int(item["total_rank_certificate"]["certified_rank_lower_bound"]),
            -min(item["post_split_nagao"]["block_score_units_1e12"]),
            -int(item["post_split_nagao"]["total_score_units_1e12"]),
            int(item["t_projective_bits"]),
            item["t"],
        )
    )

    result = {
        "schema": "elkies-k3.elkies-2026-rank9-paired-base-search.v1",
        "status": "PASS_BOUNDED_EXACT_SPLIT_BISECTION_SEARCH",
        "inputs": {
            display_path(Path(__file__).resolve()): digest(Path(__file__).resolve()),
            display_path(args.base): digest(args.base),
            display_path(args.model): digest(args.model),
            display_path(args.sections): digest(args.sections),
            display_path(args.bisections): digest(args.bisections),
            display_path(FINITE_HELPER): digest(FINITE_HELPER),
            display_path(NAGAO_HELPER): digest(NAGAO_HELPER),
        },
        "base_lattice": {
            "third_quotient_model_a1_a2_a3_a4_a6": list(THIRD_MODEL),
            "pari_ellrank_interval": rank_interval,
            "pari_ellrank_sha_information": int(rank_result[2]),
            "saturation_prime_bound_strict": args.saturation_bound,
            "saturation_basis": [point_text(point) for point in saturated],
            "old_basis_in_saturation_basis_columns": [
                [int(value) for value in row] for row in old_in_saturated.rows()
            ],
            "exact_old_to_saturation_index": relation_index,
            "old_maximum_coordinate_bits": max(map(coordinate_bits, old_points)),
            "saturation_maximum_coordinate_bits": max(map(coordinate_bits, saturated)),
            "canonical_height_precision_decimal_digits": 100,
            "canonical_height_gram": matrix_text(saturated_gram),
            "canonical_height_determinant": f"{saturated_gram.det():.70}",
            "lll_transform_rows": [[int(value) for value in row] for row in lll_transform.rows()],
            "lll_basis": [point_text(point) for point in lll_basis],
            "lll_reduced_canonical_height_gram": matrix_text(lll_gram),
            "lll_diagonal": [f"{value:.60}" for value in lll_gram.diagonal()],
            "interpretation": (
                "PARI reports rank interval [9,9]. ellsaturation through primes below the "
                f"declared bound returns a smaller basis, but the exact transition has index {relation_index}; "
                "this is a basis improvement, not an added point or a new rank claim."
            ),
        },
        "enumeration": {
            "canonical_height_bound": args.height_bound,
            "primary_candidate_order": "t_projective_bits, then t_total_bits, then canonical_height",
            "canonical_height_is_not_the_specialization_order": True,
            "requested_t_finalists": args.t_finalists,
            "pareto_frontier_size": len(pareto),
            "pareto_frontier": pareto,
            **enumeration_summary,
            "specialized_distinct_t_count": len(finalists),
        },
        "specializations": specialization_results,
        "nagao_policy": {
            "applied_after_explicit_split_count": True,
            "prime_blocks": [list(block) for block in prime_blocks],
            "rejected_primes": list(rejected_primes),
            "sort_order": (
                "other split count; certified rank lower bound; weakest Nagao block; "
                "total Nagao score; t complexity"
            ),
        },
        "runtime_seconds": perf_counter() - started,
        "reproducing_command": shlex.join(sys.argv),
        "proof_boundary": (
            "Every listed point, square test, and finite-quotient lower bound is exact. The base-height "
            "bound and t-finalist cap are bounded search choices. No missing split, rank upper bound, "
            "or exact specialized rank is inferred. Nagao scores are heuristic and are evaluated only "
            "after the explicit split count."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    best = specialization_results[0]
    print(
        "ELKIES2026R9PAIRSEARCH|"
        f"height_bound={args.height_bound}|enumerated_t={enumeration_summary['distinct_finite_t_values']}|"
        f"specialized_t={len(finalists)}|best_other_splits={best['other_split_bisection_count']}|"
        f"best_rank_lower_bound={best['total_rank_certificate']['certified_rank_lower_bound']}|"
        f"status={result['status']}|output={display_path(args.output)}"
    )


if __name__ == "__main__":
    main()
