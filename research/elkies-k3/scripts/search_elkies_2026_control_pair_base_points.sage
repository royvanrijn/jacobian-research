#!/usr/bin/env sage-python
"""Materialize and specialize points from the t=3/8 control-selected pair bases.

status: ACTIVE_SEARCH
claim: bounded exact pair-base point enumeration and split-bisection specialization
inputs: control pair catalogue/rank ledger, R17 sections, complete bisection batch
outputs: artifacts/generated-results/elkies-2026-control-pair-base-point-search.json
supersedes: none; turns certified positive pair-base rank into explicit parameters

Each certified Jacobian basis is canonical-height LLL-reduced.  The declared
coefficient box is then enumerated on the pair-base elliptic curve, transported
through the exact pointed-quartic inverse, and ordered separately by the height
of t.  All 39,120 bisections are sieved and exactly tested at every distinct t.
With ``--certify all``, the generic sections and both signs of every split
bisection are materialized and the resulting rank is certified by finite
quotients.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import importlib.util
from itertools import product
import json
from pathlib import Path
import shlex
import sys
from time import perf_counter

from sage.all import EllipticCurve, PolynomialRing, QQ, RealField, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
CATALOGUE = ROOT / "artifacts/generated-results/elkies-2026-control-pair-base-catalogue.json"
RANK_LEDGER = ROOT / "artifacts/generated-results/elkies-2026-control-pair-base-rank-ledger.json"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
BISECTIONS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
SEARCH_HELPER = ROOT / "elkies-k3/scripts/search_elkies_2026_rank9_paired_base.sage"
SIEVE_HELPER = ROOT / "elkies-k3/scripts/sieve_elkies_2026_rank9_paired_base.sage"
FINITE_HELPER = ROOT / "elliptic-curves/cas/elliptic_candidate_record.py"
NAGAO_HELPER = ROOT / "elkies-k3/scripts/search_h92_q12o5867_rootless_nagao.py"
OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-control-pair-base-point-search.json"
T0 = QQ(3) / 8


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_loader(name, SourceFileLoader(name, str(path)))
    if specification is None or specification.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def projective_bits(value) -> int:
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


def point_digest(points) -> str:
    payload = "\n".join(
        "0" if point.is_zero() else f"{rational_text(point[0])}|{rational_text(point[1])}"
        for point in points
    ) + "\n"
    return sha256(payload.encode()).hexdigest()


def linear_combination(curve, points, coefficients):
    return sum((ZZ(coefficient) * point for coefficient, point in zip(coefficients, points)), curve(0))


def direct_independence_certificate(search, helper, model, points, prime_bound):
    fraction_points = tuple(search.fraction_point(point) for point in points)
    coefficients = helper.model_from_record(model)
    original_curve_check = helper.is_on_weierstrass_curve
    if any(not original_curve_check(coefficients, point) for point in fraction_points):
        raise ArithmeticError("a selected independence point missed the specialized curve")
    # The shared helper defensively repeats the same rational on-curve test in
    # every reduction-prime signature and again during replay.  On these very
    # large specializations that dominates the finite-field calculation.  We
    # perform it once above, then suppress only the redundant copies while the
    # exact modular signatures are built and replayed.
    helper.is_on_weierstrass_curve = lambda _model, _point: True
    try:
        attempts = []
        for relation_prime in (2, 3):
            certificate = helper.build_finite_quotient_certificate(
                model,
                fraction_points,
                relation_prime=relation_prime,
                prime_bound=prime_bound,
            )
            attempts.append(
                {
                    "relation_prime": relation_prime,
                    "combined_rank_over_relation_field": int(
                        certificate["combined_rank_over_relation_field"]
                    ),
                }
            )
            if certificate["certified_independent"]:
                helper.verify_finite_quotient_certificate(model, fraction_points, certificate)
                return {
                    "certified_rank_lower_bound": len(points),
                    "independent_point_indices_zero_based": list(range(len(points))),
                    "attempts": attempts,
                    "successful_certificate": certificate,
                }
        return {
            "certified_rank_lower_bound": 0,
            "independent_point_indices_zero_based": [],
            "attempts": attempts,
            "successful_certificate": None,
        }
    finally:
        helper.is_on_weierstrass_curve = original_curve_check


def pair_transport(row, generators):
    t_ring = PolynomialRing(QQ, "t")
    s_ring = PolynomialRing(QQ, "s")
    s = s_ring.gen()
    left_q = t_ring([QQ(value) for value in row["q_coefficients_low_to_high"][0]])
    right_q = t_ring([QQ(value) for value in row["q_coefficients_low_to_high"][1]])
    base = next(point for point in row["common_control_points"] if point["t"] == "3/8")
    u0 = QQ(base["u"])
    v0 = QQ(base["v"])
    leading = QQ(left_q[2])
    derivative = QQ(left_q.derivative()(T0))
    denominator = 1 - leading * s**2
    t_of_s = T0 + (derivative * s**2 - 2 * u0 * s) / denominator
    quartic = s_ring(denominator**2 * right_q(t_of_s))
    e, d, c, b, a = (QQ(quartic[index]) for index in range(5))
    paired_curve = EllipticCurve(
        QQ,
        [
            d / v0,
            c - d**2 / (4 * v0**2),
            2 * v0 * b,
            -4 * v0**2 * a,
            a * (d**2 - 4 * v0**2 * c),
        ],
    )
    invariant_i = QQ(row["binary_quartic_I"])
    invariant_j = QQ(row["binary_quartic_J"])
    third_curve = EllipticCurve(QQ, [0, 0, 0, -27 * invariant_i, -27 * invariant_j])
    minimal_curve = EllipticCurve(QQ, row["minimal_jacobian_a1_a2_a3_a4_a6"])
    minimal_to_third = minimal_curve.isomorphism_to(third_curve)
    isogenies = [
        morphism
        for morphism in third_curve.isogenies_prime_degree(2)
        if morphism.codomain().is_isomorphic(paired_curve)
    ]
    if len(isogenies) != 1:
        raise ArithmeticError(f"pair {row['pair_key']} does not have a unique expected 2-isogeny")
    isogeny = isogenies[0]
    codomain_to_paired = (
        None
        if isogeny.codomain() == paired_curve
        else isogeny.codomain().isomorphism_to(paired_curve)
    )

    def to_pair(point):
        image = isogeny(minimal_to_third(point))
        return image if codomain_to_paired is None else codomain_to_paired(image)

    def to_cover(point):
        image = to_pair(point)
        if image.is_zero() or image[1] == 0:
            return None
        x_value, y_value = image[:2]
        s_value = (4 * v0**2 * (x_value + c) - d**2) / (2 * v0 * y_value)
        if s_value == 0 or denominator(s_value) == 0:
            return None
        quartic_y = (x_value * s_value**2 - d * s_value) / (2 * v0) - v0
        if quartic_y**2 != quartic(s_value):
            raise ArithmeticError("pointed-quartic inverse failed")
        t_value = t_of_s(s_value)
        u_value = u0 + (t_value - T0) / s_value
        v_value = quartic_y / denominator(s_value)
        if u_value**2 != left_q(t_value) or v_value**2 != right_q(t_value):
            raise ArithmeticError("pair-cover inverse failed")
        return t_value, u_value, v_value

    return minimal_curve, paired_curve, to_pair, to_cover


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", type=Path, default=CATALOGUE)
    parser.add_argument("--rank-ledger", type=Path, default=RANK_LEDGER)
    parser.add_argument("--model", type=Path, default=MODEL)
    parser.add_argument("--sections", type=Path, default=SECTIONS)
    parser.add_argument("--bisections", type=Path, default=BISECTIONS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--coefficient-radius", type=int, default=1)
    parser.add_argument("--certify", choices=("none", "all"), default="none")
    parser.add_argument("--certificate-prime-bound", type=int, default=300)
    parser.add_argument("--checkpoint-every", type=int, default=500)
    parser.add_argument(
        "--nagao-primes",
        default="19,41,43,61,71,73,79,83;89,107,113,127,131,137,139,151;157,163,167,173,179,181,191,193,197",
    )
    args = parser.parse_args()
    if args.coefficient_radius < 1:
        parser.error("coefficient-radius must be positive")
    started = perf_counter()

    search = load_module("elkies_control_pair_search_helper", SEARCH_HELPER)
    sieve = load_module("elkies_control_pair_sieve_helper", SIEVE_HELPER)
    finite_helper = load_module("elkies_control_pair_finite_helper", FINITE_HELPER)
    nagao = load_module("elkies_control_pair_nagao_helper", NAGAO_HELPER)
    catalogue = json.loads(args.catalogue.read_text())
    rank_ledger = json.loads(args.rank_ledger.read_text())
    model_document = json.loads(args.model.read_text())
    section_document = json.loads(args.sections.read_text())
    bisection_document = json.loads(args.bisections.read_text())
    if catalogue.get("schema") != "elkies-k3.elkies-2026-control-pair-base-catalogue.v1":
        raise ValueError("unexpected pair catalogue")
    if rank_ledger.get("schema") != "elkies-k3.elkies-2026-control-pair-base-rank-ledger.v1":
        raise ValueError("unexpected pair rank ledger")
    rows = {row["pair_key"]: row for row in catalogue["pairs"]}

    pari.set_real_precision(200)
    real_field = RealField(300)
    parameter_sources = {}
    pair_summaries = []
    exceptional = 0
    for result in sorted(rank_ledger["results"].values(), key=lambda item: item["arithmetic_complexity_rank"]):
        rank = int(result["certified_rank_lower_bound"])
        if rank == 0:
            continue
        row = rows[result["pair_key"]]
        generators_raw = result["generators"]
        minimal_curve = EllipticCurve(QQ, row["minimal_jacobian_a1_a2_a3_a4_a6"])
        generators = [minimal_curve(QQ(point[0]), QQ(point[1])) for point in generators_raw]
        if len(generators) != rank:
            raise ArithmeticError("certified rank and generator count differ")
        height_pari = pari(minimal_curve).ellheightmatrix(pari(generators))
        height_gram = matrix(
            real_field,
            rank,
            rank,
            [[real_field(str(height_pari[row_index, column_index])) for column_index in range(rank)] for row_index in range(rank)],
        )
        scale = ZZ(10) ** 8
        integral_gram = matrix(ZZ, [[ZZ((scale * value).round()) for value in line] for line in height_gram.rows()])
        lll_transform = matrix(ZZ, pari(integral_gram).qflllgram())
        if abs(lll_transform.det()) != 1:
            raise ArithmeticError("pair-base LLL transform is not unimodular")
        lll_points = [linear_combination(minimal_curve, generators, lll_transform.row(index)) for index in range(rank)]
        lll_gram = lll_transform * height_gram * lll_transform.transpose()
        _minimal, _paired, _to_pair, to_cover = pair_transport(row, generators)
        enumerated = 0
        mapped = 0
        for coefficients_tuple in product(range(-args.coefficient_radius, args.coefficient_radius + 1), repeat=rank):
            if not any(coefficients_tuple):
                continue
            enumerated += 1
            coefficients = vector(ZZ, coefficients_tuple)
            point = linear_combination(minimal_curve, lll_points, coefficients)
            cover = to_cover(point)
            if cover is None:
                exceptional += 1
                continue
            mapped += 1
            t_value, u_value, v_value = cover
            key = rational_text(t_value)
            source = {
                "pair_key": row["pair_key"],
                "orbit_masks": row["orbit_masks"],
                "lll_coefficients": list(map(int, coefficients_tuple)),
                "canonical_height": f"{coefficients * lll_gram * coefficients:.50}",
                "u": rational_text(u_value),
                "v": rational_text(v_value),
            }
            parameter_sources.setdefault(key, []).append(source)
        pair_summaries.append(
            {
                "pair_key": row["pair_key"],
                "certified_base_rank_lower_bound": rank,
                "lll_transform_rows": [[int(value) for value in line] for line in lll_transform.rows()],
                "lll_height_diagonal": [f"{value:.40}" for value in lll_gram.diagonal()],
                "enumerated_nonzero_vectors": enumerated,
                "mapped_cover_points": mapped,
            }
        )

    all_bisections = bisection_document["bisections"]
    quadratics = [
        sieve.normalized_homogeneous_quadratic(record["residual_chord"]["q_coefficients"])
        for record in all_bisections
    ]
    tables, all_bits = sieve.build_residue_bitsets(quadratics, sieve.DEFAULT_PRIMES)
    a_coefficients = tuple(QQ(value) for value in model_document["A_coefficients_low_to_high"])
    b_coefficients = tuple(QQ(value) for value in model_document["B_coefficients_low_to_high"])
    fibres = []
    exact_tests = 0
    certified_counts = {}
    for ordinal, (t_text, sources) in enumerate(
        sorted(parameter_sources.items(), key=lambda item: (projective_bits(QQ(item[0])), item[0])),
        start=1,
    ):
        t_value = QQ(t_text)
        numerator = ZZ(t_value.numerator())
        denominator = ZZ(t_value.denominator())
        survivor_indices = tuple(sieve.sieve_indices(numerator, denominator, tables, all_bits))
        split = []
        collisions = []
        for index in survivor_indices:
            exact_tests += 1
            exact = sieve.exact_square_value(*quadratics[index], numerator, denominator)
            if exact is None:
                continue
            q_value, root = exact
            record = all_bisections[index]
            entry = {
                "lattice_orbit_mask": int(record["lattice_orbit_mask"]),
                "orbit_hex": f"0x{int(record['lattice_orbit_mask']):05x}",
                "q_value": rational_text(q_value),
                "canonical_positive_square_root": rational_text(root),
            }
            (collisions if q_value == 0 else split).append((entry, record, root))
        split_masks = {entry["lattice_orbit_mask"] for entry, _record, _root in split}
        expected_masks = {mask for source in sources for mask in source["orbit_masks"]}
        if not expected_masks <= split_masks:
            raise ArithmeticError("a materialized pair point failed its defining square tests")
        fibre = {
            "t": t_text,
            "t_projective_bits": projective_bits(t_value),
            "source_count": len(sources),
            "sources": sources,
            "modular_survivor_count": len(survivor_indices),
            "split_bisection_count": len(split),
            "split_masks": sorted(split_masks),
            "branch_collisions": [entry for entry, _record, _root in collisions],
        }
        if args.certify == "all":
            coefficient_a = search.evaluate(a_coefficients, t_value)
            coefficient_b = search.evaluate(b_coefficients, t_value)
            curve = EllipticCurve(QQ, [coefficient_a, coefficient_b])
            generic_points = search.reconstruct_sections(section_document, t_value, curve)
            points = list(generic_points)
            certificate_points = list(generic_points)
            split_records = []
            for entry, record, root in split:
                positive = search.lifted_point(record, t_value, root, curve)
                negative = search.lifted_point(record, t_value, -root, curve)
                points.extend((positive, negative))
                certificate_points.append(positive)
                split_records.append(
                    {
                        **entry,
                        "positive_point": search.point_text(positive),
                        "negative_point": search.point_text(negative),
                    }
                )
            unique_points = []
            seen = set()
            for point in points:
                key = ("zero",) if point.is_zero() else (point[0], point[1])
                if key not in seen:
                    seen.add(key)
                    unique_points.append(point)
            model = tuple(search.to_fraction(value) for value in (0, 0, 0, coefficient_a, coefficient_b))
            certificate = direct_independence_certificate(
                search,
                finite_helper,
                model,
                certificate_points,
                args.certificate_prime_bound,
            )
            if certificate["successful_certificate"] is None:
                certificate = sieve.compact_certificate(
                    search, finite_helper, model, unique_points, args.certificate_prime_bound
                )
            certified_rank = int(certificate["certified_rank_lower_bound"])
            certified_counts[str(certified_rank)] = certified_counts.get(str(certified_rank), 0) + 1
            fibre.update(
                {
                    "specialized_model_a1_a2_a3_a4_a6": [
                        rational_text(value) for value in (0, 0, 0, coefficient_a, coefficient_b)
                    ],
                    "generic_point_count": 17,
                    "split_bisections": split_records,
                    "materialized_unique_point_count": len(unique_points),
                    "materialized_point_sequence_sha256": point_digest(unique_points),
                    "total_rank_certificate": certificate,
                }
            )
        fibres.append(fibre)
        if ordinal % args.checkpoint_every == 0:
            print(
                "ELKIES2026CONTROLPAIRPOINTS|"
                f"progress={ordinal}/{len(parameter_sources)}|max_splits={max(item['split_bisection_count'] for item in fibres)}|"
                f"exact_tests={exact_tests}",
                flush=True,
            )

    prime_blocks = nagao.parse_prime_blocks(args.nagao_primes)
    if len(prime_blocks) < 3:
        raise ArithmeticError("Nagao ranking requires at least three disjoint prime blocks")
    family = nagao.load_family_model(args.model)
    table_blocks, rejected_primes = nagao.build_residue_tables(family, prime_blocks)
    for fibre in fibres:
        t_value = QQ(fibre["t"])
        score = nagao.Candidate(
            int(t_value.numerator()),
            int(t_value.denominator()),
            max(abs(int(t_value.numerator())), int(t_value.denominator())),
        )
        for block in table_blocks:
            score = nagao.score_block(score, block, {})
        fibre["post_split_nagao"] = nagao.candidate_record(score)
    fibres.sort(
        key=lambda item: (
            -item["split_bisection_count"],
            -(
                item.get("total_rank_certificate", {}).get("certified_rank_lower_bound", 0)
                if item.get("total_rank_certificate")
                else 0
            ),
            -min(item["post_split_nagao"]["block_score_units_1e12"]),
            -int(item["post_split_nagao"]["total_score_units_1e12"]),
            item["t_projective_bits"],
            item["t"],
        )
    )
    result = {
        "schema": "elkies-k3.elkies-2026-control-pair-base-point-search.v1",
        "status": "PASS_BOUNDED_EXACT_CONTROL_PAIR_BASE_POINT_SEARCH",
        "inputs": {
            display_path(Path(__file__).resolve()): digest(Path(__file__).resolve()),
            display_path(args.catalogue): digest(args.catalogue),
            display_path(args.rank_ledger): digest(args.rank_ledger),
            display_path(args.model): digest(args.model),
            display_path(args.sections): digest(args.sections),
            display_path(args.bisections): digest(args.bisections),
            display_path(SEARCH_HELPER): digest(SEARCH_HELPER),
            display_path(SIEVE_HELPER): digest(SIEVE_HELPER),
            display_path(FINITE_HELPER): digest(FINITE_HELPER),
            display_path(NAGAO_HELPER): digest(NAGAO_HELPER),
        },
        "bounds": {
            "lll_coefficient_radius": args.coefficient_radius,
            "coefficient_box": f"[-{args.coefficient_radius},{args.coefficient_radius}]^r minus zero",
            "certification_policy": args.certify,
            "certificate_prime_bound": args.certificate_prime_bound,
            "sieve_primes": list(sieve.DEFAULT_PRIMES),
        },
        "summary": {
            "positive_rank_pair_count": len(pair_summaries),
            "enumerated_nonzero_vector_count": sum(row["enumerated_nonzero_vectors"] for row in pair_summaries),
            "mapped_pair_cover_point_count": sum(row["mapped_cover_points"] for row in pair_summaries),
            "exceptional_map_point_count": exceptional,
            "distinct_finite_t_count": len(parameter_sources),
            "exact_square_test_count": exact_tests,
            "maximum_split_bisection_count": max((row["split_bisection_count"] for row in fibres), default=0),
            "fibres_with_more_than_two_splits": sum(row["split_bisection_count"] > 2 for row in fibres),
            "certified_rank_lower_bound_counts": dict(sorted(certified_counts.items(), key=lambda item: int(item[0]))),
        },
        "pair_lattices": pair_summaries,
        "fibres": fibres,
        "nagao_policy": {
            "applied_after_explicit_split_count_and_rank_certification": True,
            "prime_blocks": [list(block) for block in prime_blocks],
            "rejected_primes": list(rejected_primes),
            "sort_order": (
                "explicit split count; certified rank lower bound; weakest Nagao block; "
                "total Nagao score; t projective height"
            ),
        },
        "runtime_seconds": perf_counter() - started,
        "reproducing_command": shlex.join(sys.argv),
        "proof_boundary": (
            "All point transports and square tests are exact inside the declared LLL coefficient boxes. "
            "The boxes are bounded searches, not exhaustive Mordell-Weil lattices. Rank lower bounds are "
            "asserted only when the corresponding finite-quotient certificate is stored."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        "ELKIES2026CONTROLPAIRPOINTS|"
        f"pairs={len(pair_summaries)}|mapped={result['summary']['mapped_pair_cover_point_count']}|"
        f"distinct_t={len(parameter_sources)}|max_splits={result['summary']['maximum_split_bisection_count']}|"
        f"status={result['status']}|output={display_path(args.output)}"
    )


if __name__ == "__main__":
    main()
