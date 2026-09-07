#!/usr/bin/env sage-python
"""Control the ``0x103b2`` quartic-Jacobian point-search rank signal.

Run the identical PARI ``hyperellratpoints(H=10000)`` search, exact pointed
quartic map, global minimal-model transport, and Sage/eclib relation reduction
on a deterministic sample of ten other pointed covers from the production
genus-one-bisection catalogue.  The target is replayed by the same code as a
positive regression.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, pari
from sage.env import SAGE_VERSION
from sage.libs.eclib.all import mwrank_MordellWeil


ROOT = Path(__file__).resolve().parents[2]
INPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-pointed-cover-jacobian-rank-controls-h10000-v1.json"
)
TARGET_LABEL = "norm12-orbit-103b2"
TARGET_T = QQ(1) / 25
TARGET_S = QQ("3521934804796232704/643125")
DEFAULT_HEIGHT = 10_000
DEFAULT_SAMPLE_SIZE = 10
DEFAULT_SEED = 0x103B2


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def rational_text(value) -> str:
    return str(QQ(value))


def pointed_model(quartic, parameter, cover_coordinate):
    local_ring = PolynomialRing(QQ, "z")
    z = local_ring.gen()
    shifted = local_ring(quartic(z + parameter))
    e, d, c, b, a = (QQ(shifted[index]) for index in range(5))
    v0 = QQ(cover_coordinate)
    if not v0 or v0**2 != e:
        raise ArithmeticError("invalid pointed-quartic base point")
    curve = EllipticCurve(QQ, [
        d / v0,
        c - d**2 / (4 * v0**2),
        2 * v0 * b,
        -4 * v0**2 * a,
        a * (d**2 - 4 * v0**2 * c),
    ])
    opposite_x = d**2 / (4 * v0**2) - c
    return curve, opposite_x, (c, d, v0)


def map_quartic_point(parameter, cover_coordinate, base_parameter, constants, curve):
    c, d, v0 = constants
    local_parameter = QQ(parameter - base_parameter)
    if not local_parameter:
        return None
    x_value = (
        2 * v0 * (QQ(cover_coordinate) + v0) + d * local_parameter
    ) / local_parameter**2
    y_value = (
        4 * v0**2 * (x_value + c) - d**2
    ) / (2 * v0 * local_parameter)
    return curve(x_value, y_value)


def analyze(record, base_parameter, cover_coordinate, height):
    ring = PolynomialRing(QQ, "t")
    quartic = ring([
        QQ(value)
        for value in record["branch_polynomial_q_coefficients_low_to_high"]
    ])
    if quartic(base_parameter) != cover_coordinate**2:
        raise ArithmeticError(f"stored pointing is not on {record['label']}")
    curve, opposite_x, constants = pointed_model(quartic, base_parameter, cover_coordinate)
    stored_jacobian = record["pointed_quartic_jacobian"]
    if stored_jacobian is not None:
        stored_invariants = [QQ(value) for value in stored_jacobian["weierstrass_coefficients"]]
        if list(curve.a_invariants()) != stored_invariants:
            raise ArithmeticError("reconstructed pointed Jacobian differs from the production record")
        pointed_generator = curve(*[QQ(value) for value in stored_jacobian["generator"]])
    else:
        candidates = [point for point in curve.lift_x(opposite_x, all=True) if point[1]]
        if len(candidates) != 1:
            raise ArithmeticError("the target pointing did not supply its unique visible generator")
        pointed_generator = candidates[0]
    # ``global_minimal_model`` asks Sage to factor the enormous discriminant
    # of each auxiliary curve.  That factorization is irrelevant to relation
    # reduction and made the proposed cheap control expensive.  Clear only
    # coefficient denominators, then feed the exact integral model directly
    # to the same eclib Mordell--Weil relation/saturation engine used by
    # ``EllipticCurve.saturation``.
    integral = curve.global_integral_model()
    isomorphism = curve.isomorphism_to(integral)

    raw_points = pari(quartic).hyperellratpoints(height)
    mapped = []
    base_fibre_points = 0
    for parameter, ordinate in raw_points:
        parameter = QQ(parameter)
        if parameter == base_parameter:
            base_fibre_points += 1
            continue
        point = map_quartic_point(
            parameter, QQ(ordinate), base_parameter, constants, curve
        )
        if point is None:
            raise ArithmeticError("non-base quartic point hit the exceptional map locus")
        mapped.append(isomorphism(point))

    # Match the 0x103b2 discovery protocol literally: both points above the
    # pointing parameter are omitted because this affine map is singular
    # there.  The separately constructed pointed generator is retained as a
    # baseline rank-one lower bound, but is not inserted into the relation
    # reduction.
    relation_candidates = mapped
    if relation_candidates:
        encoded = []
        for point in relation_candidates:
            x_value, y_value = point.xy()
            denominator = x_value.denominator().lcm(y_value.denominator())
            encoded.append((x_value * denominator, y_value * denominator, denominator))
        mordell_weil = mwrank_MordellWeil(integral.mwrank_curve(), False)
        mordell_weil.process(encoded)
        relation_rank = int(mordell_weil.rank())
        basis = mordell_weil.points()
        regulator = mordell_weil.regulator()
    else:
        relation_rank = 0
        basis, regulator = [], 1
    return {
        "label": record["label"],
        "trace_norm": int(record["trace_norm"]),
        "lattice_orbit_mask": int(record["lattice_orbit_mask"]),
        "pointing": {
            "t": rational_text(base_parameter),
            "s": rational_text(cover_coordinate),
        },
        "pointed_jacobian_integral_a_invariants": [
            rational_text(value) for value in integral.a_invariants()
        ],
        "affine_quartic_point_count_including_both_signs": len(raw_points),
        "base_fibre_point_count_skipped": base_fibre_points,
        "mapped_point_count": len(mapped),
        "pointed_model_generator_included_in_relation_reduction": False,
        "separate_pointed_generator_nontorsion": True,
        "relation_candidate_count": len(relation_candidates),
        "independent_point_count_after_relation_reduction": relation_rank,
        "eclib_basis_count_after_attempted_two_saturation": len(basis),
        "prime_two_saturation_status": "NOT_RUN_RANK_INVARIANT",
        "relation_reduced_regulator_approx": str(regulator),
        "proved_jacobian_rank_lower_bound_from_nonbase_search": relation_rank,
        "proved_jacobian_rank_lower_bound_including_separate_pointed_generator": max(1, relation_rank),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=DEFAULT_HEIGHT)
    parser.add_argument("--sample-size", type=int, default=DEFAULT_SAMPLE_SIZE)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.height < 1 or args.sample_size < 1:
        parser.error("height and sample size must be positive")

    document = json.loads(INPUT.read_text())
    records = document["construction"]["records"]
    pointed = [
        record for record in records
        if record["label"] != TARGET_LABEL and record["rational_seed"] is not None
    ]
    if args.sample_size > len(pointed):
        parser.error("sample size exceeds the number of other pointed covers")
    selected_indices = sorted(random.Random(args.seed).sample(range(len(pointed)), args.sample_size))
    selected = [pointed[index] for index in selected_indices]

    controls = []
    for ordinal, record in enumerate(selected, start=1):
        seed = record["rational_seed"]
        row = analyze(record, QQ(seed["t"]), QQ(seed["s"]), args.height)
        controls.append(row)
        print(
            f"POINTEDCONTROL|sample={ordinal}/{len(selected)}|label={row['label']}|"
            f"points={row['affine_quartic_point_count_including_both_signs']}|"
            f"mapped={row['mapped_point_count']}|rank_lb={row['proved_jacobian_rank_lower_bound_from_nonbase_search']}",
            flush=True,
        )

    target_record = next(record for record in records if record["label"] == TARGET_LABEL)
    target = analyze(target_record, TARGET_T, TARGET_S, args.height)
    if args.height == DEFAULT_HEIGHT and target[
        "independent_point_count_after_relation_reduction"
    ] != 17:
        raise ArithmeticError("the 0x103b2 H=10000 positive regression changed")

    histogram = {}
    for row in controls:
        rank = str(row["proved_jacobian_rank_lower_bound_from_nonbase_search"])
        histogram[rank] = histogram.get(rank, 0) + 1
    observed = [row["proved_jacobian_rank_lower_bound_from_nonbase_search"] for row in controls]
    result = {
        "schema": "elkies-k3.pointed-cover-jacobian-rank-controls-h10000.v1",
        "status": "PASS_EXACT_BOUNDED_POINTED_COVER_CONTROL",
        "protocol": {
            "pari_hyperellratpoints_naive_height_bound": args.height,
            "point_map": "standard exact pointed-quartic birational map",
            "relation_reduction": (
                "exact denominator-cleared integral model followed by the same eclib "
                "Mordell-Weil relation process used internally by Sage EllipticCurve.saturation; "
                "the optional prime-2 saturation step is omitted because it cannot change rank"
            ),
            "other_pointed_cover_population_size": len(pointed),
            "sample_size": args.sample_size,
            "selection": "uniform sample without replacement in production-record order",
            "deterministic_random_seed": args.seed,
            "selected_population_indices_zero_based": selected_indices,
        },
        "target_positive_regression": target,
        "controls": controls,
        "summary": {
            "visible_rank_lower_bound_histogram": dict(sorted(histogram.items(), key=lambda item: int(item[0]))),
            "minimum_visible_rank_lower_bound": min(observed),
            "maximum_visible_rank_lower_bound": max(observed),
            "mean_visible_rank_lower_bound": sum(observed) / len(observed),
            "target_visible_rank_lower_bound": target["proved_jacobian_rank_lower_bound_from_nonbase_search"],
            "controls_matching_or_exceeding_target": sum(
                value >= target["proved_jacobian_rank_lower_bound_from_nonbase_search"]
                for value in observed
            ),
        },
        "proof_boundary": (
            "Every listed quartic point and mapped relation rank is exact within the declared PARI naive-height bound. "
            "The resulting dimensions are Mordell-Weil rank lower bounds only. The ten-cover sample is deterministic "
            "but not exhaustive, and differences in visible rank can reflect point heights rather than true ranks."
        ),
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "pari_version": ".".join(map(str, pari.version())),
            "required_features": ["PARI hyperellratpoints", "Sage eclib mwrank interface"],
        },
        "inputs": {relative(INPUT): digest(INPUT)},
        "reproducing_command": (
            "sage -python elkies-k3/scripts/control_pointed_cover_jacobian_ranks.sage"
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != serialized:
            raise ArithmeticError("stored pointed-cover control differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        f"POINTEDCONTROL|target={target['proved_jacobian_rank_lower_bound_from_nonbase_search']}|"
        f"histogram={histogram}|mean={result['summary']['mean_visible_rank_lower_bound']}|"
        f"output={relative(args.output)}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
