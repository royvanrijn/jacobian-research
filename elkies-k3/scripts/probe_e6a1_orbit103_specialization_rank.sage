#!/usr/bin/env sage-python
"""Search orbit-103 finalist fibres and certify gains beyond Q_plus,Q_minus.

This is the exact promotion gate for
``search_e6a1_orbit103_specializations.py``.  Each finalist already contains
an integral short model and the exact specializations of ``Q_plus,Q_minus``.
eclib/mwrank is seeded with those two points and runs a bounded point search.
Every returned affine point is checked in exact arithmetic.  A point is
counted as a new direction only when the repository's finite-quotient
certificate proves the enlarged ordered set independent.

Failure to find or certify a point is only a bounded negative experiment.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, lcm
from pathlib import Path
import shutil
import subprocess
import sys
from time import perf_counter

from sage.all import EllipticCurve, QQ, ZZ, pari, version as sage_version
from sage.libs.eclib.interface import mwrank_EllipticCurve, mwrank_MordellWeil


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from elliptic_candidate_record import (  # noqa: E402
    build_finite_quotient_certificate,
    is_on_weierstrass_curve,
    verify_finite_quotient_certificate,
)


DEFAULT_INPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-e6a1-orbit103-specialization-search-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-e6a1-orbit103-specialization-rank-probes-v1.json"
)
Q = Fraction


def parse_point(record):
    return Q(record[0]), Q(record[1])


def point_record(point):
    return [str(point[0]), str(point[1])]


def to_mwrank_triple(point):
    x_coordinate, y_coordinate = point
    denominator = lcm(x_coordinate.denominator, y_coordinate.denominator)
    projective_x = int(x_coordinate * denominator)
    projective_y = int(y_coordinate * denominator)
    common = gcd(gcd(abs(projective_x), abs(projective_y)), denominator)
    return [projective_x // common, projective_y // common, denominator // common]


def from_mwrank_triple(triple):
    if len(triple) != 3:
        raise ValueError("an mwrank point needs three projective coordinates")
    projective_x, projective_y, denominator = map(int, triple)
    if denominator == 0:
        return None
    return Q(projective_x, denominator), Q(projective_y, denominator)


def sign_key(point):
    x_coordinate, y_coordinate = point
    direct = (
        (x_coordinate.numerator, x_coordinate.denominator),
        (y_coordinate.numerator, y_coordinate.denominator),
    )
    negative = (
        direct[0],
        (-y_coordinate.numerator, y_coordinate.denominator),
    )
    return min(direct, negative)


def sage_point_tuple(point):
    return (
        Q(int(point[0].numerator()), int(point[0].denominator())),
        Q(int(point[1].numerator()), int(point[1].denominator())),
    )


def certificate(model, points, relation_prime, prime_bound):
    result = build_finite_quotient_certificate(
        model,
        points,
        relation_prime=relation_prime,
        prime_bound=prime_bound,
    )
    verify_finite_quotient_certificate(model, points, result)
    return result


def external_gp_ellrank(minimal_model, minimal_known, effort, timeout):
    """Run the system GP as a fallback for Sage/PARI ``small_norm`` failures."""

    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("system gp is unavailable for the ellrank fallback")
    model_literal = ",".join(str(value) for value in minimal_model)
    points_literal = ",".join(
        f"[{point[0]},{point[1]}]" for point in minimal_known
    )
    program = (
        f"E=ellinit([{model_literal}]);\n"
        f"print(ellrank(E,{effort},[{points_literal}]));\n"
        "quit\n"
    )
    completed = subprocess.run(
        [executable, "-f", "-q", "-s", "256M"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if not lines:
        raise ArithmeticError(f"system gp returned no ellrank result: {completed.stderr}")
    return pari(lines[-1]), executable


def probe(finalist, args):
    exact = finalist["exact_specialization"]
    model = tuple(Q(value) for value in exact["short_integral_model"])
    if any(value.denominator != 1 for value in model):
        raise ArithmeticError("orbit-103 search model is not integral")
    known = tuple(parse_point(point) for point in exact["known_points_Q_plus_Q_minus"])
    if len(known) != 2 or any(not is_on_weierstrass_curve(model, point) for point in known):
        raise ArithmeticError("Q_plus,Q_minus failed exact specialization")
    baseline_certificate = certificate(
        model, known, args.relation_prime, args.reduction_prime_bound
    )
    if not baseline_certificate["certified_independent"]:
        raise ArithmeticError("Q_plus,Q_minus lost finite-quotient independence")

    started = perf_counter()
    raw_curve = EllipticCurve(QQ, [QQ(value) for value in model])
    minimal_curve = raw_curve.global_minimal_model()
    to_minimal = raw_curve.isomorphism_to(minimal_curve)
    from_minimal = ~to_minimal
    minimal_model = tuple(
        Q(int(value.numerator()), int(value.denominator()))
        for value in minimal_curve.ainvs()
    )
    minimal_known = tuple(
        sage_point_tuple(to_minimal(raw_curve([QQ(point[0]), QQ(point[1])])))
        for point in known
    )
    if any(value.denominator != 1 for value in minimal_model):
        raise ArithmeticError("global minimal model is unexpectedly nonintegral")
    curve = mwrank_EllipticCurve([ZZ(value.numerator) for value in minimal_model])
    subgroup = mwrank_MordellWeil(curve, verbose=False, pp=1, maxr=args.max_rank)
    growth = []
    for label, point in zip(("Q_plus", "Q_minus"), minimal_known):
        before = len(subgroup.points())
        subgroup.process([to_mwrank_triple(point)], saturation_bound=0)
        after = len(subgroup.points())
        growth.append({"label": label, "before": before, "after": after})
    baseline_eclib_rank = len(subgroup.points())
    if baseline_eclib_rank != 2:
        raise ArithmeticError("eclib did not replay the rank-two baseline")

    search_started = perf_counter()
    if args.height > 0:
        subgroup.search(args.height, verbose=False)
    search_seconds = perf_counter() - search_started
    projective_points = tuple(tuple(map(int, point)) for point in subgroup.points())
    minimal_affine_points = tuple(
        point
        for triple in projective_points
        if (point := from_mwrank_triple(triple)) is not None
    )
    if any(
        not is_on_weierstrass_curve(minimal_model, point)
        for point in minimal_affine_points
    ):
        raise ArithmeticError("eclib returned a point off the exact minimal model")
    affine_points = tuple(
        sage_point_tuple(
            from_minimal(minimal_curve([QQ(point[0]), QQ(point[1])]))
        )
        for point in minimal_affine_points
    )
    if any(not is_on_weierstrass_curve(model, point) for point in affine_points):
        raise ArithmeticError("a point failed transport back to the exact short model")

    pari_started = perf_counter()
    pari_primary_error = None
    if args.force_system_gp:
        pari_result, gp_executable = external_gp_ellrank(
            minimal_model, minimal_known, args.pari_effort, args.gp_timeout
        )
        pari_backend = f"system_gp_fallback:{gp_executable}"
    else:
        pari_backend = "sage_pari_2.17"
        try:
            pari_context = minimal_curve.pari_curve().ellrankinit()
            pari_result = pari_context.ellrank(
                args.pari_effort,
                pari([[QQ(point[0]), QQ(point[1])] for point in minimal_known]),
            )
        except Exception as error:
            pari_primary_error = f"{type(error).__name__}: {error}"
            pari_result, gp_executable = external_gp_ellrank(
                minimal_model, minimal_known, args.pari_effort, args.gp_timeout
            )
            pari_backend = f"system_gp_fallback:{gp_executable}"
    pari_seconds = perf_counter() - pari_started
    pari_minimal_points = tuple(
        (Q(str(point[0])), Q(str(point[1]))) for point in pari_result[3]
    )
    if any(
        not is_on_weierstrass_curve(minimal_model, point)
        for point in pari_minimal_points
    ):
        raise ArithmeticError("PARI returned a point off the exact minimal model")
    pari_raw_points = tuple(
        sage_point_tuple(
            from_minimal(minimal_curve([QQ(point[0]), QQ(point[1])]))
        )
        for point in pari_minimal_points
    )
    if any(not is_on_weierstrass_curve(model, point) for point in pari_raw_points):
        raise ArithmeticError("a PARI point failed transport to the short model")

    known_keys = {sign_key(point) for point in known}
    candidates = []
    seen = set(known_keys)
    for point in (*affine_points, *pari_raw_points):
        key = sign_key(point)
        if key in seen:
            continue
        seen.add(key)
        candidates.append(point)
    candidates.sort(
        key=lambda point: max(
            abs(point[0].numerator).bit_length(),
            point[0].denominator.bit_length(),
            abs(point[1].numerator).bit_length(),
            point[1].denominator.bit_length(),
        )
    )

    selected = list(known)
    accepted = []
    rejected = []
    final_certificate = baseline_certificate
    for point in candidates:
        proposed = (*selected, point)
        proposed_certificate = certificate(
            model, proposed, args.relation_prime, args.reduction_prime_bound
        )
        if proposed_certificate["certified_independent"]:
            selected.append(point)
            accepted.append(point)
            final_certificate = proposed_certificate
        else:
            rejected.append(point)

    return {
        "k": finalist["k"],
        "r": finalist["r"],
        "k_projective_pair": finalist["k_projective_pair"],
        "r_projective_pair": finalist["r_projective_pair"],
        "nagao": {
            "selection_weakest_block": finalist["selection_weakest_block"],
            "confirmation_weakest_block": finalist["confirmation_weakest_block"],
        },
        "short_integral_model": [str(value) for value in model],
        "global_minimal_model": [str(value) for value in minimal_model],
        "minimalizing_isomorphism": str(to_minimal),
        "baseline_points": [point_record(point) for point in known],
        "baseline_points_on_minimal_model": [
            point_record(point) for point in minimal_known
        ],
        "baseline_finite_quotient_certificate": baseline_certificate,
        "eclib": {
            "height_limit": args.height,
            "max_rank": args.max_rank,
            "baseline_growth": growth,
            "rank_after_search": len(projective_points),
            "projective_points": [list(point) for point in projective_points],
            "search_seconds": search_seconds,
        },
        "pari_ellrank": {
            "backend": pari_backend,
            "effort": args.pari_effort,
            "primary_error": pari_primary_error,
            "reported_lower_bound": int(pari_result[0]),
            "reported_upper_bound": int(pari_result[1]),
            "sha_pairing_rank": int(pari_result[2]),
            "returned_points": len(pari_minimal_points),
            "seconds": pari_seconds,
            "use_in_claim": "returned points only; interval is diagnostic",
        },
        "exact_candidates_beyond_baseline": [point_record(point) for point in candidates],
        "finite_quotient_independent_new_points": [point_record(point) for point in accepted],
        "finite_quotient_rejected_points": [point_record(point) for point in rejected],
        "certified_rank_lower_bound": len(selected),
        "certified_gain_beyond_Q_plus_Q_minus": len(accepted),
        "combined_finite_quotient_certificate": final_certificate,
        "all_returned_affine_points_exactly_verified": True,
        "runtime_seconds": perf_counter() - started,
        "claim_boundary": (
            "The displayed lower bound is unconditional when the combined certificate "
            "is certified independent. Rejected candidates and failure to find points "
            "give no rank upper bound."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--candidate-limit", type=int, default=20)
    parser.add_argument("--candidate-start", type=int, default=1)
    parser.add_argument(
        "--lane",
        choices=("finalists", "small_coefficient_finalists"),
        default="finalists",
    )
    parser.add_argument("--height", type=float, default=15.0)
    parser.add_argument("--pari-effort", type=int, default=0)
    parser.add_argument("--gp-timeout", type=float, default=30.0)
    parser.add_argument("--force-system-gp", action="store_true")
    parser.add_argument("--max-rank", type=int, default=20)
    # Mod 2 can collapse the images of Q_plus,Q_minus on otherwise rank-two
    # fibres; mod 3 is the cheap stable baseline in the orbit-103 smoke suite.
    parser.add_argument("--relation-prime", type=int, default=3)
    parser.add_argument("--reduction-prime-bound", type=int, default=500)
    args = parser.parse_args()
    if min(
        args.candidate_start,
        args.candidate_limit,
        args.max_rank,
        args.reduction_prime_bound,
    ) <= 0:
        raise SystemExit("candidate, rank, and finite-quotient bounds must be positive")
    if args.height < 0 or args.pari_effort < 0 or args.gp_timeout <= 0:
        raise SystemExit("search height and PARI effort must be nonnegative; GP timeout positive")
    source_bytes = args.input.read_bytes()
    source = json.loads(source_bytes)
    if source.get("status") != "PASS_BOUNDED_TWO_STAGE_NAGAO_SEARCH":
        raise SystemExit("input is not an orbit-103 two-stage search artifact")

    rows = []
    errors = []
    started = perf_counter()
    first_index = args.candidate_start - 1
    last_index = first_index + args.candidate_limit
    for index, finalist in enumerate(
        source[args.lane][first_index:last_index], start=args.candidate_start
    ):
        print(
            f"E6A1O103PROBE|index={index}|k={finalist['k']}|r={finalist['r']}|stage=start",
            flush=True,
        )
        try:
            row = probe(finalist, args)
        except Exception as error:
            errors.append(
                {
                    "index": index,
                    "k": finalist["k"],
                    "r": finalist["r"],
                    "error_type": type(error).__name__,
                    "error": str(error),
                }
            )
            print(
                f"E6A1O103PROBE|index={index}|stage=error|type={type(error).__name__}|message={error}",
                flush=True,
            )
            continue
        rows.append(row)
        print(
            f"E6A1O103PROBE|index={index}|stage=complete|"
            f"eclib_rank={row['eclib']['rank_after_search']}|"
            f"certified_rank={row['certified_rank_lower_bound']}|"
            f"gain={row['certified_gain_beyond_Q_plus_Q_minus']}|"
            f"seconds={row['runtime_seconds']:.3f}",
            flush=True,
        )

    best_rank = max((row["certified_rank_lower_bound"] for row in rows), default=2)
    payload = {
        "schema": "elkies-k3.e6a1-orbit103-specialization-rank-probes.v1",
        "status": "PASS_BOUNDED_EXACT_POINT_AND_FINITE_QUOTIENT_PROBES",
        "input": {
            "path": str(args.input),
            "sha256": sha256(source_bytes).hexdigest(),
        },
        "sage_version": str(sage_version()),
        "bounds": {
            "candidate_limit": args.candidate_limit,
            "candidate_start": args.candidate_start,
            "candidate_lane": args.lane,
            "mwrank_logarithmic_height": args.height,
            "mwrank_max_rank": args.max_rank,
            "pari_ellrank_effort": args.pari_effort,
            "system_gp_timeout_seconds": args.gp_timeout,
            "force_system_gp": args.force_system_gp,
            "finite_quotient_relation_prime": args.relation_prime,
            "finite_quotient_reduction_prime_bound": args.reduction_prime_bound,
        },
        "best_certified_rank_lower_bound": best_rank,
        "fibres_of_certified_rank_at_least_8": sum(
            row["certified_rank_lower_bound"] >= 8 for row in rows
        ),
        "records": rows,
        "errors": errors,
        "runtime_seconds": perf_counter() - started,
        "proof_boundary": (
            "Exact model and point identities plus a full finite-quotient certificate "
            "prove each displayed rank lower bound. The bounded eclib searches do not "
            "supply upper bounds, and a Nagao score alone is never rank evidence."
        ),
        "reproducing_command": sys.argv,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"E6A1O103PROBE|records={len(rows)}|errors={len(errors)}|"
        f"best_rank={best_rank}|seconds={payload['runtime_seconds']:.3f}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
