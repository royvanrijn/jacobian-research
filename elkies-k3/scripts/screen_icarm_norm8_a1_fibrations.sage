#!/usr/bin/env sage-python
"""Screen the complete alternate-Q80 norm-eight A1/MW16 layer for one ICARM curve.

This is the target-generic form of the frozen curve-398 screen.  The source is
the equation-explicit ``norm12-orbit-11952`` rootless/MW17 fibration.  Every
minimum-norm-eight trace ``w`` gives the old-degree-two isotropic class

    D_w = (2, 2, w) = O + P_w,

whose Jacobian is an A1/MW16 pencil when the finite-pole chord chart is
available.  At one declared prime this program tests the necessary projective
j-equation for a pinned ICARM target.  A no-root result excludes that
trace/target pair over QQ; survivors remain UNKNOWN until exact factorization.

The program is checkpointed by prime and accepts the survivor list from a
previous checkpoint.  The target equation is read from the exact 573-curve
snapshot embedded in the certified norm-twelve database sweep, so no live
database access is involved.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
import json
from pathlib import Path
import sys
import time

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
DEFAULT_MODEL = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
DEFAULT_TABLE = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
)
DEFAULT_TARGET_SNAPSHOT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-icarm-database-sweep-v2.json"
)
CHORD_SOURCE = SCRIPTS / "construct_elkies_2026_bisections.sage"
FROZEN_SCREEN_SOURCE = SCRIPTS / "screen_icarm_curve398_norm8_a1_fibrations.sage"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def load_source(name: str, path: Path):
    loader = SourceFileLoader(name, str(path))
    spec = spec_from_loader(name, loader)
    if spec is None:
        raise ImportError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    loader.exec_module(module)
    return module


def normalized_ainvs(values) -> tuple[Fraction, ...]:
    if len(values) == 2:
        values = [0, 0, 0, *values]
    if len(values) != 5:
        raise ValueError("target curve has neither two nor five a-invariants")
    return tuple(Fraction(str(value)) for value in values)


def target_record(snapshot: dict, curve_id: int) -> dict:
    records = snapshot.get("snapshot", {}).get("curves")
    if not isinstance(records, list):
        raise ValueError("target snapshot has no pinned curve equation list")
    matches = [record for record in records if int(record.get("id", -1)) == curve_id]
    if len(matches) != 1:
        raise ValueError(f"target snapshot has {len(matches)} rows for curve {curve_id}")
    return matches[0]


def invariants(values) -> tuple[tuple[Fraction, ...], Fraction, Fraction, Fraction]:
    a1, a2, a3, a4, a6 = normalized_ainvs(values)
    b2 = a1**2 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    delta = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    if not delta or c4**3 - c6**2 != 1728 * delta:
        raise ArithmeticError("target curve is singular or has inconsistent invariants")
    return (a1, a2, a3, a4, a6), c4, c6, delta


def reduce_rational(value, field):
    value = Fraction(str(value))
    return field(value.numerator) / field(value.denominator)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--target-snapshot", type=Path, default=DEFAULT_TARGET_SNAPSHOT)
    parser.add_argument("--curve-id", type=int, required=True)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--start-rank", type=int, default=1)
    parser.add_argument("--limit", type=int)
    parser.add_argument(
        "--candidate-ranks",
        type=Path,
        help="optional JSON integer list or earlier screen/compact survivor ledger",
    )
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.curve_id <= 0:
        parser.error("--curve-id must be positive")
    if arguments.start_rank <= 0:
        parser.error("--start-rank must be positive")
    if arguments.limit is not None and arguments.limit <= 0:
        parser.error("--limit must be positive")

    prime = int(arguments.prime)
    if prime in (2, 3) or not ZZ(prime).is_prime():
        parser.error("--prime must be a prime greater than three")
    field = GF(prime)

    model_path = arguments.model.resolve()
    table_path = arguments.table.resolve()
    snapshot_path = arguments.target_snapshot.resolve()
    output_path = arguments.output.resolve()
    model = json.loads(model_path.read_text())
    if model.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
        raise ValueError("expected a certified direct rootless equation")
    if model["sections"].get("status") != "PASS_EXACT_SATURATED_RANK17_BASIS":
        raise ValueError("rootless equation does not carry a saturated MW17 basis")

    helpers = load_source("icarm_norm8_a1_frozen_helpers", FROZEN_SCREEN_SOURCE)
    chord = load_source("icarm_norm8_a1_chord", CHORD_SOURCE)
    rows = helpers.load_rows(table_path)
    all_ranks = set(range(arguments.start_rank, len(rows) + 1))
    if arguments.limit is not None:
        all_ranks = set(
            range(
                arguments.start_rank,
                min(len(rows) + 1, arguments.start_rank + arguments.limit),
            )
        )
    candidate_input_hash = None
    if arguments.candidate_ranks is not None:
        candidate_path = arguments.candidate_ranks.resolve()
        candidate_input_hash = digest(candidate_path)
        candidate_document = json.loads(candidate_path.read_text())
        requested = (
            candidate_document
            if isinstance(candidate_document, list)
            else candidate_document.get("search", {}).get("survivor_priority_ranks")
        )
        if not isinstance(requested, list) or any(not isinstance(value, int) for value in requested):
            raise ValueError("candidate ledger has no integer survivor_priority_ranks list")
        candidate_target = candidate_document.get("target") if isinstance(candidate_document, dict) else None
        prior_curve_id = candidate_target.get("curve_id") if isinstance(candidate_target, dict) else None
        if prior_curve_id is not None and int(prior_curve_id) != arguments.curve_id:
            raise ValueError("candidate ledger belongs to a different ICARM curve")
        all_ranks &= set(requested)
    selected_rows = [rows[rank - 1] for rank in sorted(all_ranks)]

    snapshot = json.loads(snapshot_path.read_text())
    pinned_target = target_record(snapshot, arguments.curve_id)
    ainvs, c4, c6, delta = invariants(pinned_target["ainvs"])
    target_a = -27 * reduce_rational(c4, field)
    target_b = -54 * reduce_rational(c6, field)
    if not target_a or not target_b or reduce_rational(delta, field) == 0:
        raise ArithmeticError("chosen prime is not a usable target-reduction prime")

    ring = PolynomialRing(field, "u")
    function_field = ring.fraction_field()
    weierstrass = model["weierstrass_model"]
    old_a = ring(
        [reduce_rational(value, field) for value in weierstrass["A_coefficients_low_to_high"]]
    )
    old_b = ring(
        [reduce_rational(value, field) for value in weierstrass["B_coefficients_low_to_high"]]
    )
    old_curve = EllipticCurve(function_field, [old_a, old_b])
    basis = tuple(
        old_curve(
            helpers.polynomial_from_record(record["X"], ring, field),
            helpers.polynomial_from_record(record["Y"], ring, field),
        )
        for record in model["sections"]["records"]
    )
    if len(basis) != 17:
        raise ArithmeticError("rootless equation basis no longer has rank seventeen")

    lambda_ring = PolynomialRing(field, "lambda")
    lambda_variable = lambda_ring.gen()
    bivariate_ring = PolynomialRing(lambda_ring, "u")

    started = time.monotonic()
    records = []
    survivors = []
    status_histogram: dict[str, int] = {}
    for position, row in enumerate(selected_rows, start=1):
        rank = int(row["priority_rank"])
        coordinates = helpers.parse_vector(row["section_basis_w"])
        if len(coordinates) != 17:
            raise ArithmeticError(f"priority rank {rank} has a non-rank-17 trace word")
        record = {
            "priority_rank": rank,
            "orbit_mask": int(row["orbit_mask"]),
            "orbit_hex": row["orbit_hex"],
            "section_basis_w": list(coordinates),
        }
        try:
            trace = sum(
                (
                    coefficient * point
                    for coefficient, point in zip(coordinates, basis)
                    if coefficient
                ),
                old_curve(0),
            )
            if trace.is_zero():
                raise ArithmeticError("trace reduces to zero")
            frame = chord.trace_chord_frame(trace[0], trace[1], ring)
            h, nx, ny, m0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
            if h.degree() != 2:
                record.update(
                    {
                        "status": "UNAVAILABLE_POLE_AT_INFINITY_CHART",
                        "finite_pole_degree": int(h.degree()),
                    }
                )
                survivors.append(rank)
            else:
                h_symbolic = bivariate_ring(h)
                nx_symbolic = bivariate_ring(nx)
                ny_symbolic = bivariate_ring(ny)
                m_symbolic = bivariate_ring(m0) + lambda_variable * h_symbolic**2
                numerator = (
                    m_symbolic**4
                    - 6 * m_symbolic**2 * nx_symbolic
                    - 8 * m_symbolic * ny_symbolic
                    - 3 * nx_symbolic**2
                    - 4 * bivariate_ring(old_a) * h_symbolic**4
                )
                quartic, remainder = numerator.quo_rem(h_symbolic**6)
                if remainder or quartic.degree() != 4:
                    raise ArithmeticError("residual chord did not produce a binary quartic")
                invariant_i, invariant_j = helpers.binary_quartic_invariants(quartic, lambda_ring)
                comparison, roots, error = helpers.projective_target_roots(
                    invariant_i, invariant_j, target_a, target_b, field
                )
                if error is not None:
                    record.update({"status": error, "projective_roots": None})
                    survivors.append(rank)
                else:
                    status = (
                        "SURVIVES_MODULAR_J_TEST"
                        if roots
                        else "PASS_MODULAR_NO_TARGET_PARAMETER"
                    )
                    record.update(
                        {
                            "status": status,
                            "finite_pole_degree": int(h.degree()),
                            "invariant_i_degree": int(invariant_i.degree()),
                            "invariant_j_degree": int(invariant_j.degree()),
                            "comparison_degree": int(comparison.degree()),
                            "projective_roots": roots,
                        }
                    )
                    if roots:
                        survivors.append(rank)
        except (ArithmeticError, ValueError, ZeroDivisionError) as error:
            record.update({"status": "UNAVAILABLE_BAD_REDUCTION", "reason": str(error)})
            survivors.append(rank)
        status_histogram[record["status"]] = status_histogram.get(record["status"], 0) + 1
        records.append(record)
        if position % 1000 == 0:
            print(
                f"ICARMA1SCREEN|curve={arguments.curve_id}|prime={prime}"
                f"|done={position}/{len(selected_rows)}|survivors={len(survivors)}"
                f"|seconds={time.monotonic()-started:.1f}",
                flush=True,
            )

    payload = {
        "schema": "elkies-k3.icarm-norm8-a1-modular-screen.v1",
        "status": "PASS_COMPLETE_DECLARED_CHUNK_MODULAR_SCREEN",
        "source_fibration": model.get("divisor", {}).get("label", model_path.stem),
        "target": {
            "curve_id": arguments.curve_id,
            "label": f"ICARM curve {arguments.curve_id}",
            "snapshot_rank_lower_bound": int(pinned_target["snapshot_rank_lower_bound"]),
            "generalized_weierstrass_coefficients": [str(value) for value in ainvs],
            "c4": str(c4),
            "c6": str(c6),
            "discriminant": str(delta),
        },
        "search": {
            "prime": prime,
            "priority_table_class_count": len(rows),
            "start_rank": arguments.start_rank,
            "limit": arguments.limit,
            "candidate_rank_filter": (
                None if arguments.candidate_ranks is None else relative(arguments.candidate_ranks)
            ),
            "candidate_rank_filter_sha256": candidate_input_hash,
            "processed_count": len(records),
            "processed_priority_ranks": [record["priority_rank"] for record in records],
            "status_histogram": dict(sorted(status_histogram.items())),
            "survivor_count": len(survivors),
            "survivor_priority_ranks": survivors,
            "wall_seconds": time.monotonic() - started,
        },
        "records": records,
        "proof_boundary": (
            "Every PASS_MODULAR_NO_TARGET_PARAMETER row has no parameter in P1(F_p), "
            "so that trace cannot have a rational characteristic-zero target parameter. "
            "Survivors and unavailable reductions remain UNKNOWN. A complete no-hit excludes "
            "only this committed old-degree-two norm-eight layer on source chart 11952."
        ),
        "inputs": {
            relative(path): digest(path)
            for path in (
                model_path,
                table_path,
                snapshot_path,
                CHORD_SOURCE,
                FROZEN_SCREEN_SOURCE,
            )
        },
        "software": {"sage_version": SAGE_VERSION},
        "reproducing_command": (
            "sage -python elkies-k3/scripts/screen_icarm_norm8_a1_fibrations.sage "
            f"--curve-id {arguments.curve_id} --prime {prime} "
            f"--model {relative(model_path)} --table {relative(table_path)} "
            f"--target-snapshot {relative(snapshot_path)} --start-rank {arguments.start_rank} "
            + ("" if arguments.limit is None else f"--limit {arguments.limit} ")
            + (
                ""
                if arguments.candidate_ranks is None
                else f"--candidate-ranks {relative(arguments.candidate_ranks)} "
            )
            + f"--output {relative(output_path)}"
        ),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"ICARMA1SCREEN|curve={arguments.curve_id}|prime={prime}|processed={len(records)}"
        f"|excluded={len(records)-len(survivors)}|survivors={len(survivors)}"
        f"|output={relative(output_path)}|status={payload['status']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
