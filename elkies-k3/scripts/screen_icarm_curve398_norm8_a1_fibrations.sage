#!/usr/bin/env sage-python
"""Screen one complete old-degree-two A1/MW16 layer for ICARM curve 398.

For a rootless determinant-948 fibration with saturated Mordell--Weil
lattice ``M``, every minimum-norm-eight parity class represented by ``w``
defines the isotropic class

    D_w = (2, 2, w),  D_w^2 = 0,  D_w.F = 2,  D_w.O = 0.

The old zero is therefore vertical in the new genus-one pencil.  This is the
first complete lattice layer in which the hidden I2/III, generic-MW16
fibration used for ICARM curve 398 can occur.

For each committed trace this program constructs the residual-chord quartic

    q(u, lambda)

over a declared finite field, computes its binary-quartic Jacobian invariants,
and tests the exact necessary condition

    j(Jac(q_lambda)) = j(curve 398)

for every ``lambda`` in P1(F_p).  A trace with no projective root is excluded
from having a rational curve-398 parameter in characteristic zero.  Survivors
remain UNKNOWN until an exact QQ compilation and rational-root factorization.

The program is deliberately chunked.  Each output pins the model, complete
priority table, target source, prime, and processed half-open rank interval.
It is a modular exclusion certificate, not a fibration or section recovery.
"""

from __future__ import annotations

import argparse
import csv
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
TARGET_SOURCE = ROOT / "elliptic-curves/cas/icarm_curve398.py"
CHORD_SOURCE = SCRIPTS / "construct_elkies_2026_bisections.sage"


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


def parse_rational(value) -> Fraction:
    return Fraction(str(value))


def reduce_rational(value, field):
    value = parse_rational(value)
    return field(value.numerator) / field(value.denominator)


def polynomial_from_record(record, ring, field):
    numerator = ring(
        [reduce_rational(value, field) for value in record["numerator_coefficients_low_to_high"]]
    )
    denominator = ring(
        [reduce_rational(value, field) for value in record["denominator_coefficients_low_to_high"]]
    )
    return ring.fraction_field()(numerator) / ring.fraction_field()(denominator)


def parse_vector(text: str) -> tuple[int, ...]:
    return tuple(int(value) for value in text.split())


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as source:
        rows = list(csv.DictReader(source, delimiter="\t"))
    if not rows or "section_basis_w" not in rows[0]:
        raise ValueError("norm-eight priority table has no section_basis_w column")
    for expected_rank, row in enumerate(rows, start=1):
        if int(row["priority_rank"]) != expected_rank:
            raise ValueError("norm-eight priority table is not rank-contiguous")
    return rows


def binary_quartic_invariants(quartic, coefficient_ring):
    coefficients = [coefficient_ring(quartic[index]) for index in range(5)]
    e, d, c, b, a = coefficients
    invariant_i = 12 * a * e - 3 * b * d + c**2
    invariant_j = (
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d**2
        - 27 * b**2 * e
        - 2 * c**3
    )
    return invariant_i, invariant_j


def projective_target_roots(invariant_i, invariant_j, target_a, target_b, field):
    # The standard binary-quartic Jacobian is
    # y^2=x^3-27*I*x-27*J.  Cross multiplication by the target short model
    # avoids division at j=0, 1728, and at finite bad pencil parameters.
    child_a = -27 * invariant_i
    child_b = -27 * invariant_j
    comparison = child_a**3 * target_b**2 - target_a**3 * child_b**2
    if not comparison:
        return None, None, "IDENTICALLY_ZERO_COMPARISON"
    finite_roots = tuple(int(root) for root in field if comparison(root) == 0)

    # At lambda=infinity compare the rational ratio A(lambda)^3/B(lambda)^2.
    # A finite nonzero target ratio can occur there only when the two degrees
    # agree; this test is invariant under a common nonminimal gauge factor.
    numerator = child_a**3
    denominator = child_b**2
    infinity = bool(
        numerator
        and denominator
        and numerator.degree() == denominator.degree()
        and numerator.leading_coefficient() * target_b**2
        == target_a**3 * denominator.leading_coefficient()
    )
    roots = list(finite_roots)
    if infinity:
        roots.append("infinity")
    return comparison, roots, None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--start-rank", type=int, default=1)
    parser.add_argument("--limit", type=int)
    parser.add_argument(
        "--candidate-ranks",
        type=Path,
        help="optional JSON list of one-based priority ranks retained by an earlier prime",
    )
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
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
    output_path = arguments.output.resolve()
    model = json.loads(model_path.read_text())
    if model.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
        raise ValueError("expected a certified direct rootless equation")
    if model["sections"].get("status") != "PASS_EXACT_SATURATED_RANK17_BASIS":
        raise ValueError("rootless equation does not carry a saturated MW17 basis")

    rows = load_rows(table_path)
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
        requested_document = json.loads(candidate_path.read_text())
        requested = (
            requested_document
            if isinstance(requested_document, list)
            else requested_document.get("search", {}).get("survivor_priority_ranks")
        )
        if not isinstance(requested, list) or any(not isinstance(value, int) for value in requested):
            raise ValueError(
                "--candidate-ranks must contain a JSON integer list or an earlier screen ledger"
            )
        all_ranks &= set(requested)
    selected_rows = [rows[rank - 1] for rank in sorted(all_ranks)]

    target = load_source("icarm_curve398_norm8_target", TARGET_SOURCE)
    a1, a2, a3, a4, a6 = tuple(parse_rational(value) for value in target.GENERAL_WEIERSTRASS_COEFFICIENTS)
    b2 = a1**2 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3**2 + 4 * a6
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    target_a = -27 * reduce_rational(c4, field)
    target_b = -54 * reduce_rational(c6, field)
    if not target_a or not target_b:
        raise ArithmeticError("chosen prime degenerates the target short j-ratio")

    ring = PolynomialRing(field, "u")
    function_field = ring.fraction_field()
    u = ring.gen()
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
            polynomial_from_record(record["X"], ring, field),
            polynomial_from_record(record["Y"], ring, field),
        )
        for record in model["sections"]["records"]
    )
    if len(basis) != 17:
        raise ArithmeticError("rootless equation basis no longer has rank seventeen")

    chord = load_source("icarm_curve398_norm8_chord", CHORD_SOURCE)
    lambda_ring = PolynomialRing(field, "lambda")
    lambda_variable = lambda_ring.gen()
    bivariate_ring = PolynomialRing(lambda_ring, "u")

    started = time.monotonic()
    records = []
    survivors = []
    status_histogram: dict[str, int] = {}
    for position, row in enumerate(selected_rows, start=1):
        rank = int(row["priority_rank"])
        coordinates = parse_vector(row["section_basis_w"])
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
                invariant_i, invariant_j = binary_quartic_invariants(quartic, lambda_ring)
                comparison, roots, error = projective_target_roots(
                    invariant_i, invariant_j, target_a, target_b, field
                )
                if error is not None:
                    record.update({"status": error, "projective_roots": None})
                    survivors.append(rank)
                else:
                    status = (
                        "SURVIVES_MODULAR_J_TEST"
                        if roots
                        else "PASS_MODULAR_NO_CURVE398_PARAMETER"
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
            # Bad reduction cannot exclude a characteristic-zero candidate.
            record.update(
                {
                    "status": "UNAVAILABLE_BAD_REDUCTION",
                    "reason": str(error),
                }
            )
            survivors.append(rank)
        status_histogram[record["status"]] = status_histogram.get(record["status"], 0) + 1
        records.append(record)
        if position % 1000 == 0:
            print(
                f"CURVE398A1SCREEN|prime={prime}|done={position}/{len(selected_rows)}"
                f"|survivors={len(survivors)}|seconds={time.monotonic()-started:.1f}",
                flush=True,
            )

    payload = {
        "schema": "elkies-k3.icarm-curve398-norm8-a1-modular-screen.v1",
        "status": "PASS_COMPLETE_DECLARED_CHUNK_MODULAR_SCREEN",
        "source_fibration": model.get("divisor", {}).get("label", model_path.stem),
        "target": {
            "label": "ICARM curve 398",
            "generalized_weierstrass_coefficients": [str(value) for value in (a1, a2, a3, a4, a6)],
            "c4": str(c4),
            "c6": str(c6),
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
            "Every PASS_MODULAR_NO_CURVE398_PARAMETER row has no parameter in P1(F_p), "
            "so that trace cannot have a rational characteristic-zero curve-398 parameter. "
            "Survivors and unavailable reductions remain UNKNOWN. Even a complete no-hit "
            "excludes only this committed old-degree-two norm-eight layer on this source chart."
        ),
        "inputs": {
            relative(path): digest(path)
            for path in (model_path, table_path, TARGET_SOURCE, CHORD_SOURCE)
        },
        "software": {"sage_version": SAGE_VERSION},
        "reproducing_command": (
            "sage -python elkies-k3/scripts/screen_icarm_curve398_norm8_a1_fibrations.sage "
            f"--model {relative(model_path)} --table {relative(table_path)} "
            f"--prime {prime} --start-rank {arguments.start_rank} "
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
        f"CURVE398A1SCREEN|prime={prime}|processed={len(records)}"
        f"|excluded={len(records)-len(survivors)}|survivors={len(survivors)}"
        f"|output={relative(output_path)}|status={payload['status']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
