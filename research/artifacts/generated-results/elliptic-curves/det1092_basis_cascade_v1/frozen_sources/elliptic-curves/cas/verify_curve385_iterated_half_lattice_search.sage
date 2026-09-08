#!/usr/bin/env sage -python
"""Verify the frozen blind curve-385 iteration against the public 29 points.

The blind ledger is read and hashed before this program opens the public
fixture.  Mutual exact integral coordinate matrices then prove equality of the
blindly discovered rank-29 subgroup and the displayed public rank-29 subgroup.
This verifier cannot influence chart selection or the frozen point search.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys
from typing import Any, Sequence

from sage.all import EllipticCurve, Matrix, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
BLIND = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "curve385_iterated_half_lattice_blind_v1.json"
)
PUBLIC = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "curve385_iterated_half_lattice_verification_v1.json"
)

EXPECTED_BLIND_SHA256 = "356001898f738f607d984e081663a015825e11de0c606d35055af156eb2d7502"
EXPECTED_PUBLIC_SHA256 = "9a2675ab48cc37111d1f4050bd1797fc84c98b7839668d292d11406efe7a9eaa"
DIMENSION = 29

sys.path[:0] = [str(ROOT / "elliptic-curves"), str(CAS)]

from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)
from search_nagao_u135_alternate_covers import relation_proposals  # noqa: E402


Point = tuple[Fraction, Fraction]


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def point_record(point: Point) -> dict[str, str]:
    return {"x": str(point[0]), "y": str(point[1])}


def read_point(record: dict[str, str]) -> Point:
    return Fraction(record["x"]), Fraction(record["y"])


def short_data(ainvs, points):
    a1, a2, a3, a4, a6 = (QQ(value) for value in ainvs)
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    c4 = b2 * b2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    model = (QQ(0), QQ(0), QQ(0), -27 * c4, -54 * c6)
    short_points = tuple(
        (
            36 * QQ(point[0]) + 3 * b2,
            108 * (2 * QQ(point[1]) + a1 * QQ(point[0]) + a3),
        )
        for point in points
    )
    return tuple(Fraction(str(value)) for value in model), tuple(
        (Fraction(str(point[0])), Fraction(str(point[1]))) for point in short_points
    )


def exact_relations(
    model,
    basis: Sequence[Point],
    points: Sequence[Point],
    *,
    chunk_size: int,
    timeout_seconds: float,
    stack_bytes: int,
):
    answer = []
    for start in range(0, len(points), chunk_size):
        answer.extend(
            relation_proposals(
                model,
                basis,
                points[start : start + chunk_size],
                timeout=timeout_seconds,
                stack_bytes=stack_bytes,
            )
        )
    if not all(exact for unused_relation, exact in answer):
        failures = [index for index, (unused_relation, exact) in enumerate(answer) if not exact]
        raise ArithmeticError(f"exact integral relations failed at indices {failures}")
    return tuple(tuple(map(int, relation)) for relation, unused_exact in answer)


def finite_certificate(model, points):
    signatures = find_mod2_reduction_certificate(model, points, prime_bound=1_000)
    rank = combined_mod2_rank(signatures, len(points))
    return {
        "prime_bound": 1_000,
        "column_count": len(points),
        "combined_rank": rank,
        "signatures": [
            {
                "prime": row.prime,
                "group_order": row.group_order,
                "doubled_subgroup_order": row.doubled_subgroup_order,
                "quotient_dimension": row.quotient_dimension,
                "rows": [list(vector) for vector in row.rows],
            }
            for row in signatures
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blind", type=Path, default=BLIND)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--relation-chunk-size", type=int, default=64)
    parser.add_argument("--relation-timeout-seconds", type=float, default=180.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    args = parser.parse_args()
    if args.relation_chunk_size <= 0 or not 0 < args.relation_timeout_seconds <= 300:
        raise SystemExit("invalid relation budget")

    # This read-and-hash is the phase boundary.  Do not move public imports or
    # reads above it.
    blind_bytes = args.blind.read_bytes()
    blind_hash = sha256(blind_bytes).hexdigest()
    if blind_hash != EXPECTED_BLIND_SHA256:
        raise ArithmeticError("the frozen blind iterative ledger changed")
    blind = json.loads(blind_bytes)
    if blind.get("status") != "STOPPED_AT_DECLARED_LIFT_LIMIT":
        raise ArithmeticError("the blind run did not end at its declared lift safeguard")
    if blind["blindness_boundary"]["public_rank29_fixture_loaded"] is not False:
        raise ArithmeticError("the blind ledger crossed its public-fixture boundary")
    if len(blind["current_basis"]) != DIMENSION:
        raise ArithmeticError("the blind ledger did not reach rank 29")
    blind_model = tuple(Fraction(value) for value in blind["curve"]["short_model"])
    blind_basis = tuple(read_point(row) for row in blind["current_basis"])
    generic = tuple(read_point(row) for row in blind["curve"]["generic_points"])
    discoveries = tuple(read_point(row["point"]) for row in blind["discoveries"])

    # The public fixture is first read after the blind bytes have been frozen
    # and authenticated above.
    if digest(PUBLIC) != EXPECTED_PUBLIC_SHA256:
        raise ArithmeticError("the public curve fixture changed")
    public = json.loads(PUBLIC.read_text())
    record = next(row for row in public["records"] if row["id"] == 385)
    public_model, public_basis = short_data(record["ainvs"], record["points"])
    if public_model != blind_model:
        raise ArithmeticError("the blind and public short models differ")
    if public_basis[:17] != generic:
        raise ArithmeticError("the public first seventeen points differ from the blind generic subgroup")
    curve = EllipticCurve(QQ, list(public_model))
    if any(curve(point) == curve(0) for point in public_basis + blind_basis):
        raise ArithmeticError("a displayed basis contains the point at infinity")

    blind_certificate = finite_certificate(blind_model, blind_basis)
    public_certificate = finite_certificate(public_model, public_basis)
    if blind_certificate["combined_rank"] != DIMENSION:
        raise ArithmeticError("the blind subgroup failed exact rank-29 certification")
    if public_certificate["combined_rank"] != DIMENSION:
        raise ArithmeticError("the public subgroup failed exact rank-29 certification")

    blind_in_public = exact_relations(
        blind_model,
        public_basis,
        blind_basis,
        chunk_size=args.relation_chunk_size,
        timeout_seconds=args.relation_timeout_seconds,
        stack_bytes=args.stack_bytes,
    )
    public_in_blind = exact_relations(
        blind_model,
        blind_basis,
        public_basis,
        chunk_size=args.relation_chunk_size,
        timeout_seconds=args.relation_timeout_seconds,
        stack_bytes=args.stack_bytes,
    )
    blind_to_public = Matrix(ZZ, blind_in_public)
    public_to_blind = Matrix(ZZ, public_in_blind)
    identity = Matrix.identity(ZZ, DIMENSION)
    if blind_to_public * public_to_blind != identity:
        raise ArithmeticError("blind-to-public coordinate matrices are not inverse")
    if public_to_blind * blind_to_public != identity:
        raise ArithmeticError("public-to-blind coordinate matrices are not inverse")
    if abs(blind_to_public.det()) != 1 or abs(public_to_blind.det()) != 1:
        raise ArithmeticError("mutual integral containment is not unimodular")

    discovery_relations = exact_relations(
        blind_model,
        blind_basis,
        discoveries,
        chunk_size=args.relation_chunk_size,
        timeout_seconds=args.relation_timeout_seconds,
        stack_bytes=args.stack_bytes,
    )
    initial_new = blind["initial_transition"]["selected_new_direction_count"]
    iteration_new = sum(
        row["new_independent_direction_count"] for row in blind["iterations"]
    )
    if initial_new != 3 or iteration_new != 9:
        raise ArithmeticError("the blind 3+9 transition inventory changed")

    payload = {
        "schema": "elliptic-curves.curve385-iterated-half-lattice-verification.v1",
        "status": "PASS_BLIND_M29_EQUALS_DISPLAYED_PUBLIC_M29",
        "phase_boundary": {
            "blind_artifact_sha256_before_public_fixture_import": blind_hash,
            "public_fixture_loaded_only_after_blind_hash_check": True,
            "public_points_could_not_influence_search_or_selection": True,
        },
        "transition": {
            "generic_rank": 17,
            "initial_deep43_blind_gain": initial_new,
            "first_lift_round_blind_gain": iteration_new,
            "final_blind_rank": 29,
            "displayed_public_rank": 29,
            "displayed_public_quotient_dimension": 12,
            "blind_quotient_dimension_recovered": 12,
            "unrecovered_displayed_public_quotient_dimension": 0,
        },
        "subgroup_equality": {
            "blind_basis_in_public_basis_rows": [list(row) for row in blind_in_public],
            "public_basis_in_blind_basis_rows": [list(row) for row in public_in_blind],
            "coordinate_matrices_are_mutual_inverses": True,
            "blind_to_public_determinant": int(blind_to_public.det()),
            "public_to_blind_determinant": int(public_to_blind.det()),
        },
        "independence": {
            "blind_rank29": blind_certificate,
            "public_rank29": public_certificate,
        },
        "all_blind_discoveries": {
            "point_count": len(discoveries),
            "all_have_exact_integral_relations_in_final_blind_basis": True,
            "coordinate_rows_sha256": canonical_hash([list(row) for row in discovery_relations]),
        },
        "bounded_search_audit": {
            "searched_chart_count": sum(
                row["searched_new_chart_count"] for row in blind["iterations"]
            ),
            "bounded_complete_count": sum(
                row["bounded_complete_count"] for row in blind["iterations"]
            ),
            "timeout_count": sum(row["timeout_count"] for row in blind["iterations"]),
            "pari_failure_count": sum(
                row["pari_failure_count"] for row in blind["iterations"]
            ),
            "blind_terminal_status": blind["status"],
            "next_blind_lift_count": blind["stop"]["next_nonzero_lift_count"],
        },
        "input_hashes": {
            relative(args.blind): blind_hash,
            relative(PUBLIC): digest(PUBLIC),
            relative(Path(__file__)): digest(Path(__file__)),
        },
        "generation": {"python": platform.python_version()},
        "claim_boundary": [
            "The blind and displayed public rank-29 subgroups are exactly equal as integral subgroups.",
            "This proves recovery of all twelve displayed public quotient directions while preserving the search-phase blind boundary.",
            "The public 29-point subgroup is a rank lower bound; equality with the unknown full Mordell--Weil group is not proved.",
            "The blind search stopped at its declared lift limit and did not establish bounded-search stability at rank 29.",
        ],
        "reproducing_command": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elliptic-curves/cas/verify_curve385_iterated_half_lattice_search.sage"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "C385ITERVERIFY|status=PASS|blind_rank=29|public_rank=29|"
        f"det={blind_to_public.det()}|discoveries={len(discoveries)}|"
        f"output={relative(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
