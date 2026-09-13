#!/usr/bin/env python3
"""Check retained native MW17 coordinates without repeating height discovery.

Run with Sage Python. Choose --transports-only (no cover inventories) or
--with-covers (the complete inventories, about 1.18 GB for all three charts).
--curve may be repeated to replay a bounded subset. Neither mode regenerates
the original certificate. The original producer and its dependent hashes stay
unchanged; its special_fibre function runs in an isolated namespace with exact
retained-coordinate verification replacing numerical coordinate recovery.

This is a replay of the selected native component, not all eight checkers in
the calibration note. Historical height-rounding metadata is not remeasured.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import runpy
import tempfile
from time import monotonic
from types import FunctionType


ROOT = Path(__file__).resolve().parents[2]
PRODUCER = ROOT / "elkies-k3/scripts/certify_r17_norm12_native_icarm_quotient_audit.sage"
CERTIFICATE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-native-icarm-quotient-audit-v1.json"
PRODUCER_SHA256 = "878981317bcc71f72aabfd5e88ce3051a629a55bfd793156305d764adf44516c"
CERTIFICATE_SHA256 = "582898590dcfdd8f36b7564d8df65dff951b0b304eddb6ae8825140eb1a92d9a"
CURVE_IDS = (12, 395, 363, 364, 378, 393, 404)


def load_pinned(path, expected):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise ArithmeticError(f"pinned input changed: {path}")
    return json.loads(data)


def integer_matrix(rows, height, width):
    """Reject truncation, bools and nonintegral coercions before group addition."""
    if (not isinstance(rows, list) or len(rows) != height
            or any(not isinstance(row, list) or len(row) != width for row in rows)
            or any(type(value) is not int for row in rows for value in row)):
        raise ArithmeticError("retained coordinates must be an integer matrix of the exact dimensions")
    return rows


def coordinate_rows(record, *, with_covers):
    count = record["displayed_point_count"]
    if type(count) is not int or count < 17:
        raise ArithmeticError("invalid displayed point count")
    generic = integer_matrix(record["specialized_generic_subgroup"][
        "coordinate_matrix_rows_in_ordered_public_points"], count, 17)
    rows = [row[:] for row in generic]
    if with_covers:
        splits = record["alternate_q80_cover_audit"]["splits"]
        labels = [split["label"] for split in splits]
        if len(labels) != len(set(labels)):
            raise ArithmeticError("duplicate retained split label")
        columns = integer_matrix([split["plus_point_coordinates_in_ordered_public_points"]
                                  for split in splits], len(splits), count)
        for index, row in enumerate(rows):
            row.extend(column[index] for column in columns)
    return rows


def verify_integer_relations(curve, public_points, target_points, rows):
    integer_matrix(rows, len(public_points), len(target_points))
    for column, target in enumerate(target_points):
        actual = curve(0)
        for row, point in zip(rows, public_points):
            actual += row[column] * point
        if actual != target:
            raise ArithmeticError(f"retained exact relation {column + 1} failed")


def verify_retained_independence(curve, points, record, helper):
    """Verify the original finite-quotient witnesses; do not select new primes."""
    from mod2_reduction_independence import (
        mod2_reduction_signature,
        short_curve_has_no_rational_2_torsion_modular_certificate,
    )

    primes = record["certificate_primes"]
    if (not primes or any(type(prime) is not int or not 3 <= prime <= 500 for prime in primes)
            or len(primes) != len(set(primes))
            or primes != [row["prime"] for row in record["signatures"]]):
        raise ArithmeticError("invalid retained independence prime roster")
    coefficients = [0, 0, 0, helper["python_fraction"](curve.a4()), helper["python_fraction"](curve.a6())]
    affine = [(helper["python_fraction"](point[0]), helper["python_fraction"](point[1])) for point in points]
    signatures = []
    for prime, expected in zip(primes, record["signatures"]):
        signature = mod2_reduction_signature(coefficients, affine, prime)
        actual = {"prime": signature.prime, "group_order": signature.group_order,
                  "doubled_subgroup_order": signature.doubled_subgroup_order,
                  "quotient_dimension": signature.quotient_dimension,
                  "rows": [list(row) for row in signature.rows]}
        if actual != expected:
            raise ArithmeticError(f"retained independence signature differs at prime {prime}")
        signatures.append(signature)
    rank = helper["combined_mod2_rank"](signatures, len(points))
    if rank != len(points) or record["combined_exact_rank_over_F2"] != rank:
        raise ArithmeticError("retained finite quotients do not prove displayed independence")
    prime = record["two_torsion_certificate_prime"]
    if (type(prime) is not int or not 3 <= prime <= 500
            or not short_curve_has_no_rational_2_torsion_modular_certificate(coefficients, prime)):
        raise ArithmeticError("retained rational two-torsion exclusion failed")
    return record


def replay_fibre(helper, config, hit, public, direct, covers, ring, record):
    """Reuse the pinned proof body; isolate discovery from witness checking."""
    with_covers = covers is not None
    if (record["curve_id"] != int(public["id"])
            or record["curve_id"] != int(hit["curve_id"])
            or record["native_chart"] != config["source_chart"]):
        raise ArithmeticError("retained fibre/chart identity changed")
    rows = coordinate_rows(record, with_covers=with_covers)
    original = helper["special_fibre"]
    namespace = dict(original.__globals__)

    def retained_coordinates(curve, public_points, target_points):
        verify_integer_relations(curve, public_points, target_points, rows)
        return namespace["matrix"](namespace["ZZ"], rows)

    def retained_split_order(inventory, parameter, polynomial_ring):
        splits, digest = helper["evaluate_cover_splits"](inventory, parameter, polynomial_ring)
        expected = [split["label"] for split in record["alternate_q80_cover_audit"]["splits"]]
        if [split["cover"]["label"] for split in splits] != expected:
            raise ArithmeticError("retained split labels/order differ from exhaustive evaluation")
        return splits, digest

    namespace["recover_coordinates"] = retained_coordinates
    namespace["evaluate_cover_splits"] = retained_split_order
    namespace["finite_reduction_certificate"] = lambda curve, points: verify_retained_independence(
        curve, points, record["public_point_independence"], helper)
    replay = FunctionType(original.__code__, namespace, original.__name__, original.__defaults__)
    result = replay(config, hit, public, direct, covers, ring)
    # The pinned generation record's height_recovery_separation_gate is retained
    # historical metadata, not a hypothesis or a fresh numerical measurement.
    if with_covers:
        actual, expected = result, record
    else:
        actual = {key: value for key, value in result.items() if key != "alternate_q80_cover_audit"}
        expected = {key: value for key, value in record.items() if key != "alternate_q80_cover_audit"}
        if result["alternate_q80_cover_audit"]["status"] != "NOT_RUN_NO_FROZEN_NATIVE_COVER_INVENTORY":
            raise ArithmeticError("transport-only replay unexpectedly evaluated covers")
    if actual != expected:
        changed = [key for key in set(actual) | set(expected) if actual.get(key) != expected.get(key)]
        raise ArithmeticError(f"curve {record['curve_id']} differs in {sorted(changed)}")
    return result


def selected_curves(requested):
    if requested is None:
        return CURVE_IDS
    if not requested or len(set(requested)) != len(requested) or not set(requested) <= set(CURVE_IDS):
        raise ValueError("select distinct curve ids from the seven-fibre certificate")
    return tuple(curve_id for curve_id in CURVE_IDS if curve_id in requested)


def run(*, with_covers, curve_ids, checkpoint=None):
    start = monotonic()
    if tuple(curve_ids) != selected_curves(list(curve_ids)):
        raise ValueError("curve selection must use certificate order")
    certificate = load_pinned(CERTIFICATE, CERTIFICATE_SHA256)
    if hashlib.sha256(PRODUCER.read_bytes()).hexdigest() != PRODUCER_SHA256:
        raise ArithmeticError("original producer changed; review compatibility before replay")
    helper = runpy.run_path(str(PRODUCER))
    if helper["SAGE_VERSION"] != certificate["software_assumptions"]["sage_version"]:
        raise ArithmeticError("Sage version differs from the pinned Smith-coordinate convention")
    if [row["curve_id"] for row in certificate["fibres"]] != list(CURVE_IDS):
        raise ArithmeticError("original fibre roster changed")
    records = {row["curve_id"]: row for row in certificate["fibres"]}
    inputs = {str(CERTIFICATE.relative_to(ROOT)): CERTIFICATE_SHA256,
              str(PRODUCER.relative_to(ROOT)): PRODUCER_SHA256}

    def read_input(path):
        key = str(path.relative_to(ROOT))
        result = load_pinned(path, certificate["inputs"][key])
        inputs[key] = certificate["inputs"][key]
        return result

    sweep, public, local, curve12 = (read_input(helper[key]) for key in ("SWEEP", "PUBLIC", "LOCAL", "CURVE12"))
    if local["status"] != "PASS_EXACT_LOCAL_FINGERPRINTS_FOR_ALL_69_RECOGNIZED_FIBRES":
        raise ArithmeticError("local fingerprint input changed")
    if (curve12["status"] != "PROVED_CURVE12_NATIVE_ALTERNATE_Q80_AND_DISPLAYED_QUOTIENT"
            or curve12["displayed_exceptional_quotient"]["free_basis_modulo_specialized_generic"]
            != helper["CURVE12_PREFERRED_BASIS"]):
        raise ArithmeticError("curve12 quotient prerequisite changed")
    hits = {int(row["curve_id"]): row for row in sweep["rational_j_hits_and_twists"]}
    publics = {int(row["id"]): row for row in public["records"]}
    ring = helper["PolynomialRing"](helper["QQ"], "u")
    completed = []
    report = {
        "schema": "elkies-k3.native-quotient-retained-replay.v1",
        "status": "INCOMPLETE_CHECKPOINT",
        "requested_curve_ids": list(curve_ids),
        "mode": "with-covers" if with_covers else "transports-only",
        "fibres": completed,
        "inputs": inputs,
        "replayer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "sage_version": helper["SAGE_VERSION"],
        "historical_metadata_not_remeasured": ["180-bit height recovery and its 2^-100 rounding gate"],
        "not_replayed": [
            "other calibration checkers: public projection, local features, norm8 incidence, highest-rank transports, dataset and holdout",
            "upstream proofs of the supplied saturated generic bases and complete cover inventories",
        ] + ([] if with_covers else ["all native cover inventories, their exhaustive split digests and visible spans"]),
        "boundary": "Shared-producer exact replay of selected displayed-subgroup components; no full Mordell-Weil group, rank upper bound, new discovery or assurance upgrade.",
    }
    if checkpoint:
        checkpoint(report)
    for config in helper["CHARTS"]:
        selected = [curve_id for curve_id in config["curve_ids"] if curve_id in curve_ids]
        if not selected:
            continue
        direct = read_input(config["direct"])
        if direct["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS":
            raise ArithmeticError("generic basis prerequisite changed")
        covers = read_input(config["covers"]) if with_covers else None
        if covers is not None and len(covers["bisections"]) != int(covers["construction"]["record_count"]):
            raise ArithmeticError("cover inventory incomplete")
        for curve_id in selected:
            print(f"R17RETAINED|curve={curve_id}|stage=exact_replay|heights=NOT_RUN", flush=True)
            before = monotonic()
            result = replay_fibre(helper, config, hits[curve_id], publics[curve_id],
                                  direct, covers, ring, records[curve_id])
            completed.append({"curve_id": curve_id, "displayed_quotient": result["displayed_exceptional_quotient"]["quotient"],
                              "generic_relations_verified": 17,
                              "split_relations_verified": len(result["alternate_q80_cover_audit"]["splits"]),
                              "covers_evaluated": result["alternate_q80_cover_audit"]["covers_evaluated"],
                              "elapsed_seconds": round(monotonic() - before, 3)})
            print(f"R17RETAINED|curve={curve_id}|status=PASS|elapsed={completed[-1]['elapsed_seconds']}", flush=True)
            report["elapsed_seconds"] = round(monotonic() - start, 3)
            if checkpoint:
                checkpoint(report)
    if [row["curve_id"] for row in completed] != list(curve_ids):
        raise ArithmeticError("requested fibre coverage incomplete")
    report["status"] = "PASS_SELECTED_NATIVE_COMPONENTS"
    report["elapsed_seconds"] = round(monotonic() - start, 3)
    if checkpoint:
        checkpoint(report)
    return report


def write_checkpoint(path, report):
    """Replace a reserved local receipt atomically after each completed fibre."""
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--transports-only", action="store_true")
    mode.add_argument("--with-covers", action="store_true")
    parser.add_argument("--curve", type=int, action="append", choices=CURVE_IDS)
    parser.add_argument("--receipt", type=Path, help="New local JSON checkpoint, updated after each completed fibre")
    args = parser.parse_args()
    curve_ids = selected_curves(args.curve)
    receipt = args.receipt.resolve() if args.receipt else None
    if receipt is not None:
        if not receipt.is_relative_to((ROOT / "artifacts/local").resolve()) or receipt.exists():
            parser.error("receipt must be a new file under research/artifacts/local/")
        receipt.parent.mkdir(parents=True, exist_ok=True)
        with receipt.open("x") as handle:
            handle.write(json.dumps({"status": "INCOMPLETE_CHECKPOINT", "requested_curve_ids": list(curve_ids), "fibres": []}) + "\n")
    result = run(with_covers=args.with_covers, curve_ids=curve_ids,
                 checkpoint=(lambda report: write_checkpoint(receipt, report)) if receipt else None)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(encoded, end="")


if __name__ == "__main__":
    main()
