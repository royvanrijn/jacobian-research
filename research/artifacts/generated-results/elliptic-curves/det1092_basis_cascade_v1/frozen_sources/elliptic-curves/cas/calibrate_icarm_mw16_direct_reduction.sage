#!/usr/bin/env sage-python
"""Calibrate the direct-reduction quartic backend on MW16 jump controls.

This launcher replays the existing complement-blind curve-398 and curve-400
adaptive searches while replacing only their pointed-quartic search backend.
Raw checkpoints stay under ``artifacts/local``.  A compact generated summary
is emitted only when curve 398 again recovers fourteen quotient directions,
curve 400 again recovers twelve, and every declared chart completes.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
DIRECT = ROOT / "elliptic-curves/cas/half_lattice_direct_reduction.py"
CURVE398 = ROOT / "elliptic-curves/cas/run_curve398_mw16_adaptive_half_lattice_search.sage"
CURVE400 = ROOT / "elliptic-curves/cas/run_icarm_mw16_curve400_adaptive_calibration.sage"
RAW_DIR = ROOT / "artifacts/local/elliptic-curves/mw16-direct-reduction-controls"
RAW398 = RAW_DIR / "curve398-direct-reduction-v1.json"
RAW400 = RAW_DIR / "curve400-direct-reduction-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_direct_reduction_calibration_v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def inject_backend(script_module, direct_module) -> None:
    """Replace only the legacy module's engine during a control replay."""

    real_loader = script_module.SourceFileLoader
    legacy_path = script_module.LEGACY.resolve()

    class Proxy:
        def __init__(self, loader):
            self.loader = loader

        def load_module(self):
            module = self.loader.load_module()
            module.engine = direct_module
            return module

    def patched_loader(name, path):
        loader = real_loader(name, str(path))
        return Proxy(loader) if Path(path).resolve() == legacy_path else loader

    script_module.SourceFileLoader = patched_loader


def run_control(name: str, script: Path, output: Path, direct_module, args) -> None:
    module = SourceFileLoader(f"mw16_direct_{name}", str(script)).load_module()
    inject_backend(module, direct_module)
    control_argv = [
        str(script),
        "--output",
        str(output),
        "--height-bound",
        str(args.height_bound),
        "--timeout-seconds",
        str(args.timeout_seconds),
        "--stack-bytes",
        str(args.stack_bytes),
        "--relation-chunk-size",
        str(args.relation_chunk_size),
        "--relation-timeout-seconds",
        str(args.relation_timeout_seconds),
    ]
    previous = sys.argv
    try:
        sys.argv = control_argv
        module.main()
    finally:
        sys.argv = previous
    payload = json.loads(output.read_text())
    payload["direct_reduction_experiment"] = {
        **direct_module.provenance(),
        "launcher": relative(Path(__file__)),
        "launcher_sha256": digest(Path(__file__)),
        "control_script": relative(script),
        "control_script_sha256": digest(script),
        "actual_backend_overrides_legacy_engine": True,
    }
    payload["reproducing_command"] = (
        "sage -python elliptic-curves/cas/"
        f"calibrate_icarm_mw16_direct_reduction.sage --control {name}"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def nested_search_records(value):
    if isinstance(value, dict):
        if isinstance(value.get("search"), dict) and "status" in value["search"]:
            yield value["search"]
        for item in value.values():
            yield from nested_search_records(item)
    elif isinstance(value, list):
        for item in value:
            yield from nested_search_records(item)


def summarize(args, direct_module) -> dict:
    if not RAW398.is_file() or not RAW400.is_file():
        raise ArithmeticError("both direct-reduction control checkpoints are required")
    curve398 = json.loads(RAW398.read_text())
    curve400 = json.loads(RAW400.read_text())
    if curve398.get("status") != "STOPPED_AT_DECLARED_LIFT_LIMIT":
        raise ArithmeticError("curve-398 direct-reduction replay did not reach its frozen limit")
    if len(curve398.get("current_basis", [])) != 30:
        raise ArithmeticError("curve-398 direct-reduction replay did not recover M30")
    if curve400.get("status") != "PASS_COMPLETE_CURVE400_ADAPTIVE_CALIBRATION":
        raise ArithmeticError("curve-400 direct-reduction replay did not complete")
    if int(curve400.get("exact_quotient_rank_recovered_total", -1)) != 12:
        raise ArithmeticError("curve-400 direct-reduction replay did not recover twelve directions")

    rows398 = list(nested_search_records(curve398))
    rows400 = list(nested_search_records(curve400))
    expected_counts = {"398": 384, "400": 128}
    if len(rows398) != expected_counts["398"] or len(rows400) != expected_counts["400"]:
        raise ArithmeticError("direct-reduction control chart census changed")
    for curve_id, rows in (("398", rows398), ("400", rows400)):
        if any(row.get("backend") != direct_module.BACKEND_NAME for row in rows):
            raise ArithmeticError(f"curve {curve_id} contains a non-direct search record")
        if any(row.get("status") != "bounded_search_complete" for row in rows):
            raise ArithmeticError(f"curve {curve_id} contains an incomplete chart")
        if any(row.get("hyperellminimalmodel_called") is not False for row in rows):
            raise ArithmeticError(f"curve {curve_id} called quartic minimalization")

    return {
        "schema": "elliptic-curves.icarm-mw16-direct-reduction-calibration.v1",
        "status": "PASS_EXACT_DIRECT_REDUCTION_MW16_CALIBRATION",
        "backend": direct_module.provenance(),
        "declared_budget": {
            "height_bound_each_quartic": args.height_bound,
            "timeout_seconds_each_quartic": args.timeout_seconds,
            "stack_bytes_each_quartic": args.stack_bytes,
            "unrestricted_point_search": False,
        },
        "controls": {
            "curve398": {
                "chart_count": len(rows398),
                "completed_chart_count": len(rows398),
                "generic_rank": 16,
                "exact_quotient_rank_recovered": 14,
                "basis_rank_after": 30,
                "raw_checkpoint": relative(RAW398),
                "raw_checkpoint_sha256": digest(RAW398),
            },
            "curve400": {
                "chart_count": len(rows400),
                "completed_chart_count": len(rows400),
                "generic_rank": 16,
                "exact_quotient_rank_recovered": 12,
                "basis_rank_after": 28,
                "raw_checkpoint": relative(RAW400),
                "raw_checkpoint_sha256": digest(RAW400),
            },
        },
        "inputs": {
            relative(path): digest(path)
            for path in (Path(__file__), DIRECT, CURVE398, CURVE400)
        },
        "claim_boundary": [
            "Both searches retain their original complement-blind input boundaries.",
            "All returned points are mapped back and checked by exact rational arithmetic.",
            "All 512 declared quartic charts completed; no timeout is counted as a miss.",
            "The replay certifies only the recovered subgroups M30 and M28.",
            "It makes no rank upper bound, global minimality, or Selmer claim.",
        ],
        "reproducing_command": (
            "sage -python elliptic-curves/cas/"
            "calibrate_icarm_mw16_direct_reduction.sage --control all"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--control", choices=("all", "curve398", "curve400", "summarize"), default="all"
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--height-bound", type=int, default=100_000)
    parser.add_argument("--timeout-seconds", type=float, default=15.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--relation-chunk-size", type=int, default=64)
    parser.add_argument("--relation-timeout-seconds", type=float, default=180.0)
    args = parser.parse_args()
    if args.height_bound <= 0 or not 0 < args.timeout_seconds <= 60:
        raise SystemExit("invalid quartic-search budget")

    direct = SourceFileLoader("mw16_direct_calibration_backend", str(DIRECT)).load_module()
    if args.control in ("all", "curve398"):
        run_control("curve398", CURVE398, RAW398, direct, args)
    if args.control in ("all", "curve400"):
        run_control("curve400", CURVE400, RAW400, direct, args)
    if args.control in ("all", "summarize"):
        payload = summarize(args, direct)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(
            "MW16DIRECTCAL|curve398_gain=14|curve400_gain=12|charts=512|"
            f"output={relative(args.output)}|status={payload['status']}",
            flush=True,
        )


if __name__ == "__main__":
    main()
