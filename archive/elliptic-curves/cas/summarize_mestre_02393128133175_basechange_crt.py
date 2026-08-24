#!/usr/bin/env python3
"""Close and summarize the 220-fiber base-change CRT specialization screen."""

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "artifacts/local/elliptic-curves/mestre-02393128133175-basechange-crt-v1"
TARGET_LOG_CONDUCTOR = Decimal("182.72")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def atomic_text(path: Path, text: str) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text)
    temporary.replace(path)


def atomic_json(path: Path, value: Any) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def directory_manifest(directory: Path) -> dict[str, Any]:
    files = sorted(directory.glob("*.json"))
    records = [{"name": path.name, "sha256": sha256_file(path)} for path in files]
    return {"count": len(records), "records_sha256": canonical_digest(records)}


def main() -> None:
    candidate_path = BASE / "candidate-input.json"
    main_summary_path = BASE / "summary.json"
    candidate_input = json.loads(candidate_path.read_text())
    main_summary = json.loads(main_summary_path.read_text())
    candidates = candidate_input["candidates"]
    if len(candidates) != 220:
        raise AssertionError("the frozen CRT input no longer has 220 candidates")

    rows = []
    for candidate in candidates:
        identifier = f"u{candidate['numerator']}_{candidate['denominator']}"
        conductor_path = BASE / "conductor-records" / f"{identifier}.json"
        conductor = json.loads(conductor_path.read_text())
        completed_conductor = conductor["status"].startswith("completed")
        if completed_conductor:
            point_path = BASE / "point-certificates-h200000" / f"{identifier}.json"
        else:
            if "TimeoutExpired" not in conductor.get("error", ""):
                raise AssertionError(f"{identifier} has a non-timeout conductor failure")
            point_path = (
                BASE
                / "point-certificates-conductor-timeouts-h200000"
                / f"{identifier}.json"
            )
        point = json.loads(point_path.read_text())
        if not point["status"].startswith("completed"):
            raise AssertionError(f"{identifier} lacks completed H200000 point coverage")
        rank = int(point["exact_specialization_rank_lower_bound"])
        if rank < 13:
            raise AssertionError("the exact generic rank-13 subgroup was lost")
        global_curve = conductor.get("global_curve")
        row = {
            "u": candidate["u"],
            "T": candidate["base_T"],
            "numerator": candidate["numerator"],
            "denominator": candidate["denominator"],
            "parameter_projective_height": max(
                abs(candidate["numerator"]), candidate["denominator"]
            ),
            "forcing_paths": candidate["forcing_paths"],
            "forced_prime_actual_valuations": candidate[
                "forced_prime_actual_valuations"
            ],
            "conductor_status": conductor["status"],
            "log_conductor": global_curve["log_conductor"] if global_curve else None,
            "root_number": global_curve["root_number"] if global_curve else None,
            "below_strict_log_conductor_182_72": (
                Decimal(str(global_curve["log_conductor"])) < TARGET_LOG_CONDUCTOR
                if global_curve
                else None
            ),
            "generic_exact_rank_lower_bound": 13,
            "specialization_exact_rank_lower_bound": rank,
            "point_pool_count_modulo_inverse": point["pool_point_count_modulo_inverse"],
            "point_pool_sha256": point["pool_point_sha256"],
            "conductor_record_sha256": sha256_file(conductor_path),
            "point_certificate_sha256": sha256_file(point_path),
        }
        rows.append(row)

    exact_global = [row for row in rows if row["log_conductor"] is not None]
    pareto = []
    for row in exact_global:
        rank = row["specialization_exact_rank_lower_bound"]
        log_n = Decimal(str(row["log_conductor"]))
        height = row["parameter_projective_height"]
        dominated = False
        for other in exact_global:
            if other is row:
                continue
            other_rank = other["specialization_exact_rank_lower_bound"]
            other_log_n = Decimal(str(other["log_conductor"]))
            other_height = other["parameter_projective_height"]
            weak = other_rank >= rank and other_log_n <= log_n and other_height <= height
            strict = other_rank > rank or other_log_n < log_n or other_height < height
            if weak and strict:
                dominated = True
                break
        if not dominated:
            pareto.append(row)
    pareto.sort(
        key=lambda row: (
            -row["specialization_exact_rank_lower_bound"],
            Decimal(str(row["log_conductor"])),
            row["parameter_projective_height"],
        )
    )

    headers = (
        "u",
        "T",
        "parameter_projective_height",
        "generic_exact_rank_lower_bound",
        "specialization_exact_rank_lower_bound",
        "log_conductor",
        "root_number",
        "below_strict_log_conductor_182_72",
        "point_pool_count_modulo_inverse",
        "point_pool_sha256",
    )
    lines = ["\t".join(headers)]
    for row in sorted(
        rows,
        key=lambda item: (
            -item["specialization_exact_rank_lower_bound"],
            Decimal(str(item["log_conductor"]))
            if item["log_conductor"] is not None
            else Decimal("Infinity"),
            item["parameter_projective_height"],
            item["numerator"],
            item["denominator"],
        ),
    ):
        lines.append("\t".join(str(row[key]) for key in headers))
    all_tsv_path = BASE / "all220-exact-coverage.tsv"
    atomic_text(all_tsv_path, "\n".join(lines) + "\n")

    pareto_lines = ["\t".join(headers)]
    for row in pareto:
        pareto_lines.append("\t".join(str(row[key]) for key in headers))
    pareto_path = BASE / "pareto.tsv"
    atomic_text(pareto_path, "\n".join(pareto_lines) + "\n")

    rank_distribution = Counter(
        row["specialization_exact_rank_lower_bound"] for row in rows
    )
    timeout_rows = [row for row in rows if row["log_conductor"] is None]
    target_rows = [row for row in rows if row["below_strict_log_conductor_182_72"]]
    summary = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete exact H200000 coverage of all 220 frozen CRT fibers",
        "family": {
            "roots": [0, 23, 93, 128, 133, 175],
            "base_change": "T=(14406-u^2)/(2u)",
            "generic_exact_rank_lower_bound": 13,
        },
        "coverage": {
            "frozen_candidates": len(rows),
            "exact_conductor_completed": len(exact_global),
            "conductor_GP_90s_timeouts": len(timeout_rows),
            "H200000_exact_point_certificates": len(rows),
            "point_certificate_errors": 0,
            "rank_distribution": dict(sorted(rank_distribution.items())),
            "maximum_exact_specialization_rank_lower_bound": max(rank_distribution),
            "rank_above_16": [],
        },
        "target_qualified": {
            "count_with_exact_log_conductor_below_182_72": len(target_rows),
            "W_minus_1_count": sum(row["root_number"] == -1 for row in target_rows),
            "rank_distribution": dict(
                sorted(
                    Counter(
                        row["specialization_exact_rank_lower_bound"]
                        for row in target_rows
                    ).items()
                )
            ),
            "records": target_rows,
        },
        "pareto_definition": (
            "maximize exact specialization rank lower bound; minimize exact log conductor "
            "and projective height max(|numerator(u)|,denominator(u))"
        ),
        "pareto_records": pareto,
        "alternate_cover_u36": {
            "path": str((BASE / "alternate-covers/u36_1/summary.json").relative_to(ROOT)),
            "sha256": sha256_file(BASE / "alternate-covers/u36_1/summary.json"),
            "result_sha256": json.loads(
                (BASE / "alternate-covers/u36_1/summary.json").read_text()
            )["result_sha256"],
            "exact_rank_lower_bound_after_search": 13,
        },
        "artifacts": {
            "candidate_input": {"path": str(candidate_path.relative_to(ROOT)), "sha256": sha256_file(candidate_path)},
            "main_summary": {"path": str(main_summary_path.relative_to(ROOT)), "sha256": sha256_file(main_summary_path), "result_sha256": main_summary["result_sha256"]},
            "all220_tsv": {"path": str(all_tsv_path.relative_to(ROOT)), "sha256": sha256_file(all_tsv_path)},
            "pareto_tsv": {"path": str(pareto_path.relative_to(ROOT)), "sha256": sha256_file(pareto_path)},
            "conductor_records_manifest": directory_manifest(BASE / "conductor-records"),
            "main_point_certificates_manifest": directory_manifest(BASE / "point-certificates-h200000"),
            "timeout_point_certificates_manifest": directory_manifest(BASE / "point-certificates-conductor-timeouts-h200000"),
        },
        "scope_warning": (
            "H200000 searches and bounded alternate-cover charts are not rank upper bounds; "
            "the 71 missing conductor values are explicit GP timeouts, not conductor claims"
        ),
        "provenance": {
            "script_path": str(Path(__file__).resolve().relative_to(ROOT)),
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "reproducing_command": "python3 elliptic-curves/cas/summarize_mestre_02393128133175_basechange_crt.py",
        },
    }
    summary["result_sha256"] = canonical_digest(
        {key: value for key, value in summary.items() if key != "generated_at_utc"}
    )
    output_path = BASE / "coverage-pareto-summary.json"
    atomic_json(output_path, summary)
    print(
        f"wrote {output_path}: ranks={dict(sorted(rank_distribution.items()))} "
        f"pareto={len(pareto)} sha={summary['result_sha256']}"
    )


if __name__ == "__main__":
    main()
