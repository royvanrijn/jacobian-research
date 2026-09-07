#!/usr/bin/env python3
# <!-- status-consumer: EC-K3-ELKIES-2026-RESIDUAL-SELMER-GATE f7a8c94736f1b44f -->
"""Run a strictly bounded PARI 2-cover search on one q12o5867 fibre.

The input is an exact q12o5867 specialization artifact containing the global
minimal model and the independently certified ordered rank-17 subgroup.  A
fresh PARI/GP process first calls ``ell2cover``.  If that succeeds, every
returned quartic is searched with ``hyperellratpoints`` and every image is
mapped back by PARI's exact covering map.  The Python supervisor checks those
images exactly and measures escape from the rank-17 subgroup with the
repository's finite-reduction quotient implementation.

This is a bounded point-search backend.  It is disabled unless an
unconditional completed residual 2-Selmer artifact for the same minimal curve
has already passed the rank-32 room gate.  A timeout, resource stop, empty
cover search, or zero observed quotient gain is not a rank upper bound.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from typing import Any, Sequence


REPOSITORY = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = REPOSITORY / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
SCRIPTS = ELLIPTIC_ROOT / "scripts"
sys.path.insert(0, str(ELLIPTIC_ROOT))
sys.path.insert(0, str(CAS))
sys.path.insert(0, str(SCRIPTS))
from research_runtime.supervisor import Limits, capture_record

from elliptic_candidate_record import is_on_weierstrass_curve  # noqa: E402
from elkies_residual_selmer_gate import require_gate_for_specialization  # noqa: E402
from probe_q12o5867_mwrank import (  # noqa: E402
    exact_escape_records,
    parse_point,
    sign_key,
)


Q = Fraction
PROTOCOL = "Q12P2"
INPUT_STATUS = "PASS_EXACT_Q12O5867_SPECIALIZED_GENERIC_RANK17_LOWER_BOUND"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def gp_program(
    model: Sequence[int], *, stack_bytes: int, search_height: int
) -> str:
    """Return the exact GP program used by the owned subprocess."""

    if len(model) != 5:
        raise ValueError("a generalized Weierstrass model needs five coefficients")
    if stack_bytes < 1 or search_height < 1:
        raise ValueError("stack bytes and search height must be positive")
    coefficients = ",".join(str(int(value)) for value in model)
    return f'''default(parisizemax,{stack_bytes});
E=ellinit([{coefficients}]);
print("{PROTOCOL}|stage=input|status=complete|pari_version=",version());
fail(stage,err)={{
  print("{PROTOCOL}|stage=",stage,"|status=error|pari_error=",err);
  quit(2);
}};
print("{PROTOCOL}|stage=ell2cover|status=start");
gettime();
iferr(C=ell2cover(E),ERR,fail("ell2cover",ERR));
print("{PROTOCOL}|stage=ell2cover|status=complete|milliseconds=",gettime(),"|cover_count=",#C);
emitpoint(i,M,Q)={{
  my(xx=Q[1],yy=Q[2],ex,ey);
  ex=subst(subst(M[1],x,xx),y,yy);
  ey=subst(subst(M[2],x,xx),y,yy);
  if(ellisoncurve(E,[ex,ey]),print("{PROTOCOL}|stage=candidate|cover_index=",i,"|cover_x=",xx,"|cover_y=",yy,"|curve_x=",ex,"|curve_y=",ey));
}};
searchcover(i)={{
  my(R=C[i][1],M=C[i][2],H);
  print("{PROTOCOL}|stage=cover|status=start|cover_index=",i);
  print("{PROTOCOL}|cover_index=",i,"|quartic=",R);
  print("{PROTOCOL}|cover_index=",i,"|map_x=",M[1]);
  print("{PROTOCOL}|cover_index=",i,"|map_y=",M[2]);
  gettime();
  iferr(H=hyperellratpoints(R,{search_height}),ERR,fail("cover_search",ERR));
  print("{PROTOCOL}|stage=cover_search|status=complete|cover_index=",i,"|milliseconds=",gettime(),"|raw_point_count=",#H);
  for(j=1,#H,emitpoint(i,M,H[j]));
  print("{PROTOCOL}|stage=cover|status=complete|cover_index=",i);
}};
for(i=1,#C,searchcover(i));
print("{PROTOCOL}|stage=all|status=complete");
quit
'''


def _marker_fields(line: str) -> dict[str, str] | None:
    if not line.startswith(f"{PROTOCOL}|"):
        return None
    fields: dict[str, str] = {}
    for item in line.rstrip("\n").split("|")[1:]:
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        fields[key] = value
    return fields


def parse_gp_output(stdout: str) -> dict[str, Any]:
    """Parse protocol markers without interpreting ordinary PARI output."""

    cover_count: int | None = None
    completed_cover_indices: set[int] = set()
    covers: dict[int, dict[str, Any]] = {}
    candidates: list[dict[str, Any]] = []
    pari_version: str | None = None
    ell2cover_completed = False
    all_completed = False
    last_stage: str | None = None
    for line in stdout.splitlines():
        fields = _marker_fields(line)
        if fields is None:
            continue
        last_stage = fields.get("stage", last_stage)
        pari_version = fields.get("pari_version", pari_version)
        if fields.get("stage") == "ell2cover" and fields.get("status") == "complete":
            ell2cover_completed = True
            cover_count = int(fields["cover_count"])
        if fields.get("stage") == "cover" and fields.get("status") == "complete":
            completed_cover_indices.add(int(fields["cover_index"]))
        if "cover_index" in fields and fields.get("stage") is None:
            cover_index = int(fields["cover_index"])
            cover = covers.setdefault(cover_index, {"cover_index": cover_index})
            for key in ("quartic", "map_x", "map_y"):
                if key in fields:
                    cover[key] = fields[key]
        if fields.get("stage") == "cover_search" and fields.get("status") == "complete":
            cover_index = int(fields["cover_index"])
            cover = covers.setdefault(cover_index, {"cover_index": cover_index})
            cover["search_milliseconds"] = int(fields["milliseconds"])
            cover["raw_point_count"] = int(fields["raw_point_count"])
        if fields.get("stage") == "candidate":
            candidates.append(
                {
                    "cover_index": int(fields["cover_index"]),
                    "cover_point": [fields["cover_x"], fields["cover_y"]],
                    "curve_point": [fields["curve_x"], fields["curve_y"]],
                }
            )
        if fields.get("stage") == "all" and fields.get("status") == "complete":
            all_completed = True
    return {
        "pari_version": pari_version,
        "ell2cover_completed": ell2cover_completed,
        "cover_count": cover_count,
        "completed_cover_indices": sorted(completed_cover_indices),
        "covers": [covers[index] for index in sorted(covers)],
        "raw_candidate_images": candidates,
        "all_completed": all_completed,
        "last_stage": last_stage,
    }






def load_specialization(path: Path) -> tuple[dict[str, Any], tuple[int, ...]]:
    artifact = json.loads(path.read_text())
    if artifact.get("status") != INPUT_STATUS:
        raise ValueError("input is not an exact q12o5867 rank-17 specialization")
    model_values = tuple(Q(value) for value in artifact["global_minimal_specialization"]["model"])
    if len(model_values) != 5 or any(value.denominator != 1 for value in model_values):
        raise ValueError("the global minimal specialization is not an integral model")
    model = tuple(value.numerator for value in model_values)
    points = tuple(
        parse_point(record)
        for record in artifact["global_minimal_specialization"]["points"]
    )
    if len(points) != 17 or len(set(points)) != 17:
        raise ValueError("the specialization does not contain 17 distinct baseline points")
    if any(not is_on_weierstrass_curve(model, point) for point in points):
        raise ArithmeticError("a baseline point misses the global minimal model")
    return artifact, model


def exact_candidate_records(
    artifact: dict[str, Any],
    model: Sequence[int],
    parsed: dict[str, Any],
    relation_primes: Sequence[int],
    reduction_prime_bound: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    baseline = tuple(
        parse_point(record)
        for record in artifact["global_minimal_specialization"]["points"]
    )
    seen = {sign_key(model, point) for point in baseline}
    exact_records: list[dict[str, Any]] = []
    novel_points = []
    for raw in parsed["raw_candidate_images"]:
        point = parse_point(raw["curve_point"])
        if not is_on_weierstrass_curve(model, point):
            raise ArithmeticError("a parsed PARI cover image misses the minimal model")
        key = sign_key(model, point)
        novel = key not in seen
        if novel:
            seen.add(key)
            novel_points.append(point)
        exact_records.append({**raw, "exact_on_curve": True, "novel_modulo_sign": novel})
    escape = exact_escape_records(
        artifact, novel_points, relation_primes, reduction_prime_bound
    )
    return exact_records, escape


def run(args: argparse.Namespace) -> int:
    artifact, model = load_specialization(args.input)
    descent_gate = require_gate_for_specialization(
        args.residual_selmer_gate,
        artifact,
        requested_search_limits={
            "search_height": args.search_height,
            "wall_seconds": args.timeout,
            "stack_bytes": args.stack_bytes,
            "rss_limit_bytes": args.rss_limit_bytes,
        },
    )
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    program = gp_program(
        model, stack_bytes=args.stack_bytes, search_height=args.search_height
    )
    supervision=capture_record([executable,'-f','-q','-s',str(args.stack_bytes)],input_text=program,
        limits=Limits(args.timeout,args.rss_limit_bytes,pari_stack_bytes=args.stack_bytes))
    stdout,stderr=supervision['stdout'],supervision['stderr']
    peak_rss=supervision['peak_observed_rss_bytes'];wall_seconds=supervision['wall_seconds']
    outcome='running' if supervision['outcome']=='completed' else supervision['outcome']
    parsed = parse_gp_output(stdout)
    fatal = any("***" in line and "Warning:" not in line for line in stderr.splitlines())
    if outcome == "running":
        outcome = (
            "completed"
            if supervision["returncode"] == 0 and parsed["all_completed"] and not fatal
            else "pari_failure"
        )
    exact_records, escape = exact_candidate_records(
        artifact,
        model,
        parsed,
        args.relation_primes,
        args.reduction_prime_bound,
    )
    status = (
        "PASS_BOUNDED_PARI_TWO_COVER_SEARCH"
        if outcome == "completed"
        else "BOUNDED_PARI_TWO_COVER_PROBE_NO_COMPLETE_RESULT"
    )
    result = {
        "schema": "elliptic-curves.q12o5867-pari-two-cover-probe.v1",
        "status": status,
        "mathematical_status": "bounded_exact_point_search_not_rank_upper_bound",
        "input": {
            "specialization_artifact": str(args.input.resolve()),
            "sha256": sha256_file(args.input),
            "parameter": artifact["parameter"],
            "global_minimal_model": [str(value) for value in model],
            "certified_baseline_rank_lower_bound": 17,
        },
        "backend": {
            "engine": "PARI/GP ell2cover + hyperellratpoints",
            "executable": executable,
            "pari_version": parsed["pari_version"],
            "python_version": platform.python_version(),
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "gp_program_sha256": hashlib.sha256(program.encode()).hexdigest(),
        },
        "bounds": {
            "cover_naive_height": args.search_height,
            "wall_timeout_seconds": args.timeout,
            "pari_stack_bytes": args.stack_bytes,
            "rss_limit_bytes": args.rss_limit_bytes,
            "relation_primes": list(args.relation_primes),
            "reduction_prime_bound": args.reduction_prime_bound,
        },
        "process": {
            "outcome": outcome,
            "returncode": supervision["returncode"],
            "wall_seconds": wall_seconds,
            "peak_observed_rss_bytes": peak_rss,
            "last_stage": parsed["last_stage"],
            "raw_stdout": stdout,
            "raw_stderr": stderr,
        },
        "cover_search": {
            "ell2cover_completed": parsed["ell2cover_completed"],
            "cover_count": parsed["cover_count"],
            "completed_cover_indices": parsed["completed_cover_indices"],
            "covers": parsed["covers"],
            "exact_candidate_images": exact_records,
            "novel_image_count_modulo_sign": sum(
                bool(record["novel_modulo_sign"]) for record in exact_records
            ),
        },
        "finite_quotient_escape": escape,
        "residual_selmer_gate": {
            "path": str(args.residual_selmer_gate.resolve()),
            "sha256": sha256_file(args.residual_selmer_gate),
            "status": descent_gate["status"],
        },
        "promotion_threshold": 15,
        "promoted": escape["maximum_marginal_dimension"] >= 15,
        "claim_boundary": [
            "Completion constructs PARI's everywhere locally soluble 2-cover basis and performs a bounded rational-point search on each quartic.",
            "A timeout, resource stop, empty search, or zero observed quotient gain is not a rank upper bound.",
            "A residual 2-cover without a rational point may represent a Tate-Shafarevich class.",
            "No specialization is promoted unless exact finite-quotient gain reaches 15.",
        ],
        "reproducing_command": (
            f".venv/bin/python {Path(__file__).relative_to(REPOSITORY)} "
            f"--input {args.input} --output {args.output} "
            f"--residual-selmer-gate {args.residual_selmer_gate} "
            f"--search-height {args.search_height} --timeout {args.timeout} "
            f"--stack-bytes {args.stack_bytes} --rss-limit-bytes {args.rss_limit_bytes}"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.overwrite else "x"
    with args.output.open(mode) as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"status={status}")
    print(f"outcome={outcome}")
    print(f"last_stage={parsed['last_stage']}")
    print(f"ell2cover_completed={parsed['ell2cover_completed']}")
    print(f"cover_count={parsed['cover_count']}")
    print(f"novel_image_count={result['cover_search']['novel_image_count_modulo_sign']}")
    print(f"maximum_marginal_dimension={escape['maximum_marginal_dimension']}")
    print(f"output={args.output.resolve()}")
    return 0 if outcome in {"completed", "strict_wall_timeout", "strict_rss_limit"} else 1


def parse_relation_primes(text: str) -> tuple[int, ...]:
    try:
        result = tuple(int(value) for value in text.split(",") if value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("relation primes must be integers") from error
    if not result or any(value < 2 for value in result):
        raise argparse.ArgumentTypeError("relation primes must be >=2")
    return result


def main() -> None:
    sys.set_int_max_str_digits(0)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--residual-selmer-gate",
        type=Path,
        required=True,
        help="passing unconditional gate for this exact minimal curve",
    )
    parser.add_argument("--search-height", type=int, default=1_000_000)
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--stack-bytes", type=int, default=2_000_000_000)
    parser.add_argument("--rss-limit-bytes", type=int, default=4_000_000_000)
    parser.add_argument("--relation-primes", type=parse_relation_primes, default=(2, 3, 5))
    parser.add_argument("--reduction-prime-bound", type=int, default=500)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.search_height < 1:
        parser.error("--search-height must be positive")
    if args.timeout <= 0 or args.timeout > 86400:
        parser.error("--timeout must lie in (0,86400]")
    if min(args.stack_bytes, args.rss_limit_bytes) < 64_000_000:
        parser.error("stack and RSS bounds must each be at least 64MB")
    if args.reduction_prime_bound < 3:
        parser.error("--reduction-prime-bound must be at least 3")
    raise SystemExit(run(args))


if __name__ == "__main__":
    main()
