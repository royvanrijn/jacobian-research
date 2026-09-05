#!/usr/bin/env python3
"""Run Simon's exact 2-Selmer local conditions in the GP owning a BNF.

PARI binary checkpoints are safest when reloaded by the same GP build that
created them.  This supervisor therefore runs the Simon routines directly in
that GP process, rather than importing a PARI 2.19 checkpoint into Sage's
older libpari.  It records the complete Selmer dimension, algebraic basis
representatives, and leave-one-place-out local-condition ranks.  Known-point
alignment remains a separate exact step unless the total dimension equals the
independently certified rank lower bound.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import time


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
import sys

sys.path.insert(0, str(CAS))
from run_elkies_2026_relative_2selmer_checkpointed import (  # noqa: E402
    SIMON_GP_FUNCTION,
)


from research_runtime.supervisor import Limits, capture, run, preserve_previous
from research_runtime.store import checkpoint as write_checkpoint


SCHEMA = "elliptic-curves.elkies-2026-pari219-selmer-from-bnf.v1"
PROTOCOL = "ELKIESR17PARI219SELMER"
SUMMARY_RE = re.compile(
    rf"^{PROTOCOL}\|stage=summary\|selmer=(?P<selmer>\d+)"
    r"\|sclass=(?P<sclass>\d+)\|norm=(?P<norm>\d+)"
    r"\|local_rank=(?P<local_rank>\d+)\|bad=(?P<bad>.*)$"
)
DELETE_RE = re.compile(
    rf"^{PROTOCOL}\|stage=delete_one\|place=(?P<place>-?\d+)"
    r"\|allowed=(?P<allowed>\d+)\|alone=(?P<alone>\d+)"
    r"\|omitted=(?P<omitted>\d+)\|rank=(?P<rank>\d+)$"
)
BASIS_RE = re.compile(
    rf"^{PROTOCOL}\|stage=basis\|index=(?P<index>\d+)\|alpha=(?P<alpha>.*)$"
)
NORM_BASIS_RE = re.compile(
    rf"^{PROTOCOL}\|stage=norm_basis\|basis=(?P<basis>.*)$"
)
LOCAL_BASIS_RE = re.compile(
    rf"^{PROTOCOL}\|stage=local_basis\|place=(?P<place>-?\d+)"
    r"\|basis=(?P<basis>.*)$"
)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def gp_quote(path: Path) -> str:
    return str(path).replace("\\", "\\\\").replace('"', '\\"')


def parse_model(value: str) -> list[str]:
    try:
        model = json.loads(value)
    except json.JSONDecodeError as error:
        raise argparse.ArgumentTypeError("model must be a JSON array") from error
    if not isinstance(model, list) or len(model) != 5:
        raise argparse.ArgumentTypeError("model must contain five coefficients")
    return [str(coefficient) for coefficient in model]


def pari219_compatible_simon_source(path: Path) -> str:
    """Apply the documented PARI 2.19 interval syntax to pinned Simon code."""

    source = path.read_text()
    replacements = {
        "polsturm(pol,-a,a)": "polsturm(pol,[-a,a])",
        "polsturm(pol,a,c)": "polsturm(pol,[a,c])",
        "polsturm(elt2,rootapprox[1],rootapprox[2])": (
            "polsturm(elt2,[rootapprox[1],rootapprox[2]])"
        ),
    }
    for old, new in replacements.items():
        source = source.replace(old, new)
    return source


def gp_version(gp: Path) -> str:
    result = capture([str(gp), "-q", "-f"], input_text="print(version())\n",
                     limits=Limits(10, 256_000_000))
    return result.stdout.strip()


def parse_log(
    log: str,
) -> tuple[
    dict[str, object] | None,
    list[dict[str, object]],
    list[dict[str, object]],
    str | None,
    list[dict[str, object]],
]:
    summary = None
    deleted = []
    basis = []
    norm_basis = None
    local_bases = []
    for line in log.splitlines():
        match = SUMMARY_RE.match(line)
        if match:
            bad = match.group("bad")
            summary = {
                "two_selmer_dimension": int(match.group("selmer")),
                "global_s_squareclass_dimension": int(match.group("sclass")),
                "global_norm_square_subspace_dimension": int(match.group("norm")),
                "full_local_condition_matrix_rank_on_norm_subspace": int(
                    match.group("local_rank")
                ),
                "bad_rational_primes_gp": bad,
            }
            continue
        match = DELETE_RE.match(line)
        if match:
            place = int(match.group("place"))
            deleted.append(
                {
                    "place": "infinity" if place == -1 else str(place),
                    "allowed_subspace_dimension_in_global_s_squareclasses": int(
                        match.group("allowed")
                    ),
                    "norm_subspace_intersection_dimension_for_this_place_alone": int(
                        match.group("alone")
                    ),
                    "selmer_candidate_dimension_after_deleting_this_place": int(
                        match.group("omitted")
                    ),
                    "matrix_rank_after_deleting_this_place": int(match.group("rank")),
                }
            )
            continue
        match = BASIS_RE.match(line)
        if match:
            basis.append(
                {
                    "index": int(match.group("index")),
                    "field_squareclass_representative": match.group("alpha"),
                }
            )
            continue
        match = NORM_BASIS_RE.match(line)
        if match:
            norm_basis = match.group("basis")
            continue
        match = LOCAL_BASIS_RE.match(line)
        if match:
            place = int(match.group("place"))
            local_bases.append(
                {
                    "place": "infinity" if place == -1 else str(place),
                    "allowed_subspace_basis_columns_gp": match.group("basis"),
                }
            )
    return summary, deleted, basis, norm_basis, local_bases


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--gp", type=Path, required=True)
    parser.add_argument("--bnf-checkpoint", type=Path, required=True)
    parser.add_argument("--simon-directory", type=Path, required=True)
    parser.add_argument("--transformed-model", type=parse_model, required=True)
    parser.add_argument("--curve-theta", required=True)
    parser.add_argument("--known-generic-rank", type=int, default=17)
    parser.add_argument("--known-total-rank-lower-bound", type=int, default=29)
    parser.add_argument("--rigid-dimension", type=int, default=2)
    parser.add_argument("--timeout-seconds", type=float, default=1800.0)
    parser.add_argument("--rss-bytes", type=int, default=5_000_000_000)
    parser.add_argument("--stack-bytes", type=int, default=4_000_000_000)
    parser.add_argument("--random-seed", type=int, default=20260904)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.timeout_seconds <= 0 or args.stack_bytes <= 0 or args.rss_bytes <= 0:
        parser.error("resource limits must be positive")
    if not 0 <= args.known_generic_rank <= args.known_total_rank_lower_bound:
        parser.error("known rank bounds are inconsistent")
    if args.rigid_dimension < 0:
        parser.error("the rigid dimension must be nonnegative")

    gp = args.gp.resolve()
    checkpoint = args.bnf_checkpoint.resolve()
    simon = args.simon_directory.resolve()
    for path in (gp, checkpoint):
        if not path.is_file():
            raise FileNotFoundError(path)
    simon_files = [simon / name for name in ("ellQ.gp", "ell.gp", "qfsolve.gp", "resultant3.gp")]
    for path in simon_files:
        if not path.is_file():
            raise FileNotFoundError(path)
    for path in (args.log, args.output):
        if path.exists() and not args.overwrite:
            raise FileExistsError(path)
        path.parent.mkdir(parents=True, exist_ok=True)

    model = "[" + ",".join(args.transformed_model) + "]"
    simon_source = "\n".join(
        pari219_compatible_simon_source(path) for path in simon_files
    )
    definitions = SIMON_GP_FUNCTION.replace(
        "/* ELKIES_R17_GP_DEFINITION_SPLIT */", ""
    )
    program = f'''default(nbthreads,1);
setrand({args.random_seed});
{simon_source}
DEBUGLEVEL_ell=0;LIMBIGPRIME=0;LIM1=0;LIM3=0;LIMTRIV=0;
{definitions}
b=read("{gp_quote(checkpoint)}");
if(!bnfcertify(b),error("reloaded BNF failed certification"));
curve=ellinit({model});ctheta={args.curve_theta};
r=ell2selmer_basis_gen(curve,b,1,ctheta);
print("{PROTOCOL}|stage=summary|selmer=",#r[2],"|sclass=",r[5],"|norm=",r[6],"|local_rank=",r[7],"|bad=",r[3]);
print("{PROTOCOL}|stage=norm_basis|basis=",lift(r[10]));
for(i=1,#r[8],print("{PROTOCOL}|stage=local_basis|place=",r[9][i],"|basis=",lift(r[8][i])));
for(i=1,#r[4],a=r[4][i];print("{PROTOCOL}|stage=delete_one|place=",a[1],"|allowed=",a[2],"|alone=",a[3],"|omitted=",a[4],"|rank=",a[5]));
for(i=1,#r[2],a=Mod(1,b.pol);for(j=1,#r[1],if(r[2][j,i],a*=r[1][j]));print("{PROTOCOL}|stage=basis|index=",i,"|alpha=",lift(a)));
print("{PROTOCOL}|stage=complete|status=PASS");
'''

    source_path = args.output.with_suffix(".gp")
    preserve_previous(source_path)
    source_path.write_text(program)
    preserve_previous(args.output)
    supervision = run([str(gp), "-q", "-f", "-s", str(args.stack_bytes)],
        input_text=program, log_path=args.log,
        checkpoint_path=args.output.with_suffix(".supervisor.json"),
        limits=Limits(args.timeout_seconds, args.rss_bytes, pari_stack_bytes=args.stack_bytes))
    elapsed = supervision["wall_seconds"]
    outcome = "running" if supervision["outcome"] == "completed" else supervision["outcome"]
    log_text = args.log.read_text(errors="replace")
    summary, deleted, basis, norm_basis, local_bases = parse_log(log_text)
    complete = (
        outcome == "running"
        and supervision["returncode"] == 0
        and f"{PROTOCOL}|stage=complete|status=PASS" in log_text
        and summary is not None
        and len(basis) == summary["two_selmer_dimension"]
        and norm_basis is not None
        and len(local_bases) == len(deleted)
    )
    if outcome == "running":
        outcome = "completed_unconditional_two_selmer" if complete else "backend_failure"

    arithmetic = None
    if complete and summary is not None:
        total = int(summary["two_selmer_dimension"])
        residual = total - args.known_generic_rank
        if total < args.known_total_rank_lower_bound:
            raise ArithmeticError("Selmer dimension contradicts the certified rank lower bound")
        full_rank = int(summary["full_local_condition_matrix_rank_on_norm_subspace"])
        for row in deleted:
            row["rank_drop_from_full_matrix"] = full_rank - int(
                row["matrix_rank_after_deleting_this_place"]
            )
        arithmetic = {
            **summary,
            "residual_dimension_modulo_known_generic_mw": residual,
            "dimension_after_quotienting_rigid_plane": residual
            - args.rigid_dimension,
            "additional_dimension_beyond_known_total_rank_lower_bound": total
            - args.known_total_rank_lower_bound,
            "known_points_exhaust_two_selmer": total
            == args.known_total_rank_lower_bound,
            "rank_nullity_verified": full_rank + total
            == int(summary["global_norm_square_subspace_dimension"]),
            "global_norm_square_subspace_basis_columns_gp": norm_basis,
            "local_allowed_subspaces": local_bases,
            "delete_one_place": deleted,
            "basis": basis,
        }

    result = {
        "schema": SCHEMA,
        "status": outcome,
        "case_id": args.case_id,
        "arithmetic": arithmetic,
        "claim_boundary": [
            "Only completed_unconditional_two_selmer is a complete Selmer result.",
            "A timeout or backend failure supplies no Selmer or rank upper bound.",
            "If extra Selmer dimensions exist, this record does not align them with the known-point quotient basis.",
        ],
        "input": {
            "transformed_model": args.transformed_model,
            "curve_theta": args.curve_theta,
            "known_generic_rank": args.known_generic_rank,
            "known_total_rank_lower_bound": args.known_total_rank_lower_bound,
            "rigid_dimension": args.rigid_dimension,
            "bnf_checkpoint": str(checkpoint),
            "bnf_checkpoint_sha256": file_sha256(checkpoint),
            "simon_files": {
                str(path): file_sha256(path) for path in simon_files
            },
            "random_seed": args.random_seed,
            "timeout_seconds": args.timeout_seconds,
            "stack_bytes": args.stack_bytes,
            "rss_bytes": args.rss_bytes,
        },
        "backend": {"gp": str(gp), "version": gp_version(gp)},
        "measurement": {
            "wall_seconds": elapsed,
            "returncode": supervision["returncode"],
        },
        "supervision": supervision,
        "source": str(source_path),
        "source_sha256": file_sha256(source_path),
        "log": str(args.log),
        "log_sha256": file_sha256(args.log),
    }
    write_checkpoint(args.output, result)
    print(
        f"{PROTOCOL}|case={args.case_id}|status={outcome}|seconds={elapsed:.3f}"
        f"|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
