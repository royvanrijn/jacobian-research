#!/usr/bin/env python3
"""Run checkpointed norm-eight A1/MW16 screens for ordered ICARM targets.

For each target this driver runs the reusable Sage modular screen through a
fixed prime chain, stopping only when no survivor remains or the declared
chain is exhausted.  It then writes a compact standalone exclusion ledger and
factors every survivor over QQ.  Target order is preserved, making
``--curve-ids 302,273,542,548`` the canonical first campaign.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
SCREEN = ROOT / "elkies-k3/scripts/screen_icarm_norm8_a1_fibrations.sage"
EXACT = ROOT / "elkies-k3/scripts/certify_icarm_norm8_a1_survivors.sage"
MODEL = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
TABLE = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
)
TARGET_SNAPSHOT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-icarm-database-sweep-v2.json"
)
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
PRIMES = (
    1009,
    1013,
    1019,
    1021,
    1031,
    1033,
    1039,
    1049,
    1051,
    1061,
    1063,
    1069,
    1087,
    1091,
    1093,
    1097,
    1103,
    1109,
    1117,
    1123,
    1129,
    1151,
    1153,
    1163,
    1171,
    1181,
    1187,
    1193,
    1201,
    1213,
    1217,
    1223,
    1229,
    1231,
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def parse_ids(text: str) -> list[int]:
    values = [int(value) for value in text.split(",") if value.strip()]
    if not values or any(value <= 0 for value in values) or len(values) != len(set(values)):
        raise argparse.ArgumentTypeError("curve ids must be a nonempty comma-separated unique list")
    return values


def parse_primes(text: str) -> list[int]:
    values = [int(value) for value in text.split(",") if value.strip()]
    if not values or any(value <= 3 for value in values) or len(values) != len(set(values)):
        raise argparse.ArgumentTypeError("primes must be a nonempty comma-separated unique list")
    return values


def normalized_ainvs(values) -> tuple[Fraction, ...]:
    if len(values) == 2:
        values = [0, 0, 0, *values]
    if len(values) != 5:
        raise ValueError("target curve has neither two nor five a-invariants")
    return tuple(Fraction(value) for value in values)


def target_invariants(record):
    a1, a2, a3, a4, a6 = normalized_ainvs(record["ainvs"])
    b2 = a1**2 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    delta = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    return c4, c6, delta


def nonzero_mod(value: Fraction, prime: int) -> bool:
    return value.denominator % prime != 0 and value.numerator % prime != 0


def run(command: list[str]) -> None:
    subprocess.run(command, cwd=ROOT, check=True)


def compact_campaign(
    *,
    curve_id: int,
    target: dict,
    checkpoints: list[Path],
    output: Path,
    requested_primes: list[int],
) -> dict:
    if not checkpoints:
        raise ArithmeticError("campaign produced no usable-prime checkpoints")
    class_count = None
    prior_path = None
    survivors = None
    witness_prime = None
    summaries = []
    for path in checkpoints:
        document = json.loads(path.read_text())
        if document.get("status") != "PASS_COMPLETE_DECLARED_CHUNK_MODULAR_SCREEN":
            raise ArithmeticError(f"incomplete checkpoint {path}")
        if int(document["target"]["curve_id"]) != curve_id:
            raise ArithmeticError(f"checkpoint {path} belongs to a different target")
        search = document["search"]
        if class_count is None:
            class_count = int(search["priority_table_class_count"])
            if search["processed_priority_ranks"] != list(range(1, class_count + 1)):
                raise ArithmeticError("first checkpoint is not the complete priority table")
            witness_prime = [None] * class_count
        elif search["processed_priority_ranks"] != survivors:
            raise ArithmeticError(f"checkpoint chain is discontinuous at {path}")
        if prior_path is not None and search["candidate_rank_filter_sha256"] != digest(prior_path):
            raise ArithmeticError(f"predecessor hash changed at {path}")
        excluded = [
            int(record["priority_rank"])
            for record in document["records"]
            if record["status"] == "PASS_MODULAR_NO_TARGET_PARAMETER"
        ]
        prime = int(search["prime"])
        for rank in excluded:
            if witness_prime[rank - 1] is not None:
                raise ArithmeticError("a priority rank received two first-witness primes")
            witness_prime[rank - 1] = prime
        survivors = list(map(int, search["survivor_priority_ranks"]))
        summaries.append(
            {
                "prime": prime,
                "processed_count": int(search["processed_count"]),
                "excluded_count": len(excluded),
                "survivor_count": len(survivors),
                "checkpoint": relative(path),
                "checkpoint_sha256": digest(path),
            }
        )
        prior_path = path

    assert class_count is not None and witness_prime is not None and survivors is not None
    excluded_count = sum(value is not None for value in witness_prime)
    if excluded_count + len(survivors) != class_count:
        raise ArithmeticError("first-witness vector and survivors do not partition the atlas")
    if any(witness_prime[rank - 1] is not None for rank in survivors):
        raise ArithmeticError("survivor also has an exclusion witness")
    status = (
        "PASS_COMPLETE_MODULAR_EXCLUSION_OF_COMMITTED_LAYER"
        if not survivors
        else "PASS_COMPLETE_MODULAR_SCREEN_WITH_SURVIVORS"
    )
    payload = {
        "schema": "elkies-k3.icarm-norm8-a1-modular-screen-compact.v1",
        "status": status,
        "source_chart": "norm12-orbit-11952",
        "target": {
            "curve_id": curve_id,
            "label": f"ICARM curve {curve_id}",
            "snapshot_rank_lower_bound": int(target["snapshot_rank_lower_bound"]),
            "generalized_weierstrass_coefficients": list(target["ainvs"]),
        },
        "search": {
            "priority_table_class_count": class_count,
            "requested_prime_chain": requested_primes,
            "usable_prime_chain": [row["prime"] for row in summaries],
            "checkpoint_summaries": summaries,
            "excluded_count": excluded_count,
            "survivor_count": len(survivors),
            "survivor_priority_ranks": survivors,
            "first_exclusion_prime_by_priority_rank": witness_prime,
        },
        "proof_boundary": (
            "Each non-null first-exclusion entry names a good prime at which the exact "
            "projective target j-equation has no root. Survivors remain UNKNOWN until "
            "characteristic-zero factorization. Completeness applies only to the committed "
            "63,917 old-degree-two norm-eight classes on source chart 11952."
        ),
        "inputs": {
            relative(path): digest(path)
            for path in (Path(__file__).resolve(), MODEL, TABLE, TARGET_SNAPSHOT, SCREEN)
        },
        "reproducing_command": (
            "python3 elkies-k3/scripts/run_icarm_norm8_a1_atlas.py "
            f"--curve-ids {curve_id} --primes {','.join(map(str, requested_primes))}"
        ),
        "software_assumptions": {"python": sys.version.split()[0], "sage": "10.9"},
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curve-ids", type=parse_ids, required=True)
    parser.add_argument("--primes", type=parse_primes, default=list(PRIMES))
    parser.add_argument("--max-primes", type=int)
    parser.add_argument("--local-directory", type=Path, default=LOCAL)
    parser.add_argument("--generated-directory", type=Path, default=GENERATED)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--skip-exact", action="store_true")
    args = parser.parse_args()
    if args.max_primes is not None and args.max_primes <= 0:
        parser.error("--max-primes must be positive")

    requested_primes = list(args.primes)
    if args.max_primes is not None:
        requested_primes = requested_primes[: args.max_primes]
    local = args.local_directory.resolve()
    generated = args.generated_directory.resolve()
    local.mkdir(parents=True, exist_ok=True)
    generated.mkdir(parents=True, exist_ok=True)

    snapshot = json.loads(TARGET_SNAPSHOT.read_text())
    target_by_id = {
        int(record["id"]): record for record in snapshot["snapshot"]["curves"]
    }
    summaries = []
    for curve_id in args.curve_ids:
        if curve_id not in target_by_id:
            raise ValueError(f"curve {curve_id} is absent from the pinned target snapshot")
        target = target_by_id[curve_id]
        c4, c6, delta = target_invariants(target)
        usable_primes = [
            prime
            for prime in requested_primes
            if nonzero_mod(c4, prime) and nonzero_mod(c6, prime) and nonzero_mod(delta, prime)
        ]
        if not usable_primes:
            raise ArithmeticError(f"curve {curve_id} has no usable prime in the declared chain")

        checkpoints = []
        prior = None
        survivors = None
        for prime in usable_primes:
            output = local / f"icarm-curve{curve_id}-11952-norm8-a1-mod{prime}.json"
            command = [
                "sage",
                "-python",
                str(SCREEN),
                "--curve-id",
                str(curve_id),
                "--prime",
                str(prime),
                "--output",
                str(output),
            ]
            if prior is not None:
                command.extend(["--candidate-ranks", str(prior)])
            reuse = False
            if args.resume and output.is_file():
                document = json.loads(output.read_text())
                search = document.get("search", {})
                reuse = (
                    document.get("status") == "PASS_COMPLETE_DECLARED_CHUNK_MODULAR_SCREEN"
                    and int(document.get("target", {}).get("curve_id", -1)) == curve_id
                    and int(search.get("prime", -1)) == prime
                    and (
                        (prior is None and search.get("candidate_rank_filter_sha256") is None)
                        or (
                            prior is not None
                            and search.get("candidate_rank_filter_sha256") == digest(prior)
                        )
                    )
                )
            if not reuse:
                run(command)
            document = json.loads(output.read_text())
            survivors = list(map(int, document["search"]["survivor_priority_ranks"]))
            checkpoints.append(output)
            prior = output
            if not survivors:
                break

        compact_output = (
            generated
            / f"elkies-k3-icarm-curve{curve_id}-11952-norm8-a1-modular-screen-v1.json"
        )
        compact = compact_campaign(
            curve_id=curve_id,
            target=target,
            checkpoints=checkpoints,
            output=compact_output,
            requested_primes=requested_primes,
        )
        exact_output = (
            generated
            / f"elkies-k3-icarm-curve{curve_id}-11952-norm8-a1-exact-survivors-v1.json"
        )
        exact_status = "SKIPPED"
        qq_hits = None
        if not args.skip_exact:
            run(
                [
                    "sage",
                    "-python",
                    str(EXACT),
                    "--curve-id",
                    str(curve_id),
                    "--candidate-ranks",
                    str(compact_output),
                    "--output",
                    str(exact_output),
                ]
            )
            exact = json.loads(exact_output.read_text())
            exact_status = exact["status"]
            qq_hits = int(exact["qq_isomorphic_candidate_count"])
        summary = {
            "curve_id": curve_id,
            "rank_lower_bound": int(target["snapshot_rank_lower_bound"]),
            "usable_primes_run": len(checkpoints),
            "modular_survivors": int(compact["search"]["survivor_count"]),
            "exact_status": exact_status,
            "qq_isomorphic_candidates": qq_hits,
            "compact_output": relative(compact_output),
            "exact_output": None if args.skip_exact else relative(exact_output),
        }
        summaries.append(summary)
        print(
            f"ICARMA1CAMPAIGN|curve={curve_id}|primes={len(checkpoints)}"
            f"|modular_survivors={summary['modular_survivors']}"
            f"|qq_isomorphic={qq_hits}|status={exact_status}",
            flush=True,
        )

    print("ICARMA1CAMPAIGN|summary=" + json.dumps(summaries, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
