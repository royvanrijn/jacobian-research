#!/usr/bin/env python3
"""Build a resumable exact bad-place Kummer ledger for Elkies rank 28.

Each finite or infinite place is handled by a fresh, owned Sage worker.  A
completed block is written immediately and replayed on later runs; a timeout
therefore identifies one specific place and does not discard completed work.
The ledger concerns the images of the certified generic seventeen points.  It
does not compute an S-class group, enumerate the ambient K(S,2), or assert a
Selmer upper bound.
"""

from __future__ import annotations

from research_runtime.pari_context import prepared_prime_ideals, prepared_factor
from research_runtime.supervisor import Limits, capture_record

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
DEFAULT_SPECIALIZATION = (
    ROOT
    / "artifacts/local/elliptic-curves/q12o5867-specializations"
    / "q12o5867-specialization-m9529_5471.json"
)
DEFAULT_CACHE = (
    ROOT
    / "artifacts/local/elliptic-curves/q12o5867-bnf-free"
    / "elkies-rank28-bad-place-blocks"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_rank28_bad_place_kummer_ledger_v1.json"
)
DEFAULT_SAGE = Path("/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python")
SCHEMA = "elliptic-curves.elkies-2026-rank28-bad-place-kummer-ledger.v1"
BLOCK_SCHEMA = "elliptic-curves.elkies-2026-rank28-local-kummer-block.v1"
PROTOCOL = "ELKIESR28LOCAL"

# Complete factorization of the monic 2-division cubic discriminant.  The
# final two factors were discovered through FactorDB and then independently
# proved prime by Sage; the proof is replayed by this program before workers
# are launched.
DISCRIMINANT_FACTORIZATION = (
    (2, 23),
    (3, 6),
    (5, 6),
    (7, 4),
    (11, 2),
    (13, 4),
    (17, 5),
    (19, 3),
    (48463, 1),
    (20650099, 1),
    (315574902691581877528345013999136728634663121, 1),
    (
        376018840263193489397987439236873583997122096511452343225772113000611087671413,
        1,
    ),
)
ODD_BAD_PRIMES = tuple(prime for prime, _ in DISCRIMINANT_FACTORIZATION if prime != 2)

sys.path[:0] = [str(ELLIPTIC_ROOT), str(CAS)]

from build_q12o5867_bnf_free_signature import (  # noqa: E402
    evaluate_cubic,
    monic_cubic_coefficients,
    point_on_monic_cubic,
)
from elliptic_candidate_record import is_on_weierstrass_curve  # noqa: E402
from run_fermigier_rank20_auxiliary_fingerprints import (  # noqa: E402
    f2_rank,
    prime_local_rows,
    qpari,
    two_adic_coords,
)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def factorization_product() -> int:
    answer = 1
    for prime, exponent in DISCRIMINANT_FACTORIZATION:
        answer *= prime**exponent
    return answer


def load_input(path: Path) -> dict[str, Any]:
    artifact = json.loads(path.read_text())
    if artifact.get("status") != "PASS_EXACT_Q12O5867_SPECIALIZED_GENERIC_RANK17_LOWER_BOUND":
        raise ValueError("input is not an exact q12o5867 rank-17 specialization")
    parameter = artifact.get("parameter", {}).get("affine_value")
    if parameter != "-9529/5471":
        raise ValueError("input is not the published rank-28 fibre")
    minimal = artifact["global_minimal_specialization"]
    model = tuple(minimal["model"])
    points = tuple(tuple(point) for point in minimal["points"])
    if len(points) != 17 or len(set(points)) != 17:
        raise ValueError("input does not contain seventeen distinct generic points")
    if any(not is_on_weierstrass_curve(model, point) for point in points):
        raise ArithmeticError("a serialized generic point misses the minimal model")
    coefficients = monic_cubic_coefficients(model)
    transformed = tuple(point_on_monic_cubic(model, point) for point in points)
    if any(z * z != evaluate_cubic(coefficients, x) for x, z in transformed):
        raise ArithmeticError("the monic-cubic point transport failed")
    discriminant = factorization_product()
    # The closed formula is evaluated independently of the pinned factors.
    d, c, b, one = coefficients
    if one != 1:
        raise ArithmeticError("the descent cubic is not monic")
    computed = b * b * c * c - 4 * c**3 - 4 * b**3 * d - 27 * d * d + 18 * b * c * d
    if computed != discriminant:
        raise ArithmeticError("the pinned discriminant factorization has the wrong product")
    return {
        "artifact": artifact,
        "model": model,
        "points": points,
        "coefficients": coefficients,
        "discriminant": discriminant,
        "sha256": file_sha256(path),
    }


def _worker_setup(path: Path):
    from fractions import Fraction
    from sage.all import PolynomialRing, QQ, ZZ, pari

    source = load_input(path)
    coefficients = source["coefficients"]
    ring = PolynomialRing(QQ, "t")
    t = ring.gen()
    polynomial = sum(ZZ(value) * t**index for index, value in enumerate(coefficients))
    ramified_primes = [prime for prime, _ in DISCRIMINANT_FACTORIZATION]
    print(
        f"{PROTOCOL}|stage=nfinit|status=start|factorization_supplied=true",
        flush=True,
    )
    started = time.monotonic()
    # Supplying every divisor of the polynomial discriminant prevents PARI
    # from repeating the hard 168-digit factorization inside nfinit.
    from research_runtime.pari_context import prepared_nf
    nf = prepared_nf(pari(polynomial), ramified_primes)
    print(
        f"{PROTOCOL}|stage=nfinit|status=complete|seconds={time.monotonic()-started:.6f}",
        flush=True,
    )
    print(f"{PROTOCOL}|stage=nfcertify|status=start", flush=True)
    certify_started = time.monotonic()
    obstructions = list(pari.nfcertify(nf))
    if obstructions:
        raise ArithmeticError(
            f"factor-supplied nfinit did not certify the maximal order: {obstructions}"
        )
    print(
        f"{PROTOCOL}|stage=nfcertify|status=complete|seconds={time.monotonic()-certify_started:.6f}",
        flush=True,
    )
    theta = pari(f"Mod(t,{polynomial})")
    transformed = tuple(
        point_on_monic_cubic(source["model"], point) for point in source["points"]
    )
    alphas = [qpari(pari, Fraction(x)) - theta for x, _ in transformed]
    return source, polynomial, nf, alphas, transformed, pari


def worker_proof(path: Path) -> dict[str, Any]:
    from sage.all import ZZ

    source = load_input(path)
    rows = []
    for prime, exponent in DISCRIMINANT_FACTORIZATION:
        started = time.monotonic()
        proved = bool(ZZ(prime).is_prime(proof=True))
        rows.append(
            {
                "prime": str(prime),
                "exponent": exponent,
                "proved_prime": proved,
                "proof_seconds": time.monotonic() - started,
            }
        )
        if not proved:
            raise ArithmeticError(f"factor {prime} failed its primality proof")
    return {
        "schema": BLOCK_SCHEMA,
        "kind": "factorization_proof",
        "specialization_sha256": source["sha256"],
        "discriminant": str(source["discriminant"]),
        "factorization": rows,
        "product_verified": True,
        "all_factors_proved_prime": True,
    }


def worker_odd(path: Path, prime: int) -> dict[str, Any]:
    if prime not in ODD_BAD_PRIMES:
        raise ValueError("requested prime is not in the complete odd bad-place set")
    source, _polynomial, nf, alphas, _transformed, pari = _worker_setup(path)
    print(f"{PROTOCOL}|stage=odd_place|prime={prime}|status=start", flush=True)
    started = time.monotonic()
    rows, places = prime_local_rows(pari, nf, alphas, prime)
    return {
        "schema": BLOCK_SCHEMA,
        "kind": "odd_bad_prime",
        "rational_prime": str(prime),
        "specialization_sha256": source["sha256"],
        "seconds": time.monotonic() - started,
        "prime_ideals": list(places),
        "coordinate_order_per_prime_ideal": [
            "valuation_parity",
            "unit_residue_squareclass",
        ],
        "rows_for_generic_P1_through_P17": rows,
    }


def worker_two(path: Path) -> dict[str, Any]:
    source, _polynomial, nf, alphas, transformed, pari = _worker_setup(path)
    print(f"{PROTOCOL}|stage=two_adic|status=start", flush=True)
    started = time.monotonic()
    two_primes = list(prepared_prime_ideals(nf, 2))
    basis, origins, rows = two_adic_coords(pari, nf, two_primes, alphas)
    return {
        "schema": BLOCK_SCHEMA,
        "kind": "two_adic",
        "rational_prime": "2",
        "specialization_sha256": source["sha256"],
        "seconds": time.monotonic() - started,
        "prime_ideals": [str(prime) for prime in two_primes],
        "basis_dimension": len(basis),
        "basis_origins_one_based": [index + 1 for index in origins],
        "basis_generators": [
            [str(transformed[index][0]), "-1", "0"] for index in origins
        ],
        "rows_for_generic_P1_through_P17": rows,
    }


def worker_real(path: Path) -> dict[str, Any]:
    from fractions import Fraction
    from sage.all import AA, PolynomialRing, QQ, ZZ

    source = load_input(path)
    ring = PolynomialRing(QQ, "t")
    t = ring.gen()
    polynomial = sum(ZZ(value) * t**index for index, value in enumerate(source["coefficients"]))
    roots = list(polynomial.roots(AA, multiplicities=False))
    transformed = tuple(
        point_on_monic_cubic(source["model"], point) for point in source["points"]
    )
    rows = [
        [1 if QQ(Fraction(x)) - root < 0 else 0 for root in roots]
        for x, _ in transformed
    ]
    return {
        "schema": BLOCK_SCHEMA,
        "kind": "real_places",
        "specialization_sha256": source["sha256"],
        "real_root_count": len(roots),
        "root_order": "increasing_real_root",
        "rows_for_generic_P1_through_P17": rows,
    }


def emit_worker(args: argparse.Namespace) -> None:
    if args.worker_kind == "proof":
        result = worker_proof(args.specialization)
    elif args.worker_kind == "odd":
        if args.worker_prime is None:
            raise ValueError("odd worker requires --worker-prime")
        result = worker_odd(args.specialization, args.worker_prime)
    elif args.worker_kind == "two":
        result = worker_two(args.specialization)
    elif args.worker_kind == "real":
        result = worker_real(args.specialization)
    else:
        raise ValueError("unknown worker kind")
    print(f"{PROTOCOL}_JSON=" + json.dumps(result, sort_keys=True), flush=True)


def owned_worker(
    *,
    sage_python: Path,
    specialization: Path,
    kind: str,
    prime: int | None,
    timeout: float,
) -> dict[str, Any]:
    command = [
        str(sage_python),
        str(Path(__file__).resolve()),
        "--specialization",
        str(specialization.resolve()),
        "--worker-kind",
        kind,
    ]
    if prime is not None:
        command.extend(("--worker-prime", str(prime)))
    record=capture_record(command,limits=Limits(timeout,1_073_741_824))
    outcome=record['outcome']
    stdout,stderr=record['stdout'],record['stderr']
    marker = f"{PROTOCOL}_JSON="
    payloads = [line[len(marker) :] for line in stdout.splitlines() if line.startswith(marker)]
    block = json.loads(payloads[0]) if outcome == "completed" and len(payloads) == 1 else None
    if outcome == "completed" and (record["returncode"] != 0 or block is None):
        outcome = "backend_failure"
    return {
        "outcome": outcome,
        "returncode": record["returncode"],
        "wall_seconds": record["wall_seconds"],
        "command": command,
        "stdout": stdout,
        "stderr": stderr,
        "block": block,
    }


def cache_name(kind: str, prime: int | None) -> str:
    return f"odd-{prime}.json" if kind == "odd" else f"{kind}.json"


def valid_cached(path: Path, *, kind: str, prime: int | None, source_sha: str):
    if not path.is_file():
        return None
    block = json.loads(path.read_text())
    if block.get("schema") != BLOCK_SCHEMA or block.get("specialization_sha256") != source_sha:
        return None
    expected_kind = "odd_bad_prime" if kind == "odd" else {
        "proof": "factorization_proof",
        "two": "two_adic",
        "real": "real_places",
    }[kind]
    if block.get("kind") != expected_kind:
        return None
    if prime is not None and block.get("rational_prime") != str(prime):
        return None
    return block


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--specialization", type=Path, default=DEFAULT_SPECIALIZATION)
    parser.add_argument("--cache-directory", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--sage-python", type=Path, default=DEFAULT_SAGE)
    parser.add_argument("--per-place-timeout", type=float, default=120.0)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--worker-kind", choices=("proof", "odd", "two", "real"), help=argparse.SUPPRESS
    )
    parser.add_argument("--worker-prime", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    sys.set_int_max_str_digits(0)
    if args.worker_kind:
        emit_worker(args)
        return
    if args.per_place_timeout <= 0:
        parser.error("--per-place-timeout must be positive")
    sage_python = shutil.which(str(args.sage_python))
    if sage_python is None:
        raise SystemExit(f"Sage Python is unavailable: {args.sage_python}")
    source = load_input(args.specialization)
    args.cache_directory.mkdir(parents=True, exist_ok=True)

    tasks = [("proof", None), *(("odd", prime) for prime in ODD_BAD_PRIMES), ("two", None), ("real", None)]
    records = []
    blocks = []
    for kind, prime in tasks:
        cache_path = args.cache_directory / cache_name(kind, prime)
        cached = valid_cached(
            cache_path, kind=kind, prime=prime, source_sha=source["sha256"]
        )
        if cached is not None:
            records.append(
                {
                    "kind": kind,
                    "prime": str(prime) if prime is not None else None,
                    "outcome": "replayed_cache",
                    "cache_path": str(cache_path.resolve()),
                    "cache_sha256": file_sha256(cache_path),
                }
            )
            blocks.append(cached)
            continue
        result = owned_worker(
            sage_python=Path(sage_python),
            specialization=args.specialization,
            kind=kind,
            prime=prime,
            timeout=args.per_place_timeout,
        )
        block = result.pop("block")
        if block is not None:
            cache_path.write_text(json.dumps(block, indent=2, sort_keys=True) + "\n")
            result["cache_path"] = str(cache_path.resolve())
            result["cache_sha256"] = file_sha256(cache_path)
            blocks.append(block)
        records.append(
            {
                "kind": kind,
                "prime": str(prime) if prime is not None else None,
                **result,
            }
        )
        print(
            f"{PROTOCOL}|kind={kind}|prime={prime}|outcome={result['outcome']}",
            flush=True,
        )

    proof = next((block for block in blocks if block["kind"] == "factorization_proof"), None)
    local_blocks = [block for block in blocks if block["kind"] != "factorization_proof"]
    expected_local = len(ODD_BAD_PRIMES) + 2
    complete = proof is not None and len(local_blocks) == expected_local
    concatenated = [[] for _ in range(17)]
    if complete:
        for block in local_blocks:
            rows = block["rows_for_generic_P1_through_P17"]
            if len(rows) != 17:
                raise ArithmeticError("a completed local block has the wrong row count")
            concatenated = [concatenated[index] + rows[index] for index in range(17)]
    document = {
        "schema": SCHEMA,
        "status": (
            "COMPLETE_ALL_BAD_PLACE_KUMMER_IMAGES_NOT_A_SELMER_BOUND"
            if complete
            else "PARTIAL_BAD_PLACE_KUMMER_IMAGES_NOT_A_SELMER_BOUND"
        ),
        "parameter": "-9529/5471",
        "specialization": {
            "path": str(args.specialization.resolve()),
            "sha256": source["sha256"],
            "global_minimal_model": [str(value) for value in source["model"]],
            "generic_point_count": 17,
        },
        "descent_cubic_coefficients_ascending": [
            str(value) for value in source["coefficients"]
        ],
        "descent_cubic_discriminant": str(source["discriminant"]),
        "factorization": [
            {"prime": str(prime), "exponent": exponent}
            for prime, exponent in DISCRIMINANT_FACTORIZATION
        ],
        "factorization_product_verified": True,
        "factor_primality_proof_completed": proof is not None,
        "expected_local_block_count": expected_local,
        "completed_local_block_count": len(local_blocks),
        "all_bad_place_blocks_completed": complete,
        "combined_known_kummer_rank": f2_rank(concatenated) if complete else None,
        "task_records": records,
        "completed_blocks": blocks,
        "claim_boundary": [
            "Every completed block is an exact local image of the certified generic seventeen points.",
            "A timeout or missing block is computational incompleteness, not a local obstruction.",
            "Even a complete known-image ledger is not an ambient K(S,2), S-class-group, local-solubility, Selmer, or rank upper-bound certificate.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.overwrite else "x"
    with args.output.open(mode) as handle:
        json.dump(document, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(
        f"{PROTOCOL}|status={document['status']}|"
        f"completed={len(local_blocks)}/{expected_local}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
