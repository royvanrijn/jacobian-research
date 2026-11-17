#!/usr/bin/env python3
"""Discover modular quotient identities for the q-period input minors.

Let I=(f_1,...,f_16) be the primitive maximal-minor ideal and let
G=(g_1,...,g_21) be the rebuilt rational candidate basis.  At each requested
good prime this script performs Singular's ordered division

    f_i = sum_j q_{j,i} g_j + r_i

and records every coefficient of every q_{j,i}.  The ordered basis and term
order make this a reproducible modular quotient convention.  This is a
reconstruction checkpoint, not a characteristic-zero containment proof:
that proof requires CRT/rational reconstruction followed by exact integer
polynomial identity checks.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from math import gcd, isqrt
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from verify_degree_five_qper_fitting_basis import (
    DEFAULT_BASIS,
    ROOT,
    frozen_basis_lines,
    primitive_minor_lines,
)

try:
    from flint import fmpz
except ImportError:  # pragma: no cover - dependency is present in the repo venv
    fmpz = None


DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree_five_qper_input_quotients_modular.json"
)
DEFAULT_PRIMES = (31991, 32003, 65521)
LEGACY_FORMAT = "degree-five-qper-input-quotients-modular-v1"
FORMAT = "degree-five-qper-input-quotients-crt-v2"


def modular_program(prime: int, basis_path: Path) -> str:
    lines = [
        f"ring r={prime},(a,tau),dp;",
        "short=0;",
        "ideal G;",
        *frozen_basis_lines(basis_path),
        'attrib(G,"isSB",1);',
        "int quotient_timer=timer;",
        *primitive_minor_lines(),
        "list input_division=division(primitive_minors,G);",
        "matrix input_quotients=input_division[1];",
        "ideal input_remainders=input_division[2];",
        (
            "matrix input_residual=matrix(primitive_minors)"
            "-matrix(G)*input_quotients;"
        ),
        "ideal input_residual_ideal=ideal(input_residual);",
        (
            'print("REMAINDER_TERMS="'
            "+string(size(input_remainders)));"
        ),
        (
            'print("RESIDUAL_TERMS="'
            "+string(size(input_residual_ideal)));"
        ),
        (
            'print("QUOTIENT_ROWS="'
            "+string(nrows(input_quotients)));"
        ),
        (
            'print("QUOTIENT_COLS="'
            "+string(ncols(input_quotients)));"
        ),
        "int input_index;",
        "int basis_index;",
        "poly quotient_entry;",
        "intvec quotient_exponent;",
        "for(input_index=1;input_index<=16;input_index++)",
        "{",
        "  for(basis_index=1;basis_index<=21;basis_index++)",
        "  {",
        (
            "    quotient_entry="
            "input_quotients[basis_index,input_index];"
        ),
        "    while(quotient_entry!=0)",
        "    {",
        "      quotient_exponent=leadexp(quotient_entry);",
        (
            '      print("QTERM="+string(input_index)'
            '+","+string(basis_index)'
            '+","+string(quotient_exponent[1])'
            '+","+string(quotient_exponent[2])'
            '+","+string(leadcoef(quotient_entry)));'
        ),
        "      quotient_entry=quotient_entry-lead(quotient_entry);",
        "    }",
        "  }",
        "}",
        (
            'print("QUOTIENT_SECONDS="'
            "+string(timer-quotient_timer));"
        ),
        "quit;",
        "",
    ]
    return "\n".join(lines)


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def primes_after(start: int, count: int) -> list[int]:
    primes = []
    candidate = start + 1
    while len(primes) < count:
        if is_prime(candidate):
            primes.append(candidate)
        candidate += 1
    return primes


def parse_output(prime: int, output: str) -> dict:
    scalars: dict[str, int] = {}
    quotients: list[list[list[list[int]]]] = [
        [[] for _ in range(21)]
        for _ in range(16)
    ]
    for line in output.splitlines():
        if line.startswith("QTERM="):
            values = [int(value) for value in line[6:].split(",")]
            if len(values) != 5:
                raise RuntimeError(f"malformed quotient term: {line}")
            input_index, basis_index, a_exp, tau_exp, coefficient = values
            quotients[input_index - 1][basis_index - 1].append(
                [a_exp, tau_exp, coefficient % prime]
            )
        elif "=" in line:
            name, value = line.split("=", 1)
            if name in {
                "REMAINDER_TERMS",
                "RESIDUAL_TERMS",
                "QUOTIENT_ROWS",
                "QUOTIENT_COLS",
                "QUOTIENT_SECONDS",
            }:
                scalars[name] = int(value)
    expected = {
        "REMAINDER_TERMS": 0,
        "RESIDUAL_TERMS": 0,
        "QUOTIENT_ROWS": 21,
        "QUOTIENT_COLS": 16,
    }
    for name, value in expected.items():
        if scalars.get(name) != value:
            raise RuntimeError(
                f"prime {prime}: expected {name}={value}, "
                f"got {scalars.get(name)}\n{output}"
            )
    support = [
        [
            [[a_exp, tau_exp] for a_exp, tau_exp, _ in entry]
            for entry in target
        ]
        for target in quotients
    ]
    support_encoding = json.dumps(
        support,
        separators=(",", ":"),
        sort_keys=True,
    ).encode()
    return {
        "prime": prime,
        "seconds": scalars["QUOTIENT_SECONDS"],
        "nonzero_entries": sum(
            bool(entry)
            for target in quotients
            for entry in target
        ),
        "term_count": sum(
            len(entry)
            for target in quotients
            for entry in target
        ),
        "support_sha256": hashlib.sha256(support_encoding).hexdigest(),
        "quotients": quotients,
    }


def run_prime(
    singular: str,
    prime: int,
    basis_path: Path,
    timeout: int,
) -> dict:
    with tempfile.TemporaryDirectory(
        prefix=f"degree-five-qper-input-{prime}-",
    ) as directory:
        program = Path(directory) / "quotients.sing"
        program.write_text(modular_program(prime, basis_path))
        try:
            result = subprocess.run(
                [singular, "-q", str(program)],
                check=True,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as error:
            raise RuntimeError(
                f"prime {prime} exceeded {timeout} seconds"
            ) from error
    output = result.stdout + result.stderr
    errors = [
        line
        for line in output.splitlines()
        if line.lstrip().startswith("?")
    ]
    if errors:
        raise RuntimeError(
            f"Singular failed at prime {prime}:\n"
            + "\n".join(errors[:20])
        )
    return parse_output(prime, output)


def rational_reconstruct(residue: int, modulus: int) -> tuple[int, int] | None:
    """Return the balanced rational reconstruction, when it exists."""
    if fmpz is None:
        bound = isqrt(modulus // 2)
        old_remainder, remainder = modulus, residue
        old_coefficient, coefficient = 0, 1
    else:
        old_remainder = fmpz(modulus)
        remainder = fmpz(residue)
        bound = (old_remainder // 2).isqrt()
        old_coefficient, coefficient = fmpz(0), fmpz(1)
    while remainder > bound:
        quotient = old_remainder // remainder
        old_remainder, remainder = (
            remainder,
            old_remainder - quotient * remainder,
        )
        old_coefficient, coefficient = (
            coefficient,
            old_coefficient - quotient * coefficient,
        )
    numerator, denominator = remainder, coefficient
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    if (
        denominator == 0
        or denominator > bound
        or abs(numerator) > bound
        or gcd(numerator, denominator) != 1
        or (numerator - residue * denominator) % modulus
    ):
        return None
    return int(numerator), int(denominator)


def flatten_record(record: dict) -> tuple[list[list[int]], list[int]]:
    support = []
    coefficients = []
    for input_index, target in enumerate(record["quotients"], start=1):
        for basis_index, entry in enumerate(target, start=1):
            for a_exp, tau_exp, coefficient in entry:
                support.append(
                    [input_index, basis_index, a_exp, tau_exp]
                )
                coefficients.append(coefficient)
    return support, coefficients


def merge_image(state: dict, prime: int, coefficients: list[int]) -> None:
    modulus = state["modulus"]
    inverse = pow(modulus, -1, prime)
    for index, coefficient in enumerate(coefficients):
        previous = state["residues"][index]
        correction = ((coefficient - previous) * inverse) % prime
        state["residues"][index] = previous + modulus * correction
    state["modulus"] *= prime
    state["crt_primes"].append(prime)


def empty_state() -> dict:
    return {
        "support_sha256": None,
        "support": [],
        "modulus": 1,
        "residues": [],
        "crt_primes": [],
        "holdout": None,
        "records": [],
    }


def accept_modular_record(state: dict, record: dict) -> str:
    support, coefficients = flatten_record(record)
    if state["support_sha256"] is None:
        state["support_sha256"] = record["support_sha256"]
        state["support"] = support
        state["residues"] = [0] * len(support)
    status = "good"
    if (
        record["support_sha256"] != state["support_sha256"]
        or support != state["support"]
    ):
        status = "unlucky"
    else:
        if state["holdout"] is not None:
            merge_image(
                state,
                state["holdout"]["prime"],
                state["holdout"]["coefficients"],
            )
        state["holdout"] = {
            "prime": record["prime"],
            "coefficients": coefficients,
        }
    state["records"].append(
        {
            key: record[key]
            for key in (
                "prime",
                "seconds",
                "nonzero_entries",
                "term_count",
                "support_sha256",
            )
        }
    )
    return status


def migrate_legacy(checkpoint: dict) -> dict:
    records = checkpoint["records"]
    support_counts = Counter(
        record["support_sha256"]
        for record in records
    )
    dominant_support, _ = support_counts.most_common(1)[0]
    ordered = [
        record
        for record in records
        if record["support_sha256"] == dominant_support
    ] + [
        record
        for record in records
        if record["support_sha256"] != dominant_support
    ]
    state = empty_state()
    for record in ordered:
        accept_modular_record(state, record)
    # Preserve the original record order for the audit ledger.
    summaries = {
        record["prime"]: record
        for record in state["records"]
    }
    state["records"] = [
        summaries[record["prime"]]
        for record in records
    ]
    return state


def load_state(checkpoint: dict) -> dict:
    if checkpoint["format"] == LEGACY_FORMAT:
        return migrate_legacy(checkpoint)
    if checkpoint["format"] != FORMAT:
        raise ValueError(f"unrecognized format {checkpoint['format']}")
    return {
        "support_sha256": checkpoint["support_sha256"],
        "support": checkpoint["support"],
        "modulus": int(checkpoint["crt"]["modulus"], 16),
        "residues": [
            int(value, 16)
            for value in checkpoint["crt"]["residues"]
        ],
        "crt_primes": checkpoint["crt"]["primes"],
        "holdout": checkpoint["holdout"],
        "records": checkpoint["records"],
    }


def state_payload(state: dict, basis: str) -> dict:
    return {
        "format": FORMAT,
        "basis": basis,
        "ordered_division": "Singular division(primitive_minors,G)",
        "support_sha256": state["support_sha256"],
        "support": state["support"],
        "crt": {
            "modulus": hex(state["modulus"]),
            "primes": state["crt_primes"],
            "residues": [
                hex(value)
                for value in state["residues"]
            ],
        },
        "holdout": state["holdout"],
        "records": state["records"],
    }


def reconstruction_diagnostic(state: dict) -> dict | None:
    """Rationally reconstruct the CRT state and test at the holdout prime."""
    holdout = state["holdout"]
    if not state["crt_primes"] or holdout is None:
        return None

    candidates = 0
    validated = 0
    validated_by_input: Counter[int] = Counter()
    maximum_numerator_bits = 0
    maximum_denominator_bits = 0
    holdout_prime = holdout["prime"]
    for index, residue in enumerate(state["residues"]):
        fraction = rational_reconstruct(residue, state["modulus"])
        if fraction is None:
            continue
        candidates += 1
        numerator, denominator = fraction
        if denominator % holdout_prime == 0:
            continue
        expected = (
            numerator * pow(denominator, -1, holdout_prime)
        ) % holdout_prime
        if expected == holdout["coefficients"][index]:
            validated += 1
            validated_by_input[state["support"][index][0]] += 1
            maximum_numerator_bits = max(
                maximum_numerator_bits,
                abs(numerator).bit_length(),
            )
            maximum_denominator_bits = max(
                maximum_denominator_bits,
                denominator.bit_length(),
            )
    return {
        "build_primes": len(state["crt_primes"]),
        "build_modulus_bits": state["modulus"].bit_length(),
        "holdout_prime": holdout_prime,
        "candidates": candidates,
        "validated": validated,
        "validated_by_input": dict(sorted(validated_by_input.items())),
        "maximum_numerator_bits": maximum_numerator_bits,
        "maximum_denominator_bits": maximum_denominator_bits,
        "total": len(state["residues"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--basis",
        type=Path,
        default=DEFAULT_BASIS,
        help="rebuilt rational basis artifact",
    )
    parser.add_argument(
        "--primes",
        type=int,
        nargs="+",
        help="good primes used for modular quotient discovery",
    )
    parser.add_argument(
        "--prime-start",
        type=int,
        default=1_000_000_000,
        help="generate primes strictly above this value",
    )
    parser.add_argument(
        "--prime-count",
        type=int,
        default=0,
        help="number of consecutive generated primes to process",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="checkpoint JSON containing modular quotient coefficients",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=600,
        help="hard wall-clock limit per prime",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="parallel Singular workers",
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=8,
        help="write the CRT checkpoint after this many completed primes",
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help="ignore an existing output checkpoint",
    )
    parser.add_argument(
        "--skip-diagnostic",
        action="store_true",
        help="skip balanced reconstruction against the held-out prime",
    )
    args = parser.parse_args()
    if args.timeout < 1:
        parser.error("--timeout must be positive")
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.checkpoint_every < 1:
        parser.error("--checkpoint-every must be positive")
    if args.primes is not None and args.prime_count:
        parser.error("--primes and --prime-count are mutually exclusive")
    if args.prime_count < 0:
        parser.error("--prime-count cannot be negative")
    if args.prime_count:
        requested_primes = primes_after(
            args.prime_start,
            args.prime_count,
        )
    elif args.primes is not None:
        requested_primes = args.primes
    else:
        requested_primes = list(DEFAULT_PRIMES)
    if len(set(requested_primes)) != len(requested_primes):
        parser.error("--primes must be distinct")
    nonprimes = [
        prime
        for prime in requested_primes
        if not is_prime(prime)
    ]
    if nonprimes:
        parser.error(
            "--primes contains nonprime values: "
            + ",".join(str(value) for value in nonprimes)
        )

    basis_path = args.basis
    if not basis_path.is_absolute():
        basis_path = ROOT / basis_path
    if not basis_path.exists():
        raise SystemExit(f"missing rebuilt basis: {basis_path}")
    frozen_basis_lines(basis_path)
    output_path = args.output
    if not output_path.is_absolute():
        output_path = ROOT / output_path
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required")

    expected_basis = str(basis_path.relative_to(ROOT))
    state = empty_state()
    if output_path.exists() and not args.no_resume:
        saved_checkpoint = json.loads(output_path.read_text())
        if saved_checkpoint.get("basis") != expected_basis:
            raise SystemExit(
                f"checkpoint basis mismatch in {output_path}"
            )
        try:
            state = load_state(saved_checkpoint)
        except ValueError as error:
            raise SystemExit(
                f"{error} in {output_path}"
            ) from error
        print(
            f"RESUMED_PRIMES={len(state['records'])} FROM={output_path}",
            flush=True,
        )
    completed_primes = {
        record["prime"]
        for record in state["records"]
    }
    pending_primes = []
    for prime in requested_primes:
        if prime in completed_primes:
            print(f"SKIP_PRIME={prime}", flush=True)
        else:
            pending_primes.append(prime)

    def checkpoint() -> None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(
                state_payload(state, expected_basis),
                separators=(",", ":"),
                sort_keys=True,
            )
            + "\n"
        )
        print(f"CHECKPOINT={output_path}", flush=True)

    def accept_record(record: dict) -> None:
        status = accept_modular_record(state, record)
        prime = record["prime"]
        print(
            f"PRIME={prime} "
            f"SECONDS={record['seconds']} "
            f"NONZERO_ENTRIES={record['nonzero_entries']} "
            f"TERMS={record['term_count']} "
            f"STATUS={status} "
            f"SUPPORT_SHA256={record['support_sha256']}",
            flush=True,
        )

    if args.jobs == 1:
        completed_since_checkpoint = 0
        for prime in pending_primes:
            print(f"START_PRIME={prime}", flush=True)
            accept_record(
                run_prime(
                    singular,
                    prime,
                    basis_path,
                    args.timeout,
                )
            )
            completed_since_checkpoint += 1
            if completed_since_checkpoint >= args.checkpoint_every:
                checkpoint()
                completed_since_checkpoint = 0
        if completed_since_checkpoint:
            checkpoint()
    else:
        print(
            f"PARALLEL_PRIMES={len(pending_primes)} JOBS={args.jobs}",
            flush=True,
        )
        completed_since_checkpoint = 0
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = [
                (
                    prime,
                    executor.submit(
                        run_prime,
                        singular,
                        prime,
                        basis_path,
                        args.timeout,
                    ),
                )
                for prime in pending_primes
            ]
            for prime, future in futures:
                try:
                    accept_record(future.result())
                except Exception as error:
                    print(f"FAILED_PRIME={prime}: {error}", flush=True)
                    for _, pending in futures:
                        pending.cancel()
                    if completed_since_checkpoint:
                        checkpoint()
                    raise
                completed_since_checkpoint += 1
                if completed_since_checkpoint >= args.checkpoint_every:
                    checkpoint()
                    completed_since_checkpoint = 0
        if completed_since_checkpoint:
            checkpoint()

    if not pending_primes and output_path.exists():
        # Normalize older pretty-printed checkpoints to compact JSON.
        checkpoint()

    support_counts = Counter(
        record["support_sha256"]
        for record in state["records"]
    )
    dominant_support = state["support_sha256"]
    good_records = [
        record
        for record in state["records"]
        if record["support_sha256"] == dominant_support
    ]
    unlucky_primes = [
        record["prime"]
        for record in state["records"]
        if record["support_sha256"] != dominant_support
    ]
    print(f"SUPPORT_PATTERNS={len(support_counts)}")
    print(f"DOMINANT_SUPPORT_SHA256={dominant_support}")
    print(f"RECONSTRUCTION_GOOD_PRIMES={len(good_records)}")
    print(
        "UNLUCKY_PRIMES="
        + (
            ",".join(str(prime) for prime in unlucky_primes)
            if unlucky_primes
            else "none"
        )
    )
    diagnostic = (
        None
        if args.skip_diagnostic
        else reconstruction_diagnostic(state)
    )
    if diagnostic is not None:
        print(
            "RECONSTRUCTION_DIAGNOSTIC="
            f"build_primes:{diagnostic['build_primes']},"
            f"modulus_bits:{diagnostic['build_modulus_bits']},"
            f"holdout:{diagnostic['holdout_prime']},"
            f"candidates:{diagnostic['candidates']},"
            f"validated:{diagnostic['validated']},"
            f"total:{diagnostic['total']},"
            f"max_num_bits:{diagnostic['maximum_numerator_bits']},"
            f"max_den_bits:{diagnostic['maximum_denominator_bits']},"
            "by_input:"
            + "|".join(
                f"{index}={count}"
                for index, count in diagnostic[
                    "validated_by_input"
                ].items()
            )
        )


if __name__ == "__main__":
    main()
