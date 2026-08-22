#!/usr/bin/env python3

from __future__ import annotations

import argparse
import ast
import glob
import json
from math import factorial
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from icarm_curve273 import short_coefficients
from curve273_full_ideal_chain import SUPPORTS, build_relations, prime_ideal


PROTOCOL = "R30RELPOOL"

S_RATIONAL = {
    2,
    3,
    5,
    7,
    13,
    31,
    41,
    47,
    53,
    67,
    379,
    4349,
    25721454817,
    97018222656318846556561979214040553412450110580812087282349817173780902099339117104673990259247421230916714670243202937,
}


INPUT_RE = re.compile(
    r"R30(?:SPECIALQ|MULTIQ)\|stage=input\|(.*)"
)

FACTOR_RE = re.compile(
    r"R30(?:SPECIALQ|MULTIQ)\|factorbest=\d+"
    r"\|delta_k=(-?\d+)"
    r"\|m=(-?\d+)"
    r".*?"
    r"\|cofactor_bits=(\d+)"
    r"\|largest_factor_bits=(\d+)"
    r"\|factor_count=(\d+)"
    r"\|factorization=([^|]+)"
)


CRT_BEST_RE = re.compile(
    r"R30CRT\|cyclebest=\d+\|(.*)"
)


def sage_q(QQ, ZZ, value):
    n = value.numerator
    d = value.denominator

    if callable(n):
        n = n()
    if callable(d):
        d = d()

    return QQ(ZZ(n)) / QQ(ZZ(d))


def parse_factorization(text):
    result = []

    for item in text.split("*"):
        item = item.strip()

        if "^" in item:
            p, e = item.split("^", 1)
            result.append((int(p), int(e)))
        else:
            result.append((int(item), 1))

    return tuple(result)


def parse_key_values(text):
    result = {}

    for part in text.split("|"):
        if "=" not in part:
            continue

        key, value = part.split("=", 1)
        result[key] = value

    return result


def parse_ideal_labels(text):
    if not text or text == "none":
        return ()

    return tuple(
        tuple(map(int, item.split(":", 1)))
        for item in text.split(",")
    )


def ideal_key(P):
    q = int(P.smallest_integer())

    try:
        hnf = str(P.pari_hnf())
    except Exception:
        hnf = str(P)

    return q, hnf


def insert_packed(pivots, row):
    v = int(row)

    while v:
        pivot = v.bit_length() - 1

        if pivot in pivots:
            v ^= pivots[pivot]
        else:
            pivots[pivot] = v
            return True

    return False


def xor_sets(left, right):
    return left.symmetric_difference(right)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--factor-base-bound",
        type=int,
        default=1_000_000,
    )

    ap.add_argument(
        "--max-largest-factor-bits",
        type=int,
        default=52,
    )

    ap.add_argument(
        "--max-relations-per-log",
        type=int,
        default=50,
    )

    ap.add_argument(
        "--glob",
        action="append",
        default=[],
    )

    ap.add_argument(
        "--include-full-ideal-chain",
        action="store_true",
        help=(
            "add the exact pinned arbitrary-ideal/CRT relations "
            "to the sparse pool"
        ),
    )

    ap.add_argument(
        "--include-crt-cycle-logs",
        action="store_true",
        help=(
            "add exact candidate relations from local R30CRT cycle logs"
        ),
    )

    ap.add_argument(
        "--max-crt-relations-per-log",
        type=int,
        default=25,
    )

    ap.add_argument(
        "--include-ideal-lattice-logs",
        action="store_true",
        help="add exact R30IDEAL best rows from bounded lattice-search logs",
    )

    ap.add_argument(
        "--ideal-glob",
        action="append",
        default=[],
        help="glob for R30IDEAL logs (repeatable)",
    )

    ap.add_argument(
        "--max-ideal-relations-per-log",
        type=int,
        default=20,
    )

    ap.add_argument(
        "--write-principal-relations",
        type=Path,
        help=(
            "write the exact principal generators and packed prime-ideal "
            "parity rows as a local BNF-free relation ledger"
        ),
    )

    ap.add_argument(
        "--complete-factor-base",
        action="store_true",
        help=(
            "materialize every prime ideal above rational primes through the "
            "factor-base bound, enabling a conditional Minkowski generation "
            "certificate when that bound is large enough"
        ),
    )
    ap.add_argument(
        "--write-large-prime-target-plan",
        type=Path,
        help=(
            "write a machine-readable next-target plan from the exact sparse "
            "large-prime incidence matrix"
        ),
    )
    ap.add_argument(
        "--target-plan-count",
        type=int,
        default=20,
        help="number of ranked unresolved large-prime ideals in a target plan",
    )

    args = ap.parse_args()
    if args.target_plan_count < 1:
        raise SystemExit("--target-plan-count must be positive")

    from sage.all import (
        QQ,
        RealField,
        ZZ,
        NumberField,
        PolynomialRing,
        prime_range,
    )

    patterns = args.glob or [
        "artifacts/local/elliptic-curves/specialq-*-factor.log",
        "artifacts/local/elliptic-curves/multiq-*.log",
    ]

    files = sorted(
        {
            Path(name)
            for pattern in patterns
            for name in glob.glob(pattern)
        }
    )

    crt_files = []

    if args.include_crt_cycle_logs:
        crt_files = sorted(
            Path(name)
            for name in glob.glob(
                "artifacts/local/elliptic-curves/crt-cycle-*.log"
            )
        )

    ideal_files = []

    if args.include_ideal_lattice_logs:
        ideal_patterns = args.ideal_glob or [
            "artifacts/local/elliptic-curves/curve273-*.log",
            "artifacts/local/elliptic-curves/ideal-lattice-*.log",
        ]
        ideal_files = sorted(
            {
                Path(name)
                for pattern in ideal_patterns
                for name in glob.glob(pattern)
            }
        )

    if (
        not files
        and not crt_files
        and not ideal_files
        and not args.include_full_ideal_chain
    ):
        raise SystemExit("no special-q logs found")

    print(
        f"{PROTOCOL}|stage=input"
        f"|files={len(files)}"
        f"|crt_files={len(crt_files)}"
        f"|ideal_files={len(ideal_files)}"
        f"|fb_bound={args.factor_base_bound}"
        f"|max_lp_bits={args.max_largest_factor_bits}",
        flush=True,
    )

    # --------------------------------------------------------
    # Curve / cubic field
    # --------------------------------------------------------

    coeffs = short_coefficients()

    A = ZZ(
        sage_q(QQ, ZZ, coeffs[3])
    )

    B = ZZ(
        sage_q(QQ, ZZ, coeffs[4])
    )

    R = PolynomialRing(QQ, "x")
    x = R.gen()

    f = x**3 + A*x + B

    K = NumberField(f, "theta")
    theta = K.gen()

    # --------------------------------------------------------
    # Ordinary rational factor base
    # --------------------------------------------------------

    trial_primes = tuple(
        int(p)
        for p in prime_range(
            2,
            args.factor_base_bound + 1,
        )
    )

    prime_decomposition_cache = {}

    def primes_above(q):
        q = int(q)

        if q not in prime_decomposition_cache:
            prime_decomposition_cache[q] = tuple(
                K.primes_above(ZZ(q))
            )

        return prime_decomposition_cache[q]

    # --------------------------------------------------------
    # Parse candidate m values from all logs.
    # --------------------------------------------------------

    records = []

    for path in files:
        lines = path.read_text().splitlines()

        special = None

        for line in lines:
            match = INPUT_RE.match(line)

            if not match:
                continue

            kv = parse_key_values(
                match.group(1)
            )

            if "special" in kv:
                special = tuple(
                    int(q)
                    for q in kv["special"].split(",")
                )

            elif "q" in kv:
                special = (
                    int(kv["q"]),
                )

            break

        if not special:
            print(
                f"{PROTOCOL}|stage=parse"
                f"|file={path}"
                f"|status=NO_SPECIAL_Q",
                flush=True,
            )
            continue

        accepted_here = 0

        for line in lines:
            match = FACTOR_RE.search(line)

            if not match:
                continue

            (
                delta_k,
                m,
                cofactor_bits,
                largest_bits,
                factor_count,
                factorization,
            ) = match.groups()

            delta_k = int(delta_k)

            # The seed replay contributes no new relation.
            if delta_k == 0:
                continue

            largest_bits = int(
                largest_bits
            )

            if (
                largest_bits
                > args.max_largest_factor_bits
            ):
                continue

            residual = parse_factorization(
                factorization
            )

            records.append(
                {
                    "file": str(path),
                    "special": special,
                    "m": int(m),
                    "delta_k": delta_k,
                    "cofactor_bits": int(cofactor_bits),
                    "largest_factor_bits": largest_bits,
                    "residual": residual,
                }
            )

            accepted_here += 1

            if (
                accepted_here
                >= args.max_relations_per_log
            ):
                break

        print(
            f"{PROTOCOL}|stage=parse"
            f"|file={path.name}"
            f"|special={','.join(map(str,special))}"
            f"|accepted={accepted_here}",
            flush=True,
        )

    # Deduplicate by m.
    by_m = {}

    for record in records:
        prior = by_m.get(
            record["m"]
        )

        if (
            prior is None
            or record["largest_factor_bits"]
            < prior["largest_factor_bits"]
        ):
            by_m[record["m"]] = record

    records = list(
        by_m.values()
    )

    records.sort(
        key=lambda r: (
            r["largest_factor_bits"],
            r["cofactor_bits"],
            abs(r["delta_k"]),
        )
    )

    print(
        f"{PROTOCOL}|stage=parse_summary"
        f"|unique_m={len(records)}",
        flush=True,
    )

    crt_records = []

    for path in crt_files:
        accepted_here = 0

        for line in path.read_text().splitlines():
            match = CRT_BEST_RE.match(line)

            if not match:
                continue

            kv = parse_key_values(match.group(1))

            if "m" not in kv or "candidate" not in kv:
                continue

            declared = parse_ideal_labels(kv["candidate"])

            if not declared:
                continue

            crt_records.append(
                {
                    "file": str(path),
                    "m": int(kv["m"]),
                    "declared": declared,
                }
            )

            accepted_here += 1

            if accepted_here >= args.max_crt_relations_per_log:
                break

        print(
            f"{PROTOCOL}|stage=parse_crt"
            f"|file={path.name}"
            f"|accepted={accepted_here}",
            flush=True,
        )

    print(
        f"{PROTOCOL}|stage=parse_crt_summary"
        f"|relations={len(crt_records)}",
        flush=True,
    )

    # --------------------------------------------------------
    # Sparse column registries.
    #
    # FB columns:
    #   q <= 1e6 OR q in S.
    #
    # LP columns:
    #   everything else.
    #
    # Columns identify PRIME IDEALS, never merely rational p.
    # --------------------------------------------------------

    fb_index = {}
    lp_index = {}

    fb_meta = {}
    fb_primes = {}
    lp_meta = {}

    lp_occurrences = defaultdict(list)

    def get_fb_column(P, q):
        key = ideal_key(P)

        if key not in fb_index:
            index = len(fb_index)
            fb_index[key] = index

            fb_meta[index] = {
                "q": int(q),
                "key": key,
                "S": int(q) in S_RATIONAL,
            }
            fb_primes[index] = P

        return fb_index[key]

    def get_lp_column(P, q, witness, residue):
        key = ideal_key(P)

        if key not in lp_index:
            index = len(lp_index)
            lp_index[key] = index

            lp_meta[index] = {
                "q": int(q),
                "key": key,
                "residue": int(residue),
            }

        index = lp_index[key]

        lp_occurrences[index].append(
            str(witness)
        )

        return index

    # A factor base that contains every prime ideal of norm at most the
    # Minkowski bound generates the ordinary ideal class group.  The existing
    # relation rows then give a one-sided (and therefore safe) mod-2 upper
    # bound after S-primes are killed: missing principal relations can only
    # make the reported quotient larger.  This is intentionally opt-in since
    # the relevant bound can be much larger than an efficient experimental
    # factor base.
    real_field = RealField(256)
    degree = K.degree()
    _, complex_places = K.signature()
    minkowski_raw = (
        (real_field(4) / real_field.pi()) ** complex_places
        * real_field(factorial(degree))
        / real_field(degree) ** degree
        * real_field(abs(K.discriminant())).sqrt()
    )
    minkowski_bound = ZZ(minkowski_raw.ceil())

    if args.complete_factor_base:
        for q in trial_primes:
            for P in primes_above(q):
                get_fb_column(P, q)
        for p in S_RATIONAL:
            for P in primes_above(p):
                get_fb_column(P, p)

    factor_base_generates = (
        args.complete_factor_base and args.factor_base_bound >= minkowski_bound
    )
    print(
        f"{PROTOCOL}|stage=minkowski_factor_base"
        f"|field_degree={degree}|complex_places={complex_places}"
        f"|bound={minkowski_bound}|fb_bound={args.factor_base_bound}"
        f"|materialized={args.complete_factor_base}"
        f"|generates_class_group={factor_base_generates}",
        flush=True,
    )

    # --------------------------------------------------------
    # Exact ideal rows.
    # --------------------------------------------------------

    exact_rows = []
    alpha_fingerprints = set()

    def alpha_fingerprint(alpha):
        return tuple(
            (
                int(QQ(c).numerator()),
                int(QQ(c).denominator()),
            )
            for c in alpha.list()
        )

    def principal_generator(alpha):
        """Stable power-basis display for an exact principal generator.

        The sparse LP matrix is only an incidence device.  Keeping this
        expression with the row means an LP dependency can be replayed as a
        product of actual principal generators, rather than remaining an
        anonymous valuation-vector dependency.
        """

        return "(" + ",".join(str(QQ(c)) for c in alpha.list()) + ")"

    t0 = time.monotonic()

    for row_number, record in enumerate(
        records,
        1,
    ):
        m = ZZ(record["m"])

        N = abs(
            ZZ(f(m))
        )

        if N == 0:
            continue

        co = N

        rational_support = []

        # Small rational factor base.
        for p in trial_primes:
            if co % p:
                continue

            exponent = 0

            while co % p == 0:
                co //= p
                exponent += 1

            rational_support.append(
                (p, exponent)
            )

            if co == 1:
                break

        # Strip every S-prime outside the ordinary bound.
        for p0 in S_RATIONAL:
            p = ZZ(p0)

            if p <= args.factor_base_bound:
                continue

            if co % p:
                continue

            exponent = 0

            while co % p == 0:
                co //= p
                exponent += 1

            rational_support.append(
                (int(p), exponent)
            )

        # The remaining part must be explained completely by:
        #
        #   forced special q's
        #   + logged residual factorization.
        expected = ZZ(1)

        for q in record["special"]:
            expected *= ZZ(q)

        for p, exponent in record["residual"]:
            pz = ZZ(p)

            if not pz.is_prime(
                proof=True
            ):
                raise RuntimeError(
                    f"logged factor {p} is not prime"
                )

            expected *= (
                pz ** exponent
            )

        if co != expected:
            raise RuntimeError(
                "log reconstruction mismatch:\n"
                f"m={m}\n"
                f"remaining={co}\n"
                f"expected={expected}\n"
                f"file={record['file']}"
            )

        # Add the >FB support.
        for q in record["special"]:
            rational_support.append(
                (int(q), 1)
            )

        rational_support.extend(
            record["residual"]
        )

        alpha = K(m) - theta
        fingerprint = alpha_fingerprint(alpha)

        if fingerprint in alpha_fingerprints:
            continue

        alpha_fingerprints.add(fingerprint)

        fb_row = 0
        lp_row = 0

        ideal_support = []

        # Exact prime-ideal valuations.
        for q, rational_exponent in rational_support:

            for P in primes_above(q):
                valuation = int(
                    alpha.valuation(P)
                )

                if valuation == 0:
                    continue

                ideal_support.append(
                    (
                        int(q),
                        ideal_key(P),
                        valuation,
                    )
                )

                if valuation & 1:

                    if (
                        q <= args.factor_base_bound
                        or q in S_RATIONAL
                    ):
                        c = get_fb_column(
                            P,
                            q,
                        )

                        fb_row ^= (
                            1 << c
                        )

                    else:
                        c = get_lp_column(
                            P,
                            q,
                            int(m),
                            int(m % q),
                        )

                        lp_row ^= (
                            1 << c
                        )

        # Norm parity sanity.
        reconstructed_parity = defaultdict(int)

        for q, key, valuation in ideal_support:
            # Norm(P)=q^f.  For every q in f(m), relevant
            # support here should reproduce rational valuation.
            P = next(
                P
                for P in primes_above(q)
                if ideal_key(P) == key
            )

            degree = int(
                P.residue_class_degree()
            )

            reconstructed_parity[q] += (
                valuation * degree
            )

        for q, exponent in rational_support:
            if (
                reconstructed_parity[q]
                != exponent
            ):
                raise RuntimeError(
                    f"norm valuation mismatch q={q} "
                    f"got={reconstructed_parity[q]} "
                    f"expected={exponent}"
                )

        exact_rows.append(
            {
                "m": int(m),
                "fb": fb_row,
                "lp": lp_row,
                "file": record["file"],
                "generator": principal_generator(alpha),
            }
        )

        if (
            row_number % 10 == 0
            or row_number == len(records)
        ):
            print(
                f"{PROTOCOL}|stage=exact_rows"
                f"|processed={row_number}/{len(records)}"
                f"|fb_columns={len(fb_index)}"
                f"|lp_columns={len(lp_index)}"
                f"|seconds={time.monotonic()-t0:.3f}",
                flush=True,
            )

    def append_declared_relation(alpha, declared, witness, source):
        """Verify and append a relation with all non-FB ideals declared."""

        fingerprint = alpha_fingerprint(alpha)

        if fingerprint in alpha_fingerprints:
            return False

        if not alpha.is_integral():
            raise RuntimeError(f"nonintegral declared relation {witness}")

        norm = abs(ZZ(alpha.norm()))
        co = norm
        rational_support = []

        for p in trial_primes:
            if co % p:
                continue

            exponent = 0

            while co % p == 0:
                co //= p
                exponent += 1

            rational_support.append((p, exponent))

            if co == 1:
                break

        for p0 in S_RATIONAL:
            p = ZZ(p0)

            if p <= args.factor_base_bound or co % p:
                continue

            exponent = 0

            while co % p == 0:
                co //= p
                exponent += 1

            rational_support.append((int(p), exponent))

        expected = ZZ(1)

        for q, residue in declared:
            qz = ZZ(q)

            if not qz.is_prime(proof=True):
                raise RuntimeError(f"declared factor {q} is not prime")

            if q > args.factor_base_bound and q not in S_RATIONAL:
                expected *= qz

        if co != expected:
            raise RuntimeError(
                "declared relation reconstruction mismatch:\n"
                f"witness={witness}\nremaining={co}\nexpected={expected}"
            )

        rational_support.extend(
            (int(q), 1)
            for q, residue in declared
            if q > args.factor_base_bound and q not in S_RATIONAL
        )

        # Several declared degree-one ideals may lie above the same rational
        # prime (notably the conjugate-pair replacement rows).  The norm audit
        # is rational-prime based, so aggregate those exponents before
        # traversing every prime ideal above q.
        aggregated_support = defaultdict(int)
        for q, exponent in rational_support:
            aggregated_support[int(q)] += int(exponent)
        rational_support = sorted(aggregated_support.items())

        declared_residues = {}

        for label in declared:
            P = prime_ideal(K, theta, label)

            if not P.is_prime() or P.norm() != label[0]:
                raise RuntimeError(f"bad declared prime ideal {label}")

            if int(alpha.valuation(P)) != 1:
                raise RuntimeError(
                    f"declared valuation is not one: {witness} {label}"
                )

            declared_residues[ideal_key(P)] = int(label[1])

        fb_row = 0
        lp_row = 0
        reconstructed = defaultdict(int)

        for q, rational_exponent in rational_support:
            for P in primes_above(q):
                valuation = int(alpha.valuation(P))

                if valuation == 0:
                    continue

                key = ideal_key(P)
                reconstructed[q] += valuation * int(
                    P.residue_class_degree()
                )

                if (valuation & 1) == 0:
                    continue

                if q <= args.factor_base_bound or q in S_RATIONAL:
                    column = get_fb_column(P, q)
                    fb_row ^= 1 << column
                else:
                    if key not in declared_residues:
                        raise RuntimeError(
                            "undeclared large-prime ideal "
                            f"witness={witness} q={q}"
                        )

                    column = get_lp_column(
                        P,
                        q,
                        witness,
                        declared_residues[key],
                    )
                    lp_row ^= 1 << column

        for q, exponent in rational_support:
            if reconstructed[q] != exponent:
                raise RuntimeError(
                    "declared norm valuation mismatch "
                    f"witness={witness} q={q} "
                    f"got={reconstructed[q]} expected={exponent}"
                )

        alpha_fingerprints.add(fingerprint)
        exact_rows.append(
            {
                "m": str(witness),
                "fb": fb_row,
                "lp": lp_row,
                "file": str(source),
                "generator": principal_generator(alpha),
            }
        )
        return True

    if args.include_full_ideal_chain:
        chain_relations = build_relations(K, theta)
        added = 0

        for relation in chain_relations:
            added += append_declared_relation(
                relation["alpha"],
                relation["declared"],
                relation["name"],
                "full-ideal-chain",
            )

        print(
            f"{PROTOCOL}|stage=chain_rows"
            f"|added={added}"
            f"|duplicates={len(chain_relations)-added}"
            f"|fb_columns={len(fb_index)}"
            f"|lp_columns={len(lp_index)}",
            flush=True,
        )

    if args.include_crt_cycle_logs:
        added = 0

        for record in crt_records:
            added += append_declared_relation(
                K(ZZ(record["m"])) - theta,
                record["declared"],
                record["m"],
                record["file"],
            )

        print(
            f"{PROTOCOL}|stage=crt_rows"
            f"|added={added}"
            f"|duplicates={len(crt_records)-added}"
            f"|fb_columns={len(fb_index)}"
            f"|lp_columns={len(lp_index)}",
            flush=True,
        )

    if args.include_ideal_lattice_logs:
        added = 0
        duplicates = 0

        for path in ideal_files:
            lines = path.read_text().splitlines()
            input_line = next(
                (line for line in lines if line.startswith("R30IDEAL|stage=input|")),
                None,
            )
            if input_line is None:
                print(
                    f"{PROTOCOL}|stage=parse_ideal|file={path.name}|status=NO_INPUT",
                    flush=True,
                )
                continue

            input_kv = parse_key_values(input_line)
            targets = parse_ideal_labels(input_kv.get("targets", ""))
            if not targets:
                raise RuntimeError(f"missing ideal targets in {path}")

            accepted_here = 0
            for line_number, line in enumerate(lines, 1):
                if not line.startswith("R30IDEAL|best|"):
                    continue

                kv = parse_key_values(line)
                if not {"twist", "coordinates", "labels"} <= kv.keys():
                    continue

                label_records = ast.literal_eval(kv["labels"])
                residual = []
                valid = True
                for q, roots, exponent in label_records:
                    if int(exponent) != 1 or len(roots) != 1:
                        valid = False
                        break
                    residual.append((int(q), int(roots[0])))
                if not valid:
                    continue

                twist_labels = ()
                if kv["twist"] != "none":
                    twist_labels = parse_ideal_labels(
                        kv["twist"].replace("+", ",")
                    )

                source_ideal = K.ideal(1)
                for label in targets + twist_labels:
                    source_ideal *= prime_ideal(K, theta, label)
                basis = tuple(source_ideal.basis())
                coordinates = tuple(
                    ZZ(value) for value in kv["coordinates"].split(",")
                )
                if len(basis) != 3 or len(coordinates) != 3:
                    raise RuntimeError(f"bad ideal row dimensions in {path}:{line_number}")
                alpha = sum(
                    (coordinates[index] * basis[index] for index in range(3)),
                    K(0),
                )
                witness = f"{path.name}:{line_number}:{kv['twist']}"
                was_added = append_declared_relation(
                    alpha,
                    targets + tuple(residual),
                    witness,
                    path,
                )
                added += int(was_added)
                duplicates += int(not was_added)
                accepted_here += 1
                if accepted_here >= args.max_ideal_relations_per_log:
                    break

            print(
                f"{PROTOCOL}|stage=parse_ideal|file={path.name}"
                f"|accepted={accepted_here}",
                flush=True,
            )

        print(
            f"{PROTOCOL}|stage=ideal_rows|added={added}|duplicates={duplicates}"
            f"|fb_columns={len(fb_index)}|lp_columns={len(lp_index)}",
            flush=True,
        )

    if args.write_principal_relations:
        # This is an exact, local checkpoint for the subsequent quotient and
        # certification lane.  It deliberately stores generators and prime
        # ideals, rather than only a rank summary, so every later dependency
        # can be replayed in the cubic field.
        ledger = {
            "schema": "elliptic-curves.bnf-free-principal-relation-ledger.v1",
            "status": "exact_relation_rows_not_class_group_completion",
            "factor_base_bound": args.factor_base_bound,
            "fb_columns": [
                {
                    "index": index,
                    "q": meta["q"],
                    "S": meta["S"],
                    "prime_ideal_key": list(meta["key"]),
                }
                for index, meta in sorted(fb_meta.items())
            ],
            "lp_columns": [
                {
                    "index": index,
                    "q": meta["q"],
                    "residue": meta["residue"],
                    "prime_ideal_key": list(meta["key"]),
                }
                for index, meta in sorted(lp_meta.items())
            ],
            "relations": [
                {
                    "label": relation["m"],
                    "generator_power_basis": relation["generator"],
                    "fb_parity_mask_hex": hex(relation["fb"]),
                    "lp_parity_mask_hex": hex(relation["lp"]),
                    "source": relation["file"],
                }
                for relation in exact_rows
            ],
        }
        args.write_principal_relations.parent.mkdir(parents=True, exist_ok=True)
        args.write_principal_relations.write_text(
            json.dumps(ledger, indent=2, sort_keys=True) + "\n"
        )
        print(
            f"{PROTOCOL}|stage=write_principal_relations"
            f"|path={args.write_principal_relations}|relations={len(exact_rows)}"
            f"|status=EXACT_ROWS_NOT_COMPLETENESS_CERTIFICATE",
            flush=True,
        )

    # --------------------------------------------------------
    # Gaussian elimination on LARGE-PRIME columns only.
    #
    # Every dependency here produces an exact principal
    # relation supported entirely in the ordinary FB + S.
    # --------------------------------------------------------

    pivots = {}

    pure_fb_rows = []

    trivial_cycles = 0

    for relation_index, relation in enumerate(
        exact_rows
    ):
        lp = relation["lp"]
        fb = relation["fb"]

        provenance = {
            relation_index
        }

        while lp:
            pivot = lp.bit_length() - 1

            if pivot not in pivots:
                pivots[pivot] = (
                    lp,
                    fb,
                    provenance,
                )
                break

            (
                pivot_lp,
                pivot_fb,
                pivot_provenance,
            ) = pivots[pivot]

            lp ^= pivot_lp
            fb ^= pivot_fb

            provenance = xor_sets(
                provenance,
                pivot_provenance,
            )

        if lp == 0:

            if fb:
                pure_fb_rows.append(
                    (
                        fb,
                        provenance,
                    )
                )
            else:
                trivial_cycles += 1

    # Independent rank among obtained pure-FB rows.
    fb_pivots = {}

    independent_pure = []

    for fb, provenance in pure_fb_rows:

        if insert_packed(
            fb_pivots,
            fb,
        ):
            independent_pure.append(
                (
                    fb,
                    provenance,
                )
            )

    lp_rank = len(pivots)

    lp_nullity = (
        len(exact_rows)
        - lp_rank
    )

    # Quotient the exact pure-FB relations by all Selmer-support prime ideals.
    # When ``factor_base_generates`` holds this is an unconditional upper bound
    # on Cl(O_K[S^{-1}])[2], not a rank-stabilization heuristic.
    s_class_pivots = dict(fb_pivots)
    for index, meta in fb_meta.items():
        if meta["S"]:
            insert_packed(s_class_pivots, 1 << index)
    s_class_dimension_upper_bound = len(fb_index) - len(s_class_pivots)
    class_bound_status = (
        "CERTIFIED_MINKOWSKI_FACTOR_BASE"
        if factor_base_generates
        else "UNCERTIFIED_FACTOR_BASE_DOES_NOT_REACH_MINKOWSKI_BOUND"
    )

    print(
        f"{PROTOCOL}|stage=elimination"
        f"|relations={len(exact_rows)}"
        f"|lp_columns={len(lp_index)}"
        f"|lp_rank={lp_rank}"
        f"|lp_nullity={lp_nullity}"
        f"|pure_fb_relations={len(pure_fb_rows)}"
        f"|independent_pure_fb={len(independent_pure)}"
        f"|trivial_cycles={trivial_cycles}",
        flush=True,
    )
    dimension_field = (
        f"dimension_upper_bound={s_class_dimension_upper_bound}"
        if factor_base_generates
        else f"dimension_factor_base_model={s_class_dimension_upper_bound}"
    )
    print(
        f"{PROTOCOL}|stage=s_class_mod2_upper_bound"
        f"|status={class_bound_status}"
        f"|{dimension_field}"
        f"|relation_rank={len(fb_pivots)}"
        f"|S_columns={sum(meta['S'] for meta in fb_meta.values())}"
        f"|interpretation=IDEAL_CLASS_QUOTIENT_ONLY_NOT_SELMER",
        flush=True,
    )

    if args.write_principal_relations:
        # Extend the early exact-row checkpoint with the LP-eliminated rows
        # that the generic BNF-free S-class audit consumes.  The old fields
        # remain for backwards-compatible chain diagnostics; these standard
        # fields make each closed row replayable as a product of its stored
        # principal generators.
        ledger = json.loads(args.write_principal_relations.read_text())
        ledger.update(
            {
                "field_polynomial": str(f),
                "defining_polynomial_ascending": [str(B), str(A), "0", "1"],
                "field_discriminant": str(K.discriminant()),
                "generator_coordinate_order": ["1", "theta", "theta^2"],
                "factor_base_completion": {
                    "all_prime_ideals_above_rational_primes_through": (
                        args.factor_base_bound
                    ),
                    "materialized_complete_factor_base": bool(
                        args.complete_factor_base
                    ),
                    "extra_declared_S_rational_primes": sorted(S_RATIONAL),
                },
                "selmer_rational_primes": sorted(S_RATIONAL),
                "factor_base": [
                    {
                        "hnf": str(fb_primes[index].pari_hnf()),
                        "norm": int(fb_primes[index].norm()),
                        "residue_degree": int(
                            fb_primes[index].residue_class_degree()
                        ),
                        "rational_prime": int(fb_primes[index].smallest_integer()),
                    }
                    for index in sorted(fb_meta)
                ],
                "S_columns": [
                    index for index, meta in sorted(fb_meta.items()) if meta["S"]
                ],
                "generators": [
                    {
                        "power_basis": relation["generator"].strip("()").split(
                            ","
                        ),
                        "source": relation["file"],
                        "label": str(relation["m"]),
                    }
                    for relation in exact_rows
                ],
                "closed_relations": [
                    {
                        "fb_parity_mask_hex": hex(fb),
                        "generator_indices": sorted(provenance),
                        "kind": "lp_eliminated_relation",
                    }
                    for fb, provenance in pure_fb_rows
                ],
            }
        )
        args.write_principal_relations.write_text(
            json.dumps(ledger, indent=2, sort_keys=True) + "\n"
        )
        print(
            f"{PROTOCOL}|stage=write_bnf_free_closures"
            f"|path={args.write_principal_relations}"
            f"|generators={len(exact_rows)}|closed_relations={len(pure_fb_rows)}",
            flush=True,
        )

    # Solve the exact sparse system against the certified chain endpoint.
    # This proves whether the logged pool genuinely connects the original
    # special-q component to that four-ideal residual, without claiming an
    # LP-free cycle.
    if args.include_full_ideal_chain:
        target_mask = 0
        missing = []

        for label in SUPPORTS[-1]:
            P = prime_ideal(K, theta, label)
            key = ideal_key(P)

            if key not in lp_index:
                missing.append(label)
                continue

            target_mask ^= 1 << lp_index[key]

        target_lp = target_mask
        target_fb = 0
        target_provenance = set()

        while target_lp:
            pivot = target_lp.bit_length() - 1

            if pivot not in pivots:
                break

            pivot_lp, pivot_fb, pivot_provenance = pivots[pivot]
            target_lp ^= pivot_lp
            target_fb ^= pivot_fb
            target_provenance = xor_sets(
                target_provenance,
                pivot_provenance,
            )

        if missing:
            target_status = "MISSING_COLUMNS"
        elif target_lp:
            target_status = "NOT_IN_SPAN"
        else:
            target_status = "IN_SPAN"

        target_witnesses = ",".join(
            str(exact_rows[index]["m"])
            for index in sorted(target_provenance)
        )

        print(
            f"{PROTOCOL}|stage=chain_endpoint_span"
            f"|status={target_status}"
            f"|interpretation=REACHABILITY_NOT_LP_CLOSURE"
            f"|target_lp={target_mask.bit_count()}"
            f"|unresolved_lp={target_lp.bit_count()}"
            f"|combined_rows={len(target_provenance)}"
            f"|fb_support={target_fb.bit_count()}"
            f"|witnesses={target_witnesses}",
            flush=True,
        )

    # --------------------------------------------------------
    # Report exact FB-only relations.
    # --------------------------------------------------------

    for index, (
        fb,
        provenance,
    ) in enumerate(
        independent_pure[:20],
        1,
    ):
        support = []

        for column, meta in fb_meta.items():
            if (
                fb >> column
            ) & 1:
                support.append(
                    (
                        meta["q"],
                        int(meta["S"]),
                    )
                )

        ms = [
            exact_rows[i]["m"]
            for i in sorted(provenance)
        ]
        generators = [
            exact_rows[i]["generator"]
            for i in sorted(provenance)
        ]

        print(
            f"{PROTOCOL}|pure={index}"
            f"|combined_rows={len(provenance)}"
            f"|fb_support={len(support)}"
            f"|S_support={sum(flag for q,flag in support)}"
            f"|m={','.join(map(str,ms))}"
            + "|generator_product="
            + "*".join(generators)
            + f"|rational_support="
            + ",".join(
                (
                    f"{q}{'S' if is_s else ''}"
                    for q, is_s in support
                )
            ),
            flush=True,
        )

    # --------------------------------------------------------
    # Pairwise cycle engineering.
    #
    # XOR two exact principal-ideal relations.  Any shared LP
    # columns cancel.  Rank pairs by the size/height of the LP
    # support that remains.
    #
    # This does NOT create new mathematical information; it
    # identifies the cheapest partial relations we should try
    # to close with CRT special-q sieving.
    # --------------------------------------------------------

    pair_candidates = []

    for left in range(len(exact_rows)):
        L = exact_rows[left]["lp"]

        for right in range(left + 1, len(exact_rows)):
            R = exact_rows[right]["lp"]

            shared = L & R

            shared_count = shared.bit_count()

            # One shared prime is usually just a star edge.
            # Prefer pairs cancelling at least two LP ideals.
            if shared_count < 2:
                continue

            remaining = L ^ R

            remaining_columns = [
                column
                for column in range(len(lp_meta))
                if (remaining >> column) & 1
            ]

            remaining_count = len(
                remaining_columns
            )

            remaining_info = []

            for column in remaining_columns:
                meta = lp_meta[column]

                q = int(meta["q"])
                residue = int(meta["residue"])

                remaining_info.append(
                    (
                        q.bit_length(),
                        q,
                        residue,
                        column,
                    )
                )

            remaining_info.sort()

            max_bits = (
                max(
                    bits
                    for bits, q, residue, column
                    in remaining_info
                )
                if remaining_info
                else 0
            )

            total_bits = sum(
                bits
                for bits, q, residue, column
                in remaining_info
            )

            # Also track shared ideals for diagnostics.
            shared_info = []

            for column in range(len(lp_meta)):
                if not (
                    (shared >> column) & 1
                ):
                    continue

                meta = lp_meta[column]

                q = int(meta["q"])
                shared_info.append(
                    (
                        q,
                        int(meta["residue"]),
                    )
                )

            score = (
                max_bits,
                total_bits,
                remaining_count,
                -shared_count,
                left,
                right,
            )

            pair_candidates.append(
                (
                    score,
                    left,
                    right,
                    shared_count,
                    shared_info,
                    remaining_info,
                )
            )

    pair_candidates.sort(
        key=lambda item: item[0]
    )

    print(
        f"{PROTOCOL}|stage=pair_search"
        f"|candidate_pairs={len(pair_candidates)}",
        flush=True,
    )

    for rank, item in enumerate(
        pair_candidates[:30],
        1,
    ):
        (
            score,
            left,
            right,
            shared_count,
            shared_info,
            remaining_info,
        ) = item

        remaining_count = len(
            remaining_info
        )

        max_bits = (
            max(
                bits
                for bits, q, residue, column
                in remaining_info
            )
            if remaining_info
            else 0
        )

        shared_text = ",".join(
            f"{q}:{residue}"
            for q, residue
            in shared_info
        )

        remaining_text = ",".join(
            f"{q}:{residue}"
            for bits, q, residue, column
            in remaining_info
        )

        print(
            f"{PROTOCOL}|pairbest={rank}"
            f"|row1={left}"
            f"|row2={right}"
            f"|m1={exact_rows[left]['m']}"
            f"|m2={exact_rows[right]['m']}"
            f"|shared={shared_count}"
            f"|remaining={remaining_count}"
            f"|max_bits={max_bits}"
            f"|shared_ideals={shared_text}"
            f"|remaining_ideals={remaining_text}",
            flush=True,
        )

    # --------------------------------------------------------
    # Triple combinations around the strongest hub.
    #
    # Only test triples when the union shares at least two
    # columns pairwise.  With 115 rows this is still tiny.
    # --------------------------------------------------------

    triple_candidates = []

    hub_rows = defaultdict(list)

    for row_index, relation in enumerate(
        exact_rows
    ):
        mask = relation["lp"]

        for column in range(len(lp_meta)):
            if (mask >> column) & 1:
                hub_rows[column].append(
                    row_index
                )

    # Only hubs occurring frequently enough to matter.
    hubs = sorted(
        (
            (len(rows), column, rows)
            for column, rows
            in hub_rows.items()
            if len(rows) >= 3
        ),
        reverse=True,
    )

    for incidence, hub, rows in hubs[:10]:

        # Limit combinatorics to the strongest rows at a hub.
        rows = rows[:40]

        for ai in range(len(rows)):
            i = rows[ai]

            for aj in range(ai + 1, len(rows)):
                j = rows[aj]

                for ak in range(aj + 1, len(rows)):
                    k = rows[ak]

                    remaining = (
                        exact_rows[i]["lp"]
                        ^ exact_rows[j]["lp"]
                        ^ exact_rows[k]["lp"]
                    )

                    columns = [
                        column
                        for column in range(len(lp_meta))
                        if (remaining >> column) & 1
                    ]

                    info = []

                    for column in columns:
                        meta = lp_meta[column]

                        q = int(meta["q"])
                        info.append(
                            (
                                q.bit_length(),
                                q,
                                int(meta["residue"]),
                                column,
                            )
                        )

                    info.sort()

                    max_bits = (
                        max(
                            x[0]
                            for x in info
                        )
                        if info
                        else 0
                    )

                    total_bits = sum(
                        x[0]
                        for x in info
                    )

                    triple_candidates.append(
                        (
                        (
                            max_bits,
                            total_bits,
                            len(info),
                            i,
                            j,
                            k,
                            ),
                            hub,
                            i,
                            j,
                            k,
                            info,
                        )
                    )

    triple_candidates.sort(
        key=lambda item: item[0]
    )

    print(
        f"{PROTOCOL}|stage=triple_search"
        f"|candidate_triples={len(triple_candidates)}",
        flush=True,
    )

    for rank, item in enumerate(
        triple_candidates[:20],
        1,
    ):
        score, hub, i, j, k, info = item

        max_bits = (
            max(x[0] for x in info)
            if info
            else 0
        )

        remaining_text = ",".join(
            f"{q}:{residue}"
            for bits, q, residue, column
            in info
        )

        hub_q = int(
            lp_meta[hub]["q"]
        )

        print(
            f"{PROTOCOL}|triplebest={rank}"
            f"|hub_q={hub_q}"
            f"|rows={i},{j},{k}"
            f"|m={exact_rows[i]['m']},"
            f"{exact_rows[j]['m']},"
            f"{exact_rows[k]['m']}"
            f"|remaining={len(info)}"
            f"|max_bits={max_bits}"
            f"|remaining_ideals={remaining_text}",
            flush=True,
        )

    # --------------------------------------------------------
    # Which unresolved LP ideals should get another
    # special-q search?
    #
    # Prefer:
    #  1. actual pivot columns;
    #  2. smaller rational q;
    #  3. columns already seen multiple times.
    # --------------------------------------------------------

    pivot_columns = set(
        pivots
    )

    candidates = []

    for column, meta in lp_meta.items():

        occurrences = list(
            dict.fromkeys(
                lp_occurrences[column]
            )
        )

        candidates.append(
            (
                0 if column in pivot_columns else 1,
                -len(occurrences),
                int(meta["q"]).bit_length(),
                int(meta["q"]),
                column,
                occurrences,
            )
        )

    candidates.sort()

    target_plan = []
    for rank, item in enumerate(
        candidates[: args.target_plan_count],
        1,
    ):
        (
            pivot_penalty,
            neg_count,
            qbits,
            q,
            column,
            occurrences,
        ) = item

        print(
            f"{PROTOCOL}|target={rank}"
            f"|q={q}"
            f"|q_bits={qbits}"
            f"|root={lp_meta[column]['residue']}"
            f"|incidence={-neg_count}"
            f"|pivot={int(column in pivot_columns)}"
            f"|witness={occurrences[0]}"
            f"|ideal_key={lp_meta[column]['key'][1]}",
            flush=True,
        )
        target_plan.append(
            {
                "rank": rank,
                "rational_prime": int(q),
                "residue": int(lp_meta[column]["residue"]),
                "prime_ideal_hnf": lp_meta[column]["key"][1],
                "occurrence_count": -neg_count,
                "is_elimination_pivot": bool(column in pivot_columns),
                # Some exact relation sources use numeric rows while the
                # ideal-lattice source uses stable symbolic labels (for
                # example ``I9``).  A target plan is an audit trail, not a
                # row-number API, so preserve the source identifier exactly.
                "first_witness_relation_id": str(occurrences[0]),
            }
        )

    if args.write_large_prime_target_plan:
        plan = {
            "schema": "elliptic-curves.bnf-free-large-prime-target-plan.v1",
            "status": "EXACT_SPARSE_INCIDENT_TARGETS_NOT_A_SELMER_CERTIFICATE",
            "factor_base_bound": args.factor_base_bound,
            "exact_row_count": len(exact_rows),
            "large_prime_column_count": len(lp_index),
            "large_prime_rank": lp_rank,
            "large_prime_nullity": lp_nullity,
            "targets": target_plan,
            "reproduction": {
                "search_script": "elliptic-curves/cas/search_curve273_ideal_lattice_relations.sage",
                "target_syntax": "rational_prime:residue",
                "selection": (
                    "elimination pivots first, then higher incidence, then "
                    "smaller rational prime"
                ),
            },
        }
        args.write_large_prime_target_plan.parent.mkdir(parents=True, exist_ok=True)
        args.write_large_prime_target_plan.write_text(
            json.dumps(plan, indent=2, sort_keys=True) + "\n"
        )
        print(
            f"{PROTOCOL}|stage=write_target_plan"
            f"|path={args.write_large_prime_target_plan}|targets={len(target_plan)}"
            "|status=EXACT_SPARSE_INCIDENT_TARGETS",
            flush=True,
        )

    print(
        f"{PROTOCOL}|stage=done"
        f"|exact_rows={len(exact_rows)}"
        f"|pure_fb={len(independent_pure)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
