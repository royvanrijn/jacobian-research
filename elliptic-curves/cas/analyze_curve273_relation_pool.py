#!/usr/bin/env python3

from __future__ import annotations

import argparse
import glob
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from icarm_curve273 import short_coefficients


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

    args = ap.parse_args()

    from sage.all import (
        QQ,
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

    if not files:
        raise SystemExit("no special-q logs found")

    print(
        f"{PROTOCOL}|stage=input"
        f"|files={len(files)}"
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

        return fb_index[key]

    def get_lp_column(P, q, m):
        key = ideal_key(P)

        if key not in lp_index:
            index = len(lp_index)
            lp_index[key] = index

            lp_meta[index] = {
                "q": int(q),
                "key": key,
                "first_m": int(m),
            }

        index = lp_index[key]

        lp_occurrences[index].append(
            int(m)
        )

        return index

    # --------------------------------------------------------
    # Exact ideal rows.
    # --------------------------------------------------------

    exact_rows = []

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
                            m,
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

        print(
            f"{PROTOCOL}|pure={index}"
            f"|combined_rows={len(provenance)}"
            f"|fb_support={len(support)}"
            f"|S_support={sum(flag for q,flag in support)}"
            f"|m={','.join(map(str,ms))}"
            f"|rational_support="
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
                witness_m = int(
                    meta["first_m"]
                )

                residue = (
                    witness_m % q
                )

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
                witness_m = int(
                    meta["first_m"]
                )

                shared_info.append(
                    (
                        q,
                        witness_m % q,
                    )
                )

            score = (
                remaining_count,
                max_bits,
                total_bits,
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
                        witness_m = int(
                            meta["first_m"]
                        )

                        info.append(
                            (
                                q.bit_length(),
                                q,
                                witness_m % q,
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
                                len(info),
                                max_bits,
                                total_bits,
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

    for rank, item in enumerate(
        candidates[:20],
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
            f"|incidence={-neg_count}"
            f"|pivot={int(column in pivot_columns)}"
            f"|seed_m={occurrences[0]}"
            f"|ideal_key={lp_meta[column]['key'][1]}",
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
