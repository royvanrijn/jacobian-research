#!/usr/bin/env sage -python
"""Simultaneously reconstruct the q4/orbit230 projective section.

status: CONSTRUCTION EXPERIMENT, promoted only by literal QQ substitution

Use small overlapping coefficient windows rather than one large LLL lattice.
Each short lattice vector supplies a candidate common denominator; centered
residues extend it to the whole X or Y polynomial.  Candidate pairs are first
filtered at auxiliary good primes and only survivors are tested over QQ.  No
Groebner basis is used.
"""

import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
SURFACE = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--hensel-source", required=True)
parser.add_argument("--checkpoint-index", type=int, default=-1)
parser.add_argument("--combination-bound", type=int, default=3)
parser.add_argument("--short-rank", type=int, default=3)
parser.add_argument(
    "--output",
    default="artifacts/local/elkies-k3/q24-2a5-p230-simultaneous-reconstruction.json",
)
args = parser.parse_args()
assert args.combination_bound >= 1 and 1 <= args.short_rank <= 4


def resolved_path(value):
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


started = time.monotonic()
source_path = resolved_path(args.hensel_source)
output = resolved_path(args.output)
source = json.loads(source_path.read_text())
surface = json.loads(SURFACE.read_text())
assert source["section"] == "p230"
assert tuple(source["projective_chart"]["degrees_X_Y_Z"]) == (8, 12, 2)
fixed_Z = source["projective_chart"]["fixed_Z_QQ_candidate"]
assert fixed_Z is not None
checkpoint = source["hensel"]["checkpoints"][args.checkpoint_index]
p = ZZ(source["prime"])
digits = int(checkpoint["p_adic_digits"])
modulus = p**digits
residues = [ZZ(value) for value in checkpoint["residues_0_to_p_power_minus_1"]]
assert len(residues) == 22


def balanced(value):
    value = ZZ(value) % modulus
    return value - modulus if value > modulus // 2 else value


def window_lattice(block_residues, indices):
    selected = [block_residues[index] for index in indices]
    size = len(selected)
    basis = matrix(ZZ, size + 1, size + 1)
    basis.set_block(0, 0, modulus * matrix.identity(ZZ, size))
    basis[size] = selected + [1]
    return basis.LLL(delta=0.99)


def extended_candidate(block_residues, denominator):
    if not denominator:
        return None
    values = tuple(QQ(balanced(denominator * residue)) / QQ(denominator)
                   for residue in block_residues)
    max_bits = max(
        max(abs(value.numerator()).nbits(), value.denominator().nbits())
        for value in values
    )
    return values, max_bits


def candidates_from_windows(label, block_residues, windows):
    records = {}
    window_stats = []
    for indices in windows:
        before = time.monotonic()
        reduced = window_lattice(block_residues, indices)
        short_rows = reduced.rows()[:args.short_rank]
        bounds = range(-args.combination_bound, args.combination_bound + 1)
        tested = 0
        accepted = 0
        for coefficients in itertools.product(bounds, repeat=len(short_rows)):
            if not any(coefficients):
                continue
            tested += 1
            row = sum(
                (coefficient * short_rows[index]
                 for index, coefficient in enumerate(coefficients)),
                vector(ZZ, [0] * reduced.ncols()),
            )
            result = extended_candidate(block_residues, row[-1])
            if result is None:
                continue
            values, max_bits = result
            # Random centered vectors sit at essentially the full modulus bit
            # length.  Retain a generous 95% threshold for genuine shared
            # denominators and near multiples.
            if 20 * max_bits >= 19 * modulus.nbits():
                continue
            key = tuple(values)
            old = records.get(key)
            if old is None or max_bits < old["maximum_rational_bit_height"]:
                records[key] = {
                    "values": values,
                    "maximum_rational_bit_height": int(max_bits),
                    "window": list(indices),
                    "combination": list(coefficients),
                    "denominator_bit_height_before_reduction": int(abs(row[-1]).nbits()),
                }
                accepted += 1
        window_stats.append({
            "indices": list(indices),
            "lll_seconds": round(time.monotonic() - before, 6),
            "short_row_maximum_bits": [
                int(max(abs(value).nbits() for value in row)) for row in short_rows
            ],
            "combinations_tested": tested,
            "new_candidates_accepted": accepted,
        })
        print(
            f"A5A5P230SIMREC|block={label}|window={indices}|"
            f"candidates={len(records)}|elapsed={time.monotonic()-started:.3f}",
            flush=True,
        )
    ordered = sorted(records.values(), key=lambda record: record["maximum_rational_bit_height"])
    return ordered[:40], window_stats


x_candidates, x_window_stats = candidates_from_windows(
    "X", residues[:9], ((0, 1, 2, 3, 4), (2, 3, 4, 5, 6), (4, 5, 6, 7, 8)),
)
y_candidates, y_window_stats = candidates_from_windows(
    "Y", residues[9:22],
    ((0, 1, 2, 3, 4, 5), (3, 4, 5, 6, 7, 8), (7, 8, 9, 10, 11, 12)),
)

A_QQ = [QQ(value) for value in surface["child"]["minimal_A_coefficients_low_to_high"]]
B_QQ = [QQ(value) for value in surface["child"]["minimal_B_coefficients_low_to_high"]]
Z_QQ = [QQ(value) for value in fixed_Z]
auxiliary_primes = (1000003, 1000033, 1000037)


def reduce_polynomial(values, prime, name):
    field = GF(prime)
    ring = PolynomialRing(field, name)
    return ring([
        field(value.numerator()) / field(value.denominator()) for value in values
    ])


survivors = []
for x_index, x_record in enumerate(x_candidates):
    for y_index, y_record in enumerate(y_candidates):
        passed = True
        for auxiliary_prime in auxiliary_primes:
            try:
                X_F = reduce_polynomial(x_record["values"], auxiliary_prime, "T")
                Y_F = reduce_polynomial(y_record["values"], auxiliary_prime, "T")
                A_F = reduce_polynomial(A_QQ, auxiliary_prime, "T")
                B_F = reduce_polynomial(B_QQ, auxiliary_prime, "T")
                Z_F = reduce_polynomial(Z_QQ, auxiliary_prime, "T")
            except ZeroDivisionError:
                passed = False
                break
            if Y_F**2 != X_F**3 + A_F * X_F * Z_F**4 + B_F * Z_F**6:
                passed = False
                break
        if passed:
            survivors.append((x_index, y_index))

exact_section = None
RQQ = PolynomialRing(QQ, "T")
A_exact, B_exact, Z_exact = RQQ(A_QQ), RQQ(B_QQ), RQQ(Z_QQ)
for x_index, y_index in survivors:
    X_exact = RQQ(x_candidates[x_index]["values"])
    Y_exact = RQQ(y_candidates[y_index]["values"])
    if Y_exact**2 == X_exact**3 + A_exact * X_exact * Z_exact**4 + B_exact * Z_exact**6:
        exact_section = {
            "X_coefficients_low_to_high": [str(value) for value in X_exact.list()],
            "Y_coefficients_low_to_high": [str(value) for value in Y_exact.list()],
            "Z_coefficients_low_to_high": [str(value) for value in Z_exact.list()],
            "exact_Weierstrass_identity": True,
            "candidate_indices_X_Y": [x_index, y_index],
        }
        break

status = (
    "PASS_EXACT_QQ_P230_SIMULTANEOUS_RECONSTRUCTION"
    if exact_section is not None
    else "PASS_SIMULTANEOUS_RECONSTRUCTION_DIAGNOSTIC_NO_EXACT_PAIR"
)


def serializable_records(records):
    return [{
        **{key: value for key, value in record.items() if key != "values"},
        "values": [str(value) for value in record["values"]],
    } for record in records]


payload = {
    "schema": "elkies-k3.q24-2a5-p230-simultaneous-reconstruction.v1",
    "status": status,
    "software": "SageMath 10.9 (conda-forge pinned repository environment)",
    "prime": int(p),
    "p_adic_digits": digits,
    "modulus_bits": int(modulus.nbits()),
    "method": {
        "windowed_common_denominator_LLL": True,
        "combination_bound": args.combination_bound,
        "short_rank": args.short_rank,
        "random_vector_rejection_threshold": "maximum bit height < 0.95 modulus bits",
        "auxiliary_primes": list(auxiliary_primes),
    },
    "X": {
        "window_stats": x_window_stats,
        "candidates": serializable_records(x_candidates),
    },
    "Y": {
        "window_stats": y_window_stats,
        "candidates": serializable_records(y_candidates),
    },
    "candidate_pairs_surviving_auxiliary_primes": [list(pair) for pair in survivors],
    "exact_QQ_section": exact_section,
    "large_Groebner_required": False,
    "proof_boundary": (
        "Only a candidate passing literal QQ Weierstrass substitution is exact. "
        "Short vectors and auxiliary-prime filters alone are construction aids."
    ),
    "inputs": {
        "paths": [str(SURFACE), str(source_path)],
        "sha256": {str(path): sha256(path) for path in (SURFACE, source_path)},
    },
    "elapsed_seconds": round(time.monotonic() - started, 6),
}
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5A5P230SIMREC|X={}|Y={}|survivors={}|exact={}|status={}|output={}".format(
        len(x_candidates), len(y_candidates), len(survivors),
        exact_section is not None, status, output,
    ),
    flush=True,
)
