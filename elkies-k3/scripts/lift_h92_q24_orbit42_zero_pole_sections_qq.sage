#!/usr/bin/env sage -python
"""Hensel-lift the rational q24 D12 zero-pole sections to QQ.

The p=53 scan has 22 modular sections.  Exactly 18 degree-(4,6) sections pass
the four-digit Hensel obstruction test; all degree-(3,4) points fail by p^4.
This script lifts one representative of each of the nine surviving signed
x-coordinate pairs with its degree-adapted nonsingular
coefficient Jacobian, rationally reconstructs the result, and verifies the
section identity over QQ.  Signed partners are recovered by y -> -y.

The p-adic residues are checkpointed.  This establishes the 18 rational
identity-class sections on the D12 parent only; the missing nontrivial class,
MW matching, and the orbit42/A11 pencil remain downstream.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, QQ, ZZ, PolynomialRing, Zp, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
DEFAULT_SEED = LOCAL / "q24-orbit42-zero-pole-model-mod-53.json"
DEFAULT_MODEL = LOCAL / "q24-orbit42-zero-pole-seeds-mod-43.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--seed", type=Path, default=DEFAULT_SEED)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--precision", type=int, default=65536)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

seed_path = args.seed.resolve()
model_path = args.model.resolve()
seed = json.loads(seed_path.read_text())
model_cache = json.loads(model_path.read_text())
assert seed["status"] == "PASS_Q42_ZERO_POLE_MODEL_MODP_SCAN"
assert model_cache["status"] == "PASS_Q42_ZERO_POLE_SMALLPRIME_SEEDS"

p = ZZ(seed["prime"])
if p != 53:
    raise ValueError("the pinned corrected zero-pole seed prime is 53")
precision = int(args.precision)
if precision < 2:
    raise ValueError("precision must be at least two p-adic digits")

F = GF(p)
RQ = PolynomialRing(QQ, "u")
u = RQ.gen()
RF = PolynomialRing(F, "u")
uf = RF.gen()

exact_model = model_cache["exact_model"]
A = RQ([QQ(v) for v in exact_model["A_coefficients_low_to_high"]])
B = RQ([QQ(v) for v in exact_model["B_coefficients_low_to_high"]])
if (A.degree(), B.degree()) != (6, 9):
    raise ArithmeticError("cached exact D12 model has wrong degrees")


def red_q(value):
    value = QQ(value)
    denominator = ZZ(value.denominator())
    if denominator % p == 0:
        raise ZeroDivisionError("model denominator is divisible by seed prime")
    return F(ZZ(value.numerator())) / F(denominator)


Af = RF([red_q(v) for v in A.list()])
Bf = RF([red_q(v) for v in B.list()])


def polynomial_from_seed(record, key):
    return RF([F(v) for v in record[key]])


records = []
for record in seed["sections"]:
    if int(record.get("hensel_depth", 1)) < 4:
        continue
    x_value = polynomial_from_seed(record, "x_coefficients_low_to_high")
    y_value = polynomial_from_seed(record, "y_coefficients_low_to_high")
    if y_value**2 != x_value**3 + Af * x_value + Bf:
        raise ArithmeticError("modular seed misses the exact model reduction")
    records.append((record, x_value, y_value))
if len(records) != 18:
    raise ArithmeticError(
        f"expected 18 depth-four modular sections, found {len(records)}"
    )

# Pair P and -P by their common x coordinate.
by_x = {}
for record, x_value, y_value in records:
    by_x.setdefault(tuple(x_value.list()), []).append((record, x_value, y_value))
if len(by_x) != 9 or sorted(len(group) for group in by_x.values()) != [2] * 9:
    raise ArithmeticError("p=53 depth-four seeds do not form nine signed pairs")

pairs = []
for key in sorted(by_x, key=str):
    group = by_x[key]
    left, right = group
    if left[1] != right[1] or left[2] != -right[2]:
        raise ArithmeticError("modular signed section pair is inconsistent")
    pairs.append(min(group, key=lambda row: tuple(map(int, row[2].list()))))


def padded(poly, length, zero):
    return [poly[i] if i <= poly.degree() else zero for i in range(length)]


def identity_and_jacobian(ring, x_values, y_values, dx, dy, equation_degree):
    uu = ring.gen()
    base = ring.base_ring()
    xx = sum(base(x_values[i]) * uu**i for i in range(dx + 1))
    yy = sum(base(y_values[i]) * uu**i for i in range(dy + 1))
    aa = ring(A)
    bb = ring(B)
    residual = yy**2 - xx**3 - aa * xx - bb
    derivatives = []
    for i in range(dx + 1):
        derivatives.append((-3 * xx**2 - aa) * uu**i)
    for i in range(dy + 1):
        derivatives.append(2 * yy * uu**i)
    zero = base.zero()
    residual_vector = vector(
        base, padded(residual, equation_degree + 1, zero)
    )
    jacobian = matrix(base, [
        padded(derivative, equation_degree + 1, zero)
        for derivative in derivatives
    ]).transpose()
    return residual_vector, jacobian


checkpoint_path = LOCAL / f"q24-orbit42-zero-pole-hensel-p{p}.json"
checkpoint = None
if checkpoint_path.exists():
    try:
        candidate = json.loads(checkpoint_path.read_text())
        if candidate.get("prime") == int(p):
            checkpoint = candidate
    except Exception:
        checkpoint = None

padic = Zp(p, prec=precision)
RP = PolynomialRing(padic, "u")
modulus = p**precision


def valuation_floor(values):
    valuations = [value.valuation() for value in values if value]
    return precision if not valuations else min(valuations)


def reconstruct(value):
    try:
        return QQ(ZZ(value.lift()).rational_reconstruction(modulus))
    except (ArithmeticError, ValueError):
        return None


lifted = []
checkpoint_records = []

for pair_index, (seed_record, x_seed, y_seed) in enumerate(pairs):
    dx = int(x_seed.degree())
    dy = int(y_seed.degree())
    if (dx, dy) not in ((3, 4), (4, 6)):
        raise ArithmeticError(f"unexpected zero-pole degree pair {(dx, dy)}")
    equation_degree = max(2 * dy, 3 * dx, A.degree() + dx, B.degree())
    x0 = padded(x_seed, dx + 1, F.zero())
    y0 = padded(y_seed, dy + 1, F.zero())
    residual_f, jacobian_f = identity_and_jacobian(
        RF, x0, y0, dx, dy, equation_degree
    )
    if residual_f:
        raise ArithmeticError("nonzero residual at modular seed")
    unknowns = dx + dy + 2
    rank = int(jacobian_f.rank())
    if rank != unknowns:
        raise ArithmeticError(
            f"degree-adapted seed {pair_index} rank={rank}, unknowns={unknowns}"
        )
    rows = tuple(map(int, jacobian_f.transpose().pivots()))
    square_f = jacobian_f.matrix_from_rows(rows)
    if not square_f.is_invertible():
        raise ArithmeticError("selected modular Hensel minor is singular")

    values = [padic(ZZ(v)) for v in x0 + y0]
    resumed = 1
    if checkpoint is not None:
        old = next(
            (
                row for row in checkpoint.get("records", [])
                if int(row.get("pair_index", -1)) == pair_index
                and row.get("degrees") == [dx, dy]
            ),
            None,
        )
        if old is not None and len(old.get("residues", [])) == len(values):
            old_precision = int(checkpoint.get("precision", 1))
            if 1 < old_precision < precision:
                values = [padic(ZZ(v)) for v in old["residues"]]
                resumed = old_precision

    for iteration in range(1, 2 * precision + 4):
        residual_p, jacobian_p = identity_and_jacobian(
            RP,
            values[:dx + 1],
            values[dx + 1:],
            dx,
            dy,
            equation_degree,
        )
        valuation = valuation_floor(residual_p)
        print(
            "Q42ZPQQ_HENSEL|"
            f"pair={pair_index}|degrees={dx},{dy}|iteration={iteration}|"
            f"valuation={valuation}|target={precision}",
            flush=True,
        )
        if valuation >= precision:
            break
        correction = jacobian_p.matrix_from_rows(rows).solve_right(
            -vector(padic, [residual_p[row] for row in rows])
        )
        correction_valuations = [
            delta.valuation() for delta in correction if delta
        ]
        print(
            "Q42ZPQQ_CORRECTION|"
            f"pair={pair_index}|iteration={iteration}|"
            f"valuation={min(correction_valuations) if correction_valuations else precision}",
            flush=True,
        )
        values = [value + delta for value, delta in zip(values, correction)]
    else:
        raise ArithmeticError("p-adic Newton iteration did not converge")

    rationals = [reconstruct(value) for value in values]
    complete = not any(value is None for value in rationals)
    checkpoint_records.append({
        "pair_index": pair_index,
        "degrees": [dx, dy],
        "jacobian_rank": rank,
        "selected_rows": list(rows),
        "seed_index": int(seed_record["index"]),
        "resumed_precision": resumed,
        "residues": [str(ZZ(value.lift())) for value in values],
        "reconstruction_complete": complete,
    })
    if not complete:
        continue

    xq = RQ(rationals[:dx + 1])
    yq = RQ(rationals[dx + 1:])
    if yq**2 != xq**3 + A * xq + B:
        raise ArithmeticError("rationally reconstructed section identity failed")
    if RF([red_q(v) for v in xq.list()]) != x_seed:
        raise ArithmeticError("exact x section misses modular seed")
    if RF([red_q(v) for v in yq.list()]) != y_seed:
        raise ArithmeticError("exact y section misses modular seed")

    lifted.append((pair_index, xq, yq, dx, dy))
    print(
        "Q42ZPQQ_SECTION|"
        f"pair={pair_index}|degrees={dx},{dy}|identity=1|modp=1|status=PASS",
        flush=True,
    )

checkpoint_payload = {
    "schema": "elkies-k3.h3-q24-orbit42-zero-pole-hensel.v1",
    "status": "PASS_Q42_ZERO_POLE_PADIC_RESIDUES",
    "prime": int(p),
    "precision": precision,
    "records": checkpoint_records,
}
checkpoint_path.write_text(json.dumps(checkpoint_payload, indent=2, sort_keys=True) + "\n")

if len(lifted) != len(pairs):
    print(
        "Q42ZPQQ_NEEDS_MORE_PRECISION|"
        f"precision={precision}|complete={len(lifted)}/{len(pairs)}|"
        f"checkpoint={checkpoint_path}|suggested={precision*2}|"
        "status=NEEDS_MORE_PRECISION",
        flush=True,
    )
    raise SystemExit(2)


def qlist(poly):
    return [str(v) for v in poly.list()]


section_records = []
for pair_index, xq, yq, dx, dy in lifted:
    for sign in (1, -1):
        section_records.append({
            "pair_index": pair_index,
            "sign": sign,
            "degrees": [dx, dy],
            "x_coefficients_low_to_high": qlist(xq),
            "y_coefficients_low_to_high": qlist(sign * yq),
        })

payload = {
    "schema": "elkies-k3.h3-q24-orbit42-rational-zero-pole-sections-qq.v1",
    "status": "PASS_EXACT_Q42_RATIONAL_ZERO_POLE_SECTIONS_QQ",
    "inputs": {
        "seed": str(seed_path.relative_to(ROOT)),
        "model": str(model_path.relative_to(ROOT)),
        "checkpoint": str(checkpoint_path.relative_to(ROOT)),
    },
    "prime": int(p),
    "precision": precision,
    "signed_pair_count": len(lifted),
    "section_count": len(section_records),
    "sections": section_records,
    "verification": {
        "degree_adapted_nonsingular_hensel": True,
        "exact_section_identities": True,
        "modular_regression": True,
    },
    "proof_boundary": (
        "Eighteen explicit characteristic-zero rational zero-pole sections on "
        "the exact q24-derived D12 parent. The nontrivial discriminant-class "
        "pair, MW matching, orbit42 section/pencil, and A11 child remain open."
    ),
}

OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q24-orbit42-rational-zero-pole-sections-qq.json"
)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q42ZPQQ_RESULT|pairs=9|sections=18|"
    "status=PASS_EXACT_Q42_RATIONAL_ZERO_POLE_SECTIONS_QQ",
    flush=True,
)
