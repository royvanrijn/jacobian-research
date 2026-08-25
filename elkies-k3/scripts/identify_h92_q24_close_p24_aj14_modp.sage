#!/usr/bin/env sage -python
"""Match degree-14 close_P24 AJ traces against the exact D12 zero-pole shell."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, is_prime


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument(
    "--trace",
    type=Path,
    default=LOCAL / "q24-a11-close-p24-quintic-modp.json",
)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-close-p24-aj14-shell-match-modp.json",
)
args = parser.parse_args()

p = ZZ(args.prime)
if not is_prime(p) or p in (2, 3):
    raise SystemExit("--prime must be a prime other than 2 or 3")

TRACE = args.trace.resolve()
A11 = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
ZERO = LOCAL / "q24-orbit42-rational-zero-pole-sections-qq.json"
SPINOR = LOCAL / "q24-orbit42-spinor-zero-pole-sections-qq.json"
Q24 = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
INPUTS = (TRACE, A11, ZERO, SPINOR, Q24)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

trace = json.loads(TRACE.read_text())
a11 = json.loads(A11.read_text())
zero = json.loads(ZERO.read_text())
spinor = json.loads(SPINOR.read_text())
q24 = json.loads(Q24.read_text())
traces = trace["q24_stage"]["degree14_AJ_traces"]
assert len(traces) == 2 and all(int(row["degree"]) == 14 for row in traces)

F = GF(p)
UQ = PolynomialRing(QQ, "u")
VQ = PolynomialRing(QQ, "V")
KVQ = VQ.fraction_field()


def red(value):
    value = QQ(value)
    if value.denominator() % p == 0:
        raise ZeroDivisionError(f"bad denominator modulo {p}")
    return F(value.numerator()) / F(value.denominator())


def evaluate_qq_rational(text, tau):
    value = KVQ(str(text))
    numerator = sum(red(value.numerator()[i]) * tau**i for i in range(value.numerator().degree() + 1))
    denominator = sum(red(value.denominator()[i]) * tau**i for i in range(value.denominator().degree() + 1))
    if not denominator:
        raise ZeroDivisionError("coordinate change denominator vanished")
    return numerator / denominator


tau_values = {int(row["tau"]) for row in traces}
assert len(tau_values) == 1
tau = F(next(iter(tau_values)))
u_value = evaluate_qq_rational(a11["coordinate_change"]["u_of_V"], tau)
x_scale = evaluate_qq_rational(a11["coordinate_change"]["x_scale"], tau)
y_scale = evaluate_qq_rational(a11["coordinate_change"]["y_scale"], tau)

shell_points = []
for kind, rows in (("identity", zero["sections"]), ("spinor", spinor["sections"])):
    for index, row in enumerate(rows):
        x_coefficients = [red(value) for value in row["x_coefficients_low_to_high"]]
        y_coefficients = [red(value) for value in row["y_coefficients_low_to_high"]]
        x = x_scale * sum(value * u_value**i for i, value in enumerate(x_coefficients))
        y = y_scale * sum(value * u_value**i for i, value in enumerate(y_coefficients))
        shell_points.append((kind, index, int(x), int(y)))

rows = []
for item in traces:
    matches = [
        {"shell_kind": kind, "section_index": index, "same_y_sign": int(item["AJ_y"]) == y}
        for kind, index, x, y in shell_points
        if int(item["AJ_x"]) == x and int(item["AJ_y"]) in (y, int(-F(y)))
    ]
    rows.append({**item, "shell_matches": matches})

A_coefficients = [red(value) for value in q24["child"]["minimal_A_coefficients_low_to_high"]]
B_coefficients = [red(value) for value in q24["child"]["minimal_B_coefficients_low_to_high"]]
A_value = sum(value * tau**i for i, value in enumerate(A_coefficients))
B_value = sum(value * tau**i for i, value in enumerate(B_coefficients))
curve = EllipticCurve(F, [0, 0, 0, A_value, B_value])
trace_points = {
    row["branch"]: curve(F(row["AJ_x"]), F(row["AJ_y"])) for row in rows
}


def point_shell_matches(point):
    if point.is_zero():
        return [{"shell_kind": "zero", "section_index": None, "same_y_sign": True}]
    px, py = map(int, point.xy())
    return [
        {"shell_kind": kind, "section_index": index, "same_y_sign": py == y}
        for kind, index, x, y in shell_points
        if px == x and py in (y, int(-F(y)))
    ]


combination_matches = []
for plus_coefficient in range(-3, 4):
    for minus_coefficient in range(-3, 4):
        if not plus_coefficient and not minus_coefficient:
            continue
        value = (
            plus_coefficient * trace_points["pole_plus"]
            + minus_coefficient * trace_points["pole_minus"]
        )
        matches = point_shell_matches(value)
        if matches:
            combination_matches.append(
                {
                    "pole_plus_coefficient": plus_coefficient,
                    "pole_minus_coefficient": minus_coefficient,
                    "shell_matches": matches,
                }
            )

involution_matches = []
trace_sum = trace_points["pole_plus"] + trace_points["pole_minus"]
for kind, index, x, y in shell_points:
    shell_point = curve(F(x), F(y))
    for multiplier in (14, -14):
        if trace_sum == multiplier * shell_point:
            involution_matches.append(
                {
                    "relation": "AJ(pole_plus)+AJ(pole_minus)=m*R",
                    "multiplier": multiplier,
                    "shell_kind": kind,
                    "section_index": index,
                }
            )

payload = {
    "schema": "elkies-k3.h3-q24-close-p24-aj14-shell-match-modp.v1",
    "status": "PASS_Q24_CLOSE_P24_AJ14_SHELL_MATCH_AUDIT_MODP",
    "prime": int(p),
    "tau": int(tau),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "traces": rows,
    "small_combination_shell_matches": combination_matches,
    "quartic_involution_trace_matches": involution_matches,
    "proof_boundary": (
        "Exact arithmetic over one finite field compares the two degree-14 AJ "
        "branches with all twenty exact zero-pole shell points. A nonmatch is not "
        "a characteristic-zero MW identification."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q24CLOSEAJMATCH|prime={}|tau={}|matches={}|status={}".format(
        p,
        int(tau),
        ";".join(f"{row['branch']}:{len(row['shell_matches'])}" for row in rows),
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
