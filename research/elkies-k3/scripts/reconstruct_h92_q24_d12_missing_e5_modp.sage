#!/usr/bin/env sage -python
"""Recover the low-pole missing D12 section before characteristic-zero lift.

The degree-14 ``close_P24`` Abel--Jacobi trace has pole order 104 when it is
reconstructed by itself.  The exact selected-marking lattice word is

    E5 = AJ(close_P24) - 7*S2 - Spinor,

where abstract identity-shell vector S2 maps to equation-shell point 1 and
the positive spinor maps to equation point 0.  Apply this word independently
on every sampled smooth fibre, then reconstruct only the expected P.O=4
section.  Both AJ and spinor orientations are tested explicitly.  Everything
here is finite-field elliptic-curve arithmetic and univariate interpolation;
no multivariate elimination or Groebner basis is used.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, is_prime


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--samples",
    type=Path,
    default=LOCAL / "q24-close-p24-aj14-plus-700samples-mod100003.json",
)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-d12-missing-e5-section-mod100003.json",
)
args = parser.parse_args()

SAMPLES = args.samples.resolve()
Q24 = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
A11 = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
ZERO = LOCAL / "q24-orbit42-rational-zero-pole-sections-qq.json"
SPINOR = LOCAL / "q24-orbit42-spinor-zero-pole-sections-qq.json"
MATCHING = LOCAL / "q24-orbit42-identity-halving-qq.json"
MARKING = LOCAL / "q24-a11-equation-marking-orbit64-mod100003.json"
ROUTE = GENERATED / "elkies-k3-h3-q24-d12-degree14-aj-missing-parent-route.json"
INPUTS = (SAMPLES, Q24, A11, ZERO, SPINOR, MATCHING, MARKING, ROUTE)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

samples_artifact = json.loads(SAMPLES.read_text())
q24 = json.loads(Q24.read_text())
a11 = json.loads(A11.read_text())
zero = json.loads(ZERO.read_text())
spinor = json.loads(SPINOR.read_text())
matching = json.loads(MATCHING.read_text())
marking = json.loads(MARKING.read_text())
route = json.loads(ROUTE.read_text())

assert zero["status"] == "PASS_EXACT_Q42_RATIONAL_ZERO_POLE_SECTIONS_QQ"
assert spinor["status"] == "PASS_EXACT_Q42_SPINOR_ZERO_POLE_SECTIONS_QQ"
assert matching["status"] == "Q42_IDENTITY_HALVING_HAS_NO_A11_CHORD"
assert marking["status"] == "PASS_Q42_A11_EQUATION_MARKING_ORBIT64_MOD100003"
assert route["status"] == "PASS_EXACT_Q24_D12_DEGREE14_AJ_MISSING_PARENT_ROUTE"

p = ZZ(samples_artifact["prime"])
if not is_prime(p) or p in (2, 3):
    raise ArithmeticError("sample characteristic is not good")
if p != marking["good_reduction_prime"]:
    raise ArithmeticError("samples do not use the pinned equation marking prime")
F = GF(p)
R = PolynomialRing(F, "V")
V = R.gen()
VQ = PolynomialRing(QQ, "V")
KQ = VQ.fraction_field()


def red(value):
    value = QQ(value)
    if value.denominator() % p == 0:
        raise ZeroDivisionError(f"bad denominator modulo {p}")
    return F(value.numerator()) / F(value.denominator())


def qq_rational(text):
    return KQ(str(text))


def evaluate_rational(value, tau):
    value = KQ(value)
    numerator = sum(red(value.numerator()[i]) * tau**i for i in range(value.numerator().degree() + 1))
    denominator = sum(red(value.denominator()[i]) * tau**i for i in range(value.denominator().degree() + 1))
    if not denominator:
        raise ZeroDivisionError("coordinate-change denominator vanished")
    return numerator / denominator


def evaluate_coefficients(values, argument):
    return sum(red(value) * argument**i for i, value in enumerate(values))


# The orbit64 marking selected mapping 7.  Its abstract identity vector 2 is
# the S2 occurring in the shortest bridge word, and maps to equation point 1.
mapping_index = int(marking["selected"]["mapping_index"])
assert mapping_index == 7
abstract_s2_index = 2
equation_s2_index = int(
    matching["matching"]["mappings_abstract_to_equation"][mapping_index][abstract_s2_index]
)
assert equation_s2_index == 1
word = route["parent_bridge_word"]
assert word["identity_shell_coefficients"][abstract_s2_index] == -7
assert sum(abs(value) for value in word["identity_shell_coefficients"]) == 7
assert word["spinor_coefficient_for_vector_1_0_0_0_0"] == -1

s2_row = zero["sections"][equation_s2_index]
spinor_rows = spinor["sections"]
assert len(spinor_rows) == 2 and spinor_rows[0]["sign"] == 1 and spinor_rows[1]["sign"] == -1

u_of_V = qq_rational(a11["coordinate_change"]["u_of_V"])
x_scale = qq_rational(a11["coordinate_change"]["x_scale"])
y_scale = qq_rational(a11["coordinate_change"]["y_scale"])
A = R([red(value) for value in q24["child"]["minimal_A_coefficients_low_to_high"]])
B = R([red(value) for value in q24["child"]["minimal_B_coefficients_low_to_high"]])

raw_samples = samples_artifact["q24_stage"]["degree14_AJ_traces"]
if {row["branch"] for row in raw_samples} != {"pole_plus"}:
    raise ArithmeticError("expected a single pole_plus sample stream")
if len({int(row["tau"]) for row in raw_samples}) != len(raw_samples):
    raise ArithmeticError("duplicate sample abscissa")


def shell_point(row, tau, curve):
    u_value = evaluate_rational(u_of_V, tau)
    x = evaluate_rational(x_scale, tau) * evaluate_coefficients(
        row["x_coefficients_low_to_high"], u_value
    )
    y = evaluate_rational(y_scale, tau) * evaluate_coefficients(
        row["y_coefficients_low_to_high"], u_value
    )
    return curve(x, y)


combined = {(aj_sign, spinor_index): [] for aj_sign in (1, -1) for spinor_index in (0, 1)}
skipped = 0
for row in raw_samples:
    tau = F(row["tau"])
    curve = EllipticCurve(F, [0, 0, 0, A(tau), B(tau)])
    try:
        aj = curve(F(row["AJ_x"]), F(row["AJ_y"]))
        s2 = shell_point(s2_row, tau, curve)
        spins = [shell_point(item, tau, curve) for item in spinor_rows]
    except (ArithmeticError, ZeroDivisionError, ValueError):
        skipped += 1
        continue
    for aj_sign in (1, -1):
        for spinor_index, spin_point in enumerate(spins):
            point = aj_sign * aj - 7 * s2 - spin_point
            if point.is_zero():
                raise ArithmeticError("candidate bridge word specialized to zero")
            x, y = point.xy()
            combined[(aj_sign, spinor_index)].append((tau, F(x), F(y)))


def interpolation_polynomial(values):
    interpolation = R.zero()
    modulus = R.one()
    for parameter, value in values:
        scale = modulus(parameter)
        if not scale:
            raise ArithmeticError("duplicate interpolation parameter")
        interpolation += ((value - interpolation(parameter)) / scale) * modulus
        modulus *= V - parameter
    interpolation %= modulus
    return interpolation, modulus


def pade_sequence(values):
    interpolation, modulus = interpolation_polynomial(values)
    r0, r1 = modulus, interpolation
    t0, t1 = R.zero(), R.one()
    candidates = []
    while r1:
        numerator, denominator = r1, t1
        common = numerator.gcd(denominator)
        numerator //= common
        denominator //= common
        if denominator:
            scale = denominator.leading_coefficient()
            numerator /= scale
            denominator /= scale
            if all(
                denominator(parameter)
                and numerator(parameter) == value * denominator(parameter)
                for parameter, value in values
            ):
                candidates.append((numerator, denominator))
        quotient, r2 = r0.quo_rem(r1)
        t2 = t0 - quotient * t1
        r0, r1 = r1, r2
        t0, t1 = t1, t2
    return candidates


def reconstruct_section(samples):
    x_candidates = []
    for numerator, denominator in pade_sequence([(t, x) for t, x, unused in samples]):
        if denominator.degree() != 8 or not denominator.is_square():
            continue
        Z = denominator.sqrt().monic()
        denominator = Z**2
        scale = denominator.leading_coefficient()
        numerator /= scale
        if Z.degree() == 4 and numerator.degree() <= 12:
            x_candidates.append((numerator, Z))
    results = []
    for X, Z in x_candidates:
        for numerator, denominator in pade_sequence([(t, y) for t, unused, y in samples]):
            scale = denominator.leading_coefficient()
            numerator /= scale
            denominator /= scale
            if denominator != Z**3 or numerator.degree() > 18:
                continue
            Y = numerator
            if Y**2 != X**3 + A * X * Z**4 + B * Z**6:
                continue
            results.append((X, Y, Z))
    unique = {}
    for X, Y, Z in results:
        unique[(tuple(X.list()), tuple(Y.list()), tuple(Z.list()))] = (X, Y, Z)
    return list(unique.values())


orientation_records = []
solutions = []
for (aj_sign, spinor_index), values in combined.items():
    reconstructed = reconstruct_section(values)
    orientation_records.append(
        {
            "AJ_sign": aj_sign,
            "spinor_equation_index": spinor_index,
            "usable_samples": len(values),
            "P_dot_O_4_reconstruction_count": len(reconstructed),
        }
    )
    for X, Y, Z in reconstructed:
        solutions.append((aj_sign, spinor_index, X, Y, Z))

if len(solutions) != 1:
    raise ArithmeticError(f"expected one oriented P.O=4 reconstruction, got {len(solutions)}")
aj_sign, spinor_index, X, Y, Z = solutions[0]

payload = {
    "schema": "elkies-k3.h3-q24-d12-missing-e5-section-modp.v1",
    "status": "PASS_Q24_D12_MISSING_E5_SECTION_RECONSTRUCTION_MODP",
    "prime": int(p),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "bridge_word": {
        "formula": "E5=oriented_AJ(close_P24)-7*S2-Spinor",
        "abstract_identity_index": abstract_s2_index,
        "equation_identity_index": equation_s2_index,
        "selected_mapping_index": mapping_index,
        "AJ_sign": aj_sign,
        "spinor_equation_index": spinor_index,
    },
    "orientation_audit": orientation_records,
    "sample_count": len(raw_samples),
    "skipped_bad_shell_specializations": skipped,
    "section": {
        "X_coefficients_low_to_high": [int(value) for value in X.list()],
        "Y_coefficients_low_to_high": [int(value) for value in Y.list()],
        "Z_coefficients_low_to_high": [int(value) for value in Z.list()],
        "degrees_X_Y_Z": [int(X.degree()), int(Y.degree()), int(Z.degree())],
        "P_dot_O": 4,
        "D12_MW_Abel_Jacobi": [0, 0, 0, 0, 1],
        "exact_modp_weierstrass_identity": True,
    },
    "method": "fibrewise exact shell combination followed by univariate Newton/Pade reconstruction",
    "large_Groebner_required": False,
    "proof_boundary": (
        "Exact over the pinned good-reduction field with an exact characteristic-zero shell word. "
        "The displayed section still requires characteristic-zero coefficient lifting and literal QQ verification."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q24E5MODP|prime={}|AJ_sign={}|spinor={}|samples={}|degrees={},{},{}|PO=4|status={}".format(
        p,
        aj_sign,
        spinor_index,
        len(combined[(aj_sign, spinor_index)]),
        X.degree(),
        Y.degree(),
        Z.degree(),
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
