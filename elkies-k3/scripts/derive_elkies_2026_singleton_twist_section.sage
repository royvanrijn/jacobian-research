#!/usr/bin/env sage-python
"""Descend one certified bisection lift to its quadratic-twist section.

For a record with ``u^2=q(t)`` and lifted point ``P``, form

    R = P - sigma(P).

Then ``x(R)`` is invariant and ``y(R)=u*y_1(t)`` is anti-invariant.  On the
integral short twist

    Y^2 = X^3 + A*q^2*X + B*q^3

the descended point is ``(X,Y)=(q*x(R),q^2*y_1)``.  The script verifies all
identities exactly and also records the square-equivalent integral-q model
used by the modular polynomial-section exporter.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import lcm
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, sage_eval


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BISECTIONS = (
    ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
)
DEFAULT_MODEL = (
    ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
)


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def polynomial(ring, values):
    return ring([QQ(value) for value in values])


def polynomial_coefficients(value):
    return [str(coefficient) for coefficient in value.list()]


def require_polynomial(value, ring, label):
    value = value.numerator() / value.denominator()
    if value.denominator().degree() != 0:
        raise ArithmeticError(f"{label} is not a polynomial")
    return ring(value)


parser = argparse.ArgumentParser(description=__doc__)
target = parser.add_mutually_exclusive_group(required=True)
target.add_argument("--mask", type=int)
target.add_argument("--direct-label")
parser.add_argument("--prime", type=int, help="optionally identify its exporter block")
parser.add_argument("--bisections", type=Path, default=DEFAULT_BISECTIONS)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--export", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

bisections = json.loads(args.bisections.read_text())
if bisections.get("schema") != "elkies-k3.bisection-extension-input.v1":
    raise ValueError("unexpected bisection batch schema")
record = next(
    (
        item
        for item in bisections["bisections"]
        if (
            int(item["lattice_orbit_mask"]) == args.mask
            if args.mask is not None
            else item["label"] == args.direct_label
        )
    ),
    None,
)
if record is None:
    raise ValueError(
        f"unknown bisection target {args.mask if args.mask is not None else args.direct_label}"
    )
mask = int(record["lattice_orbit_mask"])
direct = args.direct_label is not None

model = json.loads(args.model.read_text())
if model.get("status") == "PASS_TRANSCRIBED_PUBLISHED_R17_MODEL":
    model_coefficients = model
elif model.get("status") == "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    model_coefficients = model["weierstrass_model"]
else:
    raise ValueError("expected the certified compact published R17 model")

R = PolynomialRing(QQ, "t")
t = R.gen()
K = R.fraction_field()
A = polynomial(R, model_coefficients["A_coefficients_low_to_high"])
B = polynomial(R, model_coefficients["B_coefficients_low_to_high"])
q = polynomial(
    R,
    record["branch"]["numerator_coefficients"]
    if direct
    else record["residual_chord"]["q_coefficients"],
)
if q.degree() != 2 or q.gcd(q.derivative()).degree() != 0:
    raise ArithmeticError("expected a squarefree quadratic twist polynomial")

extension_ring = PolynomialRing(K, "z")
z = extension_ring.gen()
L = K.extension(z**2 - q, "u")
u = L.gen()
E = EllipticCurve(L, [A, B])
lift = record["lifted_section"]
x0 = polynomial(R, lift["x0_coefficients"])
x1 = polynomial(R, lift["x1_coefficients"])
y0 = polynomial(R, lift["y0_coefficients"])
y1 = polynomial(R, lift["y1_coefficients"])
P = E(L(x0) + L(x1) * u, L(y0) + L(y1) * u)
sigma_P = E(L(x0) - L(x1) * u, L(y0) - L(y1) * u)
anti = P - sigma_P

x_coefficients = anti[0].list()
y_coefficients = anti[1].list()
if len(x_coefficients) > 1 and any(x_coefficients[1:]):
    raise ArithmeticError("x(P-sigma(P)) is not invariant")
if not y_coefficients or y_coefficients[0] != 0 or any(y_coefficients[2:]):
    raise ArithmeticError("y(P-sigma(P)) is not purely anti-invariant")
x_anti = K(x_coefficients[0])
y_anti_u = K(y_coefficients[1])

X = require_polynomial(K(q) * x_anti, R, "twist X")
Y = require_polynomial(K(q) ** 2 * y_anti_u, R, "twist Y")
if Y**2 != X**3 + A * q**2 * X + B * q**3:
    raise ArithmeticError("descended twist section fails the exact curve identity")
if X.degree() > 6 or Y.degree() > 9:
    raise ArithmeticError("descended section exceeds the chi=3 polynomial bounds")

# This is exactly the square-equivalent q normalization used by the census
# and modular exporter: if D clears q's denominators, q_int=D^2*q.
q_denominator = lcm(*(coefficient.denominator() for coefficient in q.list()))
q_integral = R(q_denominator**2 * q)
X_integral_q = R(q_denominator**2 * X)
Y_integral_q = R(q_denominator**3 * Y)
if Y_integral_q**2 != (
    X_integral_q**3
    + A * q_integral**2 * X_integral_q
    + B * q_integral**3
):
    raise ArithmeticError("section scaling to the integral-q model failed")

modular_identification = None
if args.prime is not None:
    if args.export is None:
        args.export = (
            ROOT
            / "artifacts/local/elkies-k3/twist-polynomial-sections"
            / (
                f"direct-singleton-{args.direct_label}/p{args.prime}/export.json"
                if direct
                else f"singleton-{args.mask}/p{args.prime}/export.json"
            )
        )
    export = json.loads(args.export.read_text())
    expected_kind = "direct_singleton" if direct else "singleton"
    expected_key = args.direct_label if direct else str(args.mask)
    if (
        export["candidate"]["kind"] != expected_kind
        or export["candidate"]["key"] != expected_key
    ):
        raise ValueError("modular export does not match the requested singleton")
    if int(export["prime"]) != args.prime:
        raise ValueError("modular export prime mismatch")
    prime = int(args.prime)
    if any(coefficient.denominator() % prime == 0 for coefficient in X_integral_q.list() + Y_integral_q.list()):
        raise ArithmeticError("known section is not integral at the requested prime")
    field = GF(prime)
    modular_ring = PolynomialRing(field, "s")
    s = modular_ring.gen()
    chart = export["infinity_fibre"]["chart"]
    if chart == "original_infinity":
        transformed_X = modular_ring([field(value) for value in X_integral_q])
        transformed_Y = modular_ring([field(value) for value in Y_integral_q])
    else:
        prefix = "original_t="
        if not chart.startswith(prefix) or " via " not in chart:
            raise ValueError(f"unrecognized exporter chart {chart}")
        chart_parameter = int(chart[len(prefix) :].split(" via ", 1)[0])
        chart_parameter_field = field(chart_parameter)
        transformed_X = modular_ring(
            sum(
                field(X_integral_q[index])
                * (chart_parameter_field * s + 1) ** index
                * s ** (6 - index)
                for index in range(X_integral_q.degree() + 1)
            )
        )
        transformed_Y = modular_ring(
            sum(
                field(Y_integral_q[index])
                * (chart_parameter_field * s + 1) ** index
                * s ** (9 - index)
                for index in range(Y_integral_q.degree() + 1)
            )
        )
    leading_x = transformed_X[6]
    leading_y = transformed_Y[9]
    matches = [
        int(system["block_index"])
        for system in export["systems"]
        if system["leading_x_y"] == [int(leading_x), int(leading_y)]
    ]
    negative_matches = [
        int(system["block_index"])
        for system in export["systems"]
        if system["leading_x_y"] == [int(leading_x), int(-leading_y) % prime]
    ]
    if len(matches) != 1 or len(negative_matches) != 1:
        raise ArithmeticError("known section did not identify unique modular blocks")
    system_record = export["systems"][matches[0]]
    system_path = ROOT / system_record["path"]
    system_lines = system_path.read_text().splitlines()
    system_names = tuple(system_lines[0].split(","))
    if int(system_lines[1]) != prime or system_names != tuple(
        f"x{index}" for index in range(5, -1, -1)
    ):
        raise ValueError("unexpected exported msolve system header")
    coefficient_ring = PolynomialRing(field, names=system_names, order="degrevlex")
    local_variables = coefficient_ring.gens_dict()
    substitutions = {
        local_variables[f"x{index}"]: transformed_X[index]
        for index in range(6)
    }
    equation_text = "\n".join(system_lines[2:])
    equations = [
        coefficient_ring(
            sage_eval(piece.replace("^", "**"), locals=local_variables)
        )
        for piece in equation_text.split(",\n")
    ]
    if len(equations) != 9 or any(equation.subs(substitutions) for equation in equations):
        raise ArithmeticError("known section does not solve its exported modular system")
    modular_identification = {
        "prime": prime,
        "export": relative(args.export),
        "export_sha256": digest(args.export),
        "leading_x_y": [int(leading_x), int(leading_y)],
        "section_block": matches[0],
        "negative_section_block": negative_matches[0],
        "section_X_coefficients_low_to_high": [
            int(transformed_X[index]) for index in range(7)
        ],
        "section_Y_coefficients_low_to_high": [
            int(transformed_Y[index]) for index in range(10)
        ],
        "msolve_X_variables_x5_through_x0": [
            int(transformed_X[index]) for index in range(5, -1, -1)
        ],
        "exported_system": relative(system_path),
        "exported_system_sha256": digest(system_path),
        "exported_equation_residuals_zero": True,
    }

payload = {
    "schema": "elkies-k3.elkies-2026-singleton-twist-section.v1",
    "status": "PASS_EXACT_DESCENDED_SINGLETON_TWIST_SECTION",
    "label": record.get("label"),
    "mask": mask,
    "orbit_hex": f"0x{mask:05x}",
    "construction": "R=P-sigma(P); (X,Y)=(q*x(R),q^2*coefficient_u(y(R)))",
    "twist_model": "Y^2=X^3+A*q^2*X+B*q^3",
    "q_coefficients_low_to_high": polynomial_coefficients(q),
    "X_coefficients_low_to_high": polynomial_coefficients(X),
    "Y_coefficients_low_to_high": polynomial_coefficients(Y),
    "degrees_X_Y": [int(X.degree()), int(Y.degree())],
    "integral_q_model": {
        "scale_D": int(q_denominator),
        "q_coefficients_low_to_high": polynomial_coefficients(q_integral),
        "X_coefficients_low_to_high": polynomial_coefficients(X_integral_q),
        "Y_coefficients_low_to_high": polynomial_coefficients(Y_integral_q),
    },
    "height_on_quadratic_cover": {
        "value": 12,
        "source": "certified bisection lift intersections",
        "formula": "2*(P.sigma(P)-P.P)=2*(2-(-4))=12",
    },
    "modular_identification": modular_identification,
    "proof_boundary": (
        "This certifies the already-known non-torsion section of E^q(Q(t)); "
        "it does not bound the full twist rank."
    ),
    "inputs": {
        relative(args.bisections): digest(args.bisections),
        relative(args.model): digest(args.model),
    },
}
if args.output is None:
    args.output = (
        ROOT
        / "artifacts/generated-results"
        / (
            f"elkies-k3-r17-norm12-direct-singleton-twist-section-{args.direct_label}.json"
            if direct
            else f"elkies-2026-singleton-twist-section-mask-{args.mask}.json"
        )
    )
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"ELKIES2026TWISTSECTION|mask={mask}|degrees={X.degree()},{Y.degree()}"
    f"|modular_block={None if modular_identification is None else modular_identification['section_block']}"
    f"|output={args.output}|status=PASS_EXACT",
    flush=True,
)
