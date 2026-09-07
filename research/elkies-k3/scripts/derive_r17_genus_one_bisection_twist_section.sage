#!/usr/bin/env sage-python
"""Descend a certified genus-one bisection lift to its quartic twist.

For a construction record with ``s^2=q(t)`` and lifted point ``P``, this
forms ``T=P-sigma(P)``.  Its x-coordinate is invariant and its y-coordinate
is ``s`` times a rational function.  The corresponding point on

    Y^2 = X^3 + A*q^2*X + B*q^3

is ``(X,Y)=(q*x(T),q^2*(y(T)/s))``.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import lcm
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, matrix, sage_eval


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONSTRUCTIONS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
)
DEFAULT_MODEL = (
    ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def polynomial(ring, values):
    return ring([QQ(value) for value in values])


def p_valuation(value: int, prime: int) -> int:
    if value == 0:
        return 10**9
    result = 0
    value = abs(int(value))
    while value % prime == 0:
        result += 1
        value //= prime
    return result


def coefficient_strings(value):
    return [str(coefficient) for coefficient in value.list()]


def rational_function_record(value):
    value = value.numerator() / value.denominator()
    return {
        "numerator_coefficients_low_to_high": coefficient_strings(value.numerator()),
        "denominator_coefficients_low_to_high": coefficient_strings(value.denominator()),
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--label", default="norm12-orbit-103b2")
parser.add_argument("--constructions", type=Path, default=DEFAULT_CONSTRUCTIONS)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--prime", type=int)
parser.add_argument("--export", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

construction_path = args.constructions.resolve()
model_path = args.model.resolve()
constructions = json.loads(construction_path.read_text())
records = constructions["construction"]["records"]
record = next((item for item in records if item["label"] == args.label), None)
if record is None:
    raise ValueError(f"unknown construction label {args.label}")
model = json.loads(model_path.read_text())
if model.get("status") != "PASS_TRANSCRIBED_PUBLISHED_R17_MODEL":
    raise ValueError("expected the certified compact published R17 model")

ring = PolynomialRing(QQ, "t")
fraction_field = ring.fraction_field()
A = polynomial(ring, model["A_coefficients_low_to_high"])
B = polynomial(ring, model["B_coefficients_low_to_high"])
q = polynomial(ring, record["branch_polynomial_q_coefficients_low_to_high"])

extension_ring = PolynomialRing(fraction_field, "z")
z = extension_ring.gen()
cover_field = fraction_field.extension(z**2 - q, "s")
s = cover_field.gen()
curve = EllipticCurve(cover_field, [A, B])
lift = record["lifted_section"]
x0 = polynomial(ring, lift["x0_coefficients_low_to_high"])
x1 = polynomial(ring, lift["x1_coefficients_low_to_high"])
y0 = polynomial(ring, lift["y0_coefficients_low_to_high"])
y1 = polynomial(ring, lift["y1_coefficients_low_to_high"])
P = curve(cover_field(x0) + cover_field(x1) * s, cover_field(y0) + cover_field(y1) * s)
sigma_P = curve(
    cover_field(x0) - cover_field(x1) * s,
    cover_field(y0) - cover_field(y1) * s,
)
anti = P - sigma_P
x_coefficients = anti[0].list()
y_coefficients = anti[1].list()
if len(x_coefficients) > 1 and any(x_coefficients[1:]):
    raise ArithmeticError("x(P-sigma(P)) is not invariant")
if not y_coefficients or y_coefficients[0] != 0 or any(y_coefficients[2:]):
    raise ArithmeticError("y(P-sigma(P)) is not purely anti-invariant")

X = fraction_field(q) * fraction_field(x_coefficients[0])
Y = fraction_field(q) ** 2 * fraction_field(y_coefficients[1])
if Y**2 != X**3 + fraction_field(A * q**2) * X + fraction_field(B * q**3):
    raise ArithmeticError("descended twist section fails the exact curve identity")

modular_identification = None
if args.prime is not None:
    if X.denominator().degree() or Y.denominator().degree():
        raise ArithmeticError("modular P.O=0 identification requires polynomial coordinates")
    if args.export is None:
        args.export = (
            ROOT
            / "artifacts/local/elkies-k3/twist-polynomial-sections"
            / f"genus-one-{args.label}/p{args.prime}/export.json"
        )
    export_path = args.export.resolve()
    export = json.loads(export_path.read_text())
    if export["candidate"]["kind"] != "genus_one" or export["candidate"]["key"] != args.label:
        raise ValueError("modular export does not match the requested construction")
    prime = int(args.prime)
    if int(export["prime"]) != prime:
        raise ValueError("modular export prime mismatch")
    field = GF(prime)
    modular_ring = PolynomialRing(field, "u")
    u = modular_ring.gen()
    q_denominator = lcm(*(coefficient.denominator() for coefficient in q.list()))
    integral_q = [int(q_denominator**2 * coefficient) for coefficient in q.list()]
    common_valuation = min(p_valuation(value, prime) for value in integral_q)
    if common_valuation % 2:
        raise ArithmeticError("integral twist model has odd p-adic content")
    # The exporter removes common even p-adic content from D^2*q.  If its
    # model is q' = r^2*q, then the corresponding coordinates are
    # (X',Y')=(r^2*X,r^3*Y).  Using D directly would spuriously reduce the
    # section to zero at primes dividing D.
    twist_scale = QQ(q_denominator) / prime ** (common_valuation // 2)
    polynomial_X = ring(twist_scale**2 * X)
    polynomial_Y = ring(twist_scale**3 * Y)
    chart = export["infinity_fibre"]["chart"]
    if chart == "original_infinity":
        transformed_X = modular_ring(polynomial_X)
        transformed_Y = modular_ring(polynomial_Y)
    else:
        prefix = "original_t="
        if not chart.startswith(prefix) or " via " not in chart:
            raise ValueError(f"unrecognized exporter chart {chart}")
        chart_parameter = field(int(chart[len(prefix) :].split(" via ", 1)[0]))
        transformed_X = modular_ring(
            sum(
                field(polynomial_X[index])
                * (chart_parameter * u + 1) ** index
                * u ** (8 - index)
                for index in range(polynomial_X.degree() + 1)
            )
        )
        transformed_Y = modular_ring(
            sum(
                field(polynomial_Y[index])
                * (chart_parameter * u + 1) ** index
                * u ** (12 - index)
                for index in range(polynomial_Y.degree() + 1)
            )
        )
    matches = [
        int(system["block_index"])
        for system in export["systems"]
        if system["leading_x_y"] == [int(transformed_X[8]), int(transformed_Y[12])]
    ]
    if len(matches) != 1:
        raise ArithmeticError("known section did not identify a unique modular block")
    system_record = export["systems"][matches[0]]
    system_path = ROOT / system_record["path"]
    lines = system_path.read_text().splitlines()
    names = tuple(lines[0].split(","))
    coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
    local_variables = coefficient_ring.gens_dict()
    substitutions = {
        local_variables[f"x{index}"]: transformed_X[index] for index in range(8)
    }
    equations = [
        coefficient_ring(sage_eval(piece.replace("^", "**"), locals=local_variables))
        for piece in "\n".join(lines[2:]).split(",\n")
    ]
    if len(equations) != 12 or any(equation.subs(substitutions) for equation in equations):
        raise ArithmeticError("known section does not solve its exported modular system")
    jacobian = matrix(
        coefficient_ring,
        [[equation.derivative(variable) for variable in coefficient_ring.gens()] for equation in equations],
    )
    tangent_rank = jacobian.subs(substitutions).rank()
    modular_identification = {
        "prime": prime,
        "export": relative(export_path),
        "export_sha256": digest(export_path),
        "section_block": matches[0],
        "leading_x_y": [int(transformed_X[8]), int(transformed_Y[12])],
        "section_X_coefficients_low_to_high": [int(transformed_X[index]) for index in range(9)],
        "section_Y_coefficients_low_to_high": [int(transformed_Y[index]) for index in range(13)],
        "exported_system": relative(system_path),
        "exported_system_sha256": digest(system_path),
        "exported_equation_residuals_zero": True,
        "jacobian_rank_at_known_section": int(tangent_rank),
        "known_section_is_reduced_isolated_point": tangent_rank == len(names),
    }

payload = {
    "schema": "elkies-k3.r17-genus-one-bisection-twist-section.v1",
    "status": "PASS_EXACT_DESCENDED_TWIST_SECTION",
    "label": args.label,
    "orbit_hex": f"0x{int(record['lattice_orbit_mask']):05x}",
    "construction": "T=P-sigma(P); (X,Y)=(q*x(T),q^2*coefficient_s(y(T)))",
    "twist_model": "Y^2=X^3+A*q^2*X+B*q^3",
    "q_coefficients_low_to_high": coefficient_strings(q),
    "X": rational_function_record(X),
    "Y": rational_function_record(Y),
    "degrees": {
        "X_numerator_denominator": [int(X.numerator().degree()), int(X.denominator().degree())],
        "Y_numerator_denominator": [int(Y.numerator().degree()), int(Y.denominator().degree())],
    },
    "height_on_twist": {
        "value": 8,
        "source": "the certified cover height is 16 and heights double under the degree-two base change",
    },
    "modular_identification": modular_identification,
    "proof_boundary": "This certifies one non-torsion section; it does not bound the full twist rank.",
    "inputs": {
        relative(construction_path): digest(construction_path),
        relative(model_path): digest(model_path),
    },
}
if args.output is None:
    args.output = (
        ROOT
        / "artifacts/generated-results"
        / f"elkies-k3-{args.label}-twist-section-v1.json"
    )
output_path = args.output if args.output.is_absolute() else ROOT / args.output
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"R17GENUS1TWIST|label={args.label}|degrees={payload['degrees']}"
    f"|output={relative(output_path)}|status=PASS_EXACT",
    flush=True,
)
