#!/usr/bin/env sage-python
"""Export the full monic-quadratic trace-zero P.O=0 R17 schemes modulo p.

On the native 074d9 short model, search

    q=u^2+q1*u+q0,
    Y^2=X^3+q^2*A*X+q^3*B,
    deg(X)<=6, deg(Y)<=9.

Monicity gives rational points above infinity on the quadratic base. For each
leading coefficient x6 and each square root y9 of the leading equation, the
remaining coefficient equations keep y0,...,y8 as variables and are saturated
by the discriminant of q. The union of the exported blocks is the complete
finite monic-q P.O=0 chart at the selected prime.
"""

# status: ACTIVE_COMPILER
# claim: export the complete finite monic-q polynomial P.O=0 coefficient schemes
# inputs: artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json
# outputs: caller-selected local msolve systems and compact export manifest

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MODEL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    if not ZZ(args.prime).is_prime() or args.prime < 5:
        parser.error("--prime must be a prime at least five")

    model_path = args.model.resolve()
    document = json.loads(model_path.read_text())
    if document.get("status") != "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
        raise ValueError("expected the exact R17 lineage certificate")
    equation = document["representative"]
    if equation.get("chart") != "norm12-orbit-074d9":
        raise ValueError("the exact lineage representative changed")

    prime = int(args.prime)
    field = GF(prime)

    def coefficient_mod(value):
        value = QQ(value)
        denominator = int(value.denominator()) % prime
        if denominator == 0:
            raise ZeroDivisionError(f"model coefficient denominator vanishes at p={prime}")
        return field(int(value.numerator()) % prime) / field(denominator)

    names = (
        ("q0", "q1")
        + tuple(f"x{index}" for index in range(6))
        + tuple(f"y{index}" for index in range(9))
        + ("v",)
    )
    coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
    variables = coefficient_ring.gens_dict()
    t_ring = PolynomialRing(coefficient_ring, "u")
    u = t_ring.gen()
    A = t_ring([coefficient_mod(value) for value in equation["A_coefficients_low_to_high"]])
    B = t_ring([coefficient_mod(value) for value in equation["B_coefficients_low_to_high"]])
    if A.degree() != 8 or B.degree() != 12:
        raise ValueError("the native model degree drops at the selected prime")

    q0, q1 = variables["q0"], variables["q1"]
    q = u**2 + q1 * u + q0
    x_tail = sum(variables[f"x{index}"] * u**index for index in range(6))
    discriminant = q1**2 - 4 * q0
    output_dir = args.output_dir if args.output_dir.is_absolute() else ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    blocks = []
    for x6 in field:
        leading = x6**3 + coefficient_mod(equation["A_coefficients_low_to_high"][8]) * x6 + coefficient_mod(equation["B_coefficients_low_to_high"][12])
        # Include leading=0: it is the legitimate degree-drop block y9=0.
        if not leading.is_square():
            continue
        root = leading.sqrt()
        y9_values = sorted({int(root), int(-root)})
        for y9_integer in y9_values:
            y9 = field(y9_integer)
            X = x_tail + x6 * u**6
            right = X**3 + q**2 * A * X + q**3 * B
            Y = y9 * u**9 + sum(
                variables[f"y{index}"] * u**index for index in range(9)
            )
            residual = right - Y**2
            equations = [coefficient_ring(residual[index]) for index in range(18)]
            equations.append(variables["v"] * discriminant - 1)
            if residual[18]:
                raise ArithmeticError("the selected leading block missed degree 18")
            text = ",".join(names) + f"\n{prime}\n"
            text += ",\n".join(str(value) for value in equations) + "\n"
            path = output_dir / f"x6-{int(x6):03d}-y9-{y9_integer:03d}.ms"
            path.write_text(text)
            blocks.append(
                {
                    "block_index": len(blocks),
                    "leading_x6_y9": [int(x6), y9_integer],
                    "path": relative(path),
                    "sha256": sha256(text.encode()).hexdigest(),
                }
            )

    summary = args.summary if args.summary.is_absolute() else ROOT / args.summary
    summary.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "elkies-k3.r17-rational-quadratic-tracezero-po0-msolve-export.v1",
        "status": "PASS_COMPLETE_MODP_MONIC_QUADRATIC_PO0_SCHEMES_EXPORTED",
        "prime": prime,
        "coordinate": "norm12-orbit-074d9 native u",
        "family": "q=u^2+q1*u+q0",
        "variable_names": list(names),
        "equations_per_block": 19,
        "block_count": len(blocks),
        "blocks": blocks,
        "inputs": {relative(model_path): digest(model_path)},
        "proof_boundary": (
            "The blocks exhaust monic squarefree quadratic q and polynomial P.O=0 "
            "short-twist sections with deg(X)<=6 and deg(Y)<=9 modulo the selected "
            "good prime. Solving, extracting rational finite-field points, intersecting "
            "across primes, characteristic-zero lifting, and height independence are "
            "separate gates. Nonmonic rational covers whose rational point is not above "
            "infinity require an anchored coordinate chart."
        ),
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "export_r17_rational_quadratic_tracezero_po0_msolve.sage "
            f"--model {relative(model_path)} --prime {prime} "
            f"--output-dir {relative(output_dir)} --summary {relative(summary)}"
        ),
    }
    summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "R17MONICQPO0EXPORT"
        f"|p={prime}|blocks={len(blocks)}|variables={len(names)}"
        f"|equations=19|output={relative(summary)}|status={payload['status']}"
    )


if __name__ == "__main__":
    main()
