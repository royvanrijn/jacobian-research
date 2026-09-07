#!/usr/bin/env sage-python
"""Exhaust a finite-field two-parameter genus-2 normalization export.

For every exported norm-six trace and every ``(l0,l1)`` over the displayed
prime field, factor the exact branch sextic.  Retain nonsplit rational
normalizations, namely specializations whose squarefree part has degree two.
This is complete for the finite affine parameter chart, but it is only a
modular discovery sieve for characteristic-zero pencils.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def relative(path):
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def mod_p(value, field):
    value = QQ(value)
    if int(value.denominator()) % field.characteristic() == 0:
        raise ZeroDivisionError("prime divides a rational coefficient denominator")
    return field(value.numerator()) / field(value.denominator())


def evaluate_mod_p(polynomial, l0_value, l1_value, field):
    return sum(
        mod_p(coefficient, field)
        * l0_value**exponents[0]
        * l1_value**exponents[1]
        for exponents, coefficient in polynomial.dict().items()
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--export", type=Path, required=True)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

export_path = args.export.resolve()
export = json.loads(export_path.read_text())
if export.get("schema") != "elkies-k3.r17-norm12-direct-genus2-normalization-msolve-export.v1":
    raise ValueError("unexpected normalization export schema")
prime = int(export["prime"])
field = GF(prime)
parameter_ring = PolynomialRing(QQ, names=("l0", "l1"))
l0, l1 = parameter_ring.gens()
u_ring = PolynomialRing(field, "u")

survivors = []
profile_histogram = Counter()
for system in export["systems"]:
    q_coefficients = [
        parameter_ring(text)
        for text in system["q_coefficients_in_QQ_l0_l1_low_to_high"]
    ]
    for l0_value in field:
        for l1_value in field:
            coefficients = []
            for coefficient in q_coefficients:
                coefficients.append(
                    evaluate_mod_p(coefficient, l0_value, l1_value, field)
                )
            q = u_ring(coefficients)
            if not q:
                profile_histogram["zero_polynomial"] += 1
                continue
            factorization = q.factor()
            square_part = u_ring.one()
            reduced = u_ring(factorization.unit())
            for factor, exponent in factorization:
                square_part *= factor ** (int(exponent) // 2)
                if int(exponent) % 2:
                    reduced *= factor
            if square_part**2 * reduced != q:
                raise ArithmeticError("finite-field squareclass decomposition failed")
            profile = f"q{q.degree()}_square{square_part.degree()}_reduced{reduced.degree()}"
            profile_histogram[profile] += 1
            if reduced.degree() != 2:
                continue
            survivors.append(
                {
                    "trace_index": int(system["trace_index"]),
                    "translation_orbit_mask": int(system["translation_orbit_mask"]),
                    "basis_coordinates": system["basis_coordinates"],
                    "l0_l1": [int(l0_value), int(l1_value)],
                    "branch_coefficients_low_to_high": [int(value) for value in q],
                    "removed_square_factor_coefficients_low_to_high": [
                        int(value) for value in square_part
                    ],
                    "reduced_quadratic_coefficients_low_to_high": [
                        int(value) for value in reduced
                    ],
                    "factorization": [
                        {
                            "coefficients_low_to_high": [int(value) for value in factor],
                            "multiplicity": int(exponent),
                        }
                        for factor, exponent in factorization
                    ],
                }
            )

if args.output is None:
    output_path = export_path.with_name("factor-screen.json")
else:
    output_path = args.output if args.output.is_absolute() else ROOT / args.output
output_path.parent.mkdir(parents=True, exist_ok=True)
payload = {
    "schema": "elkies-k3.r17-norm12-direct-genus2-normalization-modp-screen.v1",
    "status": (
        "PASS_COMPLETE_MODP_NONSPLIT_RATIONAL_NORMALIZATION_SURVIVORS"
        if survivors
        else "PASS_COMPLETE_MODP_NO_NONSPLIT_RATIONAL_NORMALIZATION"
    ),
    "proof_boundary": (
        "Every affine (l0,l1) pair over the displayed prime field was factored "
        "for every exported finite-chart trace. Survivors have squarefree branch "
        "degree two modulo p. This is not a characteristic-zero existence or "
        "nonexistence result, and traces skipped by the exporter are out of scope."
    ),
    "prime": prime,
    "trace_count": len(export["systems"]),
    "parameter_pair_count_per_trace": prime**2,
    "total_specialization_count": len(export["systems"]) * prime**2,
    "survivor_count": len(survivors),
    "survivor_trace_count": len({item["trace_index"] for item in survivors}),
    "factor_profile_histogram": dict(sorted(profile_histogram.items())),
    "survivors": survivors,
    "inputs": {relative(export_path): digest(export_path)},
    "reproducing_command": (
        "sage -python "
        "elkies-k3/scripts/screen_r17_norm12_direct_genus2_normalizations_modp.sage "
        f"--export {relative(export_path)}"
    ),
}
output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"R17GENUS2SCREEN|p={prime}|traces={len(export['systems'])}"
    f"|specializations={payload['total_specialization_count']}"
    f"|survivors={len(survivors)}|survivor_traces={payload['survivor_trace_count']}"
    f"|output={relative(output_path)}|status={payload['status']}",
    flush=True,
)
