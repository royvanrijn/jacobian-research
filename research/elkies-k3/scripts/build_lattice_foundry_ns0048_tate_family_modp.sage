#!/usr/bin/env sage-python
"""Build the section-translated NS0048 Tate-family system modulo p."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, binomial, matrix


ROOT = Path(__file__).resolve().parents[2]
FIBRES = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-source-ansatz-mod11-suffix600k-v1.json"
SECTIONS = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-pole0-sections-xonly-mod11-suffix600k-v1.json"
DEFAULT_OUTPUT_DIR = ROOT / "artifacts/local/elkies-k3/ns0048-tate-family-modp"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=11)
parser.add_argument("--fixed-lambda", type=int)
parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

field = GF(args.prime)
fixed_lambda = None if args.fixed_lambda is None else field(args.fixed_lambda)
names = [f"a2_{index}" for index in range(4, -1, -1)]
names += [f"r_{index}" for index in range(4, -1, -1)]
names += [f"s_{index}" for index in range(6, -1, -1)]
if fixed_lambda is None:
    names += ["lambda"]
coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
g = coefficient_ring.gens_dict()
outer = PolynomialRing(coefficient_ring, "t")
t = outer.gen()
lambda_value = coefficient_ring(fixed_lambda) if fixed_lambda is not None else g["lambda"]
support_product = (t - 1) * (t - lambda_value)
a2 = outer([g[f"a2_{index}"] for index in range(5)])
r = outer([g[f"r_{index}"] for index in range(5)])
s = outer([g[f"s_{index}"] for index in range(7)])
a3 = support_product * r
a4 = support_product * s

short_A = a4 - a2**2 / 3
short_B = a3**2 / 4 + 2 * a2**3 / 27 - a2 * a4 / 3
b2 = 4 * a2
b4 = 2 * a4
b6 = a3**2
b8 = a2 * a3**2 - a4**2
discriminant = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
# Both a3 and a4 are divisible by the two marked multiplicative supports.
# Removing this forced square before taking jets avoids expanding thousands of
# terms which cancel identically.  On the intended open locus lambda != 0,1,
# order five of this quotient at zero is equivalent to I5, while order five at
# one is equivalent to the five extra conditions upgrading the forced I2 to
# I7.
discriminant_quotient = (
    -16 * a2**3 * r**2
    + 16 * a2**2 * s**2
    - 64 * support_product * s**3
    - 27 * support_product**2 * r**4
    + 72 * a2 * support_product * s * r**2
)
if discriminant != support_product**2 * discriminant_quotient:
    raise ArithmeticError("Tate discriminant quotient identity failed")


def truncate(poly, precision):
    return outer(poly.list()[:precision])


def multiply_truncated(left, right, precision):
    return truncate(left * right, precision)


def power_truncated(poly, exponent, precision):
    answer = outer.one()
    for unused in range(exponent):
        answer = multiply_truncated(answer, poly, precision)
    return answer


def discriminant_truncated(local_a2, local_a3, local_a4, precision):
    local_b2 = truncate(4 * local_a2, precision)
    local_b4 = truncate(2 * local_a4, precision)
    local_b6 = multiply_truncated(local_a3, local_a3, precision)
    local_b8 = truncate(
        multiply_truncated(local_a2, local_b6, precision)
        - multiply_truncated(local_a4, local_a4, precision),
        precision,
    )
    return truncate(
        -multiply_truncated(
            multiply_truncated(local_b2, local_b2, precision),
            local_b8,
            precision,
        )
        - 8 * power_truncated(local_b4, 3, precision)
        - 27 * multiply_truncated(local_b6, local_b6, precision)
        + 9
        * multiply_truncated(
            multiply_truncated(local_b2, local_b4, precision),
            local_b6,
            precision,
        ),
        precision,
    )


def shifted_truncated(poly, point, precision):
    return outer(
        [
            sum(
                coefficient_ring(poly[index])
                * coefficient_ring(binomial(index, jet))
                * coefficient_ring(point) ** (index - jet)
                for index in range(jet, len(poly.list()))
            )
            for jet in range(precision)
        ]
    )


delta_quotient_zero = truncate(discriminant_quotient, 5)
delta_quotient_one = shifted_truncated(
    discriminant_quotient, coefficient_ring.one(), 5
)

equations = [coefficient_ring(short_A[0] + 3)]
equations += [coefficient_ring(short_A[index]) for index in (7, 8)]
equations += [coefficient_ring(short_B[index]) for index in range(10, 13)]
equations += [coefficient_ring(delta_quotient_zero[index]) for index in range(5)]
equations += [coefficient_ring(delta_quotient_one[index]) for index in range(5)]
equations += [coefficient_ring(discriminant_quotient[14])]
if len(equations) != 17 or any(not equation for equation in equations):
    raise ArithmeticError("unexpected reduced Tate-family equation accounting")

modular_point = None
jacobian_rank = None
if args.prime == 11 and fixed_lambda in (None, field(10)):
    values_by_name = {
        **{f"a2_{index}": value for index, value in enumerate([9, 2, 0, 4, 1])},
        **{f"r_{index}": value for index, value in enumerate([8, 5, 5, 3, 6])},
        **{f"s_{index}": value for index, value in enumerate([9, 9, 2, 9, 10, 10, 4])},
        "lambda": 10,
    }
    modular_point = [field(values_by_name[name]) for name in names]
    if any(equation(*modular_point) for equation in equations):
        raise ArithmeticError("pinned marked model does not satisfy reduced Tate system")
    jacobian = matrix(
        field,
        [
            [equation.derivative(variable)(*modular_point) for variable in coefficient_ring.gens()]
            for equation in equations
        ],
    )
    jacobian_rank = int(jacobian.rank())
    expected_rank = 17
    if jacobian_rank != expected_rank:
        raise ArithmeticError("unexpected reduced Tate-family Jacobian rank")

output_dir = args.output_dir.resolve()
output_dir.mkdir(parents=True, exist_ok=True)
stem = f"p{args.prime}" + (f"-lambda{int(fixed_lambda)}" if fixed_lambda is not None else "-curve")
msolve_path = output_dir / f"{stem}.ms"
metadata_path = output_dir / f"{stem}.json"
msolve_text = ",".join(names) + "\n" + str(args.prime) + "\n"
msolve_text += ",\n".join(str(equation).replace("**", "^") for equation in equations) + "\n"

metadata = {
    "schema": "elkies-k3.lattice-foundry-ns0048-tate-family-modp-system.v1",
    "status": "PASS_EXACT_REDUCED_TATE_FAMILY_SYSTEM",
    "prime": args.prime,
    "fixed_lambda": None if fixed_lambda is None else int(fixed_lambda),
    "inputs": {
        "fibre_artifact": relative(FIBRES),
        "fibre_artifact_sha256": hashlib.sha256(FIBRES.read_bytes()).hexdigest(),
        "section_artifact": relative(SECTIONS),
        "section_artifact_sha256": hashlib.sha256(SECTIONS.read_bytes()).hexdigest(),
    },
    "coordinates": {
        "equation": "y^2+a3*y=x^3+a2*x^2+a4*x",
        "a3": "(t-1)(t-lambda) r(t), deg r<=4",
        "a4": "(t-1)(t-lambda) s(t), deg s<=6",
        "short_A": "a4-a2^2/3",
        "short_B": "a3^2/4+2a2^3/27-a2*a4/3",
        "discriminant_quotient": "Delta/((t-1)(t-lambda))^2",
    },
    "system": {
        "variables": names,
        "variable_count": len(names),
        "equation_count": len(equations),
        "equation_total_degrees": [int(equation.degree()) for equation in equations],
        "equation_term_counts": [len(equation.monomials()) for equation in equations],
        "msolve_input": relative(msolve_path),
        "msolve_sha256": hashlib.sha256(msolve_text.encode()).hexdigest(),
    },
    "pinned_point": (
        None
        if modular_point is None
        else {
            "values_in_variable_order": [int(value) for value in modular_point],
            "jacobian_rank": jacobian_rank,
            "tangent_dimension": len(names) - jacobian_rank,
        }
    ),
    "proof_boundary": (
        "The reduced system exactly imposes the section-translated Tate form, "
        "short-model degree bounds, I5 order at zero, the five additional I7 "
        "conditions at one, and I1* leading cancellation, after exactly removing "
        "the forced square of the two marked multiplicative supports. Exact open fibre "
        "orders, squarefree residual, characteristic-zero lift, and corridor "
        "still require separate verification."
    ),
}
metadata_text = json.dumps(metadata, indent=2, sort_keys=True) + "\n"
if args.check:
    if msolve_path.read_text() != msolve_text or metadata_path.read_text() != metadata_text:
        raise SystemExit("NS0048 Tate-family system is stale")
else:
    msolve_path.write_text(msolve_text)
    metadata_path.write_text(metadata_text)

print(
    "FOUNDRYNS0048TATE|"
    f"p={args.prime}|fixed_lambda={'' if fixed_lambda is None else int(fixed_lambda)}|"
    f"vars={len(names)}|eqs={len(equations)}|rank={'' if jacobian_rank is None else jacobian_rank}|status=PASS",
    flush=True,
)
print(f"MSOLVE_INPUT|{msolve_path}", flush=True)
print(f"OUTPUT|{metadata_path}", flush=True)
