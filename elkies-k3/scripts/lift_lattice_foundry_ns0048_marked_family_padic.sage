#!/usr/bin/env sage-python
"""Hensel-lift the smooth GF(11) NS0048 marked source with a1 fixed."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, QQ, ZZ, PolynomialRing, binomial, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FIBRES = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-source-ansatz-mod11-suffix600k-v1.json"
DEFAULT_SECTIONS = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-pole0-sections-xonly-mod11-suffix600k-v1.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-marked-family-padic-lift-p11-v1.json"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def order_at(poly, point):
    if not poly:
        return None
    shifted = poly(poly.parent().gen() + point)
    return min(index for index, value in enumerate(shifted.list()) if value)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fibres", type=Path, default=DEFAULT_FIBRES)
parser.add_argument("--sections", type=Path, default=DEFAULT_SECTIONS)
parser.add_argument("--fibre-example", type=int, default=13)
parser.add_argument("--section-index", type=int, default=0)
parser.add_argument("--precision", type=int, default=80)
parser.add_argument("--fixed-parameter", choices=[f"a{index}" for index in range(1, 7)] + [f"x{index}" for index in range(5)] + [f"y{index}" for index in range(7)] + ["lambda"], default="a1")
parser.add_argument("--fixed-value", type=int)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

fibre_path = args.fibres.resolve()
section_path = args.sections.resolve()
fibres = json.loads(fibre_path.read_text())
sections = json.loads(section_path.read_text())
if sections["input"]["artifact"] != relative(fibre_path):
    raise ValueError("section artifact does not consume the supplied fibres")
prime = ZZ(fibres["prime"])
if prime != 11 or sections["prime"] != 11:
    raise ValueError("the pinned lift is the good-prime p=11 chart")
model = next(row for row in sections["models"] if row["example_index"] == args.fibre_example)
section = model["solutions"][args.section_index]
fibre = fibres["examples"][args.fibre_example]

names = [f"a{index}" for index in range(1, 7)]
names += [f"x{index}" for index in range(5)]
names += [f"y{index}" for index in range(7)]
names += ["lambda"]
coefficient_ring = PolynomialRing(ZZ, names=names, order="degrevlex")
g = coefficient_ring.gens_dict()
outer = PolynomialRing(coefficient_ring, "t")
t = outer.gen()
A = outer([-3] + [g[f"a{index}"] for index in range(1, 7)])
X = outer([g[f"x{index}"] for index in range(5)])
Y = outer([g[f"y{index}"] for index in range(7)])
B = Y**2 - X**3 - A * X
discriminant_core = 4 * A**3 + 27 * B**2
lambda_symbol = g["lambda"]
fibre_equations = [coefficient_ring(B[index]) for index in range(10, 13)]
fibre_equations += [coefficient_ring(discriminant_core[index]) for index in range(5)]
fibre_equations += [
    sum(
        coefficient_ring(discriminant_core[index])
        * coefficient_ring(binomial(index, jet))
        for index in range(jet, 19)
    )
    for jet in range(7)
]
fibre_equations += [
    sum(
        coefficient_ring(discriminant_core[index])
        * coefficient_ring(binomial(index, jet))
        * lambda_symbol ** (index - jet)
        for index in range(jet, 19)
    )
    for jet in range(2)
]
fibre_equations += [coefficient_ring(discriminant_core[18])]
component_equations = [
    coefficient_ring(Y(1)),
    coefficient_ring(3 * X(1) ** 2 + A(1)),
    coefficient_ring(Y(lambda_symbol)),
    coefficient_ring(3 * X(lambda_symbol) ** 2 + A(lambda_symbol)),
]
equations = fibre_equations + component_equations
if len(fibre_equations) != 18 or len(equations) != 22:
    raise ArithmeticError("unexpected NS0048 family-system accounting")

initial_values = [ZZ(value) for value in fibre["A_coefficients_low_to_high"][1:7]]
initial_values += [ZZ(value) for value in section["X_coefficients_low_to_high"]]
initial_values += [ZZ(value) for value in section["Y_coefficients_low_to_high"]]
initial_values += [ZZ(fibre["lambda"])]
fixed_index = names.index(args.fixed_parameter)
fixed_residue = initial_values[fixed_index] % prime
fixed_value = (
    ZZ(args.fixed_value)
    if args.fixed_value is not None
    else ZZ(fixed_residue if fixed_residue <= prime // 2 else fixed_residue - prime)
)
if fixed_value % prime != fixed_residue:
    raise ValueError("--fixed-value must reduce to the pinned modular coordinate")
unknown_indices = [index for index in range(len(names)) if index != fixed_index]

field = GF(prime)
modular_values = [field(value) for value in initial_values]
jacobian_modp = matrix(
    field,
    [
        [
            field(equation.derivative(coefficient_ring.gen(index))(*initial_values))
            for index in unknown_indices
        ]
        for equation in equations
    ],
)
if jacobian_modp.rank() != 18:
    raise ArithmeticError("fixed-parameter Jacobian is not full column rank")

values = list(initial_values)
values[fixed_index] = fixed_value
modulus = prime
for unused_digit in range(1, args.precision):
    equation_values = [ZZ(equation(*values)) for equation in equations]
    if any(value % modulus for value in equation_values):
        raise ArithmeticError("current Hensel point lost its certified precision")
    right = vector(field, [field(-(value // modulus)) for value in equation_values])
    correction = jacobian_modp.solve_right(right)
    for index, delta in zip(unknown_indices, correction):
        values[index] += modulus * ZZ(int(delta))
    values[fixed_index] = fixed_value
    modulus *= prime
    if any(ZZ(equation(*values)) % modulus for equation in equations):
        raise ArithmeticError("overdetermined Hensel correction failed")

reconstructed = []
reconstruction_error = None
for value in values:
    try:
        reconstructed.append(QQ(ZZ(value % modulus).rational_reconstruction(modulus)))
    except ArithmeticError as error:
        reconstruction_error = str(error)
        reconstructed = []
        break

exact_equations_zero = bool(reconstructed) and all(
    QQ(equation(*reconstructed)) == 0 for equation in equations
)
exact_model = None
if exact_equations_zero:
    rational_ring = PolynomialRing(QQ, "t")
    tr = rational_ring.gen()
    Ar = rational_ring([-3] + reconstructed[:6])
    Xr = rational_ring(reconstructed[6:11])
    Yr = rational_ring(reconstructed[11:18])
    lambdar = reconstructed[18]
    Br = Yr**2 - Xr**3 - Ar * Xr
    Dr = 4 * Ar**3 + 27 * Br**2
    orders = (
        order_at(Dr, QQ(0)),
        order_at(Dr, QQ(1)),
        order_at(Dr, lambdar),
        24 - Dr.degree(),
    )
    divisor = tr**5 * (tr - 1) ** 7 * (tr - lambdar) ** 2
    residual, remainder = Dr.quo_rem(divisor)
    if (
        Br.degree() > 9
        or orders != (5, 7, 2, 7)
        or remainder
        or residual.degree() != 3
        or residual.gcd(residual.derivative()).degree() != 0
    ):
        raise ArithmeticError("reconstructed rational model fails exact fibre audit")
    exact_model = {
        "A_coefficients_low_to_high": [str(value) for value in Ar],
        "B_coefficients_low_to_high": [str(value) for value in Br],
        "X_coefficients_low_to_high": [str(value) for value in Xr],
        "Y_coefficients_low_to_high": [str(value) for value in Yr],
        "lambda": str(lambdar),
        "discriminant_orders_I5_I7_I2_I1star": list(orders),
        "residual_discriminant_coefficients_low_to_high": [str(value) for value in residual],
    }

output = {
    "schema": "elkies-k3.lattice-foundry-ns0048-marked-family-padic-lift.v1",
    "status": (
        "PASS_EXACT_RATIONAL_MARKED_SOURCE_MODEL"
        if exact_model is not None
        else "PASS_PADIC_LIFT_NO_EXACT_RATIONAL_RECONSTRUCTION"
    ),
    "prime": int(prime),
    "precision": args.precision,
    "modulus": str(modulus),
    "fixed_parameter": {"name": args.fixed_parameter, "value": str(fixed_value)},
    "inputs": {
        "fibre_artifact": relative(fibre_path),
        "fibre_artifact_sha256": hashlib.sha256(fibre_path.read_bytes()).hexdigest(),
        "fibre_example": args.fibre_example,
        "section_artifact": relative(section_path),
        "section_artifact_sha256": hashlib.sha256(section_path.read_bytes()).hexdigest(),
        "section_index": args.section_index,
    },
    "system": {
        "variables": names,
        "equations": len(equations),
        "unknowns_after_fixing_parameter": len(unknown_indices),
        "fixed_parameter_jacobian_rank_mod_p": int(jacobian_modp.rank()),
    },
    "lifted_residues_in_variable_order": [str(value % modulus) for value in values],
    "rational_reconstruction": {
        "all_coordinates_reconstructed": bool(reconstructed),
        "exact_equations_zero": exact_equations_zero,
        "error": reconstruction_error,
        "values_in_variable_order": [str(value) for value in reconstructed],
    },
    "exact_model": exact_model,
    "proof_boundary": {
        "proved": (
            f"The displayed modular marked point lifts uniquely with {args.fixed_parameter} fixed "
            "through the stated p-adic precision against all 22 equations."
            + (
                " Rational reconstruction gives the exact audited model recorded here."
                if exact_model is not None
                else " No exact rational model is inferred from the failed reconstruction."
            )
        ),
        "not_proved": (
            "A p-adic point without exact reconstruction is not a rational source; "
            "a rational model still requires exact equations and fibre/component audits."
            if exact_model is None
            else "The physical elliptic-neighbour corridor to NS0048-F001 remains open."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/lift_lattice_foundry_ns0048_marked_family_padic.sage"
        + (f" --fixed-parameter {args.fixed_parameter}" if args.fixed_parameter != "a1" else "")
        + (f" --fixed-value {fixed_value}" if args.fixed_value is not None else "")
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
output_path = args.output.resolve()
if args.check:
    if output_path.read_text() != serialized:
        raise SystemExit("NS0048 p-adic lift artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "FOUNDRYNS0048PADICLIFT|"
    f"p={prime}|precision={args.precision}|rank={jacobian_modp.rank()}|"
    f"rational={int(exact_model is not None)}|status=PASS",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
