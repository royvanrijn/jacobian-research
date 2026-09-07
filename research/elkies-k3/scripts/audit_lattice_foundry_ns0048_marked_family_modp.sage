#!/usr/bin/env sage-python
"""Audit the marked NS0048 family Jacobian at modular section survivors."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, binomial, matrix


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FIBRES = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-source-ansatz-mod7.json"
DEFAULT_SECTIONS = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-pole0-sections-xonly-mod7.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0048-marked-family-jacobian-mod7.json"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fibres", type=Path, default=DEFAULT_FIBRES)
parser.add_argument("--sections", type=Path, default=DEFAULT_SECTIONS)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

fibre_path = args.fibres.resolve()
section_path = args.sections.resolve()
fibres = json.loads(fibre_path.read_text())
sections = json.loads(section_path.read_text())
if fibres["schema"] != "elkies-k3.lattice-foundry-ns0048-source-ansatz-modp.v1":
    raise ValueError("unexpected NS0048 fibre schema")
if sections["schema"] != "elkies-k3.lattice-foundry-ns0048-pole0-sections-xonly-modp-scan.v1":
    raise ValueError("unexpected NS0048 section schema")
if sections["input"]["artifact"] != relative(fibre_path):
    raise ValueError("section artifact does not consume the supplied fibre artifact")
if sections["input"]["sha256"] != hashlib.sha256(fibre_path.read_bytes()).hexdigest():
    raise ValueError("section artifact has a stale fibre hash")
prime = int(fibres["prime"])
if int(sections["prime"]) != prime:
    raise ValueError("prime mismatch")
field = GF(prime)
twist = field(sections["quadratic_twist"])

names = [f"a{index}" for index in range(1, 7)]
names += [f"x{index}" for index in range(5)]
names += [f"y{index}" for index in range(7)]
names += ["lambda"]
coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
g = coefficient_ring.gens_dict()
outer = PolynomialRing(coefficient_ring, "t")
t = outer.gen()

records = []
for model in sections["models"]:
    fibre = fibres["examples"][model["example_index"]]
    A_values = [twist**2 * field(value) for value in fibre["A_coefficients_low_to_high"]]
    B_values = [twist**3 * field(value) for value in fibre["B_coefficients_low_to_high"]]
    for section_index, section in enumerate(model["solutions"]):
        X_values = [field(value) for value in section["X_coefficients_low_to_high"]]
        Y_values = [field(value) for value in section["Y_coefficients_low_to_high"]]
        lambda_value = field(fibre["lambda"])
        values = A_values[1:7] + X_values + Y_values + [lambda_value]

        A = outer([coefficient_ring(A_values[0])] + [g[f"a{index}"] for index in range(1, 7)])
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
        if len(fibre_equations) != 18 or len(equations) != 22 or any(equation(*values) for equation in equations):
            raise ArithmeticError("marked point does not satisfy the family system")

        specialized_B = [field(coefficient(*values)) for coefficient in B.list()]
        specialized_B += [field.zero()] * (10 - len(specialized_B))
        if specialized_B[:10] != B_values:
            raise ArithmeticError("section-derived B does not match the fibre model")
        jacobian = matrix(
            field,
            [
                [equation.derivative(variable)(*values) for variable in coefficient_ring.gens()]
                for equation in equations
            ],
        )
        tangent_basis = jacobian.right_kernel().basis()
        records.append(
            {
                "fibre_example_index": model["example_index"],
                "section_index": section_index,
                "lambda": int(lambda_value),
                "equation_count": len(equations),
                "variable_count": len(names),
                "jacobian_rank": int(jacobian.rank()),
                "tangent_dimension": len(tangent_basis),
                "tangent_basis": [[int(value) for value in row] for row in tangent_basis],
                "smooth_one_dimensional_at_point": jacobian.rank() == 18,
            }
        )

if not records:
    raise ValueError("section artifact contains no marked survivors")
all_smooth = all(record["smooth_one_dimensional_at_point"] for record in records)
output = {
    "schema": "elkies-k3.lattice-foundry-ns0048-marked-family-jacobian-modp.v1",
    "status": (
        "PASS_SMOOTH_ONE_DIMENSIONAL_MARKED_MODULAR_LOCUS"
        if all_smooth
        else "PASS_MARKED_MODULAR_POINTS_WITH_EXCESS_TANGENT_DIMENSION"
    ),
    "prime": prime,
    "inputs": {
        "fibre_artifact": relative(fibre_path),
        "fibre_artifact_sha256": hashlib.sha256(fibre_path.read_bytes()).hexdigest(),
        "section_artifact": relative(section_path),
        "section_artifact_sha256": hashlib.sha256(section_path.read_bytes()).hexdigest(),
    },
    "family_system": {
        "variables": names,
        "variable_count": len(names),
        "fibre_equation_count": 18,
        "component_incidence_equation_count": 4,
        "equation_count": 22,
        "expected_jacobian_rank": 18,
        "definition": (
            "B=Y^2-X^3-A*X; B_10=B_11=B_12=0; discriminant orders at least "
            "5,7,2 at 0,1,lambda; leading degree-18 discriminant coefficient zero; "
            "Y=0 and 3X^2+A=0 at the I7 and I2 supports."
        ),
    },
    "records": records,
    "proof_boundary": {
        "proved": (
            "The displayed modular marked points satisfy the 22-equation "
            "section-built family system, and their exact Jacobian ranks and "
            "Zariski tangent dimensions are recorded."
        ),
        "not_proved": (
            "Smoothness at one modular point does not by itself give a rational "
            "characteristic-zero family, a lift, or a neighbour corridor."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/audit_lattice_foundry_ns0048_marked_family_modp.sage"
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
output_path = args.output.resolve()
if args.check:
    if output_path.read_text() != serialized:
        raise SystemExit("NS0048 marked-family Jacobian artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "FOUNDRYNS0048MARKEDJACOBIAN|"
    f"p={prime}|points={len(records)}|ranks={','.join(str(row['jacobian_rank']) for row in records)}|"
    f"status={'SMOOTH' if all_smooth else 'EXCESS_TANGENT'}",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
