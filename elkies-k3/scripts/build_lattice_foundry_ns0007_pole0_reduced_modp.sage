#!/usr/bin/env sage-python
"""Build the globally reduced NS0007 pole-zero source system modulo p.

The split-node parametrization at infinity determines all coefficients of
``a3`` and the top seven coefficients of ``a4``.  The marked-section node
conditions at zero and one determine the two remaining ``a4`` coefficients.
This eliminates the global ``a3,a4`` variables from the section-first chart.

For fixed cross-ratio lambda the exact-I2 chart has 19 variables and 19
equations: four free nonconstant coefficients of a2, two seven-jets for the
split I7 nodes, and one Rabinowitsch inverse.  Solutions still require an
independent exact-fibre and squarefree-residual audit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ, is_prime


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-c-v1.json"
POLES = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-rank1-section-poles-v1.json"
DEFAULT_OUTPUT_DIR = ROOT / "artifacts/local/elkies-k3/ns0007-pole0-reduced-modp"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=7)
parser.add_argument("--lambda-value", type=int, default=2)
parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
parser.add_argument(
    "--split-a2-node-constants-over-base-field",
    action="store_true",
    help=(
        "append x^p-x for the four a2 coefficients and the two constant "
        "I7 node jets, encoding all base-field slices in one ideal"
    ),
)
parser.add_argument(
    "--fixed-a2-node-case",
    help=(
        "six comma-separated base-field values for "
        "a2_4,a2_3,a2_2,a2_1,si_0,sl_0; append exact linear equations"
    ),
)
parser.add_argument(
    "--compact-factored-msolve",
    action="store_true",
    help=(
        "retain the final I4 and exact-I2 equations in factored input syntax "
        "instead of expanding them into tens of thousands of monomials"
    ),
)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

prime = ZZ(args.prime)
if not is_prime(prime) or prime in (2, 3):
    raise SystemExit("--prime must be a prime other than 2 or 3")
field = GF(prime)
lambda_value = field(args.lambda_value)
if lambda_value in (0, 1):
    raise SystemExit("--lambda-value must differ from zero and one")

source_payload = json.loads(SOURCE.read_text())
source = next(
    row["source"]
    for row in source_payload["sources"]
    if row["ns_id"] == "NS0007" and row["source_id"] == "NS0007-S025"
)
assert source["root_type"] == "A1+A3+2A6"
assert source["root_lattice_primitive"] and source["torsion"] == 1
assert source["mw_height_gram"] == [["11/4"]]
pole_payload = json.loads(POLES.read_text())
pole_row = next(
    row
    for row in pole_payload["sources"]
    if row["source_artifact"] == relative(SOURCE)
    and row["source_id"] == "NS0007-S025"
)
assert pole_row["minimum_section_pole_order"] == 0

names = [f"sl_{index}" for index in range(6, -1, -1)]
names += [f"si_{index}" for index in range(6, -1, -1)]
names += [f"a2_{index}" for index in range(4, 0, -1)]
names += ["u0"]
coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
g = coefficient_ring.gens_dict()
polynomial_ring = PolynomialRing(coefficient_ring, "t")
t = polynomial_ring.gen()

a2 = coefficient_ring(3) + sum(g[f"a2_{index}"] * t**index for index in range(1, 5))
a2_infinity = polynomial_ring([a2[4 - index] for index in range(5)])
si = sum(g[f"si_{index}"] * t**index for index in range(7))
sl = sum(g[f"sl_{index}"] * t**index for index in range(7))

# Work modulo t^7 in the infinity coordinate.  Reversing the resulting jets
# recovers the global coefficients in degrees 6..0 and 8..2.
a3_infinity = polynomial_ring(si * (a2_infinity - si**2) % t**7)
a4_infinity = polynomial_ring(
    ((3 * si**2 + a2_infinity) * (a2_infinity - si**2) / 4) % t**7
)
a3 = polynomial_ring([a3_infinity[6 - index] for index in range(7)])
a4_high = {8 - index: coefficient_ring(a4_infinity[index]) for index in range(7)}
a4_coefficients = [coefficient_ring.zero()] * 9
for degree, value in a4_high.items():
    a4_coefficients[degree] = value
a4_coefficients[1] = -sum(a4_coefficients[2:])
a4 = polynomial_ring(a4_coefficients)
if a4[0] or a4(1):
    raise ArithmeticError("failed to build the two marked-section node values into a4")

a2_lambda = a2(t + coefficient_ring(lambda_value))
a3_lambda = a3(t + coefficient_ring(lambda_value))
a4_lambda = a4(t + coefficient_ring(lambda_value))
lambda_relations = (
    a3_lambda - sl * (a2_lambda - sl**2),
    4 * a4_lambda - (3 * sl**2 + a2_lambda) * (a2_lambda - sl**2),
)

equations = [coefficient_ring(a3[0]), coefficient_ring(a3(1))]
for relation in lambda_relations:
    equations += [coefficient_ring(relation[index]) for index in range(7)]

# Local discriminant expansion at the marked node.  If a3=s*p, a4=s*q,
# a2=r, then Delta/s^2 starts with -16*r^2*(r*p^2-q^2).  Vanishing of its
# first two coefficients gives the displayed I4 equations after using the
# first relation in the second.  This avoids expanding the global degree-24
# discriminant after the infinity-jet substitution.
a2_one = a2(t + 1)
a3_one = a3(t + 1)
a4_one = a4(t + 1)
r0, r1 = coefficient_ring(a2_one[0]), coefficient_ring(a2_one[1])
p0, p1 = coefficient_ring(a3_one[1]), coefficient_ring(a3_one[2])
q0, q1 = coefficient_ring(a4_one[1]), coefficient_ring(a4_one[2])
k1 = r1 * p0**2 + 2 * r0 * p0 * p1 - 2 * q0 * q1
equations += [q0**2 - r0 * p0**2, q0**3 - 2 * r0**2 * k1]

# At zero a2(0)=3 is already nonzero.  The exact-I2 open is therefore the
# nonvanishing of K0=a2(0)*(a3/t)^2-(a4/t)^2.
exact_i2_open = coefficient_ring(3 * a3[1] ** 2 - a4[1] ** 2)
equations += [g["u0"] * exact_i2_open - 1]
base_field_split_variables = [
    "a2_4",
    "a2_3",
    "a2_2",
    "a2_1",
    "si_0",
    "sl_0",
]
if args.split_a2_node_constants_over_base_field and args.fixed_a2_node_case:
    raise SystemExit("the base-field split and one fixed case are mutually exclusive")
fixed_case_values = None
if args.fixed_a2_node_case:
    fixed_case_tokens = args.fixed_a2_node_case.split(",")
    if len(fixed_case_tokens) != len(base_field_split_variables):
        raise SystemExit("--fixed-a2-node-case requires six comma-separated values")
    fixed_case_values = [field(int(token)) for token in fixed_case_tokens]
    equations += [
        g[name] - value
        for name, value in zip(base_field_split_variables, fixed_case_values)
    ]
if args.split_a2_node_constants_over_base_field:
    equations += [
        g[name] ** int(prime) - g[name] for name in base_field_split_variables
    ]
expected_equations = len(names) + (
    len(base_field_split_variables)
    if args.split_a2_node_constants_over_base_field or fixed_case_values is not None
    else 0
)
if len(equations) != expected_equations or any(not equation for equation in equations):
    raise ArithmeticError("unexpected reduced equation accounting")

output_dir = args.output_dir.resolve()
output_dir.mkdir(parents=True, exist_ok=True)
stem = f"p{prime}-lambda{int(lambda_value)}" + (
    "-baseFieldA2Node0"
    if args.split_a2_node_constants_over_base_field
    else ""
)
if fixed_case_values is not None:
    stem += "-case-" + "-".join(str(int(value)) for value in fixed_case_values)
if args.compact_factored_msolve:
    stem += "-compact"
msolve_path = output_dir / f"{stem}.ms"
metadata_path = output_dir / f"{stem}.json"
msolve_text = ",".join(names) + "\n" + str(prime) + "\n"
equation_texts = [str(equation).replace("**", "^") for equation in equations]
if args.compact_factored_msolve:
    equation_texts[16:19] = [
        f"({q0})^2-({r0})*({p0})^2",
        f"({q0})^3-2*({r0})^2*({k1})",
        f"u0*({exact_i2_open})-1",
    ]
    equation_texts = [text.replace("**", "^") for text in equation_texts]
msolve_text += ",\n".join(equation_texts) + "\n"

metadata = {
    "schema": "elkies-k3.lattice-foundry-ns0007-pole0-reduced-modp-system.v1",
    "status": "PASS_EXACT_GLOBALLY_REDUCED_SECTION_FIRST_SYSTEM",
    "input": {
        "source_artifact": relative(SOURCE),
        "source_artifact_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "source_id": "NS0007-S025",
        "source_gram_sha256": source["gram_sha256"],
        "section_pole_artifact": relative(POLES),
        "section_pole_artifact_sha256": hashlib.sha256(POLES.read_bytes()).hexdigest(),
    },
    "prime": int(prime),
    "lambda": int(lambda_value),
    "source_profile": {
        "root_type": "A1+A3+2A6",
        "semistable_fibres": "I2+I4+2I7+4I1",
        "mw_height": "11/4",
        "minimum_section_pole_order": 0,
        "component_corrections_at_0_1_lambda_infinity": ["1/2", "3/4", "0", "0"],
    },
    "ansatz": {
        "equation": "y^2+a3*y=x^3+a2*x^2+a4*x",
        "marked_section": "P=(0,0)",
        "a2": str(a2),
        "a3_reconstructed_from_infinity_node_jet": str(a3),
        "a4_reconstructed_from_infinity_node_jet": str(a4),
        "lambda_node_parameter": str(sl),
        "infinity_node_parameter": str(si),
        "local_I4_equations": [str(equations[-3]), str(equations[-2])],
        "exact_I2_open": str(exact_i2_open),
    },
    "system": {
        "variables": names,
        "variable_count": len(names),
        "equation_count": len(equations),
        "condition_blocks": {
            "a3_node_values_at_0_and_1": 2,
            "split_I7_at_lambda": 14,
            "additional_I4_at_1": 2,
            "exact_I2_rabinowitsch": 1,
            "base_field_splitting_polynomials": (
                len(base_field_split_variables)
                if args.split_a2_node_constants_over_base_field
                else 0
            ),
        },
        "base_field_split_variables": (
            base_field_split_variables
            if args.split_a2_node_constants_over_base_field
            else []
        ),
        "fixed_a2_node_case": (
            {
                name: int(value)
                for name, value in zip(
                    base_field_split_variables, fixed_case_values
                )
            }
            if fixed_case_values is not None
            else None
        ),
        "equation_total_degrees": [int(equation.degree()) for equation in equations],
        "equation_term_counts": [len(equation.monomials()) for equation in equations],
        "msolve_polynomial_encoding": (
            "FACTORED_FINAL_I4_AND_EXACT_I2"
            if args.compact_factored_msolve
            else "FULLY_EXPANDED"
        ),
        "msolve_input": relative(msolve_path),
        "msolve_sha256": hashlib.sha256(msolve_text.encode()).hexdigest(),
    },
    "proof_boundary": (
        "The infinity split-node relations reconstruct a3 and a4 exactly, and "
        "the displayed square system imposes the marked nodes, split lambda I7, "
        "additional I4 order, and exact I2 open condition. Solutions still need "
        "independent full Kodaira-order, residual-squarefreeness, NS-marking, and "
        "characteristic-zero lift checks."
    ),
}
metadata_text = json.dumps(metadata, indent=2, sort_keys=True) + "\n"
if args.check:
    if msolve_path.read_text() != msolve_text:
        raise SystemExit("reduced NS0007 msolve system is stale")
    if metadata_path.read_text() != metadata_text:
        raise SystemExit("reduced NS0007 metadata is stale")
else:
    msolve_path.write_text(msolve_text)
    metadata_path.write_text(metadata_text)

print(
    "FOUNDRYNS0007POLE0REDUCED|"
    f"p={prime}|lambda={int(lambda_value)}|vars={len(names)}|eqs={len(equations)}|"
    f"max_degree={max(equation.degree() for equation in equations)}|status=PASS",
    flush=True,
)
print(f"MSOLVE_INPUT|{msolve_path}", flush=True)
print(f"OUTPUT|{metadata_path}", flush=True)
