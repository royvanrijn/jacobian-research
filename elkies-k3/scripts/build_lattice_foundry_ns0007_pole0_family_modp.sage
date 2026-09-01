#!/usr/bin/env sage-python
"""Build the section-first NS0007 pole-zero source system modulo p.

Translate the exact pole-zero MW generator to ``P=(0,0)`` and use

    y^2 + a3(t)y = x^3 + a2(t)x^2 + a4(t)x.

The component corrections ``1/2,3/4,0,0`` place P at the nodes of the I2
and I4 fibres and on the identity components of both I7 fibres.  Normalize
the supports to ``0,1,lambda,infinity`` and set

    a3=t(t-1)c3,  a4=t(t-1)c4,  a2(0)=3.

After dividing the discriminant by ``t^2(t-1)^2``, two conditions at one,
seven at lambda, and seven at infinity give a square system of 16 equations
in 16 variables for fixed lambda.  This is an exact finite-field discovery
chart, not a characteristic-zero family or a rational-lift certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ, is_prime


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-c-v1.json"
)
DEFAULT_POLES = (
    ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-rank1-section-poles-v1.json"
)
DEFAULT_OUTPUT_DIR = ROOT / "artifacts/local/elkies-k3/ns0007-pole0-family-modp"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
parser.add_argument("--source-id", default="NS0007-S025")
parser.add_argument("--section-poles", type=Path, default=DEFAULT_POLES)
parser.add_argument("--prime", type=int, default=7)
parser.add_argument("--lambda-value", type=int, default=2)
parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
parser.add_argument(
    "--saturate-exact-I2",
    action="store_true",
    help="add u0*H(0)-1 to remove higher-order/globally singular components",
)
parser.add_argument(
    "--parametrize-I7-nodes",
    action="store_true",
    help="replace dense discriminant I7 conditions by sparse split-node jets",
)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

prime = ZZ(args.prime)
if not is_prime(prime) or prime in (2, 3):
    raise SystemExit("--prime must be a prime other than 2 or 3")
field = GF(prime)
lambda_value = field(args.lambda_value)
if lambda_value in (field.zero(), field.one()):
    raise SystemExit("--lambda-value must differ from 0 and 1 modulo p")

source_path = args.source.resolve()
pole_path = args.section_poles.resolve()
source_payload = json.loads(source_path.read_text())
source_entry = next(
    row
    for row in source_payload["sources"]
    if row["ns_id"] == "NS0007" and row["source_id"] == args.source_id
)
source = source_entry["source"]
assert source["root_type"] == "A1+A3+2A6"
assert source["root_lattice_primitive"] and source["torsion"] == 1
assert source["mw_height_gram"] == [["11/4"]]
pole_payload = json.loads(pole_path.read_text())
pole_row = next(
    row
    for row in pole_payload["sources"]
    if row["source_artifact"] == relative(source_path)
    and row["source_id"] == args.source_id
)
assert pole_row["minimum_section_pole_order"] == 0

# High coefficients first is generally favorable for the seven triangular
# infinity equations in degree-reverse-lexicographic order.
base_names = [f"a2_{index}" for index in range(4, 0, -1)]
base_names += [f"c3_{index}" for index in range(4, -1, -1)]
base_names += [f"c4_{index}" for index in range(6, -1, -1)]
if args.parametrize_I7_nodes:
    jet_names = [f"sl_{index}" for index in range(6, -1, -1)]
    jet_names += [f"si_{index}" for index in range(6, -1, -1)]
    names = jet_names + base_names
else:
    names = base_names
if args.saturate_exact_I2:
    names += ["u0"]
coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
generators = coefficient_ring.gens_dict()
polynomial_ring = PolynomialRing(coefficient_ring, "t")
t = polynomial_ring.gen()

a2 = coefficient_ring(3) + sum(
    generators[f"a2_{index}"] * t**index for index in range(1, 5)
)
c3 = sum(generators[f"c3_{index}"] * t**index for index in range(5))
c4 = sum(generators[f"c4_{index}"] * t**index for index in range(7))
D = t * (t - 1)
a3 = D * c3
a4 = D * c4

b2 = 4 * a2
b4 = 2 * a4
b6 = a3**2
b8 = a2 * a3**2 - a4**2
discriminant = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
H, remainder = discriminant.quo_rem(D**2)
if remainder or H.degree() > 20:
    raise ArithmeticError("failed to remove the automatic I2+I2 discriminant factor")

at_one = H(t + 1)
equations = [coefficient_ring(at_one[index]) for index in range(2)]
if args.parametrize_I7_nodes:
    sl = sum(generators[f"sl_{index}"] * t**index for index in range(7))
    si = sum(generators[f"si_{index}"] * t**index for index in range(7))
    a2_lambda = a2(t + coefficient_ring(lambda_value))
    a3_lambda = a3(t + coefficient_ring(lambda_value))
    a4_lambda = a4(t + coefficient_ring(lambda_value))
    lambda_relations = (
        a3_lambda - sl * (a2_lambda - sl**2),
        4 * a4_lambda - (3 * sl**2 + a2_lambda) * (a2_lambda - sl**2),
    )
    a2_infinity = polynomial_ring([a2[4 - index] for index in range(5)])
    a3_infinity = polynomial_ring([a3[6 - index] for index in range(7)])
    a4_infinity = polynomial_ring([a4[8 - index] for index in range(9)])
    infinity_relations = (
        a3_infinity - si * (a2_infinity - si**2),
        4 * a4_infinity - (3 * si**2 + a2_infinity) * (a2_infinity - si**2),
    )
    for relation in lambda_relations + infinity_relations:
        equations += [coefficient_ring(relation[index]) for index in range(7)]
else:
    at_lambda = H(t + coefficient_ring(lambda_value))
    equations += [coefficient_ring(at_lambda[index]) for index in range(7)]
    equations += [coefficient_ring(H[index]) for index in range(14, 21)]
if args.saturate_exact_I2:
    equations.append(generators["u0"] * coefficient_ring(H[0]) - 1)
expected_equations = 30 if args.parametrize_I7_nodes else 16
expected_equations += int(args.saturate_exact_I2)
if len(equations) != expected_equations or any(not equation for equation in equations):
    raise ArithmeticError("unexpected section-first equation accounting")

output_dir = args.output_dir.resolve()
output_dir.mkdir(parents=True, exist_ok=True)
stem = f"p{prime}-lambda{int(lambda_value)}" + (
    "-nodeJets" if args.parametrize_I7_nodes else ""
) + (
    "-satI2" if args.saturate_exact_I2 else ""
)
msolve_path = output_dir / f"{stem}.ms"
metadata_path = output_dir / f"{stem}.json"
msolve_text = ",".join(names) + "\n" + str(prime) + "\n"
msolve_text += ",\n".join(str(equation).replace("**", "^") for equation in equations)
msolve_text += "\n"

metadata = {
    "schema": "elkies-k3.lattice-foundry-ns0007-pole0-family-modp-system.v1",
    "status": "PASS_EXACT_SECTION_FIRST_MODULAR_SYSTEM",
    "input": {
        "source_artifact": relative(source_path),
        "source_artifact_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "source_id": args.source_id,
        "source_gram_sha256": source["gram_sha256"],
        "section_pole_artifact": relative(pole_path),
        "section_pole_artifact_sha256": hashlib.sha256(pole_path.read_bytes()).hexdigest(),
    },
    "prime": int(prime),
    "lambda": int(lambda_value),
    "source_profile": {
        "root_type": "A1+A3+2A6",
        "semistable_fibres": "I2+I4+2I7+4I1",
        "mw_rank": 1,
        "mw_height": "11/4",
        "minimum_section_pole_order": 0,
        "component_corrections_at_0_1_lambda_infinity": ["1/2", "3/4", "0", "0"],
    },
    "ansatz": {
        "equation": "y^2+a3*y=x^3+a2*x^2+a4*x",
        "marked_section": "P=(0,0)",
        "a2": str(a2),
        "a3": str(a3),
        "a4": str(a4),
        "discriminant": str(discriminant),
        "discriminant_after_t2_tminus1_2": str(H),
    },
    "system": {
        "variables": names,
        "variable_count": len(names),
        "equation_count": len(equations),
        "condition_blocks": {
            "additional_I4_at_1": 2,
            "I7_at_lambda": 14 if args.parametrize_I7_nodes else 7,
            "I7_at_infinity": 14 if args.parametrize_I7_nodes else 7,
            "exact_I2_rabinowitsch": int(args.saturate_exact_I2),
        },
        "split_I7_node_parametrization": bool(args.parametrize_I7_nodes),
        "equation_total_degrees": [int(equation.degree()) for equation in equations],
        "equation_term_counts": [len(equation.monomials()) for equation in equations],
        "msolve_input": relative(msolve_path),
        "msolve_sha256": hashlib.sha256(msolve_text.encode()).hexdigest(),
    },
    "proof_boundary": (
        "Exact symbolic algebra builds the fixed-cross-ratio finite-field system "
        "with the pole-zero section and required discriminant divisibilities. "
        "A solution must still pass exact Kodaira-order, squarefree-residual, "
        "identity-component, NS determinant, and characteristic-zero lift gates."
    ),
}
metadata_text = json.dumps(metadata, indent=2, sort_keys=True) + "\n"
if args.check:
    if msolve_path.read_text() != msolve_text:
        raise SystemExit("NS0007 section-first msolve system is stale")
    if metadata_path.read_text() != metadata_text:
        raise SystemExit("NS0007 section-first metadata is stale")
else:
    msolve_path.write_text(msolve_text)
    metadata_path.write_text(metadata_text)

print(
    "FOUNDRYNS0007POLE0BUILD|"
    f"p={prime}|lambda={int(lambda_value)}|vars={len(names)}|eqs={len(equations)}|"
    f"max_degree={max(equation.degree() for equation in equations)}|status=PASS",
    flush=True,
)
print(f"MSOLVE_INPUT|{msolve_path}", flush=True)
print(f"OUTPUT|{metadata_path}", flush=True)
