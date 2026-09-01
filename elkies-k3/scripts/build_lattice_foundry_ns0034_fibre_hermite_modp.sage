#!/usr/bin/env sage-python
"""Build a square Hermite system for the NS0034 semistable fibre stratum.

Fix ``lambda`` and the leading coefficient of A.  Formal nodal square-root
jets h at the I4, I8, I3, and I5 supports satisfy ``h^2=-A/3`` and determine
the B jets by ``B=2h^3``.  The seven left-kernel equations of the 20-by-13
Hermite matrix then give a sparse square system.  Nonzero-root inverses remove
the additive degeneracies at one and lambda.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ, binomial, is_prime, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-a-v1.json"
POLES = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-rank1-section-poles-v1.json"
DEFAULT_OUTPUT_DIR = ROOT / "artifacts/local/elkies-k3/ns0034-fibre-hermite-modp"
ORDERS = (4, 8, 3, 5)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=7)
parser.add_argument("--lambda-value", type=int, default=2)
parser.add_argument("--A8", type=int, default=1)
parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

prime = ZZ(args.prime)
if not is_prime(prime) or prime in (2, 3):
    raise SystemExit("--prime must be a prime other than 2 or 3")
field = GF(prime)
lambda_value = field(args.lambda_value)
a8_value = field(args.A8)
if lambda_value in (0, 1) or not a8_value:
    raise SystemExit("require lambda outside {0,1} and nonzero A8")

source_payload = json.loads(SOURCE.read_text())
source = next(
    row["source"]
    for row in source_payload["sources"]
    if row["ns_id"] == "NS0034" and row["source_id"] == "NS0034-S008"
)
assert source["root_type"] == "A2+A3+A4+A7"
assert source["mw_height_gram"] == [["19/8"]]
pole_payload = json.loads(POLES.read_text())
pole_row = next(
    row
    for row in pole_payload["sources"]
    if row["source_artifact"] == relative(SOURCE)
    and row["source_id"] == "NS0034-S008"
)
assert pole_row["minimum_section_pole_order"] == 0

jet_specs = (("h0", 4, True), ("h1", 8, False), ("hl", 3, False), ("hi", 5, False))
names = []
for prefix, precision, fixed_constant in jet_specs:
    top = precision - 1
    bottom = 1 if fixed_constant else 0
    names += [f"{prefix}_{index}" for index in range(top, bottom - 1, -1)]
names += [f"A{index}" for index in range(7, 0, -1)]
names += ["v1", "vl"]
coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
g = coefficient_ring.gens_dict()
outer = PolynomialRing(coefficient_ring, "s")
s = outer.gen()

A_coefficients = [coefficient_ring(-3)]
A_coefficients += [g[f"A{index}"] for index in range(1, 8)]
A_coefficients += [coefficient_ring(a8_value)]
A = outer(A_coefficients)

def shifted_coefficients(point, precision):
    return [
        sum(
            A_coefficients[index]
            * coefficient_ring(binomial(index, jet))
            * coefficient_ring(point) ** (index - jet)
            for index in range(jet, 9)
        )
        for jet in range(precision)
    ]


A_series = (
    A_coefficients[:4],
    shifted_coefficients(field.one(), 8),
    shifted_coefficients(lambda_value, 3),
    [A_coefficients[8 - index] for index in range(5)],
)
h_series = []
root_equations = []
for (prefix, precision, fixed_constant), local_A in zip(jet_specs, A_series):
    coefficients = [coefficient_ring.one()] if fixed_constant else [g[f"{prefix}_0"]]
    coefficients += [g[f"{prefix}_{index}"] for index in range(1, precision)]
    h = outer(coefficients)
    h_series.append(h)
    square_relation = h**2 + outer(local_A) / 3
    start = 1 if fixed_constant else 0
    root_equations += [coefficient_ring(square_relation[index]) for index in range(start, precision)]

# The B interpolation conditions use the same support order as the fibre
# profile.  No B variables are needed: the 20 target jets are 2*h^3.
rows = []
for point, precision in (
    (field.zero(), 4),
    (field.one(), 8),
    (lambda_value, 3),
):
    for jet in range(precision):
        rows.append(
            [
                field(binomial(index, jet)) * point ** (index - jet)
                if index >= jet
                else field.zero()
                for index in range(13)
            ]
        )
for jet in range(5):
    rows.append([field(index == 12 - jet) for index in range(13)])
hermite = matrix(field, rows)
if hermite.nrows() != 20 or hermite.ncols() != hermite.rank() or hermite.rank() != 13:
    raise ArithmeticError("unexpected Hermite rank")
left_kernel = hermite.left_kernel().basis_matrix()
target = []
for h, unused_precision in zip(h_series, ORDERS):
    cube = 2 * h**3
    target += [coefficient_ring(cube[index]) for index in range(unused_precision)]
compatibility_vector = left_kernel.change_ring(coefficient_ring) * vector(coefficient_ring, target)
compatibility_equations = list(compatibility_vector)

equations = root_equations + compatibility_equations
equations += [g["v1"] * g["h1_0"] - 1, g["vl"] * g["hl_0"] - 1]
if len(equations) != len(names) or any(not equation for equation in equations):
    raise ArithmeticError("unexpected square-system accounting")

output_dir = args.output_dir.resolve()
output_dir.mkdir(parents=True, exist_ok=True)
stem = f"p{prime}-lambda{int(lambda_value)}-A8-{int(a8_value)}"
msolve_path = output_dir / f"{stem}.ms"
metadata_path = output_dir / f"{stem}.json"
msolve_text = ",".join(names) + "\n" + str(prime) + "\n"
msolve_text += ",\n".join(str(equation).replace("**", "^") for equation in equations) + "\n"

metadata = {
    "schema": "elkies-k3.lattice-foundry-ns0034-fibre-hermite-modp-system.v1",
    "status": "PASS_EXACT_SQUARE_HERMITE_FIBRE_SYSTEM",
    "input": {
        "source_artifact": relative(SOURCE),
        "source_artifact_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "source_id": "NS0034-S008",
        "source_gram_sha256": source["gram_sha256"],
        "section_pole_artifact": relative(POLES),
        "section_pole_artifact_sha256": hashlib.sha256(POLES.read_bytes()).hexdigest(),
    },
    "prime": int(prime),
    "lambda": int(lambda_value),
    "A8": int(a8_value),
    "normalization": "A(0)=-3; positive h0(0)=1 branch",
    "fibre_profile": "I4+I8+I3+I5+4I1",
    "system": {
        "variables": names,
        "variable_count": len(names),
        "equation_count": len(equations),
        "root_jet_equations": len(root_equations),
        "Hermite_compatibility_equations": len(compatibility_equations),
        "nonzero_root_inverse_equations": 2,
        "equation_total_degrees": [int(equation.degree()) for equation in equations],
        "equation_term_counts": [len(equation.monomials()) for equation in equations],
        "msolve_input": relative(msolve_path),
        "msolve_sha256": hashlib.sha256(msolve_text.encode()).hexdigest(),
    },
    "proof_boundary": (
        "The square system exactly parameterizes split nodal branch jets and B "
        "Hermite compatibility for the fixed (lambda,A8) chart. Solver points "
        "still require exact fibre-order, residual-squarefree, pole-zero section, "
        "NS-marking, and characteristic-zero lift verification."
    ),
}
metadata_text = json.dumps(metadata, indent=2, sort_keys=True) + "\n"
if args.check:
    if msolve_path.read_text() != msolve_text or metadata_path.read_text() != metadata_text:
        raise SystemExit("NS0034 Hermite fibre system is stale")
else:
    msolve_path.write_text(msolve_text)
    metadata_path.write_text(metadata_text)

print(
    "FOUNDRYNS0034HERMITE|"
    f"p={prime}|lambda={int(lambda_value)}|A8={int(a8_value)}|"
    f"vars={len(names)}|eqs={len(equations)}|"
    f"max_degree={max(equation.degree() for equation in equations)}|status=PASS",
    flush=True,
)
print(f"MSOLVE_INPUT|{msolve_path}", flush=True)
print(f"OUTPUT|{metadata_path}", flush=True)
