#!/usr/bin/env sage-python
"""Build an A-eliminated Hermite system for the NS0034 fibre stratum.

For fixed ``p``, ``lambda``, ``A8``, and a chosen infinity square-root
constant, the square-root jets at zero and infinity determine all seven free
coefficients of A.  Only the jets at one and lambda, the remaining endpoint
jets, and two nonzero-root inverses remain.  This gives a 20-variable square
system instead of the unreduced 28-variable system.
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
DEFAULT_OUTPUT_DIR = ROOT / "artifacts/local/elkies-k3/ns0034-fibre-hermite-reduced-modp"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=7)
parser.add_argument("--lambda-value", type=int, default=2)
parser.add_argument("--A8", type=int, default=1)
parser.add_argument("--hi0", type=int, default=3)
parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

prime = ZZ(args.prime)
if not is_prime(prime) or prime in (2, 3):
    raise SystemExit("--prime must be a prime other than 2 or 3")
field = GF(prime)
lambda_value = field(args.lambda_value)
a8_value = field(args.A8)
hi0_value = field(args.hi0)
if lambda_value in (0, 1) or not a8_value:
    raise SystemExit("require lambda outside {0,1} and nonzero A8")
if hi0_value**2 + a8_value / 3:
    raise SystemExit("--hi0 must satisfy hi0^2=-A8/3")

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

names = ["h0_3", "h0_2", "h0_1"]
names += ["hi_4", "hi_3", "hi_2", "hi_1"]
names += [f"h1_{index}" for index in range(7, -1, -1)]
names += [f"hl_{index}" for index in range(2, -1, -1)]
names += ["v1", "vl"]
coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
g = coefficient_ring.gens_dict()
outer = PolynomialRing(coefficient_ring, "s")
s = outer.gen()

h0 = outer([1, g["h0_1"], g["h0_2"], g["h0_3"]])
hi = outer(
    [hi0_value, g["hi_1"], g["hi_2"], g["hi_3"], g["hi_4"]]
)
h0_square = h0**2
hi_square = hi**2
A_coefficients = [coefficient_ring(-3)]
A_coefficients += [-3 * coefficient_ring(h0_square[index]) for index in range(1, 4)]
A_coefficients += [
    -3 * coefficient_ring(hi_square[index]) for index in range(4, 0, -1)
]
A_coefficients += [coefficient_ring(a8_value)]
if len(A_coefficients) != 9:
    raise ArithmeticError("unexpected A reconstruction")


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


h1 = outer([g[f"h1_{index}"] for index in range(8)])
hl = outer([g[f"hl_{index}"] for index in range(3)])
A1_series = shifted_coefficients(field.one(), 8)
Al_series = shifted_coefficients(lambda_value, 3)
root_equations = [
    coefficient_ring((h1**2 + outer(A1_series) / 3)[index])
    for index in range(8)
]
root_equations += [
    coefficient_ring((hl**2 + outer(Al_series) / 3)[index])
    for index in range(3)
]

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
for h, precision in ((h0, 4), (h1, 8), (hl, 3), (hi, 5)):
    cube = 2 * h**3
    target += [coefficient_ring(cube[index]) for index in range(precision)]
compatibility_equations = list(
    left_kernel.change_ring(coefficient_ring)
    * vector(coefficient_ring, target)
)

equations = root_equations + compatibility_equations
equations += [g["v1"] * g["h1_0"] - 1, g["vl"] * g["hl_0"] - 1]
if len(equations) != len(names) or any(not equation for equation in equations):
    raise ArithmeticError("unexpected reduced square-system accounting")

output_dir = args.output_dir.resolve()
output_dir.mkdir(parents=True, exist_ok=True)
stem = (
    f"p{prime}-lambda{int(lambda_value)}-A8-{int(a8_value)}-"
    f"hi0-{int(hi0_value)}"
)
msolve_path = output_dir / f"{stem}.ms"
metadata_path = output_dir / f"{stem}.json"
msolve_text = ",".join(names) + "\n" + str(prime) + "\n"
msolve_text += ",\n".join(str(equation).replace("**", "^") for equation in equations) + "\n"

metadata = {
    "schema": "elkies-k3.lattice-foundry-ns0034-fibre-hermite-reduced-modp-system.v1",
    "status": "PASS_EXACT_A_ELIMINATED_SQUARE_HERMITE_FIBRE_SYSTEM",
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
    "hi0": int(hi0_value),
    "normalization": "A(0)=-3; h0(0)=1; enumerated split hi(0) branch",
    "fibre_profile": "I4+I8+I3+I5+4I1",
    "reconstructed_A_coefficients": [str(item) for item in A_coefficients],
    "system": {
        "variables": names,
        "variable_count": len(names),
        "equation_count": len(equations),
        "root_jet_equations_at_one_and_lambda": len(root_equations),
        "Hermite_compatibility_equations": len(compatibility_equations),
        "nonzero_root_inverse_equations": 2,
        "equation_total_degrees": [int(equation.degree()) for equation in equations],
        "equation_term_counts": [len(equation.monomials()) for equation in equations],
        "msolve_input": relative(msolve_path),
        "msolve_sha256": hashlib.sha256(msolve_text.encode()).hexdigest(),
    },
    "proof_boundary": (
        "The square system exactly parameterizes the fixed split nodal branch "
        "chart after eliminating all A coefficients. Solver points still require "
        "exact fibre-order, residual-squarefree, pole-zero section, NS-marking, "
        "and characteristic-zero lift verification."
    ),
}
metadata_text = json.dumps(metadata, indent=2, sort_keys=True) + "\n"
if args.check:
    if msolve_path.read_text() != msolve_text or metadata_path.read_text() != metadata_text:
        raise SystemExit("NS0034 reduced Hermite fibre system is stale")
else:
    msolve_path.write_text(msolve_text)
    metadata_path.write_text(metadata_text)

print(
    "FOUNDRYNS0034HERMITEREDUCED|"
    f"p={prime}|lambda={int(lambda_value)}|A8={int(a8_value)}|"
    f"hi0={int(hi0_value)}|vars={len(names)}|eqs={len(equations)}|"
    f"max_degree={max(equation.degree() for equation in equations)}|status=PASS",
    flush=True,
)
print(f"MSOLVE_INPUT|{msolve_path}", flush=True)
print(f"OUTPUT|{metadata_path}", flush=True)
