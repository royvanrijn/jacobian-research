#!/usr/bin/env sage -python
"""Evaluate the degree-44 third-divisor generic ambient in the actual E7 quotient.

The 43 ordinary Weierstrass monomials and the already evaluated marked chord
are mapped to the length-363 actual H92 E7 quotient.  This is the reusable
local evaluator for the high-degree Riemann--Roch ambient; it intentionally
does not claim that this single E7 block is the full global condition matrix.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, QQ


ROOT = Path(__file__).resolve().parents[2]
AMBIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-generic-rr-ambient.json"
QUOTIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-actual-e7-quotient-block.json"
CHORD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-marked-chord-actual-e7-quotient.json"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-generic-ambient-actual-e7-quotient.json"
AMBIENT_SHA256 = "598480530e0c845fcecb1bdcdef4756dec61f668f0aed07c3c71ccc1049b6878"
QUOTIENT_SHA256 = "7848c3a506b2255fc1e42cab9ced9b72f8852e26aa703f9b23e2b2417474d2ed"
CHORD_SHA256 = "bbfdd659a38dbdbe6f0fd272970df251049ceb4f9344068de7f2f4abe9a6ed9f"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--ambient", type=Path, default=AMBIENT)
parser.add_argument("--quotient", type=Path, default=QUOTIENT)
parser.add_argument("--chord", type=Path, default=CHORD)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
for selected, default, expected in ((args.ambient, AMBIENT, AMBIENT_SHA256),
                                   (args.quotient, QUOTIENT, QUOTIENT_SHA256),
                                   (args.chord, CHORD, CHORD_SHA256)):
    if selected == default:
        assert digest(selected) == expected
ambient = json.loads(args.ambient.read_text())
quotient = json.loads(args.quotient.read_text())
chord_data = json.loads(args.chord.read_text())
assert ambient["status"] == "PASS_EXACT_Q6_THIRD_SYMBOLIC_GENERIC_RR_AMBIENT"
assert quotient["status"] == "PASS_EXACT_Q6_THIRD_ACTUAL_E7_QUOTIENT_BLOCK"
assert chord_data["status"] == "PASS_EXACT_Q6_THIRD_MARKED_CHORD_ACTUAL_E7_QUOTIENT"
assert ambient["dimension"] == 44 and quotient["quotient_dimension"] == 363

basis = [tuple(entry) for entry in quotient["quotient_basis_exponents"]]
assert [list(entry) for entry in basis] == chord_data["basis_exponents"]
chord = {entry: QQ(value) for entry, value in zip(basis, chord_data["coordinates_in_basis_order"]) if value != "0"}
assert len(chord) == chord_data["chord_support"] == 363
generators = {tuple(entry) for entry in quotient["complete_ideal_generators"]}
c2 = QQ(quotient["translated_coordinates"]["c2"])
c3 = QQ(quotient["translated_coordinates"]["c3"])


def add(left, right):
    answer = dict(left)
    for key, value in right.items():
        answer[key] = answer.get(key, QQ(0))+value
        if not answer[key]:
            del answer[key]
    return answer


def in_ideal(exponent):
    return any(all(left <= right for left, right in zip(generator, exponent)) for generator in generators)


def multiply_raw(left, right):
    answer = {}
    for (i, a, b), coefficient in left.items():
        for (j, c, d), other in right.items():
            key = (i+j, a+c, b+d)
            answer[key] = answer.get(key, QQ(0))+coefficient*other
    return {key: value for key, value in answer.items() if value}


x = {(0, 1, 0): QQ(1), (2, 0, 0): c2, (3, 0, 0): c3}
a = {(3, 0, 0): QQ(quotient["translated_coordinates"].get("A1", 0))}
# Recover the exact Weierstrass right hand side from the published relation
# coefficients implicit in the quotient's translated H92 block.  The needed
# terms are supplied directly by the source formula through a small anchor.
from importlib.machinery import SourceFileLoader
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
anchor = SourceFileLoader("h92_third_ambient_quotient_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)
rhs_y2 = multiply_raw(multiply_raw(x, x), x)
rhs_y2 = add(rhs_y2, multiply_raw({(3, 0, 0): A1, (4, 0, 0): A}, x))
rhs_y2 = add(rhs_y2, {(5, 0, 0): B1, (6, 0, 0): B, (7, 0, 0): B2})


def reduce(value):
    pending = list(value.items())
    answer = {}
    while pending:
        (i, a_power, b_power), coefficient = pending.pop()
        if not coefficient or in_ideal((i, a_power, b_power)):
            continue
        if b_power < 2:
            key = (i, a_power, b_power)
            answer[key] = answer.get(key, QQ(0))+coefficient
            continue
        for (j, c_power, d_power), other in rhs_y2.items():
            pending.append(((i+j, a_power+c_power, b_power-2+d_power), coefficient*other))
    return {key: value for key, value in answer.items() if value}


def multiply(left, right):
    return reduce(multiply_raw(left, right))


ordinary = []
ordinary_labels = []
for entry in ambient["basis"][:-1]:
    assert entry["kind"] == "monomial"
    value = {(0, 0, 0): QQ(1)}
    for unused in range(entry["x_power"]):
        value = multiply(value, x)
    for unused in range(entry["y_power"]):
        value = multiply(value, {(0, 0, 1): QQ(1)})
    ordinary.append(value)
    ordinary_labels.append({
        "kind": "monomial",
        "x_power": int(entry["x_power"]),
        "y_power": int(entry["y_power"]),
    })
assert len(ordinary) == 43
columns = ordinary+[chord]
ambient_labels = ordinary_labels+[{"kind": "marked_chord"}]

prime = GF(1009)
modular_columns = []
for value in columns:
    modular_columns.append([
        prime(coefficient.numerator())/prime(coefficient.denominator())
        for coefficient in (value.get(exponent, QQ(0)) for exponent in basis)
    ])
from sage.all import matrix
modular_matrix = matrix(prime, modular_columns).transpose()
modular_rank = modular_matrix.rank()
print("H92Q6THIRDAMBIENTE7|modular_rank={}".format(modular_rank), flush=True)
assert modular_rank >= 1
exact_matrix = matrix(QQ, [
    [value.get(exponent, QQ(0)) for value in columns]
    for exponent in basis
])
exact_rank = exact_matrix.rank()
assert exact_rank >= modular_rank

# This local kernel has a sparse exact description, so do not delegate its
# materialization to a dense rational nullspace routine.  The only zero
# columns are the high-pole ordinary monomials; rank-nullity proves that their
# coordinate vectors are the whole E7 kernel.  This is an exact local result,
# not a substitute for the still-required global E8 and smooth blocks.
zero_indices = tuple(index for index, value in enumerate(columns) if not value)
assert len(zero_indices) == 44-exact_rank == 11
exec(compile(CORE.read_text(), str(CORE), "exec"))
ambient_indices = tuple(range(len(columns)))
e7_block = quotient_condition(
    "actual_H92_E7_third_complete_ideal",
    ambient_indices,
    lambda column: tuple(columns[column].get(exponent, QQ(0)) for exponent in basis),
    tuple("T^{}*U^{}*Y^{}".format(*exponent) for exponent in basis),
    "length-363 actual H92 E7 complete quotient in translated (T,U,Y) coordinates",
)
compiler_replay = compile_resolved_conditions(
    ambient_indices, (e7_block,), complete=False
)
assert compiler_replay["condition_matrix"] == exact_matrix
assert (compiler_replay["rank"], compiler_replay["kernel_dimension"],
        compiler_replay["kernel_materialization"]) == (33, 11, "zero_columns")
kernel_basis = compiler_replay["kernel_basis"]
assert kernel_basis.rank() == len(zero_indices)
assert exact_matrix*kernel_basis.transpose() == matrix(QQ, exact_matrix.nrows(), len(zero_indices))
nonzero_indices = tuple(index for index in range(len(columns)) if index not in zero_indices)
assert matrix(QQ, [[exact_matrix[row, column] for column in nonzero_indices]
                   for row in range(exact_matrix.nrows())]).rank() == exact_rank
exact_kernel_labels = [ambient_labels[index] for index in zero_indices]

ordinary_records = []
for entry, value in zip(ambient["basis"][:-1], ordinary):
    coordinates = [str(value.get(exponent, QQ(0))) for exponent in basis]
    ordinary_records.append({
        "x_power": entry["x_power"], "y_power": entry["y_power"],
        "support": len(value),
        "coordinate_sha256": hashlib.sha256(json.dumps(coordinates, separators=(",", ":")).encode()).hexdigest(),
    })
payload = {
    "schema": "elkies-k3.h92-q6-third-generic-ambient-actual-e7-quotient.v1",
    "status": "PASS_EXACT_Q6_THIRD_GENERIC_AMBIENT_ACTUAL_E7_EVALUATION",
    "inputs": {
        "generic_ambient": {"path": str(args.ambient.relative_to(ROOT)), "sha256": digest(args.ambient)},
        "actual_e7_quotient": {"path": str(args.quotient.relative_to(ROOT)), "sha256": digest(args.quotient)},
        "marked_chord": {"path": str(args.chord.relative_to(ROOT)), "sha256": digest(args.chord)},
    },
    "quotient_dimension": 363,
    "ambient_dimension": 44,
    "ordinary_monomials": ordinary_records,
    "marked_chord_coordinate_sha256": chord_data["coordinate_sha256"],
    "rank_certificate": {"prime": 1009, "modular_rank": int(modular_rank), "exact_rank": int(exact_rank), "kernel_dimension": int(44-exact_rank)},
    "compiler_replay": {
        "core": {"path": str(CORE.relative_to(ROOT)), "sha256": digest(CORE)},
        "condition_block": e7_block["name"],
        "complete_resolved_chart_cover": compiler_replay["complete_resolved_chart_cover"],
        "kernel_materialization": compiler_replay["kernel_materialization"],
    },
    "exact_kernel": {
        "basis": exact_kernel_labels,
        "description": "The listed coordinate monomials vanish individually in the length-363 E7 quotient; exact rank-nullity proves that they span the full local kernel.",
    },
    "boundary": "This is the E7 evaluation map for the generic degree-44 ambient. A global high-degree coefficient ambient and its E8/smooth blocks still must be derived before solving for a transported section.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q6THIRDAMBIENTE7|ambient=44|quotient=363|rank={}|kernel={}|"
    "status=PASS_EXACT_Q6_THIRD_GENERIC_AMBIENT_ACTUAL_E7_EVALUATION".format(exact_rank, 44-exact_rank),
    flush=True,
)
