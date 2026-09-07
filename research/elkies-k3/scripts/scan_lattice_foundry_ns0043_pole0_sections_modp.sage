#!/usr/bin/env sage-python
"""Scan the NS0043 height-four pole-zero section chart on stored fibres.

The selected primitive source has root type ``A2+A6+A8``, MW height four,
and minimum section pole order zero.  Hence its Shioda correction sum is zero:
the marked section must specialize to the smooth identity component at the
I9, I7, and I3 fibres.  For each polynomial X of degree at most four, test
whether ``X^3+A*X+B`` is a polynomial square of degree at most twelve.  This
exhausts the full pole-zero section chart using only ``p^5`` X-polynomials.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import GF, PolynomialRing


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0011-source-ansatz-mod5.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0043-pole0-sections-mod5.json"
SOURCE = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-a-v1.json"
POLES = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-rank1-section-poles-v1.json"
SOURCE_ID = "NS0043-S005"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--examples", type=int, default=100)
parser.add_argument(
    "--quadratic-twist",
    type=int,
    default=1,
    help="replace (A,B) by (d^2 A,d^3 B) before scanning sections",
)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

input_path = args.input.resolve()
payload = json.loads(input_path.read_text())
if payload["ansatz"]["normalized_reducible_supports"] != [
    "0:I9",
    "1:I7",
    "infinity:I3",
]:
    raise ValueError("input does not have the NS0043 semistable fibre profile")
prime = int(payload["prime"])
if prime in (2, 3):
    raise ValueError("section scan requires odd characteristic other than three")

source_payload = json.loads(SOURCE.read_text())
source_row = next(
    row
    for row in source_payload["sources"]
    if row["ns_id"] == "NS0043" and row["source_id"] == SOURCE_ID
)
source = source_row["source"]
assert source["root_type"] == "A2+A6+A8"
assert source["mw_height_gram"] == [["4"]]
assert source["root_lattice_primitive"] and source["torsion"] == 1
assert source_row["determinant"] == 756
assert all(
    target["mw_rank_for_rho_19"] == 15
    for target in source_row["same_ns_high_rank_targets"]
)
pole_payload = json.loads(POLES.read_text())
pole = next(
    row
    for row in pole_payload["sources"]
    if row["source_artifact"] == relative(SOURCE)
    and row["source_id"] == SOURCE_ID
)
assert pole["source_gram_sha256"] == source["gram_sha256"]
assert pole["minimum_section_pole_order"] == 0
assert pole["mw_height"] == "4"

field = GF(prime)
twist = field(args.quadratic_twist)
if not twist:
    raise ValueError("--quadratic-twist must be nonzero")
ring = PolynomialRing(field, "t")
t = ring.gen()
records = []
stored_solutions = []
for example_index, example in enumerate(payload["examples"]):
    A = twist**2 * ring(example["A_coefficients_low_to_high"])
    B = twist**3 * ring(example["B_coefficients_low_to_high"])
    nodes = {
        "I9_at_zero": -field(3) * B(0) / (field(2) * A(0)),
        "I7_at_one": -field(3) * B(1) / (field(2) * A(1)),
        "I3_at_infinity": -field(3) * B[12] / (field(2) * A[8]),
    }
    raw_count = 0
    marked_count = 0
    local_examples = []
    raw_examples = []
    component_pattern_counts = {}
    for x_values in itertools.product(field, repeat=5):
        X = ring(x_values)
        right = X**3 + A * X + B
        if not right.is_square():
            continue
        positive_Y = right.sqrt()
        y_values = [positive_Y] if not positive_Y else [positive_Y, -positive_Y]
        for Y in y_values:
            raw_count += 1
            if Y**2 != right:
                raise ArithmeticError("polynomial square-root audit failed")
            smooth = {
                "I9_at_zero": not (X(0) == nodes["I9_at_zero"] and Y(0) == 0),
                "I7_at_one": not (X(1) == nodes["I7_at_one"] and Y(1) == 0),
                "I3_at_infinity": not (
                    X[4] == nodes["I3_at_infinity"] and Y[6] == 0
                ),
            }
            pattern = ",".join(
                f"{label}={int(value)}" for label, value in smooth.items()
            )
            component_pattern_counts[pattern] = component_pattern_counts.get(pattern, 0) + 1
            raw_diagnostic = {
                "X_coefficients_low_to_high": [int(value) for value in X],
                "Y_coefficients_low_to_high": [int(value) for value in Y],
                "smooth_identity_components": smooth,
            }
            if len(raw_examples) < args.examples:
                raw_examples.append(raw_diagnostic)
            if not all(smooth.values()):
                continue
            marked_count += 1
            diagnostic = raw_diagnostic | {
                "shioda_height": "4",
                "implied_NS_determinant": 756,
            }
            if len(local_examples) < args.examples:
                local_examples.append(diagnostic)
            if len(stored_solutions) < args.examples:
                stored_solutions.append({"fibre_example_index": example_index} | diagnostic)
    records.append(
        {
            "fibre_example_index": example_index,
            "X_polynomials_scanned": prime**5,
            "raw_polynomial_section_count": raw_count,
            "component_pattern_counts": component_pattern_counts,
            "stored_raw_examples": raw_examples,
            "marked_section_count": marked_count,
            "stored_examples": local_examples,
        }
    )

total_raw = sum(row["raw_polynomial_section_count"] for row in records)
total_marked = sum(row["marked_section_count"] for row in records)
output = {
    "schema": "elkies-k3.lattice-foundry-ns0043-pole0-sections-modp-scan.v1",
    "status": (
        "PASS_EXACT_EXHAUSTIVE_STORED_MODELS_WITH_MARKED_SECTIONS"
        if total_marked
        else "PASS_EXACT_EXHAUSTIVE_STORED_MODELS_EMPTY_SECTION_CHART"
    ),
    "prime": prime,
    "quadratic_twist": int(twist),
    "quadratic_twist_square_class": "square" if twist.is_square() else "nonsquare",
    "input": {
        "fibre_artifact": relative(input_path),
        "fibre_artifact_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "source_artifact": relative(SOURCE),
        "source_artifact_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "source_id": SOURCE_ID,
        "source_gram_sha256": source["gram_sha256"],
        "section_pole_artifact": relative(POLES),
        "section_pole_artifact_sha256": hashlib.sha256(POLES.read_bytes()).hexdigest(),
    },
    "source": {
        "root_type": "A2+A6+A8",
        "mw_rank": 1,
        "mw_height": "4",
        "minimum_section_pole_order": 0,
        "component_corrections": ["0", "0", "0"],
        "determinant": 756,
        "same_ns_high_rank_targets": source_row["same_ns_high_rank_targets"],
    },
    "scope": {
        "stored_fibre_models": len(records),
        "X_polynomials_per_model": prime**5,
        "all_Y_square_roots_retained": True,
        "fibre_ansatz_scan_exhausted": bool(payload["scan"]["exhausted"]),
        "twist_is_bijective_on_stored_fibre_models": True,
    },
    "accounting": {
        "total_X_polynomials_scanned": len(records) * prime**5,
        "raw_polynomial_sections": total_raw,
        "marked_polynomial_sections": total_marked,
        "models_with_marked_sections": sum(bool(row["marked_section_count"]) for row in records),
        "stored_marked_sections": len(stored_solutions),
    },
    "models": records,
    "stored_marked_sections": stored_solutions,
    "proof_boundary": {
        "proved": (
            "Every pole-zero polynomial section is exhausted on every stored "
            "finite-field fibre model. Retained sections satisfy the exact "
            "equation, meet the smooth identity components at I9, I7, and I3, "
            "and therefore have Shioda height four and determinant 756."
        ),
        "not_proved": (
            "A finite-field marked model is not a rational source family, a "
            "characteristic-zero lift, or a physical neighbour corridor."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/scan_lattice_foundry_ns0043_pole0_sections_modp.sage"
        + (f" --quadratic-twist {int(twist)}" if twist != 1 else "")
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
output_path = args.output.resolve()
if args.check:
    if output_path.read_text() != serialized:
        raise SystemExit("NS0043 pole-zero section scan is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "FOUNDRYNS0043POLE0SCAN|"
    f"models={len(records)}|X={len(records) * prime**5}|"
    f"raw={total_raw}|marked={total_marked}|status=PASS",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
