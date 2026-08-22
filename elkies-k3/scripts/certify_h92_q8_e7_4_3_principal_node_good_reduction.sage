#!/usr/bin/env sage -python
"""Lift a full-rank actual principal-node reduction to characteristic zero.

The modular local-normal-form evaluator works in the actual resolved local
ring and represents only the finite image of a fixed ambient in
``R/(t^17)``.  If that image has full ambient rank at a prime where every
clearing factor is a unit, then the corresponding characteristic-zero local
membership map is injective.  This is a bounded-ambient rejection result; it
does not claim a finite presentation of ``R/(t^17)`` or a full q8 cover.
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROBE = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-principal-node-local-normal-form-mod-43.json"
PROBE_SOURCE = ROOT / "elkies-k3/scripts/probe_h92_q8_e7_4_3_principal_node_local_normal_form_modp.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-principal-node-good-reduction.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--probe", type=Path, default=PROBE)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

probe = json.loads(args.probe.read_text())
assert probe["status"] == "EXPERIMENTAL_MODULAR_Q8_E7_4_3_LOCAL_NORMAL_FORM_BLOCK"
assert probe["inputs"]["checker_source"] == {
    "path": str(PROBE_SOURCE.relative_to(ROOT)), "sha256": digest(PROBE_SOURCE),
}
assert int(probe["prime"]) == 43
assert probe["local_ring"] == {
    "chart": "E7_4--E7_3",
    "order": "Singular ds local degree order at (Z,U,Y)",
    "surface_equation": probe["local_ring"]["surface_equation"],
    "principal_ideal": "(surface,(Z^3*U^2)^17)",
    "completed_interpretation": "t^17=Z^51*Y^68*unit",
}
assert probe["common_clearing"]["T"] == 17
assert probe["good_reduction"]["all_input_coefficient_denominators_nonzero"]
assert all(probe["good_reduction"]["node_unit_residues"].values())
image = probe["finite_ambient_image"]
assert image["ambient_dimension"] == 54
assert image["rank"] == image["ambient_dimension"]
assert image["kernel_dimension"] == 0

payload = {
    "schema": "elkies-k3.h92-q8-e7-4-3-principal-node-good-reduction.v1",
    "status": "PASS_EXACT_Q8_E7_4_3_PRINCIPAL_NODE_INJECTIVITY",
    "inputs": {
        "local_normal_form_probe": {
            "path": str(args.probe.relative_to(ROOT)), "sha256": digest(args.probe),
        },
        "probe_source": {
            "path": str(PROBE_SOURCE.relative_to(ROOT)), "sha256": digest(PROBE_SOURCE),
        },
    },
    "bounded_ambient": {
        "dimension": image["ambient_dimension"],
        "chart": "E7_4--E7_3",
        "principal_condition": "common-cleared numerator belongs to (t^17) in the actual resolved local ring",
    },
    "good_reduction": {
        "prime": probe["prime"],
        "local_image_rank": image["rank"],
        "local_image_kernel_dimension": image["kernel_dimension"],
        "clearing_unit_residues": probe["good_reduction"]["node_unit_residues"],
        "argument": (
            "A nonzero QQ kernel vector may be cleared to a primitive integral "
            "vector before reduction. The chart denominators and common-clearing "
            "factors are units modulo 43, so its actual local membership relation "
            "would give a nonzero vector in the displayed modular kernel. The "
            "full-rank local normal-form image therefore proves that this QQ "
            "principal-node condition is injective on the stated ambient."
        ),
    },
    "boundary": (
        "This rejects the displayed 54-column endpoint ambient at one actual "
        "resolved node already. It neither supplies a finite quotient of "
        "R/(t^17) nor evaluates the enlarged ambients, other nodes, overlap "
        "maps, smooth/E8 blocks, a q8 kernel, or a child equation."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8E743LOCALGOODREDUCTION|prime=43|ambient=54|rank=54|kernel=0|"
    "status=PASS_EXACT_Q8_E7_4_3_PRINCIPAL_NODE_INJECTIVITY",
    flush=True,
)
