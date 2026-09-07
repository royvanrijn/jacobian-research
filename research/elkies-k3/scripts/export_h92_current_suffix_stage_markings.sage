#!/usr/bin/env sage -python
"""Export every current suffix stage with all marked targets in local coordinates."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
SUFFIX = GENERATED / "elkies-k3-h3-pinned-r17-current-suffix-marking.json"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = GENERATED / "elkies-k3-h3-current-suffix-stage-markings.json"
FRAMES = GENERATED / "elkies-k3-h3-current-suffix-stage-frames"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


suffix = json.loads(SUFFIX.read_text())
assert suffix["status"] == "PASS_EXACT_PINNED_R17_CURRENT_SUFFIX_MARKING"
pinned = load_matrix(PINNED_FRAME)
g_pinned = block_diagonal_matrix(U2, -pinned)
targets_pinned = {
    name: vector(ZZ, value)
    for name, value in suffix["target_fibres_in_root_adapted_hub"].items()
}
# Rows of this matrix are the equation-A11 basis vectors in pinned-R17
# coordinates.  The current-route A11 stage has the same fibre ray, but its
# full marked basis need not equal the equation-A11 basis.
equation_in_pinned = matrix(
    ZZ, suffix["root_adapted_hub_to_equation_A11_basis"]
)
bases_pinned = {
    name: matrix(ZZ, data["basis_in_pinned_R17"])
    for name, data in suffix["current_suffix_stages"].items()
}
FRAMES.mkdir(parents=True, exist_ok=True)
records = {}
for source_name, source_data in suffix["current_suffix_stages"].items():
    source_in_pinned = bases_pinned[source_name]
    pinned_in_source = source_in_pinned.inverse().change_ring(ZZ)
    g_source = source_in_pinned * g_pinned * source_in_pinned.transpose()
    assert g_source[:2, :2] == U2 and g_source[:2, 2:] == 0
    frame = -g_source[2:, 2:]
    frame_path = FRAMES / f"{source_name}.txt"
    frame_path.write_text(
        f"# exact current suffix frame {source_name}\n"
        + "\n".join(" ".join(map(str, row)) for row in frame.rows()) + "\n"
    )
    targets_source = {
        name: value * pinned_in_source for name, value in targets_pinned.items()
    }
    assert targets_source[source_name] == vector(ZZ, [1, 0] + [0] * 17)
    stage_bases_source = {
        name: basis * pinned_in_source for name, basis in bases_pinned.items()
    }
    assert all(abs(basis.det()) == 1 for basis in stage_bases_source.values())
    equation_to_source = equation_in_pinned * pinned_in_source
    g_equation = equation_in_pinned * g_pinned * equation_in_pinned.transpose()
    assert equation_to_source * g_source * equation_to_source.transpose() == g_equation
    marking_path = GENERATED / f"elkies-k3-h3-{source_name}-marked-frame.json"
    payload = {
        "schema": "elkies-k3.h3-current-suffix-stage-marking.v1",
        "status": "PASS_EXACT_CURRENT_SUFFIX_STAGE_MARKING",
        "hub": source_name,
        "root_data": source_data["root_data"],
        "frame_output": str(frame_path.relative_to(ROOT)),
        "frame_sha256": hashlib.sha256(frame_path.read_bytes()).hexdigest(),
        "equation_A11_to_root_adapted_hub_basis": rows(equation_to_source),
        "root_adapted_hub_to_equation_A11_basis": rows(
            equation_to_source.inverse().change_ring(ZZ)
        ),
        "target_fibres_in_root_adapted_hub": {
            name: list(map(int, value)) for name, value in targets_source.items()
        },
        "current_suffix_stage_bases_in_root_adapted_hub": {
            name: rows(basis) for name, basis in stage_bases_source.items()
        },
        "source_basis_in_pinned_R17": rows(source_in_pinned),
        "pinned_R17_basis_in_source": rows(pinned_in_source),
        "proof_boundary": (
            "Exact full current-suffix marking in the selected local stage basis. "
            "All target fibres and complete stage bases use determinant-one transport."
        ),
        "inputs": {
            "paths": [str(SUFFIX.relative_to(ROOT)), str(PINNED_FRAME.relative_to(ROOT))],
            "sha256": {
                str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in (SUFFIX, PINNED_FRAME)
            },
        },
    }
    marking_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    records[source_name] = {
        "marking": str(marking_path.relative_to(ROOT)),
        "marking_sha256": hashlib.sha256(marking_path.read_bytes()).hexdigest(),
        "frame": str(frame_path.relative_to(ROOT)),
        "frame_sha256": hashlib.sha256(frame_path.read_bytes()).hexdigest(),
        "root_data": source_data["root_data"],
        "pinned_fibre": list(map(int, targets_source["pinned_R17"])),
    }

payload = {
    "schema": "elkies-k3.h3-current-suffix-stage-markings.v1",
    "status": "PASS_EXACT_CURRENT_SUFFIX_STAGE_MARKINGS",
    "records": records,
    "inputs": {
        "paths": [str(SUFFIX.relative_to(ROOT)), str(PINNED_FRAME.relative_to(ROOT))],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (SUFFIX, PINNED_FRAME)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "SUFFIXSTAGEMARK|stages={}|status={}|output={}".format(
        len(records), payload["status"], OUTPUT
    ),
    flush=True,
)
