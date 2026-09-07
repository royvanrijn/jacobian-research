#!/usr/bin/env sage
"""Build the exact handoff for the alternate Q80 final q6 divisor.

status: ACTIVE_PROOF
claim: replay the equation-open alternate final-divisor and MW17 lattice gate
inputs: pinned Q80 neighbour, chamber, nefness, and frame certificates
outputs: elkies-k3-q80-alternate-final-divisor-handoff-v1.json

This is the fail-closed equation-compilation input.  It packages the certified
nef divisor and its physical zero, the complete determinant-one NS transport,
and the rootless child frame.  It also records, explicitly, that no generic
characteristic-zero equation for the immediate A1/MW16 parent is currently
available.  Thus an equation compiler cannot silently substitute the CM24
specialization or mistake the abstract Riemann--Roch statement for explicit
resolved functions.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    Integer,
    ZZ,
    block_diagonal_matrix,
    gcd,
    load,
    matrix,
    pari,
    vector,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/"
    / "elkies-k3-q80-alternate-final-divisor-handoff-v1.json"
)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--check", action="store_true")
parser.add_argument("--write-artifact", action="store_true")
arguments = parser.parse_args()

load(str(HERE / "verify_q80_alternate_final_q6_nef.sage"))


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_matrix_text(value):
    return "\n".join(" ".join(map(str, row)) for row in value.rows()) + "\n"


def reflect(divisor, curve):
    return divisor + intersection(divisor, curve, alternate_ns) * curve


def physical_representative(divisor):
    lookup = dict(curves)
    for name, _ in reflections:
        divisor = reflect(divisor, lookup[name])
    return vector(ZZ, divisor)


rootless_child, q6_transport = neighbor(
    alternate_child, ZZ(6), ZZ(2), ZZ(3), q6_v
)
rootless_ns = block_diagonal_matrix(U, -rootless_child)
assert q6_transport * alternate_ns * q6_transport.transpose() == rootless_ns
assert abs(q6_transport.det()) == 1

# Rows of the transport express the new standard U+(-MW) basis in the old
# A1/MW16 basis.  Weyl-reduce both members of U by the same reflections used
# to turn the raw q6 class into the certified nef physical fibre.
child_fiber = vector(ZZ, [1, 0] + [0] * 17)
child_zero = vector(ZZ, [-1, 1] + [0] * 17)
raw_fiber_from_transport = vector(ZZ, child_fiber * q6_transport)
raw_zero_from_transport = vector(ZZ, child_zero * q6_transport)
assert raw_fiber_from_transport == raw
physical_fiber = physical_representative(raw_fiber_from_transport)
physical_zero = physical_representative(raw_zero_from_transport)
assert physical_fiber == reduced
assert physical_fiber * alternate_ns * physical_fiber == 0
assert physical_zero * alternate_ns * physical_zero == -2
assert intersection(physical_fiber, physical_zero, alternate_ns) == 1
assert gcd(list(alternate_ns * physical_fiber)) == 1

# Replay the complete lattice and linear-system gate.
root_result = pari(rootless_child).qfminim(2)
norm_four_result = pari(rootless_child).qfminim(4)
assert int(root_result[0]) == 0
assert int(norm_four_result[0]) == 2626
assert rootless_child.det() == 948
assert rootless_child.is_positive_definite()
assert len(rootless_child.elementary_divisors()) == 17
assert list(rootless_child.elementary_divisors())[-1] == 948
assert reduced == zero + horizontal_section - fiber
assert closest_distance == 3
assert minimum_section_pairing == 1
assert negative_bisection_impossible

q80_to_fourth = (
    fourth_child_transport
    * third_child_transport
    * second_transport
    * first_transport
)
complete_transport = q6_transport * alternate_transport * q80_to_fourth
q80_ns = block_diagonal_matrix(U, -start)
assert complete_transport * q80_ns * complete_transport.transpose() == rootless_ns
assert abs(complete_transport.det()) == 1

published_frame_path = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
published_frame = matrix(
    ZZ,
    [
        list(map(ZZ, line.split()))
        for line in published_frame_path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ],
)
assert pari(rootless_child).qfisom(pari(published_frame)) == 0

candidate_path = (
    ROOT / "artifacts/generated-results/elkies-k3-other-rank17-candidate.json"
)
candidate = json.loads(candidate_path.read_text())
pinned_coordinates = candidate["pinned_ns_coordinates"]
assert pinned_coordinates["transport_determinant"] == 1
assert vector(ZZ, pinned_coordinates["zero"]) == (
    vector(ZZ, pinned_coordinates["isotropic_mate"])
    - vector(ZZ, pinned_coordinates["fibre"])
)
assert candidate["frame"]["sha256"] == hashlib.sha256(
    canonical_matrix_text(rootless_child).encode()
).hexdigest()

source_paths = [
    Path(__file__).resolve(),
    HERE / "verify_q80_alternate_fifth_q6_rootless.sage",
    HERE / "analyze_q80_alternate_final_q6_chamber.sage",
    HERE / "verify_q80_alternate_final_q6_nef.sage",
    ROOT / "elkies-k3/RANK_MUTATION_AND_LIFT_THEOREMS.md",
    ROOT / "artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json",
    ROOT / "artifacts/generated-results/elkies-k3-other-rank17-candidate.json",
    ROOT / "artifacts/generated-results/elkies-k3-other-rank17-invariants.json",
]

payload = {
    "schema": "elkies-k3-q80-alternate-final-divisor-handoff-v1",
    "status": "PASS_LATTICE_HANDOFF_EQUATION_OPEN",
    "coordinate_convention": {
        "parent": "U + (-M_A1_MW16), row vectors",
        "parent_fiber": list(map(int, fiber)),
        "parent_zero": list(map(int, zero)),
        "transport_direction": (
            "rows express the child U+(-MW17) basis in the parent basis"
        ),
    },
    "parent": {
        "frame_type": "A1/MW16",
        "frame_determinant": int(alternate_child.det()),
        "frame_gram": rows(alternate_child),
        "simple_root": list(map(int, simple[0])),
        "simple_component": list(map(int, root_curve)),
        "affine_component": list(map(int, affine_curve)),
        "generic_characteristic_zero_equation": None,
        "equation_status": "NOT_YET_COMPILED",
    },
    "final_q6": {
        "q": 6,
        "a": 2,
        "b": 3,
        "frame_witness": list(map(int, q6_v)),
        "raw_fiber": list(map(int, raw)),
        "weyl_reflections": [[name, int(pairing)] for name, pairing in reflections],
        "physical_nef_fiber": list(map(int, reduced)),
        "raw_zero": list(map(int, raw_zero_from_transport)),
        "physical_zero": list(map(int, physical_zero)),
        "horizontal_section": list(map(int, horizontal_section)),
        "horizontal_section_in_fourth_frame": list(
            map(int, horizontal_section_in_fourth)
        ),
        "decomposition": "F_new = O_old + S - F_old",
        "intersection_checks": {
            "F_new_squared": 0,
            "F_new_primitive": True,
            "F_new_nef": True,
            "F_new_dot_O_new": 1,
            "F_new_dot_O_old": int(intersection(reduced, zero, alternate_ns)),
            "S_squared": -2,
            "S_dot_F_old": 1,
            "S_dot_O_old": 4,
            "minimum_section_pairing": int(minimum_section_pairing),
            "negative_bisection_impossible": True,
        },
    },
    "riemann_roch": {
        "dimension": 2,
        "basis_status": "ABSTRACT_ONLY_UNTIL_PARENT_EQUATION_AND_S_COORDINATES_EXIST",
        "generic_fiber_formal_basis": [
            "1",
            "(y+y(S))/(x-x(S)) after the certified vertical F_old cancellation",
        ],
        "resolved_polar_divisor_proof_status": "OPEN_AT_EQUATION_LEVEL",
        "proof": (
            "primitive nef isotropic divisor on a K3; h0=2 by the certified "
            "Riemann--Roch/elliptic-pencil theorem"
        ),
        "cm24_warning": (
            "the GF(73) five-generator ambient module and its two-dimensional "
            "kernel belong to a special CM24 shadow, not this generic parent"
        ),
    },
    "child": {
        "frame_type": "rootless/MW17",
        "frame_gram": rows(rootless_child),
        "frame_determinant": 948,
        "smith_invariants": list(map(int, rootless_child.elementary_divisors())),
        "root_count": 0,
        "norm_four_vector_count": 2626,
        "norm_four_pairs": 1313,
        "isometric_to_published_R17_frame": False,
        "q6_transport_to_parent": rows(q6_transport),
        "complete_transport_to_q80_ns": rows(complete_transport),
        "complete_transport_determinant_abs": int(abs(complete_transport.det())),
    },
    "pinned_h3_r17_ns_embedding": {
        "coordinate_status": (
            "exact lattice U embedding; use the parent-coordinate Weyl reduction "
            "above for the physical final representative"
        ),
        **pinned_coordinates,
    },
    "equation_frontier": {
        "immediate_parent_equation": "missing",
        "latest_exact_upstream_model": (
            "third-q12 fixed-u=-2 resolved genus-one pencil over a biquadratic "
            "field, descending to Q(sqrt(q1*q2))"
        ),
        "latest_exact_upstream_jacobian": "open",
        "required_route": [
            "compile exact third-q12 Jacobian",
            "compile fourth q12",
            "compile alternate q4 A1/MW16 parent",
            "substitute the final q6 pencil",
        ],
    },
    "source_sha256": {
        str(path.relative_to(ROOT)): sha256_file(path) for path in source_paths
    },
    "reproduce": (
        "sage elkies-k3/scripts/"
        "build_q80_alternate_final_divisor_handoff.sage --write-artifact"
    ),
}

encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
digest = hashlib.sha256(encoded.encode()).hexdigest()

if arguments.write_artifact:
    ARTIFACT.write_text(encoded)
    print(f"Q80ALTHANDOFF|artifact={ARTIFACT}|sha256={digest}|status=PASS_WRITE")
elif arguments.check:
    assert ARTIFACT.read_text() == encoded
    print(f"Q80ALTHANDOFF|artifact={ARTIFACT}|sha256={digest}|status=PASS_CHECK")
else:
    print(encoded, end="")
