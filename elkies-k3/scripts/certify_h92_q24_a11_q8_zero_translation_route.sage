#!/usr/bin/env sage -python
"""Certify the C10-zero translation needed by the A11 q8 compiler.

The pinned A11 lattice frame and the invariant-quartic Weierstrass equation
do not use the same zero.  The equation zero is the old physical component
C10, whose pinned A11 MW vector is

    a = (1,-2,1,-1,0,0).

Consequently the two horizontal poles of the q8 divisor
``O_pinned + P12 - 2F`` have equation-side MW vectors ``-a`` and ``P12-a``.
This script computes their exact height/component/pole profiles and certifies
the zero-compatible low-trace lattice word

    -a = S5 - 2*S7 - S17 - Qminus.

The three old-curve traces have degree-weighted coefficient sum zero; the
already-constructed A11 equation section Q has translation weight zero.
Only one degree-three and two degree-one traces are needed.  No section ansatz
or Groebner basis is used here; the signed trace points and q8 H0 plane remain
the next construction.
"""

import hashlib
import json
from pathlib import Path

from sage.all import (
    IntegralLattice,
    QQ,
    ZZ,
    block_diagonal_matrix,
    identity_matrix,
    lcm,
    matrix,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"

BRIDGE_SECTION = LOCAL / "q24-a11-bridge-m-section-marked-qq.json"
TARGET_COSET = LOCAL / "q24-a11-target-coset-bridge.json"
MARKING = LOCAL / "q24-a11-equation-marking-orbit64-mod100003.json"
PHYSICAL = LOCAL / "q24-d12-orbit42-i8star-physical-marking-qq.json"
NEIGHBOURS = LOCAL / "q24-downstream-lift/d12-c10a-zero-q6-all.json"
PARENT_FRAME = LOCAL / "q24-downstream-lift/d12-c10a-zero-frame.txt"
CHILD_FRAME = (
    LOCAL
    / "q24-downstream-lift/d12-c10a-zero-q6-frames/"
    "q6-o0064-r11-n132-d12-ad4a027cb197.txt"
)
Q8_LATTICE = GENERATED / "elkies-k3-h3-a11-q8-orbit12-lattice-certificate.json"
POINTED = GENERATED / "elkies-k3-h3-a11-pointed-opposite-mw-candidates.json"
OUTPUT = GENERATED / "elkies-k3-h3-q24-a11-q8-zero-translation-route.json"
INPUTS = (
    BRIDGE_SECTION,
    TARGET_COSET,
    MARKING,
    PHYSICAL,
    NEIGHBOURS,
    PARENT_FRAME,
    CHILD_FRAME,
    Q8_LATTICE,
    POINTED,
)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


bridge_section = json.loads(BRIDGE_SECTION.read_text())
target_coset = json.loads(TARGET_COSET.read_text())
marking = json.loads(MARKING.read_text())
physical = json.loads(PHYSICAL.read_text())
neighbours = json.loads(NEIGHBOURS.read_text())
q8_lattice = json.loads(Q8_LATTICE.read_text())
pointed = json.loads(POINTED.read_text())

assert bridge_section["status"] == "PASS_EXACT_Q24_A11_BRIDGE_M_SECTION_MARKED_QQ"
assert target_coset["status"] == "PASS_EXACT_A11_TARGET_COSET_BRIDGE"
assert marking["status"] == "PASS_Q42_A11_EQUATION_MARKING_ORBIT64_MOD100003"
assert q8_lattice["status"] == "PASS_EXACT_A11_Q8_EQUATION_COST_LATTICE_CERTIFICATE"
assert pointed["status"] == "PASS_EXACT_A11_POINTED_OPPOSITE_MW_PROFILE_ENUMERATION"

# Recover the C10 equation zero directly from the exact physical marking and
# the determinant-one D12-to-A11 transition.
selected_physical = next(
    row
    for row in physical["orientation_candidates"]
    if row["section_meets_physical_components"] == ["C10"]
)
anchor_root_index = selected_physical["abstract_to_physical"].index("C10")
assert anchor_root_index == 6
neighbor = next(row for row in neighbours["neighbors"] if int(row["orbit_index"]) == 64)
transition = block_diagonal_matrix(
    identity_matrix(ZZ, 2), matrix(ZZ, neighbor["child_root_adapted_basis"])
) * matrix(ZZ, neighbor["neighbor_basis"])
assert abs(transition.det()) == 1
anchor_parent = vector(
    ZZ,
    [0, 0] + [(-1 if index == anchor_root_index else 0) for index in range(17)],
)
anchor_child = anchor_parent * transition.inverse().change_ring(ZZ)
a = vector(ZZ, anchor_child[-6:])
assert anchor_child[1] == 1
assert a == vector(ZZ, (1, -2, 1, -1, 0, 0))

M = vector(ZZ, target_coset["selected_bridge"]["mw"])
P = vector(ZZ, target_coset["target"]["MW"])
assert M == vector(ZZ, bridge_section["section"]["pinned_lattice_MW_Abel_Jacobi"])
assert M - a == vector(ZZ, bridge_section["section"]["equation_MW_Abel_Jacobi"])
assert P == vector(ZZ, q8_lattice["selection"]["mw_projection"])

# Compute exact section profiles in the pinned A11 frame.  Translation by
# -a is a fibrewise surface automorphism, so profile(P-a) against the pinned
# zero equals the geometric intersection P.a on the equation model.
frame = load_matrix(CHILD_FRAME)
assert frame == matrix(ZZ, neighbor["child_root_adapted_frame"])
root_rank = 11
root = frame[:root_rank, :root_rank]
coupling = frame[:root_rank, root_rank:]
tail = frame[root_rank:, root_rank:]
height = tail - coupling.transpose() * root.inverse() * coupling
root_lattice = IntegralLattice(root)


def section_profile(values):
    z = vector(ZZ, values)
    h = QQ(z * height * z)
    base = vector(ZZ, [0] * root_rank + list(z))
    dual = vector(QQ, base * frame[:, :root_rank]) * root.inverse()
    iterator = root_lattice.enumerate_close_vectors(-dual)
    minimum = None
    minimizer_count = 0
    for unused in range(100000):
        shift = vector(ZZ, next(iterator))
        lifted = base + vector(ZZ, list(shift) + [0] * 6)
        norm = QQ(lifted * frame * lifted)
        if minimum is None:
            minimum = norm
        elif norm > minimum:
            break
        minimizer_count += 1
    correction = minimum - h
    pole_order = (h + correction - 4) / 2
    class_order = ZZ.one()
    for value in dual:
        class_order = lcm(class_order, ZZ(QQ(value).denominator()))
    assert pole_order in ZZ
    return {
        "MW": [int(value) for value in z],
        "height": str(h),
        "local_correction": str(correction),
        "P_dot_O": int(pole_order),
        "component_class_order": int(class_order),
        "minimum_root_lifts": minimizer_count,
    }


pinned_target_profile = section_profile(P)
pinned_zero_on_equation = section_profile(-a)
target_on_equation = section_profile(P - a)
bridge_on_equation = section_profile(M - a)
zero_minus_bridge = section_profile(-M)
residual = P - M
residual_profile = section_profile(residual)
residual_trace_on_equation = section_profile(residual - a)
assert pinned_target_profile["height"] == "13"
assert pinned_target_profile["local_correction"] == "3"
assert pinned_target_profile["P_dot_O"] == 6
assert pinned_zero_on_equation["height"] == "25/3"
assert pinned_zero_on_equation["local_correction"] == "5/3"
assert pinned_zero_on_equation["P_dot_O"] == 3
assert target_on_equation["height"] == "70/3"
assert target_on_equation["local_correction"] == "8/3"
assert target_on_equation["P_dot_O"] == 11
assert bridge_on_equation["P_dot_O"] == 8
assert zero_minus_bridge["P_dot_O"] == 5
assert residual_profile["P_dot_O"] == 0

# The pointed-opposite candidate Qminus completes the index-five identity
# shell.  Only the old-curve traces carry a zero-translation weight: the
# already-constructed A11 equation section Q has weight zero.  The following
# exact word has shell-degree weight zero and maximum trace degree three.
shell = [
    vector(ZZ, row)
    for row in target_coset["exact_identity_shell"]["MW_vectors_in_equation_order"]
]
qminus = vector(ZZ, (0, -1, 0, 0, 0, 0))
assert any(vector(ZZ, row["mw"]) == qminus for row in pointed["candidates"])
word = -qminus + shell[5] - 2 * shell[7] - shell[17]
assert word == -a
word_term_profiles = {index: section_profile(shell[index]) for index in (5, 7, 17)}
degrees = marking["equation_identity_shell_new_fibre_degrees"]
assert {index: int(degrees[index]) for index in (5, 7, 17)} == {
    5: 3,
    7: 1,
    17: 1,
}
degree_weight = 3 - 2 - 1
assert degree_weight == 0
word_term_equation_profiles = {
    index: section_profile(shell[index] - int(degrees[index]) * a)
    for index in (5, 7, 17)
}
qminus_equation_profile = section_profile(qminus)
s7_equation_mw = shell[7] - a
qminus_s7_pairing_profiles = {
    "Qminus_plus_S7trace": section_profile(qminus + s7_equation_mw),
    "Qminus_minus_S7trace": section_profile(qminus - s7_equation_mw),
}

payload = {
    "schema": "elkies-k3.h3-q24-a11-q8-zero-translation-route.v1",
    "status": "PASS_EXACT_A11_Q8_ZERO_TRANSLATION_ROUTE",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "equation_zero": {
        "physical_component": "C10",
        "abstract_parent_root_index": anchor_root_index,
        "pinned_A11_MW": [int(value) for value in a],
    },
    "q8_divisor": {
        "pinned_formula": "O_pinned + P12 - 2F",
        "pinned_target_MW": [int(value) for value in P],
        "equation_horizontal_MW": {
            "O_pinned": [int(value) for value in -a],
            "P12": [int(value) for value in P - a],
        },
        "equation_formula": "O_pinned(eq:-a) + P12(eq:P12-a) - 2F",
        "horizontal_difference": [int(value) for value in P],
    },
    "profiles": {
        "P12_against_pinned_zero": pinned_target_profile,
        "O_pinned_against_equation_zero": pinned_zero_on_equation,
        "P12_against_equation_zero": target_on_equation,
        "M_against_equation_zero": bridge_on_equation,
        "O_pinned_minus_M_against_equation_zero": zero_minus_bridge,
        "P12_minus_M_against_pinned_zero": residual_profile,
        "degree_weight_one_trace_of_P12_minus_M_against_equation_zero": residual_trace_on_equation,
    },
    "pinned_zero_low_trace_word": {
        "formula": "-a=S5-2*S7-S17-Qminus",
        "pointed_opposite_MW": [int(value) for value in qminus],
        "pointed_opposite_coefficient": -1,
        "pointed_opposite_translation_weight": 0,
        "identity_shell_coefficients": {"5": 1, "7": -2, "17": -1},
        "identity_shell_source_degrees": {"5": 3, "7": 1, "17": 1},
        "identity_shell_equation_profiles": {
            str(index): profile for index, profile in word_term_profiles.items()
        },
        "trace_sections_against_C10_zero": {
            str(index): profile for index, profile in word_term_equation_profiles.items()
        },
        "Qminus_against_C10_zero": qminus_equation_profile,
        "Qminus_S7trace_pairing_profiles": qminus_s7_pairing_profiles,
        "required_traces": [
            {"curve": "S5", "degree": 3, "space": "L(4O)"},
        ],
        "degree_one_inputs": ["S7", "S17"],
        "verified_MW_sum": [int(value) for value in word],
        "word_L1": 5,
        "maximum_trace_degree": 3,
        "old_curve_degree_weight": degree_weight,
    },
    "large_Groebner_required": False,
    "proof_boundary": (
        "Exact lattice/zero-translation and section-profile certificate. It "
        "identifies the two horizontal equation points required by the q8 "
        "divisor and a low-trace construction of the pinned-zero point. It "
        "does not identify which exact pointed-section sign is Qminus without "
        "the pinned-good-reduction pairing, lift O_pinned to QQ, construct the "
        "exact P12 equation point, "
        "the two-dimensional H0 plane, quartic, or 2A5 Jacobian."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A11Q8ZEROTRANSLATION|a=1,-2,1,-1,0,0|"
    "Opinned_equation_PO=3|P12_equation_PO=11|"
    "zero_word=S5-2S7-S17-Qminus|degree_weight=0|max_trace_degree=3|"
    f"status={payload['status']}",
    flush=True,
)
print(f"OUTPUT|{OUTPUT.resolve()}", flush=True)
