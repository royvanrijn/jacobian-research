#!/usr/bin/env sage -python
"""Refine the component-nef shell by regular Hensel liftability at p=103.

status: ACTIVE_PROOF
claim: exact shell marking refinement and the polynomial construction of P230
output: artifacts/local/elkies-k3/q24-2a5-zero-pole-hensel-refinement-p103.json

This replays only a norm-four lattice shell, 36 stored finite embeddings, and
one finite-field elliptic addition.  No Groebner basis is used.
"""

import hashlib
import json
from itertools import combinations_with_replacement
from pathlib import Path

from sage.all import (
    EllipticCurve, GF, PolynomialRing, QQ, ZZ, block_matrix,
    matrix, vector, zero_matrix,
)


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
FRAME_PATH = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frame.txt"
MARKING_PATH = LOCAL / "q24-a11-to-2a5-q8-equation-marking-qq.json"
WORD_PATH = LOCAL / "q24-2a5-q6o1307-horizontal-word.json"
Q4_WORD_PATH = LOCAL / "q24-2a5-q4o230-horizontal-word.json"
POOL_PATH = LOCAL / "q24-2a5-zero-pole-sections-p103.json"
MATCH_PATH = LOCAL / "q24-2a5-zero-pole-shell-match-p103.json"
HENSEL_PATH = LOCAL / "q24-2a5-p1229-hensel-p103.json"
OUTPUT = LOCAL / "q24-2a5-zero-pole-hensel-refinement-p103.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


frame = matrix(ZZ, [
    list(map(int, line.split()))
    for line in FRAME_PATH.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
gram = block_matrix(ZZ, [
    [matrix(ZZ, [[0, 1], [1, 0]]), zero_matrix(ZZ, 2, 17)],
    [zero_matrix(ZZ, 17, 2), -frame],
])
marking = json.loads(MARKING_PATH.read_text())
word = json.loads(WORD_PATH.read_text())
q4_word = json.loads(Q4_WORD_PATH.read_text())
pool = json.loads(POOL_PATH.read_text())
match = json.loads(MATCH_PATH.read_text())
hensel = json.loads(HENSEL_PATH.read_text())
assert marking["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_EQUATION_MARKING"
assert word["status"] == "PASS_EXACT_Q24_2A5_Q6O1307_LOW_POLE_HORIZONTAL_WORD"
assert q4_word["status"] == "PASS_EXACT_Q24_2A5_Q4O230_LOW_POLE_HORIZONTAL_WORD"
assert pool["status"] == "PASS_BOUNDED_MOD103_ZERO_POLE_SECTION_ENUMERATION"
assert match["status"] == "PASS_EXHAUSTIVE_MOD103_ZERO_POLE_SHELL_EMBEDDINGS_CANONICAL_MARKING"
assert hensel["status"] == "PASS_REGULAR_MOD103_POOL_FILTER_UNIQUE_P1229_BRANCH"

old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_components = [
    vector(ZZ, marking["physical_2A5"]["child_coordinates"][f"old_A11_component_{index}"])
    for index in range(11)
]
chains = marking["physical_2A5"]["chains"]
root_components = [old_components[index] for chain in chains for index in chain]
affine_components = [
    old_fibre - sum((old_components[index] for index in chain), vector(ZZ, 19))
    for chain in chains
]
fibre_components = root_components + affine_components

indexed_classes = match["exact_shell"]["indexed_component_nef_classes"]
assert [row["exact_index"] for row in indexed_classes] == list(range(54))
exact_sections = [vector(ZZ, row["NS_coordinates"]) for row in indexed_classes]
assert len(exact_sections) == 54

trivial = matrix(QQ, [old_fibre, old_components[9]] + root_components)
quotient_columns = trivial.right_kernel().basis_matrix().transpose()


def signature(section):
    return vector(QQ, section) * quotient_columns


p1229 = vector(ZZ, word["sections"]["q6_orbit1229"]["effective_section"])
p230 = vector(ZZ, q4_word["q4_orbit230_horizontal"]["effective_section"])
p1229_index = exact_sections.index(p1229)
assert p1229_index == match["exact_shell"]["P1229_exact_index"] == 24

lift_indices = set(map(int, hensel["hensel"]["full_lift_indices"]))
assert len(lift_indices) == 20
mappings = match["all_complete_mappings_exact_index_to_modular_index"]
p230_signature = signature(p230)
p230_triples = []
for indices in combinations_with_replacement(range(len(exact_sections)), 3):
    candidate_signature = sum(
        (signature(exact_sections[index]) for index in indices),
        vector(QQ, quotient_columns.ncols()),
    )
    if candidate_signature == p230_signature:
        p230_triples.append(indices)
assert (0, 2, 3) in p230_triples
p230_exact_triple = (0, 2, 3)

# The P1229 Hensel obstruction removes five of its six possible modular
# images.  The six complete shell embeddings with image 115 all agree on the
# P230 triple, so no residual marking convention is used here.
compatible_mapping_indices = [
    index for index, mapping in enumerate(mappings)
    if int(mapping[str(p1229_index)]) == 115
]
assert len(compatible_mapping_indices) == 6
p230_modular_images = {
    tuple(int(mappings[index][str(exact_index)]) for exact_index in p230_exact_triple)
    for index in compatible_mapping_indices
}
assert len(p230_modular_images) == 6
regular_p230_modular_images = sorted(
    triple for triple in p230_modular_images if set(triple) <= lift_indices
)
assert regular_p230_modular_images == [(114, 62, 36)]
p230_modular_triple = regular_p230_modular_images[0]
assert p230_modular_triple == (114, 62, 36)

# Add the three selected polynomial branches and normalize the double-pole
# P230 projective coordinates with monic quadratic Z.
F103 = GF(103)
RT = PolynomialRing(F103, "T")
K = RT.fraction_field()
E = EllipticCurve(K, [
    0, 0, 0,
    K(RT(pool["surface_mod_103"]["A_coefficients_low_to_high"])),
    K(RT(pool["surface_mod_103"]["B_coefficients_low_to_high"])),
])


def modular_point(index):
    row = pool["sections"][index]
    return E(
        K(RT(row["X_coefficients_low_to_high"])),
        K(RT(row["Y_coefficients_low_to_high"])),
    )


p230_point = sum((modular_point(index) for index in p230_modular_triple), E(0))
x230, y230 = p230_point.xy()
x_denominator = RT(x230.denominator())
y_denominator = RT(y230.denominator())
Z230 = x_denominator.gcd(x_denominator.derivative()).monic()
assert Z230.degree() == 2
assert x_denominator == Z230**2
assert y_denominator == Z230**3
X230 = RT(x230.numerator())
Y230 = RT(y230.numerator())
assert Y230**2 == X230**3 + RT(E.a4()) * X230 * Z230**4 + RT(E.a6()) * Z230**6


def coefficients(poly, length):
    return [int(poly[index]) for index in range(length)]


payload = {
    "schema": "elkies-k3.q24-2a5-zero-pole-hensel-refinement-p103.v1",
    "status": "PASS_EXACT_SHELL_RELATION_REGULARLY_SELECTED_P230_MOD103_CONSTRUCTION",
    "software": "SageMath 10.9 (conda-forge pinned repository environment)",
    "shell_refinement": {
        "component_nef_exact_classes": 54,
        "regular_hensel_branches_in_coefficient_chart": 20,
        "stored_marked_embeddings": len(mappings),
        "P1229_compatible_embedding_indices": compatible_mapping_indices,
        "P1229_exact_to_modular_index": [24, 115],
        "rank_deficient_modular_points_not_classified": True,
    },
    "low_pole_word_audit": {
        "P230_component_nef_triple_count": len(p230_triples),
        "selected_exact_identity_modulo_trivial_lattice": "P230 = Q0 + Q2 + Q3",
        "P230_selected_triple_exact_indices": list(p230_exact_triple),
        "P230_modular_triples_after_P1229_filter": [
            list(triple) for triple in sorted(p230_modular_images)
        ],
        "P230_unique_all_regular_modular_triple": list(p230_modular_triple),
    },
    "P230_mod103": {
        "construction": "polynomial_section_114 + polynomial_section_62 + polynomial_section_36",
        "X_coefficients_low_to_high": coefficients(X230, 9),
        "Y_coefficients_low_to_high": coefficients(Y230, 13),
        "Z_coefficients_low_to_high": coefficients(Z230, 3),
    },
    "compiler_consequence": {
        "direct_P230_or_simple_pole_solve_required": False,
        "modular_polynomial_branches_to_recover": [114, 62, 36],
        "large_Groebner_required": False,
        "next_gate": "recover polynomial branches 114, 62, and 36 over QQ, then add them",
    },
    "proof_boundary": (
        "The exact lattice proves the P230 triple relation. Among the six shell embeddings "
        "compatible with the uniquely surviving P1229 branch, exactly one image triple "
        "consists entirely of regular Hensel branches. Rank-deficient modular points remain unclassified. "
        "This is a construction aid: exact QQ equations for the three summands and P230 "
        "remain to be reconstructed and verified."
    ),
    "inputs": {
        "paths": [
            str(path.relative_to(ROOT))
            for path in (
                FRAME_PATH, MARKING_PATH, WORD_PATH, Q4_WORD_PATH,
                POOL_PATH, MATCH_PATH, HENSEL_PATH,
            )
        ],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (
                FRAME_PATH, MARKING_PATH, WORD_PATH, Q4_WORD_PATH,
                POOL_PATH, MATCH_PATH, HENSEL_PATH,
            )
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5A5SHELLREFINE|component_nef=54|regular=20|compatible_embeddings=6|"
    "P230=114+62+36|status={}|output={}".format(payload["status"], OUTPUT),
    flush=True,
)
