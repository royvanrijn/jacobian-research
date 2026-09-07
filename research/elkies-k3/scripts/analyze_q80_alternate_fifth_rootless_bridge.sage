#!/usr/bin/env sage
"""Locate the alternate two-section fifth child relative to the rootless frame.

The alternate q=4 witness produces a generic A1/MW16 child.  It shares its
CM horizontal class with a marked genus-one equation over GF(73), but their
vertical completions are not yet identified.  This checker composes the
chosen witness's integral NS transport with the first four pinned transports,
expresses the known rootless fiber/zero in that child, and transports the
distinguished norm-4, divisibility-4 Humbert-8 class ``h``.  No q6 shell
search is run.
"""

import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, gcd, matrix, pari, vector


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = ROOT / "elkies-k3/data/fibrations"
load(str(HERE / "analyze_q80_fifth_q4_chamber.sage"))


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


search_artifact = (
    ROOT / "artifacts/local/"
    "q80-fifth-q4-low-degree-neighbor-search-gf73-v1.json"
)
search_data = json.loads(search_artifact.read_text())
candidate = search_data["candidates"][0]
alternate_v = vector(ZZ, candidate["v"])
alternate_child, alternate_transport = neighbor(
    fourth_child_frame, ZZ(4), ZZ(2), ZZ(2), alternate_v
)
assert alternate_child.det() == 948
alternate_simple, alternate_positive = deterministic_simple_roots(
    alternate_child
)
assert alternate_simple.nrows() == 1
assert len(alternate_positive) == 1
assert abs((alternate_simple*alternate_child*alternate_simple.transpose()).det()) == 2
assert pari(alternate_child).qfisom(pari(fifth_child_frame)) == 0

q80_to_fourth = (
    fourth_child_transport
    * third_child_transport
    * second_transport
    * first_transport
)
q80_to_alternate = alternate_transport*q80_to_fourth
q80_ns = block_diagonal_matrix(U, -start)
alternate_ns = block_diagonal_matrix(U, -alternate_child)
assert q80_to_alternate*q80_ns*q80_to_alternate.transpose() == alternate_ns
assert abs(q80_to_alternate.det()) == 1

target = load_matrix(ROOT / "elkies-k3/data/lattice/rank17_gram.txt")
target_ns = block_diagonal_matrix(U, -target)
target_to_q80 = load_matrix(
    DATA / "kumar_q80_rootless_target_to_q80_ns_transport.txt"
)
target_to_alternate = (
    target_to_q80*q80_to_alternate.inverse()
).change_ring(ZZ)
assert abs(target_to_alternate.det()) == 1
assert (
    target_to_alternate*alternate_ns*target_to_alternate.transpose()
    == target_ns
)
rootless_fiber_in_alternate = vector(ZZ, target_to_alternate[0])
rootless_zero_in_alternate = vector(ZZ, target_to_alternate[1])
direct_q = (
    rootless_fiber_in_alternate[0]*rootless_fiber_in_alternate[1]
)
assert direct_q > 0
assert (
    vector(ZZ, rootless_fiber_in_alternate[2:])
    * alternate_child
    * vector(ZZ, rootless_fiber_in_alternate[2:])
    == 2*direct_q
)

h_target = vector(ZZ, (
    4, 4, -1, 0, -3, 0, 2, -2, 1, -2, 1, 1, 0, 1, 0, 0,
    -2, -2, 2,
))
assert h_target*target_ns*h_target == -4
h_divisibility = gcd([abs(ZZ(value)) for value in target_ns*h_target])
assert h_divisibility == 4
h_perp = matrix(ZZ, [list(h_target*target_ns)]).right_kernel_matrix()
h_perp_gram = h_perp*target_ns*h_perp.transpose()
assert abs(h_perp_gram.det()) == 237
h_smith = tuple(
    abs(ZZ(h_perp_gram.smith_form()[0][index, index]))
    for index in range(h_perp_gram.nrows())
)
assert h_smith == (1,)*17+(237,)
h_in_alternate = vector(ZZ, h_target*target_to_alternate)
assert h_in_alternate*alternate_ns*h_in_alternate == -4
assert gcd([abs(ZZ(value)) for value in alternate_ns*h_in_alternate]) == 4

print(
    "Q80ALTFIFTHBRIDGE|"
    f"alternate_v={tuple(alternate_v)}|child=A1/MW16|"
    f"target_fiber_in_alternate={tuple(rootless_fiber_in_alternate)}|"
    f"target_zero_in_alternate={tuple(rootless_zero_in_alternate)}|"
    f"direct_neighbor_q={direct_q}|transport_det={target_to_alternate.det()}|"
    "isometric_to_pinned_fifth=0|"
    "status=PASS_EXACT_ROOTLESS_TRANSPORT",
    flush=True,
)

# The first fibers inside the intrinsic H237 complement occur at q=9 in the
# rootless marking.  They are not automatically short from the alternate
# fifth marking, so transport all thirteen exact representatives before
# proposing a geometric continuation from the pair14 equation.
h8_q9_data = json.loads((
    ROOT / "artifacts/generated-results/rank17-h8-orthogonal-q9-fibers.json"
).read_text())
assert h8_q9_data["classes_up_to_sign"] == 13
alternate_q9_rows = []
for row in h8_q9_data["rows"]:
    rootless_q9_fiber = vector(ZZ, [3, 3]+row["v"])
    alternate_coordinates = vector(
        ZZ, rootless_q9_fiber*target_to_alternate
    )
    alternate_q = alternate_coordinates[0]*alternate_coordinates[1]
    assert alternate_q > 0
    assert (
        vector(ZZ, alternate_coordinates[2:])
        * alternate_child
        * vector(ZZ, alternate_coordinates[2:])
        == 2*alternate_q
    )
    alternate_q9_rows.append((alternate_q, tuple(alternate_coordinates), row))
alternate_q9_rows.sort(key=lambda item: (item[0], item[1]))
best_q9 = alternate_q9_rows[0]
print(
    "Q80ALTFIFTHBRIDGE|"
    f"h8_q9_classes={len(alternate_q9_rows)}|"
    f"minimum_q_from_alternate={best_q9[0]}|"
    f"minimum_q_target_v={tuple(best_q9[2]['v'])}|"
    f"minimum_q_target_roots={(best_q9[2]['root_rank'], best_q9[2]['roots'], best_q9[2]['rootdet'])}|"
    f"minimum_q_alternate_coordinates={best_q9[1]}|"
    "status=PASS_H8_Q9_ALTERNATE_DISTANCE",
    flush=True,
)
print(
    "Q80ALTFIFTHBRIDGE|"
    f"h_in_target={tuple(h_target)}|h_in_alternate={tuple(h_in_alternate)}|"
    f"h_square=-4|h_divisibility={h_divisibility}|"
    f"h_perp_det={abs(h_perp_gram.det())}|h_perp_smith={h_smith}|"
    "status=PASS_H8_CLASS_TRANSPORT",
    flush=True,
)
