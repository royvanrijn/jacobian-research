#!/usr/bin/env sage
"""Reduce the fifth q=4 class in the exact generic fourth-child chamber.

This continues the pinned q80-to-rootless path after the second q12 step.
It constructs the generic ``4A1/MW13`` fourth child, reduces the fifth q4
fiber against the chosen zero and every affine/simple old component, and
records its root and Mordell--Weil projections.  It is a bounded planning
certificate for using the explicit degree-one sections on the compact CM24
cubic; it does not claim an equation-level fifth Jacobian.
"""

from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, vector


HERE = Path(__file__).resolve().parent
load(str(HERE / "analyze_q80_fourth_q12_chamber.sage"))


fourth_child_frame, fourth_child_transport = neighbor(
    third_child_frame,
    ZZ(fourth["q"]),
    ZZ(fourth["a"]),
    ZZ(fourth["b"]),
    vector(ZZ, map(ZZ, fourth["v"].split(","))),
)
assert fourth_child_frame.det() == 948
fourth_child_ns = block_diagonal_matrix(U, -fourth_child_frame)
simple_fifth, positive_fifth = deterministic_simple_roots(fourth_child_frame)
assert simple_fifth.nrows() == 4
root_gram_fifth = (
    simple_fifth * fourth_child_frame * simple_fifth.transpose()
)
components_fifth = connected_components(root_gram_fifth)
assert sorted(map(len, components_fifth)) == [1, 1, 1, 1]

fiber_fifth = vector(ZZ, [1, 0] + [0] * 17)
zero_fifth = vector(ZZ, [-1, 1] + [0] * 17)
curves_fifth = [("O", zero_fifth)] + [
    (f"R{index}", vector(ZZ, [0, 0] + list(root)))
    for index, root in enumerate(simple_fifth.rows(), 1)
]
affine_fifth = highest_roots(
    fourth_child_frame, simple_fifth, positive_fifth
)
for component_index, (_, root, _) in enumerate(affine_fifth, 1):
    curves_fifth.append(
        (
            f"Theta0_{component_index}",
            fiber_fifth-vector(ZZ, [0, 0] + list(root)),
        )
    )

fifth = steps[4]
raw_fifth = vector(
    ZZ,
    [ZZ(fifth["a"]), ZZ(fifth["b"])]
    + list(map(ZZ, fifth["v"].split(","))),
)
assert raw_fifth * fourth_child_ns * raw_fifth == 0
reduced_fifth, reflection_sequence_fifth = chamber_reduce(
    raw_fifth, curves_fifth, fourth_child_ns
)
assert reduced_fifth * fourth_child_ns * reduced_fifth == 0
assert all(
    intersection(reduced_fifth, curve, fourth_child_ns) >= 0
    for _, curve in curves_fifth
)

frame_part_fifth = vector(QQ, reduced_fifth[2:])
root_coordinates_fifth = (
    frame_part_fifth
    * fourth_child_frame
    * simple_fifth.transpose()
    * root_gram_fifth.inverse()
)
mw_projection_fifth = (
    frame_part_fifth-root_coordinates_fifth*simple_fifth
)
mw_norm_fifth = (
    mw_projection_fifth
    * fourth_child_frame
    * mw_projection_fifth
)

fifth_child_frame, fifth_child_transport = neighbor(
    fourth_child_frame,
    ZZ(fifth["q"]),
    ZZ(fifth["a"]),
    ZZ(fifth["b"]),
    vector(ZZ, map(ZZ, fifth["v"].split(","))),
)
fifth_child_simple, fifth_child_positive = deterministic_simple_roots(
    fifth_child_frame
)
fifth_child_root_gram = (
    fifth_child_simple
    * fifth_child_frame
    * fifth_child_simple.transpose()
)
fifth_child_roots = (
    fifth_child_simple.nrows(),
    2*len(fifth_child_positive),
    abs(fifth_child_root_gram.det()),
)
assert fifth_child_roots == (1, 2, 2)

print(
    "Q80FIFTHQ4CHAMBER|source=4A1/MW13|q=4|ab=2,2|"
    f"raw={tuple(raw_fifth)}|"
    f"raw_D.F={intersection(raw_fifth, fiber_fifth, fourth_child_ns)}|"
    f"raw_D.O={intersection(raw_fifth, zero_fifth, fourth_child_ns)}|"
    f"reflection_count={len(reflection_sequence_fifth)}|"
    f"reduction={reflection_sequence_fifth}",
    flush=True,
)
print(
    "Q80FIFTHQ4CHAMBER|"
    f"reduced={tuple(reduced_fifth)}|"
    f"D.F={intersection(reduced_fifth, fiber_fifth, fourth_child_ns)}|"
    f"D.O={intersection(reduced_fifth, zero_fifth, fourth_child_ns)}|"
    f"component_pairings={tuple((name, intersection(reduced_fifth, curve, fourth_child_ns)) for name, curve in curves_fifth)}|"
    f"root_coordinates={tuple(root_coordinates_fifth)}|"
    f"mw_projection={tuple(mw_projection_fifth)}|"
    f"mw_norm={mw_norm_fifth}|"
    f"child_roots={fifth_child_roots}|child=A1/MW16|"
    "status=PASS_OLD_COMPONENT_CHAMBER",
    flush=True,
)
