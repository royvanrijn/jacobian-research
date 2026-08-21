#!/usr/bin/env sage
"""Reduce the alternate final q6 class in its generic A1/MW16 chamber.

This is the equation-planning gate after the selected alternate fifth-q4
lattice class.  It proves exact component-chamber nonnegativity and records
the degree and Mordell--Weil projection of the q6 divisor.  Full nefness is
certified separately by ``verify_q80_alternate_final_q6_nef.sage``.  An
equation-level pencil remains open; the pair14 model is now known to be a
different CM24-only vertical class.
"""

from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, vector


HERE = Path(__file__).resolve().parent
load(str(HERE / "analyze_q80_fifth_q4_chamber.sage"))

alternate_q4_v = vector(ZZ, (
    -9, 8, -11, 10, -4, 0, 5, 1, -6,
    6, 1, -2, -1, -1, 1, 2, 0,
))
alternate_child, alternate_transport = neighbor(
    fourth_child_frame, ZZ(4), ZZ(2), ZZ(2), alternate_q4_v
)
alternate_ns = block_diagonal_matrix(U, -alternate_child)
simple, positive = deterministic_simple_roots(alternate_child)
assert simple.nrows() == 1
root_gram = simple*alternate_child*simple.transpose()
assert root_gram.nrows() == 1 and root_gram[0, 0] == 2

fiber = vector(ZZ, [1, 0]+[0]*17)
zero = vector(ZZ, [-1, 1]+[0]*17)
root_curve = vector(ZZ, [0, 0]+list(simple[0]))
affine_root = highest_roots(alternate_child, simple, positive)[0][1]
affine_curve = fiber-vector(ZZ, [0, 0]+list(affine_root))
curves = (("O", zero), ("R1", root_curve), ("Theta0", affine_curve))

q6_v = vector(ZZ, (
    0, -2, 4, 2, -1, 2, 1, -1, 1,
    0, 1, -1, 1, 0, 0, 0, 0,
))
raw = vector(ZZ, [2, 3]+list(q6_v))
assert raw*alternate_ns*raw == 0
reduced, reflections = chamber_reduce(raw, curves, alternate_ns)
assert reduced*alternate_ns*reduced == 0
assert all(intersection(reduced, curve, alternate_ns) >= 0 for _, curve in curves)

frame_part = vector(QQ, reduced[2:])
root_coordinates = (
    frame_part*alternate_child*simple.transpose()*root_gram.inverse()
)
mw_projection = frame_part-root_coordinates*simple
mw_norm = mw_projection*alternate_child*mw_projection

print(
    "Q80ALTFINALQ6|"
    f"raw={tuple(raw)}|raw_D.F={intersection(raw, fiber, alternate_ns)}|"
    f"raw_D.O={intersection(raw, zero, alternate_ns)}|"
    f"reflections={tuple(reflections)}|reduced={tuple(reduced)}|"
    f"D.F={intersection(reduced, fiber, alternate_ns)}|"
    f"D.O={intersection(reduced, zero, alternate_ns)}|"
    f"component_pairings={tuple((name, intersection(reduced, curve, alternate_ns)) for name, curve in curves)}|"
    f"root_coordinates={tuple(root_coordinates)}|"
    f"mw_norm={mw_norm}|child=rootless/MW17|"
    "status=PASS_COMPONENT_CHAMBER",
    flush=True,
)
