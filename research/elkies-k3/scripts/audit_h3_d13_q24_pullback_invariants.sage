#!/usr/bin/env sage -python
"""
Invariant pullback audit for the first rank-growing D13 q24 neighbour.

Purpose
-------
Determine whether the q24 horizontal section P and the full q24 divisor D
are cheap on the preceding q6 (E8+E6) fibration.

Crucially, this calculation does NOT identify the dominant D13 presentation
with the later component-nef equation presentation.  That identification is
unnecessary here.

Any Weyl reflections used to move the q8 fibre between chambers are reflections
in components of the q6 fibration.  Such reflections preserve intersection
with F_q6.  Therefore P.F_q6 and D.F_q6 can be computed in the already-certified
dominant D13 presentation.

Coordinate chain checked exactly:

 pinned D13 NS
   -- stored d13_root_adapted_basis -->
 raw dominant D13 child NS
   -- stored q8 neighbor basis -->
 q6 simple-frame NS
   -- stored simple/root-MW changes -->
 raw q6-child NS
   -- stored q6 neighbor basis -->
 source H3 NS

Every Gram identity is asserted independently.

Run:
  sage -python ~/Downloads/audit_h3_d13_q24_pullback_invariants.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix, vector

Q24_WITNESS = vector(ZZ, (
    0, 5, 0, 1, 2, 1, 2, 2, 2, 2, 4, 8, 2, 0, -1, 1, 1,
))


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents" / "jacobian-research",
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
    ]
    seen = set()
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if (
            (candidate / "elkies-k3" / "scripts").is_dir()
            and (candidate / "artifacts" / "generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
PINNED_D13 = (
    ROOT / "elkies-k3/data/fibrations/"
           "h3_q6_q8_d13_mw4_root_adapted_frame.txt"
)
ORBITS = (
    ROOT / "artifacts/generated-results/"
           "elkies-k3-h3-q6-q8-orbits.json"
)
OUTPUT = (
    args.output.resolve()
    if args.output
    else ROOT / "artifacts/local/elkies-k3/"
                "d13-q24-pullback-invariants.json"
)

for path in (FRAME, PINNED_D13, ORBITS):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

source_frame = load_gram(FRAME)
pinned_d13 = load_gram(PINNED_D13)
data = json.loads(ORBITS.read_text())

assert data["status"] == "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION"
assert source_frame.det() == pinned_d13.det() == 948

U = matrix(ZZ, ((0, 1), (1, 0)))
source_ns = block_diagonal_matrix(U, -source_frame)
source_F = vector(ZZ, [1, 0] + [0] * 17)
source_O = vector(ZZ, [-1, 1] + [0] * 17)

# ---------------------------------------------------------------------------
# 1. q6 source -> raw child -> root/MW -> simple frame.
# ---------------------------------------------------------------------------

B6 = matrix(ZZ, data["q6"]["neighbor_basis_in_source_ns"])
assert B6.nrows() == B6.ncols() == 19
assert abs(B6.det()) == 1

q6_raw_ns = B6 * source_ns * B6.transpose()
assert q6_raw_ns[:2, :2] == U
assert q6_raw_ns[:2, 2:] == matrix(ZZ, 2, 17)
assert q6_raw_ns[2:, :2] == matrix(ZZ, 17, 2)
q6_raw_frame = -q6_raw_ns[2:, 2:]
assert q6_raw_frame.det() == 948

root_mw_basis = matrix(ZZ, data["q6"]["root_mw_basis_in_child"])
assert abs(root_mw_basis.det()) == 1
root_mw_frame = root_mw_basis * q6_raw_frame * root_mw_basis.transpose()
assert root_mw_frame == matrix(ZZ, data["q6"]["root_adapted_gram"])

simple_root_change = matrix(
    ZZ, data["q8"]["simple_root_change_in_root_block"]
)
assert simple_root_change.nrows() == simple_root_change.ncols() == 14
assert abs(simple_root_change.det()) == 1

simple_change = block_diagonal_matrix(
    simple_root_change, identity_matrix(ZZ, 3)
)
simple_to_raw = simple_change * root_mw_basis
assert abs(simple_to_raw.det()) == 1

simple_frame = matrix(ZZ, data["q8"]["simple_frame_gram"])
assert simple_to_raw * q6_raw_frame * simple_to_raw.transpose() == simple_frame

Bsimple = (
    block_diagonal_matrix(identity_matrix(ZZ, 2), simple_to_raw)
    * B6
)
q6_simple_ns = block_diagonal_matrix(U, -simple_frame)
assert abs(Bsimple.det()) == 1
assert Bsimple * source_ns * Bsimple.transpose() == q6_simple_ns

# The q6 fibre is untouched by positive-frame basis changes.
F_q6_simple = vector(ZZ, [1, 0] + [0] * 17)
O_q6_abstract = vector(ZZ, [-1, 1] + [0] * 17)
F_q6_source = F_q6_simple * Bsimple
assert F_q6_source == vector(ZZ, B6.row(0))
assert F_q6_source * source_ns * F_q6_source == 0
assert F_q6_source * source_ns * source_F == 2

print(
    "D13Q24INV_COORD|q6_neighbor=PASS|root_mw=PASS|"
    "simple=PASS|simple_to_source=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# 2. Dominant q8 D13 child -> pinned D13 frame.
# ---------------------------------------------------------------------------

# The classifier explicitly states that the pinned frame is mw4_hits[0].
hits = data["q8"]["d13_mw4_hits"]
assert len(hits) == 2
dominant_hit = hits[0]
assert matrix(ZZ, dominant_hit["d13_root_adapted_gram"]) == pinned_d13

B8 = matrix(ZZ, dominant_hit["neighbor_basis_in_q6_ns"])
raw_d13_frame = matrix(ZZ, dominant_hit["child_frame"])
assert B8.nrows() == B8.ncols() == 19
assert abs(B8.det()) == 1
assert (
    B8 * q6_simple_ns * B8.transpose()
    == block_diagonal_matrix(U, -raw_d13_frame)
)

A13 = matrix(ZZ, dominant_hit["d13_root_adapted_basis_in_child"])
assert A13.nrows() == A13.ncols() == 17
assert abs(A13.det()) == 1
assert A13 * raw_d13_frame * A13.transpose() == pinned_d13

Bpinned_to_simple = (
    block_diagonal_matrix(identity_matrix(ZZ, 2), A13)
    * B8
)
G_d13 = block_diagonal_matrix(U, -pinned_d13)
assert (
    Bpinned_to_simple
    * q6_simple_ns
    * Bpinned_to_simple.transpose()
    == G_d13
)

F8_simple = vector(ZZ, Bpinned_to_simple.row(0))
assert F8_simple == vector(
    ZZ, [4, 2] + dominant_hit["witness_simple_frame"]
)
assert F8_simple * q6_simple_ns * F_q6_simple == 2

print(
    "D13Q24INV_COORD|dominant_neighbor=PASS|"
    "d13_adaptation=PASS|pinned_d13=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# 3. Exact q24 section and divisor in pinned D13 coordinates.
# ---------------------------------------------------------------------------

F8_pinned = vector(ZZ, [1, 0] + [0] * 17)
O8_pinned = vector(ZZ, [-1, 1] + [0] * 17)
P_pinned = vector(ZZ, [23, 1] + list(Q24_WITNESS))
D_pinned = vector(ZZ, [12, 2] + list(Q24_WITNESS))

assert Q24_WITNESS * pinned_d13 * Q24_WITNESS == 48
assert P_pinned * G_d13 * P_pinned == -2
assert P_pinned * G_d13 * F8_pinned == 1
assert P_pinned * G_d13 * O8_pinned == 22
assert D_pinned * G_d13 * D_pinned == 0
assert D_pinned * G_d13 * F8_pinned == 2
assert D_pinned == O8_pinned + P_pinned - 10 * F8_pinned

# Height/component check, independent of the neighbour decomposition.
root = pinned_d13[:13, :13]
coupling = pinned_d13[:13, 13:]
tail = pinned_d13[13:, 13:]
height = tail - coupling.transpose() * root.inverse() * coupling
mw = vector(ZZ, Q24_WITNESS[13:])
assert mw == vector(ZZ, (0, -1, 1, 1))
section_height = mw * height * mw
assert section_height == 47

root_pairing = vector(QQ, Q24_WITNESS) * pinned_d13[:, :13]
dual = root_pairing * root.inverse()
local_correction = dual * root * dual
assert local_correction == 1
assert (section_height + local_correction - 4) / 2 == 22

print(
    "D13Q24INV_SECTION|height=47|D13_correction=1|"
    "PdotO=22|decomposition=D=O+P-10F|status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# 4. Pull P, O8, D into q6 simple coordinates.
# ---------------------------------------------------------------------------

P_simple = P_pinned * Bpinned_to_simple
O8_simple = O8_pinned * Bpinned_to_simple
D_simple = D_pinned * Bpinned_to_simple

assert P_simple * q6_simple_ns * P_simple == -2
assert P_simple * q6_simple_ns * F8_simple == 1
assert O8_simple * q6_simple_ns * O8_simple == -2
assert O8_simple * q6_simple_ns * F8_simple == 1
assert D_simple * q6_simple_ns * D_simple == 0
assert D_simple * q6_simple_ns * F8_simple == 2
assert D_simple == O8_simple + P_simple - 10 * F8_simple

# Since F_q6 is [1,0,...] in this frame, intersection with F_q6 must equal
# the second coordinate.  Check both formulations explicitly.
P_q6_degree = int(P_simple * q6_simple_ns * F_q6_simple)
O8_q6_degree = int(O8_simple * q6_simple_ns * F_q6_simple)
D_q6_degree = int(D_simple * q6_simple_ns * F_q6_simple)

assert P_q6_degree == int(P_simple[1])
assert O8_q6_degree == int(O8_simple[1])
assert D_q6_degree == int(D_simple[1])
assert D_q6_degree == O8_q6_degree + P_q6_degree - 20

P_q6_O = int(P_simple * q6_simple_ns * O_q6_abstract)
D_q6_O = int(D_simple * q6_simple_ns * O_q6_abstract)

# ---------------------------------------------------------------------------
# 5. Also map to original source H3 coordinates as a secondary diagnostic.
# ---------------------------------------------------------------------------

P_source = P_simple * Bsimple
O8_source = O8_simple * Bsimple
D_source = D_simple * Bsimple
F8_source = F8_simple * Bsimple

assert P_source * source_ns * P_source == -2
assert P_source * source_ns * F8_source == 1
assert D_source * source_ns * D_source == 0
assert D_source * source_ns * F8_source == 2
assert D_source == O8_source + P_source - 10 * F8_source

# Cross-check q6 degree after the second coordinate transformation.
assert int(P_source * source_ns * F_q6_source) == P_q6_degree
assert int(O8_source * source_ns * F_q6_source) == O8_q6_degree
assert int(D_source * source_ns * F_q6_source) == D_q6_degree

P_h3_degree = int(P_source * source_ns * source_F)
P_h3_O = int(P_source * source_ns * source_O)
D_h3_degree = int(D_source * source_ns * source_F)
D_h3_O = int(D_source * source_ns * source_O)

print(
    "D13Q24INV_PULLBACK|"
    f"P_q6_degree={P_q6_degree}|P_q6_Oabstract={P_q6_O}|"
    f"O8_q6_degree={O8_q6_degree}|"
    f"D_q6_degree={D_q6_degree}|D_q6_Oabstract={D_q6_O}|"
    f"P_h3_degree={P_h3_degree}|P_h3_O={P_h3_O}|"
    f"D_h3_degree={D_h3_degree}|D_h3_O={D_h3_O}|"
    f"P_q6_simple={','.join(map(str, P_simple))}|"
    f"D_q6_simple={','.join(map(str, D_simple))}",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h3-d13-q24-pullback-invariants.v1",
    "status": "PASS_EXACT_D13_Q24_PULLBACK_INVARIANTS",
    "scope": (
        "Invariant pullback through the certified dominant D13 presentation. "
        "No dominant/component-nef D13 identification is assumed."
    ),
    "q24": {
        "mw_coordinates": [0, -1, 1, 1],
        "height": "47",
        "D13_component_correction": "1",
        "P_dot_O_D13": 22,
        "decomposition": "D=O+P-10F",
    },
    "q6_pullback": {
        "P_degree": P_q6_degree,
        "P_abstract_zero_intersection": P_q6_O,
        "q8_zero_degree": O8_q6_degree,
        "D_degree": D_q6_degree,
        "D_abstract_zero_intersection": D_q6_O,
        "P_simple_coordinates": list(map(int, P_simple)),
        "D_simple_coordinates": list(map(int, D_simple)),
    },
    "source_h3_diagnostic": {
        "P_degree": P_h3_degree,
        "P_zero_intersection": P_h3_O,
        "D_degree": D_h3_degree,
        "D_zero_intersection": D_h3_O,
        "P_source_coordinates": list(map(int, P_source)),
        "D_source_coordinates": list(map(int, D_source)),
    },
    "invariance_reason": (
        "Any later q8 chamber reflections are in q6 fibre components and "
        "therefore preserve intersection with F_q6."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("D13Q24INV_RESULT|status=PASS_EXACT_D13_Q24_PULLBACK_INVARIANTS")
print(f"OUTPUT|{OUTPUT}")
