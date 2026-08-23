#!/usr/bin/env sage -python
"""
Pull the first rank-growing D13 q24 horizontal section back through the
actual component-nef q8 fibration to the H3/q6-child Neron--Severi marking.

The pinned D13 q24 neighbour is

    D = [12,2,w],  w=(0,5,0,1,2,1,2,2,2,2,4,8,2,0,-1,1,1).

Its MW projection is (0,-1,1,1), height 47.  The same frame lift w has
norm 48, so

    P = [23,1,w]

is the corresponding effective D13 section.  Hence exactly

    D = O + P - 10 F,

with no exceptional vertical correction.  This script transports P back
through the *component-nef* q8 fibre actually used by the repaired equation
construction, then reports its degree/intersections in the previous q6 and
original H3 fibrations.

This is a lattice transport audit only; it does not construct P's equation.

Run:
  sage -python ~/Downloads/audit_h3_d13_q24_section_pullback.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix, pari, vector, xgcd
)


REFLECTIONS = (
    1,2,4,3,5,4,2,6,5,4,3,1,7,6,5,4,2,3,4,5,6,7
)
Q24_WITNESS = vector(ZZ, (
    0,5,0,1,2,1,2,2,2,2,4,8,2,0,-1,1,1,
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
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


def bezout_mate(ns, fibre):
    current = ZZ(0)
    entries = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(ns * fibre):
        if not value:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        entries = [left * item for item in entries]
        entries[index] += right
        current = divisor
    assert abs(current) == 1
    if current == -1:
        entries = [-item for item in entries]
    mate = vector(ZZ, entries)
    mate -= (mate * ns * mate // 2) * fibre
    assert mate * ns * mate == 0 and mate * ns * fibre == 1
    return mate


def child_frame(ns, fibre):
    mate = bezout_mate(ns, fibre)
    complement = matrix(
        ZZ, [list(fibre * ns), list(mate * ns)]
    ).right_kernel_matrix()
    child = -(complement * ns * complement.transpose())
    basis = matrix(
        ZZ, [list(fibre), list(mate)] + [list(row) for row in complement.rows()]
    )
    assert abs(basis.det()) == 1
    assert (
        basis * ns * basis.transpose()
        == block_diagonal_matrix(matrix(ZZ, ((0,1),(1,0))), -child)
    )
    return child, basis


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    if not count:
        return (), matrix(ZZ, 0, gram.nrows()), (0,0,1)
    half = [vector(ZZ, column) for column in matrix(ZZ, result[2]).columns()]
    roots = tuple(half + [-root for root in half])
    basis = matrix(ZZ, [list(root) for root in roots]).row_module().basis_matrix()
    root_gram = basis * gram * basis.transpose()
    return roots, basis, (basis.rank(), count, abs(ZZ(root_gram.det())))


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    rank = data[0]
    regular = None
    for shift in range(1,1000):
        candidate = vector(ZZ, [
            (i+1)**2 + shift*(i+1) + 1 for i in range(gram.nrows())
        ])
        if all(candidate * root != 0 for root in roots):
            regular = candidate
            break
    assert regular is not None
    positive = [root for root in roots if regular * root > 0]
    positive_set = {tuple(root) for root in positive}
    simple = [
        root for root in positive
        if not any(tuple(root-left) in positive_set for left in positive)
    ]
    simple = matrix(ZZ, [list(root) for root in simple])
    assert simple.nrows() == simple.rank() == rank
    return simple, simple * gram * simple.transpose()


def d13_root_adaptation(child):
    unused, root_basis, invariants = roots_and_data(child)
    assert invariants == (13,312,4)
    simple, cartan = deterministic_simple_roots(child)
    assert cartan.det() == 4

    smith, left, right = root_basis.smith_form()
    assert smith == left * root_basis * right
    assert tuple(abs(smith[i,i]) for i in range(13)) == (1,)*13
    completion = right.inverse()
    initial = simple.stack(completion[13:])
    assert abs(initial.det()) == 1

    adapted = initial * child * initial.transpose()
    coupling = adapted[:13,13:]
    tail = adapted[13:,13:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling

    scale = ZZ(1)
    for value in height.list():
        scale = scale.lcm(ZZ(QQ(value).denominator()))
    lll = matrix(ZZ, pari((scale*height).change_ring(ZZ)).qflllgram())
    assert abs(lll.det()) == 1
    quotient_change = block_diagonal_matrix(
        identity_matrix(ZZ,13), lll.transpose()
    )
    basis = quotient_change * initial
    adapted = basis * child * basis.transpose()
    return basis, adapted


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
PINNED_D13 = ROOT / "elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt"
ORBITS = ROOT / "artifacts/generated-results/elkies-k3-h3-q6-q8-orbits.json"
OUTPUT = (
    args.output.resolve()
    if args.output
    else ROOT / "artifacts/local/elkies-k3/d13-q24-section-pullback.json"
)

source_frame = load_gram(FRAME)
pinned_d13 = load_gram(PINNED_D13)
assert source_frame.det() == pinned_d13.det() == 948

U = matrix(ZZ, ((0,1),(1,0)))
source_ns = block_diagonal_matrix(U, -source_frame)
source_F = vector(ZZ, [1,0] + [0]*17)
source_O = vector(ZZ, [-1,1] + [0]*17)

# q6 raw fibre/old-zero convention used by the physical-target compiler.
q6_F = vector(ZZ, [3,2] + [
    0,0,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,1,0
])
simple_old = tuple(
    vector(ZZ, [0,0] + [ZZ(i == node) for i in range(17)])
    for node in range(15)
)

def reflect(value, nodes):
    result = vector(ZZ, value)
    for node in nodes:
        root = simple_old[node-1]
        result += (result * source_ns * root) * root
    return result

q6_O = reflect(source_O, tuple(reversed(REFLECTIONS)))
assert q6_F * source_ns * q6_F == 0
assert q6_O * source_ns * q6_O == -2
assert q6_O * source_ns * q6_F == 1

# Build actual E6+E8 q6 fibre roots exactly as the physical-target script.
q6_orth = matrix(
    ZZ, [list(q6_F*source_ns), list((q6_O+q6_F)*source_ns)]
).right_kernel_matrix()
q6_child = -(q6_orth * source_ns * q6_orth.transpose())
qf_basis = matrix(
    ZZ, pari(q6_child).qfminim(2)[2]
).transpose().row_module().basis_matrix()
assert qf_basis.rank() == 14
physical_roots = qf_basis * q6_orth

E6_QF_INDICES = (0,1,2,3,12,13)
E6_SIMPLE_IN_QF = (
    (-1,-1,-1,-1,0,0),
    (0,0,0,1,0,0),
    (0,0,1,0,0,0),
    (0,1,0,0,0,0),
    (1,0,0,0,0,1),
    (2,1,0,0,-1,1),
)
e6_qf = matrix(ZZ, [list(qf_basis[i]) for i in E6_QF_INDICES])
e6_roots = matrix(ZZ, [
    vector(ZZ, row) * e6_qf for row in E6_SIMPLE_IN_QF
]) * q6_orth
e8_roots = physical_roots[4:12,:]
actual_roots = tuple(e6_roots.rows()) + tuple(e8_roots.rows())

# Use the STORED dominant D13 hit whose root-adapted Gram is exactly the
# repository's pinned D13 frame.  This avoids a fresh 17-dimensional qfisom
# search (which can require >1GB of PARI stack).
q8_data = json.loads(ORBITS.read_text())
assert q8_data["status"] == "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION"

dominant_hit = next(
    hit for hit in q8_data["q8"]["d13_mw4_hits"]
    if hit["mw_projection"] == [0,-2,0]
)
assert matrix(ZZ, dominant_hit["d13_root_adapted_gram"]) == pinned_d13

dominant_q8_F = vector(ZZ, dominant_hit["fiber_source_h3_ns"])
assert dominant_q8_F * source_ns * dominant_q8_F == 0
assert dominant_q8_F * source_ns * q6_F == 2

# Reduce the dominant D13 fibre directly in the ACTUAL q6 E6+E8 component
# chamber.  Uniqueness of the Weyl chamber representative makes this the
# same component-nef fibre used by the repaired q8 equation construction.
q8_F = vector(ZZ, dominant_q8_F)
reflections = []
for unused in range(500):
    pairings = [int(q8_F * source_ns * root) for root in actual_roots]
    negative = [i for i,value in enumerate(pairings) if value < 0]
    if not negative:
        break
    i = negative[0]
    pairing = pairings[i]
    q8_F += pairing * actual_roots[i]
    reflections.append((i,pairing))
else:
    raise RuntimeError("dominant-to-component-nef q8 reduction did not terminate")

assert q8_F * source_ns * q8_F == 0
assert q8_F * source_ns * q6_F == 2
assert all(q8_F * source_ns * root >= 0 for root in actual_roots)

# Build the PINNED D13 NS basis in source-H3 coordinates from stored exact
# change-of-basis data:
#
# pinned D13 NS -> raw D13 NS -> q6-child NS -> source H3 NS.
q6_neighbor_basis = matrix(
    ZZ, q8_data["q6"]["neighbor_basis_in_source_ns"]
)
q8_neighbor_basis = matrix(
    ZZ, dominant_hit["neighbor_basis_in_q6_ns"]
)
d13_adapt = matrix(
    ZZ, dominant_hit["d13_root_adapted_basis_in_child"]
)
full_d13_adapt = block_diagonal_matrix(identity_matrix(ZZ,2), d13_adapt)

pinned_basis_dominant = (
    full_d13_adapt * q8_neighbor_basis * q6_neighbor_basis
)
assert abs(pinned_basis_dominant.det()) == 1
assert (
    pinned_basis_dominant
    * source_ns
    * pinned_basis_dominant.transpose()
    == block_diagonal_matrix(U, -pinned_d13)
)
assert vector(ZZ, pinned_basis_dominant.row(0)) == dominant_q8_F

# Replay the SAME ambient Weyl reflections on every basis row.  This maps the
# dominant D13 presentation isometrically to the actual component-nef q8
# presentation while preserving the pinned D13 coordinates.
def reflect_row_in_actual_root(row, root):
    row = vector(ZZ, row)
    return row + (row * source_ns * root) * root

pinned_basis_in_source = matrix(ZZ, pinned_basis_dominant)
running_fibre = vector(ZZ, dominant_q8_F)
for index, pairing in reflections:
    root = actual_roots[index]
    assert running_fibre * source_ns * root == pairing
    running_fibre = reflect_row_in_actual_root(running_fibre, root)
    pinned_basis_in_source = matrix(
        ZZ,
        [
            list(reflect_row_in_actual_root(row, root))
            for row in pinned_basis_in_source.rows()
        ],
    )

assert running_fibre == q8_F
assert abs(pinned_basis_in_source.det()) == 1
assert (
    pinned_basis_in_source
    * source_ns
    * pinned_basis_in_source.transpose()
    == block_diagonal_matrix(U, -pinned_d13)
)
assert vector(ZZ, pinned_basis_in_source.row(0)) == q8_F

# q24 divisor and its exact horizontal section in the PINNED D13 frame.
D_pinned = vector(ZZ, [12,2] + list(Q24_WITNESS))
P_pinned = vector(ZZ, [23,1] + list(Q24_WITNESS))
O_pinned = vector(ZZ, [-1,1] + [0]*17)
F_pinned = vector(ZZ, [1,0] + [0]*17)

d13_ns = block_diagonal_matrix(U, -pinned_d13)
assert Q24_WITNESS * pinned_d13 * Q24_WITNESS == 48
assert P_pinned * d13_ns * P_pinned == -2
assert P_pinned * d13_ns * F_pinned == 1
assert P_pinned * d13_ns * O_pinned == 22
assert D_pinned == O_pinned + P_pinned - 10*F_pinned

# MW height and D13 local correction.
root = pinned_d13[:13,:13]
coupling = pinned_d13[:13,13:]
tail = pinned_d13[13:,13:]
height = tail - coupling.transpose()*root.inverse()*coupling
z = vector(ZZ, Q24_WITNESS[13:])
assert z == vector(ZZ, (0,-1,1,1))
hP = z * height * z
assert hP == 47

pairing = vector(QQ, Q24_WITNESS) * pinned_d13 * matrix(
    ZZ, [list(identity_matrix(ZZ,17).row(i)) for i in range(13)]
).transpose()
dual = pairing * root.inverse()
correction = dual * root * dual
assert correction == 1
assert (hP + correction - 4)/2 == 22

# Pull all classes back to the original H3 NS marking.
P_source = P_pinned * pinned_basis_in_source
D_source = D_pinned * pinned_basis_in_source
Oq8_source = O_pinned * pinned_basis_in_source
Fq8_source = F_pinned * pinned_basis_in_source

assert Fq8_source == q8_F
assert P_source * source_ns * P_source == -2
assert P_source * source_ns * q8_F == 1
assert P_source * source_ns * Oq8_source == 22
assert D_source == Oq8_source + P_source - 10*q8_F

# Degrees in the earlier fibrations.
q6_degree = int(P_source * source_ns * q6_F)
h3_degree = int(P_source * source_ns * source_F)
q6_zero_intersection = int(P_source * source_ns * q6_O)
h3_zero_intersection = int(P_source * source_ns * source_O)

print(
    "D13Q24SECTION|height=47|D13_correction=1|PdotO=22|"
    "decomposition=D=O+P-10F|status=PASS",
    flush=True,
)
print(
    "D13Q24PULLBACK|"
    f"q8_component_reflections={len(reflections)}|"
    f"q6_degree={q6_degree}|q6_zero_intersection={q6_zero_intersection}|"
    f"h3_degree={h3_degree}|h3_zero_intersection={h3_zero_intersection}|"
    f"source_vector={','.join(map(str,P_source))}",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h3-d13-q24-section-pullback.v1",
    "status": "PASS_EXACT_D13_Q24_SECTION_PULLBACK",
    "q24": {
        "mw_coordinates": [0,-1,1,1],
        "height": "47",
        "D13_component_correction": "1",
        "P_dot_O": 22,
        "divisor_decomposition": "D=O+P-10F",
        "section_in_pinned_D13_NS": list(map(int,P_pinned)),
    },
    "q8_component_nef": {
        "reflection_count": len(reflections),
        "fibre_source_h3_ns": list(map(int,q8_F)),
    },
    "pullback": {
        "section_source_h3_ns": list(map(int,P_source)),
        "q6_fibre_degree": q6_degree,
        "q6_zero_intersection": q6_zero_intersection,
        "original_H3_fibre_degree": h3_degree,
        "original_H3_zero_intersection": h3_zero_intersection,
    },
    "boundary": (
        "This is an exact lattice pullback of the q24 horizontal section. "
        "It does not construct its rational equation or the D12 child."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}")
