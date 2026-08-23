#!/usr/bin/env sage -python
"""
Choose an actual effective zero for the corrected H92 q8 fibration and profile
the q24 divisor relative to every q6 fibre-component zero candidate.

Key point
---------
The pinned D13 zero transported through Weyl isometries need not remain
effective in the OLD q6 fibration.  We therefore do not try to identify it
with a q6 component.

Instead:
  1. transport the pinned q24 FIBRE CLASS through the exact
       pinned-D13 -> dominant -> nef -> component-nef
     lattice isometry;
  2. enumerate actual irreducible q6 fibre components C with C.F_q8 = 1;
  3. use each C as a genuine zero section of the q8 fibration;
  4. construct that D13 frame exactly;
  5. read the q24 divisor's MW projection in that frame;
  6. recover the corresponding degree-one (-2) section P by exact D13 CVP;
  7. certify D = O + P + vertical and record P.O, q6 degree, and vertical data.

This answers which EFFECTIVE q8 zero makes the q24 equation-level hop cheapest,
without identifying the abstract pinned zero with the standard Jacobian zero.

Run:
  sage -python ~/Downloads/audit_h92_q8_q24_effective_zero_choices.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    IntegralLattice, PolynomialRing, QQ, ZZ, block_diagonal_matrix, gcd,
    identity_matrix, lcm, matrix, pari, sage_eval, vector, xgcd
)

Q24_WITNESS = vector(ZZ, (
    0, 5, 0, 1, 2, 1, 2, 2, 2, 2, 4, 8, 2, 0, -1, 1, 1,
))

Q6_REFLECTIONS = (
    1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3, 1,
    7, 6, 5, 4, 2, 3, 4, 5, 6, 7,
)

E6_QF_INDICES = (0, 1, 2, 3, 12, 13)
E6_SIMPLE_IN_QF = (
    (-1, -1, -1, -1, 0, 0),
    (0, 0, 0, 1, 0, 0),
    (0, 0, 1, 0, 0, 0),
    (0, 1, 0, 0, 0, 0),
    (1, 0, 0, 0, 0, 1),
    (2, 1, 0, 0, -1, 1),
)
E6_CARTAN = matrix(ZZ, [
    [2, -1, 0, 0, 0, 0],
    [-1, 2, -1, 0, 0, 0],
    [0, -1, 2, -1, 0, -1],
    [0, 0, -1, 2, -1, 0],
    [0, 0, 0, -1, 2, 0],
    [0, 0, -1, 0, 0, 2],
])
E8_CARTAN = matrix(ZZ, [
    [2, 0, -1, 0, 0, 0, 0, 0],
    [0, 2, 0, -1, 0, 0, 0, 0],
    [-1, 0, 2, -1, 0, 0, 0, 0],
    [0, -1, -1, 2, -1, 0, 0, 0],
    [0, 0, 0, -1, 2, -1, 0, 0],
    [0, 0, 0, 0, -1, 2, -1, 0],
    [0, 0, 0, 0, 0, -1, 2, -1],
    [0, 0, 0, 0, 0, 0, -1, 2],
])


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


def reflect(row, gram, root):
    # Force a real copy: Sage vectors may alias when wrapped by vector().
    row = vector(ZZ, list(row))
    root = vector(ZZ, list(root))
    assert root * gram * root == -2
    return row + (row * gram * root) * root


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


def child_frame_with_zero(ns, fibre, zero):
    assert fibre * ns * fibre == 0
    assert zero * ns * zero == -2
    assert zero * ns * fibre == 1
    mate = zero + fibre
    assert mate * ns * mate == 0 and mate * ns * fibre == 1
    complement = matrix(
        ZZ, [list(fibre * ns), list(mate * ns)]
    ).right_kernel_matrix()
    basis = matrix(
        ZZ, [list(fibre), list(mate)] + [list(row) for row in complement.rows()]
    )
    assert abs(basis.det()) == 1
    child = -(complement * ns * complement.transpose())
    assert (
        basis * ns * basis.transpose()
        == block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -child)
    )
    return child, basis


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    if not count:
        return (), matrix(ZZ, 0, gram.nrows()), (0, 0, 1)
    half = [vector(ZZ, col) for col in matrix(ZZ, result[2]).columns()]
    roots = tuple(half + [-r for r in half])
    rb = matrix(ZZ, [list(r) for r in roots]).row_module().basis_matrix()
    rg = rb * gram * rb.transpose()
    return roots, rb, (rb.rank(), count, abs(ZZ(rg.det())))


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    rank = data[0]
    regular = None
    for shift in range(1, 1000):
        candidate = vector(ZZ, [
            (i + 1) ** 2 + shift * (i + 1) + 1
            for i in range(gram.nrows())
        ])
        if all(candidate * root != 0 for root in roots):
            regular = candidate
            break
    assert regular is not None
    positive = [r for r in roots if regular * r > 0]
    pset = {tuple(r) for r in positive}
    simple = [
        r for r in positive
        if not any(tuple(r - left) in pset for left in positive)
    ]
    simple = matrix(ZZ, [list(r) for r in simple])
    assert simple.nrows() == simple.rank() == rank
    cartan = simple * gram * simple.transpose()
    return simple, cartan


def d13_root_adaptation(child):
    unused, root_basis, inv = roots_and_data(child)
    assert inv == (13, 312, 4), inv
    simple, cartan = deterministic_simple_roots(child)
    assert cartan.det() == 4

    smith, left, right = root_basis.smith_form()
    assert smith == left * root_basis * right
    assert tuple(abs(smith[i, i]) for i in range(13)) == (1,) * 13
    completion = right.inverse()
    initial = simple.stack(completion[13:])
    assert abs(initial.det()) == 1

    adapted = initial * child * initial.transpose()
    coupling = adapted[:13, 13:]
    tail = adapted[13:, 13:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling

    scale = ZZ(1)
    for value in height.list():
        scale = lcm(scale, ZZ(QQ(value).denominator()))
    lll = matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram())
    assert abs(lll.det()) == 1

    change = block_diagonal_matrix(identity_matrix(ZZ, 13), lll.transpose())
    basis = change * initial
    adapted = basis * child * basis.transpose()

    root = adapted[:13, :13]
    coupling = adapted[:13, 13:]
    tail = adapted[13:, 13:]
    H = tail - coupling.transpose() * root.inverse() * coupling
    return basis, adapted, H


def highest_root(cartan):
    half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
    roots = [vector(ZZ, r) for r in half]
    roots += [-vector(ZZ, r) for r in half]
    return max(
        (r for r in roots if all(v >= 0 for v in r)),
        key=lambda r: sum(r),
    )


def class_order(dual):
    order = ZZ(1)
    for value in dual:
        order = lcm(order, ZZ(QQ(value).denominator()))
    return order


def d13_min_correction(dual, root):
    order = class_order(dual)
    expected = {ZZ(1): QQ(0), ZZ(2): QQ(1), ZZ(4): QQ(13) / 4}
    assert order in expected
    correction = expected[order]
    raw = QQ(dual * root * dual)
    mod2 = lambda x: QQ(x) - 2 * (QQ(x) / 2).floor()
    assert mod2(raw) == mod2(correction)
    return order, correction


def section_for_mw(adapted, z, cvp_cap):
    root = adapted[:13, :13]
    coupling = adapted[:13, 13:]
    tail = adapted[13:, 13:]
    H = tail - coupling.transpose() * root.inverse() * coupling
    z = vector(ZZ, z)
    height = QQ(z * H * z)

    base = vector(ZZ, [0] * 13 + list(z))
    pairing = vector(QQ, base * adapted[:, :13])
    dual = pairing * root.inverse()
    order, correction = d13_min_correction(dual, root)
    target_norm = height + correction
    assert target_norm in ZZ

    lattice = IntegralLattice(root)
    iterator = lattice.enumerate_close_vectors(-dual)
    chosen = None
    for unused in range(cvp_cap):
        shift = vector(ZZ, next(iterator))
        candidate = base + vector(ZZ, list(shift) + [0] * 4)
        norm = ZZ(candidate * adapted * candidate)
        if norm == target_norm:
            chosen = candidate
            break
        if norm > target_norm:
            break
    assert chosen is not None, (tuple(z), target_norm)

    pole = ZZ((target_norm - 4) / 2)
    a = ZZ((target_norm - 2) / 2)
    section = vector(ZZ, [a, 1] + list(chosen))
    ns = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -adapted)
    assert section * ns * section == -2
    return section, height, order, correction, pole


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--cvp-cap", type=int, default=8192)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts" / "generated-results"
LOCAL = ROOT / "artifacts" / "local" / "elkies-k3"

FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
PINNED = ROOT / "elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt"
ORBITS = GEN / "elkies-k3-h3-q6-q8-orbits.json"
Q6_CHILD = GEN / "elkies-k3-h92-q6-child-jacobian.json"
TARGET = LOCAL / "q8-target-component-nef.json"
q8_child_candidates = [
    LOCAL / "q8-corrected2cover-qq-child.json",
    GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8_CHILD = next(
    (
        path for path in q8_child_candidates
        if path.exists()
        and "branch_quartic" in json.loads(path.read_text()).get("pencil", {})
    ),
    None,
)
if Q8_CHILD is None:
    raise SystemExit(
        "No complete corrected q8 child artifact containing pencil.branch_quartic"
    )

OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q8-q24-effective-zero-choices.json"
)

for path in (FRAME, PINNED, ORBITS, Q6_CHILD, TARGET, Q8_CHILD):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

source_frame = load_gram(FRAME)
pinned = load_gram(PINNED)
data = json.loads(ORBITS.read_text())
q6child = json.loads(Q6_CHILD.read_text())
target = json.loads(TARGET.read_text())
q8child = json.loads(Q8_CHILD.read_text())

assert data["status"] == "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION"
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert target["normalization"]["representative"] == "component-nef"
assert q8child["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"

U2 = matrix(ZZ, ((0, 1), (1, 0)))
source_ns = block_diagonal_matrix(U2, -source_frame)
source_F = vector(ZZ, [1, 0] + [0] * 17)
source_O = vector(ZZ, [-1, 1] + [0] * 17)

# ===========================================================================
# A. Certified q6 coordinate tower and dominant D13 pinning.
# ===========================================================================

B6 = matrix(ZZ, data["q6"]["neighbor_basis_in_source_ns"])
assert abs(B6.det()) == 1
q6_raw_ns = B6 * source_ns * B6.transpose()
q6_raw = -q6_raw_ns[2:, 2:]

root_mw = matrix(ZZ, data["q6"]["root_mw_basis_in_child"])
root_mw_frame = root_mw * q6_raw * root_mw.transpose()
assert root_mw_frame == matrix(ZZ, data["q6"]["root_adapted_gram"])

simple14 = matrix(ZZ, data["q8"]["simple_root_change_in_root_block"])
simple_change = block_diagonal_matrix(simple14, identity_matrix(ZZ, 3))
simple_to_raw = simple_change * root_mw
simple_frame = matrix(ZZ, data["q8"]["simple_frame_gram"])
assert simple_to_raw * q6_raw * simple_to_raw.transpose() == simple_frame

Bsimple = (
    block_diagonal_matrix(identity_matrix(ZZ, 2), simple_to_raw) * B6
)
q6_simple_ns = block_diagonal_matrix(U2, -simple_frame)
assert Bsimple * source_ns * Bsimple.transpose() == q6_simple_ns

nef_record = data["q8"]["nef_representative"]
nef_raw = vector(ZZ, nef_record["witness_q6_child"])
assert vector(ZZ, [4, 2] + list(nef_raw)) * B6 == vector(
    ZZ, nef_record["fiber_source_h3_ns"]
)

nef_simple = vector(
    ZZ,
    list(nef_raw * root_mw.inverse() * simple_change.inverse()),
)
assert tuple(nef_simple[-3:]) == tuple(nef_record["mw_projection"])

dominant_simple = vector(ZZ, list(nef_simple))
bridge = []
while True:
    labels = vector(ZZ, dominant_simple * simple_frame)[:14]
    negative = [i for i, value in enumerate(labels) if value < 0]
    if not negative:
        break
    i = negative[0]
    pairing = ZZ(labels[i])
    dominant_simple[i] -= pairing
    bridge.append((i, pairing))

assert tuple((i + 1, int(p)) for i, p in bridge) == tuple(
    (int(i), int(p)) for i, p in nef_record["to_dominant_reflections"]
)

matches = [
    hit for hit in data["q8"]["d13_mw4_hits"]
    if tuple(hit["witness_simple_frame"]) == tuple(dominant_simple)
]
assert len(matches) == 1
dominant_hit = matches[0]
assert matrix(ZZ, dominant_hit["d13_root_adapted_gram"]) == pinned

B8 = matrix(ZZ, dominant_hit["neighbor_basis_in_q6_ns"])
raw_d13 = matrix(ZZ, dominant_hit["child_frame"])
assert B8 * q6_simple_ns * B8.transpose() == block_diagonal_matrix(U2, -raw_d13)

A13 = matrix(ZZ, dominant_hit["d13_root_adapted_basis_in_child"])
assert A13 * raw_d13 * A13.transpose() == pinned
Bpinned_to_simple = block_diagonal_matrix(identity_matrix(ZZ, 2), A13) * B8
Gpinned = block_diagonal_matrix(U2, -pinned)
assert Bpinned_to_simple * q6_simple_ns * Bpinned_to_simple.transpose() == Gpinned

print(
    "Q24ZERO_COORD|q6=PASS|simple=PASS|bridge=PASS|dominant=PASS|pinned=PASS",
    flush=True,
)

# ===========================================================================
# B. Transport pinned q24 FIBRE CLASS to nef, then physical component-nef.
# ===========================================================================

F_pinned = vector(ZZ, [1, 0] + [0] * 17)
D_pinned = vector(ZZ, [12, 2] + list(Q24_WITNESS))
assert D_pinned * Gpinned * D_pinned == 0
assert D_pinned * Gpinned * F_pinned == 2

F_dom_simple = F_pinned * Bpinned_to_simple
D_dom_simple = D_pinned * Bpinned_to_simple

simple_root_classes = tuple(
    vector(ZZ, [0, 0] + [ZZ(j == i) for j in range(17)])
    for i in range(14)
)

F_nef_simple = vector(ZZ, list(F_dom_simple))
D_nef_simple = vector(ZZ, list(D_dom_simple))
for i, unused_pairing in reversed(bridge):
    root = simple_root_classes[i]
    F_nef_simple = reflect(F_nef_simple, q6_simple_ns, root)
    D_nef_simple = reflect(D_nef_simple, q6_simple_ns, root)

assert F_nef_simple == vector(ZZ, [4, 2] + list(nef_simple))
F_nef_source = F_nef_simple * Bsimple
D_nef_source = D_nef_simple * Bsimple
assert F_nef_source == vector(ZZ, nef_record["fiber_source_h3_ns"])
assert D_nef_source * source_ns * D_nef_source == 0
assert D_nef_source * source_ns * F_nef_source == 2

# Physical q6 E6/E8 roots.
q6_F = vector(ZZ, [3, 2] + [
    0, 0, -1, -1, -1, -1, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
])

old_simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(i == node) for i in range(17)])
    for node in range(15)
)
q6_O = vector(ZZ, list(source_O))
for node1 in reversed(Q6_REFLECTIONS):
    q6_O = reflect(q6_O, source_ns, old_simple[node1 - 1])

q6_orth = matrix(
    ZZ,
    [list(q6_F * source_ns), list((q6_O + q6_F) * source_ns)],
).right_kernel_matrix()
physical_child = -(q6_orth * source_ns * q6_orth.transpose())
qf_basis = matrix(
    ZZ, pari(physical_child).qfminim(2)[2]
).transpose().row_module().basis_matrix()

e6_qf = matrix(ZZ, [list(qf_basis[i]) for i in E6_QF_INDICES])
e6_roots = matrix(
    ZZ, [list(vector(ZZ, row) * e6_qf) for row in E6_SIMPLE_IN_QF]
) * q6_orth
physical_roots = qf_basis * q6_orth
e8_roots = physical_roots[4:12, :]

assert -e6_roots * source_ns * e6_roots.transpose() == E6_CARTAN
assert -e8_roots * source_ns * e8_roots.transpose() == E8_CARTAN

actual_roots = tuple(e6_roots.rows()) + tuple(e8_roots.rows())

F_actual = vector(ZZ, list(F_nef_source))
D_actual = vector(ZZ, list(D_nef_source))
physical_reflections = []

for unused in range(500):
    pairings = [int(F_actual * source_ns * r) for r in actual_roots]
    negative = [i for i, value in enumerate(pairings) if value < 0]
    if not negative:
        break
    i = negative[0]
    r = actual_roots[i]
    p = pairings[i]
    F_actual = reflect(F_actual, source_ns, r)
    D_actual = reflect(D_actual, source_ns, r)
    physical_reflections.append((i, p))
else:
    raise RuntimeError("physical q8 reduction did not terminate")

assert len(physical_reflections) == 102
target_F = vector(ZZ, target["selected_q8"]["source_h3_ns_vector"])
assert F_actual == target_F
assert D_actual * source_ns * D_actual == 0
assert D_actual * source_ns * F_actual == 2

print(
    "Q24ZERO_TRANSPORT|bridge_steps={}|physical_reflections=102|"
    "target=PASS|D_square=0|D_degree=2".format(len(bridge)),
    flush=True,
)

# ===========================================================================
# C. Actual q6 fibre components with q8 degree one.
# ===========================================================================

e6_highest = highest_root(E6_CARTAN)
e8_highest = highest_root(E8_CARTAN)
e6_affine = q6_F - e6_highest * e6_roots
e8_affine = q6_F - e8_highest * e8_roots

components = []
for i, row in enumerate(e6_roots.rows(), 1):
    components.append((f"IV*_E6_{i}", "IV*", vector(ZZ, row)))
components.append(("IV*_affine", "IV*", vector(ZZ, e6_affine)))
for i, row in enumerate(e8_roots.rows(), 1):
    components.append((f"II*_E8_{i}", "II*", vector(ZZ, row)))
components.append(("II*_affine", "II*", vector(ZZ, e8_affine)))

candidates = [
    (name, kind, curve)
    for name, kind, curve in components
    if curve * source_ns * F_actual == 1
]
assert candidates

print(
    "Q24ZERO_CANDIDATES|count={}|names={}".format(
        len(candidates), ",".join(name for name, unused, curve in candidates)
    ),
    flush=True,
)

# Quartic values at the two additive old-base points are a useful anchor
# diagnostic, independent of the later D13 frame choice.
old_R = PolynomialRing(QQ, "T")
old_T = old_R.gen()
UR = PolynomialRing(QQ, "U")
U = UR.gen()
UF = UR.fraction_field()
TR = PolynomialRing(UF, "T")
T = TR.gen()
quartic = TR(sage_eval(q8child["pencil"]["branch_quartic"], locals={"U": U, "T": T}))

quartic_by_kind = {}
for kind in {kind for unused, kind, unused2 in candidates}:
    factor_text = next(
        item["factor"] for item in q6child["finite_fibres"]
        if item["kodaira"] == kind
    )
    factor = old_R(sage_eval(factor_text, locals={"T": old_T}))
    assert factor.degree() == 1
    T0 = QQ(-factor[0] / factor[1])
    value = UF(quartic(UF(T0)))
    square = bool(value.is_square())
    sqrt = None if not square else (UF(0) if not value else value.sqrt())
    quartic_by_kind[kind] = (T0, value, square, sqrt)
    print(
        "Q24ZERO_QUARTIC|kind={}|T0={}|zero={}|square={}|sqrt={}".format(
            kind, T0, int(not value), int(square),
            "NA" if sqrt is None else str(sqrt),
        ),
        flush=True,
    )

# ===========================================================================
# D. Profile q24 relative to each EFFECTIVE zero.
# ===========================================================================

profiles = []
for name, kind, O in candidates:
    child, Bzero = child_frame_with_zero(source_ns, F_actual, O)
    assert roots_and_data(child)[2] == (13, 312, 4)

    A, adapted, H = d13_root_adaptation(child)
    Badapt = block_diagonal_matrix(identity_matrix(ZZ, 2), A) * Bzero
    Gadapt = block_diagonal_matrix(U2, -adapted)
    assert Badapt * source_ns * Badapt.transpose() == Gadapt

    Dcoords = vector(QQ, D_actual) * Badapt.inverse()
    assert all(value in ZZ for value in Dcoords)
    Dcoords = vector(ZZ, Dcoords)
    assert Dcoords[1] == 2

    wD_adapt = vector(ZZ, Dcoords[2:])
    z = vector(ZZ, wD_adapt[-4:])
    Pcoords, height, class_ord, correction, pole = section_for_mw(
        adapted, z, args.cvp_cap
    )

    Psource = vector(ZZ, Pcoords) * Badapt
    assert Psource * source_ns * Psource == -2
    assert Psource * source_ns * F_actual == 1

    vertical = D_actual - O - Psource
    assert vertical * source_ns * F_actual == 0

    # Express vertical exactly in F + D13-root basis for this zero choice.
    root_source = matrix(
        ZZ,
        [
            list(vector(ZZ, [0, 0] + list(A.row(i))) * Bzero)
            for i in range(13)
        ],
    )
    vertical_basis = matrix(
        QQ, [list(F_actual)] + [list(row) for row in root_source.rows()]
    )
    vertical_coeffs = vertical_basis.solve_left(vector(QQ, vertical))
    assert all(value in ZZ for value in vertical_coeffs)
    vertical_coeffs = vector(ZZ, vertical_coeffs)

    q6_degree = int(Psource * source_ns * q6_F)
    q6_Oint = int(Psource * source_ns * q6_O)
    old_h3_degree = int(Psource * source_ns * source_F)
    PdotO = int(Psource * source_ns * O)
    assert PdotO == pole

    T0, qval, qsquare, qsqrt = quartic_by_kind[kind]

    record = {
        "zero": name,
        "old_fibre_kind": kind,
        "T0": str(T0),
        "quartic_value": str(qval),
        "quartic_square": qsquare,
        "quartic_sqrt": None if qsqrt is None else str(qsqrt),
        "q24_mw_coordinates": list(map(int, z)),
        "height": str(height),
        "D13_class_order": int(class_ord),
        "D13_correction": str(correction),
        "P_dot_O": PdotO,
        "P_q6_degree": q6_degree,
        "P_q6_oldzero_intersection": q6_Oint,
        "P_H3_degree": old_h3_degree,
        "vertical_fibre_coefficient": int(vertical_coeffs[0]),
        "vertical_root_coefficients": list(map(int, vertical_coeffs[1:])),
        "zero_source_h3_ns": list(map(int, O)),
        "section_source_h3_ns": list(map(int, Psource)),
    }
    profiles.append(record)

    print(
        "Q24ZERO_PROFILE|zero={}|kind={}|mw={}|height={}|"
        "class_order={}|correction={}|PdotO={}|P_q6_degree={}|"
        "P_q6_O={}|P_h3_degree={}|vertical_F={}|"
        "quartic_square={}".format(
            name, kind, ",".join(map(str, z)), height,
            class_ord, correction, PdotO, q6_degree,
            q6_Oint, old_h3_degree, vertical_coeffs[0],
            int(qsquare),
        ),
        flush=True,
    )

profiles.sort(key=lambda r: (
    not r["quartic_square"],
    r["P_q6_degree"],
    r["P_dot_O"],
    r["zero"],
))
best = profiles[0]

print(
    "Q24ZERO_BEST|zero={}|kind={}|P_q6_degree={}|PdotO={}|"
    "quartic_square={}|mw={}".format(
        best["zero"], best["old_fibre_kind"], best["P_q6_degree"],
        best["P_dot_O"], int(best["quartic_square"]),
        ",".join(map(str, best["q24_mw_coordinates"])),
    ),
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-q8-q24-effective-zero-choices.v1",
    "status": "PASS_EXACT_Q24_EFFECTIVE_ZERO_CHOICES",
    "transport": {
        "classifier_bridge_steps": [[i + 1, int(p)] for i, p in bridge],
        "physical_reflection_count": len(physical_reflections),
        "target_endpoint_match": True,
        "q24_divisor_source_h3_ns": list(map(int, D_actual)),
    },
    "candidate_count": len(profiles),
    "profiles": profiles,
    "recommended_zero": best,
    "boundary": (
        "This chooses among actual q6 fibre components as q8 zero sections "
        "and certifies the induced q24 marked section/vertical decomposition "
        "in NS. A quartic sign and explicit Weierstrass point still require "
        "the resolved local branch at the selected component."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("Q24ZERO_RESULT|status=PASS_EXACT_Q24_EFFECTIVE_ZERO_CHOICES")
print(f"OUTPUT|{OUTPUT}")
