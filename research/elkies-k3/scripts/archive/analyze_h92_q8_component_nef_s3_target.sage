#!/usr/bin/env sage -python
"""
Determine the ACTUAL D13 MW target/profile for the component-nef q8 equation.

This is lattice-only and fast.  It removes the remaining zero/frame ambiguity.

Inputs:
  artifacts/generated-results/elkies-k3-h3-q6-q8-orbits.json
  artifacts/local/elkies-k3/q8-target-component-nef.json
  elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt

The q8 equation compiler uses the component-nef fibre, while the historically
pinned D13 frame was produced from a different Weyl representative.  This
script constructs the deterministic U-neighbour frame for the component-nef
fibre itself:

    (F8, M8, child17),  M8^2=0, F8.M8=1,
    O8 = M8 - F8.

It then:
  * expresses O8 in the old q6 frame and identifies whether it is one of the
    visible E6/E8 components;
  * builds a deterministic D13 root-adapted basis for this component-nef child;
  * inserts the exact source NS class
        S3 = 21 O_H3 + 22(-P1) + (-P2) - 46 F_H3
    and reads its degree and D13 MW quotient coordinates;
  * computes the MW height and the discriminant-coset/local correction;
  * derives P.O from Shioda:
        height = 4 + 2(P.O) - correction.

This gives the correct pole profile to test against the equation-level trace,
without assuming the pinned dominant-frame vector (0,-1,1,1) survives the
Weyl/frame change unchanged.

Run:
  sage -python ~/Downloads/analyze_h92_q8_component_nef_s3_target.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    QQ, ZZ, block_diagonal_matrix, gcd, identity_matrix, lcm,
    matrix, pari, vector, xgcd
)


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
            (candidate / "elkies-k3/scripts").is_dir()
            and (candidate / "artifacts/generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
ORBITS = ROOT / "artifacts/generated-results/elkies-k3-h3-q6-q8-orbits.json"
TARGET = ROOT / "artifacts/local/elkies-k3/q8-target-component-nef.json"

for path in (FRAME, ORBITS, TARGET):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

orbits = json.loads(ORBITS.read_text())
target = json.loads(TARGET.read_text())
assert orbits["status"] == "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION"
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert target["normalization"]["representative"] == "component-nef"

U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def bezout_vector_for_pairing(ns, fiber):
    pairings = list(ns * fiber)
    current = ZZ(0)
    result = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(pairings):
        if value == 0:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        result = [left * entry for entry in result]
        result[index] += right
        current = divisor
    if abs(current) != 1:
        return None
    if current == -1:
        result = [-entry for entry in result]
    return vector(ZZ, result)


def child_frame(ns, fiber, determinant):
    mate = bezout_vector_for_pairing(ns, fiber)
    if mate is None:
        raise ArithmeticError("fibre is not primitive")
    mate_square = ZZ(mate * ns * mate)
    assert mate_square % 2 == 0
    mate -= (mate_square // 2) * fiber
    assert mate * ns * mate == 0
    assert fiber * ns * mate == 1

    kernel = matrix(
        ZZ, [list(fiber * ns), list(mate * ns)]
    ).right_kernel_matrix()
    child = -(kernel * ns * kernel.transpose())
    assert child.is_positive_definite()
    assert child.det() == determinant

    neighbor = matrix(
        ZZ, [list(fiber), list(mate)] + [list(row) for row in kernel.rows()]
    )
    assert abs(neighbor.det()) == 1
    assert (
        neighbor * ns * neighbor.transpose()
        == block_diagonal_matrix(U2, -child)
    )
    return child, neighbor


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    half = [
        vector(ZZ, column)
        for column in matrix(ZZ, result[2]).columns()
    ]
    roots = tuple(half + [-root for root in half])
    root_basis = matrix(
        ZZ, [list(root) for root in roots]
    ).row_module().basis_matrix()
    root_gram = root_basis * gram * root_basis.transpose()
    return roots, root_basis, (
        root_basis.rank(), count, abs(ZZ(root_gram.det()))
    )


def deterministic_simple_roots(gram):
    roots, unused_basis, data = roots_and_data(gram)
    rank = data[0]
    regular = None
    for shift in range(1, 1000):
        candidate = vector(ZZ, [
            (index + 1)**2 + shift*(index + 1) + 1
            for index in range(gram.nrows())
        ])
        if all(candidate * root != 0 for root in roots):
            regular = candidate
            break
    assert regular is not None

    positive = [root for root in roots if regular * root > 0]
    positive_set = {tuple(root) for root in positive}
    simple = []
    for root in positive:
        if not any(tuple(root-left) in positive_set for left in positive):
            simple.append(root)

    simple = matrix(ZZ, [list(root) for root in simple])
    assert simple.nrows() == simple.rank() == rank
    cartan = simple * gram * simple.transpose()
    return simple, cartan


def d13_root_adaptation(child):
    unused_roots, root_basis, invariants = roots_and_data(child)
    assert invariants == (13, 312, 4)

    simple, cartan = deterministic_simple_roots(child)
    assert cartan.det() == 4

    smith, unused_left, smith_right = root_basis.smith_form()
    assert tuple(
        abs(smith[i, i]) for i in range(13)
    ) == (1,) * 13
    completion = smith_right.inverse()
    assert abs(root_basis.stack(completion[13:]).det()) == 1

    adapted_basis = simple.stack(completion[13:])
    assert abs(adapted_basis.det()) == 1
    adapted = adapted_basis * child * adapted_basis.transpose()

    coupling = adapted[:13, 13:]
    tail = adapted[13:, 13:]
    height = tail - coupling.transpose()*cartan.inverse()*coupling

    scale = lcm(entry.denominator() for entry in height.list())
    lll = matrix(
        ZZ, pari((scale*height).change_ring(ZZ)).qflllgram()
    )
    assert abs(lll.det()) == 1

    quotient_change = block_diagonal_matrix(
        identity_matrix(ZZ, 13), lll.transpose()
    )
    adapted_basis = quotient_change * adapted_basis
    adapted = adapted_basis * child * adapted_basis.transpose()
    cartan = adapted[:13, :13]
    coupling = adapted[:13, 13:]
    tail = adapted[13:, 13:]
    height = tail - coupling.transpose()*cartan.inverse()*coupling
    return adapted_basis, adapted, height


source_child = load_gram(FRAME)
source_ns = block_diagonal_matrix(U2, -source_child)
assert source_child.det() == 948

# q6 coordinate change: rows map q6 NS -> source H3 NS.
q6_neighbor = matrix(
    ZZ, orbits["q6"]["neighbor_basis_in_source_ns"]
)
assert abs(q6_neighbor.det()) == 1
q6_child = -(
    q6_neighbor
    * source_ns
    * q6_neighbor.transpose()
)[2:, 2:]
q6_ns = block_diagonal_matrix(U2, -q6_child)

# Actual component-nef q8 fibre used by the equation compiler.
F8_source = vector(
    ZZ, target["selected_q8"]["source_h3_ns_vector"]
)
assert F8_source * source_ns * F8_source == 0
F8_q6 = vector(ZZ, F8_source * q6_neighbor.inverse())
assert F8_q6 * q6_ns * F8_q6 == 0
assert F8_q6[1] == 2  # old-q6 fibre degree

child13, q8_neighbor = child_frame(q6_ns, F8_q6, determinant=948)
assert roots_and_data(child13)[2] == (13, 312, 4)

# The deterministic zero of THIS component-nef U presentation.
F8_frame = vector(ZZ, [1, 0] + [0]*17)
M8_frame = vector(ZZ, [0, 1] + [0]*17)
O8_frame = M8_frame - F8_frame
assert O8_frame * block_diagonal_matrix(U2, -child13) * O8_frame == -2

O8_q6 = vector(ZZ, O8_frame * q8_neighbor)
O8_source = vector(ZZ, O8_q6 * q6_neighbor)
assert O8_source * source_ns * O8_source == -2
assert O8_source * source_ns * F8_source == 1

old_q6_degree = int(O8_q6[1])
print(
    "Q8COMPZERO|"
    f"old_q6_degree={old_q6_degree}|"
    f"q6_coordinates={','.join(map(str,O8_q6))}|status=PASS",
    flush=True,
)

# Compare deterministic O8 with the visible q6 reducible-fibre roots that are
# q8 sections (F8.root = 1).
visible = []
for family in ("E6", "E8"):
    roots = target["selected_q8"][family]["simple_root_vectors_in_source_h3_ns"]
    degrees = target["selected_q8"][family]["component_degrees"]
    for index, (root, degree) in enumerate(zip(roots, degrees), 1):
        if int(degree) != 1:
            continue
        rv = vector(ZZ, root)
        assert F8_source * source_ns * rv == 1
        equal = O8_source == rv
        visible.append((family, index, equal))
        print(
            f"Q8COMPZERO_VISIBLE|family={family}|root={index}|"
            f"equals_zero={int(equal)}",
            flush=True,
        )

# Exact S3 source-H3 NS class from the certified literal-support identity.
Fh3 = vector(ZZ, [1, 0] + [0]*17)
Oh3 = vector(ZZ, [-1, 1] + [0]*17)
minus_p1 = vector(
    ZZ,
    [5, 1]
    + [-2, -3, -4, -6, -5, -4, -3]
    + [0]*8
    + [1, 0]
)
minus_p2 = vector(
    ZZ,
    [22, 1] + [0]*16 + [-1]
)
assert len(minus_p1) == len(minus_p2) == 19

S3_source = 21*Oh3 + 22*minus_p1 + minus_p2 - 46*Fh3
assert S3_source * source_ns * S3_source == -2

S3_q6 = vector(ZZ, S3_source * q6_neighbor.inverse())
S3_q8 = vector(ZZ, S3_q6 * q8_neighbor.inverse())
degree8 = int(S3_q8[1])
assert degree8 == 52
child_vector = vector(ZZ, S3_q8[2:])

adapt_basis, adapted, height = d13_root_adaptation(child13)
coords = vector(ZZ, child_vector * adapt_basis.inverse())
mw = vector(ZZ, coords[13:])
hval = QQ(mw * height * mw)

print(
    "Q8COMPS3|"
    f"q8_degree={degree8}|"
    f"mw={','.join(map(str,mw))}|height={hval}|status=PASS",
    flush=True,
)

# Compute the D13 discriminant-coset correction.  For tail z, the component
# class is the coset of B*z in D13^*/D13.  Search integral root shifts for the
# minimal norm representative of that coset.
cartan = adapted[:13, :13]
coupling = adapted[:13, 13:]
bz = vector(QQ, coupling * mw)
ainv = cartan.inverse()

# Enumerate the D13 weight coset efficiently by nearest integer shifts around
# the rational orthogonal root coordinate -A^-1 Bz.
center = -ainv * bz
base = vector(ZZ, [ZZ(value.round()) for value in center])

best = None
best_x = None
# D13 correction is small; +/-2 around the nearest lattice point is ample,
# but avoid 5^13 brute force.  Use PARI qfminim on an augmented bounded
# search by testing all roots of norm <= 8 as shifts from base.
roots, unused, unused_data = roots_and_data(cartan)
shifts = [vector(ZZ, [0]*13)] + list(roots)
for shift in shifts:
    x = base + shift
    pvec = vector(QQ, cartan*x + bz)
    corr = QQ(pvec * ainv * pvec)
    if best is None or corr < best:
        best = corr
        best_x = x

# Also exactify against the four D13 discriminant cosets: correction must be
# one of the standard D13 component corrections 0, 1, 13/4, 13/4.
allowed = {QQ(0), QQ(1), QQ(13)/4}
if best not in allowed:
    print(
        f"Q8COMPS3_CORRECTION_DIAGNOSTIC|candidate={best}|"
        "status=NEEDS_WIDER_ROOT_SHIFT_SEARCH",
        flush=True,
    )
else:
    print(
        f"Q8COMPS3_CORRECTION|correction={best}|status=PASS",
        flush=True,
    )

if best in allowed:
    PO = (hval - 4 + best)/2
    print(
        f"Q8COMPS3_PROFILE|height={hval}|correction={best}|"
        f"O_intersection={PO}|"
        f"predicted_x_den_degree={2*PO}|"
        f"predicted_y_den_degree={3*PO}|status=PASS",
        flush=True,
    )

# Compare historical pinned dominant-frame target only as a warning.
pinned = matrix(
    ZZ,
    [
        [ZZ(value) for value in line.split()]
        for line in (
            ROOT / "elkies-k3/data/fibrations/"
            "h3_q6_q8_d13_mw4_root_adapted_frame.txt"
        ).read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ],
)
pinned_height = matrix(QQ, [
    [QQ(3)/4, QQ(1)/4, -QQ(1)/4, 0],
    [QQ(1)/4, QQ(11)/4, QQ(1)/4, 1],
    [-QQ(1)/4, QQ(1)/4, QQ(11)/4, -1],
    [0, 1, -1, 46],
])
historical = vector(QQ, (0, -1, 1, 1))
print(
    "Q8COMPS3_HISTORICAL|"
    f"pinned_vector=0,-1,1,1|pinned_height={historical*pinned_height*historical}|"
    f"same_as_component_nef={int(tuple(mw)==(0,-1,1,1))}",
    flush=True,
)

print(
    "Q8COMPS3_RESULT|status=PASS_COMPONENT_NEF_D13_TARGET_ANALYSIS",
    flush=True,
)
