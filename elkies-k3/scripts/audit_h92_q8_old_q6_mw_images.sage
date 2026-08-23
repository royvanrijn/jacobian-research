#!/usr/bin/env sage -python
"""
Profile the three already-certified q6 MW section divisors in the anchored
q8/D13 Mordell-Weil group.

The q8 origin is the actual II*_E8_1 branch-point section certified by
derive_h92_q8_d13_branch_anchor.sage.  IV*_E6_1 is therefore the explicit
height-3/4 generator Gnew.

For each of the three effective q6 MW sections already certified by
certify_h3_q6_weyl_section_transport.sage, this script computes:

  * its degree as a multisection of q8;
  * its Pic^0 / Abel-Jacobi MW coordinates in the anchored D13 frame;
  * the exact D13 height of that rational Jacobian point;
  * the minimal P.O and effective (-2)-section lift for the same MW vector.

Finally it tests whether
    { IV*_E6_1, AJ(S1), AJ(S2), AJ(S3) }
is a unimodular MW basis of the D13 fibration.

This is purely exact NS arithmetic.  It does not yet compute the rational
Weierstrass coordinates of the Abel-Jacobi images of degree>1 multisections.

Run:
  sage -python ~/Downloads/audit_h92_q8_old_q6_mw_images.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    IntegralLattice, QQ, ZZ, block_diagonal_matrix, identity_matrix, lcm,
    matrix, pari, vector, xgcd
)

Q6_REFLECTIONS = (
    1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3, 1,
    7, 6, 5, 4, 2, 3, 4, 5, 6, 7,
)

H3_LIFTS = matrix(ZZ, [
    [-5, -4, -3, 0, 0, 0, 0, 0, 0, 0, 0, -4, 1, 0, -4, 2, -2],
    [-10, -8, -6, 0, 0, 0, 0, 0, 0, 0, 0, -8, 4, 1, -8, 5, -4],
    [-5, -4, -3, 0, 0, 0, 0, 0, 0, 0, 0, -3, 2, 0, -4, 2, -2],
])

OLD_ZERO_ROOT_SHIFTS = matrix(ZZ, [
    [5, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, -1, 4],
    [12, 10, 8, 0, 0, 0, 0, 0, 0, 0, 0, 6, -1, 9],
    [5, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 4],
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
    row = vector(ZZ, list(row))
    root = vector(ZZ, list(root))
    assert root * gram * root == -2
    return row + (row * gram * root) * root


def isotropic_mate(ns, fibre):
    current = ZZ(0)
    data = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(ns * fibre):
        if not value:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        data = [left * entry for entry in data]
        data[index] += right
        current = divisor
    assert abs(current) == 1
    if current == -1:
        data = [-entry for entry in data]
    mate = vector(ZZ, data)
    mate -= (mate * ns * mate // 2) * fibre
    assert mate * ns * mate == 0 and mate * ns * fibre == 1
    return mate


def child_frame_with_zero(ns, fibre, zero):
    assert zero * ns * zero == -2
    assert zero * ns * fibre == 1
    mate = zero + fibre
    orth = matrix(
        ZZ, [list(fibre * ns), list(mate * ns)]
    ).right_kernel_matrix()
    basis = matrix(
        ZZ, [list(fibre), list(mate)] + [list(row) for row in orth.rows()]
    )
    assert abs(basis.det()) == 1
    child = -(orth * ns * orth.transpose())
    U2 = matrix(ZZ, ((0, 1), (1, 0)))
    assert basis * ns * basis.transpose() == block_diagonal_matrix(U2, -child)
    return child, basis


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    if not count:
        return (), matrix(ZZ, 0, gram.nrows()), (0, 0, 1)
    half = [vector(ZZ, c) for c in matrix(ZZ, result[2]).columns()]
    roots = tuple(half + [-r for r in half])
    rb = matrix(ZZ, [list(r) for r in roots]).row_module().basis_matrix()
    rg = rb * gram * rb.transpose()
    return roots, rb, (rb.rank(), count, abs(ZZ(rg.det())))


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    rank = data[0]
    regular = None
    for shift in range(1, 1000):
        candidate = vector(
            ZZ,
            [(i + 1) ** 2 + shift * (i + 1) + 1 for i in range(gram.nrows())],
        )
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
    return simple, simple * gram * simple.transpose()


def d13_root_adaptation(child):
    unused, root_basis, invariants = roots_and_data(child)
    assert invariants == (13, 312, 4), invariants
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
    H = tail - coupling.transpose() * cartan.inverse() * coupling
    scale = ZZ(1)
    for value in H.list():
        scale = lcm(scale, ZZ(QQ(value).denominator()))
    lll = matrix(ZZ, pari((scale * H).change_ring(ZZ)).qflllgram())
    assert abs(lll.det()) == 1

    change = block_diagonal_matrix(identity_matrix(ZZ, 13), lll.transpose())
    basis = change * initial
    adapted = basis * child * basis.transpose()
    root = adapted[:13, :13]
    coupling = adapted[:13, 13:]
    tail = adapted[13:, 13:]
    H = tail - coupling.transpose() * root.inverse() * coupling
    assert H.det() == 237
    return basis, adapted, H


def class_order(dual):
    order = ZZ(1)
    for value in dual:
        order = lcm(order, ZZ(QQ(value).denominator()))
    return order


def d13_correction(dual, root):
    order = class_order(dual)
    expected = {ZZ(1): QQ(0), ZZ(2): QQ(1), ZZ(4): QQ(13)/4}
    assert order in expected
    correction = expected[order]
    raw = QQ(dual * root * dual)
    mod2 = lambda x: QQ(x) - 2 * (QQ(x)/2).floor()
    assert mod2(raw) == mod2(correction)
    return order, correction


def minimal_section_for_mw(adapted, z, cap):
    root = adapted[:13, :13]
    coupling = adapted[:13, 13:]
    tail = adapted[13:, 13:]
    H = tail - coupling.transpose() * root.inverse() * coupling

    z = vector(ZZ, z)
    height = QQ(z * H * z)
    base = vector(ZZ, [0]*13 + list(z))
    pairing = vector(QQ, base * adapted[:, :13])
    dual = pairing * root.inverse()
    order, correction = d13_correction(dual, root)
    target_norm = height + correction
    assert target_norm in ZZ and target_norm >= 0

    L = IntegralLattice(root)
    iterator = L.enumerate_close_vectors(-dual)
    chosen = None
    for unused in range(cap):
        shift = vector(ZZ, next(iterator))
        candidate = base + vector(ZZ, list(shift) + [0]*4)
        norm = ZZ(candidate * adapted * candidate)
        if norm == target_norm:
            chosen = candidate
            break
        if norm > target_norm:
            break
    assert chosen is not None, (tuple(z), height, correction, target_norm)

    pole = QQ(height + correction - 4) / 2
    assert pole in ZZ and pole >= 0
    a = ZZ((target_norm - 2) / 2)
    section = vector(ZZ, [a, 1] + list(chosen))
    U2 = matrix(ZZ, ((0,1),(1,0)))
    ns = block_diagonal_matrix(U2, -adapted)
    assert section * ns * section == -2
    return {
        "height": height,
        "class_order": order,
        "correction": correction,
        "pole": ZZ(pole),
        "section": section,
        "lift": chosen,
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--cvp-cap", type=int, default=8192)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
TARGET = LOCAL / "q8-target-component-nef.json"
BRANCH = LOCAL / "q8-d13-branch-anchor.json"
OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q8-old-q6-mw-images.json"
)

for path in (FRAME, TARGET, BRANCH):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

frame = load_gram(FRAME)
target = json.loads(TARGET.read_text())
branch = json.loads(BRANCH.read_text())
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert branch["status"] == "PASS_EXACT_D13_BRANCH_ANCHOR"

U2 = matrix(ZZ, ((0, 1), (1, 0)))
ns = block_diagonal_matrix(U2, -frame)
source_F = vector(ZZ, [1, 0] + [0]*17)
source_O = vector(ZZ, [-1, 1] + [0]*17)
source_simple = tuple(
    vector(ZZ, [0,0] + [ZZ(i == node) for i in range(17)])
    for node in range(15)
)

# ===========================================================================
# 1. Reconstruct the three exact effective q6 MW section classes.
# ===========================================================================

raw_fibre = vector(ZZ, [3, 2] + [
    0, 0, -1, -1, -1, -1, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
])
reflection_roots = tuple(source_simple[node - 1] for node in Q6_REFLECTIONS)

q6_F = vector(ZZ, list(raw_fibre))
for root in reflection_roots:
    q6_F = reflect(q6_F, ns, root)

raw_zero = vector(ZZ, list(source_O))
for root in reversed(reflection_roots):
    raw_zero = reflect(raw_zero, ns, root)

raw_mate = isotropic_mate(ns, raw_fibre)
raw_orth = matrix(
    ZZ, [list(raw_fibre * ns), list(raw_mate * ns)]
).right_kernel_matrix()
raw_transport = matrix(
    ZZ,
    [list(raw_fibre), list(raw_mate)] + [list(row) for row in raw_orth.rows()],
)
assert abs(raw_transport.det()) == 1
raw_child = -(raw_orth * ns * raw_orth.transpose())
raw_roots = matrix(
    ZZ, pari(raw_child).qfminim(2)[2]
).transpose().row_module().basis_matrix()
assert raw_roots.rank() == 14
raw_zero_coords = vector(ZZ, raw_zero * raw_transport.inverse())
assert raw_zero_coords[1] == 1
raw_zero_mw_lift = vector(ZZ, raw_zero_coords[2:])

q6_sections = []
for h3_lift, shift in zip(H3_LIFTS.rows(), OLD_ZERO_ROOT_SHIFTS.rows()):
    lift = (
        raw_zero_mw_lift
        + vector(ZZ, h3_lift)
        + vector(ZZ, shift) * raw_roots
    )
    norm = ZZ(lift * raw_child * lift)
    pole = norm // 2 - 2
    candidate_raw = (pole + 1) * raw_fibre + raw_mate + lift * raw_orth

    candidate = vector(ZZ, list(candidate_raw))
    for root in reflection_roots:
        candidate = reflect(candidate, ns, root)

    assert candidate * ns * candidate == -2
    assert candidate * ns * q6_F == 1
    q6_sections.append(candidate)

assert [int(S * ns * source_F) for S in q6_sections] == [10, 3, 44]

print(
    "Q8OLDMW_Q6|sections=3|source_degrees=10,3,44|"
    "q6_section_check=PASS",
    flush=True,
)

# ===========================================================================
# 2. Anchored q8/D13 frame: zero=II*_E8_1, Gnew=IV*_E6_1.
# ===========================================================================

F8 = vector(ZZ, target["selected_q8"]["source_h3_ns_vector"])
E6 = matrix(
    ZZ, target["selected_q8"]["E6"]["simple_root_vectors_in_source_h3_ns"]
)
E8 = matrix(
    ZZ, target["selected_q8"]["E8"]["simple_root_vectors_in_source_h3_ns"]
)
O8 = vector(ZZ, E8.row(0))
Gnew_curve = vector(ZZ, E6.row(0))

assert O8 * ns * O8 == -2 and O8 * ns * F8 == 1
assert Gnew_curve * ns * Gnew_curve == -2
assert Gnew_curve * ns * F8 == 1
assert Gnew_curve * ns * O8 == 0

child8, Bzero = child_frame_with_zero(ns, F8, O8)
A13, adapted, H = d13_root_adaptation(child8)
Badapt = block_diagonal_matrix(identity_matrix(ZZ, 2), A13) * Bzero
Gadapt = block_diagonal_matrix(U2, -adapted)
assert Badapt * ns * Badapt.transpose() == Gadapt

def coords(curve):
    result = vector(QQ, curve) * Badapt.inverse()
    assert all(v in ZZ for v in result)
    return vector(ZZ, result)

Gnew_coords = coords(Gnew_curve)
Gnew_z = vector(ZZ, Gnew_coords[-4:])
assert Gnew_coords[1] == 1
assert Gnew_z * H * Gnew_z == QQ(3)/4

# Branch certificate selected deterministic frame should agree on this exact
# orientation, not merely up to sign.
assert list(map(int, Gnew_z)) == branch["iv_sections"]["lattice"]["E6_1_mw"]

print(
    "Q8OLDMW_ANCHOR|zero=II*_E8_1|"
    f"Gnew_mw={','.join(map(str,Gnew_z))}|height=3/4|status=PASS",
    flush=True,
)

# ===========================================================================
# 3. Read old q6 section divisors in the new q8 Jacobian.
# ===========================================================================

records = []
mw_rows = [Gnew_z]
for i, S in enumerate(q6_sections, 1):
    c = coords(S)
    degree = ZZ(c[1])
    assert degree == S * ns * F8

    # For a degree-d divisor S on the q8 generic fibre, [S-d*O8] lies in
    # Pic^0.  Subtracting d*O8 changes only the U coordinates; the positive
    # frame tail, hence the MW quotient z, is unchanged.
    z = vector(ZZ, c[-4:])
    height = QQ(z * H * z)

    prof = minimal_section_for_mw(adapted, z, args.cvp_cap)
    Psource = prof["section"] * Badapt
    assert Psource * ns * Psource == -2
    assert Psource * ns * F8 == 1
    assert int(Psource * ns * O8) == int(prof["pole"])

    q6_degree_of_new_point = int(Psource * ns * q6_F)
    source_degree_of_new_point = int(Psource * ns * source_F)

    record = {
        "index": i,
        "q8_multisection_degree": int(degree),
        "q8_zero_intersection_of_curve": int(S * ns * O8),
        "mw_coordinates": list(map(int, z)),
        "height": str(height),
        "D13_class_order": int(prof["class_order"]),
        "D13_local_correction": str(prof["correction"]),
        "jacobian_section_P_dot_O": int(prof["pole"]),
        "jacobian_section_q6_degree": q6_degree_of_new_point,
        "jacobian_section_source_H3_degree": source_degree_of_new_point,
        "q6_section_source_h3_ns": list(map(int, S)),
        "jacobian_section_source_h3_ns": list(map(int, Psource)),
        "effective_lift_in_anchored_D13_frame": list(map(int, prof["lift"])),
    }
    records.append(record)
    mw_rows.append(z)

    print(
        "Q8OLDMW_IMAGE|"
        f"old_q6_generator={i}|q8_degree={degree}|"
        f"curve_O8={S * ns * O8}|"
        f"mw={','.join(map(str,z))}|height={height}|"
        f"class_order={prof['class_order']}|"
        f"correction={prof['correction']}|"
        f"PdotO={prof['pole']}|"
        f"P_q6_degree={q6_degree_of_new_point}|"
        f"P_h3_degree={source_degree_of_new_point}",
        flush=True,
    )

M = matrix(ZZ, [list(row) for row in mw_rows])
det = ZZ(M.det())
index = abs(det)

height_sub = M * H * M.transpose()
height_sub_det = QQ(height_sub.det())
assert height_sub_det == QQ(index**2) * H.det()

print(
    "Q8OLDMW_BASIS|"
    f"det={det}|index={index}|"
    f"sub_height_det={height_sub_det}|full_height_det={H.det()}|"
    f"unimodular={int(index == 1)}",
    flush=True,
)

# Pairings with the explicit new generator help reveal a simple basis shape.
pairings = [QQ(Gnew_z * H * vector(ZZ, r["mw_coordinates"])) for r in records]
print(
    "Q8OLDMW_PAIRINGS|Gnew_with_old={}".format(
        ",".join(map(str, pairings))
    ),
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-q8-old-q6-mw-images.v1",
    "status": "PASS_EXACT_Q8_OLD_Q6_MW_IMAGES",
    "anchor": {
        "zero": "II*_E8_1",
        "explicit_new_generator": "IV*_E6_1",
        "explicit_new_generator_mw": list(map(int, Gnew_z)),
        "explicit_new_generator_height": "3/4",
        "height_gram": [[str(v) for v in row] for row in H.rows()],
    },
    "old_q6_generators": records,
    "combined_basis": {
        "rows": [list(map(int, row)) for row in mw_rows],
        "determinant": int(det),
        "index": int(index),
        "unimodular": bool(index == 1),
        "height_gram_in_combined_basis":
            [[str(v) for v in row] for row in height_sub.rows()],
        "height_determinant": str(height_sub_det),
        "full_D13_MW_height_determinant": str(H.det()),
        "Gnew_pairings_with_old": [str(v) for v in pairings],
    },
    "boundary": (
        "For old q6 curves of q8 degree >1, this certifies their rational "
        "Pic^0/Jacobian MW class and the corresponding effective q8 section "
        "class in NS. It does not yet evaluate the Abel-Jacobi map on the "
        "q8 quartic to produce rational Weierstrass coordinates."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("Q8OLDMW_RESULT|status=PASS_EXACT_Q8_OLD_Q6_MW_IMAGES")
print(f"OUTPUT|{OUTPUT}")
