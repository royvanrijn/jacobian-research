#!/usr/bin/env sage -python
"""
Certify the equation-level q8 divisor in the exact H3 source NS frame.

Key affine datum already recovered:
    old_zero in standard q6 MW = (2,-1,0)
hence, in the old-zero q6 group,
    O_std       = (-2, 1,0)
    S_q8,std    = (-4,-1,0)  because S_q8 has standard coordinate (-2,-2,0)
    S3_std      = (-2, 1,1)  because S3 has standard coordinate (0,0,1).

All three have trivial IV* component correction because a-b is divisible by 3.
Using the exact old-zero q6 Shioda basis reconstructed by the certified Weyl
transport, build their integral section classes.

The corrected q8 compiler uses the component-nef target only for the exact
vertical E6/E8 cycle; its horizontal part is the STANDARD-Weierstrass divisor
O_std + S_q8,std.  Therefore construct

    F8_eq = O_std + S_q8,std + V_component_nef

with V read from q8-target-component-nef.json.

Certify:
  * O_std.S_q8 = 10 (known q8 marking collision degree);
  * O_std.S3 = 21 (known exact S3 denominator-root degree);
  * F8_eq^2 = 0 and F8_eq.F_q6 = 2;
  * all E6/E8 component degrees equal the physical target certificate;
  * F8_eq.S3_std = 52, matching the exact equation bridge;
  * the orthogonal child has D13 root data (13,312,4).

Run:
  sage -python ~/Downloads/certify_h92_q8_equation_ns_divisor.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    IntegralLattice, QQ, ZZ, block_diagonal_matrix, gcd, identity_matrix,
    lcm, matrix, pari, vector, xgcd
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
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
TARGET = ROOT / "artifacts/local/elkies-k3/q8-target-component-nef.json"
TRANSLATION = ROOT / "artifacts/local/elkies-k3/q6-standard-zero-translation.json"
OUTPUT = (
    args.output.resolve()
    if args.output and args.output.is_absolute()
    else ROOT / (
        args.output
        if args.output
        else Path("artifacts/local/elkies-k3/q8-s3-q24-equation-frame-closure.json")
    )
)

for path in (FRAME, TARGET, TRANSLATION):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

target = json.loads(TARGET.read_text())
translation = json.loads(TRANSLATION.read_text())
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert target["normalization"]["representative"] == "component-nef"
assert translation["status"] == "PASS_EXACT_MW_TRANSLATION_PARAMETER"
assert translation["standard_MW_coordinates"]["old_zero"] == [2,-1,0]
assert translation["translation"]["standard_group_translation_vector"] == [-2,1,0]

U = matrix(ZZ, ((0, 1), (1, 0)))
REFLECTIONS = (
    1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3,
    1, 7, 6, 5, 4, 2, 3, 4, 5, 6, 7,
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
assert H3_LIFTS.dimensions() == (3,17)
assert OLD_ZERO_ROOT_SHIFTS.dimensions() == (3,14)

EXPECTED_HEIGHT = matrix(QQ, [
    [QQ(8)/3, QQ(1)/3, -1],
    [QQ(1)/3, QQ(8)/3, 1],
    [-1, 1, 46],
])


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def reflect(value, root, gram):
    value = vector(ZZ, value)
    root = vector(ZZ, root)
    assert root * gram * root == -2
    return value + (value * gram * root) * root


def isotropic_mate(ns, fiber):
    current = ZZ(0)
    data = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(ns * fiber):
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
    mate -= (mate * ns * mate // 2) * fiber
    assert mate * ns * mate == 0 and mate * ns * fiber == 1
    return mate


def roots_and_data(gram):
    qf = pari(gram).qfminim(2)
    count = ZZ(qf[0])
    if not count:
        return (), matrix(ZZ, 0, gram.nrows()), (0,0,1)
    half = [
        vector(ZZ, col)
        for col in matrix(ZZ, qf[2]).columns()
    ]
    roots = tuple(half + [-r for r in half])
    basis = matrix(
        ZZ, [list(r) for r in roots]
    ).row_module().basis_matrix()
    rg = basis * gram * basis.transpose()
    return roots, basis, (
        basis.rank(), count, abs(ZZ(rg.det()))
    )


frame = load_gram(FRAME)
assert frame.nrows() == 17 and frame.det() == 948
ns = block_diagonal_matrix(U, -frame)
Fh3 = vector(ZZ, [1,0] + [0]*17)
Oold = vector(ZZ, [-1,1] + [0]*17)
simple = tuple(
    vector(ZZ, [0,0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)

raw_fiber = vector(
    ZZ,
    [3,2] + [
        0,0,-1,-1,-1,-1,-1,
        0,0,0,0,0,0,0,0,1,0
    ],
)
reflection_roots = tuple(simple[node-1] for node in REFLECTIONS)


def weyl_transport(value):
    result = vector(ZZ, value)
    for root in reflection_roots:
        result = reflect(result, root, ns)
    return result


def inverse_weyl_transport(value):
    result = vector(ZZ, value)
    for root in reversed(reflection_roots):
        result = reflect(result, root, ns)
    return result


F6 = weyl_transport(raw_fiber)
assert F6 * ns * F6 == 0
assert F6 * ns * Fh3 == 2

# Reconstruct the exact old-zero q6 MW basis sections.
raw_mate = isotropic_mate(ns, raw_fiber)
raw_orth = matrix(
    ZZ, [list(raw_fiber*ns), list(raw_mate*ns)]
).right_kernel_matrix()
raw_transport = matrix(
    ZZ,
    [list(raw_fiber), list(raw_mate)] + [list(r) for r in raw_orth.rows()],
)
assert abs(raw_transport.det()) == 1
raw_child = -(raw_orth * ns * raw_orth.transpose())
raw_roots = matrix(
    ZZ, pari(raw_child).qfminim(2)[2]
).transpose().row_module().basis_matrix()
assert raw_roots.dimensions() == (14,17)

raw_zero = inverse_weyl_transport(Oold)
raw_zero_coords = vector(ZZ, raw_zero * raw_transport.inverse())
raw_zero_mw_lift = vector(ZZ, raw_zero_coords[2:])

basis_sections = []
for h3_lift, shift in zip(H3_LIFTS.rows(), OLD_ZERO_ROOT_SHIFTS.rows()):
    candidate_lift = (
        raw_zero_mw_lift
        + vector(ZZ, h3_lift)
        + vector(ZZ, shift) * raw_roots
    )
    norm = ZZ(candidate_lift * raw_child * candidate_lift)
    assert norm >= 4 and norm % 2 == 0
    pole = norm//2 - 2
    candidate_raw = (
        (pole+1)*raw_fiber
        + raw_mate
        + candidate_lift*raw_orth
    )
    candidate = weyl_transport(candidate_raw)
    assert candidate * ns * candidate == -2
    assert candidate * ns * F6 == 1
    basis_sections.append(candidate)

# Exact q6 Shioda map relative to old O.
nef_mate = isotropic_mate(ns, F6)
nef_orth = matrix(
    ZZ, [list(F6*ns), list(nef_mate*ns)]
).right_kernel_matrix()
nef_child = -(nef_orth * ns * nef_orth.transpose())
root_basis_child = matrix(
    ZZ, pari(nef_child).qfminim(2)[2]
).transpose().row_module().basis_matrix()
roots_source = root_basis_child * nef_orth
root_gram = roots_source * ns * roots_source.transpose()
assert roots_source.nrows() == 14 and abs(root_gram.det()) == 3
projection = (
    matrix(QQ, identity_matrix(ZZ, 19))
    - matrix(QQ, ns)
      * matrix(QQ, roots_source.transpose())
      * matrix(QQ, root_gram.inverse())
      * matrix(QQ, roots_source)
)


def shioda(section):
    horizontal = section - Oold - (section*ns*Oold + 2)*F6
    assert horizontal*ns*F6 == 0 and horizontal*ns*Oold == 0
    return vector(QQ, horizontal) * projection


phis = [shioda(section) for section in basis_sections]
height = matrix(QQ, [
    [-left*ns*right for right in phis]
    for left in phis
])
assert height == EXPECTED_HEIGHT


def iv_correction(mw):
    mw = vector(ZZ, mw)
    return QQ(0) if (mw[0]-mw[1]) % 3 == 0 else QQ(4)/3


def section_from_old_mw(mw):
    mw = vector(ZZ, mw)
    phi = sum((ZZ(mw[i])*phis[i] for i in range(3)), vector(QQ,[0]*19))
    h = QQ(mw*EXPECTED_HEIGHT*mw)
    corr = iv_correction(mw)
    PO = (h - 4 + corr)/2
    assert PO in ZZ and PO >= 0
    if corr != 0:
        raise ArithmeticError(
            f"requested section {tuple(mw)} has nonzero IV* correction; "
            "this direct reconstruction intentionally handles only identity-component vectors"
        )
    value = vector(QQ, Oold) + QQ(PO+2)*vector(QQ,F6) + phi
    assert all(entry in ZZ for entry in value)
    value = vector(ZZ, value)
    assert value*ns*value == -2 and value*ns*F6 == 1
    assert value*ns*Oold == PO
    return value, ZZ(PO), h


# Coordinates in old-zero group.
mw_Ostd = vector(ZZ, (-2,1,0))
mw_marked = vector(ZZ, (-4,-1,0))
mw_S3std = vector(ZZ, (-2,1,1))

Ostd, Ostd_oldO, hOstd = section_from_old_mw(mw_Ostd)
Smarked, Smarked_oldO, hSmarked_old = section_from_old_mw(mw_marked)
S3std, S3_oldO, hS3_old = section_from_old_mw(mw_S3std)

# Standard-zero regressions from exact equation-level denominators.
assert Ostd * ns * Smarked == 10
assert Ostd * ns * S3std == 21
assert (
    vector(QQ, (-2,-2,0))
    * EXPECTED_HEIGHT
    * vector(QQ, (-2,-2,0))
    == 24
)
assert (
    vector(QQ, (0,0,1))
    * EXPECTED_HEIGHT
    * vector(QQ, (0,0,1))
    == 46
)

print(
    "Q8EQNS_SECTIONS|"
    f"Ostd_oldO={Ostd_oldO}|marked_oldO={Smarked_oldO}|S3_oldO={S3_oldO}|"
    "Ostd_marked=10|Ostd_S3=21|status=PASS",
    flush=True,
)

# Exact vertical cycle used by the corrected q8 compiler.
selected = target["selected_q8"]
V = ZZ(selected["vertical_fibre_coefficient"]) * F6

for family in ("E6","E8"):
    roots = selected[family]["simple_root_vectors_in_source_h3_ns"]
    cycle = selected[family]["vertical_cycle"]
    assert len(roots) == len(cycle)
    for coefficient, root in zip(cycle, roots):
        V += ZZ(coefficient) * vector(ZZ, root)

# Equation-level q8 divisor.
F8eq = Ostd + Smarked + V
assert F8eq * ns * F8eq == 0
assert F8eq * ns * F6 == 2
assert gcd(tuple(F8eq)) == 1

# The stored component degrees must replay exactly.
for family in ("E6","E8"):
    roots = selected[family]["simple_root_vectors_in_source_h3_ns"]
    expected = list(map(int, selected[family]["component_degrees"]))
    actual = [
        int(F8eq * ns * vector(ZZ, root))
        for root in roots
    ]
    assert actual == expected

degree_S3 = int(F8eq * ns * S3std)
print(
    "Q8EQNS_FIBRE|square=0|q6_degree=2|primitive=1|"
    f"S3_degree={degree_S3}|expected=52|match={int(degree_S3==52)}|status=PASS",
    flush=True,
)
assert degree_S3 == 52

# Independent child root classification.
mate8 = isotropic_mate(ns, F8eq)
orth8 = matrix(
    ZZ, [list(F8eq*ns), list(mate8*ns)]
).right_kernel_matrix()
child8 = -(orth8 * ns * orth8.transpose())
root_data = roots_and_data(child8)[2]
assert root_data == (13,312,4)

# ===========================================================================
# Exact equation-frame AJ(S3) / q24 closure.
#
# IMPORTANT: use F8eq, not the earlier physical root-target isotropic class.
# The genuine equation zero is the actual II*_E8_1 component, whose stored
# component degree against F8eq is one.
# ===========================================================================

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
    U2 = matrix(ZZ, ((0,1),(1,0)))
    assert basis * ns * basis.transpose() == block_diagonal_matrix(U2, -child)
    return child, basis


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    rank = data[0]
    regular = None
    for shift in range(1, 1000):
        candidate = vector(
            ZZ,
            [(i + 1)**2 + shift*(i + 1) + 1 for i in range(gram.nrows())],
        )
        if all(candidate * root != 0 for root in roots):
            regular = candidate
            break
    assert regular is not None
    positive = [root for root in roots if regular * root > 0]
    pset = {tuple(root) for root in positive}
    simple_roots = [
        root for root in positive
        if not any(tuple(root-left) in pset for left in positive)
    ]
    simple_roots = matrix(ZZ, [list(root) for root in simple_roots])
    assert simple_roots.nrows() == simple_roots.rank() == rank
    return simple_roots, simple_roots * gram * simple_roots.transpose()


def d13_root_adaptation(child):
    unused, root_basis, invariants = roots_and_data(child)
    assert invariants == (13,312,4), invariants
    simple_roots, cartan = deterministic_simple_roots(child)
    assert cartan.det() == 4

    smith, left, right = root_basis.smith_form()
    assert smith == left * root_basis * right
    assert tuple(abs(smith[i,i]) for i in range(13)) == (1,)*13
    completion = right.inverse()
    initial = simple_roots.stack(completion[13:])
    assert abs(initial.det()) == 1

    adapted = initial * child * initial.transpose()
    coupling = adapted[:13,13:]
    tail = adapted[13:,13:]
    H0 = tail - coupling.transpose()*cartan.inverse()*coupling

    scale = ZZ(1)
    for value in H0.list():
        scale = lcm(scale, ZZ(QQ(value).denominator()))
    lll = matrix(ZZ, pari((scale*H0).change_ring(ZZ)).qflllgram())
    assert abs(lll.det()) == 1

    change = block_diagonal_matrix(identity_matrix(ZZ,13), lll.transpose())
    basis = change * initial
    adapted = basis * child * basis.transpose()

    root = adapted[:13,:13]
    coupling = adapted[:13,13:]
    tail = adapted[13:,13:]
    H = tail - coupling.transpose()*root.inverse()*coupling
    assert H.det() == 237
    return basis, adapted, H


def class_order(dual):
    answer = ZZ(1)
    for value in dual:
        answer = lcm(answer, ZZ(QQ(value).denominator()))
    return answer


def d13_correction(dual, root):
    order = class_order(dual)
    expected = {ZZ(1): QQ(0), ZZ(2): QQ(1), ZZ(4): QQ(13)/4}
    assert order in expected
    correction = expected[order]
    raw = QQ(dual * root * dual)
    mod2 = lambda value: QQ(value) - 2*(QQ(value)/2).floor()
    assert mod2(raw) == mod2(correction)
    return order, correction


def minimal_section_for_mw(adapted, z, cap=16384):
    root = adapted[:13,:13]
    coupling = adapted[:13,13:]
    tail = adapted[13:,13:]
    H = tail - coupling.transpose()*root.inverse()*coupling

    z = vector(ZZ,z)
    height = QQ(z*H*z)
    base = vector(ZZ, [0]*13 + list(z))
    pairing = vector(QQ, base * adapted[:, :13])
    dual = pairing * root.inverse()
    order, correction = d13_correction(dual, root)
    target_norm = height + correction
    assert target_norm in ZZ and target_norm >= 4 and target_norm % 2 == 0

    lattice = IntegralLattice(root)
    iterator = lattice.enumerate_close_vectors(-dual)
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

    pole = QQ(height + correction - 4)/2
    assert pole in ZZ and pole >= 0
    pole = ZZ(pole)

    a = ZZ((target_norm - 2)/2)
    section = vector(ZZ, [a,1] + list(chosen))
    U2 = matrix(ZZ, ((0,1),(1,0)))
    child_ns = block_diagonal_matrix(U2, -adapted)
    assert section * child_ns * section == -2
    assert section * child_ns * vector(ZZ,[1,0]+[0]*17) == 1
    assert section * child_ns * vector(ZZ,[-1,1]+[0]*17) == pole

    return {
        "mw": z,
        "height": height,
        "class_order": order,
        "correction": correction,
        "pole": pole,
        "lift": chosen,
        "section": section,
    }


# The equation-level fibre has the SAME E6/E8 component degrees that were
# certified for the physical target, so II*_E8_1 is a genuine equation zero.
O8eq = vector(ZZ, selected["E8"]["simple_root_vectors_in_source_h3_ns"][0])
assert O8eq * ns * O8eq == -2
assert O8eq * ns * F8eq == 1

child_eq, Bzero_eq = child_frame_with_zero(ns, F8eq, O8eq)
assert roots_and_data(child_eq)[2] == (13,312,4)
A13eq, adapted_eq, H13eq = d13_root_adaptation(child_eq)
Badapt_eq = block_diagonal_matrix(identity_matrix(ZZ,2), A13eq) * Bzero_eq
Gadapt_eq = block_diagonal_matrix(matrix(ZZ,((0,1),(1,0))), -adapted_eq)
assert Badapt_eq * ns * Badapt_eq.transpose() == Gadapt_eq

def eq_coords(curve):
    result = vector(QQ,curve) * Badapt_eq.inverse()
    assert all(value in ZZ for value in result)
    return vector(ZZ,result)

# For a degree-d divisor C, [C-d O] has the same positive-frame tail, so this
# tail is exactly its Pic^0 / Abel-Jacobi MW vector.
cS3 = eq_coords(S3std)
assert cS3[1] == 52
zS3 = vector(ZZ, cS3[-4:])
profS3 = minimal_section_for_mw(adapted_eq, zS3)
P8eq = profS3["section"] * Badapt_eq

assert P8eq * ns * P8eq == -2
assert P8eq * ns * F8eq == 1
assert P8eq * ns * O8eq == profS3["pole"]

print(
    "Q8S3EQFRAME_AJ|"
    f"mw={','.join(map(str,zS3))}|height={profS3['height']}|"
    f"class_order={profS3['class_order']}|"
    f"correction={profS3['correction']}|PdotO={profS3['pole']}|"
    "status=PASS_EXACT_EQUATION_FRAME_AJ",
    flush=True,
)

# Cross-check the independently discovered direct Weierstrass pole profile.
DIRECT_MOD = ROOT / "artifacts/local/elkies-k3/q8-s3-direct-x-mod-100003.json"
direct_profile_match = None
if DIRECT_MOD.exists():
    dm = json.loads(DIRECT_MOD.read_text())
    assert dm["status"] == "PASS_DIRECT_ANCHORED_Q8_S3_PROFILE_DISCOVERY"
    direct_pole = ZZ(dm["weierstrass_structure"]["denominator_root_degree"])
    direct_x = (
        ZZ(dm["x"]["numerator_degree"]),
        ZZ(dm["x"]["denominator_degree"]),
    )
    direct_y = (
        ZZ(dm["weierstrass_structure"]["y_abs_numerator_degree"]),
        ZZ(dm["weierstrass_structure"]["y_denominator_degree"]),
    )
    direct_profile_match = (
        direct_pole == profS3["pole"]
        and direct_x == (2*profS3["pole"]+4, 2*profS3["pole"])
        and direct_y == (3*profS3["pole"]+6, 3*profS3["pole"])
    )
    print(
        "Q8S3EQFRAME_DIRECT_COMPARE|"
        f"lattice_PdotO={profS3['pole']}|direct_Z={direct_pole}|"
        f"direct_x={direct_x[0]}/{direct_x[1]}|"
        f"direct_y={direct_y[0]}/{direct_y[1]}|"
        f"match={int(direct_profile_match)}|status=PASS",
        flush=True,
    )
    assert direct_profile_match

# The no-root-correction q24 candidate forced by the equation-frame section.
# For D=O+P-kF, D^2=0 gives k=(P.O-2)/2.
assert (profS3["pole"] - 2) % 2 == 0
q24_twist = ZZ((profS3["pole"] - 2)//2)
D24eq = O8eq + P8eq - q24_twist*F8eq

assert D24eq * ns * D24eq == 0
assert D24eq * ns * F8eq == 2
assert gcd(tuple(D24eq)) == 1

mate24 = isotropic_mate(ns, D24eq)
orth24 = matrix(
    ZZ, [list(D24eq*ns), list(mate24*ns)]
).right_kernel_matrix()
child24 = -(orth24 * ns * orth24.transpose())
q24_root_data = roots_and_data(child24)[2]
q24_is_d12 = q24_root_data == (12,264,4)

print(
    "Q8S3EQFRAME_Q24|"
    f"twist={q24_twist}|square=0|old_q8_degree=2|"
    f"root_data={q24_root_data[0]},{q24_root_data[1]},{q24_root_data[2]}|"
    f"D12={int(q24_is_d12)}|MW_rank_if_rho19={17-q24_root_data[0]}|"
    "status=PASS" if q24_is_d12 else
    "Q8S3EQFRAME_Q24|status=NOT_D12",
    flush=True,
)
assert q24_is_d12

# Optional absolute comparison against the previously transported q24 divisor.
# A mismatch here is informative: it would mean that artifact belongs to the
# historical/physical q8 frame.  It does NOT invalidate the new equation-frame
# D12 neighbour just certified above.
Q24_OLD = ROOT / "artifacts/local/elkies-k3/q8-q24-effective-zero-choices.json"
historical_compare = None
historical_difference = None
if Q24_OLD.exists():
    old24 = json.loads(Q24_OLD.read_text())
    assert old24["status"] == "PASS_EXACT_Q24_EFFECTIVE_ZERO_CHOICES"
    Dhist = vector(ZZ, old24["transport"]["q24_divisor_source_h3_ns"])
    historical_compare = bool(Dhist == D24eq)
    delta24 = D24eq - Dhist
    historical_difference = {
        "vector": list(map(int,delta24)),
        "square": int(delta24*ns*delta24),
        "equation_q8_degree": int(Dhist*ns*F8eq),
    }
    print(
        "Q8S3EQFRAME_HISTORICAL_Q24|"
        f"same_absolute_class={int(historical_compare)}|"
        f"historical_degree_on_F8eq={Dhist*ns*F8eq}|"
        f"difference_square={delta24*ns*delta24}|status=PASS_DIAGNOSTIC",
        flush=True,
    )

# If the characteristic-zero Hensel lift has already completed, bind it into
# the closure certificate.  The direct trace computes AJ(S3), this exact lattice
# calculation identifies AJ(S3), and the nonsingular Hensel lift is the unique
# characteristic-zero lift of that modular section.
EXACT_SECTION = ROOT / "artifacts/local/elkies-k3/q8-s3-direct-section-qq.json"
exact_section_bound = False
if EXACT_SECTION.exists():
    exact = json.loads(EXACT_SECTION.read_text())
    assert exact["status"] == "PASS_EXACT_Q8_S3_DIRECT_SECTION"
    assert exact["verification"]["exact_weierstrass_identity"] is True
    assert exact["profile"]["P_dot_O_from_denominator"] == int(profS3["pole"])
    assert exact["profile"]["x_degrees"] == [
        int(2*profS3["pole"]+4), int(2*profS3["pole"])
    ]
    assert exact["profile"]["y_degrees"] == [
        int(3*profS3["pole"]+6), int(3*profS3["pole"])
    ]
    exact_section_bound = True
    print(
        "Q8S3EQFRAME_CHAR0_BIND|"
        f"PdotO={profS3['pole']}|exact_identity=1|"
        "status=PASS_EXACT_CHAR0_AJ_BINDING",
        flush=True,
    )

closure_payload = {
    "equation_zero": {
        "curve": "II*_E8_1",
        "source_h3_ns": list(map(int,O8eq)),
    },
    "AJ_S3_equation_frame": {
        "degree_of_S3": int(cS3[1]),
        "mw_coordinates_in_deterministic_equation_D13_basis":
            list(map(int,zS3)),
        "height": str(profS3["height"]),
        "D13_class_order": int(profS3["class_order"]),
        "D13_local_correction": str(profS3["correction"]),
        "P_dot_O": int(profS3["pole"]),
        "effective_section_source_h3_ns": list(map(int,P8eq)),
    },
    "direct_modular_profile_match": direct_profile_match,
    "equation_frame_q24": {
        "divisor_source_h3_ns": list(map(int,D24eq)),
        "decomposition": f"D=O+P-{q24_twist}F",
        "old_q8_degree": int(D24eq*ns*F8eq),
        "square": int(D24eq*ns*D24eq),
        "primitive": bool(gcd(tuple(D24eq)) == 1),
        "child_root_data": list(map(int,q24_root_data)),
        "child_root_lattice": "D12",
        "MW_rank_if_rho19": int(17-q24_root_data[0]),
    },
    "historical_q24_absolute_same": historical_compare,
    "historical_q24_difference": historical_difference,
    "exact_characteristic_zero_section_bound": exact_section_bound,
}

print(
    "Q8S3EQFRAME_RESULT|"
    f"height={profS3['height']}|correction={profS3['correction']}|"
    f"PdotO={profS3['pole']}|q24_twist={q24_twist}|"
    "child=D12|MW=5|"
    f"char0_bound={int(exact_section_bound)}|"
    "status=PASS_EXACT_AJ_Q24_LATTICE_CLOSURE",
    flush=True,
)

physical = vector(ZZ, selected["source_h3_ns_vector"])
assert physical * ns * physical == 0
delta = F8eq - physical
print(
    "Q8EQNS_COMPARE|"
    f"physical_S3old_degree={int(physical*ns*basis_sections[2])}|"
    f"equation_S3std_degree={degree_S3}|"
    f"difference_square={int(delta*ns*delta)}|"
    f"difference_q6_degree={int(delta*ns*F6)}|status=PASS",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-q8-equation-ns-divisor.v1",
    "status": "PASS_EXACT_Q8_EQUATION_NS_DIVISOR",
    "q6": {
        "standard_zero_old_group_MW": list(map(int,mw_Ostd)),
        "marked_section_old_group_MW": list(map(int,mw_marked)),
        "S3_standard_old_group_MW": list(map(int,mw_S3std)),
        "standard_intersections": {
            "Ostd_marked": int(Ostd*ns*Smarked),
            "Ostd_S3": int(Ostd*ns*S3std),
        },
    },
    "q8_equation_fibre": {
        "source_h3_ns_vector": list(map(int,F8eq)),
        "square": int(F8eq*ns*F8eq),
        "q6_degree": int(F8eq*ns*F6),
        "S3_degree": degree_S3,
        "root_data": list(map(int,root_data)),
    },
    "physical_target_comparison": {
        "source_h3_ns_vector": list(map(int,physical)),
        "physical_old_q6_third_degree": int(physical*ns*basis_sections[2]),
        "difference_vector": list(map(int,delta)),
        "difference_square": int(delta*ns*delta),
        "difference_q6_degree": int(delta*ns*F6),
    },
    "AJ_q24_equation_frame_closure": closure_payload,
    "boundary": (
        "This pins the equation-level q8 NS fibre class in the source H3 frame. "
        "It does not yet identify the standard zero / MW basis of the D13 "
        "Jacobian produced by the binary quartic."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q8EQNS_RESULT|S3_degree=52|root=D13|q24_child=D12|MW=5|"
    "status=PASS_EXACT_Q8_EQUATION_NS_AND_AJ_Q24_CLOSURE",
    flush=True,
)
