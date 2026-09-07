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
    QQ, ZZ, block_diagonal_matrix, gcd, matrix, pari, vector, xgcd
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
        else Path("artifacts/local/elkies-k3/q8-equation-ns-divisor.json")
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
    matrix.identity(QQ, 19)
    - matrix(QQ, ns)
      * matrix(QQ, roots_source.transpose())
      * matrix(QQ, root_gram.inverse())
      * matrix(QQ, roots_source)
)


def shioda(section):
    horizontal = section - Oold - (section*ns*Oold + 2)*F6
    assert horizontal*ns*F6 == 0 and horizontal*ns*Oold == 0
    phi = vector(QQ, horizontal) * projection
    phi -= (phi*ns*Oold) * vector(QQ, F6)
    assert phi*ns*F6 == 0 and phi*ns*Oold == 0
    return phi


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

# Exact equation-level q8 fibre.
#
# The component-nef target already stores the exact isotropic source-H3 NS
# class.  Do not interpret its local `vertical_cycle` valuation data as
# coefficients that can simply be added to Ostd+Smarked.  Instead derive the
# actual vertical difference in the q6 fibre lattice.
selected = target["selected_q8"]

# The stored component-nef target is in the RAW q6 representative:
# raw_fiber, inverse_weyl_transport(Oold).  Our reconstructed standard-zero
# sections are in the transported representative F6,Oold.  Apply the same
# Weyl isometry to the entire q8 target before comparing divisors.
physical_raw = vector(ZZ, selected["source_h3_ns_vector"])
assert physical_raw * ns * physical_raw == 0

physical = weyl_transport(physical_raw)
assert inverse_weyl_transport(physical) == physical_raw
assert physical * ns * physical == 0
assert physical * ns * F6 == 2

# The stored component-nef target is the PHYSICAL divisor before
# fibrewise translation by -old_zero.  Its horizontal pair is
#     Oold + Qphysical,
# where Qphysical has old-zero MW coordinate (-2,-2,0).
#
# Translation by Ostd=(-2,1,0) in the old group sends
#     Oold      -> Ostd
#     Qphysical -> Smarked=(-4,-1,0).
mw_Qphysical = vector(ZZ, (-2,-2,0))
Qphysical, Qphysical_oldO, hQphysical = section_from_old_mw(mw_Qphysical)
assert Qphysical_oldO == 10
assert hQphysical == 24
assert mw_Qphysical + mw_Ostd == mw_marked

physical_horizontal = Oold + Qphysical
equation_horizontal = Ostd + Smarked
assert physical_horizontal * ns * F6 == 2
assert equation_horizontal * ns * F6 == 2

E6roots = [
    weyl_transport(vector(ZZ, row))
    for row in selected["E6"]["simple_root_vectors_in_source_h3_ns"]
]
E8roots = [
    weyl_transport(vector(ZZ, row))
    for row in selected["E8"]["simple_root_vectors_in_source_h3_ns"]
]

vertical_basis = matrix(
    QQ,
    [list(F6)] + [list(r) for r in E6roots] + [list(r) for r in E8roots],
)
vertical_delta = physical - physical_horizontal

derived = vertical_basis.transpose().solve_right(vector(QQ, vertical_delta))
assert len(derived) == 15
assert all(value in ZZ for value in derived)
derived = vector(ZZ, derived)

reconstructed_delta = vector(QQ, derived) * vertical_basis
assert reconstructed_delta == vector(QQ, vertical_delta)

derived_fibre = int(derived[0])
derived_E6 = list(map(int, derived[1:7]))
derived_E8 = list(map(int, derived[7:15]))

stored_fibre = int(selected["vertical_fibre_coefficient"])
stored_E6 = list(map(int, selected["E6"]["vertical_cycle"]))
stored_E8 = list(map(int, selected["E8"]["vertical_cycle"]))

print(
    "Q8EQNS_VERTICAL|"
    f"derived_fibre={derived_fibre}|"
    f"derived_E6={','.join(map(str,derived_E6))}|"
    f"derived_E8={','.join(map(str,derived_E8))}|"
    f"stored_fibre={stored_fibre}|"
    f"stored_E6={','.join(map(str,stored_E6))}|"
    f"stored_E8={','.join(map(str,stored_E8))}|"
    f"same={int(derived_fibre==stored_fibre and derived_E6==stored_E6 and derived_E8==stored_E8)}|"
    "status=PASS_EXACT_VERTICAL_DECOMPOSITION",
    flush=True,
)

# Translation is by an identity-component q6 section:
# II* has trivial component group and mw_Ostd has trivial IV* class.
# Hence the fibre/root component classes are fixed; only the horizontal
# sections move.
assert iv_correction(mw_Ostd) == 0
F8eq = equation_horizontal + vector(ZZ, vertical_delta)
assert F8eq * ns * F8eq == 0
assert F8eq * ns * F6 == 2
assert gcd(tuple(F8eq)) == 1

# The stored component degrees must replay exactly.
for family, roots in (("E6", E6roots), ("E8", E8roots)):
    expected = list(map(int, selected[family]["component_degrees"]))
    actual = [int(F8eq * ns * root) for root in roots]
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

# Physical target round-trips to the stored raw component-nef class.
assert inverse_weyl_transport(physical) == physical_raw

# This difference is precisely the full-NS effect of translating the
# horizontal degree-two divisor from the physical old-zero presentation
# to the standard-Weierstrass presentation.
delta = F8eq - physical
assert delta * ns * F6 == 0
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
    "Q8EQNS_RESULT|S3_degree=52|root=D13|"
    "status=PASS_EXACT_Q8_EQUATION_NS_DIVISOR",
    flush=True,
)
