#!/usr/bin/env sage -python
"""
Analyze AJ_q8(S3) relative to the THREE genuine low-degree q8 sections.

This fixes two mistakes in the previous lattice diagnostic:
  1. it does not call an arbitrary isotropic mate "the zero section";
  2. it reconstructs the actual certified third q6 section class using the
     exact Weyl/root-shift recipe from certify_h3_q6_weyl_section_transport.sage.

The component-nef physical q8 target has exactly three visible old-q6 fibre
components of q8 degree one:
    E6 simple root 1,
    E6 simple root 5,
    E8 simple root 1.

The first two are the two IV* branches already sampled equation-level.  The
third is an II*/E8 component and has not yet been used as a q8 origin.

For each genuine section R, form the U-plane <F8,R+F8>, compute its D13
root/MW quotient, and read the Abel-Jacobi class of the degree-d multisection
S3 relative to R.  This gives exact:
    MW height,
    D13 discriminant coset/local correction,
    P.O,
    predicted x/y coordinate degrees.

Run:
  sage -python ~/Downloads/analyze_h92_q8_visible_zero_profiles.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    QQ, ZZ, block_diagonal_matrix, identity_matrix, lcm,
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
TARGET = ROOT / "artifacts/local/elkies-k3/q8-target-component-nef.json"

for path in (FRAME, TARGET):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

target = json.loads(TARGET.read_text())
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert target["normalization"]["representative"] == "component-nef"

U = matrix(ZZ, ((0, 1), (1, 0)))
REFLECTIONS = (
    1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3,
    1, 7, 6, 5, 4, 2, 3, 4, 5, 6, 7,
)
H3_LIFTS = matrix(ZZ, [
    [-5, -4, -3, 0, 0, 0, 0, 0, 0, 0, 0, -4, 1, 0, -4, 2, -2],
    [-10, -8, -6, 0, 0, 0, 0, 0, 0, 0, 0, -8, 4, 1, -8, 5, -4],
    [-5, -4, -3, 0, 0, 0, 0, 0, 0, 0, -3, 2, 0, -4, 2, -2],
])
OLD_ZERO_ROOT_SHIFTS = matrix(ZZ, [
    [5, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, -1, 4],
    [12, 10, 8, 0, 0, 0, 0, 0, 0, 0, 0, 6, -1, 9],
    [5, 4, 3, 0, 0, 0, 0, 0, 0, 0, 2, 0, 4],
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
    assert mate * ns * mate == 0
    assert mate * ns * fiber == 1
    return mate


def roots_and_data(gram):
    qf = pari(gram).qfminim(2)
    count = ZZ(qf[0])
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


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    rank = data[0]
    regular = None
    for shift in range(1, 1000):
        candidate = vector(ZZ, [
            (i + 1)**2 + shift*(i + 1) + 1
            for i in range(gram.nrows())
        ])
        if all(candidate * r != 0 for r in roots):
            regular = candidate
            break
    assert regular is not None
    positive = [r for r in roots if regular * r > 0]
    positive_set = {tuple(r) for r in positive}
    simple = []
    for r in positive:
        if not any(tuple(r-left) in positive_set for left in positive):
            simple.append(r)
    simple = matrix(ZZ, [list(r) for r in simple])
    assert simple.nrows() == simple.rank() == rank
    return simple, simple * gram * simple.transpose()


def d13_adaptation(child):
    unused_roots, root_basis, invariants = roots_and_data(child)
    assert invariants == (13, 312, 4)
    simple, cartan = deterministic_simple_roots(child)
    assert cartan.det() == 4

    smith, unused_left, smith_right = root_basis.smith_form()
    assert tuple(abs(smith[i, i]) for i in range(13)) == (1,) * 13
    completion = smith_right.inverse()
    adapted_basis = simple.stack(completion[13:])
    assert abs(adapted_basis.det()) == 1

    adapted = adapted_basis * child * adapted_basis.transpose()
    A = adapted[:13, :13]
    B = adapted[:13, 13:]
    C = adapted[13:, 13:]
    H = C - B.transpose()*A.inverse()*B

    scale = lcm(v.denominator() for v in H.list())
    lll = matrix(ZZ, pari((scale*H).change_ring(ZZ)).qflllgram())
    assert abs(lll.det()) == 1
    change = block_diagonal_matrix(
        identity_matrix(ZZ, 13), lll.transpose()
    )
    adapted_basis = change * adapted_basis
    adapted = adapted_basis * child * adapted_basis.transpose()
    A = adapted[:13, :13]
    B = adapted[:13, 13:]
    C = adapted[13:, 13:]
    H = C - B.transpose()*A.inverse()*B
    return adapted_basis, adapted, H


def correction_from_coset(adapted, mw):
    """D13 has discriminant group Z/4: corrections 0,1,13/4."""
    A = adapted[:13, :13]
    B = adapted[:13, 13:]
    q = A.inverse() * B * vector(QQ, mw)

    # Reduce each coordinate modulo Z only to classify the discriminant order.
    fracs = []
    for value in q:
        floor = value.floor()
        fracs.append(QQ(value - floor))

    if all(v == 0 for v in fracs):
        return QQ(0), 1
    if all((2*v).denominator() == 1 for v in fracs):
        return QQ(1), 2
    if all((4*v).denominator() == 1 for v in fracs):
        return QQ(13)/4, 4
    raise ArithmeticError(f"unexpected D13 discriminant coset: {fracs}")


# ---------------------------------------------------------------------------
# Source H3 NS and the ACTUAL certified third q6 section class.
# ---------------------------------------------------------------------------
frame = load_gram(FRAME)
ns = block_diagonal_matrix(U, -frame)
old_fiber = vector(ZZ, [1, 0] + [0]*17)
old_zero = vector(ZZ, [-1, 1] + [0]*17)
simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)

raw_fiber = vector(
    ZZ,
    [3, 2] + [
        0, 0, -1, -1, -1, -1, -1,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
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


nef_q6_fiber = weyl_transport(raw_fiber)
assert nef_q6_fiber * ns * nef_q6_fiber == 0
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

raw_zero = inverse_weyl_transport(old_zero)
raw_zero_coords = vector(ZZ, raw_zero * raw_transport.inverse())
raw_zero_mw_lift = vector(ZZ, raw_zero_coords[2:])

sections = []
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
    assert candidate * ns * nef_q6_fiber == 1
    sections.append(candidate)

S3 = sections[2]
assert S3 * ns * S3 == -2

# Regression against the known q6 data.
assert S3 * ns * old_fiber == 44
assert S3 * ns * nef_q6_fiber == 1
print(
    "Q8VISIBLE_S3|source_class=certified_weyl_transport|"
    "old_fiber_degree=44|q6_degree=1|status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# Physical component-nef q8 fibre and its three visible degree-one sections.
# ---------------------------------------------------------------------------
F8 = vector(ZZ, target["selected_q8"]["source_h3_ns_vector"])
assert F8 * ns * F8 == 0
degree8 = int(S3 * ns * F8)
print(
    f"Q8VISIBLE_FIBRE|S3_degree={degree8}|"
    "equation_bridge_degree=52|"
    f"match={int(degree8 == 52)}",
    flush=True,
)

visible = []
for family in ("E6", "E8"):
    roots = target["selected_q8"][family][
        "simple_root_vectors_in_source_h3_ns"
    ]
    degrees = target["selected_q8"][family]["component_degrees"]
    for index, (root_data, degree) in enumerate(zip(roots, degrees), 1):
        if int(degree) != 1:
            continue
        R = vector(ZZ, root_data)
        assert R * ns * R == -2
        assert F8 * ns * R == 1
        visible.append((family, index, R))

assert [(a,b) for a,b,unused in visible] == [
    ("E6", 1), ("E6", 5), ("E8", 1)
]

# If this is not 52, the physical old-zero target is NOT yet the
# equation-level translated q8 fibre.  Still report visible-zero data, but do
# not pretend it predicts the equation sampler.
if degree8 != 52:
    print(
        "Q8VISIBLE_BOUNDARY|physical_target_not_equation_fibre|"
        "reason=standard_zero_translation_unresolved|status=STOP",
        flush=True,
    )
    raise SystemExit(0)

# ---------------------------------------------------------------------------
# For each genuine section R, build the q8 D13 frame with R as zero and
# compute AJ_R(S3).
# ---------------------------------------------------------------------------
for family, index, R in visible:
    mate = R + F8
    assert mate * ns * mate == 0
    assert mate * ns * F8 == 1

    orth = matrix(
        ZZ, [list(F8*ns), list(mate*ns)]
    ).right_kernel_matrix()
    child = -(orth * ns * orth.transpose())
    assert roots_and_data(child)[2] == (13, 312, 4)

    frame8 = matrix(
        ZZ, [list(F8), list(mate)] + [list(r) for r in orth.rows()]
    )
    assert abs(frame8.det()) == 1
    coords = vector(ZZ, S3 * frame8.inverse())
    assert int(coords[1]) == 52

    # Tail is the degree-zero class modulo U and hence determines AJ_R.
    tail = vector(ZZ, coords[2:])
    adapted_basis, adapted, height = d13_adaptation(child)
    root_mw_coords = vector(ZZ, tail * adapted_basis.inverse())
    mw = vector(ZZ, root_mw_coords[13:])
    hval = QQ(mw * height * mw)

    correction, coset_order = correction_from_coset(adapted, mw)
    PO = (hval - 4 + correction)/2
    assert PO in ZZ and PO >= 0
    PO = ZZ(PO)

    xden = 2*PO
    yden = 3*PO
    xnum = xden + 4
    ynum = yden + 6

    print(
        f"Q8VISIBLE_PROFILE|zero={family}_{index}|"
        f"mw={','.join(map(str,mw))}|"
        f"height={hval}|D13_correction={correction}|"
        f"coset_order={coset_order}|O={PO}|"
        f"x={xnum}/{xden}|y={ynum}/{yden}|status=PASS",
        flush=True,
    )

print(
    "Q8VISIBLE_RESULT|status=PASS_VISIBLE_Q8_ZERO_PROFILES",
    flush=True,
)
