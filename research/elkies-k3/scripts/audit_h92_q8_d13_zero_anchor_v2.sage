#!/usr/bin/env sage -python
"""
Pin the abstract D13 zero to the actual corrected q8 quartic/Jacobian.

This script resolves the remaining origin ambiguity without assuming that
independently constructed q8 chamber representatives are already identical.

It reconstructs the classifier's nef->dominant bridge literally from the
current source constants and matrices, transports the pinned D13 zero through
that bridge and through the 102 physical E6/E8 reflections used by the
equation-level q8 compiler, identifies the resulting q8 zero among the actual
q6 fibre components, and evaluates the exact q8 quartic at that component's
fixed old-base value T=T0.

If that quartic point is rational, the script constructs the elliptic curve
birationally anchored at that point (so the pinned D13 zero maps to infinity)
and checks that it is isomorphic over QQ(U) to the certified I9*+9I1 child.

No binary-quartic covariant image is used to identify the origin.

Run:
  sage -python ~/Downloads/audit_h92_q8_d13_zero_anchor.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    EllipticCurve, PolynomialRing, QQ, ZZ, block_diagonal_matrix,
    identity_matrix, matrix, pari, sage_eval, vector, xgcd
)

# Current classifier constant, copied verbatim from
# classify_h3_q6_child_q8_orbits.sage and checked against the generated artifact.
NEF_Q8_WITNESS_Q6_RAW = vector(ZZ, (
    -5, -4, -3, 4, 5, 7, 10, 8, 6, 4, 2, -4, 2, -2, -4, 0, -2,
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
    row = vector(ZZ, row)
    root = vector(ZZ, root)
    assert root * gram * root == -2
    return row + (row * gram * root) * root


def highest_root(cartan):
    half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
    roots = [vector(ZZ, r) for r in half]
    roots += [-vector(ZZ, r) for r in half]
    positive = [r for r in roots if all(v >= 0 for v in r)]
    assert positive
    return max(positive, key=lambda r: sum(r))


def proportional_to_fibre(diff, fibre):
    """Return integer k if diff=k*fibre, else None."""
    k = None
    for d, f in zip(diff, fibre):
        if f:
            if d % f:
                return None
            here = d // f
            if k is None:
                k = here
            elif k != here:
                return None
        elif d:
            return None
    return ZZ(0) if k is None else ZZ(k)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--target", type=Path)
parser.add_argument("--q8-child", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts" / "generated-results"
LOCAL = ROOT / "artifacts" / "local" / "elkies-k3"

FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
PINNED_D13 = (
    ROOT / "elkies-k3/data/fibrations/"
           "h3_q6_q8_d13_mw4_root_adapted_frame.txt"
)
ORBITS = GEN / "elkies-k3-h3-q6-q8-orbits.json"
Q6_CHILD = GEN / "elkies-k3-h92-q6-child-jacobian.json"
TARGET = (
    args.target.resolve()
    if args.target
    else LOCAL / "q8-target-component-nef.json"
)

if args.q8_child:
    Q8_CHILD = args.q8_child.resolve()
else:
    q8_candidates = [
        GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
        LOCAL / "q8-corrected2cover-qq-child.json",
    ]
    Q8_CHILD = next((p for p in q8_candidates if p.exists()), q8_candidates[0])

OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q8-d13-zero-anchor.json"
)

for path in (FRAME, PINNED_D13, ORBITS, Q6_CHILD, TARGET, Q8_CHILD):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

source_frame = load_gram(FRAME)
pinned_d13 = load_gram(PINNED_D13)
q8data = json.loads(ORBITS.read_text())
q6child = json.loads(Q6_CHILD.read_text())
target = json.loads(TARGET.read_text())
q8child = json.loads(Q8_CHILD.read_text())

assert q8data["status"] == "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION"
assert q6child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert target["normalization"]["representative"] == "component-nef"
assert q8child["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"

U2 = matrix(ZZ, ((0, 1), (1, 0)))
source_ns = block_diagonal_matrix(U2, -source_frame)
source_O = vector(ZZ, [-1, 1] + [0] * 17)

# ===========================================================================
# 1. Certified q6 coordinate tower.
# ===========================================================================

B6 = matrix(ZZ, q8data["q6"]["neighbor_basis_in_source_ns"])
assert abs(B6.det()) == 1
q6_raw_ns = B6 * source_ns * B6.transpose()
assert q6_raw_ns[:2, :2] == U2
q6_raw = -q6_raw_ns[2:, 2:]
assert q6_raw.det() == 948

root_mw = matrix(ZZ, q8data["q6"]["root_mw_basis_in_child"])
assert abs(root_mw.det()) == 1
root_mw_frame = root_mw * q6_raw * root_mw.transpose()
assert root_mw_frame == matrix(ZZ, q8data["q6"]["root_adapted_gram"])

simple14 = matrix(ZZ, q8data["q8"]["simple_root_change_in_root_block"])
simple_change = block_diagonal_matrix(simple14, identity_matrix(ZZ, 3))
simple_to_raw = simple_change * root_mw
simple_frame = matrix(ZZ, q8data["q8"]["simple_frame_gram"])
assert simple_to_raw * q6_raw * simple_to_raw.transpose() == simple_frame

Bsimple = (
    block_diagonal_matrix(identity_matrix(ZZ, 2), simple_to_raw) * B6
)
q6_simple_ns = block_diagonal_matrix(U2, -simple_frame)
assert Bsimple * source_ns * Bsimple.transpose() == q6_simple_ns

print(
    "Q8ZERO_COORD|q6_raw=PASS|root_mw=PASS|simple=PASS|source=PASS",
    flush=True,
)

# ===========================================================================
# 2. Rebuild nef -> dominant exactly as the classifier source does.
# ===========================================================================

nef_record = q8data["q8"]["nef_representative"]
nef_witness_raw = vector(ZZ, nef_record["witness_q6_child"])
assert nef_witness_raw * q6_raw * nef_witness_raw == 16
nef_root_mw = nef_witness_raw * root_mw.inverse()
nef_simple = vector(ZZ, nef_root_mw * simple_change.inverse())
assert tuple(nef_simple[-3:]) == tuple(nef_record["mw_projection"])

# Recompute the bridge literally from the SAME artifact matrices/witness, then
# compare it to the artifact's stored bridge.  This prevents mixing a generated
# certificate with a copied constant from a different classifier revision.
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

computed_bridge = tuple((i + 1, int(p)) for i, p in bridge)
stored_bridge = tuple(
    (int(i), int(p)) for i, p in nef_record["to_dominant_reflections"]
)
assert computed_bridge == stored_bridge

dominant_matches = [
    hit for hit in q8data["q8"]["d13_mw4_hits"]
    if tuple(hit["witness_simple_frame"]) == tuple(dominant_simple)
]
assert len(dominant_matches) == 1
dominant_hit = dominant_matches[0]
assert dominant_hit["mw_projection"] == nef_record["mw_projection"]
assert matrix(ZZ, dominant_hit["d13_root_adapted_gram"]) == pinned_d13

# Verify the same four reflections on full NS classes and map the roots to source.
simple_root_classes = tuple(
    vector(ZZ, [0, 0] + [ZZ(j == i) for j in range(17)])
    for i in range(14)
)
simple_root_source = tuple(root * Bsimple for root in simple_root_classes)
for r in simple_root_source:
    assert r * source_ns * r == -2

F_nef_simple = vector(ZZ, [4, 2] + list(nef_simple))
F_dom_simple = vector(ZZ, [4, 2] + list(dominant_simple))

# Direct raw-q6 -> source identity is the artifact's native construction.
F_nef_source_direct = vector(ZZ, [4, 2] + list(nef_witness_raw)) * B6
assert F_nef_source_direct == vector(ZZ, nef_record["fiber_source_h3_ns"])

# Independent simple-frame route must give the same class.
F_nef_source = F_nef_simple * Bsimple
F_dom_source = F_dom_simple * Bsimple
assert F_nef_source == F_nef_source_direct
assert F_dom_source == vector(ZZ, dominant_hit["fiber_source_h3_ns"])

check = vector(ZZ, F_nef_source)
for i, pairing in bridge:
    # Positive-frame label is pairing; NS intersection is -pairing.
    assert check * source_ns * simple_root_source[i] == -pairing
    check = reflect(check, source_ns, simple_root_source[i])
assert check == F_dom_source

# ===========================================================================
# 3. Pinned D13 zero -> dominant -> nef source coordinates.
# ===========================================================================

B8dom = matrix(ZZ, dominant_hit["neighbor_basis_in_q6_ns"])
raw_d13 = matrix(ZZ, dominant_hit["child_frame"])
assert B8dom * q6_simple_ns * B8dom.transpose() == block_diagonal_matrix(U2, -raw_d13)

A13 = matrix(ZZ, dominant_hit["d13_root_adapted_basis_in_child"])
assert A13 * raw_d13 * A13.transpose() == pinned_d13
Bpinned_to_simple = (
    block_diagonal_matrix(identity_matrix(ZZ, 2), A13) * B8dom
)
Gd13 = block_diagonal_matrix(U2, -pinned_d13)
assert Bpinned_to_simple * q6_simple_ns * Bpinned_to_simple.transpose() == Gd13

O8_pinned = vector(ZZ, [-1, 1] + [0] * 17)
O_dom_simple = O8_pinned * Bpinned_to_simple
O_dom_source = O_dom_simple * Bsimple

assert O_dom_source * source_ns * O_dom_source == -2
assert O_dom_source * source_ns * F_dom_source == 1

# Invert the four bridge reflections on the zero.
O_nef_source = vector(ZZ, O_dom_source)
F_back = vector(ZZ, F_dom_source)
for i, unused_pairing in reversed(bridge):
    root = simple_root_source[i]
    O_nef_source = reflect(O_nef_source, source_ns, root)
    F_back = reflect(F_back, source_ns, root)

assert F_back == F_nef_source
assert O_nef_source * source_ns * O_nef_source == -2
assert O_nef_source * source_ns * F_nef_source == 1

print(
    "Q8ZERO_BRIDGE|nef_to_dominant={}|forward=PASS|inverse_zero=PASS".format(len(bridge)),
    flush=True,
)

# ===========================================================================
# 4. Reconstruct the physical E6+E8 roots and replay the target's reduction.
# ===========================================================================

q6_F_source = vector(ZZ, [3, 2] + [
    0, 0, -1, -1, -1, -1, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
])

old_simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(i == node) for i in range(17)])
    for node in range(15)
)

q6_O_source = vector(ZZ, source_O)
for node1 in reversed(Q6_REFLECTIONS):
    q6_O_source = reflect(q6_O_source, source_ns, old_simple[node1 - 1])

assert q6_F_source * source_ns * q6_F_source == 0
assert q6_O_source * source_ns * q6_O_source == -2
assert q6_O_source * source_ns * q6_F_source == 1

q6_orth = matrix(
    ZZ,
    [
        list(q6_F_source * source_ns),
        list((q6_O_source + q6_F_source) * source_ns),
    ],
).right_kernel_matrix()
physical_child = -(q6_orth * source_ns * q6_orth.transpose())
qf_basis = matrix(
    ZZ, pari(physical_child).qfminim(2)[2]
).transpose().row_module().basis_matrix()
assert qf_basis.rank() == 14

e6_qf = matrix(ZZ, [list(qf_basis[i]) for i in E6_QF_INDICES])
e6_roots = matrix(
    ZZ, [list(vector(ZZ, row) * e6_qf) for row in E6_SIMPLE_IN_QF]
) * q6_orth
physical_roots = qf_basis * q6_orth
e8_roots = physical_roots[4:12, :]

assert -e6_roots * source_ns * e6_roots.transpose() == E6_CARTAN
assert -e8_roots * source_ns * e8_roots.transpose() == E8_CARTAN
assert e6_roots * source_ns * e8_roots.transpose() == matrix(ZZ, 6, 8)

actual_roots = tuple(e6_roots.rows()) + tuple(e8_roots.rows())

F_actual = vector(ZZ, F_nef_source)
O_actual = vector(ZZ, O_nef_source)
physical_reflections = []

for unused in range(500):
    pairings = [int(F_actual * source_ns * r) for r in actual_roots]
    negative = [i for i, value in enumerate(pairings) if value < 0]
    if not negative:
        break
    i = negative[0]
    p = pairings[i]
    r = actual_roots[i]
    F_actual = reflect(F_actual, source_ns, r)
    O_actual = reflect(O_actual, source_ns, r)
    physical_reflections.append((i, p))
else:
    raise RuntimeError("physical q8 reduction did not terminate")

assert len(physical_reflections) == 102
assert O_actual * source_ns * O_actual == -2
assert O_actual * source_ns * F_actual == 1

target_F = vector(ZZ, target["selected_q8"]["source_h3_ns_vector"])
endpoint_match = (F_actual == target_F)

print(
    "Q8ZERO_PHYSICAL|reflections={}|target_endpoint_match={}|"
    "O_square={}|O_degree={}".format(
        len(physical_reflections), int(endpoint_match),
        O_actual * source_ns * O_actual,
        O_actual * source_ns * F_actual,
    ),
    flush=True,
)

if not endpoint_match:
    diff = F_actual - target_F
    print(
        "Q8ZERO_PHYSICAL_MISMATCH|"
        f"diff_square={diff * source_ns * diff}|"
        f"diff_q6_degree={diff * source_ns * q6_F_source}|"
        f"direct={','.join(map(str, F_actual))}|"
        f"target={','.join(map(str, target_F))}",
        flush=True,
    )
    raise RuntimeError(
        "reconstructed classifier/physical path does not land on the q8 target artifact"
    )

# ===========================================================================
# 5. Identify the actual q8 zero among q6 fibre components.
# ===========================================================================

e6_highest = highest_root(E6_CARTAN)
e8_highest = highest_root(E8_CARTAN)
e6_affine = q6_F_source - e6_highest * e6_roots
e8_affine = q6_F_source - e8_highest * e8_roots

component_candidates = []
for i, row in enumerate(e6_roots.rows(), 1):
    component_candidates.append((f"IV*_E6_{i}", "IV*", vector(ZZ, row)))
component_candidates.append(("IV*_affine", "IV*", vector(ZZ, e6_affine)))
for i, row in enumerate(e8_roots.rows(), 1):
    component_candidates.append((f"II*_E8_{i}", "II*", vector(ZZ, row)))
component_candidates.append(("II*_affine", "II*", vector(ZZ, e8_affine)))

exact_matches = [
    (name, kind, comp)
    for name, kind, comp in component_candidates
    if O_actual == comp
]
fibre_shift_matches = []
for name, kind, comp in component_candidates:
    k = proportional_to_fibre(O_actual - comp, q6_F_source)
    if k is not None:
        fibre_shift_matches.append((name, kind, int(k), comp))

print(
    "Q8ZERO_COMPONENT|exact_matches={}|shift_matches={}".format(
        ",".join(name for name, _, _ in exact_matches) or "none",
        ",".join(f"{name}:{k}" for name, _, k, _ in fibre_shift_matches) or "none",
    ),
    flush=True,
)

if len(exact_matches) == 1:
    component_name, component_kind, unused = exact_matches[0]
    fibre_shift = 0
elif len(fibre_shift_matches) == 1:
    component_name, component_kind, fibre_shift, unused = fibre_shift_matches[0]
else:
    raise RuntimeError(
        "could not uniquely identify the transported D13 zero as a q6 fibre component"
    )

# A genuine irreducible zero section should be the component itself, not C+kF.
assert fibre_shift == 0

# ===========================================================================
# 6. Evaluate the exact q8 quartic at that fixed q6-base value.
# ===========================================================================

old_R = PolynomialRing(QQ, "T")
old_T = old_R.gen()
q6_factor_text = next(
    item["factor"] for item in q6child["finite_fibres"]
    if item["kodaira"] == component_kind
)
q6_factor = old_R(sage_eval(q6_factor_text, locals={"T": old_T}))
assert q6_factor.degree() == 1
T0 = QQ(-q6_factor[0] / q6_factor[1])

UR = PolynomialRing(QQ, "U")
U = UR.gen()
UF = UR.fraction_field()
TR = PolynomialRing(UF, "T")
T = TR.gen()

quartic_text = q8child["pencil"]["branch_quartic"]
quartic = TR(sage_eval(quartic_text, locals={"U": U, "T": T}))
assert quartic.degree() == 4

q_at_zero = UF(quartic(UF(T0)))
quartic_point_rational = bool(q_at_zero.is_square())

if quartic_point_rational:
    W0 = UF(0) if not q_at_zero else q_at_zero.sqrt()
else:
    W0 = None

print(
    "Q8ZERO_QUARTIC|component={}|kind={}|T0={}|value_zero={}|"
    "is_square={}|sqrt={}".format(
        component_name, component_kind, T0,
        int(not q_at_zero), int(quartic_point_rational),
        "NA" if W0 is None else str(W0),
    ),
    flush=True,
)

if not quartic_point_rational:
    raise RuntimeError("transported zero does not give a rational point on the q8 quartic")

# ===========================================================================
# 7. Build a Weierstrass model anchored at this quartic point.
# ===========================================================================

rR = PolynomialRing(UF, "r")
r = rR.gen()
shifted = rR(quartic(r + UF(T0)))
assert shifted.degree() == 4

e = UF(shifted[0])
d = UF(shifted[1])
c = UF(shifted[2])
b = UF(shifted[3])
a = UF(shifted[4])
assert e == q_at_zero

anchor_type = None
anchor_curves = []

if not e:
    # Branch-point case.  With d != 0:
    #   X=d/r, Y=d*w/r^2
    # gives Y^2=X^3+cX^2+bd X+a d^2 and sends (r,w)=(0,0) to infinity.
    assert d
    Eanchor = EllipticCurve(UF, [0, c, 0, b*d, a*d**2])
    anchor_type = "branch_point_q0"
    anchor_curves.append(("zero", UF(0), Eanchor))
else:
    # Nonbranch rational point.  For each sign q with q^2=e, Washington
    # Theorem 2.17 gives an anchored generalized Weierstrass model sending
    # (r,w)=(0,q) to infinity.  The correct sign is a separate local-chart
    # question if both are possible.
    q = W0
    for sign_name, qsign in (("plus", q), ("minus", -q)):
        a1 = d/qsign
        a2 = c - (d/(2*qsign))**2
        a3 = 2*qsign*b
        a4 = -4*qsign**2*a
        a6 = a2*a4
        Eanchor = EllipticCurve(UF, [a1, a2, a3, a4, a6])
        anchor_curves.append((sign_name, qsign, Eanchor))
    anchor_type = "nonbranch_two_signs"

# Canonical certified D13 short model.
child_data = q8child["child"]
Amin = UR([QQ(v) for v in child_data["minimal_A_coefficients_low_to_high"]])
Bmin = UR([QQ(v) for v in child_data["minimal_B_coefficients_low_to_high"]])
Ecanon = EllipticCurve(UF, [0, 0, 0, UF(Amin), UF(Bmin)])

anchor_results = []
for sign_name, qsign, Eanchor in anchor_curves:
    j_equal = (Eanchor.j_invariant() == Ecanon.j_invariant())
    isomorphic = False
    maps = None
    error = None
    try:
        iso = Eanchor.isomorphism_to(Ecanon)
        isomorphic = True
        maps = tuple(str(value) for value in iso.rational_maps())
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}"

    anchor_results.append({
        "sign": sign_name,
        "q": str(qsign),
        "j_equal": bool(j_equal),
        "isomorphic": bool(isomorphic),
        "isomorphism_rational_maps": maps,
        "error": error,
    })
    print(
        "Q8ZERO_ANCHOR|sign={}|j_equal={}|isomorphic={}|maps={}|error={}".format(
            sign_name, int(j_equal), int(isomorphic),
            "NA" if maps is None else ";".join(maps),
            "none" if error is None else error,
        ),
        flush=True,
    )

assert all(item["j_equal"] for item in anchor_results)
if anchor_type == "branch_point_q0":
    assert anchor_results[0]["isomorphic"]

status = (
    "PASS_EXACT_D13_ZERO_ANCHOR"
    if anchor_type == "branch_point_q0" and anchor_results[0]["isomorphic"]
    else "PASS_RATIONAL_ZERO_ANCHOR_SIGN_REMAINS"
)

payload = {
    "schema": "elkies-k3.h92-q8-d13-zero-anchor.v1",
    "status": status,
    "classifier_bridge": {
        "nef_witness_q6_raw": list(map(int, nef_witness_raw)),
        "nef_simple": list(map(int, nef_simple)),
        "dominant_simple": list(map(int, dominant_simple)),
        "reflections": [[int(i+1), int(p)] for i, p in bridge],
        "dominant_mw_projection": dominant_hit["mw_projection"],
    },
    "physical_target": {
        "reflection_count": len(physical_reflections),
        "endpoint_matches_target_artifact": bool(endpoint_match),
        "fibre_source_h3_ns": list(map(int, F_actual)),
    },
    "pinned_zero": {
        "source_h3_ns_after_physical_transport": list(map(int, O_actual)),
        "q6_component": component_name,
        "q6_fibre_kind": component_kind,
        "q6_base_value_T0": str(T0),
        "fibre_shift": fibre_shift,
    },
    "quartic_anchor": {
        "quartic_value_at_T0": str(q_at_zero),
        "rational": quartic_point_rational,
        "sqrt": None if W0 is None else str(W0),
        "type": anchor_type,
        "anchored_models": anchor_results,
    },
    "boundary": (
        "A branch-point PASS pins the abstract D13 zero to the standard "
        "Jacobian through an explicit anchored Weierstrass isomorphism. "
        "If the quartic value is a nonzero square, the two possible signs "
        "remain until a local resolved-component sign is selected."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"Q8ZERO_RESULT|status={status}")
print(f"OUTPUT|{OUTPUT}")
