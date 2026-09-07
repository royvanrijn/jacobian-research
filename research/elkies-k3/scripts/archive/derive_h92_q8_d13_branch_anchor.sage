#!/usr/bin/env sage -python
"""
Anchor the corrected H92 q8 quartic at the actual II*_E8_1 branch point.

The effective-zero audit proved:
  * II*_E8_1 is a genuine q8 section;
  * its old-base value T=T_II is a branch point of the exact q8 quartic;
  * IV*_E6_1 and IV*_E6_5 are the other two old-fibre component sections.

Choosing II*_E8_1 as zero removes all sign/origin ambiguity.  If

    W^2 = a r^4 + b r^3 + c r^2 + d r,   r=T-T_II,

then

    X0 = d/r,
    Y0 = d W/r^2

gives

    Y0^2 = X0^3 + c X0^2 + b d X0 + a d^2,

and sends the branch point (r,W)=(0,0) to infinity.

This script:
  1. certifies that branch-point transform exactly;
  2. constructs the anchored elliptic curve over QQ(U);
  3. obtains an exact Sage Weierstrass isomorphism to the certified minimal
     I9*+9I1 D13 child;
  4. records explicit quartic -> canonical-D13 rational maps;
  5. maps the two IV* points W=+/-sqrt(f(T_IV)) to canonical D13 points;
  6. independently computes the MW coordinates of the actual E6_1/E6_5
     component sections relative to E8_1 zero in a deterministic D13 frame.

Expected structural result:
    E6_5 = -E6_1 in MW,
    height(E6_1)=3/4,
so this gives one explicit primitive D13 MW generator, up to the harmless
choice of which IV* branch is called plus.

Run:
  sage -python ~/Downloads/derive_h92_q8_d13_branch_anchor.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    EllipticCurve, PolynomialRing, QQ, ZZ, block_diagonal_matrix,
    identity_matrix, lcm, matrix, pari, sage_eval, vector
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


def child_frame_with_zero(ns, fibre, zero):
    assert fibre * ns * fibre == 0
    assert zero * ns * zero == -2
    assert zero * ns * fibre == 1
    mate = zero + fibre
    assert mate * ns * mate == 0
    assert mate * ns * fibre == 1
    complement = matrix(
        ZZ, [list(fibre * ns), list(mate * ns)]
    ).right_kernel_matrix()
    basis = matrix(
        ZZ,
        [list(fibre), list(mate)] + [list(row) for row in complement.rows()],
    )
    assert abs(basis.det()) == 1
    child = -(complement * ns * complement.transpose())
    U2 = matrix(ZZ, ((0, 1), (1, 0)))
    assert (
        basis * ns * basis.transpose()
        == block_diagonal_matrix(U2, -child)
    )
    return child, basis


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    if not count:
        return (), matrix(ZZ, 0, gram.nrows()), (0, 0, 1)
    half = [vector(ZZ, column) for column in matrix(ZZ, result[2]).columns()]
    roots = tuple(half + [-root for root in half])
    root_basis = matrix(
        ZZ, [list(root) for root in roots]
    ).row_module().basis_matrix()
    root_gram = root_basis * gram * root_basis.transpose()
    return roots, root_basis, (
        root_basis.rank(),
        count,
        abs(ZZ(root_gram.det())),
    )


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    rank = data[0]
    regular = None
    for shift in range(1, 1000):
        candidate = vector(
            ZZ,
            [(i + 1) ** 2 + shift * (i + 1) + 1
             for i in range(gram.nrows())],
        )
        if all(candidate * root != 0 for root in roots):
            regular = candidate
            break
    assert regular is not None
    positive = [root for root in roots if regular * root > 0]
    positive_set = {tuple(root) for root in positive}
    simple = [
        root
        for root in positive
        if not any(tuple(root - left) in positive_set for left in positive)
    ]
    simple = matrix(ZZ, [list(root) for root in simple])
    assert simple.nrows() == simple.rank() == rank
    cartan = simple * gram * simple.transpose()
    return simple, cartan


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

    quotient_change = block_diagonal_matrix(
        identity_matrix(ZZ, 13), lll.transpose()
    )
    basis = quotient_change * initial
    adapted = basis * child * basis.transpose()

    root = adapted[:13, :13]
    coupling = adapted[:13, 13:]
    tail = adapted[13:, 13:]
    H = tail - coupling.transpose() * root.inverse() * coupling
    return basis, adapted, H


def rat_degrees(value, ring):
    value = value.parent()(value)
    return (
        int(ring(value.numerator()).degree()),
        int(ring(value.denominator()).degree()),
    )


def point_payload(point, ring):
    if point.is_zero():
        return {"zero": True}
    x, y = point.xy()
    return {
        "zero": False,
        "x": str(x),
        "y": str(y),
        "x_degrees": list(rat_degrees(x, ring)),
        "y_degrees": list(rat_degrees(y, ring)),
    }


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
Q6_CHILD = GEN / "elkies-k3-h92-q6-child-jacobian.json"

TARGET = (
    args.target.resolve()
    if args.target
    else LOCAL / "q8-target-component-nef.json"
)

if args.q8_child:
    Q8_CHILD = args.q8_child.resolve()
else:
    candidates = [
        GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
        LOCAL / "q8-corrected2cover-qq-child.json",
    ]
    Q8_CHILD = next((path for path in candidates if path.exists()), candidates[0])

OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q8-d13-branch-anchor.json"
)

for path in (FRAME, Q6_CHILD, TARGET, Q8_CHILD):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

source_frame = load_gram(FRAME)
q6child = json.loads(Q6_CHILD.read_text())
target = json.loads(TARGET.read_text())
q8child = json.loads(Q8_CHILD.read_text())

assert q6child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert target["normalization"]["representative"] == "component-nef"
assert q8child["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"

U2 = matrix(ZZ, ((0, 1), (1, 0)))
source_ns = block_diagonal_matrix(U2, -source_frame)

F8 = vector(ZZ, target["selected_q8"]["source_h3_ns_vector"])
E6 = matrix(
    ZZ,
    target["selected_q8"]["E6"]["simple_root_vectors_in_source_h3_ns"],
)
E8 = matrix(
    ZZ,
    target["selected_q8"]["E8"]["simple_root_vectors_in_source_h3_ns"],
)

O = vector(ZZ, E8.row(0))  # II*_E8_1
Piv1 = vector(ZZ, E6.row(0))  # IV*_E6_1
Piv5 = vector(ZZ, E6.row(4))  # IV*_E6_5

for curve in (O, Piv1, Piv5):
    assert curve * source_ns * curve == -2
    assert curve * source_ns * F8 == 1

assert Piv1 * source_ns * O == 0
assert Piv5 * source_ns * O == 0

# ---------------------------------------------------------------------------
# 1. Deterministic D13 lattice frame for the effective branch-point zero.
# ---------------------------------------------------------------------------

child, Bzero = child_frame_with_zero(source_ns, F8, O)
assert roots_and_data(child)[2] == (13, 312, 4)
A13, adapted, H = d13_root_adaptation(child)

Badapt = block_diagonal_matrix(identity_matrix(ZZ, 2), A13) * Bzero
Gadapt = block_diagonal_matrix(U2, -adapted)
assert Badapt * source_ns * Badapt.transpose() == Gadapt

def coordinates(curve):
    result = vector(QQ, curve) * Badapt.inverse()
    assert all(value in ZZ for value in result)
    return vector(ZZ, result)

c1 = coordinates(Piv1)
c5 = coordinates(Piv5)
assert c1[1] == c5[1] == 1

z1 = vector(ZZ, c1[-4:])
z5 = vector(ZZ, c5[-4:])
h1 = QQ(z1 * H * z1)
h5 = QQ(z5 * H * z5)

assert h1 == h5 == QQ(3) / 4
assert z5 == -z1

# With P.O=0, Shioda gives local correction 4-height = 13/4.
assert QQ(4) - h1 == QQ(13) / 4

print(
    "Q8BRANCH_LATTICE|zero=II*_E8_1|"
    f"E6_1_mw={','.join(map(str,z1))}|"
    f"E6_5_mw={','.join(map(str,z5))}|"
    "inverse=1|height=3/4|PdotO=0|correction=13/4|status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# 2. Exact q8 quartic and branch-point transformation.
# ---------------------------------------------------------------------------

oldR = PolynomialRing(QQ, "Told")
Told = oldR.gen()

def old_fibre_value(kind):
    text = next(
        item["factor"]
        for item in q6child["finite_fibres"]
        if item["kodaira"] == kind
    )
    factor = oldR(sage_eval(text, locals={"T": Told, "Told": Told}))
    assert factor.degree() == 1
    return QQ(-factor[0] / factor[1])

Tii = old_fibre_value("II*")
Tiv = old_fibre_value("IV*")
assert Tii != Tiv

UR = PolynomialRing(QQ, "U")
U = UR.gen()
UF = UR.fraction_field()

TR = PolynomialRing(UF, "T")
T = TR.gen()
quartic = TR(
    sage_eval(q8child["pencil"]["branch_quartic"], locals={"U": U, "T": T})
)
assert quartic.degree() == 4

q_ii = UF(quartic(UF(Tii)))
q_iv = UF(quartic(UF(Tiv)))
assert not q_ii
assert q_iv and q_iv.is_square()
Wiv = q_iv.sqrt()

# Shift T=Tii+r.
rR = PolynomialRing(UF, "r")
r = rR.gen()
shifted = rR(quartic(r + UF(Tii)))
assert shifted[0] == 0
assert shifted.degree() == 4

d = UF(shifted[1])
c = UF(shifted[2])
b = UF(shifted[3])
a = UF(shifted[4])
assert d

Eanchor = EllipticCurve(UF, [0, c, 0, b*d, a*d**2])

# Verify the quartic transform algebraically in a polynomial ring in r,W.
# It is enough to verify the coefficient identity:
# Y0^2 = X0^3+cX0^2+bdX0+ad^2 after substitution X0=d/r.
assert shifted == rR(d*r + c*r**2 + b*r**3 + a*r**4)

child_data = q8child["child"]
Amin = UR([QQ(v) for v in child_data["minimal_A_coefficients_low_to_high"]])
Bmin = UR([QQ(v) for v in child_data["minimal_B_coefficients_low_to_high"]])
Ecanon = EllipticCurve(UF, [0, 0, 0, UF(Amin), UF(Bmin)])

assert Eanchor.j_invariant() == Ecanon.j_invariant()
assert Eanchor.is_isomorphic(Ecanon)

iso = Eanchor.isomorphism_to(Ecanon)
u_iso, r_iso, s_iso, t_iso = iso.tuple()
map_x, map_y = iso.rational_maps()

# Explicit combined quartic -> canonical formulas:
#   Xa = d/(T-Tii)
#   Ya = d*W/(T-Tii)^2
# followed by (u,r,s,t).
# Store as strings rather than attempting to construct a two-variable field.
combined_x = (
    f"(({d})/(T-({Tii})) - ({r_iso}))/({u_iso})^2"
)
combined_y = (
    f"(({d})*W/(T-({Tii}))^2"
    f" - ({s_iso})*((({d})/(T-({Tii})))-({r_iso}))"
    f" - ({t_iso}))/({u_iso})^3"
)

print(
    "Q8BRANCH_ANCHOR|T_II={}|T_IV={}|II_branch=1|IV_square=1|"
    "j_equal=1|isomorphic=1|urst=({},{},{},{})|status=PASS".format(
        Tii, Tiv, u_iso, r_iso, s_iso, t_iso
    ),
    flush=True,
)

# ---------------------------------------------------------------------------
# 3. Map the two IV* quartic points to the canonical D13 equation.
# ---------------------------------------------------------------------------

deltaT = UF(Tiv - Tii)
Xa = d / deltaT
Ya_plus = d * Wiv / deltaT**2

Pplus_anchor = Eanchor(Xa, Ya_plus)
Pminus_anchor = -Pplus_anchor
assert Pminus_anchor.xy()[0] == Pplus_anchor.xy()[0]

Pplus = iso(Pplus_anchor)
Pminus = iso(Pminus_anchor)

assert Pplus in Ecanon and Pminus in Ecanon
assert Pminus == -Pplus
assert not Pplus.is_zero()

xp, yp = Pplus.xy()
xm, ym = Pminus.xy()
assert xp == xm
assert ym == -yp  # canonical child is short Weierstrass

print(
    "Q8BRANCH_POINTS|IV_plus_minus_inverse=1|"
    f"x_deg={rat_degrees(xp,UR)[0]}/{rat_degrees(xp,UR)[1]}|"
    f"y_deg={rat_degrees(yp,UR)[0]}/{rat_degrees(yp,UR)[1]}|"
    f"lattice_mw_plus_choice={','.join(map(str,z1))}|"
    f"lattice_mw_minus_choice={','.join(map(str,z5))}|status=PASS",
    flush=True,
)

# Group-law sanity checks.
assert 2 * Pminus == -(2 * Pplus)
assert Pplus + Pminus == Ecanon(0)

payload = {
    "schema": "elkies-k3.h92-q8-d13-branch-anchor.v1",
    "status": "PASS_EXACT_D13_BRANCH_ANCHOR",
    "zero": {
        "curve": "II*_E8_1",
        "old_base_T": str(Tii),
        "quartic_point": [str(Tii), "0"],
        "interpretation": (
            "actual q6 fibre component, q8 section, chosen as q8 zero; "
            "the quartic branch point maps to infinity"
        ),
    },
    "iv_sections": {
        "old_base_T": str(Tiv),
        "quartic_value": str(q_iv),
        "quartic_sqrt": str(Wiv),
        "lattice": {
            "E6_1_mw": list(map(int, z1)),
            "E6_5_mw": list(map(int, z5)),
            "inverse": True,
            "height": "3/4",
            "P_dot_O": 0,
            "D13_local_correction": "13/4",
            "height_gram": [[str(v) for v in row] for row in H.rows()],
        },
        "canonical_D13_points": {
            "plus": point_payload(Pplus, UR),
            "minus": point_payload(Pminus, UR),
            "sign_assignment": (
                "The algebraic plus branch is assigned to E6_1 and minus to "
                "E6_5; swapping both is the harmless choice of generator sign."
            ),
        },
    },
    "quartic_to_anchor": {
        "shift": f"r=T-({Tii})",
        "shifted_coefficients": {
            "a_r4": str(a),
            "b_r3": str(b),
            "c_r2": str(c),
            "d_r1": str(d),
        },
        "X0": f"({d})/(T-({Tii}))",
        "Y0": f"({d})*W/(T-({Tii}))^2",
        "weierstrass": {
            "a1": "0",
            "a2": str(c),
            "a3": "0",
            "a4": str(b*d),
            "a6": str(a*d**2),
        },
    },
    "anchor_to_canonical": {
        "urst": [str(u_iso), str(r_iso), str(s_iso), str(t_iso)],
        "sage_rational_maps": [str(map_x), str(map_y)],
        "combined_quartic_to_canonical": {
            "x": combined_x,
            "y": combined_y,
        },
    },
    "canonical_child": {
        "minimal_A_coefficients_low_to_high":
            child_data["minimal_A_coefficients_low_to_high"],
        "minimal_B_coefficients_low_to_high":
            child_data["minimal_B_coefficients_low_to_high"],
        "fibre_signature": "I9* + 9 I1",
        "root_lattice": "D13",
    },
    "boundary": (
        "This pins the q8 origin and one primitive height-3/4 MW generator "
        "explicitly on the canonical D13 equation. It does not yet recover "
        "the remaining three MW generators or execute the q24 neighbour."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("Q8BRANCH_RESULT|status=PASS_EXACT_D13_BRANCH_ANCHOR")
print(f"OUTPUT|{OUTPUT}")
