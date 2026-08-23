#!/usr/bin/env sage -python
"""
Recover the standard-q6 MW coordinate of the transported old zero.

Exact standard-Weierstrass relations from the q8 marking:
    Pmap = E7_7 - old_zero     has MW coordinate (-2,0,0)
    Qmap = E7_7 - affine_E7   has MW coordinate (0,-2,0)
    S3                           MW coordinate (0,0,1)

Hence if z is the standard-group MW coordinate of old_zero,
    E7_7     = z + (-2,0,0)
    affine_E7= z + (-2,2,0).

For old_zero, E7_7 and affine_E7 the displayed rational coordinates determine
P.O exactly from the common square/cube denominator.  On E8+E6 the only local
height correction is at IV*, and is either 0 (identity component) or 4/3
(nonidentity multiplicity-one component).  The certified difference-section heights pin the E6 component-group map,
so the correction of any integral MW vector is determined modulo 3; then solve
for the unique integral z in the saturated height lattice

    [[8/3,1/3,-1],
     [1/3,8/3,1],
     [-1,1,46]].

Run:
  sage -python ~/Downloads/recover_h92_q6_standard_zero_translation.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, matrix, vector


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
parser.add_argument("--search-bound", type=int, default=40)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
CHILD = GEN / "elkies-k3-h92-q6-child-jacobian.json"
ZERO = GEN / "elkies-k3-h92-q6-child-zero-section.json"
COMP = GEN / "elkies-k3-h92-q6-child-e7-infinity-sections.json"
BRIDGE = LOCAL / "q6-third-to-q8-bridge.json"
OUTPUT = (
    args.output.resolve()
    if args.output and args.output.is_absolute()
    else ROOT / (
        args.output
        if args.output
        else Path("artifacts/local/elkies-k3/q6-standard-zero-translation.json")
    )
)

for path in (CHILD, ZERO, COMP, BRIDGE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

child = json.loads(CHILD.read_text())
zero = json.loads(ZERO.read_text())
components = json.loads(COMP.read_text())
bridge = json.loads(BRIDGE.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero["status"] == "PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"
assert components["status"] == "PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT"
assert bridge["status"] == "PASS_EXACT_Q6_THIRD_TO_Q8_DEGREE52"

R = PolynomialRing(QQ, "T")
T = R.gen()
K = R.fraction_field()


def poly(values):
    return R([QQ(v) for v in values])


def rat(data, np, dp):
    return K(poly(data[np])) / K(poly(data[dp]))


A = poly(child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B = poly(child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])

zdata = zero["section"]
old_zero = E(
    rat(
        zdata,
        "x_numerator_coefficients_low_to_high",
        "x_denominator_coefficients_low_to_high",
    ),
    rat(
        zdata,
        "y_numerator_coefficients_low_to_high",
        "y_denominator_coefficients_low_to_high",
    ),
)

points = {}
point_records = {}
for entry in components["sections"]:
    P = E(
        rat(
            entry,
            "x_numerator_coefficients_low_to_high",
            "x_denominator_coefficients_low_to_high",
        ),
        rat(
            entry,
            "y_numerator_coefficients_low_to_high",
            "y_denominator_coefficients_low_to_high",
        ),
    )
    points[entry["sign"]] = P
    point_records[entry["sign"]] = entry

affine_sign = components["source"]["affine_E7_sign"]
e77_sign = components["source"]["E7_7_sign"]
affine = points[affine_sign]
e77 = points[e77_sign]

# Exact group-law regressions from the corrected q8 marking.
Pmap = e77 - old_zero
Qmap = e77 - affine
assert not Pmap.is_zero() and not Qmap.is_zero()

H = matrix(QQ, [
    [QQ(8)/3, QQ(1)/3, -1],
    [QQ(1)/3, QQ(8)/3, 1],
    [-1, 1, 46],
])
assert vector(QQ, (-2,0,0))*H*vector(QQ, (-2,0,0)) == QQ(32)/3
assert vector(QQ, (0,-2,0))*H*vector(QQ, (0,-2,0)) == QQ(32)/3

iv_record = next(
    item for item in child["finite_fibres"] if item["kodaira"] == "IV*"
)
iv = R(iv_record["factor"]).monic()
assert iv.degree() == 1
u0 = -iv[0]/iv[1]


def monic_power_root(value, exponent):
    value = R(value)
    root = R.one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= factor.monic()**(multiplicity//exponent)
    return root.monic()


def section_pole_profile(P, label):
    x, y = P.xy()
    dx = R(x.denominator())
    dy = R(y.denominator())
    Z = monic_power_root(dx, 2)
    assert Z == monic_power_root(dy, 3)
    assert dx // Z**2 in QQ
    assert dy // Z**3 in QQ
    PO = ZZ(Z.degree())

    ox = int(R(x.numerator()).valuation(iv) - dx.valuation(iv))
    oy = int(R(y.numerator()).valuation(iv) - dy.valuation(iv))

    print(
        f"Q6ZEROTRANS_POLE|section={label}|O={PO}|"
        f"iv_orders={ox},{oy}|status=PASS",
        flush=True,
    )
    return {"O": PO, "iv_orders": [ox, oy]}


profiles = {
    "old_zero": section_pole_profile(old_zero, "old_zero"),
    "E7_7": section_pole_profile(e77, "E7_7"),
    "affine_E7": section_pole_profile(affine, "affine_E7"),
}

# The difference sections have certified height 32/3. Their exact pole
# degrees force correction 4/3 and P.O=4.
pmap_profile = section_pole_profile(Pmap, "Pmap=E7_7-old_zero")
qmap_profile = section_pole_profile(Qmap, "Qmap=E7_7-affine_E7")
pmap_correction = QQ(4 + 2*pmap_profile["O"]) - QQ(32)/3
qmap_correction = QQ(4 + 2*qmap_profile["O"]) - QQ(32)/3
assert pmap_profile["O"] == qmap_profile["O"] == 4
assert pmap_correction == qmap_correction == QQ(4)/3

# Pmap=-2 e1 and Qmap=-2 e2 are nonzero IV* component classes, while their
# sum (-2,-2,0) is the certified identity-component q8 marking. Therefore
# e1,e2 have opposite nonzero classes in E6^*/E6 = Z/3. The third basis
# direction has correction zero (height 46, P.O=21).
def iv_correction_from_mw(v):
    v = vector(ZZ, v)
    component = (v[0] - v[1]) % 3
    return QQ(0) if component == 0 else QQ(4)/3


assert iv_correction_from_mw((-2,0,0)) == QQ(4)/3
assert iv_correction_from_mw((0,-2,0)) == QQ(4)/3
assert iv_correction_from_mw((-2,-2,0)) == 0
assert iv_correction_from_mw((0,0,1)) == 0

# If z=coord(old_zero), then
#   coord(E7_7)      = z + (-2,0,0)
#   coord(affine_E7) = z + (-2,2,0).
candidates = []
bound = ZZ(args.search_bound)
for a in range(-bound, bound + 1):
    for b in range(-bound, bound + 1):
        for c in range(-bound, bound + 1):
            z = vector(ZZ, (a,b,c))
            e77v = z + vector(ZZ, (-2,0,0))
            affv = z + vector(ZZ, (-2,2,0))

            hz_expected = (
                QQ(4 + 2*profiles["old_zero"]["O"])
                - iv_correction_from_mw(z)
            )
            he77_expected = (
                QQ(4 + 2*profiles["E7_7"]["O"])
                - iv_correction_from_mw(e77v)
            )
            haff_expected = (
                QQ(4 + 2*profiles["affine_E7"]["O"])
                - iv_correction_from_mw(affv)
            )

            if vector(QQ,z)*H*vector(QQ,z) != hz_expected:
                continue
            if vector(QQ,e77v)*H*vector(QQ,e77v) != he77_expected:
                continue
            if vector(QQ,affv)*H*vector(QQ,affv) != haff_expected:
                continue
            candidates.append(vector(ZZ, (a,b,c)))

print(
    f"Q6ZEROTRANS_SOLVE|bound={bound}|candidates={len(candidates)}|"
    f"values={';'.join(','.join(map(str,v)) for v in candidates) if candidates else 'none'}",
    flush=True,
)

if len(candidates) != 1:
    raise ArithmeticError(
        f"standard-q6 old-zero MW coordinate not unique: {len(candidates)} candidates"
    )

z = candidates[0]

# More exact lattice/group consistency.
e77v = z + vector(ZZ, (-2,0,0))
affv = z + vector(ZZ, (-2,2,0))
hz = vector(QQ,z)*H*vector(QQ,z)
he77 = vector(QQ,e77v)*H*vector(QQ,e77v)
haff = vector(QQ,affv)*H*vector(QQ,affv)
assert hz == QQ(4 + 2*profiles["old_zero"]["O"]) - iv_correction_from_mw(z)
assert he77 == QQ(4 + 2*profiles["E7_7"]["O"]) - iv_correction_from_mw(e77v)
assert haff == QQ(4 + 2*profiles["affine_E7"]["O"]) - iv_correction_from_mw(affv)

# In the old-zero group, standard zero has coordinate -z.  The surface
# translation sending old_zero -> standard_zero is translation by -z in the
# standard Weierstrass group.
translation = -z

print(
    "Q6ZEROTRANS_SOLVED_PROFILE|"
    f"old_height={hz}|old_correction={iv_correction_from_mw(z)}|"
    f"E7_7_height={he77}|affine_height={haff}|status=PASS",
    flush=True,
)

print(
    "Q6ZEROTRANS_RESULT|"
    f"old_zero_standard_MW={','.join(map(str,z))}|"
    f"translate_old_to_standard={','.join(map(str,translation))}|"
    "status=PASS_EXACT_MW_TRANSLATION_PARAMETER",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-q6-standard-zero-translation.v1",
    "status": "PASS_EXACT_MW_TRANSLATION_PARAMETER",
    "height_gram": [[str(v) for v in row] for row in H.rows()],
    "standard_MW_coordinates": {
        "old_zero": list(map(int, z)),
        "E7_7": list(map(int, z + vector(ZZ, (-2,0,0)))),
        "affine_E7": list(map(int, z + vector(ZZ, (-2,2,0)))),
        "Pmap": [-2,0,0],
        "Qmap": [0,-2,0],
        "S3": [0,0,1],
    },
    "profiles": {
        "old_zero": {
            "O_intersection": int(profiles["old_zero"]["O"]),
            "IVstar_orders": profiles["old_zero"]["iv_orders"],
            "IVstar_correction": str(iv_correction_from_mw(z)),
            "height": str(hz),
        },
        "E7_7": {
            "O_intersection": int(profiles["E7_7"]["O"]),
            "IVstar_orders": profiles["E7_7"]["iv_orders"],
            "IVstar_correction": str(iv_correction_from_mw(e77v)),
            "height": str(he77),
        },
        "affine_E7": {
            "O_intersection": int(profiles["affine_E7"]["O"]),
            "IVstar_orders": profiles["affine_E7"]["iv_orders"],
            "IVstar_correction": str(iv_correction_from_mw(affv)),
            "height": str(haff),
        },
    },
    "translation": {
        "surface_automorphism": "fiberwise translation sending old_zero to standard Weierstrass zero",
        "standard_group_translation_vector": list(map(int, translation)),
        "old_group_coordinate_of_standard_zero": list(map(int, translation)),
    },
    "boundary": (
        "This certifies the MW translation parameter only.  Applying the "
        "translation to the full Neron-Severi lattice, including the IV* "
        "component permutation, is the next certificate."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
