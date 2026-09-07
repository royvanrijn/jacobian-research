#!/usr/bin/env sage -python
"""
Construct the current cheapest H92 q24 bridge as an EXACT q6 rational section.

Recovered exhaustive bridge data:
    raw-q6 MW word n = (0,-1,1)
    q8 degree          = 46
    D13 AJ             = (-1,-1,0,1)
    q24 correction     = +3*G1 - G3

The important point is to identify the raw-q6 MW basis with the older
H3-lift q6 basis *exactly*, rather than assuming their coordinates agree.

Once that change of basis is certified, convert the selected old-zero-group
word to the standard q6 Weierstrass group using the already-certified old-zero
translation.  The resulting standard MW vector is then built from the exact
known q6 points

    Pmap = E7_7 - old_zero = (-2,0,0)
    Qmap = E7_7 - affine   = (0,-2,0)
    S3                     = (0,0,1).

Finally replay the repaired q8 RR pencil and prove that its restriction to the
constructed section has degree 46.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    EllipticCurve, PolynomialRing, QQ, ZZ, block_diagonal_matrix,
    identity_matrix, matrix, pari, vector
)


Q6_REFLECTIONS = (
    1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3, 1,
    7, 6, 5, 4, 2, 3, 4, 5, 6, 7,
)

H3_LIFTS = matrix(ZZ, [
    [-5,-4,-3,0,0,0,0,0,0,0,0,-4,1,0,-4,2,-2],
    [-10,-8,-6,0,0,0,0,0,0,0,0,-8,4,1,-8,5,-4],
    [-5,-4,-3,0,0,0,0,0,0,0,0,-3,2,0,-4,2,-2],
])

OLD_ZERO_ROOT_SHIFTS = matrix(ZZ, [
    [5,4,3,0,0,0,0,0,0,0,0,3,-1,4],
    [12,10,8,0,0,0,0,0,0,0,0,6,-1,9],
    [5,4,3,0,0,0,0,0,0,0,0,2,0,4],
])

EXPECTED_H = matrix(QQ, [
    [QQ(8)/3, QQ(1)/3, -1],
    [QQ(1)/3, QQ(8)/3, 1],
    [-1, 1, 46],
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
    for c in candidates:
        try:
            c = c.resolve()
        except Exception:
            continue
        if c in seen:
            continue
        seen.add(c)
        if (
            (c / "elkies-k3/scripts").is_dir()
            and (c / "artifacts/generated-results").is_dir()
        ):
            return c
    raise SystemExit("Could not locate jacobian-research")


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


def polynomial(ring, values):
    return ring([QQ(v) for v in values])


def rational(field, ring, data, nk, dk):
    return field(polynomial(ring, data[nk])) / field(polynomial(ring, data[dk]))


def rf_record(value, ring):
    value = value.parent()(value)
    return {
        "numerator_coefficients_low_to_high": [
            str(v) for v in ring(value.numerator()).list()
        ],
        "denominator_coefficients_low_to_high": [
            str(v) for v in ring(value.denominator()).list()
        ],
    }


def monic_power_root(value, exponent):
    root = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= factor.monic() ** (multiplicity // exponent)
    return root.monic()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"

FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
ORBITS = GEN / "elkies-k3-h3-q6-q8-orbits.json"
CHILD = GEN / "elkies-k3-h92-q6-child-jacobian.json"
ZERO = GEN / "elkies-k3-h92-q6-child-zero-section.json"
COMP = GEN / "elkies-k3-h92-q6-child-e7-infinity-sections.json"
S3BRIDGE = LOCAL / "q6-third-to-q8-bridge.json"
TRANSLATION = LOCAL / "q6-standard-zero-translation.json"
CHEAP = LOCAL / "q24-cheapest-bridge-current.json"

q8_candidates = [
    LOCAL / "q8-corrected2cover-qq-child.json",
    GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8 = next(
    (
        p for p in q8_candidates
        if p.exists()
        and "rr" in json.loads(p.read_text())
        and "kernel_polynomials" in json.loads(p.read_text()).get("rr", {})
        and "marking" in json.loads(p.read_text())
    ),
    None,
)
if Q8 is None:
    raise SystemExit("No complete corrected q8 child artifact found")

OUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q24-cheapest-bridge-q6-exact.json"
)

for p in (FRAME, ORBITS, CHILD, ZERO, COMP, S3BRIDGE, TRANSLATION, CHEAP, Q8):
    if not p.exists():
        raise SystemExit(f"Missing prerequisite: {p}")

orbits = json.loads(ORBITS.read_text())
child = json.loads(CHILD.read_text())
zero = json.loads(ZERO.read_text())
components = json.loads(COMP.read_text())
s3bridge = json.loads(S3BRIDGE.read_text())
translation = json.loads(TRANSLATION.read_text())
cheap = json.loads(CHEAP.read_text())
q8 = json.loads(Q8.read_text())

assert orbits["status"] == "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION"
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero["status"] == "PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"
assert components["status"] == "PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT"
assert s3bridge["status"] == "PASS_EXACT_Q6_THIRD_TO_Q8_DEGREE52"
assert translation["status"] == "PASS_EXACT_MW_TRANSLATION_PARAMETER"
assert cheap["status"] == "PASS_PINNED_CHEAPEST_Q24_BRIDGE_BELOW_47"
assert q8["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"

selected = cheap["selected"]
assert selected["kind"] == "raw_q6"
n_raw = vector(ZZ, selected["word"])
assert n_raw == vector(ZZ, (0, -1, 1))
assert int(selected["q8_degree"]) == 46
assert selected["d13_mw"] == [-1, -1, 0, 1]
assert int(selected["q24_correction_G1"]) == 3
assert int(selected["q24_correction_G3"]) == -1

# ===========================================================================
# 1. Exact change of q6 MW basis:
#    older H3-lift basis -> stored raw root/MW basis.
# ===========================================================================

frame = load_gram(FRAME)
U2 = matrix(ZZ, ((0,1),(1,0)))
source_ns = block_diagonal_matrix(U2, -frame)
source_O = vector(ZZ, [-1,1] + [0]*17)
source_simple = tuple(
    vector(ZZ, [0,0] + [ZZ(i == node) for i in range(17)])
    for node in range(15)
)
reflection_roots = tuple(source_simple[node-1] for node in Q6_REFLECTIONS)

F6_raw = vector(ZZ, [3,2] + [
    0,0,-1,-1,-1,-1,-1,
    0,0,0,0,0,0,0,0,1,0,
])
O6_raw = vector(ZZ, source_O)
for root in reversed(reflection_roots):
    O6_raw = reflect(O6_raw, source_ns, root)

# Stored raw child frame.
B6 = matrix(ZZ, orbits["q6"]["neighbor_basis_in_source_ns"])
raw_ns = B6 * source_ns * B6.transpose()
raw_child = -raw_ns[2:,2:]
assert raw_child.det() == 948

# Regression: the explicit raw fibre is the first stored neighbour row.
assert vector(ZZ, B6.row(0)) == F6_raw

root_mw = matrix(ZZ, orbits["q6"]["root_mw_basis_in_child"])
assert abs(root_mw.det()) == 1
root_mw_frame = root_mw * raw_child * root_mw.transpose()
assert root_mw_frame == matrix(ZZ, orbits["q6"]["root_adapted_gram"])
root_block = root_mw_frame[:14,:14]
coupling = root_mw_frame[:14,14:]
tail = root_mw_frame[14:,14:]
H_raw = tail - coupling.transpose()*root_block.inverse()*coupling
assert H_raw == EXPECTED_H

raw_roots = matrix(
    ZZ, pari(raw_child).qfminim(2)[2]
).transpose().row_module().basis_matrix()
assert raw_roots.nrows() == 14

# The H3-lift generator difference from the old zero is exactly
# h3_lift + root_shift. Express each such direction in the stored root/MW basis.
old_to_raw_rows = []
for h3_lift, shift in zip(H3_LIFTS.rows(), OLD_ZERO_ROOT_SHIFTS.rows()):
    direction = vector(ZZ, h3_lift) + vector(ZZ, shift) * raw_roots
    c = vector(QQ, direction) * root_mw.inverse()
    assert all(v in ZZ for v in c)
    c = vector(ZZ, c)
    old_to_raw_rows.append(vector(ZZ, c[-3:]))

M_old_to_raw = matrix(ZZ, [list(v) for v in old_to_raw_rows])
assert abs(M_old_to_raw.det()) == 1
assert M_old_to_raw * EXPECTED_H * M_old_to_raw.transpose() == EXPECTED_H

n_old_q = vector(QQ, n_raw) * M_old_to_raw.inverse()
assert all(v in ZZ for v in n_old_q)
n_old = vector(ZZ, n_old_q)
assert n_old * M_old_to_raw == n_raw

# Old-group -> standard-Weierstrass group:
# an old-group word w represents standard coordinate w + z(old_zero).
z_oldzero_std = vector(
    ZZ, translation["standard_MW_coordinates"]["old_zero"]
)
assert z_oldzero_std == vector(ZZ, (2,-1,0))
n_std = n_old + z_oldzero_std

print(
    "Q24W_BASIS|"
    f"old_to_raw={';'.join(','.join(map(str,row)) for row in M_old_to_raw.rows())}|"
    f"raw_word={','.join(map(str,n_raw))}|"
    f"old_word={','.join(map(str,n_old))}|"
    f"standard_mw={','.join(map(str,n_std))}|"
    "status=PASS_EXACT_BASIS_CHANGE",
    flush=True,
)

# ===========================================================================
# 2. Build that exact standard-q6 point from known rational points.
# ===========================================================================

R = PolynomialRing(QQ, "T")
T = R.gen()
K = R.fraction_field()
model = child["minimal_short_weierstrass"]
A = polynomial(R, model["A_coefficients_low_to_high"])
Bcurve = polynomial(R, model["B_coefficients_low_to_high"])
E = EllipticCurve(K, [0,0,0,K(A),K(Bcurve)])

zdata = zero["section"]
old_zero_point = E(
    rational(K,R,zdata,
             "x_numerator_coefficients_low_to_high",
             "x_denominator_coefficients_low_to_high"),
    rational(K,R,zdata,
             "y_numerator_coefficients_low_to_high",
             "y_denominator_coefficients_low_to_high"),
)

points = {}
for entry in components["sections"]:
    points[entry["sign"]] = E(
        rational(K,R,entry,
                 "x_numerator_coefficients_low_to_high",
                 "x_denominator_coefficients_low_to_high"),
        rational(K,R,entry,
                 "y_numerator_coefficients_low_to_high",
                 "y_denominator_coefficients_low_to_high"),
    )
affine = points[components["source"]["affine_E7_sign"]]
e77 = points[components["source"]["E7_7_sign"]]

Pmap = e77 - old_zero_point
Qmap = e77 - affine

s3data = s3bridge["third_section_canonical_q6"]
S3 = E(
    rational(K,R,s3data["x"],
             "numerator_coefficients_low_to_high",
             "denominator_coefficients_low_to_high"),
    rational(K,R,s3data["y"],
             "numerator_coefficients_low_to_high",
             "denominator_coefficients_low_to_high"),
)

assert translation["standard_MW_coordinates"]["Pmap"] == [-2,0,0]
assert translation["standard_MW_coordinates"]["Qmap"] == [0,-2,0]
assert translation["standard_MW_coordinates"]["S3"] == [0,0,1]
assert s3bridge["lattice_target"]["q6_MW_word"] == [0,0,1]

# The exact points available span all vectors with even first/second coordinate.
a,b,c = map(ZZ, n_std)
if a % 2 or b % 2:
    raise ArithmeticError(
        f"selected standard q6 word {tuple(n_std)} is not expressible "
        "from exact Pmap,Qmap,S3 without primitive halving"
    )
kp = -a // 2
kq = -b // 2
ks = c

Wphysical = kp*Pmap + kq*Qmap + ks*S3
assert Wphysical in E and not Wphysical.is_zero()

wpx, wpy = Wphysical.xy()
assert wpy**2 == wpx**3 + K(A)*wpx + K(Bcurve)

physical_formula = f"{kp}*Pmap+{kq}*Qmap+{ks}*S3"
pretty_physical_formula = (
    "-Pmap+Qmap+S3"
    if (kp,kq,ks) == (-1,1,1)
    else physical_formula
)

# CRITICAL repaired-route step:
# the raw bridge search used the PHYSICAL/component-nef q8 fibre.
# The repaired RR pencil uses the equation q8 fibre, obtained by the exact
# q6 translation sending old_zero to the standard zero. On the standard
# q6 Weierstrass group this is translation by -old_zero_point.
Weq = Wphysical - old_zero_point
assert Weq in E and not Weq.is_zero()

# In the expected basis orientation this simplifies to E7_7 + S3.
assert Weq == e77 + S3

wx, wy = Weq.xy()
assert wy**2 == wx**3 + K(A)*wx + K(Bcurve)

formula = "translate_old_to_standard(Wphysical)"
pretty_formula = "E7_7+S3"

def degrees(value):
    value = K(value)
    return (
        int(R(value.numerator()).degree()),
        int(R(value.denominator()).degree()),
    )

print(
    "Q24W_Q6_PHYSICAL|"
    f"raw_word={','.join(map(str,n_raw))}|"
    f"old_word={','.join(map(str,n_old))}|"
    f"standard_mw_before_translation={','.join(map(str,n_std))}|"
    f"formula={pretty_physical_formula}|"
    f"x={degrees(wpx)[0]}/{degrees(wpx)[1]}|"
    f"y={degrees(wpy)[0]}/{degrees(wpy)[1]}|"
    "identity=PASS|status=PASS_EXACT_PHYSICAL_BRIDGE",
    flush=True,
)
print(
    "Q24W_Q6_EQUATION|"
    "translation=-old_zero|standard_mw_after_translation=0,-1,1|"
    f"formula={pretty_formula}|"
    f"x={degrees(wx)[0]}/{degrees(wx)[1]}|"
    f"y={degrees(wy)[0]}/{degrees(wy)[1]}|"
    "identity=PASS|status=PASS_EXACT_EQUATION_BRIDGE",
    flush=True,
)

# ===========================================================================
# 3. Replay the repaired q8 pencil on W and certify degree 46.
# ===========================================================================

# Use the serialized corrected marked section, and independently check it equals
# the actual Pmap+Qmap we just reconstructed.
mdata = q8["marking"]["section"]
sx = rational(
    K,R,mdata,
    "x_numerator_coefficients_low_to_high",
    "x_denominator_coefficients_low_to_high",
)
sy = rational(
    K,R,mdata,
    "y_numerator_coefficients_low_to_high",
    "y_denominator_coefficients_low_to_high",
)
Smark = E(sx,sy)
assert Smark == Pmap + Qmap

nx,dx = R(sx.numerator()), R(sx.denominator())
ny,dy = R(sy.numerator()), R(sy.denominator())
h = monic_power_root(dx,2)
assert h == monic_power_root(dy,3)
assert h.degree() == 10

ii = R(next(item for item in child["finite_fibres"] if item["kodaira"]=="II*")["factor"]).monic()
iv = R(next(item for item in child["finite_fibres"] if item["kodaira"]=="IV*")["factor"]).monic()
M = (ii**2 * iv**2).monic()

normalizer = (ny * dx * (h*dy).inverse_mod(nx)).mod(nx)
assert (normalizer*h*dy - ny*dx) % nx == 0
p_fun = -sy/sx
rho = (normalizer * nx.inverse_mod(M)).mod(M)

pairs = []
for entry in q8["rr"]["kernel_polynomials"]:
    sp = R(entry["s"])
    tp = R(entry["t"])
    Bcoef = K(sp)/K(h)
    Acoef = (
        -K(sp)*p_fun/K(h)
        - K(sp)*K(normalizer)/K(nx)
        + K(sp*rho)
        + K(tp*M)
    )
    pairs.append((Acoef,Bcoef))
assert len(pairs) == 2
(A0,B0),(A1,B1) = pairs

mW = (wy + sy)/(wx - sx)
UW = K((A1+B1*mW)/(A0+B0*mW))
u_num = R(UW.numerator())
u_den = R(UW.denominator())
common = u_num.gcd(u_den)
assert common in QQ
q8_degree = max(u_num.degree(),u_den.degree())
assert q8_degree == 46, q8_degree

print(
    "Q24W_Q8|"
    f"degree={q8_degree}|numdeg={u_num.degree()}|dendeg={u_den.degree()}|"
    f"expected_AJ={','.join(map(str,selected['d13_mw']))}|"
    "q24_add=3*G1-G3|status=PASS_EXACT_DEGREE46_BRIDGE",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-q24-cheapest-q6-bridge-exact.v1",
    "status": "PASS_EXACT_Q24_CHEAPEST_Q6_BRIDGE_DEGREE46",
    "basis_certificate": {
        "raw_search_word": list(map(int,n_raw)),
        "old_H3_lift_to_raw_MW_matrix_rows": [
            list(map(int,row)) for row in M_old_to_raw.rows()
        ],
        "old_zero_group_word": list(map(int,n_old)),
        "old_zero_standard_MW": list(map(int,z_oldzero_std)),
        "standard_q6_MW": list(map(int,n_std)),
    },
    "physical_bridge_q6_point": {
        "formula": pretty_physical_formula,
        "standard_MW_before_translation": list(map(int,n_std)),
        "coefficients_on_Pmap_Qmap_S3": [int(kp),int(kq),int(ks)],
        "x": rf_record(wpx,R),
        "y": rf_record(wpy,R),
        "weierstrass_identity": True,
        "x_degrees": list(degrees(wpx)),
        "y_degrees": list(degrees(wpy)),
    },
    "equation_bridge_q6_point": {
        "translation": "-old_zero_point",
        "formula": pretty_formula,
        "standard_MW_after_translation": [0,-1,1],
        "exact_identity": "W_eq = W_physical - old_zero = E7_7 + S3",
        "x": rf_record(wx,R),
        "y": rf_record(wy,R),
        "weierstrass_identity": True,
        "x_degrees": list(degrees(wx)),
        "y_degrees": list(degrees(wy)),
    },
    # Deliberately retain aliases matching the old S3 bridge shape. These
    # are the EQUATION-translated bridge coordinates, suitable for the repaired
    # q8 pencil and direct branch-zero trace.
    "third_section_canonical_q6": {
        "x": rf_record(wx,R),
        "y": rf_record(wy,R),
        "weierstrass_identity": True,
    },
    "q8_parameter_on_third": {
        "degree": int(q8_degree),
        "numerator_degree": int(u_num.degree()),
        "denominator_degree": int(u_den.degree()),
        "numerator_coefficients_low_to_high": [str(v) for v in u_num.list()],
        "denominator_coefficients_low_to_high": [str(v) for v in u_den.list()],
    },
    "q8_bridge": {
        "degree": int(q8_degree),
        "numerator_degree": int(u_num.degree()),
        "denominator_degree": int(u_den.degree()),
        "numerator_coefficients_low_to_high": [str(v) for v in u_num.list()],
        "denominator_coefficients_low_to_high": [str(v) for v in u_den.list()],
    },
    "D13_target": {
        "AJ_of_bridge": list(map(int,selected["d13_mw"])),
        "q24": [2,-1,-1,1],
        "correction": {
            "G1": int(selected["q24_correction_G1"]),
            "G3": int(selected["q24_correction_G3"]),
            "formula": "q24 = AJ(W) + 3*G1 - G3",
        },
    },
    "inputs": {
        "q8_child": str(Q8.relative_to(ROOT)),
        "q8_child_sha256": hashlib.sha256(Q8.read_bytes()).hexdigest(),
        "cheapest_bridge": str(CHEAP.relative_to(ROOT)),
        "cheapest_bridge_sha256": hashlib.sha256(CHEAP.read_bytes()).hexdigest(),
        "s3_bridge": str(S3BRIDGE.relative_to(ROOT)),
        "s3_bridge_sha256": hashlib.sha256(S3BRIDGE.read_bytes()).hexdigest(),
    },
    "next": (
        "Run the direct II*_E8_1 branch-zero Abel-Jacobi trace on this degree-46 "
        "multisection using L(47O), then add the already-explicit canonical "
        "D13 points 3*G1-G3. No q8 S3 Hensel lift is required."
    ),
}

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24W_RESULT|"
    f"physical_standard_mw={','.join(map(str,n_std))}|"
    "equation_standard_mw=0,-1,1|"
    f"formula={pretty_formula}|degree=46|"
    "status=PASS_EXACT_Q24_CHEAPEST_Q6_BRIDGE_DEGREE46",
    flush=True,
)
