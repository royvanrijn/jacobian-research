#!/usr/bin/env sage -python
"""Recover the R17-directed q24/orbit42 section by exact MW halving.

The exact polynomial D12 model uses the R3 zero.  Its eighteen rational
identity-class zero-pole sections form the complete height-four identity
shell.  Match that shell to the exact R3-zero MW lattice, form twice the
orbit42 target, and halve it in the function-field elliptic curve.

The modular shell matching is only a permutation/seed aid.  Any emitted QQ
point must be replayed with exact rational group law and exact duplication.
"""

import argparse
import json
from pathlib import Path

import networkx as nx
from networkx.algorithms import isomorphism as nxiso

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
ZERO = LOCAL / "q24-orbit42-rational-zero-pole-sections-qq.json"
MODEL = LOCAL / "q24-orbit42-zero-pole-seeds-mod-43.json"
AUDIT = LOCAL / "q24-orbit42-identity-halving-audit.json"
SPINOR = LOCAL / "q24-orbit42-spinor-zero-profiles.json"
POINTED = LOCAL / "q24-downstream-lift/pointed-d12-a11-profile-p100003.json"
EASY = LOCAL / "q24-downstream-lift/a11-easy-spinor-shift-p100003.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=100003)
parser.add_argument("--matching-only", action="store_true")
parser.add_argument("--output", type=Path)
args = parser.parse_args()
p = ZZ(args.prime)
F = GF(p)

for path in (ZERO, MODEL, AUDIT, SPINOR, POINTED, EASY):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

zero = json.loads(ZERO.read_text())
model = json.loads(MODEL.read_text())
audit = json.loads(AUDIT.read_text())
spinor = json.loads(SPINOR.read_text())
pointed = json.loads(POINTED.read_text())
easy = json.loads(EASY.read_text())

assert zero["status"] == "PASS_EXACT_Q42_RATIONAL_ZERO_POLE_SECTIONS_QQ"
assert model["status"] == "PASS_Q42_ZERO_POLE_SMALLPRIME_SEEDS"
assert audit["status"] == "PASS_Q42_ORBIT42_IDENTITY_HALVING_LATTICE_GATE"
assert spinor["status"] == "PASS_Q24_ORBIT42_EXACT_SPINOR_ZERO_PROFILES"
assert pointed["status"] == "PASS_Q24_POINTED_D12_A11_PROFILE"

exact_model = model["exact_model"]
RQ = PolynomialRing(QQ, "u")
u = RQ.gen()
KQ = RQ.fraction_field()
Aq = RQ([QQ(v) for v in exact_model["A_coefficients_low_to_high"]])
Bq = RQ([QQ(v) for v in exact_model["B_coefficients_low_to_high"]])
EQ = EllipticCurve(KQ, [0, 0, 0, KQ(Aq), KQ(Bq)])

RF = PolynomialRing(F, "u")
uf = RF.gen()
KF = RF.fraction_field()


def red_q(value):
    value = QQ(value)
    den = ZZ(value.denominator())
    if den % p == 0:
        raise ZeroDivisionError(f"denominator divisible by matching prime {p}")
    return F(ZZ(value.numerator())) / F(den)


def red_poly(poly):
    return RF([red_q(v) for v in RQ(poly).list()])


Af = red_poly(Aq)
Bf = red_poly(Bq)
EF = EllipticCurve(KF, [0, 0, 0, KF(Af), KF(Bf)])

equation_points_q = []
equation_points_f = []
for row in zero["sections"]:
    xq = RQ([QQ(v) for v in row["x_coefficients_low_to_high"]])
    yq = RQ([QQ(v) for v in row["y_coefficients_low_to_high"]])
    if yq**2 != xq**3 + Aq * xq + Bq:
        raise ArithmeticError("exact zero-pole point misses exact D12 model")
    equation_points_q.append(EQ(KQ(xq), KQ(yq)))
    equation_points_f.append(EF(KF(red_poly(xq)), KF(red_poly(yq))))

if len(equation_points_q) != 18:
    raise ArithmeticError("expected eighteen exact equation points")

# The zero-pole lift uses a different exact base/gauge normalization from the
# pointed p=100003 profiler.  Recover the affine base/gauge transport below;
# it will anchor the MW-shell isometry by the known A0-R3 point.

# -------------------------------------------------------------------------
# 1. Abstract R3-zero height-four identity shell.
# -------------------------------------------------------------------------
r3_profile = next(row for row in spinor["profiles"] if row["zero"] == "R3")
G = matrix(ZZ, r3_profile["frame"])
root = G[:12, :12]
coupling = G[:12, 12:]
tail = G[12:, 12:]
H = tail - coupling.transpose() * root.inverse() * coupling

r3_audit = audit["exact_model_R3_zero"]
abstract = [vector(ZZ, row) for row in r3_audit["identity_vectors"]]
target = vector(ZZ, r3_audit["target_mw"])
double_target = vector(ZZ, r3_audit["double_target_mw"])
if len(abstract) != 18 or QQ(target * H * target) != 19:
    raise ArithmeticError("R3 identity/target lattice gate is inconsistent")

coeff_by_vector = {
    tuple(row["mw"]): ZZ(row["coefficient"])
    for row in r3_audit["terms"]
}
if sum(
    (coeff_by_vector.get(tuple(z), ZZ.zero()) * z for z in abstract),
    vector(ZZ, [0] * 5),
) != double_target:
    raise ArithmeticError("stored R3 doubling expression is inconsistent")

# -------------------------------------------------------------------------
# 2. Pairing-colored complete graph over GF(p).
#
# For distinct identity-component P.O=0 sections on a K3,
#     <P,Q> = 2 - P.Q.
# The finite intersection multiplicity is deg gcd(xP-xQ,yP-yQ).  The complete
# colored graph is matched to the exact MW Gram graph.  The transported
# A0-R3 point anchors both its vector and its sign.
# -------------------------------------------------------------------------
abstract_pair = matrix(ZZ, 18, 18)
equation_pair = matrix(ZZ, 18, 18)
for i in range(18):
    for j in range(18):
        abstract_pair[i, j] = ZZ(abstract[i] * H * abstract[j])
        if i == j:
            equation_pair[i, j] = 4
            continue
        xi, yi = equation_points_f[i].xy()
        xj, yj = equation_points_f[j].xy()
        xd = RF(xi - xj)
        yd = RF(yi - yj)
        intersection = xd.gcd(yd).degree()
        equation_pair[i, j] = ZZ(2 - intersection)

if sorted(abstract_pair.list()) != sorted(equation_pair.list()):
    raise ArithmeticError("equation shell pairing multiset misses abstract shell")

# Known Q=A0-R3 in the pointed R3-zero group law.
q_mw = vector(ZZ, pointed["R3_zero_lattice_marking"]["explicit_A0_minus_R3_mw"])
if tuple(q_mw) not in {tuple(z) for z in abstract}:
    raise ArithmeticError("known A0-R3 MW vector is absent from identity shell")
abstract_anchor = next(i for i, z in enumerate(abstract) if z == q_mw)

if p != 100003:
    raise ValueError("the pointed-model marking anchor is pinned at p=100003")

pointed_model = easy["I8star_infinity"]
Apoint = RF(pointed_model["A"]["coefficients_low_to_high"])
Bpoint = RF(pointed_model["B"]["coefficients_low_to_high"])

Siso = PolynomialRing(F, names=("s", "c", "k4", "k6"))
siso, ciso, k4iso, k6iso = Siso.gens()
Uiso = PolynomialRing(Siso, "U")
U0 = Uiso.gen()
transport_ideal = Siso.ideal(
    list(Uiso(Af) - k4iso * Uiso(Apoint)(siso * U0 + ciso))
    + list(Uiso(Bf) - k6iso * Uiso(Bpoint)(siso * U0 + ciso))
    + [k6iso**2 - k4iso**3]
)
transport_solutions = transport_ideal.variety()
if len(transport_solutions) != 1:
    raise ArithmeticError(
        f"expected one pointed-to-zero-model transport, got {len(transport_solutions)}"
    )
transport = transport_solutions[0]
base_scale = F(transport[siso])
base_shift = F(transport[ciso])
k4 = F(transport[k4iso])
k6 = F(transport[k6iso])
x_scale = k6 / k4
if x_scale**2 != k4 or x_scale**3 != k6:
    raise ArithmeticError("pointed-to-zero-model x gauge is inconsistent")

qxp = RF(pointed_model["known_Q_x"]["coefficients_low_to_high"])
qyp = RF(pointed_model["known_Q_y"]["coefficients_low_to_high"])
qx = RF(x_scale * qxp(base_scale * uf + base_shift))

y_root = k6.sqrt() if k6.is_square() else None
if y_root is None:
    raise ArithmeticError("pointed-to-zero-model y gauge has no square root")
y_scales = [y_root, -y_root]
anchor_hits = []
for y_scale in y_scales:
    qy = RF(y_scale * qyp(base_scale * uf + base_shift))
    for i, point in enumerate(equation_points_f):
        if RF(point.xy()[0]) == qx and RF(point.xy()[1]) == qy:
            anchor_hits.append((i, y_scale))
if len(anchor_hits) != 2:
    raise ArithmeticError(
        f"transported A0-R3 point matches {len(anchor_hits)} exact shell points"
    )
equation_anchor_set = {index for index, unused_scale in anchor_hits}
if len(equation_anchor_set) != 2:
    raise ArithmeticError("the two y gauges do not give a signed anchor pair")
abstract_anchor_set = {
    i for i, z in enumerate(abstract) if z == q_mw or z == -q_mw
}
if len(abstract_anchor_set) != 2:
    raise ArithmeticError("abstract A0-R3 signed anchor pair is incomplete")

print(
    "Q42HALF_NORMALIZATION|"
    f"base_scale={int(base_scale)}|base_shift={int(base_shift)}|"
    f"x_scale={int(x_scale)}|y_scale=SIGNED_PAIR|"
    f"anchor_pair={','.join(map(str, sorted(abstract_anchor_set)))}->"
    f"{','.join(map(str, sorted(equation_anchor_set)))}|"
    "status=PASS_POINTED_ANCHOR_PAIR",
    flush=True,
)

GA = nx.Graph()
GE = nx.Graph()
for i in range(18):
    GA.add_node(i, anchor=(i in abstract_anchor_set))
    GE.add_node(i, anchor=(i in equation_anchor_set))
for i in range(18):
    for j in range(i + 1, 18):
        GA.add_edge(i, j, pairing=int(abstract_pair[i, j]))
        GE.add_edge(i, j, pairing=int(equation_pair[i, j]))

matcher = nxiso.GraphMatcher(
    GA,
    GE,
    node_match=nxiso.categorical_node_match("anchor", False),
    edge_match=nxiso.categorical_edge_match("pairing", None),
)

mappings = []
mapping_keys = set()
for mapping in matcher.isomorphisms_iter():
    key = tuple(mapping[i] for i in range(18))
    if key in mapping_keys:
        continue
    mapping_keys.add(key)
    mappings.append(key)

if not mappings:
    raise ArithmeticError("no anchor-preserving identity-shell marking")

print(
    "Q42HALF_MATCH|"
    f"prime={p}|abstract=18|equation=18|anchor=SIGNED_PAIR|"
    f"mappings={len(mappings)}|status=PASS_MODULAR_SHELL_MARKING",
    flush=True,
)

matching_payload = {
    "prime": int(p),
    "abstract_anchor_index": abstract_anchor,
    "equation_anchor_indices": sorted(equation_anchor_set),
    "known_anchor_mw": list(map(int, q_mw)),
    "pointed_to_zero_model": {
        "base_scale": int(base_scale),
        "base_shift": int(base_shift),
        "x_scale": int(x_scale),
        "y_scales": sorted(int(value) for value in y_scales),
    },
    "mapping_count": len(mappings),
    "mappings_abstract_to_equation": [list(map(int, row)) for row in mappings],
    "abstract_pairing": [list(map(int, row)) for row in abstract_pair.rows()],
    "equation_pairing_mod_p": [list(map(int, row)) for row in equation_pair.rows()],
}

if args.matching_only:
    payload = {
        "schema": "elkies-k3.h3-q24-orbit42-identity-halving-qq.v1",
        "status": "PASS_Q42_IDENTITY_SHELL_MATCHING_ONLY",
        "matching": matching_payload,
        "proof_boundary": (
            "The exact rational points are verified over QQ, and their reduction "
            "is matched to the exact R3-zero MW shell modulo p. The doubled "
            "target has not yet been formed or halved."
        ),
    }
    output = (
        args.output.resolve()
        if args.output
        else LOCAL / "q24-orbit42-identity-halving-qq.json"
    )
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"OUTPUT|{output}", flush=True)
    print(
        "Q42HALF_QQ_RESULT|stage=MATCHING|"
        f"mappings={len(mappings)}|status=PASS_Q42_IDENTITY_SHELL_MATCHING_ONLY",
        flush=True,
    )
    raise SystemExit(0)

# -------------------------------------------------------------------------
# 3. Form every marking-compatible doubled target modulo p and halve it.
# -------------------------------------------------------------------------
def point_key(point):
    if point.is_zero():
        return ("ZERO",)
    xx, yy = point.xy()
    return (str(xx), str(yy))


double_candidates = {}
for mapping_index, mapping in enumerate(mappings):
    point = EF(0)
    used = []
    for abstract_index, z in enumerate(abstract):
        coefficient = coeff_by_vector.get(tuple(z), ZZ.zero())
        if not coefficient:
            continue
        equation_index = mapping[abstract_index]
        point += int(coefficient) * equation_points_f[equation_index]
        used.append((abstract_index, equation_index, int(coefficient)))
    if point.is_zero():
        raise ArithmeticError("doubled orbit42 candidate collapsed to zero")
    double_candidates.setdefault(point_key(point), {
        "point": point,
        "mapping_indices": [],
        "used": used,
    })["mapping_indices"].append(mapping_index)

XR = PolynomialRing(KF, "Xh")
Xh = XR.gen()


def power_root(poly, exponent):
    poly = RF(poly)
    if not poly or poly.degree() <= 0:
        return RF.one()
    lc = poly.leading_coefficient()
    monic = poly / lc
    answer = RF.one()
    for factor, multiplicity in monic.factor():
        if int(multiplicity) % exponent:
            return None
        answer *= factor.monic()**(int(multiplicity) // exponent)
    return answer.monic()


def pole_degree(point):
    if point.is_zero():
        return -1
    xx, yy = point.xy()
    zx = power_root(RF(xx.denominator()), 2)
    zy = power_root(RF(yy.denominator()), 3)
    if zx is None or zy is None or zx != zy:
        raise ArithmeticError("halved point has incompatible section denominators")
    return int(zx.degree())


halves = {}
for double_index, row in enumerate(double_candidates.values()):
    Q2 = row["point"]
    qx, qy = Q2.xy()
    quartic = (
        Xh**4
        - 4 * qx * Xh**3
        - 2 * KF(Af) * Xh**2
        - (4 * KF(Af) * qx + 8 * KF(Bf)) * Xh
        + KF(Af)**2
        - 4 * KF(Bf) * qx
    )
    linear = []
    half_factorization = quartic.factor()
    half_profile = []
    for factor, multiplicity in half_factorization:
        half_profile.append([int(factor.degree()), int(multiplicity)])
        if factor.degree() == 1:
            linear.extend([factor] * int(multiplicity))

    print(
        "Q42HALF_MOD_DOUBLE|"
        f"index={double_index}|mappings={len(row['mapping_indices'])}|"
        f"factors={half_profile}|linear_halves={len(linear)}|"
        "status=PASS_DUPLICATION_FACTOR",
        flush=True,
    )

    for factor in linear:
        xp = KF(-factor[0] / factor[1])
        numerator = (
            (3 * xp**2 + KF(Af)) * (xp - qx)
            - 2 * (xp**3 + KF(Af) * xp + KF(Bf))
        )
        if not qy:
            continue
        yp = numerator / (2 * qy)
        try:
            P = EF(xp, yp)
        except (ArithmeticError, TypeError, ValueError):
            continue
        if 2 * P != Q2:
            continue
        pd = pole_degree(P)
        key = point_key(P)
        halves.setdefault(key, {
            "point": P,
            "pole_degree": pd,
            "double_indices": [],
            "mapping_indices": [],
        })["double_indices"].append(double_index)
        halves[key]["mapping_indices"].extend(row["mapping_indices"])

degree9 = [row for row in halves.values() if row["pole_degree"] == 9]
relative_sections = {}
for row in degree9:
    for mapping_index in sorted(set(row["mapping_indices"])):
        mapping = mappings[mapping_index]
        Qanchor = equation_points_f[mapping[abstract_anchor]]
        # D42 has old-fibre degree two.  Changing the zero from R3 to A0=Q
        # changes its Abel--Jacobi point by D-2Q, hence subtract TWO copies
        # of Q.  Subtracting one copy is the degree-one section formula and
        # gives the wrong height nine class.
        relative = row["point"] - 2 * Qanchor
        relative_degree = pole_degree(relative)
        key = point_key(relative)
        relative_sections.setdefault(key, {
            "point": relative,
            "pole_degree": relative_degree,
            "mapping_indices": [],
        })["mapping_indices"].append(mapping_index)

degree3 = [
    row for row in relative_sections.values()
    if row["pole_degree"] == 3
]

# Compile the surviving relative sections by the same direct chord/Yun method
# as compile_h92_q24_orbit42_a11_chord_modp.sage.
Rm = PolynomialRing(F, "m")
m0 = Rm.gen()
Km = Rm.fraction_field()
RUm = PolynomialRing(Km, "u")
m = Km(m0)


def projectivize(point):
    xx, yy = point.xy()
    nx, dx = RF(xx.numerator()), RF(xx.denominator())
    ny, dy = RF(yy.numerator()), RF(yy.denominator())
    lcx, lcy = dx.leading_coefficient(), dy.leading_coefficient()
    nx, dx = nx / lcx, dx / lcx
    ny, dy = ny / lcy, dy / lcy
    zx, zy = power_root(dx, 2), power_root(dy, 3)
    if zx is None or zy is None or zx != zy:
        raise ArithmeticError("relative section does not projectivize")
    if ny**2 != nx**3 + Af * nx * zx**4 + Bf * zx**6:
        raise ArithmeticError("projectivized relative section identity failed")
    return RF(nx), RF(ny), RF(zx)


def lift_poly(poly):
    return RUm([Km(c) for c in RF(poly).list()])


def binary_invariants(quartic):
    coeff = [Km.zero()] * 5
    for i, value in enumerate(RUm(quartic).list()):
        if i < 5:
            coeff[i] = Km(value)
    e, d, c, b, a = coeff
    I = 12*a*e - 3*b*d + c*c
    J = 72*a*c*e + 9*b*c*d - 27*a*d*d - 27*b*b*e - 2*c*c*c
    return I, J


def as_m_poly(value):
    value = Km(value)
    denominator = Rm(value.denominator())
    if denominator.degree() > 0:
        return None
    return Rm(value.numerator()) / F(denominator[0])


chord_results = []
for section_index, row in enumerate(degree3):
    X, Y, Z = projectivize(row["point"])
    if Z.degree() != 3:
        raise ArithmeticError("relative section projective Z degree is not three")
    XR, YR, ZR, AR = map(lift_poly, (X, Y, Z, Af))
    raw = RUm(
        m**4 * ZR**4
        - 6 * XR * m**2 * ZR**2
        - 8 * YR * m * ZR
        - 3 * XR**2
        - 4 * AR * ZR**4
    )
    factorization = raw.factor()
    reduced = RUm(factorization.unit())
    square = RUm.one()
    factor_profile = []
    for factor, multiplicity in factorization:
        multiplicity = int(multiplicity)
        factor_profile.append([int(factor.degree()), multiplicity])
        square *= factor**(multiplicity // 2)
        if multiplicity % 2:
            reduced *= factor
    reduced = RUm(reduced)

    has_a11 = False
    fibre_profile = []
    jacA = jacB = None
    if reduced.degree() in (3, 4):
        I, J = binary_invariants(reduced)
        jacA, jacB = -27 * I, -27 * J
        discriminant = -16 * (4 * jacA**3 + 27 * jacB**2)
        Ap, Bp, Dp = map(as_m_poly, (jacA, jacB, discriminant))
        if Ap is not None and Bp is not None and Dp is not None:
            Ap, Bp, Dp = Rm(Ap), Rm(Bp), Rm(Dp)
            for factor, multiplicity in Dp.factor():
                multiplicity = int(multiplicity)
                va, ta = 0, Ap
                while ta and ta % factor == 0:
                    va += 1
                    ta //= factor
                vb, tb = 0, Bp
                while tb and tb % factor == 0:
                    vb += 1
                    tb //= factor
                fibre_profile.append({
                    "factor": [int(x) for x in factor.list()],
                    "degree": int(factor.degree()),
                    "orders": [va, vb, multiplicity],
                })
                if multiplicity == 12 and va == 0 and vb == 0:
                    has_a11 = True
            infinity_orders = [
                8 - Ap.degree(), 12 - Bp.degree(), 24 - Dp.degree()
            ]
            fibre_profile.append({
                "factor": "infinity",
                "orders": list(map(int, infinity_orders)),
            })
            if infinity_orders == [0, 0, 12]:
                has_a11 = True

    chord_results.append({
        "section_index": section_index,
        "mapping_indices": sorted(set(row["mapping_indices"])),
        "X": [int(v) for v in X.list()],
        "Y": [int(v) for v in Y.list()],
        "Z": [int(v) for v in Z.list()],
        "raw_degree": int(raw.degree()),
        "factor_profile": factor_profile,
        "square_factor_degree": int(square.degree()),
        "quartic_degree": int(reduced.degree()),
        "quartic": [str(v) for v in reduced.list()],
        "jacobian_A": None if jacA is None else str(jacA),
        "jacobian_B": None if jacB is None else str(jacB),
        "fibre_profile": fibre_profile,
        "has_I12_A11": has_a11,
    })
    print(
        "Q42HALF_MOD_CHORD|"
        f"section={section_index}|raw={raw.degree()}|quartic={reduced.degree()}|"
        f"square={square.degree()}|A11={int(has_a11)}|"
        f"status={'PASS_A11_MODP' if has_a11 else 'COMPILED_NON_A11'}",
        flush=True,
    )

a11_chords = [row for row in chord_results if row["has_I12_A11"]]
print(
    "Q42HALF_MOD_RESULT|"
    f"double_candidates={len(double_candidates)}|halves={len(halves)}|"
    f"PdotO9={len(degree9)}|relative={len(relative_sections)}|"
    f"PdotO3={len(degree3)}|A11={len(a11_chords)}|"
    "status=PASS_MODULAR_HALVING_CENSUS",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h3-q24-orbit42-identity-halving-qq.v1",
    "status": (
        "PASS_Q42_MODULAR_IDENTITY_HALVING_A11_CANDIDATES"
        if a11_chords else "Q42_IDENTITY_HALVING_HAS_NO_A11_CHORD"
    ),
    "matching": matching_payload,
    "modular_halving": {
        "double_candidate_count": len(double_candidates),
        "rational_half_count": len(halves),
        "P_dot_O_9_count": len(degree9),
        "relative_section_count": len(relative_sections),
        "relative_P_dot_O_3_count": len(degree3),
        "A11_chord_count": len(a11_chords),
        "candidates": [
            {
                "index": index,
                "pole_degree": row["pole_degree"],
                "double_indices": sorted(set(row["double_indices"])),
                "mapping_indices": sorted(set(row["mapping_indices"])),
                "x": str(row["point"].xy()[0]),
                "y": str(row["point"].xy()[1]),
            }
            for index, row in enumerate(halves.values())
        ],
        "relative_candidates": [
            {
                "index": index,
                "pole_degree": row["pole_degree"],
                "mapping_indices": sorted(set(row["mapping_indices"])),
                "x": str(row["point"].xy()[0]),
                "y": str(row["point"].xy()[1]),
            }
            for index, row in enumerate(relative_sections.values())
        ],
        "chord_results": chord_results,
    },
    "next": (
        "Replay the surviving mapping(s) with exact QQ group law, factor the "
        "exact duplication quartic, and require exact P.O=9 plus the orbit42 "
        "A11 child regression."
    ),
    "proof_boundary": (
        "Exact input points and exact lattice shell, with modular shell marking "
        "and modular duplication/halving only. No characteristic-zero orbit42 "
        "section or A11 child is claimed."
    ),
}
output = (
    args.output.resolve()
    if args.output
    else LOCAL / "q24-orbit42-identity-halving-qq.json"
)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{output}", flush=True)
print(
    "Q42HALF_QQ_RESULT|stage=MODULAR_HALVING|"
    f"PdotO9={len(degree9)}|PdotO3={len(degree3)}|A11={len(a11_chords)}|"
    f"status={payload['status']}",
    flush=True,
)
