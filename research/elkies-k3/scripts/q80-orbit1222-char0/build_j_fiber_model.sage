#!/usr/bin/env sage
"""
Build the exact characteristic-zero short Weierstrass model for Q80 orbit 1222.

Requires the already-certified permanent data file:
    q80_char0_orbit1222_mu_critical_factors.sage

If
    j(V) = 1728*C^3 / (C^3 - mu*S^2),
then choosing
    c4 = mu*C,
    c6 = mu^2*S
gives
    c4^3 - c6^2 = mu^3*(C^3-mu*S^2).

Thus no sqrt(mu) is needed:
    Delta = mu^3*(C^3-mu*S^2)/1728
    A = -c4/48
    B = -c6/864.

The model is entirely over QQ(sqrt(-3)).
"""

import re
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
ORBIT_DATA = (
    REPO_ROOT
    / "elkies-k3"
    / "data"
    / "fibrations"
    / "q80-orbit1222-char0"
)
SOURCE = ORBIT_DATA / "q80_char0_orbit1222_mu_critical_factors.sage"
OUT = ORBIT_DATA / "q80_char0_orbit1222_weierstrass.sage"
NOTE = ORBIT_DATA / "Q80_CHAR0_ORBIT1222_WEIERSTRASS.md"

if not SOURCE.exists():
    raise SystemExit(f"missing exact mu source: {SOURCE}")

load(str(SOURCE))

# SOURCE defines K=QQ(sqrt(-3)), j, R=K[V], V, C, S, mu.
c4 = mu*C
c6 = mu**2*S
Delta = mu**3*(C**3-mu*S**2)/1728

assert c4**3-c6**2 == 1728*Delta

A = -c4/48
B = -c6/864

assert -48*A == c4
assert -864*B == c6
assert -16*(4*A**3+27*B**2) == Delta

assert A.degree() == 8
assert B.degree() == 12
assert Delta.degree() == 24

print(
    "Q80CHAR0MODEL|"
    f"A_degree={A.degree()}|B_degree={B.degree()}|Delta_degree={Delta.degree()}|"
    "field=QQ(sqrt(-3))|stage=EXACT_SHORT_MODEL",
    flush=True,
)

# Extract the fiber-support polynomials WITHOUT full factorization.
# For Delta = P7^7 * P2^2 * P1 (up to scalar):
#   G1 = gcd(Delta,Delta')  = P7^6 * P2
#   G2 = gcd(G1,G1')       = P7^5
#   G3 = gcd(G2,G2')       = P7^4
G1 = Delta.gcd(Delta.derivative()).monic()
G2 = G1.gcd(G1.derivative()).monic()
G3 = G2.gcd(G2.derivative()).monic()

assert (G1.degree(), G2.degree(), G3.degree()) == (15, 10, 8)

P7, r = G2.quo_rem(G3)
assert r == 0
P7 = P7.monic()

P7P2, r = G1.quo_rem(G2)
assert r == 0
P7P2 = P7P2.monic()

P2, r = P7P2.quo_rem(P7)
assert r == 0
P2 = P2.monic()

repeated = P7**7 * P2**2
P1raw, r = Delta.quo_rem(repeated)
assert r == 0
P1 = P1raw.monic()

assert (P7.degree(), P2.degree(), P1.degree()) == (2, 3, 4)
assert Delta.monic() == P7**7 * P2**2 * P1

print(
    "Q80CHAR0MODEL|"
    f"I7_support_degree={P7.degree()}|I2_support_degree={P2.degree()}|"
    f"I1_support_degree={P1.degree()}|"
    "fibers=2I7+3I2+4I1|stage=PASS_EXACT_FIBER_SUPPORT",
    flush=True,
)


def parse_expected_kernel(path):
    text = path.read_text()
    match = re.search(
        r"expected_kernel\s*=\s*vector\(\s*finite\s*,\s*\[(.*?)\]\s*,?\s*\)",
        text,
        re.S,
    )
    if not match:
        raise RuntimeError(f"could not parse expected_kernel from {path}")
    values = [int(x) for x in re.findall(r"-?\d+", match.group(1))]
    if len(values) != 50:
        raise RuntimeError(f"{path}: parsed {len(values)} values, expected 50")
    return values


def reduce_poly(poly, jroot):
    F = GF(73)
    RV = PolynomialRing(F, "V")
    jr = F(jroot)
    return RV([
        F(QQ(K(c)[0])) + jr*F(QQ(K(c)[1]))
        for c in poly.list()
    ])


def validate_kernel(jroot, expected):
    F = GF(73)
    RV = PolynomialRing(F, "V")
    Np = RV([F(x) for x in expected[:25]])
    Dp = RV([F(x) for x in expected[25:]])

    Ne = reduce_poly(c4**3, jroot)
    De = reduce_poly(Delta, jroot)

    # Models may differ by a nonzero global scalar; compare rational j-maps.
    assert Ne*Dp == De*Np
    return Ne.degree(), De.degree()


orig = parse_expected_kernel(
    REPO_ROOT / "elkies-k3/scripts/reconstruct_q80_third_q12_jacobian_gf73.sage"
)
conj = parse_expected_kernel(
    REPO_ROOT / "elkies-k3/scripts/analyze_q80_third_q12_galois_descent_gf73.sage"
)

deg_plus = validate_kernel(17, orig)
deg_minus = validate_kernel(56, conj)

print(
    "Q80CHAR0MODEL|"
    f"p73_plus_degrees={deg_plus[0]},{deg_plus[1]}|"
    f"p73_minus_degrees={deg_minus[0]},{deg_minus[1]}|"
    "status=PASS_P73_BOTH_GALOIS_JMAPS",
    flush=True,
)

# Write a compact exact model file.
lines = [
    "#!/usr/bin/env sage",
    "from sage.all import PolynomialRing, QuadraticField",
    'K = QuadraticField(-3, "j")',
    "j = K.gen()",
    'R = PolynomialRing(K, "V")',
    "V = R.gen()",
    f"A = {A}",
    f"B = {B}",
    "c4 = -48*A",
    "c6 = -864*B",
    "Delta = -16*(4*A^3+27*B^2)",
    f"P7 = {P7}",
    f"P2 = {P2}",
    f"P1 = {P1}",
    "assert Delta.monic() == P7^7*P2^2*P1",
    'print("Q80ORBIT1222WEIERSTRASS|fibers=2I7+3I2+4I1|field=QQ(sqrt(-3))|status=PASS_EXACT_MODEL")',
]
OUT.write_text("\n".join(lines) + "\n")

NOTE.write_text(
    "# Q80 orbit 1222 — exact characteristic-zero Weierstrass model\n\n"
    "Status: **PASS_EXACT_MODEL**\n\n"
    "The exact critical-value computation determines `mu`, while certified "
    "monic polynomials `C` (degree 8) and `S` (degree 12) determine the "
    "j-map.  No square-root extension is needed: choose\n\n"
    "    c4 = mu*C\n"
    "    c6 = mu^2*S\n"
    "    Delta = mu^3*(C^3-mu*S^2)/1728.\n\n"
    "Then `c4^3-c6^2=1728*Delta` identically, so the short model is\n\n"
    "    y^2 = x^3 + A(V)*x + B(V)\n"
    "    A = -c4/48\n"
    "    B = -c6/864.\n\n"
    "The model is defined over `QQ(sqrt(-3))`.\n\n"
    "Exact fiber-support polynomials were extracted from the gcd tower, "
    "without full discriminant factorization:\n\n"
    "    Delta ~ P7^7 * P2^2 * P1\n\n"
    f"with degrees `(deg P7, deg P2, deg P1)=({P7.degree()}, {P2.degree()}, {P1.degree()})`, "
    "hence `2 I7 + 3 I2 + 4 I1` over the algebraic closure.\n\n"
    "Both independent p=73 Galois j-maps validate by exact cross-multiplication.\n"
)

print(
    f"Q80CHAR0MODEL|model={OUT}|note={NOTE}|status=PASS_EXACT_MODEL",
    flush=True,
)
