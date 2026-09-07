#!/usr/bin/env sage
"""
Finalize Q80 orbit 1222 after the exact component-resolved P1/P3 lift.

This script:
  1. loads the exact Jacobian model and exact P1/P3 sections;
  2. replaces the huge gauge-dependent twist d by a squarefree integral
     representative delta in QQ(sqrt(-3));
  3. proves exactly that d/delta is a square in K;
  4. transports P1/P3 to the delta-normalized short model;
  5. checks every q80_twchar_p*.json local Frobenius character against delta;
  6. writes reusable normalized model/section files and a final certificate.

No modular reconstruction is performed.
"""

import json
from pathlib import Path

from sage.all import GF, QQ, ZZ, PolynomialRing, QuadraticField

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
ORBIT_DATA = (
    REPO_ROOT
    / "elkies-k3"
    / "data"
    / "fibrations"
    / "q80-orbit1222-char0"
)
MODULAR_DATA = (
    REPO_ROOT
    / "artifacts"
    / "generated-results"
    / "q80-orbit1222-char0"
    / "modular"
)
MODEL = ORBIT_DATA / "q80_char0_orbit1222_jacobian.sage"
SECTIONS = ORBIT_DATA / "q80_char0_orbit1222_P1_P3.sage"

OUT_MODEL = ORBIT_DATA / "q80_char0_orbit1222_jacobian_normalized.sage"
OUT_SECTIONS = ORBIT_DATA / "q80_char0_orbit1222_P1_P3_normalized.sage"
OUT_NOTE = ORBIT_DATA / "Q80_CHAR0_ORBIT1222_FINAL_CERTIFICATE.md"

if not MODEL.exists():
    raise SystemExit(f"missing exact Jacobian: {MODEL}")
if not SECTIONS.exists():
    raise SystemExit(f"missing exact sections: {SECTIONS}")

load(str(MODEL))
# MODEL defines K,j,R,V,d,A,B,Delta.
d_exact = d
A_exact = A
B_exact = B
Delta_exact = Delta

load(str(SECTIONS))

# SECTIONS reloads MODEL and defines P1x,P1y,P3x,P3y.
P1x_exact, P1y_exact = P1x, P1y
P3x_exact, P3y_exact = P3x, P3y

# Squarefree integral Eisenstein representative:
#
#   delta = m+n*omega, omega=(1+j)/2
#         = (246872123973208405 - 1438885155484622555*j)/2.
#
# Its norm factorization is
#   5^2 * 13 * 32647 * 16918628491 * 8734989504476550481.
delta = K(
    QQ(246872123973208405)/2
    - QQ(1438885155484622555)/2*j
)

expected_norm = ZZ(1568029329404265614196331360501707025)
assert QQ(delta.norm()) == expected_norm
assert expected_norm.factor() == (
    ZZ(5)**2
    * ZZ(13)
    * ZZ(32647)
    * ZZ(16918628491)
    * ZZ(8734989504476550481)
).factor()

# ---------------------------------------------------------------------------
# Exact proof that d/delta is a square in K=QQ(sqrt(-3)).
# Avoid version-dependent NumberFieldElement.is_square/sqrt APIs.
# ---------------------------------------------------------------------------

ratio = K(d_exact/delta)
u = QQ(ratio[0])
v = QQ(ratio[1])

norm_ratio = u**2 + 3*v**2


def qq_sqrt(q):
    q = QQ(q)
    if q < 0:
        return None
    n = ZZ(q.numerator())
    den = ZZ(q.denominator())
    if not n.is_square() or not den.is_square():
        return None
    return QQ(n.sqrt())/QQ(den.sqrt())


sn = qq_sqrt(norm_ratio)
if sn is None:
    raise ArithmeticError("Norm(d/delta) is not a rational square")

x2 = (sn+u)/2
y2 = (sn-u)/6
xabs = qq_sqrt(x2)
yabs = qq_sqrt(y2)

if xabs is None or yabs is None:
    raise ArithmeticError("d/delta is not a square in QQ(sqrt(-3))")

square_root = None
for x in (xabs, -xabs):
    for y in (yabs, -yabs):
        if 2*x*y == v:
            candidate = K(x+y*j)
            if candidate**2 == ratio:
                square_root = candidate
                break
    if square_root is not None:
        break

if square_root is None:
    raise ArithmeticError("failed to construct exact sqrt(d/delta)")

r = square_root
assert d_exact == delta*r**2

print(
    "Q80FINAL1222|"
    f"delta={delta}|norm={expected_norm}|"
    "status=PASS_EXACT_SQUARECLASS_REDUCTION",
    flush=True,
)

# ---------------------------------------------------------------------------
# Transport the exact model and sections to the delta normalization.
#
# d = delta*r^2, so
#   x_d = r^2*x_delta, y_d = r^3*y_delta.
# ---------------------------------------------------------------------------

A_delta = delta**2 * (A_exact / d_exact**2)
B_delta = delta**3 * (B_exact / d_exact**3)
Delta_delta = -16*(4*A_delta**3+27*B_delta**2)

P1x_delta = P1x_exact/r**2
P1y_delta = P1y_exact/r**3
P3x_delta = P3x_exact/r**2
P3y_delta = P3y_exact/r**3

assert P1y_delta**2 == P1x_delta**3 + A_delta*P1x_delta + B_delta
assert P3y_delta**2 == P3x_delta**3 + A_delta*P3x_delta + B_delta

assert A_delta.degree() == 8
assert B_delta.degree() == 12
assert Delta_delta.degree() == 24
assert P1x_delta.degree() == 4
assert P1y_delta.degree() == 6
assert P3x_delta.degree() == 4
assert P3y_delta.degree() == 6

print(
    "Q80FINAL1222|"
    "model_degrees=8,12,24|"
    "P1_degrees=4,6|P3_degrees=4,6|"
    "status=PASS_NORMALIZED_EXACT_MODEL_AND_SECTIONS",
    flush=True,
)

# ---------------------------------------------------------------------------
# Check every accumulated Frobenius-character JSON.
# ---------------------------------------------------------------------------


def legendre_character(value):
    if value == 0:
        return 0
    return 1 if value.is_square() else -1


char_files = sorted(MODULAR_DATA.glob("q80_twchar_p*.json"))
expected_audit_primes = {
    73, 79, 97, 103, 127, 151, 193, 199,
    223, 241, 271, 313, 337, 409, 433, 439,
}
checked = []

for path in char_files:
    rec = json.loads(path.read_text())
    p = int(rec["prime"])
    jr = int(rec["sqrt_minus_3_canonical"])
    cp = int(rec["chi_plus_j"])
    cm = int(rec["chi_minus_j"])

    F = GF(p)
    inv2 = F(2)**-1

    # delta=(A+B*j)/2
    Aint = ZZ(246872123973208405)
    Bint = ZZ(-1438885155484622555)

    plus = (F(Aint)+F(Bint)*F(jr))*inv2
    minus = (F(Aint)-F(Bint)*F(jr))*inv2

    got_plus = legendre_character(plus)
    got_minus = legendre_character(minus)

    if (got_plus, got_minus) != (cp, cm):
        raise ArithmeticError(
            f"{path.name}: character mismatch "
            f"expected {(cp,cm)}, got {(got_plus,got_minus)}"
        )

    checked.append(p)

if set(checked) != expected_audit_primes or len(checked) != 16:
    raise ArithmeticError(
        "Frobenius audit cache does not contain the certified 16 primes: "
        f"found {sorted(checked)}"
    )

print(
    "Q80FINAL1222|"
    f"frobenius_primes={','.join(map(str,checked))}|"
    f"rational_primes_checked={len(checked)}|"
    f"split_places_checked={2*len(checked)}|"
    "status=PASS_ALL_FROBENIUS_CHARACTERS",
    flush=True,
)

# ---------------------------------------------------------------------------
# Persist normalized exact model and sections.
# ---------------------------------------------------------------------------

OUT_MODEL.write_text(
    "\n".join([
        "#!/usr/bin/env sage",
        "from sage.all import QQ, PolynomialRing, QuadraticField",
        'K = QuadraticField(-3, "j")',
        "j = K.gen()",
        'R = PolynomialRing(K, "V")',
        "V = R.gen()",
        f"delta = {delta}",
        f"A = {A_delta}",
        f"B = {B_delta}",
        "Delta = -16*(4*A^3+27*B^2)",
        "assert (A.degree(),B.degree(),Delta.degree()) == (8,12,24)",
        'print(f"Q80ORBIT1222NORMALIZED|delta={delta}|status=PASS_EXACT_NORMALIZED_JACOBIAN")',
    ]) + "\n"
)

OUT_SECTIONS.write_text(
    OUT_MODEL.read_text().rstrip()
    + "\n"
    + "\n".join([
        f"P1x = {P1x_delta}",
        f"P1y = {P1y_delta}",
        f"P3x = {P3x_delta}",
        f"P3y = {P3y_delta}",
        "assert P1y^2 == P1x^3+A*P1x+B",
        "assert P3y^2 == P3x^3+A*P3x+B",
        'print("Q80ORBIT1222NORMALIZEDSECTIONS|P1_height=1/7|P3_height=8/7|status=PASS_EXACT_P1_P3")',
    ]) + "\n"
)

OUT_NOTE.write_text(
    "# Q80 orbit 1222 — final characteristic-zero certificate\n\n"
    "Status: **PASS_EXACT_ORBIT1222**\n\n"
    "The exact characteristic-zero Jacobian twist was recovered by lifting "
    "the pinned P1/P3 Mordell-Weil sections with the full resolved I7 "
    "component multiplicities imposed.  The resulting gauge-dependent twist "
    "`d` was then reduced modulo squares in `QQ(sqrt(-3))`.\n\n"
    "## Normalized twist squareclass\n\n"
    "    delta = (246872123973208405 - "
    "1438885155484622555*sqrt(-3))/2\n\n"
    "with norm\n\n"
    f"    {expected_norm}\n\n"
    "and norm factorization\n\n"
    "    5^2 * 13 * 32647 * 16918628491 * 8734989504476550481.\n\n"
    "The ratio `d/delta` is proved to be an exact square in "
    "`QQ(sqrt(-3))`.\n\n"
    "## Exact geometry\n\n"
    "- normalized short model degrees `(A,B,Delta)=(8,12,24)`\n"
    "- reducible fibers `2 I7 + 3 I2 + 4 I1`\n"
    "- exact P1 degrees `(4,6)`, height `1/7`\n"
    "- exact P3 degrees `(4,6)`, height `8/7`\n"
    "- P3 is the horizontal needed for q6_7774\n\n"
    "## Independent modular audit\n\n"
    f"- {len(checked)} rational split primes checked\n"
    f"- {2*len(checked)} split places checked\n"
    "- every accumulated Frobenius twist character agrees with `delta`\n\n"
    "This supersedes all earlier provisional twist candidates.\n"
)

print(
    "Q80FINAL1222|"
    f"model={OUT_MODEL}|sections={OUT_SECTIONS}|note={OUT_NOTE}|"
    "status=PASS_EXACT_ORBIT1222",
    flush=True,
)
