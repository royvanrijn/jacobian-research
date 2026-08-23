#!/usr/bin/env sage -python
"""
Rebase the q6 child and its p-adic third section to the intrinsic II*/IV* gauge.

Let L_II(T), L_IV(T) be the two linear reducible-fibre factors.  Use

    s = L_II(T) / L_IV(T),

so II* goes to s=0 and IV* to s=infinity.

If
    T = (a*s+b)/(c*s+d),
then for a K3 short Weierstrass model
    y^2 = x^3 + A(T)x + B(T)
the homogeneous rebase is
    x_new = (c*s+d)^4 x_old,
    y_new = (c*s+d)^6 y_old,

and therefore
    A_new = (c*s+d)^8 A(T(s)),
    B_new = (c*s+d)^12 B(T(s)).

The section profile transforms as binary forms with NO extra denominators:
    Z_21 -> degree-21 binary transform,
    X_46 -> degree-46 binary transform,
    Y_69 -> degree-69 binary transform.
After making Z monic, X and Y are rescaled by the corresponding square/cube.

This script uses the existing p-adic Hensel residues, so it requires no new
Newton solve.  It reports rational-reconstruction stability in the new gauge
at 512 and 1024 digits and exactifies immediately if possible.

Run:
  sage -python ~/Downloads/rebase_h92_q6_third_iistar_ivstar.sage \
    artifacts/local/elkies-k3/q6-third-hensel-p1024.json
"""

import argparse
import json
from math import gcd
from pathlib import Path

from sage.all import QQ, ZZ, PolynomialRing, Zp


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
parser.add_argument("input", type=Path)
parser.add_argument("--repo", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
INPUT = args.input
if not INPUT.is_absolute():
    INPUT = ROOT / INPUT
INPUT = INPUT.resolve()
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
OUTPUT = (
    args.output.resolve()
    if args.output and args.output.is_absolute()
    else ROOT / (
        args.output
        if args.output
        else Path("artifacts/local/elkies-k3/q6-third-iistar-ivstar-gauge.json")
    )
)

hensel = json.loads(INPUT.read_text())
child = json.loads(CHILD.read_text())
assert hensel["schema"] == "elkies-k3.h92-q6-third-hensel-lift.v1"
assert hensel["status"] == "PASS_Q6_THIRD_HENSEL"
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"

p = ZZ(hensel["prime"])
precision = int(hensel["precision"])
if precision < 1024:
    raise ValueError("use the p1024 artifact for the gauge comparison")
M = p**precision

QT = PolynomialRing(QQ, "T")
T = QT.gen()
QS = PolynomialRing(QQ, "s")
s = QS.gen()

model = child["minimal_short_weierstrass"]
Aold = QT([QQ(v) for v in model["A_coefficients_low_to_high"]])
Bold = QT([QQ(v) for v in model["B_coefficients_low_to_high"]])

ii = next(item for item in child["finite_fibres"] if item["kodaira"] == "II*")
iv = next(item for item in child["finite_fibres"] if item["kodaira"] == "IV*")
Lii = QT(ii["factor"])
Liv = QT(iv["factor"])
assert Lii.degree() == Liv.degree() == 1
assert tuple(ii["minimal_orders"])[:2] == (4, 5)
assert tuple(iv["minimal_orders"])[:2] == (3, 4)

# s = Lii/Liv.  Solve for T:
# (Liv[1]*s - Lii[1]) T = Lii[0] - Liv[0]*s.
a = -QQ(Liv[0])
b = QQ(Lii[0])
c = QQ(Liv[1])
d = -QQ(Lii[1])
det = a*d - b*c
if not det:
    raise ArithmeticError("II*/IV* fibre factors do not define a PGL2 rebase")

print(
    "Q6REBASE_MAP|"
    f"LII={Lii}|LIV={Liv}|"
    f"T=({a}*s+({b}))/({c}*s+({d}))|status=PASS",
    flush=True,
)


def binary_transform_QQ(poly, degree):
    result = QS(0)
    for i in range(degree + 1):
        coeff = QQ(poly[i]) if i <= poly.degree() else QQ(0)
        if coeff:
            result += coeff*(a*s+b)**i*(c*s+d)**(degree-i)
    return QS(result)


Anew = binary_transform_QQ(Aold, 8)
Bnew = binary_transform_QQ(Bold, 12)

# Structural gauge checks: II* at zero, IV* at infinity.
assert Anew.valuation() >= 4
assert Bnew.valuation() >= 5
assert Anew.degree() <= 5
assert Bnew.degree() <= 8

Acore = Anew // s**4
Bcore = Bnew // s**5
print(
    "Q6REBASE_MODEL|"
    f"A_degree={Anew.degree()}|A_valuation={Anew.valuation()}|"
    f"Acore_degree={Acore.degree()}|"
    f"B_degree={Bnew.degree()}|B_valuation={Bnew.valuation()}|"
    f"Bcore_degree={Bcore.degree()}|status=PASS",
    flush=True,
)


def digits_int(value):
    value = abs(ZZ(value))
    return 1 if not value else len(str(value))


def rational_height_summary(poly):
    nums = [digits_int(QQ(v).numerator()) for v in poly.list() if v]
    dens = [digits_int(QQ(v).denominator()) for v in poly.list() if v]
    return (
        max(nums) if nums else 1,
        max(dens) if dens else 1,
    )


oldA_h = rational_height_summary(Aold)
oldB_h = rational_height_summary(Bold)
newA_h = rational_height_summary(Anew)
newB_h = rational_height_summary(Bnew)
print(
    "Q6REBASE_HEIGHTS|"
    f"oldA_numden={oldA_h[0]},{oldA_h[1]}|"
    f"oldB_numden={oldB_h[0]},{oldB_h[1]}|"
    f"newA_numden={newA_h[0]},{newA_h[1]}|"
    f"newB_numden={newB_h[0]},{newB_h[1]}",
    flush=True,
)

residues = [ZZ(v) for v in hensel["residues"]]
assert len(residues) == 138
Zold_int = residues[:21] + [ZZ(1)]
Xold_int = residues[21:68]
Yold_int = residues[68:]
assert len(Zold_int) == 22 and len(Xold_int) == 47 and len(Yold_int) == 70


def transform_at_precision(k):
    padic = Zp(p, prec=k)
    PR = PolynomialRing(padic, "s")
    sp = PR.gen()
    ap, bp, cp, dp = map(padic, (a, b, c, d))

    def transform(values, degree):
        result = PR(0)
        for i in range(degree + 1):
            coeff = padic(ZZ(values[i]))
            if coeff:
                result += coeff*(ap*sp+bp)**i*(cp*sp+dp)**(degree-i)
        return result

    zraw = transform(Zold_int, 21)
    xraw = transform(Xold_int, 46)
    yraw = transform(Yold_int, 69)

    lc = zraw[21]
    if not lc or lc.valuation() != 0:
        raise ArithmeticError(
            "rebased Z leading coefficient is not a p-adic unit; "
            "the section has a special pole at IV*"
        )

    z = zraw/lc
    x = xraw/lc**2
    y = yraw/lc**3
    assert z[21] == 1

    modulus = p**k

    def lifts(poly, degree):
        return [ZZ(poly[i].lift()) % modulus for i in range(degree + 1)]

    return lifts(z, 20), lifts(x, 46), lifts(y, 69)


def rr(value, k):
    modulus = p**k
    try:
        return ZZ(value % modulus).rational_reconstruction(modulus)
    except (ArithmeticError, ValueError):
        return None


data_by_level = {}
for k in (512, 1024):
    zvals, xvals, yvals = transform_at_precision(k)
    vals = zvals + xvals + yvals
    rec = [rr(v, k) for v in vals]
    data_by_level[k] = (vals, rec)
    print(
        "Q6REBASE_RR|"
        f"precision={k}|recovered={sum(v is not None for v in rec)}/138|"
        f"Z={sum(v is not None for v in rec[:21])}/21|"
        f"X={sum(v is not None for v in rec[21:68])}/47|"
        f"Y={sum(v is not None for v in rec[68:])}/70",
        flush=True,
    )

r512 = data_by_level[512][1]
r1024 = data_by_level[1024][1]
stable = [
    i for i, (u, v) in enumerate(zip(r512, r1024))
    if u is not None and u == v
]
print(
    "Q6REBASE_STABLE|512->1024|"
    f"stable={len(stable)}/138|"
    f"Z={sum(i < 21 for i in stable)}/21|"
    f"X={sum(21 <= i < 68 for i in stable)}/47|"
    f"Y={sum(i >= 68 for i in stable)}/70",
    flush=True,
)

complete = all(v is not None for v in r1024)
exact = False
if complete:
    Z = QS(list(r1024[:21]) + [QQ(1)])
    X = QS(list(r1024[21:68]))
    Y = QS(list(r1024[68:]))
    exact = (
        Z.degree() == 21
        and X.degree() == 46
        and Y.degree() == 69
        and Y**2 == X**3 + Anew*X*Z**4 + Bnew*Z**6
    )
    print(
        f"Q6REBASE_EXACT|complete=1|identity={int(exact)}",
        flush=True,
    )

payload = {
    "schema": "elkies-k3.h92-q6-third-iistar-ivstar-gauge.v1",
    "status": (
        "PASS_EXACT_Q6_THIRD_REBASED"
        if exact
        else "PASS_Q6_THIRD_REBASE_DIAGNOSTIC"
    ),
    "source_hensel": str(INPUT.relative_to(ROOT)),
    "prime": int(p),
    "precision": precision,
    "base_change": {
        "definition": "s=L_II(T)/L_IV(T)",
        "L_II": str(Lii),
        "L_IV": str(Liv),
        "T_numerator": [str(b), str(a)],
        "T_denominator": [str(d), str(c)],
    },
    "model": {
        "A": [str(v) for v in Anew.list()],
        "B": [str(v) for v in Bnew.list()],
        "Acore": [str(v) for v in Acore.list()],
        "Bcore": [str(v) for v in Bcore.list()],
    },
    "rr_1024": {
        "recovered": int(sum(v is not None for v in r1024)),
        "Z": int(sum(v is not None for v in r1024[:21])),
        "X": int(sum(v is not None for v in r1024[21:68])),
        "Y": int(sum(v is not None for v in r1024[68:])),
    },
    "stable_512_1024": {
        "total": len(stable),
        "Z": int(sum(i < 21 for i in stable)),
        "X": int(sum(21 <= i < 68 for i in stable)),
        "Y": int(sum(i >= 68 for i in stable)),
    },
}

if exact:
    payload["Z"] = [str(v) for v in Z.list()]
    payload["X"] = [str(v) for v in X.list()]
    payload["Y"] = [str(v) for v in Y.list()]
    payload["exact_identity"] = True

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
print(f"Q6REBASE_RESULT|status={payload['status']}", flush=True)
