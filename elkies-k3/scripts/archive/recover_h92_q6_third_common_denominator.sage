#!/usr/bin/env sage -python
"""
Try exact q6-third recovery from the common denominator lattice of X(T).

At p^1024 ordinary balanced reconstruction already recovers:
  Z: 21/21 coefficients
  X: 40/47 coefficients
  Y: 36/70 coefficients.

This script exploits the 40 exact X coefficients.

For several denominator anchors D built from already reconstructed
coefficients, replace each missing X_i by the unique centered integer n_i
satisfying

    n_i == D * residue_i  (mod p^precision),

and set X_i=n_i/D.  If D is a multiple of the true denominator lattice and
the integer numerator fits inside (-M/2,M/2), this is the exact coefficient.

Each candidate is accepted ONLY if

    X^3 + A*X*Z^4 + B*Z^6

is an exact square in QQ[T].  Its square root then recovers all 70 Y
coefficients at once.  Known p-adic/reconstructed Y coefficients are replayed
as sign and congruence checks.

If no denominator anchor works, the script reports exact denominator sizes
and missing X indices for the next low-dimensional simultaneous-LLL step.

Run:
  sage -python ~/Downloads/recover_h92_q6_third_common_denominator.sage \
    artifacts/local/elkies-k3/q6-third-hensel-p1024.json
"""

import argparse
import json
from math import gcd
from pathlib import Path

from sage.all import QQ, ZZ, PolynomialRing, lcm


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
        else Path("artifacts/local/elkies-k3/q6-third-exact-rational.json")
    )
)

data = json.loads(INPUT.read_text())
child = json.loads(CHILD.read_text())
assert data["schema"] == "elkies-k3.h92-q6-third-hensel-lift.v1"
assert data["status"] == "PASS_Q6_THIRD_HENSEL"
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"

p = ZZ(data["prime"])
precision = int(data["precision"])
M = p**precision
residues = [ZZ(value) % M for value in data["residues"]]
assert len(residues) == 138

# Unknown-vector block layout in the Hensel system.
Zres = residues[:21]
Xres = residues[21:68]
Yres = residues[68:138]
assert len(Zres) == 21 and len(Xres) == 47 and len(Yres) == 70


def digits(value):
    value = abs(ZZ(value))
    return 1 if not value else len(str(value))


def balanced(residue):
    try:
        return ZZ(residue).rational_reconstruction(M)
    except (ArithmeticError, ValueError):
        return None


Zknown = [balanced(value) for value in Zres]
Xknown = [balanced(value) for value in Xres]
Yknown = [balanced(value) for value in Yres]

assert all(value is not None for value in Zknown)
missing_x = [i for i, value in enumerate(Xknown) if value is None]
missing_y = [i for i, value in enumerate(Yknown) if value is None]

print(
    "Q6THIRDCOMMON|"
    f"Z={sum(v is not None for v in Zknown)}/21|"
    f"X={sum(v is not None for v in Xknown)}/47|"
    f"Y={sum(v is not None for v in Yknown)}/70|"
    f"missing_X={','.join(map(str, missing_x))}|"
    f"missing_Y_count={len(missing_y)}",
    flush=True,
)


def denominator_lcm(values):
    answer = ZZ(1)
    for value in values:
        if value is not None:
            answer = lcm(answer, ZZ(value.denominator()))
    return answer


Lz = denominator_lcm(Zknown)
Lx = denominator_lcm(Xknown)
Ly = denominator_lcm(Yknown)

Acoeff = [
    QQ(value)
    for value in child["minimal_short_weierstrass"]["A_coefficients_low_to_high"]
]
Bcoeff = [
    QQ(value)
    for value in child["minimal_short_weierstrass"]["B_coefficients_low_to_high"]
]
Lchild = denominator_lcm(Acoeff + Bcoeff)

anchors = []
def add_anchor(name, value):
    value = ZZ(value)
    if value <= 0 or value % p == 0:
        return
    if all(existing != value for unused, existing in anchors):
        anchors.append((name, value))

add_anchor("Lx", Lx)
add_anchor("lcm_XY", lcm(Lx, Ly))
add_anchor("lcm_XYZchild", lcm(lcm(Lx, Ly), lcm(Lz, Lchild)))
# Extra exponent coverage without an uncontrolled huge product.
add_anchor("lcm_XY_times_X_gcdsafe", lcm(Lx*Lx, Ly))
add_anchor("lcm_XY_times_Y_gcdsafe", lcm(Lx, Ly*Ly))

print(
    "Q6THIRDCOMMON_DENOMS|"
    f"Lz_digits={digits(Lz)}|Lx_digits={digits(Lx)}|"
    f"Ly_digits={digits(Ly)}|Lchild_digits={digits(Lchild)}|"
    f"modulus_digits={digits(M)}",
    flush=True,
)

TR = PolynomialRing(QQ, "T")
T = TR.gen()
Z = TR(list(Zknown) + [QQ(1)])
assert Z.degree() == 21 and Z.leading_coefficient() == 1
Apoly = TR(Acoeff)
Bpoly = TR(Bcoeff)


def centered(value):
    value = ZZ(value) % M
    return value if value <= M//2 else value-M


def residue_of(q):
    q = QQ(q)
    den = ZZ(q.denominator())
    if den % p == 0:
        raise ZeroDivisionError("candidate denominator divisible by p")
    return (ZZ(q.numerator()) * den.inverse_mod(M)) % M


def y_matches_padic(Y):
    # Up to global sign. Check all 70 residues directly.
    coeff = list(Y.list()) + [QQ(0)] * max(0, 70-len(Y.list()))
    coeff = coeff[:70]
    direct = all(residue_of(coeff[i]) == Yres[i] for i in range(70))
    negative = all(residue_of(-coeff[i]) == Yres[i] for i in range(70))
    if direct:
        return Y, "direct"
    if negative:
        return -Y, "negated"
    return None, None


success = None
for name, D in anchors:
    # D must be invertible mod p.
    assert gcd(int(D), int(p)) == 1
    Xcoeff = []
    missing_num_digits = []
    for index, residue in enumerate(Xres):
        if Xknown[index] is not None:
            Xcoeff.append(QQ(Xknown[index]))
            continue
        numerator = centered((D % M) * residue)
        missing_num_digits.append(digits(numerator))
        Xcoeff.append(QQ(numerator) / QQ(D))

    X = TR(Xcoeff)
    R = X**3 + Apoly*X*Z**4 + Bpoly*Z**6
    square = R.is_square()

    print(
        f"Q6THIRDCOMMON_TRY|anchor={name}|D_digits={digits(D)}|"
        f"missing_num_digits="
        f"{min(missing_num_digits) if missing_num_digits else 0}.."
        f"{max(missing_num_digits) if missing_num_digits else 0}|"
        f"rhs_square={int(bool(square))}",
        flush=True,
    )

    if not square:
        continue

    Ycandidate = R.sqrt()
    Yexact, sign = y_matches_padic(Ycandidate)
    if Yexact is None:
        print(
            f"Q6THIRDCOMMON_TRY|anchor={name}|square=1|"
            "padic_Y_match=0",
            flush=True,
        )
        continue

    assert Yexact.degree() == 69
    assert Yexact**2 == R

    # Replay all p-adic residues for X as well.
    if any(residue_of(Xcoeff[i]) != Xres[i] for i in range(47)):
        raise ArithmeticError("exact-square X candidate missed its p-adic seed")

    success = (name, D, X, Yexact, sign)
    break


if success is None:
    print(
        "Q6THIRDCOMMON_RESULT|status=NEEDS_LOW_DIMENSION_SIMULTANEOUS_LLL|"
        f"missing_X={len(missing_x)}|missing_Y={len(missing_y)}",
        flush=True,
    )
    raise SystemExit(0)

name, D, X, Y, sign = success
identity = Y**2 - X**3 - Apoly*X*Z**4 - Bpoly*Z**6
assert identity == 0

all_values = list(Z.list()) + list(X.list()) + list(Y.list())
max_num = max(digits(value.numerator()) for value in all_values)
max_den = max(digits(value.denominator()) for value in all_values)

print(
    f"Q6THIRDCOMMON_EXACT|anchor={name}|Y_sign={sign}|"
    f"Z_degree={Z.degree()}|X_degree={X.degree()}|Y_degree={Y.degree()}|"
    f"max_num_digits={max_num}|max_den_digits={max_den}|identity=PASS",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-q6-third-exact-rational.v1",
    "status": "PASS_EXACT_Q6_THIRD_RATIONAL_SECTION",
    "source_hensel": str(INPUT.relative_to(ROOT)),
    "prime": int(p),
    "precision": precision,
    "reconstruction": {
        "method": "known-X common denominator plus exact square-root verification",
        "denominator_anchor": name,
        "anchor_decimal_digits": digits(D),
        "Y_sign_against_padic_seed": sign,
    },
    "profile": {
        "Z_degree": int(Z.degree()),
        "X_degree": int(X.degree()),
        "Y_degree": int(Y.degree()),
    },
    "max_numerator_decimal_digits": int(max_num),
    "max_denominator_decimal_digits": int(max_den),
    "Z": [str(value) for value in list(Z.list())],
    "X": [str(value) for value in list(X.list())],
    "Y": [str(value) for value in list(Y.list())],
    "exact_identity": "Y^2=X^3+A*X*Z^4+B*Z^6",
    "boundary": (
        "This certifies the explicit characteristic-zero third q6 section. "
        "Its exact Shioda/component profile and downstream q8 Abel-Jacobi "
        "transport remain separate certificates."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q6THIRDCOMMON_RESULT|status=PASS_EXACT_Q6_THIRD_RATIONAL_SECTION",
    flush=True,
)
