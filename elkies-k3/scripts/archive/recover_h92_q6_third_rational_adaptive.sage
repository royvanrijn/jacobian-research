#!/usr/bin/env sage -python
"""
Adaptive rational reconstruction for the q6-third p-adic Hensel solution.

Sage's ordinary rational_reconstruction uses essentially balanced bounds.
At p^1024 this can still fail for a perfectly rational coefficient if its
numerator and denominator sizes are strongly asymmetric.

This script:
  1. reports ordinary balanced reconstruction progress;
  2. tries provably unique generalized rational reconstruction in many
     unbalanced rectangles |num|<=A, 0<den<=B with 2*A*B < modulus;
  3. requires a unique candidate per coefficient;
  4. if all 138 coefficients recover, verifies exactly over QQ[T]:
         Y^2 = X^3 + A*X*Z^4 + B*Z^6
     and writes a compact exact artifact.

If coefficients remain unresolved at p^1024, do NOT silently accept a guess;
the output identifies them for the next simultaneous-LLL stage.

Run:
  sage -python ~/Downloads/recover_h92_q6_third_rational_adaptive.sage \
    artifacts/local/elkies-k3/q6-third-hensel-p1024.json
"""

import argparse
import json
from math import gcd
from pathlib import Path

from sage.all import QQ, ZZ, PolynomialRing


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
parser.add_argument(
    "--output",
    type=Path,
    help="default: artifacts/local/elkies-k3/q6-third-exact-rational.json",
)
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

blocks = (
    ("Z", 0, 21),
    ("X", 21, 68),
    ("Y", 68, 138),
)


def label(index):
    for name, lo, hi in blocks:
        if lo <= index < hi:
            return f"{name}[{index-lo}]"
    raise AssertionError(index)


def decimal_digits(value):
    value = abs(ZZ(value))
    return 1 if not value else len(str(value))


def bounded_rational_reconstruction(x, modulus, numerator_bound, denominator_bound):
    """
    Generalized rational reconstruction in a prescribed uniqueness rectangle.

    Return a/b satisfying
        a == b*x (mod modulus),
        |a| <= numerator_bound,
        0 < b <= denominator_bound,
    provided 2*A*B < modulus and the Euclidean reconstruction lands there.
    """
    modulus = ZZ(modulus)
    x = ZZ(x) % modulus
    A = ZZ(numerator_bound)
    B = ZZ(denominator_bound)
    if A <= 0 or B <= 0 or 2*A*B >= modulus:
        raise ValueError("bounds do not define a unique reconstruction rectangle")
    if x == 0:
        return QQ(0)

    r0, r1 = modulus, x
    t0, t1 = ZZ(0), ZZ(1)

    while abs(r1) > A:
        if r1 == 0:
            return None
        q = r0 // r1
        r0, r1 = r1, r0 - q*r1
        t0, t1 = t1, t0 - q*t1

    a, b = ZZ(r1), ZZ(t1)
    if b < 0:
        a, b = -a, -b

    if not (0 < b <= B and abs(a) <= A):
        return None
    if gcd(int(abs(a)), int(b)) != 1:
        return None
    if gcd(int(b), int(modulus)) != 1:
        return None
    if (a - b*x) % modulus:
        return None
    return QQ(a) / QQ(b)


# First show ordinary balanced progress at useful precision levels.
levels = sorted(set(
    [k for k in (64, 128, 256, 512, 1024) if k <= precision] + [precision]
))
for k in levels:
    modulus = p**k
    recovered = []
    max_num = max_den = 0
    for i, residue in enumerate(residues):
        try:
            q = ZZ(residue % modulus).rational_reconstruction(modulus)
        except (ArithmeticError, ValueError):
            continue
        recovered.append(i)
        max_num = max(max_num, decimal_digits(q.numerator()))
        max_den = max(max_den, decimal_digits(q.denominator()))
    parts = []
    for name, lo, hi in blocks:
        count = sum(lo <= index < hi for index in recovered)
        parts.append(f"{name}={count}/{hi-lo}")
    print(
        f"Q6THIRDADAPT_BALANCED|precision={k}|recovered={len(recovered)}/138|"
        + "|".join(parts)
        + f"|max_num_digits={max_num}|max_den_digits={max_den}",
        flush=True,
    )

# Adaptive bit splits.  Leave several safety bits so 2AB < M strictly.
mbits = int(M.nbits())
budget = mbits - 4
# Dense around balanced, but also strongly asymmetric.
fractions = [
    0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20, 0.225, 0.25,
    0.275, 0.30, 0.325, 0.35, 0.375, 0.40, 0.425, 0.45, 0.475,
    0.50,
    0.525, 0.55, 0.575, 0.60, 0.625, 0.65, 0.675, 0.70, 0.725,
    0.75, 0.775, 0.80, 0.825, 0.85, 0.875, 0.90, 0.925, 0.95,
]

candidates = []
ambiguous = {}
for index, residue in enumerate(residues):
    found = {}

    # Include Sage's balanced result if available.
    try:
        q = ZZ(residue).rational_reconstruction(M)
        found[str(q)] = q
    except (ArithmeticError, ValueError):
        pass

    for fraction in fractions:
        nbits = max(1, int(budget*fraction))
        dbits = max(1, budget-nbits)
        A = ZZ(1) << nbits
        B = ZZ(1) << dbits
        try:
            q = bounded_rational_reconstruction(residue, M, A, B)
        except ValueError:
            continue
        if q is not None:
            found[str(q)] = q

    values = list(found.values())
    if len(values) == 1:
        candidates.append(values[0])
    else:
        candidates.append(None)
        if len(values) > 1:
            ambiguous[index] = values

recovered_indices = [i for i, q in enumerate(candidates) if q is not None]
parts = []
for name, lo, hi in blocks:
    count = sum(lo <= index < hi for index in recovered_indices)
    parts.append(f"{name}={count}/{hi-lo}")
print(
    f"Q6THIRDADAPT|recovered={len(recovered_indices)}/138|"
    + "|".join(parts)
    + f"|ambiguous={len(ambiguous)}",
    flush=True,
)

if ambiguous:
    for index, values in sorted(ambiguous.items())[:10]:
        print(
            f"Q6THIRDADAPT_AMBIGUOUS|{label(index)}|"
            + ";".join(map(str, values)),
            flush=True,
        )

failed = [
    label(index)
    for index, value in enumerate(candidates)
    if value is None
]
print(
    "Q6THIRDADAPT_FAILED|"
    + ("none" if not failed else ",".join(failed)),
    flush=True,
)

if failed or ambiguous:
    print(
        "Q6THIRDADAPT_RESULT|status=NEEDS_SIMULTANEOUS_RECONSTRUCTION|"
        f"recovered={len(recovered_indices)}/138",
        flush=True,
    )
    raise SystemExit(0)

# Rebuild exact section and verify globally over QQ[T].
Zcoeff = list(candidates[:21]) + [QQ(1)]
Xcoeff = list(candidates[21:68])
Ycoeff = list(candidates[68:138])
assert len(Zcoeff) == 22 and len(Xcoeff) == 47 and len(Ycoeff) == 70

TR = PolynomialRing(QQ, "T")
T = TR.gen()
Z = TR(Zcoeff)
X = TR(Xcoeff)
Y = TR(Ycoeff)
A = TR([
    QQ(value)
    for value in child["minimal_short_weierstrass"]["A_coefficients_low_to_high"]
])
B = TR([
    QQ(value)
    for value in child["minimal_short_weierstrass"]["B_coefficients_low_to_high"]
])

assert Z.degree() == 21 and Z.leading_coefficient() == 1
assert X.degree() == 46
assert Y.degree() == 69
identity = Y**2 - X**3 - A*X*Z**4 - B*Z**6
assert identity == 0

max_num = max(
    decimal_digits(value.numerator())
    for value in Zcoeff + Xcoeff + Ycoeff
)
max_den = max(
    decimal_digits(value.denominator())
    for value in Zcoeff + Xcoeff + Ycoeff
)

payload = {
    "schema": "elkies-k3.h92-q6-third-exact-rational.v1",
    "status": "PASS_EXACT_Q6_THIRD_RATIONAL_SECTION",
    "source_hensel": str(INPUT.relative_to(ROOT)),
    "prime": int(p),
    "precision": precision,
    "profile": {
        "Z_degree": 21,
        "X_degree": 46,
        "Y_degree": 69,
    },
    "max_numerator_decimal_digits": int(max_num),
    "max_denominator_decimal_digits": int(max_den),
    "Z": [str(value) for value in Zcoeff],
    "X": [str(value) for value in Xcoeff],
    "Y": [str(value) for value in Ycoeff],
    "exact_identity": "Y^2=X^3+A*X*Z^4+B*Z^6",
    "boundary": (
        "This certifies the explicit characteristic-zero q6 third section on "
        "the pinned child Weierstrass model. Shioda/component and downstream "
        "q8 Abel-Jacobi checks remain separate."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    f"Q6THIRDADAPT_EXACT|max_num_digits={max_num}|max_den_digits={max_den}|"
    "identity=PASS",
    flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q6THIRDADAPT_RESULT|status=PASS_EXACT_Q6_THIRD_RATIONAL_SECTION",
    flush=True,
)
