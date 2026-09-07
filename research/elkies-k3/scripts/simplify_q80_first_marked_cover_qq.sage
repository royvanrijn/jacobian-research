#!/usr/bin/env sage
"""Simplify the first Q80 marking cover and run a bounded rational-point search.

The global coefficient parameter used by the Q80 reconstruction makes the
degree-six marking cover look enormous.  Two rational branch points give a
canonical Möbius chart in which the remaining quadratics have the exact
shape

    s^2+s+1/3,       s^2+2s+2.

After one rational square rescaling the cover is therefore

    Y^2 = 2*s*(s^2+s+1/3)*(s^2+2*s+2).

This script proves the birational identity over Q and optionally asks PARI's
``hyperellratpoints`` for all affine points in one declared height box.  A
completed bounded search is evidence, not a proof of the full rational point
set of this genus-two curve.
"""

import argparse
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path

from sage.all import PolynomialRing, QQ, gcd, lcm


ROOT = Path(__file__).resolve().parents[2]
INPUT = (
    ROOT / "artifacts/generated-results/q80-slope-8-87-first-marked-cover-qq.json"
)
OUTPUT = (
    ROOT / "artifacts/generated-results/q80-first-marked-cover-simplified-qq.json"
)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=INPUT)
parser.add_argument("--output", type=Path, default=OUTPUT)
parser.add_argument("--search-height", type=int, default=100_000_000)
parser.add_argument("--gp-timeout", type=int, default=60)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
if args.search_height < 1 or args.gp_timeout < 1:
    raise ValueError("search height and timeout must be positive")


def sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def sha256(path):
    return sha256_bytes(path.read_bytes())


def rational_record(value):
    value = QQ(value)
    return {
        "value": str(value),
        "numerator_bits": int(abs(value.numerator()).nbits()),
        "denominator_bits": int(value.denominator().nbits()),
    }


started = time.monotonic()
source_raw = args.input.read_bytes()
source = json.loads(source_raw)
if source.get("schema") != "q80-first-marked-cover-qq-v1":
    raise ValueError("unexpected first-marked-cover schema")
if source.get("status") != "PASS_EXACT_FIRST_MARKED_COVER":
    raise ValueError("the source marked-cover certificate has not passed")

ring = PolynomialRing(QQ, "x")
x = ring.gen()
cover = ring(source["squarefree_cover"]["polynomial"].replace("t", "x"))
if cover.degree() != 6 or cover.gcd(cover.derivative()) != 1:
    raise ArithmeticError("source cover is not a squarefree sextic")

factorization = tuple(cover.factor())
if sorted((factor.degree(), int(exponent)) for factor, exponent in factorization) != [
    (1, 1),
    (1, 1),
    (2, 1),
    (2, 1),
]:
    raise ArithmeticError("source cover does not have the expected 1,1,2,2 split")
branch_roots = [
    -factor[0] / factor[1]
    for factor, _ in factorization
    if factor.degree() == 1
]
r0, r_infinity = branch_roots

# Put r0 at z=0 and r_infinity at z=infinity:
# z=(t-r0)/(t-r_infinity), t=(r_infinity*z-r0)/(z-1).
z = x
mobius_numerator = r_infinity * z - r0
mobius_denominator = z - 1
transformed = ring(
    sum(
        cover[index]
        * mobius_numerator**index
        * mobius_denominator ** (6 - index)
        for index in range(7)
    )
)
if transformed.degree() != 5 or transformed[0] != 0:
    raise ArithmeticError("Möbius normalization did not move one branch to infinity")

transformed_factors = tuple(transformed.factor())
quadratics = sorted(
    (factor / factor[2] for factor, _ in transformed_factors if factor.degree() == 2),
    key=lambda factor: factor[1],
)
if len(quadratics) != 2:
    raise ArithmeticError("normalized cover does not have two quadratic factors")
first_quadratic, second_quadratic = quadratics
a = QQ(first_quadratic[1])
if not a:
    raise ArithmeticError("quadratic scaling coefficient vanished")
if first_quadratic != z**2 + a * z + a**2 / 3:
    raise ArithmeticError("first normalized quadratic has the wrong exact shape")
if second_quadratic != z**2 + 2 * a * z + 2 * a**2:
    raise ArithmeticError("second normalized quadratic has the wrong exact shape")

leading_unit = QQ(transformed.leading_coefficient())
if transformed != leading_unit * z * first_quadratic * second_quadratic:
    raise ArithmeticError("normalized factorization identity failed")

# Substituting z=a*s contributes a^5.  Its squareclass is exactly 2.
cover_scale_squared = QQ(leading_unit * a**5 / 2)
if not cover_scale_squared.is_square():
    raise ArithmeticError("normalized cover scale does not have squareclass 2")
cover_scale = QQ(cover_scale_squared.sqrt())
s = x
simplified = ring(2 * s * (s**2 + s + QQ(1) / 3) * (s**2 + 2 * s + 2))
substituted = ring(transformed(a * s))
if substituted != cover_scale**2 * simplified:
    raise ArithmeticError("simplified genus-two cover identity failed")
if simplified.gcd(simplified.derivative()) != 1 or simplified.degree() != 5:
    raise ArithmeticError("simplified cover is not a squarefree quintic")

# CM24 is t=0 in the source local parameter.
cm24_s = QQ((r0 / r_infinity) / a)
if cm24_s != -1:
    raise ArithmeticError("CM24 did not map to s=-1")
if simplified(cm24_s) != -QQ(2) / 3:
    raise ArithmeticError("CM24 simplified squareclass changed")

gp = shutil.which("gp")
if gp is None:
    raise RuntimeError("PARI/GP executable 'gp' is required")
gp_polynomial = "+".join(
    f"({coefficient.numerator()}/{coefficient.denominator()})*x^{index}"
    for index, coefficient in enumerate(simplified.list())
)
program = (
    "x='x;Q="
    + gp_polynomial
    + f";gettime();R=hyperellratpoints(Q,{args.search_height});"
    + 'print("PARI_MILLISECONDS|",gettime());'
    + 'for(i=1,#R,print(R[i][1],"\\t",R[i][2]));\n'
)
search_started = time.monotonic()
try:
    completed = subprocess.run(
        [gp, "-q", "-f"],
        input=program,
        text=True,
        capture_output=True,
        timeout=args.gp_timeout,
        check=False,
    )
except subprocess.TimeoutExpired as error:
    raise RuntimeError(
        f"PARI search exceeded the declared {args.gp_timeout}-second timeout"
    ) from error
if completed.returncode != 0 or "***" in completed.stdout + completed.stderr:
    raise RuntimeError("PARI hyperellratpoints failed: " + completed.stderr[-1000:])
lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
if not lines or not lines[0].startswith("PARI_MILLISECONDS|"):
    raise RuntimeError("PARI search did not report its runtime")
pari_milliseconds = int(lines[0].split("|", 1)[1])
affine_points = []
for line in lines[1:]:
    fields = line.split("\t")
    if len(fields) != 2:
        raise RuntimeError(f"unexpected PARI point row: {line!r}")
    point = [str(QQ(fields[0])), str(QQ(fields[1]))]
    if QQ(point[1]) ** 2 != simplified(QQ(point[0])):
        raise ArithmeticError("PARI point failed exact substitution")
    affine_points.append(point)

status = "PASS_EXACT_SIMPLIFIED_COVER_AND_BOUNDED_SEARCH"
output = {
    "schema": "q80-first-marked-cover-simplified-qq-v1",
    "status": status,
    "input": {
        "path": str(args.input.resolve().relative_to(ROOT)),
        "sha256": sha256(args.input),
    },
    "exact_birational_simplification": {
        "source": "w^2=f(t)",
        "first_rational_branch_t": rational_record(r0),
        "second_rational_branch_t": rational_record(r_infinity),
        "mobius_coordinate": "z=(t-r0)/(t-r_infinity)",
        "scaled_coordinate": "z=a*s",
        "a": rational_record(a),
        "inverse_parameter_map": "t=(a*s*r_infinity-r0)/(a*s-1)",
        "ordinate_map": "Y=(a*s-1)^3*w/cover_scale",
        "cover_scale": rational_record(cover_scale),
        "simplified_equation": "Y^2=2*s*(s^2+s+1/3)*(s^2+2*s+2)",
        "simplified_coefficients_low_to_high": list(map(str, simplified.list())),
        "squarefree": True,
        "degree": 5,
        "genus": 2,
        "cm24_s": str(cm24_s),
        "cm24_rhs": str(simplified(cm24_s)),
    },
    "bounded_rational_point_search": {
        "engine": "PARI/GP hyperellratpoints",
        "height_bound": int(args.search_height),
        "timeout_seconds": int(args.gp_timeout),
        "completed": True,
        "program_sha256": sha256_bytes(program.encode()),
        "pari_milliseconds": pari_milliseconds,
        "wall_seconds": time.monotonic() - search_started,
        "affine_points": affine_points,
        "affine_point_count": len(affine_points),
        "rational_point_at_infinity": True,
        "nonbranch_affine_points_found": sum(QQ(row[1]) != 0 for row in affine_points),
    },
    "claim_boundary": {
        "proved": [
            "exact Q-birational simplification of the first marking cover",
            "CM24 maps to s=-1 with nonsquare rational right side -2/3",
            "the displayed bounded PARI search completed",
        ],
        "not_proved": [
            "the complete rational point set of the genus-two cover",
            "nonexistence of a rational first marking outside the search bound",
            "the third-q12 horizontal or child equation",
        ],
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/simplify_q80_first_marked_cover_qq.sage "
        f"--search-height {args.search_height} --gp-timeout {args.gp_timeout}"
    ),
    "runtime_seconds": time.monotonic() - started,
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists():
        raise SystemExit(f"missing simplified-cover artifact: {args.output}")
    existing = json.loads(args.output.read_text())
    for payload in (existing, output):
        payload.pop("runtime_seconds", None)
        search = payload["bounded_rational_point_search"]
        search.pop("wall_seconds", None)
        search.pop("pari_milliseconds", None)
    if existing != output:
        raise SystemExit("stale simplified-cover artifact")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)

print(
    "Q80MARKEDCOVERSIMPLE|equation=2s(s2+s+1/3)(s2+2s+2)|"
    f"height={args.search_height}|affine_points={len(affine_points)}|"
    f"nonbranch={sum(QQ(row[1]) != 0 for row in affine_points)}|"
    f"status={status}",
    flush=True,
)
