#!/usr/bin/env sage -python
"""
Hensel-lift the recovered q6 third MW generator from its modular trace seed.

The modular trace/halving reconstruction gives the exact pole profile

    x = X/Z^2,  deg(X,Z) = (46,21),
    y = Y/Z^3,  deg(Y,Z) = (69,21),

on the certified H92 q6 child

    y^2 = x^3 + A(T)*x + B(T).

With Z monic, the cleared identity

    Y^2 - X^3 - A*X*Z^4 - B*Z^6 = 0

has 138 unknown coefficients:
  * 21 coefficients of Z below its fixed monic T^21 term,
  * 47 coefficients of X,
  * 70 coefficients of Y.

There are 139 coefficient equations (degrees 0..138).  At a nonsingular
modular seed the Jacobian has rank 138; selecting a full-rank 138-row minor
gives a square p-adic Newton system.

The script supports restartable p-adic lifting via --seed, matching the
successful original H92 P2 reconstruction workflow.

First run:
  sage -python ~/Downloads/lift_h92_q6_third_hensel.sage --rank-only

Then:
  sage -python ~/Downloads/lift_h92_q6_third_hensel.sage \
    --precision 32 \
    --output artifacts/local/elkies-k3/q6-third-hensel-p32.json

If incomplete, continue:
  sage -python ~/Downloads/lift_h92_q6_third_hensel.sage \
    --precision 64 \
    --seed artifacts/local/elkies-k3/q6-third-hensel-p32.json \
    --output artifacts/local/elkies-k3/q6-third-hensel-p64.json
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, QQ, ZZ, PolynomialRing, Zp, matrix, vector


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
parser.add_argument(
    "--input",
    type=Path,
    help=(
        "modular trace seed; default is "
        "artifacts/local/elkies-k3/q6-third-trace-samples-mod-100003.json"
    ),
)
parser.add_argument("--precision", type=int, default=32, help="target p-adic digits")
parser.add_argument("--output", type=Path)
parser.add_argument("--seed", type=Path, help="previous p-adic residue output")
parser.add_argument("--rank-only", action="store_true")
args = parser.parse_args()

ROOT = locate_repo(args.repo)
INPUT = (
    args.input.resolve()
    if args.input
    else ROOT / "artifacts/local/elkies-k3/q6-third-trace-samples-mod-100003.json"
)
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"

for path in (INPUT, CHILD):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

record = json.loads(INPUT.read_text())
child = json.loads(CHILD.read_text())

assert record["status"] == "PASS_MODULAR_Q6_THIRD_TRACE_SECTION"
assert record["schema"] == "elkies-k3.h92-q6-third-trace-samples-modp.v1"
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"

prime = ZZ(record["prime"])
assert prime == 100003
finite = GF(prime)
finite_ring = PolynomialRing(finite, "T")
Tf = finite_ring.gen()

profile = record["profile"]
assert (
    int(profile["Z_degree"]),
    int(profile["X_degree"]),
    int(profile["Y_degree"]),
    int(profile["x_denominator_degree"]),
) == (21, 46, 69, 42)

Zraw = [ZZ(value) for value in record["Z"]]
Xraw = [ZZ(value) for value in record["X"]]
Yraw = [ZZ(value) for value in record["Y"]]

# Normalize/pad exact declared degrees.
if len(Zraw) > 22 or len(Xraw) > 47 or len(Yraw) > 70:
    raise ValueError("modular seed exceeds the certified degree profile")
Zraw += [ZZ(0)] * (22 - len(Zraw))
Xraw += [ZZ(0)] * (47 - len(Xraw))
Yraw += [ZZ(0)] * (70 - len(Yraw))
assert len(Zraw) == 22 and len(Xraw) == 47 and len(Yraw) == 70
assert finite(Zraw[21]) == 1  # monic Z.

zf = finite_ring([finite(value) for value in Zraw])
xf = finite_ring([finite(value) for value in Xraw])
yf = finite_ring([finite(value) for value in Yraw])
assert zf.degree() == 21 and zf.leading_coefficient() == 1
assert xf.degree() == 46
assert yf.degree() == 69

A_QQ = [
    QQ(value)
    for value in child["minimal_short_weierstrass"]["A_coefficients_low_to_high"]
]
B_QQ = [
    QQ(value)
    for value in child["minimal_short_weierstrass"]["B_coefficients_low_to_high"]
]
assert len(A_QQ) <= 9 and len(B_QQ) <= 13


def polynomial_coefficients(value, degree):
    actual_degree = value.degree()
    return [
        value[index] if actual_degree >= 0 and index <= actual_degree
        else value.base_ring()(0)
        for index in range(degree + 1)
    ]


def child_coefficients(ring):
    T = ring.gen()
    a = sum(ring(value) * T**index for index, value in enumerate(A_QQ))
    b = sum(ring(value) * T**index for index, value in enumerate(B_QQ))
    return a, b


def identity_and_jacobian(ring, z_values, x_values, y_values):
    """
    Cleared section identity and its 139x138 coefficient Jacobian.

    z_values contains coefficients T^0..T^20; T^21 is fixed to 1.
    """
    T = ring.gen()
    z = sum(ring(z_values[index]) * T**index for index in range(21)) + T**21
    x = sum(ring(x_values[index]) * T**index for index in range(47))
    y = sum(ring(y_values[index]) * T**index for index in range(70))
    a, b = child_coefficients(ring)

    identity = y**2 - x**3 - a*x*z**4 - b*z**6

    derivatives = []
    # d/d z_i
    for index in range(21):
        derivatives.append(
            (-4*a*x*z**3 - 6*b*z**5) * T**index
        )
    # d/d x_i
    for index in range(47):
        derivatives.append(
            (-3*x**2 - a*z**4) * T**index
        )
    # d/d y_i
    for index in range(70):
        derivatives.append(
            2*y*T**index
        )
    assert len(derivatives) == 138

    return (
        vector(
            ring.base_ring(),
            polynomial_coefficients(identity, 138),
        ),
        matrix(
            ring.base_ring(),
            [
                [
                    polynomial_coefficients(derivative, 138)[row]
                    for derivative in derivatives
                ]
                for row in range(139)
            ],
        ),
    )


z0 = [finite(value) for value in Zraw[:21]]
x0 = [finite(value) for value in Xraw]
y0 = [finite(value) for value in Yraw]

residue, jacobian = identity_and_jacobian(finite_ring, z0, x0, y0)
assert not residue

rank = jacobian.rank()
print(
    f"Q6THIRDHENSEL|prime={prime}|jacobian_rank={rank}|"
    f"dimensions={jacobian.dimensions()}",
    flush=True,
)
assert rank == 138

# Pivot columns of J^T are independent rows of J.
rows = tuple(int(row) for row in jacobian.transpose().pivots())
assert len(rows) == 138
square = jacobian.matrix_from_rows(rows)
assert square.nrows() == square.ncols() == 138
assert square.is_invertible()

print(
    f"Q6THIRDHENSEL|selected_rows={len(rows)}|"
    "status=PASS_NONSINGULAR_MODULAR_SEED",
    flush=True,
)

if args.rank_only:
    raise SystemExit(0)

target = int(args.precision)
if target < 2:
    raise ValueError("precision must be at least 2")

padic = Zp(prime, prec=target)
padic_ring = PolynomialRing(padic, "T")

# Start from the integer representatives of the finite-field seed.  Coercing
# GF(p) elements directly to Zp would retain only O(p) precision metadata.
coefficients = [
    padic(ZZ(value))
    for value in Zraw[:21] + Xraw + Yraw
]
assert len(coefficients) == 138

seed_precision = 1
if args.seed:
    seed_path = args.seed.resolve()
    seed = json.loads(seed_path.read_text())
    if (
        seed.get("schema") != "elkies-k3.h92-q6-third-hensel-lift.v1"
        or ZZ(seed.get("prime")) != prime
        or "residues" not in seed
    ):
        raise ValueError("seed is not a compatible q6-third Hensel artifact")
    seed_precision = int(seed["precision"])
    if seed_precision >= target:
        raise ValueError("seed precision must be smaller than target precision")
    if len(seed["residues"]) != 138:
        raise ValueError("seed residue count is inconsistent")
    coefficients = [padic(ZZ(value)) for value in seed["residues"]]
    print(
        f"Q6THIRDHENSEL|seed_precision={seed_precision}|"
        f"target_precision={target}",
        flush=True,
    )


def valuation_floor(values):
    finite_values = [value.valuation() for value in values if value]
    return target if not finite_values else min(finite_values)


for iteration in range(1, 2*target + 3):
    residual, derivative = identity_and_jacobian(
        padic_ring,
        coefficients[:21],
        coefficients[21:68],
        coefficients[68:],
    )
    valuation = valuation_floor(residual)
    print(
        f"Q6THIRDHENSEL|iteration={iteration}|"
        f"residual_valuation={valuation}",
        flush=True,
    )
    if valuation >= target:
        break

    selected_derivative = derivative.matrix_from_rows(rows)
    correction = selected_derivative.solve_right(
        -vector(padic, [residual[row] for row in rows])
    )
    coefficients = [
        value + delta
        for value, delta in zip(coefficients, correction, strict=True)
    ]
else:
    raise ArithmeticError("p-adic Newton iteration did not reach target precision")

# Recheck all 139 equations, not only the selected square subsystem.
final_residual, unused_derivative = identity_and_jacobian(
    padic_ring,
    coefficients[:21],
    coefficients[21:68],
    coefficients[68:],
)
final_valuation = valuation_floor(final_residual)
assert final_valuation >= target

modulus = prime**target


def reconstruct(values):
    answer = []
    for value in values:
        try:
            answer.append(
                str(ZZ(value.lift()).rational_reconstruction(modulus))
            )
        except (ArithmeticError, ValueError):
            answer.append(None)
    return answer


payload = {
    "schema": "elkies-k3.h92-q6-third-hensel-lift.v1",
    "status": "PASS_Q6_THIRD_HENSEL",
    "prime": int(prime),
    "precision": target,
    "seed_precision": int(seed_precision),
    "jacobian_rank": int(rank),
    "selected_rows": list(rows),
    "final_residual_valuation": int(final_valuation),
    "complete": False,
    "profile": {
        "Z_degree": 21,
        "X_degree": 46,
        "Y_degree": 69,
    },
    "Z": reconstruct(coefficients[:21]) + ["1"],
    "X": reconstruct(coefficients[21:68]),
    "Y": reconstruct(coefficients[68:]),
    "residues": [str(ZZ(value.lift())) for value in coefficients],
    "source_modular_seed": str(INPUT.relative_to(ROOT)),
    "child_model": str(CHILD.relative_to(ROOT)),
}
payload["complete"] = not any(
    value is None
    for part in ("Z", "X", "Y")
    for value in payload[part]
)

if args.output:
    output = args.output
    if not output.is_absolute():
        output = ROOT / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"OUTPUT|{output}", flush=True)

print(
    f"Q6THIRDHENSEL|precision={target}|"
    f"complete={int(payload['complete'])}|"
    f"final_residual_valuation={final_valuation}|status=PASS",
    flush=True,
)
