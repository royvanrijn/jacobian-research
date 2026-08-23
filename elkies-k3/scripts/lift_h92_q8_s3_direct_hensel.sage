#!/usr/bin/env sage -python
"""
Hensel-lift the DIRECT branch-anchored H92 q8 S3 Abel-Jacobi section.

The direct modular trace/profile discovery certified over GF(100003):

    x = X/Z^2,  deg(X,Z) = (52,24),
    y = Y/Z^3,  deg(Y,Z) = (78,24),

on the certified minimal D13 q8 child

    y^2 = x^3 + A(U)*x + B(U).

With Z monic, the cleared identity

    Y^2 - X^3 - A*X*Z^4 - B*Z^6 = 0

has 156 unknown coefficients:
  * 24 coefficients of Z below fixed monic U^24,
  * 53 coefficients of X,
  * 79 coefficients of Y.

There are 157 coefficient equations (degrees 0..156).  If the modular
Jacobian has rank 156, a full-rank 156-row minor gives an ordinary p-adic
Newton/Hensel system.

The modular Y seed is one of the two signs; either lifts the same x-section
up to inversion.

This script is restartable with --seed.  If rational reconstruction becomes
complete, it substitutes the candidate back into the exact QQ(U) D13 equation.
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
        "direct modular profile seed; default "
        "artifacts/local/elkies-k3/q8-s3-direct-x-mod-100003.json"
    ),
)
parser.add_argument("--precision", type=int, default=32, help="target p-adic digits")
parser.add_argument("--output", type=Path)
parser.add_argument("--seed", type=Path, help="previous q8 S3 Hensel artifact")
parser.add_argument("--rank-only", action="store_true")
parser.add_argument("--exact-output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"
GEN = ROOT / "artifacts/generated-results"

INPUT = (
    args.input.resolve()
    if args.input
    else LOCAL / "q8-s3-direct-x-mod-100003.json"
)

q8_candidates = [
    LOCAL / "q8-corrected2cover-qq-child.json",
    GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
CHILD = next(
    (
        path for path in q8_candidates
        if path.exists()
        and "child" in json.loads(path.read_text())
        and "minimal_A_coefficients_low_to_high"
            in json.loads(path.read_text()).get("child", {})
    ),
    None,
)
if CHILD is None:
    raise SystemExit("No complete corrected q8 D13 child artifact found")

if not INPUT.exists():
    raise SystemExit(f"Missing prerequisite: {INPUT}")

record = json.loads(INPUT.read_text())
child = json.loads(CHILD.read_text())

assert record["schema"] == "elkies-k3.h92-q8-s3-direct-profile-discovery-modp.v1"
assert record["status"] == "PASS_DIRECT_ANCHORED_Q8_S3_PROFILE_DISCOVERY"
assert child["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"

prime = ZZ(record["prime"])
assert prime == 100003

xrec = record["x"]
wrec = record["weierstrass_structure"]
assert int(xrec["numerator_degree"]) == 52
assert int(xrec["denominator_degree"]) == 48
assert bool(wrec["denominator_is_square"])
assert int(wrec["denominator_root_degree"]) == 24
assert bool(wrec["curve_rhs_numerator_is_square"])
assert int(wrec["y_abs_numerator_degree"]) == 78
assert int(wrec["y_denominator_degree"]) == 72

Zraw = [ZZ(v) for v in wrec["denominator_root_coefficients_low_to_high"]]
Xraw = [ZZ(v) for v in xrec["numerator_coefficients_low_to_high"]]
Yraw = [ZZ(v) for v in wrec["y_abs_numerator_coefficients_low_to_high"]]

if len(Zraw) > 25 or len(Xraw) > 53 or len(Yraw) > 79:
    raise ValueError("modular seed exceeds discovered q8 degree profile")
Zraw += [ZZ(0)] * (25 - len(Zraw))
Xraw += [ZZ(0)] * (53 - len(Xraw))
Yraw += [ZZ(0)] * (79 - len(Yraw))
assert len(Zraw) == 25 and len(Xraw) == 53 and len(Yraw) == 79

finite = GF(prime)
finite_ring = PolynomialRing(finite, "U")
Uf = finite_ring.gen()

assert finite(Zraw[24]) == 1
zf = finite_ring([finite(v) for v in Zraw])
xf = finite_ring([finite(v) for v in Xraw])
yf = finite_ring([finite(v) for v in Yraw])
assert zf.degree() == 24 and zf.leading_coefficient() == 1
assert xf.degree() == 52
assert yf.degree() == 78

A_QQ = [QQ(v) for v in child["child"]["minimal_A_coefficients_low_to_high"]]
B_QQ = [QQ(v) for v in child["child"]["minimal_B_coefficients_low_to_high"]]
assert len(A_QQ) <= 9 and len(B_QQ) <= 13


def polynomial_coefficients(value, degree):
    actual_degree = value.degree()
    return [
        value[index] if actual_degree >= 0 and index <= actual_degree
        else value.base_ring()(0)
        for index in range(degree + 1)
    ]


def child_coefficients(ring):
    U = ring.gen()
    a = sum(ring(v) * U**i for i, v in enumerate(A_QQ))
    b = sum(ring(v) * U**i for i, v in enumerate(B_QQ))
    return a, b


def identity_and_jacobian(ring, z_values, x_values, y_values):
    """
    Cleared section identity and its 157x156 coefficient Jacobian.

    z_values are U^0..U^23; U^24 is fixed to 1.
    """
    U = ring.gen()
    z = sum(ring(z_values[i]) * U**i for i in range(24)) + U**24
    x = sum(ring(x_values[i]) * U**i for i in range(53))
    y = sum(ring(y_values[i]) * U**i for i in range(79))
    a, b = child_coefficients(ring)

    identity = y**2 - x**3 - a*x*z**4 - b*z**6

    derivatives = []
    for i in range(24):
        derivatives.append((-4*a*x*z**3 - 6*b*z**5) * U**i)
    for i in range(53):
        derivatives.append((-3*x**2 - a*z**4) * U**i)
    for i in range(79):
        derivatives.append(2*y * U**i)
    assert len(derivatives) == 156

    coeffs = lambda poly: polynomial_coefficients(poly, 156)
    return (
        vector(ring.base_ring(), coeffs(identity)),
        matrix(
            ring.base_ring(),
            [
                [coeffs(derivative)[row] for derivative in derivatives]
                for row in range(157)
            ],
        ),
    )


z0 = [finite(v) for v in Zraw[:24]]
x0 = [finite(v) for v in Xraw]
y0 = [finite(v) for v in Yraw]

residue, jacobian = identity_and_jacobian(finite_ring, z0, x0, y0)
assert not residue

rank = jacobian.rank()
print(
    f"Q8S3HENSEL|prime={prime}|jacobian_rank={rank}|"
    f"dimensions={jacobian.dimensions()}",
    flush=True,
)

if rank != 156:
    print(
        f"Q8S3HENSEL_RESULT|rank={rank}|expected=156|"
        "status=SINGULAR_MODULAR_SEED",
        flush=True,
    )
    raise SystemExit(2)

# Pivot columns of J^T are independent rows of J.
rows = tuple(int(row) for row in jacobian.transpose().pivots())
assert len(rows) == 156
square = jacobian.matrix_from_rows(rows)
assert square.nrows() == square.ncols() == 156
assert square.is_invertible()

print(
    f"Q8S3HENSEL|selected_rows={len(rows)}|"
    "status=PASS_NONSINGULAR_MODULAR_SEED",
    flush=True,
)

if args.rank_only:
    print("Q8S3HENSEL_RESULT|status=PASS_RANK_ONLY", flush=True)
    raise SystemExit(0)

target = int(args.precision)
if target < 2:
    raise ValueError("precision must be at least 2")

padic = Zp(prime, prec=target)
padic_ring = PolynomialRing(padic, "U")

# Start from integer representatives of the finite-field solution.
coefficients = [
    padic(ZZ(v))
    for v in Zraw[:24] + Xraw + Yraw
]
assert len(coefficients) == 156

seed_precision = 1
if args.seed:
    seed_path = args.seed.resolve()
    seed = json.loads(seed_path.read_text())
    if (
        seed.get("schema") != "elkies-k3.h92-q8-s3-direct-hensel-lift.v1"
        or ZZ(seed.get("prime")) != prime
        or "residues" not in seed
    ):
        raise ValueError("seed is not a compatible q8 S3 direct Hensel artifact")
    seed_precision = int(seed["precision"])
    if seed_precision >= target:
        raise ValueError("seed precision must be smaller than target precision")
    if len(seed["residues"]) != 156:
        raise ValueError("seed residue count is inconsistent")
    coefficients = [padic(ZZ(v)) for v in seed["residues"]]
    print(
        f"Q8S3HENSEL|seed_precision={seed_precision}|"
        f"target_precision={target}",
        flush=True,
    )


def valuation_floor(values):
    vals = [value.valuation() for value in values if value]
    return target if not vals else min(vals)


for iteration in range(1, 2*target + 3):
    residual, derivative = identity_and_jacobian(
        padic_ring,
        coefficients[:24],
        coefficients[24:77],
        coefficients[77:],
    )
    valuation = valuation_floor(residual)
    print(
        f"Q8S3HENSEL|iteration={iteration}|"
        f"residual_valuation={valuation}",
        flush=True,
    )
    if valuation >= target:
        break

    selected = derivative.matrix_from_rows(rows)
    correction = selected.solve_right(
        -vector(padic, [residual[row] for row in rows])
    )
    coefficients = [
        value + delta
        for value, delta in zip(coefficients, correction, strict=True)
    ]
else:
    raise ArithmeticError("p-adic Newton iteration did not reach target precision")

# Verify all 157 equations, not only the selected minor.
final_residual, unused = identity_and_jacobian(
    padic_ring,
    coefficients[:24],
    coefficients[24:77],
    coefficients[77:],
)
final_valuation = valuation_floor(final_residual)
assert final_valuation >= target

modulus = prime**target


def reconstruct(values):
    out = []
    for value in values:
        try:
            out.append(ZZ(value.lift()).rational_reconstruction(modulus))
        except (ArithmeticError, ValueError):
            out.append(None)
    return out


Zq = reconstruct(coefficients[:24]) + [QQ(1)]
Xq = reconstruct(coefficients[24:77])
Yq = reconstruct(coefficients[77:])

complete = not any(v is None for part in (Zq, Xq, Yq) for v in part)
exact_verified = False
exact_payload = None

if complete:
    RQ = PolynomialRing(QQ, "U")
    UQ = RQ.gen()
    Z = RQ(Zq)
    X = RQ(Xq)
    Y = RQ(Yq)
    A = RQ(A_QQ)
    B = RQ(B_QQ)

    exact_identity = Y**2 - X**3 - A*X*Z**4 - B*Z**6
    exact_verified = not exact_identity

    print(
        f"Q8S3HENSEL_EXACT|complete=1|"
        f"identity={int(exact_verified)}|"
        f"degrees={X.degree()}/{(Z**2).degree()},{Y.degree()}/{(Z**3).degree()}|"
        f"status={'PASS' if exact_verified else 'CANDIDATE_NOT_EXACT'}",
        flush=True,
    )

    if exact_verified:
        assert Z.degree() == 24 and Z.leading_coefficient() == 1
        assert X.degree() == 52
        assert Y.degree() == 78
        assert (Z**2).degree() == 48
        assert (Z**3).degree() == 72

        exact_payload = {
            "schema": "elkies-k3.h92-q8-s3-direct-section-qq.v1",
            "status": "PASS_EXACT_Q8_S3_DIRECT_SECTION",
            "zero": "II*_E8_1_branch_anchor",
            "profile": {
                "P_dot_O_from_denominator": 24,
                "Z_degree": 24,
                "X_degree": 52,
                "Y_degree": 78,
                "x_degrees": [52, 48],
                "y_degrees": [78, 72],
            },
            "section": {
                "Z_coefficients_low_to_high": [str(v) for v in Z.list()],
                "X_coefficients_low_to_high": [str(v) for v in X.list()],
                "Y_coefficients_low_to_high": [str(v) for v in Y.list()],
                "x": "X/Z^2",
                "y": "Y/Z^3",
            },
            "verification": {
                "exact_weierstrass_identity": True,
                "equation": "Y^2=X^3+A*X*Z^4+B*Z^6",
            },
            "source_modular_seed": str(INPUT.relative_to(ROOT)),
            "child_model": str(CHILD.relative_to(ROOT)),
            "hensel_prime": int(prime),
            "hensel_precision": target,
            "boundary": (
                "This certifies the exact rational Weierstrass section in the "
                "II*_E8_1 branch-anchored q8 D13 equation. Its exact D13 MW "
                "coordinate/local component class and the q24 neighbour equation "
                "are separate downstream certificates. Y is fixed only up to "
                "global section inversion."
            ),
        }

        exact_output = (
            args.exact_output.resolve()
            if args.exact_output
            else LOCAL / "q8-s3-direct-section-qq.json"
        )
        exact_output.parent.mkdir(parents=True, exist_ok=True)
        exact_output.write_text(
            json.dumps(exact_payload, indent=2, sort_keys=True) + "\n"
        )
        print(f"EXACT_OUTPUT|{exact_output}", flush=True)


payload = {
    "schema": "elkies-k3.h92-q8-s3-direct-hensel-lift.v1",
    "status": "PASS_Q8_S3_DIRECT_HENSEL",
    "prime": int(prime),
    "precision": target,
    "seed_precision": int(seed_precision),
    "jacobian_rank": int(rank),
    "selected_rows": list(rows),
    "final_residual_valuation": int(final_valuation),
    "complete": bool(complete),
    "exact_verified": bool(exact_verified),
    "profile": {
        "Z_degree": 24,
        "X_degree": 52,
        "Y_degree": 78,
        "x_denominator_degree": 48,
        "y_denominator_degree": 72,
    },
    "Z": [None if v is None else str(v) for v in Zq],
    "X": [None if v is None else str(v) for v in Xq],
    "Y": [None if v is None else str(v) for v in Yq],
    "residues": [str(ZZ(value.lift())) for value in coefficients],
    "source_modular_seed": str(INPUT.relative_to(ROOT)),
    "child_model": str(CHILD.relative_to(ROOT)),
}

if args.output:
    output = args.output
    if not output.is_absolute():
        output = ROOT / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"OUTPUT|{output}", flush=True)

print(
    f"Q8S3HENSEL_RESULT|precision={target}|"
    f"complete={int(complete)}|exact={int(exact_verified)}|"
    f"final_residual_valuation={final_valuation}|status=PASS",
    flush=True,
)
