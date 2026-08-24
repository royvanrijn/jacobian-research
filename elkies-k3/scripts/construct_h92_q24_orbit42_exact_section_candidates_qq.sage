#!/usr/bin/env sage -python
"""Construct the four exact orbit42 section candidates on the R3-zero model.

The exact identity shell and the exact opposite spinor pair leave four
height-pairing-compatible extensions of the pinned mod-100003 shell marking.
For each extension this script performs the group law over QQ(u), changes the
zero from R3 to A0, projectivizes the result, and verifies the characteristic-
zero Weierstrass identity.  The four reductions are required to equal the
previous modular degree-three candidate set.

This is an exact equation-point construction with a modular marking boundary.
The resolved RR kernel must still select the physical orientation and certify
the A11 child.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents/jacobian-research",
        home / "jacobian-research",
        home / "src/jacobian-research",
        home / "git/jacobian-research",
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
        if (candidate / "elkies-k3/scripts").is_dir():
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=100003)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q24-orbit42-zero-pole-seeds-mod-43.json"
IDENTITY = LOCAL / "q24-orbit42-rational-zero-pole-sections-qq.json"
SPINOR = LOCAL / "q24-orbit42-spinor-zero-pole-sections-qq.json"
MATCHING = LOCAL / "q24-orbit42-identity-halving-qq.json"
OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q24-orbit42-exact-section-candidates-qq.json"
)

for path in (MODEL, IDENTITY, SPINOR, MATCHING):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

model_artifact = json.loads(MODEL.read_text())
identity_artifact = json.loads(IDENTITY.read_text())
spinor_artifact = json.loads(SPINOR.read_text())
matching_artifact = json.loads(MATCHING.read_text())
if identity_artifact.get("status") != "PASS_EXACT_Q42_RATIONAL_ZERO_POLE_SECTIONS_QQ":
    raise ArithmeticError("identity-shell prerequisite is not passing")
if spinor_artifact.get("status") != "PASS_EXACT_Q42_SPINOR_ZERO_POLE_SECTIONS_QQ":
    raise ArithmeticError("spinor prerequisite is not passing")
if matching_artifact.get("status") != "Q42_IDENTITY_HALVING_HAS_NO_A11_CHORD":
    raise ArithmeticError("pinned modular shell-matching artifact is not passing")

exact_model = model_artifact["exact_model"]
R = PolynomialRing(QQ, "u")
u = R.gen()
K = R.fraction_field()
A = R([QQ(value) for value in exact_model["A_coefficients_low_to_high"]])
B = R([QQ(value) for value in exact_model["B_coefficients_low_to_high"]])
E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])


def exact_point(row):
    x = R([QQ(value) for value in row["x_coefficients_low_to_high"]])
    y = R([QQ(value) for value in row["y_coefficients_low_to_high"]])
    if y**2 != x**3 + A*x + B:
        raise ArithmeticError("input exact point misses the exact D12 model")
    return E(K(x), K(y))


identity_points = [exact_point(row) for row in identity_artifact["sections"]]
spinor_points = [exact_point(row) for row in spinor_artifact["sections"]]
if len(identity_points) != 18 or len(spinor_points) != 2:
    raise ArithmeticError("expected an 18+2 exact zero-pole shell")
if spinor_points[1] != -spinor_points[0]:
    raise ArithmeticError("exact spinor points are not opposite")

matching = matching_artifact["matching"]
mappings = matching["mappings_abstract_to_equation"]
anchor = int(matching["abstract_anchor_index"])
if len(mappings) != 8:
    raise ArithmeticError("expected eight pinned identity-shell mappings")

# These are precisely the four extensions selected by the spinor-to-identity
# height pairings.  Abstractly, with the R3 zero,
#   target = -spinor + identity[1] - identity[5] + identity[8].
# Changing the zero to A0 subtracts twice the anchored A0-R3 section because
# the neighbour has old-fibre degree two.
extensions = ((0, 0), (1, 0), (6, 1), (7, 1))


def power_root(poly, exponent):
    poly = R(poly)
    if not poly or poly.degree() <= 0:
        return R.one()
    leading = QQ(poly.leading_coefficient())
    monic = R(poly / leading)
    answer = R.one()
    for factor, multiplicity in monic.factor():
        if int(multiplicity) % exponent:
            return None
        answer *= factor.monic()**(int(multiplicity) // exponent)
    return answer.monic()


def projectivize(point):
    x, y = point.xy()
    nx, dx = R(x.numerator()), R(x.denominator())
    ny, dy = R(y.numerator()), R(y.denominator())
    nx, dx = nx / dx.leading_coefficient(), dx / dx.leading_coefficient()
    ny, dy = ny / dy.leading_coefficient(), dy / dy.leading_coefficient()
    zx = power_root(dx, 2)
    zy = power_root(dy, 3)
    if zx is None or zy is None or zx != zy:
        raise ArithmeticError("exact target has incompatible section denominators")
    X = R(K(x) * K(zx**2))
    Y = R(K(y) * K(zx**3))
    if Y**2 != X**3 + A*X*zx**4 + B*zx**6:
        raise ArithmeticError("projectivized exact target identity failed")
    return X, Y, zx


p = ZZ(args.prime)
F = GF(p)


def reduce_q(value):
    value = QQ(value)
    denominator = ZZ(value.denominator())
    if denominator % p == 0:
        raise ZeroDivisionError(f"candidate denominator divisible by {p}")
    return F(ZZ(value.numerator())) / F(denominator)


def reduced_coefficients(poly):
    return [int(reduce_q(value)) for value in R(poly).list()]


expected_modular = {
    (
        tuple(map(int, row["X"])),
        tuple(map(int, row["Y"])),
        tuple(map(int, row["Z"])),
    )
    for row in matching_artifact["modular_halving"]["chord_results"]
}
if len(expected_modular) != 4:
    raise ArithmeticError("expected four pinned modular projective candidates")

candidate_records = []
actual_modular = set()
for mapping_index, spinor_index in extensions:
    mapping = mappings[mapping_index]
    target = (
        -spinor_points[spinor_index]
        + identity_points[mapping[1]]
        - identity_points[mapping[5]]
        + identity_points[mapping[8]]
        - 2*identity_points[mapping[anchor]]
    )
    X, Y, Z = projectivize(target)
    if Z.degree() != 3 or Z.leading_coefficient() != 1:
        raise ArithmeticError("exact orbit42 candidate does not have P.O=3")
    modular = (
        tuple(reduced_coefficients(X)),
        tuple(reduced_coefficients(Y)),
        tuple(reduced_coefficients(Z)),
    )
    actual_modular.add(modular)
    max_bits = max(
        max(abs(ZZ(QQ(value).numerator())).nbits(), abs(ZZ(QQ(value).denominator())).nbits())
        for poly in (X, Y, Z)
        for value in poly
    )
    candidate_records.append({
        "mapping_index": mapping_index,
        "spinor_index": spinor_index,
        "degrees": [int(X.degree()), int(Y.degree()), int(Z.degree())],
        "P_dot_O": 3,
        "max_coefficient_bits": int(max_bits),
        "X_coefficients_low_to_high": [str(value) for value in X.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y.list()],
        "Z_coefficients_low_to_high": [str(value) for value in Z.list()],
        "mod_100003": {
            "X_coefficients_low_to_high": list(modular[0]),
            "Y_coefficients_low_to_high": list(modular[1]),
            "Z_coefficients_low_to_high": list(modular[2]),
        },
    })
    print(
        "Q42CANDQQ|stage=TARGET|"
        f"mapping={mapping_index}|spinor={spinor_index}|"
        f"degrees={X.degree()},{Y.degree()},{Z.degree()}|"
        f"max_bits={max_bits}|identity=PASS|status=PASS",
        flush=True,
    )

if actual_modular != expected_modular:
    raise ArithmeticError("exact candidate set misses pinned mod-100003 census")

payload = {
    "schema": "elkies-k3.h3-q24-orbit42-exact-section-candidates-qq.v1",
    "status": "PASS_EXACT_Q42_ORBIT42_SECTION_CANDIDATES_QQ",
    "inputs": {
        "model": str(MODEL.relative_to(ROOT)),
        "identity_sections": str(IDENTITY.relative_to(ROOT)),
        "spinor_sections": str(SPINOR.relative_to(ROOT)),
        "modular_matching": str(MATCHING.relative_to(ROOT)),
        "input_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (MODEL, IDENTITY, SPINOR, MATCHING)
        },
    },
    "marking": {
        "boundary": "pinned_mod_100003_identity_shell_plus_exact_height_extensions",
        "extensions": [list(row) for row in extensions],
        "abstract_formula": "target_R3=-spinor+identity[1]-identity[5]+identity[8]",
        "zero_change": "target_A0=target_R3-2*(A0-R3)",
    },
    "candidates": candidate_records,
    "verification": {
        "candidate_count": len(candidate_records),
        "exact_weierstrass_identities": True,
        "all_P_dot_O": 3,
        "mod_100003_candidate_set_regression": True,
    },
    "next_required": "Q42_RESOLVED_RR_TRIVIALIZATION",
    "proof_boundary": (
        "The four displayed QQ(u) points and their projective Weierstrass "
        "identities are exact. Their identification as the four surviving "
        "orbit42 representatives uses the pinned mod-100003 shell marking. "
        "No resolved RR kernel, physical-orientation selection, quartic, or "
        "A11 child is claimed."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q42CANDQQ_RESULT|targets=4|PdotO=3|"
    "next=Q42_RESOLVED_RR_TRIVIALIZATION|"
    "status=PASS_EXACT_Q42_ORBIT42_SECTION_CANDIDATES_QQ",
    flush=True,
)
