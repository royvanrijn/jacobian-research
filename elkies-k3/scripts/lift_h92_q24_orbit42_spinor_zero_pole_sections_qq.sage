#!/usr/bin/env sage -python
"""Lift the two missing orbit42 spinor-class zero-pole sections over QQ.

The exact R3-zero D12 model already has eighteen certified rational
identity-class zero-pole sections.  The remaining abstract height shell has
two opposite spinor-class vectors.  Their pinned mod-100003 residue lies on
the cancellation branch

    deg(x) = 3,  deg(y) = 0.

This script solves that branch directly over QQ.  It is an exact section
certificate and a construction aid for the later resolved RR compiler; it
does not construct the orbit42 pencil or the A11 child.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


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
PROFILE = LOCAL / "q24-orbit42-spinor-zero-profiles.json"
OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q24-orbit42-spinor-zero-pole-sections-qq.json"
)

for path in (MODEL, PROFILE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

model_artifact = json.loads(MODEL.read_text())
profile_artifact = json.loads(PROFILE.read_text())
if model_artifact.get("status") != "PASS_Q42_ZERO_POLE_SMALLPRIME_SEEDS":
    raise ArithmeticError("zero-pole model prerequisite is not passing")
if profile_artifact.get("status") != "PASS_Q24_ORBIT42_EXACT_SPINOR_ZERO_PROFILES":
    raise ArithmeticError("spinor-profile prerequisite is not passing")

exact_model = model_artifact["exact_model"]
R = PolynomialRing(QQ, "u")
u = R.gen()
A = R([QQ(value) for value in exact_model["A_coefficients_low_to_high"]])
B = R([QQ(value) for value in exact_model["B_coefficients_low_to_high"]])
if (A.degree(), B.degree()) != (6, 9):
    raise ArithmeticError(f"unexpected exact D12 degrees {(A.degree(), B.degree())}")

p = ZZ(args.prime)
F = GF(p)


def reduce_q(value):
    value = QQ(value)
    denominator = ZZ(value.denominator())
    if denominator % p == 0:
        raise ZeroDivisionError(f"denominator divisible by regression prime {p}")
    return F(ZZ(value.numerator())) / F(denominator)


# The degree-nine coefficient must vanish.  The desired residue is the
# rational double root; the simple-root branch has no geometric section.
C = PolynomialRing(QQ, "c")
c = C.gen()
leading_equation = c**3 + A[6] * c + B[9]
roots = [
    (root, int(multiplicity))
    for root, multiplicity in leading_equation.roots()
    if reduce_q(root) == F(49445)
]
if len(roots) != 1 or roots[0][1] != 2:
    raise ArithmeticError(
        f"expected the unique pinned rational double root, found {roots}"
    )
c3 = QQ(roots[0][0])

# With y constant, coefficients u^1,...,u^8 give an overdetermined exact
# system in the remaining three x-coefficients.  Its unique QQ solution is
# replayed against every coefficient below.
S = PolynomialRing(QQ, names=("c2", "c1", "c0"), order="degrevlex")
c2, c1, c0 = S.gens()
K = S.fraction_field()
U = PolynomialRing(K, "z")
z = U.gen()
AA = U([K(value) for value in A.list()])
BB = U([K(value) for value in B.list()])
x_generic = K(c3) * z**3 + K(c2) * z**2 + K(c1) * z + K(c0)
rhs_generic = x_generic**3 + AA * x_generic + BB
equations = [S(rhs_generic[index]) for index in range(1, 9)]
ideal = S.ideal(equations)

print(
    "Q42SPINQQ|stage=SOLVE|vars=3|equations=8|status=GROEBNER_START",
    flush=True,
)
groebner_basis = ideal.groebner_basis()
if ideal.dimension() != 0 or len(groebner_basis) != 3:
    raise ArithmeticError(
        f"unexpected exact spinor ideal: dim={ideal.dimension()}, "
        f"gb={len(groebner_basis)}"
    )
solutions = ideal.variety()
if len(solutions) != 1:
    raise ArithmeticError(f"expected one exact spinor x-solution, got {len(solutions)}")
solution = solutions[0]
x_section = R([QQ(solution[c0]), QQ(solution[c1]), QQ(solution[c2]), c3])
rhs = R(x_section**3 + A * x_section + B)
if rhs.degree() > 0 or not QQ(rhs[0]).is_square():
    raise ArithmeticError("exact spinor rhs is not a rational constant square")
y_section = QQ(rhs[0]).sqrt()

expected_x_mod = [74701, 88005, 5652, 49445]
expected_y_mod = {45430, 54573}
x_mod = [int(reduce_q(value)) for value in x_section.list()]
y_mod = {int(reduce_q(y_section)), int(reduce_q(-y_section))}
if x_mod != expected_x_mod or y_mod != expected_y_mod:
    raise ArithmeticError(
        f"pinned residue mismatch: x={x_mod}, y={sorted(y_mod)}"
    )

sections = []
for sign in (1, -1):
    yy = QQ(sign) * y_section
    if R(yy**2) != x_section**3 + A * x_section + B:
        raise ArithmeticError("exact spinor Weierstrass identity failed")
    sections.append({
        "sign": sign,
        "x_coefficients_low_to_high": [str(value) for value in x_section.list()],
        "y_coefficients_low_to_high": [str(yy)],
        "degrees": [3, 0],
        "mod_100003": {
            "x_coefficients_low_to_high": x_mod,
            "y_coefficients_low_to_high": [int(reduce_q(yy))],
        },
    })

payload = {
    "schema": "elkies-k3.h3-q24-orbit42-spinor-zero-pole-sections-qq.v1",
    "status": "PASS_EXACT_Q42_SPINOR_ZERO_POLE_SECTIONS_QQ",
    "inputs": {
        "zero_pole_model": str(MODEL.relative_to(ROOT)),
        "zero_pole_model_sha256": hashlib.sha256(MODEL.read_bytes()).hexdigest(),
        "spinor_profiles": str(PROFILE.relative_to(ROOT)),
        "spinor_profiles_sha256": hashlib.sha256(PROFILE.read_bytes()).hexdigest(),
    },
    "model": {
        "A_degree": int(A.degree()),
        "B_degree": int(B.degree()),
        "zero": "R3",
    },
    "solve": {
        "ansatz": "deg(x)=3, deg(y)=0",
        "leading_root_multiplicity": 2,
        "equation_count": len(equations),
        "unknown_count": 3,
        "groebner_basis_length": len(groebner_basis),
        "solution_count": len(solutions),
    },
    "sections": sections,
    "verification": {
        "exact_weierstrass_identity": True,
        "opposite_pair": True,
        "mod_100003_regression": True,
    },
    "next_required": "Q42_RESOLVED_RR_TRIVIALIZATION",
    "proof_boundary": (
        "This certifies the two missing rational spinor-class zero-pole "
        "sections on the exact R3-zero D12 model. It does not select one of "
        "the surviving equation markings, construct H0(D42), compile a "
        "quartic, or certify the A11 child."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "Q42SPINQQ|stage=VERIFY|sections=2|xdeg=3|ydeg=0|"
    "identity=PASS|mod_100003=PASS|status=PASS",
    flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q42SPINQQ_RESULT|sections=2|next=Q42_RESOLVED_RR_TRIVIALIZATION|"
    "status=PASS_EXACT_Q42_SPINOR_ZERO_POLE_SECTIONS_QQ",
    flush=True,
)
