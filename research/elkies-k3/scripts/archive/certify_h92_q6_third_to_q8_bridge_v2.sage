#!/usr/bin/env sage -python
"""
Bind the newly exact q6 third section to the already-certified q8 pencil.

Input:
  artifacts/local/elkies-k3/q6-third-iistar-ivstar-gauge.json

The exact third section is stored in the sparse II*/IV* base coordinate s.
This script converts it back to the canonical q6 child T-coordinate, verifies
the original q6 Weierstrass equation exactly, reconstructs the two already-
certified corrected q8 pencil generators from their stored kernel polynomials,
and evaluates the q8 base function U on S3.

Expected exact bridge:
    q8 map degree of S3 = 52.

That is the equation-level hand-off required before taking the degree-52
Abel-Jacobi trace on the D13 child.

Run:
  sage -python ~/Downloads/certify_h92_q6_third_to_q8_bridge.sage
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ


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
    "--third",
    type=Path,
    help="exact rebased q6 third section artifact",
)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"

THIRD = (
    args.third.resolve()
    if args.third
    else LOCAL / "q6-third-iistar-ivstar-gauge.json"
)
CHILD = GEN / "elkies-k3-h92-q6-child-jacobian.json"
Q8 = GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json"
TARGET = LOCAL / "q8-target-component-nef.json"
OUTPUT = (
    args.output.resolve()
    if args.output and args.output.is_absolute()
    else ROOT / (
        args.output
        if args.output
        else Path("artifacts/local/elkies-k3/q6-third-to-q8-bridge.json")
    )
)

for path in (THIRD, CHILD, Q8, TARGET):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

third = json.loads(THIRD.read_text())
child = json.loads(CHILD.read_text())
q8 = json.loads(Q8.read_text())
target = json.loads(TARGET.read_text())

assert third["status"] == "PASS_EXACT_Q6_THIRD_REBASED"
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert q8["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert target["normalization"]["representative"] == "component-nef"

# ---------------------------------------------------------------------------
# Exact third section: sparse s-gauge -> canonical T-gauge.
# ---------------------------------------------------------------------------
SR = PolynomialRing(QQ, "s")
s = SR.gen()
SF = SR.fraction_field()

Zs = SR([QQ(v) for v in third["Z"]])
Xs = SR([QQ(v) for v in third["X"]])
Ys = SR([QQ(v) for v in third["Y"]])
assert (Zs.degree(), Xs.degree(), Ys.degree()) == (21, 46, 69)
xs = SF(Xs) / SF(Zs**2)
ys = SF(Ys) / SF(Zs**3)

base = third["base_change"]
# T(s)=(b+a*s)/(d+c*s); the exact rebase used
# x_s=(c*s+d)^4*x_T and y_s=(c*s+d)^6*y_T.
b, a = map(QQ, base["T_numerator"])
d, c = map(QQ, base["T_denominator"])
D_s = c*s + d

TR = PolynomialRing(QQ, "T")
T = TR.gen()
TF = TR.fraction_field()

Lii = TR(base["L_II"])
Liv = TR(base["L_IV"])
s_of_T = TF(Lii) / TF(Liv)

def subst_s_to_T(value):
    value = SF(value)
    return TF(value.numerator()(s_of_T)) / TF(value.denominator()(s_of_T))

D_T = TF(D_s(s_of_T))
x3 = subst_s_to_T(xs) / D_T**4
y3 = subst_s_to_T(ys) / D_T**6

model = child["minimal_short_weierstrass"]
A = TR([QQ(v) for v in model["A_coefficients_low_to_high"]])
B = TR([QQ(v) for v in model["B_coefficients_low_to_high"]])
Delta = TR([QQ(v) for v in model["Delta_coefficients_low_to_high"]])
assert y3**2 == x3**3 + TF(A)*x3 + TF(B)

def degree_map(value):
    value = TF(value)
    return max(
        TR(value.numerator()).degree(),
        TR(value.denominator()).degree(),
    )

print(
    "Q6THIRDQ8BRIDGE|stage=canonical_q6_section|"
    f"x_degree={degree_map(x3)}|y_degree={degree_map(y3)}|"
    "identity=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# Reconstruct corrected q8 marked section S and its regular q frame.
# All data below are already certified by the q8 child artifact; we replay
# only enough algebra to evaluate the pencil on S3.
# ---------------------------------------------------------------------------
E = EllipticCurve(TF, [0, 0, 0, TF(A), TF(B)])

mdata = q8["marking"]["section"]
def rational_from_data(data, np, dp):
    return TF(TR([QQ(v) for v in data[np]])) / TF(
        TR([QQ(v) for v in data[dp]])
    )

sx = rational_from_data(
    mdata,
    "x_numerator_coefficients_low_to_high",
    "x_denominator_coefficients_low_to_high",
)
sy = rational_from_data(
    mdata,
    "y_numerator_coefficients_low_to_high",
    "y_denominator_coefficients_low_to_high",
)
S = E(sx, sy)

nx, dx = TR(sx.numerator()), TR(sx.denominator())
ny, dy = TR(sy.numerator()), TR(sy.denominator())

def monic_power_root(value, exponent):
    root = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= factor.monic() ** (multiplicity // exponent)
    return root.monic()

h = monic_power_root(dx, 2)
assert h == monic_power_root(dy, 3)
assert h.degree() == 10

ii = TR(next(item for item in child["finite_fibres"] if item["kodaira"] == "II*")["factor"]).monic()
iv = TR(next(item for item in child["finite_fibres"] if item["kodaira"] == "IV*")["factor"]).monic()
M = (ii**2 * iv**2).monic()

# Corrected q-normalization; keep Dx.
normalizer = (ny * dx * (h*dy).inverse_mod(nx)).mod(nx)
assert (normalizer*h*dy - ny*dx) % nx == 0

p_fun = -sy/sx
alpha = -p_fun/TF(h) - TF(normalizer)/TF(nx)
beta = TF(T**2)/TF(h)
rho = (normalizer * nx.inverse_mod(M)).mod(M)

abase = alpha + TF(rho)
bbase = beta
tbase = TF(M)

kernel_polys = q8["rr"]["kernel_polynomials"]
assert len(kernel_polys) == 2

pairs = []
for entry in kernel_polys:
    sp = TR(entry["s"])
    tp = TR(entry["t"])
    Bcoef = TF(sp)/TF(h)
    Acoef = (
        -TF(sp)*p_fun/TF(h)
        - TF(sp)*TF(normalizer)/TF(nx)
        + TF(sp*rho)
        + TF(tp*M)
    )
    # Replay the stored regular-frame algebra carefully.
    #
    # Acoef is the actual chord A coefficient and also equals the regular
    # frame aa.  For B, however, the q8 compiler deliberately uses
    #
    #     Bcoef = s/h
    #     bb    = s*T^2/h = T^2*Bcoef
    #
    # so that bb has the required infinity regularity.  Do NOT identify
    # Bcoef with bb.
    assert Acoef == TF(sp)*abase + TF(tp)*tbase
    assert TF(sp)*bbase == TF(T**2)*Bcoef
    pairs.append((Acoef, Bcoef))

(A0, B0), (A1, B1) = pairs

# The corrected marked chord is m=(y+y(S))/(x-x(S)); it passes through -S.
m3 = (y3 + sy) / (x3 - sx)
U3 = (A1 + B1*m3) / (A0 + B0*m3)
U3 = TF(U3)

u_num = TR(U3.numerator())
u_den = TR(U3.denominator())
common = u_num.gcd(u_den)
assert common in QQ
q8_degree = max(u_num.degree(), u_den.degree())

print(
    "Q6THIRDQ8BRIDGE|"
    f"q8_degree={q8_degree}|"
    f"numerator_degree={u_num.degree()}|"
    f"denominator_degree={u_den.degree()}|"
    "status=PASS" if q8_degree == 52 else
    "status=UNEXPECTED_DEGREE",
    flush=True,
)
assert q8_degree == 52

# Exact section profile in canonical q6 coordinate for downstream use.
def rf_data(value):
    value = TF(value)
    return {
        "numerator_coefficients_low_to_high": [
            str(v) for v in TR(value.numerator()).list()
        ],
        "denominator_coefficients_low_to_high": [
            str(v) for v in TR(value.denominator()).list()
        ],
    }

payload = {
    "schema": "elkies-k3.h92-q6-third-to-q8-bridge.v1",
    "status": "PASS_EXACT_Q6_THIRD_TO_Q8_DEGREE52",
    "inputs": {
        "third_rebased": {
            "path": str(THIRD.relative_to(ROOT)),
            "sha256": hashlib.sha256(THIRD.read_bytes()).hexdigest(),
        },
        "q6_child": {
            "path": str(CHILD.relative_to(ROOT)),
            "sha256": hashlib.sha256(CHILD.read_bytes()).hexdigest(),
        },
        "q8_child": {
            "path": str(Q8.relative_to(ROOT)),
            "sha256": hashlib.sha256(Q8.read_bytes()).hexdigest(),
        },
        "q8_target": {
            "path": str(TARGET.relative_to(ROOT)),
            "sha256": hashlib.sha256(TARGET.read_bytes()).hexdigest(),
        },
    },
    "third_section_canonical_q6": {
        "x": rf_data(x3),
        "y": rf_data(y3),
        "weierstrass_identity": True,
    },
    "q8_parameter_on_third": {
        "degree": q8_degree,
        "numerator_coefficients_low_to_high": [str(v) for v in u_num.list()],
        "denominator_coefficients_low_to_high": [str(v) for v in u_den.list()],
    },
    "lattice_target": {
        "q6_MW_word": [0, 0, 1],
        "expected_D13_AJ": [0, -1, 1, 1],
        "expected_q8_degree": 52,
    },
    "boundary": (
        "This binds the exact characteristic-zero q6 third section to the "
        "certified q8 pencil and proves its degree-52 multisection profile. "
        "The D13 Abel-Jacobi trace remains the next equation-level gate."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q6THIRDQ8BRIDGE_RESULT|"
    "q8_degree=52|expected_D13_AJ=0,-1,1,1|"
    "status=PASS_EXACT_Q6_THIRD_TO_Q8_DEGREE52",
    flush=True,
)
