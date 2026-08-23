#!/usr/bin/env sage -python
"""
Derive a globally coherent q8 quartic-branch and D13 scaling orientation.

The modular q8 trace smoke test showed that choosing square roots independently
at each U=tau destroys global branch labels.  This exact characteristic-zero
script reconstructs the already-certified q8 binary quartic over QQ(U) and
extracts:

  1. w_IV(U), a rational square root of the quartic at the old IV* base point;
     +w_IV and -w_IV give globally coherent labels for the two IV* q8 sections.

  2. c2(U)=c(U)^2 and c3(U)=c(U)^3, the globally coherent scaling from the
     classical covariant Jacobian
         y^2 = x^3 + (-27 I)x + (-27 J)
     to the stored globally minimal D13 child.

These remove BOTH specialization-dependent signs from later sampling.

Run:
  sage -python ~/Downloads/derive_h92_q8_global_orientation.sage
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ


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
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
Q6 = GEN / "elkies-k3-h92-q6-child-jacobian.json"
Q8 = GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json"
OUTPUT = (
    args.output.resolve()
    if args.output and args.output.is_absolute()
    else ROOT / (
        args.output
        if args.output
        else Path("artifacts/local/elkies-k3/q8-global-orientation.json")
    )
)

for path in (CORE, Q6, Q8):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

scope = {}
exec(compile(CORE.read_text(), str(CORE), "exec"), scope)
squarefree_binary_quartic = scope["squarefree_binary_quartic"]
binary_quartic_invariants = scope["binary_quartic_invariants"]

q6 = json.loads(Q6.read_text())
q8 = json.loads(Q8.read_text())
assert q6["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert q8["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"

RT = PolynomialRing(QQ, "T")
T = RT.gen()
KT = RT.fraction_field()

def polynomial(values):
    return RT([QQ(value) for value in values])

def rational_from_data(data, np, dp):
    return KT(polynomial(data[np])) / KT(polynomial(data[dp]))

A6 = polynomial(q6["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B6 = polynomial(q6["minimal_short_weierstrass"]["B_coefficients_low_to_high"])

mdata = q8["marking"]["section"]
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
assert sy**2 == sx**3 + KT(A6)*sx + KT(B6)

def monic_power_root(value, exponent):
    result = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        result *= factor.monic()**(multiplicity//exponent)
    return result.monic()

nx, dx = RT(sx.numerator()), RT(sx.denominator())
ny, dy = RT(sy.numerator()), RT(sy.denominator())
h = monic_power_root(dx, 2)
assert h == monic_power_root(dy, 3)
assert h.degree() == 10

ii = RT(next(
    item for item in q6["finite_fibres"] if item["kodaira"] == "II*"
)["factor"]).monic()
iv = RT(next(
    item for item in q6["finite_fibres"] if item["kodaira"] == "IV*"
)["factor"]).monic()
assert ii.degree() == iv.degree() == 1
t_iv = -iv[0]/iv[1]
M = (ii**2 * iv**2).monic()

normalizer = (ny*dx*(h*dy).inverse_mod(nx)).mod(nx)
assert (normalizer*h*dy - ny*dx) % nx == 0
p_fun = -sy/sx
rho = (normalizer*nx.inverse_mod(M)).mod(M)

pairs = []
for entry in q8["rr"]["kernel_polynomials"]:
    sp = RT(entry["s"])
    tp = RT(entry["t"])
    Bcoef = KT(sp)/KT(h)
    Acoef = (
        -KT(sp)*p_fun/KT(h)
        - KT(sp)*KT(normalizer)/KT(nx)
        + KT(sp*rho)
        + KT(tp*M)
    )
    pairs.append((Acoef, Bcoef))
assert len(pairs) == 2
(A0, B0), (A1, B1) = pairs

# Rebuild the exact q8 quartic over QQ(U), using exactly the certified pencil.
print(
    "Q8GLOBALORIENT|stage=generic_rebuild_start|status=RUNNING",
    flush=True,
)
RU = PolynomialRing(QQ, "U")
U = RU.gen()
KU = RU.fraction_field()
TU = PolynomialRing(KU, "T")
TT = TU.gen()
KTU = TU.fraction_field()

def lift_poly(value):
    value = RT(value)
    return TU([KU(coefficient) for coefficient in value.list()])

def lift_rat(value):
    value = KT(value)
    return KTU(lift_poly(value.numerator())) / KTU(lift_poly(value.denominator()))

m_value = -(lift_rat(A1)-KU(U)*lift_rat(A0)) / (
    lift_rat(B1)-KU(U)*lift_rat(B0)
)
sxU, syU = lift_rat(sx), lift_rat(sy)
A6U, B6U = lift_poly(A6), lift_poly(B6)

XR = PolynomialRing(KTU, "x")
x = XR.gen()
y_line = XR(m_value)*(x-XR(sxU))-XR(syU)
relation = y_line**2 - x**3 - XR(A6U)*x - XR(B6U)
quadratic, remainder = relation.quo_rem(x-XR(sxU))
assert not remainder and quadratic.degree() == 2
disc = KTU(quadratic[1]**2 - 4*quadratic[2]*quadratic[0])
assert disc

quartic, square_factor = squarefree_binary_quartic(disc, TU)
assert quartic.degree() == 4
I, J = binary_quartic_invariants(quartic)
stdA = KU(-27*I)
stdB = KU(-27*J)
print(
    "Q8GLOBALORIENT|quartic_degree=4|stage=generic_quartic|status=PASS",
    flush=True,
)

# Stored globally minimal D13 model.
Amin = RU([
    QQ(value) for value in q8["child"]["minimal_A_coefficients_low_to_high"]
])
Bmin = RU([
    QQ(value) for value in q8["child"]["minimal_B_coefficients_low_to_high"]
])
AminK, BminK = KU(Amin), KU(Bmin)

# c2=c^2 is determined without extracting a root:
#   Amin=c^4*stdA, Bmin=c^6*stdB
# hence c^2=(Bmin/stdB)/(Amin/stdA).
c2 = KU(BminK*stdA/(stdB*AminK))
assert c2**2 == AminK/stdA
assert c2**3 == BminK/stdB

def qq_sqrt(value):
    value = QQ(value)
    if value < 0:
        return None
    num = ZZ(value.numerator())
    den = ZZ(value.denominator())
    if not num.is_square() or not den.is_square():
        return None
    return QQ(num.sqrt()) / QQ(den.sqrt())

def polynomial_square_root(poly):
    poly = RU(poly)
    if not poly:
        return RU.zero()
    factors = poly.factor()
    unit_root = qq_sqrt(QQ(factors.unit()))
    if unit_root is None:
        return None
    result = RU(unit_root)
    for factor, multiplicity in factors:
        if multiplicity % 2:
            return None
        result *= factor**(multiplicity//2)
    assert result**2 == poly
    return result

def rational_square_root(value):
    value = KU(value)
    numerator = RU(value.numerator())
    denominator = RU(value.denominator())
    nr = polynomial_square_root(numerator)
    dr = polynomial_square_root(denominator)
    if nr is None or dr is None or not dr:
        return None
    result = KU(nr)/KU(dr)
    assert result**2 == value
    # Canonical global sign: numerator leading coefficient positive after
    # fraction-field normalization.  This is only a label convention.
    if QQ(result.numerator().leading_coefficient()) < 0:
        result = -result
    return result

c3 = rational_square_root(c2**3)
if c3 is None:
    raise ArithmeticError("global D13 c^3 scaling is not rational")
assert c3**2 == c2**3

# At the old IV* base point the quartic has two rational QQ(U) points.
qiv = KU(quartic(t_iv))
wiv = rational_square_root(qiv)
if wiv is None:
    raise ArithmeticError(
        "q8 quartic at old IV* is not a square in QQ(U); "
        "the assumed rational IV* section interpretation must be revisited"
    )
assert wiv**2 == qiv

def rf_record(value):
    value = KU(value)
    return {
        "numerator_coefficients_low_to_high": [
            str(v) for v in RU(value.numerator()).list()
        ],
        "denominator_coefficients_low_to_high": [
            str(v) for v in RU(value.denominator()).list()
        ],
        "numerator_degree": int(RU(value.numerator()).degree()),
        "denominator_degree": int(RU(value.denominator()).degree()),
    }


def tu_polynomial_record(poly):
    poly = TU(poly)
    return {
        "degree": int(poly.degree()),
        "coefficients_low_to_high": [
            rf_record(poly[i]) for i in range(poly.degree() + 1)
        ],
    }


def ktu_record(value):
    value = KTU(value)
    return {
        "numerator": tu_polynomial_record(value.numerator()),
        "denominator": tu_polynomial_record(value.denominator()),
    }

print(
    "Q8GLOBALORIENT|"
    f"wIV={RU(wiv.numerator()).degree()}/{RU(wiv.denominator()).degree()}|"
    f"c2={RU(c2.numerator()).degree()}/{RU(c2.denominator()).degree()}|"
    f"c3={RU(c3.numerator()).degree()}/{RU(c3.denominator()).degree()}|"
    "squares=PASS|status=PASS",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-q8-global-orientation.v2",
    "status": "PASS_EXACT_Q8_GLOBAL_ORIENTATION",
    "inputs": {
        "q6_child": {
            "path": str(Q6.relative_to(ROOT)),
            "sha256": hashlib.sha256(Q6.read_bytes()).hexdigest(),
        },
        "q8_child": {
            "path": str(Q8.relative_to(ROOT)),
            "sha256": hashlib.sha256(Q8.read_bytes()).hexdigest(),
        },
    },
    "quartic": {
        "old_IVstar_T": str(t_iv),
        "quartic_at_old_IVstar_is_square": True,
        "global_w_plus": rf_record(wiv),
        "global_w_minus": rf_record(-wiv),
        "label_convention": (
            "plus is the rational square root whose normalized numerator "
            "has positive leading coefficient; minus is its negative"
        ),
    },
    "generic_transport": {
        "quartic_in_T_over_QQ_U": tu_polynomial_record(quartic),
        "square_factor_in_QQ_U_T": ktu_record(square_factor),
        "purpose": (
            "Serialize the exact generic q8 quartic normalization so modular "
            "samplers never rebuild it or choose specialization-local square "
            "normalizations."
        ),
    },
    "minimalization": {
        "c2_equals_c_squared": rf_record(c2),
        "c3_equals_c_cubed": rf_record(c3),
        "checks": [
            "c2^2 = A_min/A_covariant",
            "c2^3 = B_min/B_covariant",
            "c3^2 = c2^3",
        ],
    },
    "boundary": (
        "This fixes coherent algebraic signs over QQ(U).  It does not identify "
        "which IV* section equals the pinned lattice zero; that is resolved "
        "later from section transport/profile data."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q8GLOBALORIENT_RESULT|wIV=RATIONAL|c2=RATIONAL|c3=RATIONAL|"
    "status=PASS_EXACT_Q8_GLOBAL_ORIENTATION",
    flush=True,
)
