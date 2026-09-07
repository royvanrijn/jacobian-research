#!/usr/bin/env sage -python
"""
Complete t-adic principal-part test on the omitted E7 identity component for
the corrected H92 q8 movable class, restricted to the current 7D survivor.

At the generic point of the identity component (the smooth locus of the
special cubic), t=0 and the old special fibre is

    y^2 = x^3,   m = y/x,

so the identity branch is x=m^2.  Hence the crude coefficient cap d<=64 is
too strong term-by-term: combinations x*m^j-m^(j+2) vanish normally.

This script solves the ACTUAL H92 Weierstrass equation on that branch in
K[[t]], K=GF(p)(m), using the exact marked chord

    y = y(P1) + m*(x-x(P1)),

with initial value x(0)=m^2.  It then computes the entire negative t-principal
part of each of the seven corrected survivors.  The common h_rev(t)^-16 is a
unit and may be omitted from the regularity test.

If this adds rank 5, the corrected RR space is two-dimensional.

Run:
  sage -python ~/Downloads/probe_h92_q8_corrected1278_e7_identity_full_jets.sage

Optional:
  --repo /path/to/jacobian-research
  --prime 43
"""

import argparse
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import (
    GF, PolynomialRing, PowerSeriesRing, QQ, ZZ, matrix
)


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
        home / "Documents" / "jacobian-research",
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
            (candidate / "elkies-k3" / "scripts").is_dir()
            and (candidate / "artifacts" / "generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=43)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts" / "generated-results"
P = int(args.prime)
if not ZZ(P).is_prime() or P in (2, 3):
    raise ValueError("prime must be odd and !=3")
finite = GF(P)

AMBIENT = GEN / "zz-h92-q8-corrected1278-ambient.json"
GLOBAL = GEN / f"zz-h92-q8-corrected1278-global-kernel-mod-{P}.json"
TRANSLATED = GEN / f"zz-h92-q8-corrected1278-two-translated-divisors-mod-{P}.json"
P1_PATH = GEN / "elkies-k3-h92-p1-lift.json"
H92 = ROOT / "artifacts" / "local" / "humbert-inputs" / "92" / "igusa92.txt"
ANCHOR = ROOT / "elkies-k3" / "scripts" / "verify_h3_noncm_q6_source_anchor.sage"
OUTPUT = GEN / f"zz-h92-q8-corrected1278-e7-identity-full-jets-mod-{P}.json"

for path in (AMBIENT, GLOBAL, TRANSLATED, P1_PATH, H92, ANCHOR):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

ambient = json.loads(AMBIENT.read_text())
global_kernel = json.loads(GLOBAL.read_text())
translated = json.loads(TRANSLATED.read_text())
p1 = json.loads(P1_PATH.read_text())

assert int(ambient["ambient_dimension"]) == 1278
assert int(global_kernel["dimensions"]["kernel"]) == 14
assert int(translated["combined"]["remaining_dimension"]) == 7
assert p1["status"] == "PASS_EXACT_H92_P1"

basis = ambient["ambient_basis"]

# Reconstruct the current 7D ambient coefficient space.
Kglobal = matrix(
    finite,
    [[finite(v) for v in row] for row in global_kernel["kernel_basis_rows"]],
)
trows = []
for record in translated["divisors"]:
    trows.extend([[finite(v) for v in row] for row in record["row_space_basis"]])
T = matrix(finite, trows).row_space().basis_matrix()
assert T.rank() == 7 and T.ncols() == 14
C7 = T.right_kernel().basis_matrix()
S7 = C7 * Kglobal
assert S7.nrows() == 7 and S7.ncols() == 1278 and S7.rank() == 7

# Coefficient field K=GF(p)(m).
MR = PolynomialRing(finite, "m")
m = MR.gen()
MF = MR.fraction_field()
mm = MF(m)

# Exact H92 surface coefficients.
anchor = SourceFileLoader(
    "h92_q8_identity_fulljets_anchor", str(ANCHOR)
).load_module()
r92, s92 = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)


def ff(value):
    value = QQ(value)
    den = finite(ZZ(value.denominator()))
    if not den:
        raise ValueError(f"prime {P} divides denominator {value}")
    return finite(ZZ(value.numerator())) / den


A1, A, B1, B, B2 = [ff(value(r92, s92)) for value in formulas]

# Reversed P1 functions in t=1/u.
TR = PolynomialRing(finite, "t")
tt = TR.gen()


def finite_poly(values):
    return TR([ff(v) for v in values])


def reversed_fraction(num_values, den_values):
    num = finite_poly(num_values)
    den = finite_poly(den_values)
    rn = sum(num[i] * tt**(num.degree()-i) for i in range(num.degree()+1))
    rd = sum(den[i] * tt**(den.degree()-i) for i in range(den.degree()+1))
    shift = den.degree() - num.degree()
    assert shift >= 0
    return tt**shift * rn, rd


xp_num, xp_den = reversed_fraction(
    p1["x_entrance_base"]["numerator_coefficients"],
    p1["x_entrance_base"]["denominator_coefficients"],
)
yp_num, yp_den = reversed_fraction(
    p1["y_entrance_base"]["numerator_coefficients"],
    p1["y_entrance_base"]["denominator_coefficients"],
)

# The ambient has d<=80, so the worst raw identity pole is t^-16.
M = max(max(0, int(e["actual_u_power"]) - 64) for e in basis)
assert M == 16
PREC = M + 24

PS = PowerSeriesRing(MF, "t", default_prec=PREC)
t = PS.gen()


def eval_t_poly(poly):
    out = PS.zero()
    for degree, coeff in poly.dict().items():
        out += MF(coeff) * t**int(degree)
    return PS(out)


xp_d = eval_t_poly(xp_den)
yp_d = eval_t_poly(yp_den)
assert xp_d[0] != 0 and yp_d[0] != 0
xp = PS(eval_t_poly(xp_num) * xp_d**(-1))
yp = PS(eval_t_poly(yp_num) * yp_d**(-1))

a_t = PS(MF(A1)*t**3 + MF(A)*t**4)
b_t = PS(
    MF(B1)*t**5 + MF(B)*t**6 + MF(B2)*t**7
)

# Newton/Hensel solve the identity branch x(0)=m^2.
x = PS(mm**2)
for _ in range(8):
    y = PS(yp + mm*(x-xp))
    F = PS(y**2 - x**3 - a_t*x - b_t)
    Fx = PS(2*mm*y - 3*x**2 - a_t)
    assert Fx[0] != 0
    x = PS(x - F * Fx**(-1))

y = PS(yp + mm*(x-xp))
check = PS(y**2 - x**3 - a_t*x - b_t)
assert all(check[i] == 0 for i in range(M + 12)), "Hensel precision insufficient"

z = PS(x - mm**2)
z_order = int(z.valuation())
print(
    f"E7IDENTITY_BRANCH|x0=m^2|z=x-m^2|z_order={z_order}|"
    f"precision={PREC}|status=PASS",
    flush=True,
)

# Powers needed by the generic basis.
x_powers = {0: PS.one(), 1: x}

# For each 7D survivor, multiply by t^M.  Terms of actual degree d contribute
# t^(M+64-d) x^a m^b.  Coefficients t^0..t^(M-1) are precisely the negative
# principal part before multiplication by t^M.
images = []
for survivor in S7.rows():
    image = PS.zero()
    for col, c in enumerate(survivor):
        if not c:
            continue
        entry = basis[col]
        d = int(entry["actual_u_power"])
        if d <= 64:
            # This starts at t^M or later and cannot affect the negative part.
            continue
        a = int(entry["x_power"])
        b = int(entry["m_power"])
        shift = M + 64 - d
        assert 0 <= shift < M
        image += MF(c) * t**shift * x_powers[a] * mm**b
    images.append(PS(sum(image[i]*t**i for i in range(M))))

# Convert MF-valued coefficients to honest GF(p)-linear rows.
rows = []
order_records = []
for t_order in range(M):
    coeffs = [MF(image[t_order]) for image in images]
    if not any(coeffs):
        continue
    common_den = MR.one()
    for value in coeffs:
        if value:
            common_den = common_den.lcm(MR(value.denominator()))
    polys = []
    max_degree = -1
    for value in coeffs:
        if not value:
            poly = MR.zero()
        else:
            den = MR(value.denominator())
            q, rem = common_den.quo_rem(den)
            assert not rem
            poly = MR(value.numerator()) * q
        polys.append(poly)
        if poly:
            max_degree = max(max_degree, int(poly.degree()))

    before = len(rows)
    for md in range(max_degree + 1):
        row = [finite(poly[md]) for poly in polys]
        if any(row):
            rows.append(row)

    block = matrix(finite, rows[before:]) if len(rows) > before else matrix(finite, 0, 7)
    order_records.append({
        "t_order_after_clearing": t_order,
        "laurent_order": t_order - M,
        "rows": len(rows)-before,
        "block_rank": int(block.rank()),
    })

Mmat = matrix(finite, rows) if rows else matrix(finite, 0, 7)
rank = int(Mmat.rank())
remaining = 7-rank
nonzero = sum(bool(image) for image in images)

C = Mmat.right_kernel().basis_matrix()
S = C * S7
assert S.nrows() == remaining and S.rank() == remaining

print(
    "Q8E7IDENTITYFULLJETS|"
    f"prime={P}|before=7|z_order={z_order}|rows={Mmat.nrows()}|"
    f"nonzero_images={nonzero}|rank={rank}|remaining={remaining}",
    flush=True,
)

for rec in order_records:
    print(
        "E7IDENTITY_JET|"
        f"laurent_order={rec['laurent_order']}|rows={rec['rows']}|"
        f"block_rank={rec['block_rank']}",
        flush=True,
    )

OUTPUT.write_text(json.dumps({
    "schema": "elkies-k3.h92-q8-corrected1278-e7-identity-full-jets-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_CORRECTED_Q8_E7_IDENTITY_FULL_JETS",
    "prime": P,
    "input_dimension": 7,
    "branch": {
        "special_fibre": "y^2=x^3",
        "identity_parameter_relation": "x=m^2 at t=0",
        "adapted_coordinate": "z=x-m^2",
        "z_normal_order": z_order,
        "hensel_precision": PREC,
    },
    "principal_part": {
        "worst_raw_order": -M,
        "rows": int(Mmat.nrows()),
        "nonzero_survivor_images": nonzero,
        "rank": rank,
        "remaining": remaining,
        "order_records": order_records,
    },
    "survivor_basis_rows": [[int(v) for v in row] for row in S.rows()],
    "boundary": (
        "This is the complete modular t-principal-part test at the generic "
        "point of the omitted E7 identity component, allowing exact "
        "x=m^2 cancellations. If remaining=2, repeat at another good prime "
        "and rebuild the global compiler with this component included."
    ),
}, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
