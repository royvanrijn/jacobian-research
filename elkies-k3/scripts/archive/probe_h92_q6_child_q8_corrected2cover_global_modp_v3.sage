#!/usr/bin/env sage -python
"""
Corrected H92 q6-child q8 global RR probe after the binary-quartic 2-cover fix.

This script starts from the corrected standard-Jacobian marked section S of
height 24 and O-intersection 10.  It does NOT use the stale degree-46 marking
or any source-side 1600/1278 q8 compiler.

Geometry used:
  * generic divisor: O + S on the q6 child;
  * smooth collision divisor h has degree 10;
  * q=(m-p)/h with m=(y+yS)/(x-xS), p=-yS/xS is a regular frame at h;
  * the component-nef E8 and E6 cycles are the negatives of the x-valuation
    cycles, so the complete finite ideals are (u^2,X,Y) at II* and IV*;
  * therefore q maps to zero in both finite quotients and a local pair C+B*q
    is allowed exactly when C is divisible by
        M=f_II^2*f_IV^2;
  * normalize q globally by qreg=q-R/Nx, where
        R*h*Dy == Ny (mod Nx),
    removing q's generic base poles;
  * in qreg coordinates the finite module is
        < qreg + rho, M >,
    rho = R/Nx mod M;
  * place the global fibre coefficient -2 at infinity, so both coefficients
    in the smooth infinity frame <1,m_inf> must have order >=2.

The infinity bounds are derived, not guessed.  The resulting small complete
polynomial ambient is intersected modulo a good prime.  With --screen-levels,
a two-dimensional kernel is also tested by branch squareclass.

Run:
  sage -python ~/Downloads/probe_h92_q6_child_q8_corrected2cover_global_modp.sage --prime 43 --screen-levels
  sage -python ~/Downloads/probe_h92_q6_child_q8_corrected2cover_global_modp.sage --prime 59 --screen-levels
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix


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
            (candidate / "elkies-k3" / "scripts").is_dir()
            and (candidate / "artifacts" / "generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


def coeff(F, value):
    value = QQ(value)
    den = F(ZZ(value.denominator()))
    if not den:
        raise ValueError("chosen prime divides an input denominator")
    return F(ZZ(value.numerator())) / den


def poly(R, F, values):
    return R([coeff(F, value) for value in values])


def degree_or_minus_one(value):
    return -1 if not value else int(value.degree())


def infinity_order(value, R):
    if not value:
        return 10**9
    return int(R(value.denominator()).degree() - R(value.numerator()).degree())


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--screen-levels", action="store_true")
parser.add_argument("--corrected-marking", type=Path)
parser.add_argument("--component-target", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

if not ZZ(args.prime).is_prime() or args.prime in (2, 3):
    raise ValueError("prime must be an odd prime different from 3")

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts" / "generated-results"
LOCAL = ROOT / "artifacts" / "local" / "elkies-k3"

CHILD = GEN / "elkies-k3-h92-q6-child-jacobian.json"
CORRECTED = (
    args.corrected_marking.resolve()
    if args.corrected_marking
    else LOCAL / "q8-marking-2cover-corrected.json"
)
TARGET = (
    args.component_target.resolve()
    if args.component_target
    else LOCAL / "q8-target-component-nef.json"
)
OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / f"q8-corrected2cover-global-mod-{args.prime}.json"
)

for path in (CHILD, CORRECTED, TARGET):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

child = json.loads(CHILD.read_text())
corrected = json.loads(CORRECTED.read_text())
target = json.loads(TARGET.read_text())

assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert corrected["status"] == "PASS_EXACT_Q8_MARKING_2COVER_CORRECTION"
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert target["normalization"]["representative"] == "component-nef"

# Recheck the component-nef target that is actually being compiled.
selected = target["selected_q8"]
assert int(selected["vertical_fibre_coefficient"]) == -2
assert list(map(int, selected["E8"]["vertical_cycle"])) == [-4,-5,-7,-10,-8,-6,-4,-2]
assert list(map(int, selected["E8"]["component_degrees"])) == [1,0,0,0,0,0,0,0]
assert int(selected["E8"]["affine_component_degree"]) == 0
assert list(map(int, selected["E6"]["vertical_cycle"])) == [-2,-3,-4,-3,-2,-2]
assert list(map(int, selected["E6"]["component_degrees"])) == [1,0,0,0,1,0]
assert int(selected["E6"]["affine_component_degree"]) == 0

# In actual resolved chart order, -cycle is exactly the x-valuation cycle.
II_MAP = (2,7,6,4,8,3,5,1)
IV_MAPS = {
    "minus": (2,4,6,5,3,1),
    "plus": (3,5,6,4,2,1),
}
II_X = (2,4,6,10,4,7,5,8)
IV_X = (2,2,2,3,3,4)

def chart_cycle(physical, mapping):
    result = [None]*len(mapping)
    for i,j in enumerate(mapping):
        result[j-1] = int(physical[i])
    return tuple(result)

ii_cycle = chart_cycle(selected["E8"]["vertical_cycle"], II_MAP)
assert tuple(-v for v in ii_cycle) == II_X
for mapping in IV_MAPS.values():
    iv_cycle = chart_cycle(selected["E6"]["vertical_cycle"], mapping)
    assert tuple(-v for v in iv_cycle) == IV_X

F = GF(args.prime)
R = PolynomialRing(F, "T")
T = R.gen()
K = R.fraction_field()

model = child["minimal_short_weierstrass"]
A = poly(R, F, model["A_coefficients_low_to_high"])
Bcurve = poly(R, F, model["B_coefficients_low_to_high"])
Delta = poly(R, F, model["Delta_coefficients_low_to_high"])
assert A.degree() <= 8 and Bcurve.degree() <= 12 and Delta.degree() <= 24

sdata = corrected["selected_q8"][
    "relative_child_section_standard_jacobian_coordinates"
]
nx = poly(R, F, sdata["x_numerator_coefficients_low_to_high"])
dx = poly(R, F, sdata["x_denominator_coefficients_low_to_high"])
ny = poly(R, F, sdata["y_numerator_coefficients_low_to_high"])
dy = poly(R, F, sdata["y_denominator_coefficients_low_to_high"])
sx, sy = K(nx)/K(dx), K(ny)/K(dy)
assert sy**2 == sx**3 + K(A)*sx + K(Bcurve)

h = poly(
    R, F,
    corrected["selected_q8"]["collision_divisor"]["coefficients_low_to_high"],
).monic()
assert h.degree() == 10
qx, rx = dx.quo_rem(h**2)
qy, ry = dy.quo_rem(h**3)
assert not rx and not ry and qx in F and qy in F
assert nx.gcd(h).degree() == 0 and ny.gcd(h).degree() == 0
assert h.gcd(Delta).degree() == 0

# Exact finite additive factors.
srcR = PolynomialRing(QQ, "T")
ii = poly(
    R, F,
    srcR(next(item for item in child["finite_fibres"]
              if item["kodaira"] == "II*")["factor"]).list(),
).monic()
iv = poly(
    R, F,
    srcR(next(item for item in child["finite_fibres"]
              if item["kodaira"] == "IV*")["factor"]).list(),
).monic()
assert ii.degree() == iv.degree() == 1 and ii.gcd(iv).degree() == 0
assert h.gcd(ii*iv).degree() == 0
M = (ii**2 * iv**2).monic()
assert M.degree() == 4
assert nx.gcd(M).degree() == 0

# q=(m-p)/h.  Cancel its generic Nx pole by qreg=q-Rnorm/Nx.
#
# IMPORTANT: m=(y+Ny/Dy)/(x-Nx/Dx), so
#
#   -p/h = yS/(h*xS) = Ny*Dx/(h*Dy*Nx).
#
# Therefore the Nx-residue of q is Ny*Dx/(h*Dy), NOT Ny/(h*Dy).
# The missing Dx factor in the old q_regular normalization leaves a degree-24
# base pole and produces the spurious branch-degree ~100 pencil.
normalizer = (ny * dx * (h*dy).inverse_mod(nx)).mod(nx)
assert (normalizer*h*dy - ny*dx) % nx == 0
assert normalizer.degree() < nx.degree()

p_fun = -sy/sx
alpha = -p_fun/K(h) - K(normalizer)/K(nx)
beta = K(T**2)/K(h)

# Regression guard: after the corrected normalization, the base coefficient
# alpha has no Nx pole.  This is the exact condition that the previous
# normalization failed.
assert R(alpha.denominator()).gcd(nx).degree() == 0

# At II*/IV*, q maps to zero in R/(u^2,X,Y), so qreg maps to
# -Rnorm/Nx.  The finite module is C+B*qreg with
# C-B*Rnorm/Nx == 0 mod M.
rho = (normalizer * nx.inverse_mod(M)).mod(M)
assert rho.degree() < 4

# Verify the corrected degree profile.  These are consequences of
# O.S=10 on the globally minimal K3.
assert (nx.degree()-dx.degree(), ny.degree()-dy.degree()) == (4,6)
q_inf_order = int(h.degree()-2)
assert q_inf_order == 8
r_over_nx_order = int(nx.degree()-normalizer.degree())
assert r_over_nx_order >= 1

# qreg = alpha + beta*m_inf.  Finite generator:
#   g1 = qreg + rho
#   g2 = M.
abase = alpha + K(rho)
bbase = beta
tbase = K(M)

abase_order = infinity_order(abase, R)
bbase_order = infinity_order(bbase, R)
tbase_order = infinity_order(tbase, R)
assert bbase_order == 8
assert tbase_order == -4

required = 2
max_s = int(bbase_order - required)
assert max_s == 6

# A t term more polar than every possible s*abase term cannot cancel.
most_polar_s_a = int(abase_order - max_s)
max_t = int(tbase_order - most_polar_s_a)
assert max_t >= 0

labels = [("s", d) for d in range(max_s+1)]
labels += [("t", d) for d in range(max_t+1)]

a_cols = []
b_cols = []
for kind, d in labels:
    if kind == "s":
        a_cols.append(K(T**d)*abase)
        b_cols.append(K(T**d)*bbase)
    else:
        a_cols.append(K(T**d)*tbase)
        b_cols.append(K.zero())

def common_den(values):
    out = R.one()
    for value in values:
        if not value:
            continue
        out = out.lcm(R(value.denominator()))
    return out.monic()

def principal_rows(values):
    den = common_den(values)
    nums = [R(v*den) if v else R.zero() for v in values]
    cutoff = int(den.degree()-required)
    top = max((degree_or_minus_one(v) for v in nums), default=-1)
    rows = [
        [
            value[d] if value and d <= value.degree() else F.zero()
            for value in nums
        ]
        for d in range(cutoff+1, top+1)
    ]
    return rows, den, cutoff, top

a_rows, a_den, a_cutoff, a_top = principal_rows(a_cols)
b_rows, b_den, b_cutoff, b_top = principal_rows(b_cols)
condition = matrix(F, a_rows+b_rows, ncols=len(labels))
kernel = condition.right_kernel_matrix()

def st_from_row(row):
    s = R.zero()
    t = R.zero()
    for i,(kind,d) in enumerate(labels):
        if kind == "s":
            s += row[i]*T**d
        else:
            t += row[i]*T**d
    return R(s), R(t)

def pair_from_row(row):
    s,t = st_from_row(row)
    a = K(s)*abase + K(t)*tbase
    b = K(s)*bbase
    assert infinity_order(a,R) >= 2
    assert not b or infinity_order(b,R) >= 2
    # In original m coordinate:
    # qreg+rho = (m-p)/h - R/Nx + rho
    Acoef = -K(s)*p_fun/K(h) - K(s)*normalizer/K(nx) + K(s*rho) + K(t*M)
    Bcoef = K(s)/K(h)
    assert a == Acoef
    assert b == K(T**2)*Bcoef
    return s,t,Acoef,Bcoef,a,b

kernel_records = []
for row in kernel.rows():
    s,t,Acoef,Bcoef,a,b = pair_from_row(row)
    kernel_records.append({
        "s_degree": degree_or_minus_one(s),
        "t_degree": degree_or_minus_one(t),
        "a_infinity_order": infinity_order(a,R),
        "b_infinity_order": None if not b else infinity_order(b,R),
        "s_coefficients_low_to_high": [int(v) for v in s.list()],
        "t_coefficients_low_to_high": [int(v) for v in t.list()],
    })

branch = None
if args.screen_levels and kernel.nrows() == 2:
    _,_,A0,B0,_,_ = pair_from_row(kernel.row(0))
    _,_,A1,B1,_,_ = pair_from_row(kernel.row(1))
    xR = PolynomialRing(K, "x")
    x = xR.gen()
    histogram = {}
    good = []
    singular = []
    records = []
    for level in F:
        den = B1-K(level)*B0
        if not den:
            singular.append(int(level))
            continue
        mvalue = -(A1-K(level)*A0)/den
        yline = xR(mvalue)*(x-xR(sx))-xR(sy)
        rel = yline**2 - x**3 - xR(A)*x - xR(Bcurve)
        quad, rem = rel.quo_rem(x-xR(sx))
        assert not rem and quad.degree() == 2
        disc = xR.base_ring()(quad[1]**2 - 4*quad[2]*quad[0])
        num, denp = R(disc.numerator()), R(disc.denominator())
        odd_degree = 0
        odd_factors = []
        for side, value in (("n",num),("d",denp)):
            for fac,mult in value.squarefree_decomposition():
                if mult % 2:
                    odd_degree += fac.degree()
                    odd_factors.append([side,int(fac.degree()),int(mult)])
        inf_branch = (denp.degree()-num.degree()) % 2
        degree = int(odd_degree+inf_branch)
        histogram[degree] = histogram.get(degree,0)+1
        if degree == 4:
            good.append(int(level))
        records.append({
            "level": int(level),
            "branch_degree": degree,
            "disc_num_degree": int(num.degree()),
            "disc_den_degree": int(denp.degree()),
            "odd_factor_degrees": odd_factors,
            "infinity_branch": int(inf_branch),
        })
    branch = {
        "histogram": {str(k):v for k,v in sorted(histogram.items())},
        "genus_one_levels": good,
        "singular_levels": singular,
        "records": records,
    }

print(
    "Q8CORR2COVER_FRAME|"
    f"prime={args.prime}|h_deg={h.degree()}|Nx_deg={nx.degree()}|"
    f"Ny_deg={ny.degree()}|R_deg={normalizer.degree()}|rho_deg={rho.degree()}|"
    f"q_inf={q_inf_order}|qreg_alpha_inf={infinity_order(alpha,R)}|"
    f"g1_a_inf={abase_order}|g1_b_inf={bbase_order}|M_inf={tbase_order}",
    flush=True,
)
print(
    "Q8CORR2COVER_GLOBAL|"
    f"prime={args.prime}|s_max={max_s}|t_max={max_t}|"
    f"ambient={len(labels)}|a_rows={len(a_rows)}|b_rows={len(b_rows)}|"
    f"rank={condition.rank()}|kernel={kernel.nrows()}|"
    f"degrees={';'.join('s{}_t{}_a{}_b{}'.format(r['s_degree'],r['t_degree'],r['a_infinity_order'],'inf' if r['b_infinity_order'] is None else r['b_infinity_order']) for r in kernel_records) or 'none'}",
    flush=True,
)
if branch is not None:
    print(
        "Q8CORR2COVER_BRANCH|"
        f"prime={args.prime}|hist={','.join('{}:{}'.format(k,v) for k,v in branch['histogram'].items())}|"
        f"good={','.join(map(str,branch['genus_one_levels'])) or 'none'}|"
        f"singular={','.join(map(str,branch['singular_levels'])) or 'none'}",
        flush=True,
    )

def jsonable(value):
    """Recursively convert Sage scalar/container values to JSON-native types."""
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    # Sage Integer / Rational / finite-field integer-like values.
    try:
        parent = value.parent()
        if parent is ZZ:
            return int(value)
        if parent is QQ:
            return str(value)
    except Exception:
        pass
    return value


status = (
    "PASS_EXPECTED_Q8_MODULAR_PENCIL"
    if kernel.nrows() == 2 and (
        branch is None
        or bool(branch["genus_one_levels"])
    )
    else "DIAGNOSTIC_UNEXPECTED_DIMENSION_OR_BRANCH"
)

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-corrected2cover-global-modp.v1",
    "status": status,
    "prime": int(args.prime),
    "marking": {
        "height": corrected["selected_q8"]["height"],
        "O_intersection": corrected["selected_q8"]["O_intersection"],
        "h_degree": int(h.degree()),
    },
    "component_nef_target": {
        "vertical_fibre_coefficient": -2,
        "finite_ideal_II": "(u^2,X,Y)",
        "finite_ideal_IV": "(u^2,X,Y)",
        "finite_reason": "-physical cycle equals x valuation cycle",
    },
    "q_regular": {
        "normalizer_identity": "R*h*Dy == Ny*Dx mod Nx",
        "Nx_pole_cancelled": True,
        "normalizer_degree": int(normalizer.degree()),
        "rho_degree": int(rho.degree()),
        "q_infinity_order": q_inf_order,
        "alpha_infinity_order": infinity_order(alpha,R),
        "beta_infinity_order": infinity_order(beta,R),
    },
    "global_intersection": {
        "required_infinity_order": required,
        "s_degree_bound": max_s,
        "t_degree_bound": max_t,
        "ambient_dimension": len(labels),
        "a_rows": len(a_rows),
        "b_rows": len(b_rows),
        "rank": int(condition.rank()),
        "kernel_dimension": int(kernel.nrows()),
        "kernel": kernel_records,
    },
    "branch_screen": branch,
    "boundary": (
        "This is a modular complete-base intersection for the corrected "
        "degree-10 q8 marking and the exact component-nef finite valuation "
        "target. Characteristic-zero kernel reconstruction and elimination "
        "remain separate even when kernel=2."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(jsonable(payload), indent=2, sort_keys=True) + "\n")
print(f"Q8CORR2COVER_RESULT|status={status}")
print(f"OUTPUT|{OUTPUT}")
