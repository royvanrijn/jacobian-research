#!/usr/bin/env sage -python
"""
Probe the translated third t-divisor visible on the actual H92
E7_2--E7_5 chart, restricted to the true 18-dimensional q8 survivor.

Motivation
----------
On this chart
    t = Z*U*(Z + c)
with c != 0 at the E7_2--E7_5 node.  The node-local calculation correctly
inverts Z+c, so it cannot see the divisor
    L := Z+c = 0.
This script tests regularity along that divisor directly.

It uses the ACTUAL global numerator degree d, not the q6^9 helper exponent
i=d-9.  After stripping the helper t^9, the reduced q8 condition along a
component on which the integral E7 twist g has valuation zero is simply the
regularity of the global q8 section.

For
    f = t^(72-d)/h_rev(t)^18 * x^a * m^b,
    m = (y-yP)/(x-xP),
write
    Lx = Nx/xP_den,  Ly = Ny/yP_den.
Multiplying every section by the common factor
    C = t^17 * h_rev^18 * Nx^9 * yP_den^9
gives polynomial columns
    N_i =
      t^(89-d) * x^a * Ny^b * Nx^(9-b)
      * xP_den^b * yP_den^(9-b).

Along L=0 on this chart, t,x,Nx,Ny each have one factor L generically,
so ord_L(C)=17+9=26.  Therefore f is regular along L=0 iff its common
cleared numerator is divisible by L^26 in the actual surface ring.

The surface is monic quadratic in Y.  We reduce exactly in the rank-two
module {1,Y} over GF(p)[L,U], truncating only L^26.  No Groebner basis is
computed.

Run:
  sage -python ~/Downloads/probe_h92_q8_true1600_translated_divisor.sage

Optional:
  --repo /path/to/jacobian-research
  --prime 43
  --chart E7_2--E7_5
"""

import argparse
import hashlib
import json
import time
from collections import defaultdict
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, sage_eval


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


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ff_coefficient(field, value):
    value = QQ(value)
    den = field(ZZ(value.denominator()))
    if not den:
        raise ValueError("prime divides an input coefficient denominator")
    return field(ZZ(value.numerator())) / den


def ff_polynomial(ring, field, values):
    return ring([ff_coefficient(field, value) for value in values])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--chart", default="E7_2--E7_5")
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts" / "generated-results"
P = int(args.prime)
if not ZZ(P).is_prime() or P in (2, 3):
    raise ValueError("prime must be an odd prime different from 3")
finite = GF(P)

P1_PATH = GEN / "elkies-k3-h92-p1-lift.json"
PULLBACKS_PATH = GEN / "elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
GLUING_PATH = GEN / "elkies-k3-h92-q8-actual-e7-gluing.json"
AMBIENT_PATH = GEN / "zz-h92-q8-true1600-ambient.json"
KERNEL_PATH = GEN / f"zz-h92-q8-true1600-global-kernel-mod-{P}.json"

for path in (P1_PATH, PULLBACKS_PATH, GLUING_PATH, AMBIENT_PATH, KERNEL_PATH):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

p1 = json.loads(P1_PATH.read_text())
pullbacks = json.loads(PULLBACKS_PATH.read_text())
gluing = json.loads(GLUING_PATH.read_text())
ambient = json.loads(AMBIENT_PATH.read_text())
kernel = json.loads(KERNEL_PATH.read_text())

assert p1["status"] == "PASS_EXACT_H92_P1"
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"
assert gluing["status"] == "PASS_EXACT_Q8_ACTUAL_E7_GLUING"
assert int(ambient["ambient_dimension"]) == 1600
assert int(kernel["prime"]) == P
assert int(kernel["dimensions"]["ambient"]) == 1600
assert int(kernel["dimensions"]["kernel"]) == 18

basis = ambient["ambient_basis"]
kernel_rows = [[finite(v) for v in row] for row in kernel["kernel_basis_rows"]]
assert len(basis) == 1600
assert len(kernel_rows) == 18
assert all(len(row) == 1600 for row in kernel_rows)

charts = {entry["name"]: entry for entry in pullbacks["charts"]}
edges = {entry["name"]: entry for entry in gluing["actual_edge_chart_gluing"]}
if args.chart not in charts:
    raise ValueError(f"unknown chart {args.chart}")
chart = charts[args.chart]
edge = edges[args.chart]

# This first probe is intentionally pinned to the translated E7_2--E7_5 chart.
# E7_3--E7_6 has another translated factor and can be added immediately if
# this one is active.
if args.chart != "E7_2--E7_5":
    raise SystemExit("This version is certified for --chart E7_2--E7_5 only.")

# Original chart over GF(p).
R = PolynomialRing(finite, names=("Z", "U", "Y"))
Z, U, Y = R.gens()

QQR = PolynomialRing(QQ, names=("Z", "U", "Y"))
Zq, Uq, Yq = QQR.gens()
loc = {"Z": Zq, "U": Uq, "Y": Yq}

surface = R(QQR(sage_eval(chart["surface_equation"], locals=loc)))
t_map = R(QQR(sage_eval(chart["old_coordinate_pullback"]["t"], locals=loc)))
x_map = R(QQR(sage_eval(chart["old_coordinate_pullback"]["x"], locals=loc)))
y_map = R(QQR(sage_eval(chart["old_coordinate_pullback"]["y"], locals=loc)))
g_map = R(QQR(sage_eval(edge["w_cartier_equation"], locals=loc)))

# Identify the node-local unit factor L in t.  The common monomial at the node
# is Z*U, and t/(Z*U) is the translated linear factor.
def common_monomial_exponents(value):
    terms = list(value.dict())
    assert terms
    return tuple(min(exp[i] for exp in terms) for i in range(3))

t_mono_exp = common_monomial_exponents(t_map)
assert t_mono_exp == (1, 1, 0), t_mono_exp
t_monomial = Z * U
L_old = R(t_map // t_monomial)
assert L_old.degree(Z) == 1
assert L_old.degree(U) == 0 and L_old.degree(Y) == 0
z_coeff = L_old.monomial_coefficient(Z)
constant = L_old(0, 0, 0)
assert z_coeff and constant
# L_old = z_coeff*Z + constant. Introduce new L equal to L_old.
assert L_old == z_coeff * Z + constant

# New chart coordinates (L,U,Y), with Z=(L-constant)/z_coeff.
S = PolynomialRing(finite, names=("L", "U", "Y"))
L, Us, Ys = S.gens()
Z_of_L = (L - finite(constant)) / finite(z_coeff)

def translate(poly):
    # Polynomial substitution preserves exact coefficients mod p.
    return S(poly(Z_of_L, Us, Ys))

surface_L = translate(surface)
t_L = translate(t_map)
x_L = translate(x_map)
y_L = translate(y_map)
g_L = translate(g_map)
assert t_L % L == 0
assert x_L % L == 0

# The translated divisor is not one of the two root components appearing in g:
# g must have L-valuation zero, i.e. it is nonzero modulo L.
assert S(g_L.subs({L: 0})) != 0

# P1 reversed functions over the E7 base t.
TBASE = PolynomialRing(finite, "t")
tt = TBASE.gen()

def reversed_fraction(numerator_values, denominator_values):
    numerator = ff_polynomial(TBASE, finite, numerator_values)
    denominator = ff_polynomial(TBASE, finite, denominator_values)
    rev_numerator = sum(
        numerator[i] * tt ** (numerator.degree() - i)
        for i in range(numerator.degree() + 1)
    )
    rev_denominator = sum(
        denominator[i] * tt ** (denominator.degree() - i)
        for i in range(denominator.degree() + 1)
    )
    shift = denominator.degree() - numerator.degree()
    assert shift >= 0
    return tt**shift * rev_numerator, rev_denominator

xp_num, xp_den = reversed_fraction(
    p1["x_entrance_base"]["numerator_coefficients"],
    p1["x_entrance_base"]["denominator_coefficients"],
)
yp_num, yp_den = reversed_fraction(
    p1["y_entrance_base"]["numerator_coefficients"],
    p1["y_entrance_base"]["denominator_coefficients"],
)

xp_n_L = S(xp_num(t_L))
xp_d_L = S(xp_den(t_L))
yp_n_L = S(yp_num(t_L))
yp_d_L = S(yp_den(t_L))
assert xp_d_L.subs({L: 0}) != 0
assert yp_d_L.subs({L: 0}) != 0

nx_L = S(x_L * xp_d_L - xp_n_L)
ny_L = S(y_L * yp_d_L - yp_n_L)

# Generic translated-divisor orders. Both differences carry exactly one L.
def l_valuation(poly):
    if not poly:
        return 10**9
    return min(int(exp[0]) for exp in poly.dict())

orders = {
    "t": l_valuation(t_L),
    "x": l_valuation(x_L),
    "y": l_valuation(y_L),
    "Nx": l_valuation(nx_L),
    "Ny": l_valuation(ny_L),
    "g": l_valuation(g_L),
}
assert orders == {"t": 1, "x": 1, "y": 1, "Nx": 1, "Ny": 1, "g": 0}, orders

# Divide the known L factors. These are regular elements of the translated ring.
t1 = S(t_L // L)
x1 = S(x_L // L)
nx1 = S(nx_L // L)
ny1 = S(ny_L // L)
assert t1.subs({L: 0}) != 0
assert x1.subs({L: 0}) != 0
assert nx1.subs({L: 0}) != 0
assert ny1.subs({L: 0}) != 0

# Surface is monic quadratic in Y:
#     Y^2 = H(L,U)
# Extract H by checking surface = Y^2 - H(L,U).
H_LU = S(Ys**2 - surface_L)
assert H_LU.degree(Ys) == 0
assert surface_L == Ys**2 - H_LU

# Rank-two coefficient ring A = GF(p)[L,U], basis {1,Y}.
A = PolynomialRing(finite, names=("L", "U"))
La, Ua = A.gens()
H = A(H_LU(La, Ua, 0))

MOD_L = 26

def trunc(poly):
    poly = A(poly)
    return A({
        exp: coeff
        for exp, coeff in poly.dict().items()
        if int(exp[0]) < MOD_L
    })

def add_pair(p, q):
    return (trunc(p[0] + q[0]), trunc(p[1] + q[1]))

def scale_pair(c, p):
    return (trunc(c * p[0]), trunc(c * p[1]))

def mul_pair(p, q):
    return (
        trunc(p[0] * q[0] + p[1] * q[1] * H),
        trunc(p[0] * q[1] + p[1] * q[0]),
    )

ONE = (A.one(), A.zero())
ZERO = (A.zero(), A.zero())

def pair_power(p, exponent):
    result = ONE
    base = p
    while exponent:
        if exponent & 1:
            result = mul_pair(result, base)
        exponent //= 2
        if exponent:
            base = mul_pair(base, base)
    return result

def to_pair(poly):
    poly = S(poly)
    even = A.zero()
    odd = A.zero()
    # General exact reduction using Y^(2q+r)=H^q Y^r.
    hp_cache = {0: A.one()}
    for (le, ue, ye), coeff in poly.dict().items():
        q, parity = divmod(int(ye), 2)
        if q not in hp_cache:
            hp_cache[q] = trunc(H**q)
        term = trunc(A({(int(le), int(ue)): coeff}) * hp_cache[q])
        if parity:
            odd = trunc(odd + term)
        else:
            even = trunc(even + term)
    return (even, odd)

# All base factors as rank-two elements.
t1p = to_pair(t1)
x1p = to_pair(x1)
nx1p = to_pair(nx1)
ny1p = to_pair(ny1)
xpd = to_pair(xp_d_L)
ypd = to_pair(yp_d_L)
assert xpd[1] == 0 and ypd[1] == 0 and t1p[1] == 0 and x1p[1] == 0 and nx1p[1] == 0

# Precompute family unit factors:
#   x1^a * ny1^b * nx1^(9-b) * xp_d^b * yp_d^(9-b).
family_factor = {}
for a in (0, 1):
    max_b = 9 if a == 0 else 7
    for b in range(max_b + 1):
        f = ONE
        for base, exponent in (
            (x1p, a),
            (ny1p, b),
            (nx1p, 9 - b),
            (xpd, b),
            (ypd, 9 - b),
        ):
            if exponent:
                f = mul_pair(f, pair_power(base, exponent))
        family_factor[(a, b)] = f

# t1 powers needed after extracting explicit L powers.
t1_powers = [ONE]
for e in range(1, 18):
    t1_powers.append(mul_pair(t1_powers[-1], t1p))

def actual_degree(entry):
    if "actual_u_power" in entry:
        return int(entry["actual_u_power"])
    return int(entry["u_power"]) + 9

# Organize only columns that can have a negative translated-divisor order.
active_by_family = defaultdict(list)
for index, entry in enumerate(basis):
    a = int(entry["x_power"])
    b = int(entry["m_power"])
    d = actual_degree(entry)
    k = int(entry["h_power"])
    assert k == 18
    # ord_L(f)=72-d+a. If >=0 this column is already regular.
    pole_order = d - 72 - a
    if pole_order <= 0:
        continue
    # In the common cleared numerator, explicit L shift is
    # 26 - pole_order = 98-d+a.
    shift = 26 - pole_order
    t_exp = 89 - d
    assert 0 <= t_exp <= 17
    assert shift == 9 + a + t_exp
    assert 0 <= shift < 26
    active_by_family[(a, b)].append((index, d, t_exp, shift, pole_order))

active_indices = sorted(index for rows in active_by_family.values() for index, *_ in rows)
print(
    "TRANSLATED_DIVISOR_START|"
    f"prime={P}|chart={args.chart}|L={L_old}|"
    f"orders=t:{orders['t']},x:{orders['x']},Nx:{orders['Nx']},Ny:{orders['Ny']},g:{orders['g']}|"
    f"threshold=L^{MOD_L}|active_columns={len(active_indices)}|survivors=18",
    flush=True,
)

# Build the 18 exact obstruction images modulo L^26.
images = []
build_start = time.perf_counter()

for survivor_index, row in enumerate(kernel_rows):
    image = ZERO
    kernel_terms = 0
    for family, entries in active_by_family.items():
        # Scalar t-polynomial piece after extracting the common family factor.
        family_sum = ZERO
        for index, d, t_exp, shift, pole_order in entries:
            c = row[index]
            if not c:
                continue
            kernel_terms += 1
            term = t1_powers[t_exp]
            # multiply explicitly by L^shift
            Lshift = A({(shift, 0): finite.one()})
            term = scale_pair(c * Lshift, term)
            family_sum = add_pair(family_sum, term)
        if family_sum != ZERO:
            image = add_pair(image, mul_pair(family_factor[family], family_sum))
    images.append(image)
    terms_out = len(image[0].dict()) + len(image[1].dict())
    print(
        "TRANSLATED_DIVISOR_SURVIVOR|"
        f"index={survivor_index+1}/18|kernel_terms={kernel_terms}|"
        f"obstruction_terms={terms_out}",
        flush=True,
    )

build_seconds = time.perf_counter() - build_start

# Coordinate the rank-two principal parts: L^ell U^u * {1,Y}, ell<26.
coords = sorted({
    (int(le), int(ue), parity)
    for even, odd in images
    for parity, poly in ((0, even), (1, odd))
    for (le, ue) in poly.dict()
})
coord_index = {key: i for i, key in enumerate(coords)}
M = matrix(finite, len(coords), 18)

for col, (even, odd) in enumerate(images):
    for parity, poly in ((0, even), (1, odd)):
        for (le, ue), coeff in poly.dict().items():
            M[coord_index[(int(le), int(ue), parity)], col] = coeff

rank_start = time.perf_counter()
rank = int(M.rank())
rank_seconds = time.perf_counter() - rank_start
remaining = 18 - rank
nonzero = sum(bool(even or odd) for even, odd in images)

OUTPUT = GEN / f"zz-h92-q8-true1600-translated-divisor-E7_2-E7_5-mod-{P}.json"
payload = {
    "schema": "elkies-k3.h92-q8-true1600-translated-divisor-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_TRUE1600_TRANSLATED_DIVISOR_IMAGE",
    "prime": int(P),
    "chart": args.chart,
    "divisor": {
        "equation_in_original_chart": str(L_old),
        "new_parameter": "L",
        "node_local_unit": True,
        "orders": orders,
        "surface_rank_two_relation": f"Y^2={H}",
        "common_clearing_L_order": MOD_L,
    },
    "normalization": {
        "actual_degree_used": True,
        "helper_shift": "i=d-9 is already stripped; no extra q6^9 t^9 is imposed",
        "regularity_order": "ord_L(f)=72-d+x_power",
        "active_column_count": int(len(active_indices)),
    },
    "restricted_image": {
        "global_survivor_dimension": 18,
        "coordinate_rows": int(len(coords)),
        "nonzero_survivor_images": int(nonzero),
        "rank": int(rank),
        "remaining_dimension": int(remaining),
    },
    "timing_seconds": {
        "build": float(build_seconds),
        "rank": float(rank_seconds),
    },
    "inputs": {
        "p1": {"path": str(P1_PATH.relative_to(ROOT)), "sha256": digest(P1_PATH)},
        "pullbacks": {"path": str(PULLBACKS_PATH.relative_to(ROOT)), "sha256": digest(PULLBACKS_PATH)},
        "gluing": {"path": str(GLUING_PATH.relative_to(ROOT)), "sha256": digest(GLUING_PATH)},
        "ambient": {"path": str(AMBIENT_PATH.relative_to(ROOT)), "sha256": digest(AMBIENT_PATH)},
        "global_kernel": {"path": str(KERNEL_PATH.relative_to(ROOT)), "sha256": digest(KERNEL_PATH)},
    },
    "boundary": (
        "This tests the translated third t-divisor visible on the actual "
        "E7_2--E7_5 chart. It is deliberately not a whole-chart Groebner "
        "membership calculation. If active, its conditions should be stacked "
        "with the already-verified root-component and node conditions; the "
        "second translated chart E7_3--E7_6 should then be audited separately."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "Q8TRUETRANSLATEDDIVISOR|"
    f"prime={P}|chart={args.chart}|global_survivor=18|"
    f"rows={len(coords)}|nonzero_images={nonzero}|"
    f"restricted_rank={rank}|remaining={remaining}|"
    f"build_seconds={build_seconds:.4f}|rank_seconds={rank_seconds:.4f}|"
    "status=EXPERIMENTAL_MODULAR_TRUE1600_TRANSLATED_DIVISOR_IMAGE",
    flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
