#!/usr/bin/env sage -python
"""
Exact higher-normal-jet test along the generic E7_5 component, restricted to
the true 18-dimensional H92 q=8 survivor, and stacked with the already-found
translated/affine-component rank-8 obstruction.

Why this is needed
------------------
The existing generic E7_5 residue compiler uses only the INITIAL Z-coefficient
of t,x and the marked chord for each starting negative normal order.  It
therefore imposes the leading residue at each order, but after such residues
cancel a linear combination can still have further negative Z-jets.

On the actual E7_2--E7_5 chart the generic E7_5 component is Z=0 and
    ord_Z(t,x,y) = (1,2,3),
    ord_Z(x-xP) = ord_Z(y-yP) = 3,
    ord_Z(g) = 6.
After the repository's required q6^9 helper shift i=d-9 has been stripped,
the corrected reduced residual order of an actual true-ambient term is
    6 + (72-d) + 2*a = 78-d+2*a.
Thus the worst remaining pole is 9.

We compute the COMPLETE Z-principal part modulo Z^9, not merely its first
residue, using the exact monic surface relation Y^2=H(Z,U).  Everything is
rank-two arithmetic over GF(p)[Z,U]; no Groebner basis is used.

A common unit clearing gives polynomial columns proportional to
    Z^(87-d+2a) * (t/Z)^(89-d)
    * (g/Z^6) * (x/Z^2)^a
    * ((y-yP)/Z^3)^b * ((x-xP)/Z^3)^(9-b)
    * xP_den^b * yP_den^(9-b).
Modulo Z^9 only 162 of the 1600 ambient columns can contribute.

Run:
  sage -python ~/Downloads/probe_h92_q8_true1600_e75_full_normal_jets.sage

Optional:
  --repo /path/to/jacobian-research
  --prime 43
"""

import argparse
import hashlib
import json
import time
from collections import defaultdict
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, sage_eval


MOD_Z = 9
CHART = "E7_2--E7_5"


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
TRANSLATED_PATH = GEN / f"zz-h92-q8-true1600-two-translated-divisors-mod-{P}.json"
OUTPUT_PATH = GEN / f"zz-h92-q8-true1600-e75-full-normal-jets-mod-{P}.json"

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
chart = charts[CHART]
edge = edges[CHART]
assert edge["w_cartier_equation"] == "Z^6*Y^5"

# Exact chart over GF(p).
R = PolynomialRing(finite, names=("Z", "U", "Y"))
Z, U, Y = R.gens()

QQR = PolynomialRing(QQ, names=("Z", "U", "Y"))
Zq, Uq, Yq = QQR.gens()
qloc = {"Z": Zq, "U": Uq, "Y": Yq}

surface = R(QQR(sage_eval(chart["surface_equation"], locals=qloc)))
t_map = R(QQR(sage_eval(chart["old_coordinate_pullback"]["t"], locals=qloc)))
x_map = R(QQR(sage_eval(chart["old_coordinate_pullback"]["x"], locals=qloc)))
y_map = R(QQR(sage_eval(chart["old_coordinate_pullback"]["y"], locals=qloc)))
g_map = R(QQR(sage_eval(edge["w_cartier_equation"], locals=qloc)))

# P1 as exact reversed t-functions.
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

xp_n = R(xp_num(t_map))
xp_d = R(xp_den(t_map))
yp_n = R(yp_num(t_map))
yp_d = R(yp_den(t_map))
assert xp_d.subs({Z: 0}) != 0
assert yp_d.subs({Z: 0}) != 0

nx = R(x_map * xp_d - xp_n)
ny = R(y_map * yp_d - yp_n)


def z_valuation(poly):
    if not poly:
        return 10**9
    return min(int(exp[0]) for exp in poly.dict())


orders = {
    "t": z_valuation(t_map),
    "x": z_valuation(x_map),
    "y": z_valuation(y_map),
    "Nx": z_valuation(nx),
    "Ny": z_valuation(ny),
    "g": z_valuation(g_map),
}
expected_orders = {"t": 1, "x": 2, "y": 3, "Nx": 3, "Ny": 3, "g": 6}
assert orders == expected_orders, orders

# Strip exact Z-powers. The residues need only be nonzero generically on E7_5;
# they need not be units at the marked point.
t1 = R(t_map // Z)
x2 = R(x_map // Z**2)
nx3 = R(nx // Z**3)
ny3 = R(ny // Z**3)
g6 = R(g_map // Z**6)
assert t1.subs({Z: 0}) != 0
assert x2.subs({Z: 0}) != 0
assert nx3.subs({Z: 0}) != 0
assert ny3.subs({Z: 0}) != 0
assert g6.subs({Z: 0}) != 0

# Surface is monic quadratic in Y: Y^2=H(Z,U).
H_ZU = R(Y**2 - surface)
assert H_ZU.degree(Y) == 0
assert surface == Y**2 - H_ZU

# Rank-two module A + A*Y over A=GF(p)[Z,U], truncated modulo Z^9.
A = PolynomialRing(finite, names=("Z", "U"))
Za, Ua = A.gens()
H = A(H_ZU(Za, Ua, 0))


def trunc(poly):
    poly = A(poly)
    return A({
        exp: coeff
        for exp, coeff in poly.dict().items()
        if int(exp[0]) < MOD_Z
    })


ZERO = (A.zero(), A.zero())
ONE = (A.one(), A.zero())


def add_pair(left, right):
    return (trunc(left[0] + right[0]), trunc(left[1] + right[1]))


def scale_pair(c, value):
    return (trunc(c * value[0]), trunc(c * value[1]))


def mul_pair(left, right):
    return (
        trunc(left[0] * right[0] + left[1] * right[1] * H),
        trunc(left[0] * right[1] + left[1] * right[0]),
    )


def pair_power(value, exponent):
    answer = ONE
    base = value
    while exponent:
        if exponent & 1:
            answer = mul_pair(answer, base)
        exponent //= 2
        if exponent:
            base = mul_pair(base, base)
    return answer


hp_cache = {0: A.one()}


def to_pair(poly):
    poly = R(poly)
    even = A.zero()
    odd = A.zero()
    for (ze, ue, ye), coeff in poly.dict().items():
        q, parity = divmod(int(ye), 2)
        if q not in hp_cache:
            hp_cache[q] = trunc(H**q)
        mon = A({(int(ze), int(ue)): coeff})
        term = trunc(mon * hp_cache[q])
        if parity:
            odd = trunc(odd + term)
        else:
            even = trunc(even + term)
    return (even, odd)


t1p = to_pair(t1)
x2p = to_pair(x2)
nx3p = to_pair(nx3)
ny3p = to_pair(ny3)
g6p = to_pair(g6)
xpdp = to_pair(xp_d)
ypdp = to_pair(yp_d)

assert t1p[1] == 0
assert x2p[1] == 0
assert nx3p[1] == 0
assert xpdp[1] == 0 and ypdp[1] == 0

# Family factor after the common generic-component unit clearing:
# nx3^9*yp_den^9*h^18*t1^17*Z^9.
family_factor = {}
for a in (0, 1):
    max_b = 9 if a == 0 else 7
    for b in range(max_b + 1):
        value = g6p
        for base_value, exponent in (
            (x2p, a),
            (ny3p, b),
            (nx3p, 9 - b),
            (xpdp, b),
            (ypdp, 9 - b),
        ):
            if exponent:
                value = mul_pair(value, pair_power(base_value, exponent))
        family_factor[(a, b)] = value

# t1^(89-d), 0<=89-d<=17.
t1_powers = [ONE]
for exponent in range(1, 18):
    t1_powers.append(mul_pair(t1_powers[-1], t1p))


def actual_degree(entry):
    if "actual_u_power" in entry:
        return int(entry["actual_u_power"])
    return int(entry["u_power"]) + 9


# Correct reduced residual order is 78-d+2a; multiply by Z^9, so the explicit
# cleared shift is 87-d+2a. Only shifts 0..8 survive modulo Z^9.
active_by_family = defaultdict(list)
for index, entry in enumerate(basis):
    a = int(entry["x_power"])
    b = int(entry["m_power"])
    d = actual_degree(entry)
    k = int(entry["h_power"])
    assert k == 18
    residual_order = 78 - d + 2 * a
    if residual_order >= 0:
        continue
    shift = MOD_Z + residual_order  # =87-d+2a
    t_exp = 89 - d
    assert shift == 87 - d + 2 * a
    assert 0 <= shift < MOD_Z
    assert 0 <= t_exp <= 17
    active_by_family[(a, b)].append(
        (index, d, residual_order, shift, t_exp)
    )

active_indices = sorted(
    index
    for entries in active_by_family.values()
    for index, *_ in entries
)
assert len(active_indices) == 162, len(active_indices)

print(
    "E75FULLJETS_START|"
    f"prime={P}|component=E7_5|chart={CHART}|normal=Z|"
    f"orders=t:{orders['t']},x:{orders['x']},Nx:{orders['Nx']},"
    f"Ny:{orders['Ny']},g:{orders['g']}|"
    f"mod=Z^{MOD_Z}|active_columns={len(active_indices)}|survivors=18",
    flush=True,
)

images = []
build_start = time.perf_counter()

for survivor_index, row in enumerate(kernel_rows):
    image = ZERO
    kernel_terms = 0
    for family, entries in active_by_family.items():
        family_sum = ZERO
        for index, d, residual_order, shift, t_exp in entries:
            c = row[index]
            if not c:
                continue
            kernel_terms += 1
            Zshift = A({(shift, 0): finite.one()})
            term = scale_pair(c * Zshift, t1_powers[t_exp])
            family_sum = add_pair(family_sum, term)
        if family_sum != ZERO:
            image = add_pair(image, mul_pair(family_factor[family], family_sum))
    images.append(image)
    terms_out = len(image[0].dict()) + len(image[1].dict())
    print(
        "E75FULLJETS_SURVIVOR|"
        f"index={survivor_index+1}/18|kernel_terms={kernel_terms}|"
        f"principal_terms={terms_out}",
        flush=True,
    )

build_seconds = time.perf_counter() - build_start

# Coordinate exact Z-principal parts Z^z U^u * {1,Y}, z<9.
coords = sorted({
    (int(ze), int(ue), parity)
    for even, odd in images
    for parity, poly in ((0, even), (1, odd))
    for (ze, ue) in poly.dict()
})
coord_index = {key: i for i, key in enumerate(coords)}
M = matrix(finite, len(coords), 18)

for col, (even, odd) in enumerate(images):
    for parity, poly in ((0, even), (1, odd)):
        for (ze, ue), coeff in poly.dict().items():
            M[coord_index[(int(ze), int(ue), parity)], col] = coeff

rank_start = time.perf_counter()
full_rank = int(M.rank())
rank_seconds = time.perf_counter() - rank_start
nonzero = sum(bool(even or odd) for even, odd in images)

# Stack with the already-computed translated/affine component obstruction.
translated_rank = None
combined_rank = None
incremental_gain = None
if TRANSLATED_PATH.exists():
    translated = json.loads(TRANSLATED_PATH.read_text())
    assert int(translated["prime"]) == P
    row_basis = translated["divisors"][0]["row_space_basis"]
    translated_rows = [[finite(v) for v in row] for row in row_basis]
    translated_rank = int(matrix(finite, translated_rows).rank()) if translated_rows else 0
    combined_rows = translated_rows + [list(row) for row in M.row_space().basis()]
    combined_rank = int(matrix(finite, combined_rows).rank()) if combined_rows else 0
    incremental_gain = combined_rank - translated_rank

    if P == 43:
        assert translated_rank == 8, translated_rank

print(
    "E75FULLJETS_RESULT|"
    f"prime={P}|component=E7_5|rows={len(coords)}|"
    f"nonzero_images={nonzero}|full_rank={full_rank}|"
    f"translated_rank={translated_rank if translated_rank is not None else 'NA'}|"
    f"gain_after_translated={incremental_gain if incremental_gain is not None else 'NA'}|"
    f"combined_rank={combined_rank if combined_rank is not None else 'NA'}|"
    f"remaining={18-combined_rank if combined_rank is not None else 18-full_rank}|"
    f"build_seconds={build_seconds:.4f}|rank_seconds={rank_seconds:.4f}",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-q8-true1600-e75-full-normal-jets-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_TRUE1600_E75_FULL_NORMAL_JETS",
    "prime": int(P),
    "component": "E7_5",
    "chart": CHART,
    "normal_parameter": "Z",
    "normalization": {
        "helper_shift": "i=d-9 already stripped",
        "correct_reduced_residual_order": "78-d+2*x_power",
        "worst_negative_order": -9,
        "principal_part_modulus": "Z^9",
        "active_ambient_columns": int(len(active_indices)),
        "method": (
            "Complete normal principal part in the exact rank-two surface "
            "module GF(p)[Z,U] + Y*GF(p)[Z,U], not only the initial residue."
        ),
    },
    "component_orders": orders,
    "restricted_image": {
        "global_survivor_dimension": 18,
        "coordinate_rows": int(len(coords)),
        "nonzero_survivor_images": int(nonzero),
        "rank": int(full_rank),
        "remaining_if_alone": int(18 - full_rank),
        "row_space_basis": [
            [int(v) for v in row]
            for row in M.row_space().basis()
        ],
    },
    "stack_with_translated": {
        "translated_artifact_present": bool(TRANSLATED_PATH.exists()),
        "translated_rank": translated_rank,
        "incremental_rank_gain": incremental_gain,
        "combined_rank": combined_rank,
        "remaining_dimension": (
            int(18 - combined_rank) if combined_rank is not None else None
        ),
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
        "translated": (
            {"path": str(TRANSLATED_PATH.relative_to(ROOT)), "sha256": digest(TRANSLATED_PATH)}
            if TRANSLATED_PATH.exists() else None
        ),
    },
    "boundary": (
        "This completes the generic normal Z-principal part on E7_5 beyond "
        "the repository's initial-residue rows and stacks it with the "
        "translated/affine-component obstruction. It does not yet prove that "
        "other E7 components have no higher normal-jet contribution, nor "
        "supply a second-prime or characteristic-zero certificate."
    ),
}
OUTPUT_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "Q8TRUEE75FULLJETS|"
    f"prime={P}|global_survivor=18|component_rank={full_rank}|"
    f"translated_rank={translated_rank if translated_rank is not None else 'NA'}|"
    f"gain={incremental_gain if incremental_gain is not None else 'NA'}|"
    f"combined_rank={combined_rank if combined_rank is not None else 'NA'}|"
    f"remaining={18-combined_rank if combined_rank is not None else 18-full_rank}|"
    "status=EXPERIMENTAL_MODULAR_TRUE1600_E75_FULL_NORMAL_JETS",
    flush=True,
)
print(f"OUTPUT|{OUTPUT_PATH}", flush=True)
