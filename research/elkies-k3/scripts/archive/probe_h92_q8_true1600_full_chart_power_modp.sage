#!/usr/bin/env sage -python
"""
Exact full resolved-chart q8 E7 ninth-power membership on the true 1600
ambient, restricted immediately to the 18-dimensional global survivor.

Default chart: E7_2--E7_5 (marked/cancellation-sensitive chart).

This tests the exact cleared condition
    (x-xP)^9 * g * f  in  t^9 * (x-xP,y-yP)^9
on the selected actual resolved chart, after localizing only at the known
P1/h denominators that are units along the E7 fibre.

Unlike the node probes, this is a WHOLE-CHART membership test. It can therefore
see finite conditions away from the edge origin, including the marked branch
and chart-cover effects.

Run:
  sage -python ~/Downloads/probe_h92_q8_true1600_full_chart_power_modp.sage

Options:
  --repo /path/to/jacobian-research
  --prime 43
  --chart E7_2--E7_5
"""

import argparse
import hashlib
import json
import time
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
assert len(basis) == 1600
kernel_rows = [[finite(v) for v in row] for row in kernel["kernel_basis_rows"]]
assert len(kernel_rows) == 18
assert all(len(row) == 1600 for row in kernel_rows)

charts = {entry["name"]: entry for entry in pullbacks["charts"]}
edges = {entry["name"]: entry for entry in gluing["actual_edge_chart_gluing"]}
if args.chart not in charts:
    raise ValueError(f"unknown chart: {args.chart}")
chart = charts[args.chart]
edge = edges[args.chart]

# Work over GF(p)[Z,U,Y,S], where S localizes the known denominator product.
ring = PolynomialRing(finite, names=("Z", "U", "Y", "S"), order="degrevlex")
Z, U, Y, S = ring.gens()

# Parse exact QQ chart data first, then reduce coefficients mod p.
qqring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Zq, Uq, Yq = qqring.gens()
locals_qq = {"Z": Zq, "U": Uq, "Y": Yq}
surface_qq = qqring(sage_eval(chart["surface_equation"], locals=locals_qq))
t_qq = qqring(sage_eval(chart["old_coordinate_pullback"]["t"], locals=locals_qq))
x_qq = qqring(sage_eval(chart["old_coordinate_pullback"]["x"], locals=locals_qq))
y_qq = qqring(sage_eval(chart["old_coordinate_pullback"]["y"], locals=locals_qq))
g_qq = qqring(sage_eval(edge["w_cartier_equation"], locals=locals_qq))
surface, t_map, x_map, y_map, g = tuple(
    ring(v) for v in (surface_qq, t_qq, x_qq, y_qq, g_qq)
)

# Exact P1 functions in the E7 base parameter t=1/u, represented as reversed
# numerator/denominator polynomials. Same convention as the repository's
# existing power-membership probe.
t_ring = PolynomialRing(finite, "t")
tt = t_ring.gen()


def reversed_fraction(numerator_values, denominator_values):
    numerator = ff_polynomial(t_ring, finite, numerator_values)
    denominator = ff_polynomial(t_ring, finite, denominator_values)
    rev_numerator = sum(
        numerator[index] * tt ** (numerator.degree() - index)
        for index in range(numerator.degree() + 1)
    )
    rev_denominator = sum(
        denominator[index] * tt ** (denominator.degree() - index)
        for index in range(denominator.degree() + 1)
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
h_u = ff_polynomial(
    t_ring, finite, p1["structured_denominator"]["Z4_coefficients"]
)
h_reverse = sum(
    h_u[index] * tt ** (h_u.degree() - index)
    for index in range(h_u.degree() + 1)
)
assert xp_den(0) and yp_den(0) and h_reverse(0)

xp_n = ring(xp_num(t_map))
xp_d = ring(xp_den(t_map))
yp_n = ring(yp_num(t_map))
yp_d = ring(yp_den(t_map))
h_r = ring(h_reverse(t_map))

# Along the E7 fibre these denominators are units.
assert xp_d(0, 0, 0, 0)
assert yp_d(0, 0, 0, 0)
assert h_r(0, 0, 0, 0)

# Nx/xp_d = x-xP1; Ny/yp_d = y-yP1.
nx = x_map * xp_d - xp_n
ny = y_map * yp_d - yp_n
assert nx and ny

# Only ambient coordinates that actually occur in the 18 survivors matter.
used_indices = [
    index for index in range(1600)
    if any(row[index] for row in kernel_rows)
]
assert used_indices

# IMPORTANT: the true1600 ambient stores helper u_power=d-9 for the generic
# E7 residue compiler, and actual_u_power=d for direct local calculations.
def actual_degree(entry):
    if "actual_u_power" in entry:
        return int(entry["actual_u_power"])
    return int(entry["u_power"]) + 9


common_h_power = max(int(basis[index]["h_power"]) for index in used_indices)
assert common_h_power == 18
common_t_clear = max(
    max(actual_degree(basis[index]) - 4 * int(basis[index]["h_power"]), 0)
    for index in used_indices
)
assert common_t_clear <= 17

# Localize only at factors known to be units near the complete E7 fibre.
unit_denominator = h_r**common_h_power * xp_d**9 * yp_d**9
assert unit_denominator(0, 0, 0, 0)

print(
    "FULLCHART_START|"
    f"prime={P}|chart={args.chart}|ambient=1600|survivors=18|"
    f"used_columns={len(used_indices)}|common_h={common_h_power}|"
    f"common_t_clear={common_t_clear}|power=9",
    flush=True,
)

# Exact localized ninth-power ideal:
#
#   (surface, S*D-1, t^(clear+9) * Nx^(9-j) Ny^j, j=0..9).
#
# The omitted xp_d/yp_d/h factors are units after localization by D.
ideal_generators = [surface, S * unit_denominator - ring.one()]
ideal_generators.extend(
    t_map ** (common_t_clear + 9) * nx ** (9 - j) * ny**j
    for j in range(10)
)
local_ideal = ring.ideal(ideal_generators)

print(
    "FULLCHART_GB_START|"
    f"generators={len(ideal_generators)}|"
    f"max_generator_degree={max(int(v.total_degree()) for v in ideal_generators)}",
    flush=True,
)
gb_start = time.perf_counter()
groebner = local_ideal.groebner_basis()
gb_seconds = time.perf_counter() - gb_start
print(
    "FULLCHART_GB_DONE|"
    f"seconds={gb_seconds:.3f}|basis_size={len(groebner)}|"
    f"max_degree={max((int(v.total_degree()) for v in groebner), default=-1)}",
    flush=True,
)

# Precompute each used ambient column's common-cleared polynomial once.
column_polynomials = {}
build_start = time.perf_counter()
for seq, index in enumerate(used_indices):
    entry = basis[index]
    a = int(entry["x_power"])
    b = int(entry["m_power"])
    d = actual_degree(entry)
    k = int(entry["h_power"])
    assert a in (0, 1) and 0 <= b <= 9 and k == 18
    exponent = common_t_clear + 4 * k - d
    assert exponent >= 0
    # Multiply Lx^9*g*f by common t-clear and by the localized unit D.
    column_polynomials[index] = (
        g
        * x_map**a
        * nx ** (9 - b)
        * ny**b
        * t_map**exponent
        * h_r ** (common_h_power - k)
        * xp_d**b
        * yp_d ** (9 - b)
    )
print(
    "FULLCHART_COLUMNS_DONE|"
    f"columns={len(column_polynomials)}|seconds={time.perf_counter()-build_start:.3f}",
    flush=True,
)

remainders = []
reduce_start = time.perf_counter()
for survivor_index, row in enumerate(kernel_rows):
    build_one = time.perf_counter()
    cleared = ring.zero()
    terms = 0
    for index in used_indices:
        c = row[index]
        if c:
            cleared += c * column_polynomials[index]
            terms += 1
    built_seconds = time.perf_counter() - build_one

    one_reduce = time.perf_counter()
    remainder = cleared.reduce(groebner)
    reduce_seconds = time.perf_counter() - one_reduce
    remainders.append(remainder)
    print(
        "FULLCHART_REDUCE|"
        f"index={survivor_index+1}/18|kernel_terms={terms}|"
        f"input_terms={len(cleared.dict())}|remainder_terms={len(remainder.dict())}|"
        f"build_seconds={built_seconds:.3f}|reduce_seconds={reduce_seconds:.3f}",
        flush=True,
    )
total_reduce_seconds = time.perf_counter() - reduce_start

# The remainder map is linear. Coordinate all 18 remainders in the monomial
# basis actually appearing and compute its exact GF(p) rank.
monomials = sorted({
    tuple(monomial)
    for remainder in remainders
    for monomial in remainder.dict()
})
monomial_index = {m: i for i, m in enumerate(monomials)}
M = matrix(finite, len(monomials), 18)
for column, remainder in enumerate(remainders):
    for monomial, coefficient in remainder.dict().items():
        M[monomial_index[tuple(monomial)], column] = coefficient

rank_start = time.perf_counter()
restricted_rank = int(M.rank())
rank_seconds = time.perf_counter() - rank_start
remaining = 18 - restricted_rank
nonzero_images = sum(bool(r) for r in remainders)

OUTPUT_PATH = GEN / (
    f"zz-h92-q8-true1600-full-chart-{args.chart.replace('--','-')}-"
    f"power-mod-{P}.json"
)
payload = {
    "schema": "elkies-k3.h92-q8-true1600-full-chart-power-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_TRUE1600_EXACT_FULL_CHART_POWER_IMAGE",
    "prime": int(P),
    "chart": args.chart,
    "inputs": {
        "p1": {"path": str(P1_PATH.relative_to(ROOT)), "sha256": digest(P1_PATH)},
        "pullbacks": {"path": str(PULLBACKS_PATH.relative_to(ROOT)), "sha256": digest(PULLBACKS_PATH)},
        "gluing": {"path": str(GLUING_PATH.relative_to(ROOT)), "sha256": digest(GLUING_PATH)},
        "ambient": {"path": str(AMBIENT_PATH.relative_to(ROOT)), "sha256": digest(AMBIENT_PATH)},
        "global_kernel": {"path": str(KERNEL_PATH.relative_to(ROOT)), "sha256": digest(KERNEL_PATH)},
    },
    "clearing": {
        "identity": "(x-xP1)^9*g*f in t^9*(x-xP1,y-yP1)^9",
        "common_h_power": int(common_h_power),
        "common_t_clear": int(common_t_clear),
        "localized_unit_denominator": "h_reverse^18*xP_den^9*yP_den^9",
        "ideal_generator_count": int(len(ideal_generators)),
    },
    "groebner": {
        "basis_size": int(len(groebner)),
        "seconds": float(gb_seconds),
    },
    "restricted_image": {
        "global_survivor_dimension": 18,
        "used_ambient_columns": int(len(used_indices)),
        "coordinate_rows": int(len(monomials)),
        "nonzero_survivor_images": int(nonzero_images),
        "rank": int(restricted_rank),
        "remaining_dimension": int(remaining),
        "sparse_columns": [
            [
                [int(monomial_index[tuple(m)]), int(c)]
                for m, c in sorted(r.dict().items())
            ]
            for r in remainders
        ],
        "coordinate_monomials": [
            {"Z": int(m[0]), "U": int(m[1]), "Y": int(m[2]), "S": int(m[3])}
            for m in monomials
        ],
    },
    "timing_seconds": {
        "groebner": float(gb_seconds),
        "all_survivor_reductions": float(total_reduce_seconds),
        "rank": float(rank_seconds),
    },
    "boundary": (
        "This is exact modular whole-chart membership in the cleared ninth-power "
        "q6 marked module on one actual resolved E7 chart, restricted to the "
        "18-dimensional E8+h+generic-E7 survivor. Other charts must still be "
        "stacked before making a complete finite-E7 or h0 claim."
    ),
}
OUTPUT_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "Q8TRUEFULLCHART|"
    f"prime={P}|chart={args.chart}|global_survivor=18|"
    f"rows={len(monomials)}|nonzero_images={nonzero_images}|"
    f"restricted_rank={restricted_rank}|remaining={remaining}|"
    f"gb_seconds={gb_seconds:.3f}|reduce_seconds={total_reduce_seconds:.3f}|"
    "status=EXPERIMENTAL_MODULAR_TRUE1600_EXACT_FULL_CHART_POWER_IMAGE",
    flush=True,
)
print(f"OUTPUT|{OUTPUT_PATH}", flush=True)
