#!/usr/bin/env sage -python
"""
Fast exact finite-node probe for the H92 q=8 true-1600 global survivor.

Instead of mapping all 1600 ambient columns into R/(t^T), this script:
  1. loads the already-computed 18-dimensional global kernel;
  2. discards ambient columns carrying t^e with e >= T (zero mod t^T);
  3. groups the remaining columns by the 18 RR families;
  4. combines them into the 18 global survivor numerators first;
  5. performs only 18 Singular local normal-form reductions.

For E7_4--E7_3 it additionally replaces the local generator
    Z^(3T) U^(2T)
by the exactly equivalent local generator
    Z^(3T) Y^(4T),
using surface = Y^2-U*H with H(0) a unit.  This usually makes the local
standard-basis computation much simpler.  Use --no-unit-swap to retain the
original generator exactly as written in the older probe.

Run:
  sage -python ~/Downloads/probe_h92_q8_true1600_restricted_node_fast.sage

Optional:
  --repo /path/to/jacobian-research
  --prime 43
  --node E7_4--E7_3
  --no-unit-swap
  --local-standard-basis lazard
  --verify-idempotent
"""

import argparse
import hashlib
import json
import time
from collections import defaultdict
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, matrix, sage_eval, singular


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
    raise SystemExit(
        "Could not locate jacobian-research. Re-run with "
        "--repo /path/to/jacobian-research"
    )


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def invert_base(rational_u):
    numerator = rational_u.numerator()
    denominator = rational_u.denominator()
    t_ring = PolynomialRing(QQ, "t")
    t = t_ring.gen()
    field = t_ring.fraction_field()
    return field(
        t ** (denominator.degree() - numerator.degree())
        * t_ring(list(reversed(numerator.list())))
        / t_ring(list(reversed(denominator.list())))
    )


def common_monomial_exponents(value):
    terms = list(value.dict())
    assert terms
    return tuple(min(exponent[index] for exponent in terms) for index in range(3))


def reduce_coefficient(value, finite):
    value = QQ(value)
    denominator = finite(value.denominator())
    if not denominator:
        raise ValueError("prime divides an input coefficient denominator")
    return finite(value.numerator()) / denominator


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path, default=None)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--node", default="E7_4--E7_3")
parser.add_argument(
    "--local-standard-basis", choices=("std", "lazard"), default="std"
)
parser.add_argument(
    "--no-unit-swap", action="store_true",
    help="Do not replace Z^(3T)U^(2T) by the equivalent Z^(3T)Y^(4T) at E7_4--E7_3."
)
parser.add_argument(
    "--verify-idempotent", action="store_true",
    help="Reduce each of the 18 remainders a second time as a regression check."
)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
P = int(args.prime)
if P <= 1:
    raise ValueError("prime must be greater than one")

GEN = ROOT / "artifacts" / "generated-results"
P1_PATH = GEN / "elkies-k3-h92-p1-lift.json"
PULLBACKS_PATH = GEN / "elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
GLUING_PATH = GEN / "elkies-k3-h92-q8-actual-e7-gluing.json"
AMBIENT_PATH = GEN / "zz-h92-q8-true1600-ambient.json"
CLEARINGS_PATH = GEN / "zz-h92-q8-true1600-node-principal-clearings.json"
KERNEL_PATH = GEN / f"zz-h92-q8-true1600-global-kernel-mod-{P}.json"
OUTPUT_PATH = GEN / (
    f"zz-h92-q8-true1600-{args.node.replace('--', '-')}-restricted-fast-mod-{P}.json"
)

for path in (P1_PATH, PULLBACKS_PATH, GLUING_PATH, AMBIENT_PATH, CLEARINGS_PATH, KERNEL_PATH):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

p1 = json.loads(P1_PATH.read_text())
pullbacks = json.loads(PULLBACKS_PATH.read_text())
gluing = json.loads(GLUING_PATH.read_text())
ambient = json.loads(AMBIENT_PATH.read_text())
clearings = json.loads(CLEARINGS_PATH.read_text())
kernel = json.loads(KERNEL_PATH.read_text())

assert p1["status"] == "PASS_EXACT_H92_P1"
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"
assert gluing["status"] == "PASS_EXACT_Q8_ACTUAL_E7_GLUING"
assert clearings["status"] == "PASS_EXACT_Q8_E7_NODE_PRINCIPAL_CLEARINGS"
assert ambient["status"] in {
    "PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT",
    "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT",
}
assert int(ambient["ambient_dimension"]) == 1600
assert int(kernel["prime"]) == P
assert int(kernel["dimensions"]["ambient"]) == 1600
assert int(kernel["dimensions"]["kernel"]) == 18

node = next((item for item in clearings["nodes"] if item["chart"] == args.node), None)
if node is None:
    raise ValueError(f"unknown cleared E7 chart: {args.node}")

T = int(clearings["common_parameters"]["T"])
K = int(clearings["common_parameters"]["K"])
assert (T, K) == (17, 18), (T, K)
assert K == max(int(entry["h_power"]) for entry in ambient["ambient_basis"])
assert T == 9 + max(
    int(entry["u_power"]) - 4 * int(entry["h_power"])
    for entry in ambient["ambient_basis"]
)

finite = GF(P)
kernel_rows = [
    [finite(value) for value in row]
    for row in kernel["kernel_basis_rows"]
]
assert len(kernel_rows) == 18
assert all(len(row) == 1600 for row in kernel_rows)

print(
    f"FASTNODE|prime={P}|node={args.node}|ambient=1600|survivors=18|T={T}|K={K}",
    flush=True,
)

# ---------------------------------------------------------------------------
# Reconstruct the exact chart data, matching the audited local probe.
# ---------------------------------------------------------------------------

u_ring = PolynomialRing(QQ, "u")
u_field = u_ring.fraction_field()
x_p = u_field(polynomial(u_ring, p1["x_entrance_base"]["numerator_coefficients"]))
x_p /= u_field(polynomial(u_ring, p1["x_entrance_base"]["denominator_coefficients"]))
y_p = u_field(polynomial(u_ring, p1["y_entrance_base"]["numerator_coefficients"]))
y_p /= u_field(polynomial(u_ring, p1["y_entrance_base"]["denominator_coefficients"]))
x_p_t = invert_base(x_p)
y_p_t = invert_base(y_p)
t_ring = x_p_t.parent()
t_formal = t_ring.gen()
r, s = x_p_t / t_formal**2, y_p_t / t_formal**3
assert r.valuation() == 0 and s.valuation() == 0

t_poly = PolynomialRing(QQ, "t")
r_num, r_den = t_poly(r.numerator()), t_poly(r.denominator())
s_num, s_den = t_poly(s.numerator()), t_poly(s.denominator())
h_reverse = t_poly(list(reversed(polynomial(
    u_ring, p1["structured_denominator"]["Z4_coefficients"]
).list())))

chart = next(item for item in pullbacks["charts"] if item["name"] == args.node)
edge = next(item for item in gluing["actual_edge_chart_gluing"] if item["name"] == args.node)

qq_ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Zq, Uq, Yq = qq_ring.gens()
locals_qq = {"Z": Zq, "U": Uq, "Y": Yq}
surface_qq = qq_ring(sage_eval(chart["surface_equation"], locals=locals_qq))
t_qq = qq_ring(sage_eval(chart["old_coordinate_pullback"]["t"], locals=locals_qq))
x_qq = qq_ring(sage_eval(chart["old_coordinate_pullback"]["x"], locals=locals_qq))
y_qq = qq_ring(sage_eval(chart["old_coordinate_pullback"]["y"], locals=locals_qq))
g_qq = qq_ring(sage_eval(edge["w_cartier_equation"], locals=locals_qq))

ring = PolynomialRing(finite, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
surface, t_value, x_value, y_value, g = tuple(
    ring(value) for value in (surface_qq, t_qq, x_qq, y_qq, g_qq)
)


def evaluate_t_polynomial(value):
    return ring(sum(
        reduce_coefficient(coefficient, finite) * t_value**degree
        for degree, coefficient in enumerate(value.list())
    ))


R_num, R_den, S_num, S_den, H_reverse = tuple(
    evaluate_t_polynomial(value)
    for value in (r_num, r_den, s_num, s_den, h_reverse)
)
numerator_x = x_value * R_den - t_value**2 * R_num
numerator_y = y_value * S_den - t_value**3 * S_num
mx, my = common_monomial_exponents(numerator_x), common_monomial_exponents(numerator_y)
m_exponents = tuple(my[index] - mx[index] for index in range(3))
assert all(exponent >= 0 for exponent in m_exponents)
mx_value, my_value = ring.monomial(*mx), ring.monomial(*my)
A, B = numerator_x // mx_value, numerator_y // my_value
assert A(0, 0, 0) and R_den(0, 0, 0) and S_den(0, 0, 0) and H_reverse(0, 0, 0)

t_exponents = common_monomial_exponents(t_value)
t_monomial = ring.monomial(*t_exponents)
assert (t_value // t_monomial)(0, 0, 0)


def truncate_t_power(value):
    """Discard monomials already divisible by t_monomial^T."""
    return ring({
        monomial: coefficient
        for monomial, coefficient in ring(value).dict().items()
        if not all(
            monomial[index] >= T * t_exponents[index]
            for index in range(3)
        )
    })


def truncated_product(left, right):
    return truncate_t_power(ring(left) * ring(right))


def truncated_power(value, exponent):
    answer = ring.one()
    value = truncate_t_power(value)
    while exponent:
        if exponent & 1:
            answer = truncated_product(answer, value)
        exponent //= 2
        if exponent:
            value = truncated_product(value, value)
    return answer


m_monomial = ring.monomial(*m_exponents)

# ---------------------------------------------------------------------------
# Key optimization 1: t^e is zero modulo t^T whenever e >= T.
# For this true ambient e = 89-d, leaving only 286/1600 columns.
# ---------------------------------------------------------------------------

active_by_family = defaultdict(list)
skipped = 0
for index, entry in enumerate(ambient["ambient_basis"]):
    a = int(entry["x_power"])
    b = int(entry["m_power"])
    i = int(entry["u_power"])
    k = int(entry["h_power"])
    assert k == K
    t_exponent = T + 4 * k - i - 9
    assert t_exponent >= 0
    if t_exponent >= T:
        skipped += 1
        continue
    active_by_family[(a, b)].append((index, t_exponent))

active = 1600 - skipped
assert skipped == 1314, skipped
assert active == 286, active
assert len(active_by_family) == 18
print(
    f"TRUNCATION|zero_columns={skipped}|active_columns={active}|families={len(active_by_family)}",
    flush=True,
)

# Precompute t^0,...,t^(T-1) only, always truncating along the way.
t_powers = [ring.one()]
for exponent in range(1, T):
    t_powers.append(truncated_product(t_powers[-1], t_value))

# Precompute the d-independent factor once per RR family.
family_factor = {}
for (a, b), terms in sorted(active_by_family.items()):
    answer = ring.one()
    for factor in (
        g,
        x_value**a,
        m_monomial**b,
        truncated_power(B, b),
        truncated_power(R_den, b),
        truncated_power(A, 9 - b),
        truncated_power(S_den, 9 - b),
        # H_reverse^(K-k)=1 in the true uniform h^-18 ambient.
    ):
        answer = truncated_product(answer, factor)
    family_factor[(a, b)] = answer

# ---------------------------------------------------------------------------
# Key optimization 2: combine ambient columns into the 18 global survivors
# before any Singular normal-form reduction.
# ---------------------------------------------------------------------------

build_start = time.perf_counter()
survivor_numerators = []
active_coefficients = 0
for survivor_index, row in enumerate(kernel_rows):
    numerator = ring.zero()
    for family, terms in active_by_family.items():
        t_polynomial = ring.zero()
        for ambient_index, t_exponent in terms:
            coefficient = row[ambient_index]
            if coefficient:
                active_coefficients += 1
                t_polynomial += coefficient * t_powers[t_exponent]
        if t_polynomial:
            numerator += truncated_product(family_factor[family], t_polynomial)
    survivor_numerators.append(truncate_t_power(numerator))
    print(
        f"BUILD_SURVIVOR|index={survivor_index+1}/18|terms={len(survivor_numerators[-1].dict())}",
        flush=True,
    )
build_seconds = time.perf_counter() - build_start
print(
    f"SURVIVOR_BUILD_DONE|seconds={build_seconds:.3f}|active_kernel_coefficients={active_coefficients}",
    flush=True,
)

# ---------------------------------------------------------------------------
# Exact local ideal and standard basis.
# ---------------------------------------------------------------------------

singular.eval(f"ring q8node={P},(Z,U,Y),ds;")
singular.eval(f"poly surface={surface};")

original_generator = t_monomial**T
principal_generator = original_generator
generator_note = "original t_monomial^T"

# At E7_4--E7_3 the audited chart has
#   t = Z^3 U^2, surface = Y^2-U*H, H(0)=1.
# Hence in the local ring U^(2T) and Y^(4T) differ by the unit H^(2T), so
# (surface, Z^(3T)U^(2T)) = (surface, Z^(3T)Y^(4T)).
# The quotient is therefore exactly the same, while the second leading
# monomial no longer contains U.
if args.node == "E7_4--E7_3" and not args.no_unit_swap:
    assert t_monomial == Z**3 * U**2
    unit_h, remainder = (Y**2 - surface).quo_rem(U)
    assert not remainder
    assert unit_h(0, 0, 0) != 0
    principal_generator = Z**(3 * T) * Y**(4 * T)
    generator_note = "unit-equivalent Z^(3T)Y^(4T) using Y^2=U*H, H(0)!=0"
    print(
        f"IDEAL_REWRITE|original=Z^{3*T}*U^{2*T}|equivalent=Z^{3*T}*Y^{4*T}|status=EXACT_LOCAL_UNIT_EQUIVALENCE",
        flush=True,
    )

singular.eval(f"ideal principal=surface,{principal_generator};")
print(
    f"STDBASIS_START|method={args.local_standard_basis}|generator={generator_note}",
    flush=True,
)
std_start = time.perf_counter()
if args.local_standard_basis == "std":
    singular.eval("ideal standard_principal=std(principal);")
else:
    singular.eval('LIB "teachstd.lib";')
    singular.eval("ideal standard_principal=localstd(principal);")
std_seconds = time.perf_counter() - std_start
print(f"STDBASIS_DONE|seconds={std_seconds:.3f}", flush=True)

# Exactly 18 reductions, not 1600.  The second idempotence reduction from the
# old audit script is optional because it doubles this cost and is not needed
# to determine the restricted rank.
remainders = []
reduce_start = time.perf_counter()
for index, numerator in enumerate(survivor_numerators):
    one_start = time.perf_counter()
    singular.eval(f"poly node_num={numerator};")
    remainder = ring(singular("reduce(node_num,standard_principal)").sage())
    if args.verify_idempotent:
        singular.eval(f"poly remainder={remainder};")
        assert ring(singular("reduce(remainder,standard_principal)").sage()) == remainder
    remainders.append(remainder)
    print(
        f"REDUCE_SURVIVOR|index={index+1}/18|terms={len(remainder.dict())}|seconds={time.perf_counter()-one_start:.3f}",
        flush=True,
    )
reduce_seconds = time.perf_counter() - reduce_start

# Build the 18-column image matrix directly from remainder monomials.
monomials = sorted({
    tuple(monomial)
    for remainder in remainders
    for monomial in remainder.dict()
})
monomial_index = {monomial: index for index, monomial in enumerate(monomials)}
restricted = matrix(finite, len(monomials), 18)
for column, remainder in enumerate(remainders):
    for monomial, coefficient in remainder.dict().items():
        restricted[monomial_index[tuple(monomial)], column] = coefficient

restricted_rank = int(restricted.rank())
remaining = 18 - restricted_rank

payload = {
    "schema": "elkies-k3.h92-q8-true1600-restricted-node-fast-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_TRUE1600_EXACT_RESTRICTED_NODE_IMAGE",
    "prime": int(P),
    "node": args.node,
    "inputs": {
        "p1": {"path": str(P1_PATH.relative_to(ROOT)), "sha256": digest(P1_PATH)},
        "pullbacks": {"path": str(PULLBACKS_PATH.relative_to(ROOT)), "sha256": digest(PULLBACKS_PATH)},
        "gluing": {"path": str(GLUING_PATH.relative_to(ROOT)), "sha256": digest(GLUING_PATH)},
        "ambient": {"path": str(AMBIENT_PATH.relative_to(ROOT)), "sha256": digest(AMBIENT_PATH)},
        "clearings": {"path": str(CLEARINGS_PATH.relative_to(ROOT)), "sha256": digest(CLEARINGS_PATH)},
        "global_kernel": {"path": str(KERNEL_PATH.relative_to(ROOT)), "sha256": digest(KERNEL_PATH)},
    },
    "optimization": {
        "ambient_dimension": 1600,
        "columns_zero_mod_tT": int(skipped),
        "active_ambient_columns": int(active),
        "global_survivor_dimension": 18,
        "normal_form_reductions": 18,
        "unit_equivalent_generator_used": bool(
            args.node == "E7_4--E7_3" and not args.no_unit_swap
        ),
        "generator_note": generator_note,
    },
    "local_quotient": {
        "T": int(T),
        "K": int(K),
        "coordinate_rows": int(len(monomials)),
        "restricted_rank": int(restricted_rank),
        "remaining_dimension": int(remaining),
        "coordinate_monomials": [
            f"Z^{a}*U^{b}*Y^{c}" for a, b, c in monomials
        ],
        "sparse_survivor_columns": [
            [
                [monomial_index[tuple(monomial)], int(coefficient)]
                for monomial, coefficient in sorted(remainder.dict().items())
            ]
            for remainder in remainders
        ],
    },
    "timing_seconds": {
        "survivor_build": build_seconds,
        "standard_basis": std_seconds,
        "reductions": reduce_seconds,
    },
    "boundary": (
        "This is the exact modular image of the already-computed 18-dimensional "
        "generic global survivor in one actual E7 node quotient.  It does not "
        "yet impose the other five nodes, marked/overlap gluing, or certify a "
        "characteristic-zero pencil."
    ),
}
OUTPUT_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "Q8TRUEFINITENODEFAST|"
    f"prime={P}|node={args.node}|global_survivor=18|"
    f"active_ambient={active}|node_rows={len(monomials)}|"
    f"restricted_rank={restricted_rank}|remaining={remaining}|"
    f"std_seconds={std_seconds:.3f}|reduce_seconds={reduce_seconds:.3f}|"
    "status=EXPERIMENTAL_MODULAR_TRUE1600_EXACT_RESTRICTED_NODE_IMAGE",
    flush=True,
)
print(f"OUTPUT|{OUTPUT_PATH}", flush=True)
