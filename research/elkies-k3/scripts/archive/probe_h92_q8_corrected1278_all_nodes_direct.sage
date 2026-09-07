#!/usr/bin/env sage -python
"""
Exact all-six finite E7 node probe on the corrected H92 q8 movable-class survivor.

This avoids Singular entirely.

For each actual resolved edge chart:
  * start from the corrected 14-dimensional E8+H+generic-E7 kernel;
  * discard ambient columns with common-clearing t exponent >= T=16;
  * combine the remaining 238 ambient columns into the 14 survivors;
  * use the actual surface relation
        Y^2 = S * H(Z,U),  H(0,0) != 0, S in {Z,U},
    so the local surface ring is a free rank-two module over
        GF(p)[Z,U]_(Z,U)
    with basis {1,Y};
  * strip the unit factor from t = Z^a U^b * unit;
  * reduce Y^(2q+r) exactly to S^q H^q Y^r;
  * test the coefficient pair modulo
        (Z^(aT) U^(bT)).

This gives the exact modular node quotient image on the 14-dimensional
corrected survivor. The known translated-divisor rowspace is stacked first,
so node gains are measured beyond the existing rank-7 affine condition.

Run:
  sage -python ~/Downloads/probe_h92_q8_corrected1278_all_nodes_direct.sage

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

from sage.all import GF, PolynomialRing, QQ, matrix, sage_eval


NODES = [
    "E7_4--E7_3",  # already checked independently; keep first as regression
    "E7_1--E7_4",
    "E7_3--E7_7",
    "E7_7--E7_2",
    "E7_3--E7_6",
    "E7_2--E7_5",  # cancellation-sensitive node
]


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
args = parser.parse_args()

ROOT = locate_repo(args.repo)
P = int(args.prime)
if P <= 1:
    raise ValueError("prime must be greater than one")

GEN = ROOT / "artifacts" / "generated-results"
P1_PATH = GEN / "elkies-k3-h92-p1-lift.json"
PULLBACKS_PATH = GEN / "elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
GLUING_PATH = GEN / "elkies-k3-h92-q8-actual-e7-gluing.json"
AMBIENT_PATH = GEN / "zz-h92-q8-corrected1278-ambient.json"
KERNEL_PATH = GEN / f"zz-h92-q8-corrected1278-global-kernel-mod-{P}.json"
TRANSLATED_PATH = GEN / f"zz-h92-q8-corrected1278-two-translated-divisors-mod-{P}.json"
OUTPUT_PATH = GEN / f"zz-h92-q8-corrected1278-all-e7-nodes-direct-mod-{P}.json"

for path in (
    P1_PATH, PULLBACKS_PATH, GLUING_PATH, AMBIENT_PATH,
    KERNEL_PATH, TRANSLATED_PATH,
):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

p1 = json.loads(P1_PATH.read_text())
pullbacks = json.loads(PULLBACKS_PATH.read_text())
gluing = json.loads(GLUING_PATH.read_text())
ambient = json.loads(AMBIENT_PATH.read_text())
kernel = json.loads(KERNEL_PATH.read_text())
translated = json.loads(TRANSLATED_PATH.read_text())

assert p1["status"] == "PASS_EXACT_H92_P1"
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"
assert gluing["status"] == "PASS_EXACT_Q8_ACTUAL_E7_GLUING"
assert ambient["status"] == "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT"
assert int(ambient["ambient_dimension"]) == 1278
assert int(kernel["prime"]) == P
assert int(kernel["dimensions"]["ambient"]) == 1278
assert int(kernel["dimensions"]["kernel"]) == 14
assert translated["status"] == "EXPERIMENTAL_MODULAR_CORRECTED1278_TWO_TRANSLATED_DIVISOR_IMAGES"
assert int(translated["prime"]) == P
assert int(translated["common"]["global_survivor_dimension"]) == 14

K = 16
T = 16
assert K == max(int(entry["h_power"]) for entry in ambient["ambient_basis"])
assert T == 8 + max(
    int(entry["u_power"]) - 4 * int(entry["h_power"])
    for entry in ambient["ambient_basis"]
)

finite = GF(P)
kernel_rows = [
    [finite(value) for value in row]
    for row in kernel["kernel_basis_rows"]
]
assert len(kernel_rows) == 14
assert all(len(row) == 1278 for row in kernel_rows)

charts = {entry["name"]: entry for entry in pullbacks["charts"]}
edges = {entry["name"]: entry for entry in gluing["actual_edge_chart_gluing"]}
assert set(NODES) == set(charts) == set(edges)

OLD_TWIST = {
    f"E7_{i+1}": value
    for i, value in enumerate((2,5,6,4,6,3,5))
}
NEW_TWIST = {
    f"E7_{i+1}": value
    for i, value in enumerate((2,6,8,5,6,4,7))
}

# P1 entrance data, shared by all six node charts.
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

qq_ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Zq, Uq, Yq = qq_ring.gens()
locals_qq = {"Z": Zq, "U": Uq, "Y": Yq}

print(
    f"CORRECTED_ALLNODESDIRECT|prime={P}|ambient=1278|survivors=14|nodes=6|T={T}|K={K}|translated_rank=7",
    flush=True,
)

node_records = []
combined_rows = []
for record in translated["divisors"]:
    combined_rows.extend(
        [[finite(v) for v in row] for row in record["row_space_basis"]]
    )
previous_combined_rank = (
    int(matrix(finite, combined_rows).rank()) if combined_rows else 0
)
assert previous_combined_rank == 7
translated_rank = previous_combined_rank

for node_index, node_name in enumerate(NODES):
    node_start = time.perf_counter()
    chart = charts[node_name]
    edge = edges[node_name]

    surface_qq = qq_ring(sage_eval(chart["surface_equation"], locals=locals_qq))
    t_qq = qq_ring(sage_eval(chart["old_coordinate_pullback"]["t"], locals=locals_qq))
    x_qq = qq_ring(sage_eval(chart["old_coordinate_pullback"]["x"], locals=locals_qq))
    y_qq = qq_ring(sage_eval(chart["old_coordinate_pullback"]["y"], locals=locals_qq))
    old_g_qq = qq_ring(sage_eval(edge["w_cartier_equation"], locals=locals_qq))
    g_qq = old_g_qq
    for component in edge["components"]:
        name = component["name"]
        delta = NEW_TWIST[name] - OLD_TWIST[name]
        assert delta >= 0
        if delta:
            equation = {"Z": Zq, "U": Uq, "Y": Yq}[component["equation"]]
            g_qq *= equation**delta

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
    assert all(exponent >= 0 for exponent in m_exponents), (node_name, m_exponents)

    mx_value = ring.monomial(*mx)
    my_value = ring.monomial(*my)
    A = numerator_x // mx_value
    B = numerator_y // my_value
    assert A(0, 0, 0), (node_name, "A not unit")
    assert R_den(0, 0, 0), (node_name, "R_den not unit")
    assert S_den(0, 0, 0), (node_name, "S_den not unit")
    assert H_reverse(0, 0, 0), (node_name, "H_reverse not unit")

    t_exponents = common_monomial_exponents(t_value)
    assert t_exponents[2] == 0
    t_monomial = ring.monomial(*t_exponents)
    t_unit = t_value // t_monomial
    assert t_unit(0, 0, 0), (node_name, "t/m is not a unit")

    Z_LIMIT = T * int(t_exponents[0])
    U_LIMIT = T * int(t_exponents[1])

    def truncate_t_power(value):
        # Since t/t_monomial is a unit, (t^T)=(t_monomial^T) locally.
        return ring({
            monomial: coefficient
            for monomial, coefficient in ring(value).dict().items()
            if not (
                int(monomial[0]) >= Z_LIMIT
                and int(monomial[1]) >= U_LIMIT
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

    # Common clearing exponent depends only on the global ambient column,
    # not on the node. e >= T vanishes in every R/(t^T).
    active_by_family = defaultdict(list)
    skipped = 0
    for ambient_index, entry in enumerate(ambient["ambient_basis"]):
        a = int(entry["x_power"])
        b = int(entry["m_power"])
        i = int(entry["u_power"])
        k = int(entry["h_power"])
        assert k == K
        t_exponent = T + 4 * k - i - 8
        assert t_exponent >= 0
        if t_exponent >= T:
            skipped += 1
            continue
        active_by_family[(a, b)].append((ambient_index, t_exponent))

    active = 1278 - skipped
    assert skipped == 1040
    assert active == 238
    assert len(active_by_family) == 16

    # Only t^0,...,t^16 can survive.
    t_powers = [ring.one()]
    for exponent in range(1, T):
        t_powers.append(truncated_product(t_powers[-1], t_value))

    family_factor = {}
    for (a, b), terms in sorted(active_by_family.items()):
        answer = ring.one()
        for factor in (
            g,
            x_value**a,
            m_monomial**b,
            truncated_power(B, b),
            truncated_power(R_den, b),
            truncated_power(A, 8 - b),
            truncated_power(S_den, 8 - b),
            # H_reverse^(K-k)=1 because all corrected ambient columns have k=16.
        ):
            answer = truncated_product(answer, factor)
        family_factor[(a, b)] = answer

    survivor_numerators = []
    active_coefficients = 0
    for survivor_index, kernel_row in enumerate(kernel_rows):
        numerator = ring.zero()
        for family, terms in active_by_family.items():
            t_polynomial = ring.zero()
            for ambient_index, t_exponent in terms:
                coefficient = kernel_row[ambient_index]
                if coefficient:
                    active_coefficients += 1
                    t_polynomial += coefficient * t_powers[t_exponent]
            if t_polynomial:
                numerator += truncated_product(family_factor[family], t_polynomial)
        survivor_numerators.append(truncate_t_power(numerator))

    # Find the certified solved coordinate S in Y^2 = S*H, H(0) != 0.
    solved_candidates = []
    for solved_name, solved_coordinate in (("Z", Z), ("U", U)):
        quotient, remainder = (Y**2 - surface).quo_rem(solved_coordinate)
        if (
            not remainder
            and quotient(0, 0, 0)
            and quotient.degree(Y) == 0
        ):
            solved_candidates.append((solved_name, solved_coordinate, quotient))
    assert len(solved_candidates) == 1, (node_name, solved_candidates)
    solved_name, solved_coordinate, unit_h = solved_candidates[0]
    assert surface == Y**2 - solved_coordinate * unit_h

    max_y_degree = max(
        (
            int(monomial[2])
            for numerator in survivor_numerators
            for monomial in numerator.dict()
        ),
        default=0,
    )
    max_h_power = max_y_degree // 2
    h_powers = [ring.one()]
    for exponent in range(1, max_h_power + 1):
        h_powers.append(h_powers[-1] * unit_h)

    def rank_two_remainder(value):
        answer = {}
        for (z_exp, u_exp, y_exp), coefficient in ring(value).dict().items():
            q, parity = divmod(int(y_exp), 2)
            hp = h_powers[q]
            for (hz, hu, hy), hcoef in hp.dict().items():
                assert hy == 0
                z2 = int(z_exp) + int(hz) + (q if solved_name == "Z" else 0)
                u2 = int(u_exp) + int(hu) + (q if solved_name == "U" else 0)
                # Exact coefficient-ring quotient by Z^Z_LIMIT U^U_LIMIT.
                if z2 >= Z_LIMIT and u2 >= U_LIMIT:
                    continue
                key = (z2, u2, parity)
                value2 = coefficient * hcoef
                if not value2:
                    continue
                new_value = answer.get(key, finite.zero()) + value2
                if new_value:
                    answer[key] = new_value
                elif key in answer:
                    del answer[key]
        return answer

    remainders = [rank_two_remainder(value) for value in survivor_numerators]
    coordinates = sorted({
        coordinate
        for remainder in remainders
        for coordinate in remainder
    })
    coordinate_index = {coordinate: index for index, coordinate in enumerate(coordinates)}

    node_matrix = matrix(finite, len(coordinates), 14)
    for column, remainder in enumerate(remainders):
        for coordinate, coefficient in remainder.items():
            node_matrix[coordinate_index[coordinate], column] = coefficient

    node_rank = int(node_matrix.rank())
    for row in node_matrix.rows():
        combined_rows.append(list(row))
    combined_rank = (
        int(matrix(finite, combined_rows).rank())
        if combined_rows else 0
    )
    gain = combined_rank - previous_combined_rank
    previous_combined_rank = combined_rank
    remaining = 14 - combined_rank
    node_seconds = time.perf_counter() - node_start

    nonzero_survivors = sum(bool(remainder) for remainder in remainders)

    record = {
        "node": node_name,
        "solved_coordinate": solved_name,
        "surface_unit_H": str(unit_h),
        "t_monomial_exponents": {
            "Z": int(t_exponents[0]),
            "U": int(t_exponents[1]),
        },
        "t_unit_at_origin": int(t_unit(0, 0, 0)),
        "coefficient_ideal": f"(Z^{Z_LIMIT}*U^{U_LIMIT})",
        "active_ambient_columns": int(active),
        "active_kernel_coefficients": int(active_coefficients),
        "quotient_coordinate_rows": int(len(coordinates)),
        "nonzero_survivor_images": int(nonzero_survivors),
        "restricted_rank": int(node_rank),
        "combined_rank_after_node": int(combined_rank),
        "new_rank_gain": int(gain),
        "rank_beyond_translated": int(combined_rank - translated_rank),
        "remaining_dimension": int(remaining),
        "seconds": float(node_seconds),
    }
    node_records.append(record)

    print(
        "NODEDIRECT|"
        f"index={node_index+1}/6|node={node_name}|solved={solved_name}|"
        f"tmon=Z^{t_exponents[0]}U^{t_exponents[1]}|"
        f"ideal=Z^{Z_LIMIT}U^{U_LIMIT}|rows={len(coordinates)}|"
        f"nonzero_survivors={nonzero_survivors}|"
        f"rank={node_rank}|gain={gain}|combined_rank={combined_rank}|"
        f"remaining={remaining}|seconds={node_seconds:.4f}",
        flush=True,
    )

payload = {
    "schema": "elkies-k3.h92-q8-corrected1278-all-e7-nodes-direct-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_CORRECTED1278_EXACT_ALL_E7_NODE_IMAGES",
    "prime": int(P),
    "inputs": {
        "p1": {"path": str(P1_PATH.relative_to(ROOT)), "sha256": digest(P1_PATH)},
        "pullbacks": {"path": str(PULLBACKS_PATH.relative_to(ROOT)), "sha256": digest(PULLBACKS_PATH)},
        "gluing": {"path": str(GLUING_PATH.relative_to(ROOT)), "sha256": digest(GLUING_PATH)},
        "ambient": {"path": str(AMBIENT_PATH.relative_to(ROOT)), "sha256": digest(AMBIENT_PATH)},
        "global_kernel": {"path": str(KERNEL_PATH.relative_to(ROOT)), "sha256": digest(KERNEL_PATH)},
        "translated": {"path": str(TRANSLATED_PATH.relative_to(ROOT)), "sha256": digest(TRANSLATED_PATH)},
    },
    "common": {
        "ambient_dimension": 1278,
        "global_survivor_dimension": 14,
        "T": int(T),
        "K": int(K),
        "columns_zero_before_surface_reduction": 1040,
        "active_ambient_columns": 238,
        "method": (
            "Exact rank-two coefficient-module reduction via each actual "
            "surface relation Y^2=S*H, stripping only node-local unit factors "
            "from t. No Singular/local Groebner reduction."
        ),
    },
    "nodes": node_records,
    "combined": {
        "translated_rank": int(translated_rank),
        "nodes_plus_translated_rank": int(previous_combined_rank),
        "node_gain_beyond_translated": int(previous_combined_rank - translated_rank),
        "remaining_dimension": int(14 - previous_combined_rank),
    },
    "boundary": (
        "This stacks the exact modular quotient images at the six corrected E7 edge "
        "nodes after the rank-7 translated-divisor condition. It does not "
        "yet impose the distinct marked smooth point -P1, sibling-chart "
        "overlap/gluing constraints beyond node-local regularity, or certify "
        "a characteristic-zero q8 pencil."
    ),
}
OUTPUT_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "Q8CORRECTEDALLFINITENODESDIRECT|"
    f"prime={P}|global_survivor=14|translated_rank={translated_rank}|nodes=6|"
    f"combined_rank={previous_combined_rank}|"
    f"node_gain={previous_combined_rank-translated_rank}|remaining={14-previous_combined_rank}|"
    "status=EXPERIMENTAL_MODULAR_CORRECTED1278_EXACT_ALL_E7_NODE_IMAGES",
    flush=True,
)
print(f"OUTPUT|{OUTPUT_PATH}", flush=True)
