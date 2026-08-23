#!/usr/bin/env sage -python
"""
Exact two-translated-divisor q8 E7 principal-part probe on the true H92
1600-column ambient, restricted to the 18-dimensional global survivor.

The two actual resolved charts with translated t-factors are

  E7_2--E7_5:  t = Z*U*(Z+second)
  E7_3--E7_6:  t = Z^2*U^2*(Z+third)

At the corresponding edge node the translated factor is a unit, so the
previous six node-local quotient tests do not see the divisor L=0.  Here each
translated factor is promoted to a local parameter L and its principal parts
are computed exactly.

For an actual global term
    f = t^(72-d)/h_rev(t)^18 * x^a * m^b
the translated divisors have generically
    ord_L(t)=ord_L(x)=ord_L(x-xP)=ord_L(y-yP)=1,
    ord_L(g)=0,
so
    ord_L(f)=72-d+a.

A common polynomial clearing has L-order 17+9=26.  We compute the image
modulo L^26 exactly using the monic surface equation Y^2=H(L,U), i.e. as a
rank-two module over GF(p)[L,U] with basis {1,Y}.  No Groebner basis is used.

The two quotient maps are stacked as independent target summands, so the
combined matrix rank is the number of independent conditions imposed by both
translated divisors on the common 18-dimensional survivor.

Run:
  sage -python ~/Downloads/probe_h92_q8_true1600_two_translated_divisors_p59.sage

Optional:
  --repo /path/to/jacobian-research
  --prime 59
"""

import argparse
import hashlib
import json
import time
from collections import defaultdict
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, sage_eval


CHARTS_TO_TEST = ("E7_2--E7_5", "E7_3--E7_6")
MOD_L = 26


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


def common_monomial_exponents(value):
    terms = list(value.dict())
    assert terms
    return tuple(min(exp[i] for exp in terms) for i in range(3))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=59)
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
OUTPUT_PATH = GEN / f"zz-h92-q8-true1600-two-translated-divisors-mod-{P}.json"

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
assert all(name in charts and name in edges for name in CHARTS_TO_TEST)

# Exact P1 reversed functions in the E7 base t.
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


def actual_degree(entry):
    if "actual_u_power" in entry:
        return int(entry["actual_u_power"])
    return int(entry["u_power"]) + 9


def compute_chart(chart_name):
    started = time.perf_counter()
    chart = charts[chart_name]
    edge = edges[chart_name]

    # Original chart.
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

    # The monomial factor is exactly the part seen by the edge-node local
    # calculation. The quotient is the translated linear factor invisible
    # there.
    t_mono_exp = common_monomial_exponents(t_map)
    assert t_mono_exp[2] == 0
    t_monomial = R.monomial(*t_mono_exp)
    L_old = R(t_map // t_monomial)
    assert L_old.degree(Z) == 1
    assert L_old.degree(U) == 0 and L_old.degree(Y) == 0
    z_coeff = L_old.monomial_coefficient(Z)
    constant = L_old(0, 0, 0)
    assert z_coeff and constant
    assert L_old == z_coeff * Z + constant

    # Translate to L=L_old.
    S = PolynomialRing(finite, names=("L", "U", "Y"))
    L, Us, Ys = S.gens()
    Z_of_L = (L - finite(constant)) / finite(z_coeff)

    def translate(poly):
        return S(poly(Z_of_L, Us, Ys))

    surface_L = translate(surface)
    t_L = translate(t_map)
    x_L = translate(x_map)
    y_L = translate(y_map)
    g_L = translate(g_map)

    xp_n_L = S(xp_num(t_L))
    xp_d_L = S(xp_den(t_L))
    yp_n_L = S(yp_num(t_L))
    yp_d_L = S(yp_den(t_L))
    assert xp_d_L.subs({L: 0}) != 0
    assert yp_d_L.subs({L: 0}) != 0

    nx_L = S(x_L * xp_d_L - xp_n_L)
    ny_L = S(y_L * yp_d_L - yp_n_L)

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
    expected_orders = {"t": 1, "x": 1, "y": 1, "Nx": 1, "Ny": 1, "g": 0}
    assert orders == expected_orders, (chart_name, orders)

    # Strip the explicit L factors. Each residue is nonzero modulo L; hence
    # it is a unit in the height-one local ring at the generic point L=0.
    t1 = S(t_L // L)
    x1 = S(x_L // L)
    nx1 = S(nx_L // L)
    ny1 = S(ny_L // L)
    for label, value in (("t1", t1), ("x1", x1), ("nx1", nx1), ("ny1", ny1)):
        assert value.subs({L: 0}) != 0, (chart_name, label)

    # Surface is monic quadratic in Y: Y^2=H(L,U).
    H_LU = S(Ys**2 - surface_L)
    assert H_LU.degree(Ys) == 0
    assert surface_L == Ys**2 - H_LU

    # Rank-two coefficient module A + A*Y, A=GF(p)[L,U], truncated mod L^26.
    A = PolynomialRing(finite, names=("L", "U"))
    La, Ua = A.gens()
    H = A(H_LU(La, Ua, 0))

    def trunc(poly):
        poly = A(poly)
        return A({
            exp: coeff
            for exp, coeff in poly.dict().items()
            if int(exp[0]) < MOD_L
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
        poly = S(poly)
        even = A.zero()
        odd = A.zero()
        for (le, ue, ye), coeff in poly.dict().items():
            q, parity = divmod(int(ye), 2)
            if q not in hp_cache:
                hp_cache[q] = trunc(H**q)
            mon = A({(int(le), int(ue)): coeff})
            term = trunc(mon * hp_cache[q])
            if parity:
                odd = trunc(odd + term)
            else:
                even = trunc(even + term)
        return (even, odd)

    t1p = to_pair(t1)
    x1p = to_pair(x1)
    nx1p = to_pair(nx1)
    ny1p = to_pair(ny1)
    xpdp = to_pair(xp_d_L)
    ypdp = to_pair(yp_d_L)

    # All except ny may carry no Y; keep generic pair arithmetic regardless.
    assert t1p[1] == 0 and x1p[1] == 0 and nx1p[1] == 0
    assert xpdp[1] == 0 and ypdp[1] == 0

    family_factor = {}
    for a in (0, 1):
        max_b = 9 if a == 0 else 7
        for b in range(max_b + 1):
            value = ONE
            for base_value, exponent in (
                (x1p, a),
                (ny1p, b),
                (nx1p, 9 - b),
                (xpdp, b),
                (ypdp, 9 - b),
            ):
                if exponent:
                    value = mul_pair(value, pair_power(base_value, exponent))
            family_factor[(a, b)] = value

    # In the common clearing
    #   C=t^17*h_rev^18*Nx^9*yp_den^9,
    # the L order is 17+9=26. The remaining t exponent is 89-d.
    t1_powers = [ONE]
    for exponent in range(1, 18):
        t1_powers.append(mul_pair(t1_powers[-1], t1p))

    active_by_family = defaultdict(list)
    for index, entry in enumerate(basis):
        a = int(entry["x_power"])
        b = int(entry["m_power"])
        d = actual_degree(entry)
        k = int(entry["h_power"])
        assert k == 18
        pole_order = d - 72 - a
        if pole_order <= 0:
            continue
        t_exp = 89 - d
        shift = MOD_L - pole_order
        assert 0 <= t_exp <= 17
        assert shift == 9 + a + t_exp
        assert 0 <= shift < MOD_L
        active_by_family[(a, b)].append((index, t_exp, shift))

    active_indices = sorted(
        index
        for entries in active_by_family.values()
        for index, _, _ in entries
    )
    assert len(active_indices) == 278

    images = []
    for row in kernel_rows:
        image = ZERO
        for family, entries in active_by_family.items():
            family_sum = ZERO
            for index, t_exp, shift in entries:
                c = row[index]
                if not c:
                    continue
                Lshift = A({(shift, 0): finite.one()})
                term = scale_pair(c * Lshift, t1_powers[t_exp])
                family_sum = add_pair(family_sum, term)
            if family_sum != ZERO:
                image = add_pair(image, mul_pair(family_factor[family], family_sum))
        images.append(image)

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

    rank = int(M.rank())
    nonzero = sum(bool(even or odd) for even, odd in images)
    elapsed = time.perf_counter() - started

    record = {
        "chart": chart_name,
        "translated_factor": str(L_old),
        "node_monomial_exponents": {
            "Z": int(t_mono_exp[0]),
            "U": int(t_mono_exp[1]),
        },
        "orders": orders,
        "active_ambient_columns": int(len(active_indices)),
        "coordinate_rows": int(len(coords)),
        "nonzero_survivor_images": int(nonzero),
        "restricted_rank": int(rank),
        "remaining_if_alone": int(18 - rank),
        "seconds": float(elapsed),
        # Store a row-space basis only; enough to reproduce the combined rank
        # without bloating the artifact with all principal-part coordinates.
        "row_space_basis": [
            [int(v) for v in row]
            for row in M.row_space().basis()
        ],
    }
    return M, record


print(
    f"TWOTRANSLATED_START|prime={P}|ambient=1600|survivors=18|threshold=L^{MOD_L}",
    flush=True,
)

combined_rows = []
records = []
combined_rank = 0

for sequence, chart_name in enumerate(CHARTS_TO_TEST, 1):
    M, record = compute_chart(chart_name)
    previous_rank = combined_rank
    combined_rows.extend([list(row) for row in M.row_space().basis()])
    combined_rank = int(matrix(finite, combined_rows).rank()) if combined_rows else 0
    gain = combined_rank - previous_rank
    record["incremental_rank_gain"] = int(gain)
    record["combined_rank_after"] = int(combined_rank)
    record["remaining_after"] = int(18 - combined_rank)
    records.append(record)

    print(
        "TRANSLATED_RESULT|"
        f"index={sequence}/2|chart={chart_name}|L={record['translated_factor']}|"
        f"rows={record['coordinate_rows']}|"
        f"nonzero_images={record['nonzero_survivor_images']}|"
        f"rank={record['restricted_rank']}|gain={gain}|"
        f"combined_rank={combined_rank}|remaining={18-combined_rank}|"
        f"seconds={record['seconds']:.4f}",
        flush=True,
    )

    # The first divisor is now a useful p43 regression from the previous run.
    if P == 43 and chart_name == "E7_2--E7_5":
        if record["restricted_rank"] != 8:
            raise SystemExit(
                "REGRESSION_MISMATCH: E7_2--E7_5 translated divisor "
                f"rank={record['restricted_rank']}, expected 8 at p43"
            )

payload = {
    "schema": "elkies-k3.h92-q8-true1600-two-translated-divisors-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_TRUE1600_TWO_TRANSLATED_DIVISOR_IMAGES",
    "prime": int(P),
    "common": {
        "ambient_dimension": 1600,
        "global_survivor_dimension": 18,
        "principal_part_threshold": MOD_L,
        "actual_degree_used": True,
        "helper_normalization": (
            "The q6^9 helper shift i=d-9 is stripped. No extra helper t^9 "
            "is imposed in these reduced-divisor regularity tests."
        ),
        "method": (
            "Exact rank-two surface-module principal parts modulo L^26; "
            "no Groebner basis."
        ),
    },
    "divisors": records,
    "combined": {
        "restricted_rank": int(combined_rank),
        "remaining_dimension": int(18 - combined_rank),
    },
    "inputs": {
        "p1": {"path": str(P1_PATH.relative_to(ROOT)), "sha256": digest(P1_PATH)},
        "pullbacks": {"path": str(PULLBACKS_PATH.relative_to(ROOT)), "sha256": digest(PULLBACKS_PATH)},
        "gluing": {"path": str(GLUING_PATH.relative_to(ROOT)), "sha256": digest(GLUING_PATH)},
        "ambient": {"path": str(AMBIENT_PATH.relative_to(ROOT)), "sha256": digest(AMBIENT_PATH)},
        "global_kernel": {"path": str(KERNEL_PATH.relative_to(ROOT)), "sha256": digest(KERNEL_PATH)},
    },
    "boundary": (
        "These are the two translated t-divisor principal-part images on the "
        "p-modular true global survivor. A combined rank 16 would match the "
        "degree checksum, but still requires a second-prime regression and "
        "a geometric audit that these are precisely the omitted finite E7 "
        "divisor conditions before claiming h0=2 in characteristic zero."
    ),
}
OUTPUT_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "Q8TRUETWOTRANSLATED|"
    f"prime={P}|global_survivor=18|"
    f"rank1={records[0]['restricted_rank']}|"
    f"rank2={records[1]['restricted_rank']}|"
    f"gain2={records[1]['incremental_rank_gain']}|"
    f"combined_rank={combined_rank}|remaining={18-combined_rank}|"
    "status=EXPERIMENTAL_MODULAR_TRUE1600_TWO_TRANSLATED_DIVISOR_IMAGES",
    flush=True,
)
print(f"OUTPUT|{OUTPUT_PATH}", flush=True)
