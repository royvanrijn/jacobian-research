#!/usr/bin/env sage -python
"""
Complete higher normal jets on the three H92 E7 components whose certified
normal parameter is Y:

    E7_1: chart E7_1--E7_4,  Z = e^2*A(e,s), U=s, Y=e
    E7_2: chart E7_2--E7_5,  Z=s, U = e^2*A(e,s), Y=e
    E7_3: chart E7_4--E7_3,  Z=s, U = e^2*A(e,s), Y=e

The existing generic-residue compiler uses only the leading solution A(0,s)
and the initial e-coefficient. Here A(e,s) is Hensel/Newton lifted in
GF(p)(s)[[e]] far enough to recover the ENTIRE negative normal principal
part of every true-global survivor.

For an actual ambient term
    u^d/h^16 * x^a*m^b
the corrected reduced E7 order is
    r = ord(g) + (64-d)*ord(t) + a*ord(x) + b*ord(m).
The q6^8 helper shift i=d-8 is already stripped.

For each component, let M=-min(r). Multiplying by a common UNIT clearing
t0^16 * nx0^8 * yp_den^8 and by e^M gives polynomial/power-series columns
with exponent e^(M+r). The coefficient series are computed exactly modulo
e^M in K=GF(p)(s).

Each coefficient in K is then cleared to polynomial equations in s over
GF(p), so the reported matrix rank is the actual GF(p)-linear rank on the
14-dimensional survivor, not merely rank over GF(p)(s).

The three row spaces are stacked on the already-derived translated/affine
rank-7 obstruction.

No Groebner basis is used.

Run:
  sage -python ~/Downloads/probe_h92_q8_corrected1278_y_root_full_jets.sage

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

from sage.all import (
    GF, PolynomialRing, PowerSeriesRing, QQ, ZZ, matrix, sage_eval
)


COMPONENTS = (
    ("E7_1", "E7_1--E7_4", "Z"),  # solve Z=e^2*A, generic s=U
    ("E7_2", "E7_2--E7_5", "U"),  # solve U=e^2*A, generic s=Z
    ("E7_3", "E7_4--E7_3", "U"),  # solve U=e^2*A, generic s=Z
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
AMBIENT_PATH = GEN / "zz-h92-q8-corrected1278-ambient.json"
KERNEL_PATH = GEN / f"zz-h92-q8-corrected1278-global-kernel-mod-{P}.json"
CONDITIONS_PATH = GEN / "zz-h92-q8-corrected1278-generic-conditions.json"
TRANSLATED_PATH = GEN / f"zz-h92-q8-corrected1278-two-translated-divisors-mod-{P}.json"
OUTPUT_PATH = GEN / f"zz-h92-q8-corrected1278-y-root-full-jets-mod-{P}.json"

for path in (
    P1_PATH, PULLBACKS_PATH, GLUING_PATH, AMBIENT_PATH,
    KERNEL_PATH, CONDITIONS_PATH
):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

p1 = json.loads(P1_PATH.read_text())
pullbacks = json.loads(PULLBACKS_PATH.read_text())
gluing = json.loads(GLUING_PATH.read_text())
ambient = json.loads(AMBIENT_PATH.read_text())
kernel = json.loads(KERNEL_PATH.read_text())
conditions = json.loads(CONDITIONS_PATH.read_text())

assert p1["status"] == "PASS_EXACT_H92_P1"
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"
assert gluing["status"] == "PASS_EXACT_Q8_ACTUAL_E7_GLUING"
assert int(ambient["ambient_dimension"]) == 1278
assert int(kernel["prime"]) == P
assert int(kernel["dimensions"]["ambient"]) == 1278
assert int(kernel["dimensions"]["kernel"]) == 14
assert conditions["status"] == "PASS_EXACT_Q8_ALL_COMPONENT_GENERIC_CONDITION_TEMPLATE"

basis = ambient["ambient_basis"]
kernel_rows = [[finite(v) for v in row] for row in kernel["kernel_basis_rows"]]
assert len(basis) == 1278
assert len(kernel_rows) == 14
assert all(len(row) == 1278 for row in kernel_rows)

charts = {entry["name"]: entry for entry in pullbacks["charts"]}
edges = {entry["name"]: entry for entry in gluing["actual_edge_chart_gluing"]}
condition_by_component = {
    entry["component"]: entry for entry in conditions["component_conditions"]
}

OLD_TWIST = {
    f"E7_{i+1}": value for i, value in enumerate((2,5,6,4,6,3,5))
}
NEW_TWIST = {
    f"E7_{i+1}": value for i, value in enumerate((2,6,8,5,6,4,7))
}

# Reversed P1 rational functions over GF(p)[t].
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
    return int(entry["u_power"]) + 8


def compute_component(component, chart_name, solved_name):
    started = time.perf_counter()
    chart = charts[chart_name]
    edge = edges[chart_name]

    # Actual chart polynomials over GF(p).
    R = PolynomialRing(finite, names=("Z", "U", "Y"))
    Z, U, Y = R.gens()

    QQR = PolynomialRing(QQ, names=("Z", "U", "Y"))
    Zq, Uq, Yq = QQR.gens()
    qloc = {"Z": Zq, "U": Uq, "Y": Yq}

    surface = R(QQR(sage_eval(chart["surface_equation"], locals=qloc)))
    t_poly = R(QQR(sage_eval(chart["old_coordinate_pullback"]["t"], locals=qloc)))
    x_poly = R(QQR(sage_eval(chart["old_coordinate_pullback"]["x"], locals=qloc)))
    y_poly = R(QQR(sage_eval(chart["old_coordinate_pullback"]["y"], locals=qloc)))
    g_poly = R(QQR(sage_eval(edge["w_cartier_equation"], locals=qloc)))
    for edge_component in edge["components"]:
        name = edge_component["name"]
        delta = NEW_TWIST[name] - OLD_TWIST[name]
        assert delta >= 0
        if delta:
            equation = {"Z": Z, "U": U, "Y": Y}[edge_component["equation"]]
            g_poly *= equation**delta

    solved_index = {"Z": 0, "U": 1}[solved_name]
    solved_var = (Z, U)[solved_index]
    generic_var = U if solved_name == "Z" else Z

    # Coefficient field K = GF(p)(s).
    B = PolynomialRing(finite, "s")
    s = B.gen()
    K = B.fraction_field()

    # Determine the corrected orders before choosing precision.
    expected = condition_by_component[component]["orders"]
    vt_expected = int(expected["t"])
    vx_expected = int(expected["x"])
    vm_expected = int(expected["m"])

    # g has Y as the root-component equation on all three charts.
    # Read its weighted Y-normal order from direct exponents:
    weights = (2, 0, 1) if solved_name == "Z" else (0, 2, 1)

    def weighted_order(poly):
        return min(
            int(exp[0]) * weights[0]
            + int(exp[1]) * weights[1]
            + int(exp[2]) * weights[2]
            for exp in poly.dict()
        )

    vg_expected = weighted_order(g_poly)

    residuals_expected = []
    for entry in basis:
        a = int(entry["x_power"])
        b = int(entry["m_power"])
        d = actual_degree(entry)
        residuals_expected.append(
            vg_expected
            + (64 - d) * vt_expected
            + a * vx_expected
            + b * vm_expected
        )
    min_residual = min(residuals_expected)
    M = max(0, -min_residual)
    assert M > 0

    # Extra precision is only needed internally for Newton and for divisions
    # by the known e-valuations. Final obstruction is modulo e^M.
    PREC = M + 20
    PS = PowerSeriesRing(K, "e", default_prec=PREC)
    e = PS.gen()

    def ps_constant(value):
        return PS(K(value))

    def eval_chart_poly(poly, z_value, u_value, y_value):
        result = PS.zero()
        for (ze, ue, ye), coeff in poly.dict().items():
            result += (
                K(coeff)
                * z_value**int(ze)
                * u_value**int(ue)
                * y_value**int(ye)
            )
        return PS(result)

    def shift_down(value, amount):
        value = PS(value)
        for i in range(amount):
            assert value[i] == 0
        return PS(sum(
            value[i + amount] * e**i
            for i in range(max(0, PREC - amount))
        ))

    # Initial A0(s): coefficient of e^2 in
    # F(e^2*A,s,e), linear in A.
    AP = PolynomialRing(K, "A")
    AA = AP.gen()
    leading = AP.zero()
    for (ze, ue, ye), coeff in surface.dict().items():
        if solved_name == "Z":
            e_weight = 2 * int(ze) + int(ye)
            if e_weight == 2:
                leading += K(coeff) * K(s)**int(ue) * AA**int(ze)
        else:
            e_weight = 2 * int(ue) + int(ye)
            if e_weight == 2:
                leading += K(coeff) * K(s)**int(ze) * AA**int(ue)
    assert leading.degree() == 1, (component, leading)
    A0 = -leading[0] / leading[1]
    A_series = PS(A0)

    d_surface = surface.derivative(solved_var)

    # Newton/Hensel lift A(e,s), where solved=e^2*A.
    hensel_start = time.perf_counter()
    for _ in range(8):
        if solved_name == "Z":
            z_value, u_value = e**2 * A_series, ps_constant(s)
        else:
            z_value, u_value = ps_constant(s), e**2 * A_series
        y_value = e

        F_full = eval_chart_poly(surface, z_value, u_value, y_value)
        F = shift_down(F_full, 2)
        deriv = eval_chart_poly(d_surface, z_value, u_value, y_value)
        assert deriv[0] != 0, (component, "nonunit implicit derivative")
        deriv_inv = PS(deriv**(-1))
        correction = PS(F * deriv_inv)
        A_series = PS(A_series - correction)

    if solved_name == "Z":
        z_value, u_value = e**2 * A_series, ps_constant(s)
    else:
        z_value, u_value = ps_constant(s), e**2 * A_series
    y_value = e
    surface_check = eval_chart_poly(surface, z_value, u_value, y_value)
    # The Hensel lift should solve the surface through every coefficient
    # relevant to the principal part.
    assert all(surface_check[i] == 0 for i in range(M + 10)), (
        component, "Hensel precision insufficient"
    )
    hensel_seconds = time.perf_counter() - hensel_start

    # Pull back t,x,y,g to K[[e]].
    t_value = eval_chart_poly(t_poly, z_value, u_value, y_value)
    x_value = eval_chart_poly(x_poly, z_value, u_value, y_value)
    y_value_old = eval_chart_poly(y_poly, z_value, u_value, y_value)
    g_value = eval_chart_poly(g_poly, z_value, u_value, y_value)

    def eval_t_poly(poly, t_series):
        result = PS.zero()
        for degree, coeff in poly.dict().items():
            result += K(coeff) * t_series**int(degree)
        return PS(result)

    xp_n = eval_t_poly(xp_num, t_value)
    xp_d = eval_t_poly(xp_den, t_value)
    yp_n = eval_t_poly(yp_num, t_value)
    yp_d = eval_t_poly(yp_den, t_value)
    assert xp_d[0] != 0 and yp_d[0] != 0

    xp_value = PS(xp_n * PS(xp_d**(-1)))
    yp_value = PS(yp_n * PS(yp_d**(-1)))
    nx = PS(x_value - xp_value)
    ny = PS(y_value_old - yp_value)

    vt = int(t_value.valuation())
    vx = int(x_value.valuation())
    vnx = int(nx.valuation())
    vny = int(ny.valuation())
    vg = int(g_value.valuation())
    vm = vny - vnx

    assert (vt, vx, vm) == (vt_expected, vx_expected, vm_expected), (
        component, (vt, vx, vm), (vt_expected, vx_expected, vm_expected)
    )
    assert vg == vg_expected, (component, vg, vg_expected)

    def unit_part(value, order):
        return shift_down(value, order)

    t0 = unit_part(t_value, vt)
    x0 = unit_part(x_value, vx)
    nx0 = unit_part(nx, vnx)
    ny0 = unit_part(ny, vny)
    g0 = unit_part(g_value, vg)

    assert t0[0] != 0
    assert x0[0] != 0
    assert nx0[0] != 0
    assert ny0[0] != 0
    assert g0[0] != 0

    # Recompute residuals from the actually lifted series orders.
    residuals = []
    for entry in basis:
        a = int(entry["x_power"])
        b = int(entry["m_power"])
        d = actual_degree(entry)
        residuals.append(
            vg + (64 - d) * vt + a * vx + b * vm
        )
    assert residuals == residuals_expected
    assert min(residuals) == min_residual

    # All unit factors are ordinary truncated K[[e]] series. The common unit
    # clearing t0^17*nx0^9*yp_d^9 avoids every negative unit exponent.
    family_factor = {}
    for a in (0, 1):
        max_b = 8 if a == 0 else 6
        for b in range(max_b + 1):
            family_factor[(a, b)] = PS(
                g0
                * x0**a
                * ny0**b
                * nx0**(8 - b)
                * xp_d**b
                * yp_d**(8 - b)
            )

    t0_powers = [PS.one()]
    for exponent in range(1, 17):
        t0_powers.append(PS(t0_powers[-1] * t0))

    active_by_family = defaultdict(list)
    for index, entry in enumerate(basis):
        r = residuals[index]
        if r >= 0:
            continue
        a = int(entry["x_power"])
        b = int(entry["m_power"])
        d = actual_degree(entry)
        shift = M + r
        t_exp = 80 - d
        assert 0 <= shift < M
        assert 0 <= t_exp <= 16
        active_by_family[(a, b)].append((index, shift, t_exp, r))

    active_count = sum(len(rows) for rows in active_by_family.values())

    # Build all 14 obstruction series modulo e^M.
    build_start = time.perf_counter()
    images = []
    for kernel_row in kernel_rows:
        image = PS.zero()
        for family, entries in active_by_family.items():
            family_sum = PS.zero()
            for index, shift, t_exp, r in entries:
                c = kernel_row[index]
                if not c:
                    continue
                family_sum += K(c) * e**shift * t0_powers[t_exp]
            if family_sum:
                image += family_factor[family] * family_sum
        # Explicitly discard coefficients e^M and above.
        image = PS(sum(image[i] * e**i for i in range(M)))
        images.append(image)
    build_seconds = time.perf_counter() - build_start

    # Convert K(s)-valued e coefficients into honest GF(p)-rows by clearing
    # one common denominator at each e order and expanding in s.
    rows = []
    Bpoly = B
    for e_order in range(M):
        coeffs = [K(image[e_order]) for image in images]
        if not any(coeffs):
            continue
        denominator = Bpoly.one()
        for value in coeffs:
            if value:
                denominator = denominator.lcm(Bpoly(value.denominator()))
        polynomial_columns = []
        max_s_degree = -1
        for value in coeffs:
            if not value:
                poly = Bpoly.zero()
            else:
                den = Bpoly(value.denominator())
                q, rem = denominator.quo_rem(den)
                assert not rem
                poly = Bpoly(value.numerator()) * q
            polynomial_columns.append(poly)
            if poly:
                max_s_degree = max(max_s_degree, int(poly.degree()))
        for s_degree in range(max_s_degree + 1):
            row = [finite(poly[s_degree]) for poly in polynomial_columns]
            if any(row):
                rows.append(row)

    Mmat = matrix(finite, rows) if rows else matrix(finite, 0, 14)
    rank = int(Mmat.rank())
    nonzero = sum(bool(image) for image in images)
    elapsed = time.perf_counter() - started

    return Mmat, {
        "component": component,
        "chart": chart_name,
        "solved_coordinate": solved_name,
        "orders": {
            "t": vt, "x": vx, "Nx": vnx, "Ny": vny,
            "m": vm, "g": vg,
        },
        "worst_residual_order": int(min_residual),
        "principal_part_modulus": f"e^{M}",
        "hensel_precision": int(PREC),
        "hensel_seconds": float(hensel_seconds),
        "active_ambient_columns": int(active_count),
        "gf43_rows": int(Mmat.nrows()),
        "nonzero_survivor_images": int(nonzero),
        "restricted_rank": int(rank),
        "row_space_basis": [
            [int(v) for v in row] for row in Mmat.row_space().basis()
        ],
        "build_seconds": float(build_seconds),
        "total_seconds": float(elapsed),
    }


print(
    f"CORRECTED_YROOTFULLJETS_START|prime={P}|ambient=1278|survivors=14|components=3",
    flush=True,
)

# Seed the combined row space with the already-derived translated/affine
# divisor obstruction.
combined_rows = []
translated_rank = 0
if TRANSLATED_PATH.exists():
    translated = json.loads(TRANSLATED_PATH.read_text())
    assert int(translated["prime"]) == P
    translated_basis = translated["divisors"][0]["row_space_basis"]
    combined_rows.extend([[finite(v) for v in row] for row in translated_basis])
    translated_rank = int(matrix(finite, combined_rows).rank()) if combined_rows else 0
    if P == 43:
        assert translated_rank == 7

combined_rank = translated_rank
records = []

for seq, (component, chart_name, solved_name) in enumerate(COMPONENTS, 1):
    Mmat, record = compute_component(component, chart_name, solved_name)
    previous = combined_rank
    combined_rows.extend([list(row) for row in Mmat.row_space().basis()])
    combined_rank = int(matrix(finite, combined_rows).rank()) if combined_rows else 0
    gain = combined_rank - previous
    record["incremental_gain"] = int(gain)
    record["combined_rank_after"] = int(combined_rank)
    record["remaining_after"] = int(14 - combined_rank)
    records.append(record)

    print(
        "YROOT_RESULT|"
        f"index={seq}/3|component={component}|chart={chart_name}|"
        f"solve={solved_name}|"
        f"orders=t:{record['orders']['t']},x:{record['orders']['x']},"
        f"m:{record['orders']['m']},g:{record['orders']['g']}|"
        f"worst={record['worst_residual_order']}|"
        f"active={record['active_ambient_columns']}|rows={record['gf43_rows']}|"
        f"nonzero_images={record['nonzero_survivor_images']}|"
        f"rank={record['restricted_rank']}|gain={gain}|"
        f"combined_rank={combined_rank}|remaining={14-combined_rank}|"
        f"hensel_seconds={record['hensel_seconds']:.4f}|"
        f"seconds={record['total_seconds']:.4f}",
        flush=True,
    )

payload = {
    "schema": "elkies-k3.h92-q8-corrected1278-y-root-full-jets-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_CORRECTED1278_Y_ROOT_FULL_NORMAL_JETS",
    "prime": int(P),
    "common": {
        "ambient_dimension": 1278,
        "global_survivor_dimension": 14,
        "translated_start_rank": int(translated_rank),
        "helper_shift": "i=d-8 already stripped",
        "method": (
            "Hensel lift of the actual implicit Y-normal branch in "
            "GF(p)(s)[[e]], complete negative principal parts, then exact "
            "coefficient clearing back to GF(p)-linear rows. No Groebner basis."
        ),
    },
    "components": records,
    "combined": {
        "restricted_rank": int(combined_rank),
        "remaining_dimension": int(14 - combined_rank),
    },
    "inputs": {
        "p1": {"path": str(P1_PATH.relative_to(ROOT)), "sha256": digest(P1_PATH)},
        "pullbacks": {"path": str(PULLBACKS_PATH.relative_to(ROOT)), "sha256": digest(PULLBACKS_PATH)},
        "gluing": {"path": str(GLUING_PATH.relative_to(ROOT)), "sha256": digest(GLUING_PATH)},
        "ambient": {"path": str(AMBIENT_PATH.relative_to(ROOT)), "sha256": digest(AMBIENT_PATH)},
        "global_kernel": {"path": str(KERNEL_PATH.relative_to(ROOT)), "sha256": digest(KERNEL_PATH)},
        "generic_conditions": {"path": str(CONDITIONS_PATH.relative_to(ROOT)), "sha256": digest(CONDITIONS_PATH)},
        "translated": (
            {"path": str(TRANSLATED_PATH.relative_to(ROOT)), "sha256": digest(TRANSLATED_PATH)}
            if TRANSLATED_PATH.exists() else None
        ),
    },
    "boundary": (
        "Together with the previous Z/U-normal sweep, this exhausts higher "
        "normal jets on all seven exceptional/root components. A combined "
        "rank below 12 would mean the remaining finite E7 information is not "
        "a missed root-component normal jet and must be sought in the marked "
        "horizontal/cover-gluing layer or in the determinant checksum."
    ),
}
OUTPUT_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "Q8CORRECTEDYROOTFULLJETS|"
    f"prime={P}|global_survivor=14|translated_rank={translated_rank}|"
    f"combined_rank={combined_rank}|remaining={14-combined_rank}|"
    "status=EXPERIMENTAL_MODULAR_CORRECTED1278_Y_ROOT_FULL_NORMAL_JETS",
    flush=True,
)
print(f"OUTPUT|{OUTPUT_PATH}", flush=True)
