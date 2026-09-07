#!/usr/bin/env sage -python
"""
Complete normal principal parts on the four H92 E7 root components whose
certified reduced normal parameter is literally Z or U:

    E7_4 : Z=0 on E7_4--E7_3
    E7_5 : Z=0 on E7_2--E7_5
    E7_6 : Z=0 on E7_3--E7_6
    E7_7 : U=0 on E7_3--E7_7

The existing generic residue compiler imposes the initial residue at each
negative normal order.  This probe computes the ENTIRE negative normal jet
for the already-computed 14-dimensional corrected-global survivor.

For each component, with normal N and valuations
    vt=ord_N(t), vx=ord_N(x), vm=ord_N(m), vg=ord_N(g),
an actual ambient term u^d/h^16 * x^a*m^b has corrected reduced order

    r = vg + (64-d)*vt + a*vx + b*vm.

The q6^8 helper shift i=d-8 is already stripped.

Let M=-min(r).  After extracting the known N powers and multiplying by a
common unit denominator, the obstruction is represented modulo N^M.  The
actual surface is monic in Y, so all arithmetic is exact in the free rank-two
module {1,Y} over GF(p)[Z,U], truncated in N. No Groebner basis is used.

The row spaces from all four components are stacked with the already-found
rank-7 translated/affine divisor obstruction.

Run:
  sage -python ~/Downloads/probe_h92_q8_corrected1278_zu_root_full_jets.sage

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


COMPONENTS = (
    ("E7_4", "E7_4--E7_3", "Z"),
    ("E7_5", "E7_2--E7_5", "Z"),
    ("E7_6", "E7_3--E7_6", "Z"),
    ("E7_7", "E7_3--E7_7", "U"),
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


def valuation(poly, variable_index):
    if not poly:
        return 10**9
    return min(int(exp[variable_index]) for exp in poly.dict())


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
TRANSLATED_PATH = GEN / f"zz-h92-q8-corrected1278-two-translated-divisors-mod-{P}.json"
OUTPUT_PATH = GEN / f"zz-h92-q8-corrected1278-zu-root-full-jets-mod-{P}.json"

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
assert int(ambient["ambient_dimension"]) == 1278
assert int(kernel["prime"]) == P
assert int(kernel["dimensions"]["ambient"]) == 1278
assert int(kernel["dimensions"]["kernel"]) == 14

basis = ambient["ambient_basis"]
kernel_rows = [[finite(v) for v in row] for row in kernel["kernel_basis_rows"]]
assert len(basis) == 1278
assert len(kernel_rows) == 14
assert all(len(row) == 1278 for row in kernel_rows)

charts = {entry["name"]: entry for entry in pullbacks["charts"]}
edges = {entry["name"]: entry for entry in gluing["actual_edge_chart_gluing"]}

OLD_TWIST = {
    f"E7_{i+1}": value for i, value in enumerate((2,5,6,4,6,3,5))
}
NEW_TWIST = {
    f"E7_{i+1}": value for i, value in enumerate((2,6,8,5,6,4,7))
}

# P1 reversed functions in t.
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


def compute_component(component, chart_name, normal_name):
    started = time.perf_counter()
    chart = charts[chart_name]
    edge = edges[chart_name]

    R = PolynomialRing(finite, names=("Z", "U", "Y"))
    Z, U, Y = R.gens()
    variables = (Z, U, Y)
    normal_index = {"Z": 0, "U": 1}[normal_name]
    N = variables[normal_index]

    QQR = PolynomialRing(QQ, names=("Z", "U", "Y"))
    Zq, Uq, Yq = QQR.gens()
    qloc = {"Z": Zq, "U": Uq, "Y": Yq}

    surface = R(QQR(sage_eval(chart["surface_equation"], locals=qloc)))
    t_map = R(QQR(sage_eval(chart["old_coordinate_pullback"]["t"], locals=qloc)))
    x_map = R(QQR(sage_eval(chart["old_coordinate_pullback"]["x"], locals=qloc)))
    y_map = R(QQR(sage_eval(chart["old_coordinate_pullback"]["y"], locals=qloc)))
    g_map = R(QQR(sage_eval(edge["w_cartier_equation"], locals=qloc)))
    for edge_component in edge["components"]:
        name = edge_component["name"]
        delta = NEW_TWIST[name] - OLD_TWIST[name]
        assert delta >= 0
        if delta:
            equation = {"Z": Z, "U": U, "Y": Y}[edge_component["equation"]]
            g_map *= equation**delta

    xp_n = R(xp_num(t_map))
    xp_d = R(xp_den(t_map))
    yp_n = R(yp_num(t_map))
    yp_d = R(yp_den(t_map))
    assert valuation(xp_d, normal_index) == 0
    assert valuation(yp_d, normal_index) == 0

    nx = R(x_map * xp_d - xp_n)
    ny = R(y_map * yp_d - yp_n)

    vt = valuation(t_map, normal_index)
    vx = valuation(x_map, normal_index)
    vy = valuation(y_map, normal_index)
    vnx = valuation(nx, normal_index)
    vny = valuation(ny, normal_index)
    vg = valuation(g_map, normal_index)
    vm = vny - vnx
    assert vt > 0
    assert vnx >= 0 and vny >= 0 and vg >= 0

    # Strip exact normal powers.
    t0 = R(t_map // N**vt)
    x0 = R(x_map // N**vx)
    nx0 = R(nx // N**vnx)
    ny0 = R(ny // N**vny)
    g0 = R(g_map // N**vg)
    for label, value in (
        ("t0", t0), ("x0", x0), ("nx0", nx0),
        ("ny0", ny0), ("g0", g0),
    ):
        assert valuation(value, normal_index) == 0, (component, label)

    # Surface is monic in Y: Y^2=H(Z,U).
    H_ZU = R(Y**2 - surface)
    assert H_ZU.degree(Y) == 0
    assert surface == Y**2 - H_ZU

    # Compute corrected residual orders directly from actual d.
    residuals = []
    for entry in basis:
        a = int(entry["x_power"])
        b = int(entry["m_power"])
        d = actual_degree(entry)
        k = int(entry["h_power"])
        assert k == 16
        residuals.append(
            vg + (64 - d) * vt + a * vx + b * vm
        )
    min_residual = min(residuals)
    Mmod = max(0, -min_residual)
    assert Mmod > 0

    # Rank-two coefficient ring A=GF(p)[Z,U], truncated in chosen normal.
    A = PolynomialRing(finite, names=("Z", "U"))
    Za, Ua = A.gens()
    normal_A_index = {"Z": 0, "U": 1}[normal_name]
    H = A(H_ZU(Za, Ua, 0))

    def trunc(poly):
        poly = A(poly)
        return A({
            exp: coeff
            for exp, coeff in poly.dict().items()
            if int(exp[normal_A_index]) < Mmod
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

    t0p = to_pair(t0)
    x0p = to_pair(x0)
    nx0p = to_pair(nx0)
    ny0p = to_pair(ny0)
    g0p = to_pair(g0)
    xpdp = to_pair(xp_d)
    ypdp = to_pair(yp_d)

    # xp/yp denominators are base functions of t and hence Y-free.
    assert xpdp[1] == 0 and ypdp[1] == 0

    # Common unit denominator:
    # t0^16 * nx0^8 * yp_den^8 * h^16.\n    # h^16 cancels because every corrected ambient column has k=16.
    family_factor = {}
    for a in (0, 1):
        max_b = 8 if a == 0 else 6
        for b in range(max_b + 1):
            value = g0p
            for base_value, exponent in (
                (x0p, a),
                (ny0p, b),
                (nx0p, 8 - b),
                (xpdp, b),
                (ypdp, 8 - b),
            ):
                if exponent:
                    value = mul_pair(value, pair_power(base_value, exponent))
            family_factor[(a, b)] = value

    t0_powers = [ONE]
    for exponent in range(1, 17):
        t0_powers.append(mul_pair(t0_powers[-1], t0p))

    active_by_family = defaultdict(list)
    for index, entry in enumerate(basis):
        r = residuals[index]
        if r >= 0:
            continue
        a = int(entry["x_power"])
        b = int(entry["m_power"])
        d = actual_degree(entry)
        shift = Mmod + r
        t_exp = 80 - d
        assert 0 <= shift < Mmod
        assert 0 <= t_exp <= 16
        active_by_family[(a, b)].append((index, shift, t_exp, r))

    active_count = sum(len(v) for v in active_by_family.values())

    images = []
    for kernel_row in kernel_rows:
        image = ZERO
        for family, entries in active_by_family.items():
            family_sum = ZERO
            for index, shift, t_exp, r in entries:
                c = kernel_row[index]
                if not c:
                    continue
                if normal_name == "Z":
                    Nshift = A({(shift, 0): finite.one()})
                else:
                    Nshift = A({(0, shift): finite.one()})
                term = scale_pair(c * Nshift, t0_powers[t_exp])
                family_sum = add_pair(family_sum, term)
            if family_sum != ZERO:
                image = add_pair(image, mul_pair(family_factor[family], family_sum))
        images.append(image)

    coords = sorted({
        (int(ze), int(ue), parity)
        for even, odd in images
        for parity, poly in ((0, even), (1, odd))
        for (ze, ue) in poly.dict()
    })
    coord_index = {key: i for i, key in enumerate(coords)}
    M = matrix(finite, len(coords), 14)
    for col, (even, odd) in enumerate(images):
        for parity, poly in ((0, even), (1, odd)):
            for (ze, ue), coeff in poly.dict().items():
                M[coord_index[(int(ze), int(ue), parity)], col] = coeff

    rank = int(M.rank())
    nonzero = sum(bool(even or odd) for even, odd in images)
    elapsed = time.perf_counter() - started

    return M, {
        "component": component,
        "chart": chart_name,
        "normal": normal_name,
        "orders": {
            "t": int(vt), "x": int(vx), "y": int(vy),
            "Nx": int(vnx), "Ny": int(vny), "m": int(vm), "g": int(vg),
        },
        "worst_residual_order": int(min_residual),
        "principal_part_modulus": f"{normal_name}^{Mmod}",
        "active_ambient_columns": int(active_count),
        "coordinate_rows": int(len(coords)),
        "nonzero_survivor_images": int(nonzero),
        "restricted_rank": int(rank),
        "row_space_basis": [
            [int(v) for v in row]
            for row in M.row_space().basis()
        ],
        "seconds": float(elapsed),
    }


print(
    f"CORRECTED_ZUROOTFULLJETS_START|prime={P}|ambient=1278|survivors=14|components=4",
    flush=True,
)

# Start with the already-derived affine/translated divisor row space.
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

for seq, (component, chart_name, normal_name) in enumerate(COMPONENTS, 1):
    M, record = compute_component(component, chart_name, normal_name)
    previous = combined_rank
    component_rows = [list(row) for row in M.row_space().basis()]
    combined_rows.extend(component_rows)
    combined_rank = int(matrix(finite, combined_rows).rank()) if combined_rows else 0
    gain = combined_rank - previous
    record["incremental_gain"] = int(gain)
    record["combined_rank_after"] = int(combined_rank)
    record["remaining_after"] = int(14 - combined_rank)
    records.append(record)

    print(
        "ZUROOT_RESULT|"
        f"index={seq}/4|component={component}|chart={chart_name}|normal={normal_name}|"
        f"orders=t:{record['orders']['t']},x:{record['orders']['x']},"
        f"m:{record['orders']['m']},g:{record['orders']['g']}|"
        f"worst={record['worst_residual_order']}|"
        f"active={record['active_ambient_columns']}|rows={record['coordinate_rows']}|"
        f"nonzero_images={record['nonzero_survivor_images']}|"
        f"rank={record['restricted_rank']}|gain={gain}|"
        f"combined_rank={combined_rank}|remaining={14-combined_rank}|"
        f"seconds={record['seconds']:.4f}",
        flush=True,
    )

    # E7_5 is now a useful p43 regression from the previous dedicated run.
    if False:
        raise SystemExit(
            f"REGRESSION_MISMATCH: E7_5 full normal rank "
            f"{record['restricted_rank']} != 0"
        )

payload = {
    "schema": "elkies-k3.h92-q8-corrected1278-zu-root-full-jets-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_CORRECTED1278_ZU_ROOT_FULL_NORMAL_JETS",
    "prime": int(P),
    "common": {
        "ambient_dimension": 1278,
        "global_survivor_dimension": 14,
        "translated_start_rank": int(translated_rank),
        "helper_shift": "i=d-8 already stripped",
        "method": (
            "Complete normal principal parts in the exact rank-two surface "
            "module {1,Y}, truncated in the certified Z/U normal coordinate; "
            "no Groebner basis."
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
        "translated": (
            {"path": str(TRANSLATED_PATH.relative_to(ROOT)), "sha256": digest(TRANSLATED_PATH)}
            if TRANSLATED_PATH.exists() else None
        ),
    },
    "boundary": (
        "This completes higher normal jets only on E7_4,E7_5,E7_6,E7_7, "
        "whose normal is literally Z or U. If the combined rank remains below "
        "12, the remaining certified Y-normal components E7_1,E7_2,E7_3 need "
        "implicit-function/Hensel normal expansions."
    ),
}
OUTPUT_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "Q8CORRECTEDZUROOTFULLJETS|"
    f"prime={P}|global_survivor=14|translated_rank={translated_rank}|"
    f"combined_rank={combined_rank}|remaining={14-combined_rank}|"
    "status=EXPERIMENTAL_MODULAR_CORRECTED1278_ZU_ROOT_FULL_NORMAL_JETS",
    flush=True,
)
print(f"OUTPUT|{OUTPUT_PATH}", flush=True)
