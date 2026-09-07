#!/usr/bin/env sage -python
"""
H92 q8 true-envelope twist-ladder audit.

For n=-1,0,+1 compute the modular space of the currently compiled q8
pushforward E(n), including:

  * exact smooth h-collision block,
  * exact reduced E8 block + literal -11F once,
  * exact seven generic E7 residue blocks,
  * the translated/affine E7 divisor, seen in both charts and stacked.

Twisting by O_P1(n) is represented with the infinity/E7 trivialization:
  * numerator degree caps become 87+n and 89+n;
  * the E7 helper exponent is i=d-9-n (not d-9);
  * translated-divisor residual orders gain n*ord_L(t), equivalently
    t-unit exponent becomes 89+n-d.

Finite H and E8 conditions are unchanged.

If the final ladder is
    h0(E(-1))=0, h0(E)=10, h0(E(1))=28,
then because rank(E)=18 the splitting type is forced to
    O^10 + O(-1)^8
and deg(E)=-8.

Run:
  sage -python ~/Downloads/probe_h92_q8_true_twist_ladder.sage

Optional:
  --repo /path/to/jacobian-research
  --prime 43
  --twists -1,0,1
"""

import argparse
import json
import shutil
import subprocess
import time
from collections import defaultdict
from pathlib import Path

from sage.all import (
    GF, PolynomialRing, QQ, ZZ, binomial, gcd, matrix, sage_eval
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


def qq_to_finite(value, finite):
    value = QQ(value)
    den = finite(ZZ(value.denominator()))
    if not den:
        raise ZeroDivisionError(f"prime divides coefficient denominator: {value}")
    return finite(ZZ(value.numerator())) / den


def ff_polynomial(ring, finite, values):
    return ring([qq_to_finite(value, finite) for value in values])


def common_monomial_exponents(value):
    terms = list(value.dict())
    assert terms
    return tuple(min(exp[i] for exp in terms) for i in range(3))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--twists", default="-1,0,1")
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts" / "generated-results"
SCRIPTS = ROOT / "elkies-k3" / "scripts"
P = int(args.prime)
finite = GF(P)
twists = [int(x.strip()) for x in args.twists.split(",") if x.strip()]
if not twists:
    raise ValueError("no twists supplied")

sage = shutil.which("sage")
if not sage:
    raise SystemExit("sage executable not found")

P1_PATH = GEN / "elkies-k3-h92-p1-lift.json"
GENERIC_PATH = GEN / "elkies-k3-h92-q8-generic-rr-ambient.json"
PULLBACKS_PATH = GEN / "elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
GLUING_PATH = GEN / "elkies-k3-h92-q8-actual-e7-gluing.json"

for path in (P1_PATH, GENERIC_PATH, PULLBACKS_PATH, GLUING_PATH):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

p1 = json.loads(P1_PATH.read_text())
generic = json.loads(GENERIC_PATH.read_text())
pullbacks = json.loads(PULLBACKS_PATH.read_text())
gluing = json.loads(GLUING_PATH.read_text())

assert p1["status"] == "PASS_EXACT_H92_P1"
assert generic["status"] == "PASS_EXACT_Q8_GENERIC_RR_AMBIENT"
assert len(generic["basis"]) == 18
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"
assert gluing["status"] == "PASS_EXACT_Q8_ACTUAL_E7_GLUING"

charts = {entry["name"]: entry for entry in pullbacks["charts"]}
edges = {entry["name"]: entry for entry in gluing["actual_edge_chart_gluing"]}
TRANSLATED_CHARTS = ("E7_2--E7_5", "E7_3--E7_6")


def run_repo_script(name, *extra):
    command = [sage, "-python", str(SCRIPTS / name), *map(str, extra)]
    print("RUN|" + " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


# Shared finite-base data.
u_ring = PolynomialRing(finite, "u")
u = u_ring.gen()
u_field = u_ring.fraction_field()

h = ff_polynomial(
    u_ring, finite, p1["structured_denominator"]["Z4_coefficients"]
)
assert h.degree() == 4 and h(0)
assert gcd(h, h.derivative()) == 1

x_p = u_field(ff_polynomial(
    u_ring, finite, p1["x_entrance_base"]["numerator_coefficients"]
)) / u_field(ff_polynomial(
    u_ring, finite, p1["x_entrance_base"]["denominator_coefficients"]
))
y_p = u_field(ff_polynomial(
    u_ring, finite, p1["y_entrance_base"]["numerator_coefficients"]
)) / u_field(ff_polynomial(
    u_ring, finite, p1["y_entrance_base"]["denominator_coefficients"]
))


def residue_mod(value, modulus):
    value = u_field(value)
    num = u_ring(value.numerator())
    den = u_ring(value.denominator())
    if gcd(den, modulus) != 1:
        raise ZeroDivisionError(f"non-unit denominator modulo {modulus}")
    return u_ring((num * den.inverse_mod(modulus)) % modulus)


# P1 reversed functions for translated E7 charts.
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


def build_ambient(twist):
    ambient_basis = []
    for family in generic["basis"]:
        a = int(family["x_power"])
        b = int(family["m_power"])
        max_d = (87 if a == 0 else 89) + twist
        if max_d < 0:
            continue
        for d in range(max_d + 1):
            ambient_basis.append({
                "kind": family["kind"],
                "x_power": a,
                "m_power": b,
                # E(n) at E7: add n*v(t), hence helper i=d-9-n.
                "u_power": d - 9 - twist,
                "h_power": 18,
                "actual_u_power": d,
                "base_twist": twist,
                "coefficient": f"u^{d}/h(u)^18",
            })
    expected = 1600 + 18 * twist
    assert len(ambient_basis) == expected, (twist, len(ambient_basis), expected)
    return ambient_basis


def build_global_matrix(twist, ambient_basis, prefix):
    N = len(ambient_basis)
    ambient_path = GEN / f"{prefix}-ambient.json"
    cond_path = GEN / f"{prefix}-generic-conditions.json"
    r13_path = GEN / f"{prefix}-e7-1-3.json"
    r56_path = GEN / f"{prefix}-e7-5-6.json"
    r47_path = GEN / f"{prefix}-e7-4-7.json"

    ambient_payload = {
        "schema": "elkies-k3.h92-q8-twist-ladder-ambient.v1",
        "status": "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT",
        "ambient_dimension": N,
        "ambient_basis": ambient_basis,
        "base_twist": twist,
        "normalization": {
            "actual_coefficient": "u^d/h(u)^18*x^a*m^b",
            "stored_u_power": "d-9-n",
            "reason": "O_P1(n) is trivialized away from E7/infinity.",
        },
    }
    ambient_path.write_text(json.dumps(ambient_payload, indent=2, sort_keys=True) + "\n")

    run_repo_script(
        "derive_h92_q8_all_component_generic_conditions.sage",
        "--ambient", ambient_path,
        "--output", cond_path,
    )
    run_repo_script(
        "derive_h92_q8_e7_1_3_generic_residue_rows.sage",
        "--conditions", cond_path,
        "--output", r13_path,
    )
    run_repo_script(
        "derive_h92_q8_e7_5_6_generic_residue_rows.sage",
        "--conditions", cond_path,
        "--output", r56_path,
    )
    run_repo_script(
        "derive_h92_q8_e7_4_7_generic_residue_rows.sage",
        "--conditions", cond_path,
        "--output", r47_path,
    )

    # Smooth H block: independent of base twist because twist is at infinity.
    H_POLE = 27
    h_modulus = h ** H_POLE
    h_dim = h_modulus.degree()
    rho = u_field(h) * y_p / x_p
    frame_coordinates = (
        [(0, j) for j in range(10)] + [(1, j) for j in range(8)]
    )
    coord_index = {v: i for i, v in enumerate(frame_coordinates)}
    H = matrix(finite, len(frame_coordinates) * h_dim, N)

    for col, entry in enumerate(ambient_basis):
        a = int(entry["x_power"])
        b = int(entry["m_power"])
        d = int(entry["actual_u_power"])
        for j in range(b + 1):
            exponent = 2 * j - b - 18 - 2 * a
            if exponent >= 0:
                continue
            value = (
                finite(binomial(b, j))
                * u ** d
                * rho ** (b - j)
                * h ** (H_POLE + exponent)
            )
            remainder = residue_mod(value, h_modulus)
            offset = coord_index[(a, j)] * h_dim
            for degree, coefficient in enumerate(remainder.list()):
                if coefficient:
                    H[offset + degree, col] += coefficient

    # E8 finite block unchanged by twist at infinity.
    E8 = matrix(finite, 378, N)
    floor_row = 0
    for col, entry in enumerate(ambient_basis):
        a = int(entry["x_power"])
        b = int(entry["m_power"])
        d = int(entry["actual_u_power"])
        floor = (11 + 2*b) if a == 0 else (15 + 2*b)
        if d < floor:
            E8[floor_row, col] = 1
            floor_row += 1
    assert floor_row == 376

    u2 = u**2
    s_e8 = u_field(u**2) * y_p / x_p
    assert residue_mod(s_e8, u2)[0]
    for col, entry in enumerate(ambient_basis):
        if int(entry["x_power"]) != 0:
            continue
        b = int(entry["m_power"])
        d = int(entry["actual_u_power"])
        floor = 11 + 2*b
        if d < floor or d-floor >= 2:
            continue
        value = u_field(u**(d-floor)) * s_e8**b / u_field(h**18)
        remainder = residue_mod(value, u2)
        for jet in range(2):
            E8[376+jet, col] = remainder[jet]

    conditions = json.loads(cond_path.read_text())
    singleton_indices = list(map(
        int, conditions["singleton_coordinate_block"]["basis_indices"]
    ))
    specs = [("singleton", idx, None) for idx in singleton_indices]
    for path in (r13_path, r56_path, r47_path):
        payload = json.loads(path.read_text())
        for component in payload["components"]:
            for row in component["non_singleton_residue_rows"]:
                specs.append((
                    component["component"],
                    int(row["residual_order"]),
                    row["entries"],
                ))

    E7 = matrix(finite, len(specs), N)
    for row_index, spec in enumerate(specs):
        if spec[0] == "singleton":
            E7[row_index, int(spec[1])] = 1
        else:
            for item in spec[2]:
                E7[row_index, int(item["basis_index"])] += qq_to_finite(
                    item["coefficient"], finite
                )

    GLOBAL = E8.stack(H).stack(E7)
    rank = int(GLOBAL.rank())
    K = GLOBAL.right_kernel().basis_matrix()
    return {
        "N": N,
        "H_rank": int(H.rank()),
        "E8_rank": int(E8.rank()),
        "E7_rank": int(E7.rank()),
        "global_rows": int(GLOBAL.nrows()),
        "global_rank": rank,
        "kernel": K,
        "generic_h0": int(K.nrows()),
    }


def translated_matrix(twist, ambient_basis, kernel_matrix, chart_name):
    """
    Exact translated-divisor principal part on the current generic kernel.
    Generalizes the previous n=0 direct rank-two calculation to E(n).
    """
    N = len(ambient_basis)
    survivor_dim = kernel_matrix.nrows()
    kernel_rows = [list(row) for row in kernel_matrix.rows()]

    chart = charts[chart_name]
    edge = edges[chart_name]

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

    t_mono_exp = common_monomial_exponents(t_map)
    assert t_mono_exp[2] == 0
    t_monomial = R.monomial(*t_mono_exp)
    L_old = R(t_map // t_monomial)
    assert L_old.degree(Z) == 1 and L_old.degree(U) == 0 and L_old.degree(Y) == 0
    z_coeff = L_old.monomial_coefficient(Z)
    constant = L_old(0,0,0)
    assert z_coeff and constant

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
    nx_L = S(x_L * xp_d_L - xp_n_L)
    ny_L = S(y_L * yp_d_L - yp_n_L)

    def lval(poly):
        return min(int(exp[0]) for exp in poly.dict()) if poly else 10**9

    orders = {
        "t": lval(t_L), "x": lval(x_L), "y": lval(y_L),
        "Nx": lval(nx_L), "Ny": lval(ny_L), "g": lval(g_L),
    }
    assert orders == {"t":1, "x":1, "y":1, "Nx":1, "Ny":1, "g":0}, (
        chart_name, orders
    )

    t1 = S(t_L // L)
    x1 = S(x_L // L)
    nx1 = S(nx_L // L)
    ny1 = S(ny_L // L)

    H_LU = S(Ys**2 - surface_L)
    assert H_LU.degree(Ys) == 0
    assert surface_L == Ys**2 - H_LU

    A = PolynomialRing(finite, names=("L", "U"))
    La, Ua = A.gens()
    Hcoef = A(H_LU(La, Ua, 0))
    MOD_L = 26

    def trunc(poly):
        poly = A(poly)
        return A({
            exp: coeff for exp, coeff in poly.dict().items()
            if int(exp[0]) < MOD_L
        })

    ZERO = (A.zero(), A.zero())
    ONE = (A.one(), A.zero())

    def add_pair(a,b):
        return (trunc(a[0]+b[0]), trunc(a[1]+b[1]))

    def scale_pair(c,a):
        return (trunc(c*a[0]), trunc(c*a[1]))

    def mul_pair(a,b):
        return (
            trunc(a[0]*b[0] + a[1]*b[1]*Hcoef),
            trunc(a[0]*b[1] + a[1]*b[0]),
        )

    def ppow(a,n):
        out = ONE
        base = a
        while n:
            if n & 1:
                out = mul_pair(out,base)
            n //= 2
            if n:
                base = mul_pair(base,base)
        return out

    hp = {0:A.one()}
    def to_pair(poly):
        even, odd = A.zero(), A.zero()
        for (le,ue,ye), coeff in S(poly).dict().items():
            q, parity = divmod(int(ye),2)
            if q not in hp:
                hp[q] = trunc(Hcoef**q)
            term = trunc(A({(int(le),int(ue)): coeff}) * hp[q])
            if parity:
                odd = trunc(odd+term)
            else:
                even = trunc(even+term)
        return even,odd

    t1p, x1p, nx1p, ny1p = map(to_pair, (t1,x1,nx1,ny1))
    xpdp, ypdp = to_pair(xp_d_L), to_pair(yp_d_L)

    family_factor = {}
    for a in (0,1):
        max_b = 9 if a == 0 else 7
        for b in range(max_b+1):
            f = ONE
            for base, exponent in (
                (x1p,a),(ny1p,b),(nx1p,9-b),(xpdp,b),(ypdp,9-b)
            ):
                if exponent:
                    f = mul_pair(f, ppow(base,exponent))
            family_factor[(a,b)] = f

    t1_powers = [ONE]
    for e in range(1,18):
        t1_powers.append(mul_pair(t1_powers[-1], t1p))

    active_by_family = defaultdict(list)
    for index, entry in enumerate(ambient_basis):
        a = int(entry["x_power"])
        b = int(entry["m_power"])
        d = int(entry["actual_u_power"])
        # E(n): ord_L = 72-d+a+n.
        pole_order = d - 72 - a - twist
        if pole_order <= 0:
            continue
        shift = MOD_L - pole_order
        t_exp = 89 + twist - d
        assert 0 <= shift < MOD_L
        assert 0 <= t_exp <= 17
        assert shift == 9 + a + t_exp
        active_by_family[(a,b)].append((index,t_exp,shift))

    images = []
    for row in kernel_rows:
        image = ZERO
        for family, entries in active_by_family.items():
            family_sum = ZERO
            for index,t_exp,shift in entries:
                c = row[index]
                if not c:
                    continue
                Lshift = A({(shift,0): finite.one()})
                family_sum = add_pair(
                    family_sum,
                    scale_pair(c*Lshift, t1_powers[t_exp])
                )
            if family_sum != ZERO:
                image = add_pair(image, mul_pair(family_factor[family], family_sum))
        images.append(image)

    coords = sorted({
        (int(le),int(ue),parity)
        for even,odd in images
        for parity,poly in ((0,even),(1,odd))
        for (le,ue) in poly.dict()
    })
    ci = {key:i for i,key in enumerate(coords)}
    M = matrix(finite, len(coords), survivor_dim)
    for col,(even,odd) in enumerate(images):
        for parity,poly in ((0,even),(1,odd)):
            for (le,ue),coeff in poly.dict().items():
                M[ci[(int(le),int(ue),parity)],col] = coeff
    return M, str(L_old)



# ---------------------------------------------------------------------------
# Focused lower-twist multiplication diagnostic.
# ---------------------------------------------------------------------------

if P != 43:
    raise SystemExit("This focused artifact currently expects the existing p43 artifacts.")

AMBIENT0_PATH = GEN / "zz-h92-q8-true1600-ambient.json"
KERNEL0_PATH = GEN / "zz-h92-q8-true1600-global-kernel-mod-43.json"
TRANS0_PATH = GEN / "zz-h92-q8-true1600-two-translated-divisors-mod-43.json"

M1_PREFIX = "zz-h92-q8-twistladder-p43-twist-m1"
M1_AMBIENT_PATH = GEN / f"{M1_PREFIX}-ambient.json"
M1_COND_PATH = GEN / f"{M1_PREFIX}-generic-conditions.json"
M1_R13_PATH = GEN / f"{M1_PREFIX}-e7-1-3.json"
M1_R56_PATH = GEN / f"{M1_PREFIX}-e7-5-6.json"
M1_R47_PATH = GEN / f"{M1_PREFIX}-e7-4-7.json"

# Recreate cheap metadata artifacts only if they have been removed.
if not M1_AMBIENT_PATH.exists():
    m1_basis_tmp = build_ambient(-1)
    M1_AMBIENT_PATH.write_text(json.dumps({
        "schema": "elkies-k3.h92-q8-twist-ladder-ambient.v1",
        "status": "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT",
        "ambient_dimension": len(m1_basis_tmp),
        "ambient_basis": m1_basis_tmp,
        "base_twist": -1,
    }, indent=2, sort_keys=True) + "\n")
if not M1_COND_PATH.exists():
    run_repo_script(
        "derive_h92_q8_all_component_generic_conditions.sage",
        "--ambient", M1_AMBIENT_PATH, "--output", M1_COND_PATH,
    )
if not M1_R13_PATH.exists():
    run_repo_script(
        "derive_h92_q8_e7_1_3_generic_residue_rows.sage",
        "--conditions", M1_COND_PATH, "--output", M1_R13_PATH,
    )
if not M1_R56_PATH.exists():
    run_repo_script(
        "derive_h92_q8_e7_5_6_generic_residue_rows.sage",
        "--conditions", M1_COND_PATH, "--output", M1_R56_PATH,
    )
if not M1_R47_PATH.exists():
    run_repo_script(
        "derive_h92_q8_e7_4_7_generic_residue_rows.sage",
        "--conditions", M1_COND_PATH, "--output", M1_R47_PATH,
    )

for path in (
    AMBIENT0_PATH, KERNEL0_PATH, TRANS0_PATH,
    M1_AMBIENT_PATH, M1_COND_PATH, M1_R13_PATH, M1_R56_PATH, M1_R47_PATH
):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

ambient0 = json.loads(AMBIENT0_PATH.read_text())
kernel0 = json.loads(KERNEL0_PATH.read_text())
trans0 = json.loads(TRANS0_PATH.read_text())
ambient_m1 = json.loads(M1_AMBIENT_PATH.read_text())

B0 = ambient0["ambient_basis"]
Bm1 = ambient_m1["ambient_basis"]
assert len(B0) == 1600 and len(Bm1) == 1582

K0 = matrix(
    finite,
    [[finite(v) for v in row] for row in kernel0["kernel_basis_rows"]],
)
assert K0.nrows() == 18 and K0.ncols() == 1600

# The two translated charts have identical row space at p43; use their
# combined row-space anyway to avoid relying on that interpretation.
trans_rows = []
for record in trans0["divisors"]:
    trans_rows.extend(
        [[finite(v) for v in row] for row in record["row_space_basis"]]
    )
A0 = matrix(finite, trans_rows)
A0_basis = A0.row_space().basis_matrix()
assert A0_basis.rank() == 8

# Coefficients in the global-18 coordinates that survive the affine divisor.
C0 = A0_basis.right_kernel().basis_matrix()
assert C0.nrows() == 10 and C0.ncols() == 18
F0 = C0 * K0
assert F0.nrows() == 10 and F0.ncols() == 1600 and F0.rank() == 10


def key(entry):
    return (
        int(entry["x_power"]),
        int(entry["m_power"]),
        int(entry["actual_u_power"]),
    )


idx0 = {key(entry): i for i, entry in enumerate(B0)}
idxm1 = {key(entry): i for i, entry in enumerate(Bm1)}
assert len(idx0) == 1600 and len(idxm1) == 1582
assert set(idxm1).issubset(idx0)

shared0 = [idx0[key(entry)] for entry in Bm1]
excluded0 = [i for i, entry in enumerate(B0) if key(entry) not in idxm1]
assert len(excluded0) == 18

# Restrict the n=0 final 10-space to the n=-1 degree envelope.
degree_rows = F0[:, excluded0].transpose()
degree_rank = int(degree_rows.rank())

# Build the n=-1 generic E7 row matrix from the already-derived residue data.
conditions_m1 = json.loads(M1_COND_PATH.read_text())
specs = [
    ("singleton", int(idx), None)
    for idx in conditions_m1["singleton_coordinate_block"]["basis_indices"]
]
for path in (M1_R13_PATH, M1_R56_PATH, M1_R47_PATH):
    payload = json.loads(path.read_text())
    for component in payload["components"]:
        for row in component["non_singleton_residue_rows"]:
            specs.append((
                component["component"],
                int(row["residual_order"]),
                row["entries"],
            ))

R_m1 = matrix(finite, len(specs), 1582)
for row_index, spec in enumerate(specs):
    if spec[0] == "singleton":
        R_m1[row_index, int(spec[1])] = 1
    else:
        for item in spec[2]:
            R_m1[row_index, int(item["basis_index"])] += qq_to_finite(
                item["coefficient"], finite
            )

# F0 rows are basis sections. Columns of Fshared^T map 10 coordinates to
# the n=-1 ambient after the high-degree coefficients are required to vanish.
Fshared = F0[:, shared0]
generic_extra = R_m1 * Fshared.transpose()
extra = degree_rows.stack(generic_extra)
Cpre = extra.right_kernel().basis_matrix()

print(
    "LOWERTWIST_PRE|"
    f"h0_n0=10|degree_rank={degree_rank}|"
    f"extra_generic_rank={generic_extra.rank()}|"
    f"pretranslated_dim={Cpre.nrows()}",
    flush=True,
)

# Candidate n=-1 generic survivors, expressed in the actual 1582 ambient.
Fm1_pre = Cpre * Fshared
assert Fm1_pre.ncols() == 1582

# Apply the exact translated/affine divisor at n=-1 directly to this small
# candidate space.  translated_matrix accepts any row-basis kernel matrix.
T1, L1 = translated_matrix(-1, Bm1, Fm1_pre, TRANSLATED_CHARTS[0])
T2, L2 = translated_matrix(-1, Bm1, Fm1_pre, TRANSLATED_CHARTS[1])
Trows = [list(row) for row in T1.row_space().basis()]
Trows += [list(row) for row in T2.row_space().basis()]
T = matrix(finite, Trows) if Trows else matrix(finite, 0, Fm1_pre.nrows())
Ct = T.right_kernel().basis_matrix()
Fm1 = Ct * Fm1_pre
hm1 = int(Fm1.nrows())
assert hm1 == 4, f"expected known h0(-1)=4, got {hm1}"

# Embed the four lower-twist sections into the n=0 ambient in the two ways
# induced by H^0(O(1))=<1,u>.
embed_one = matrix(finite, hm1, 1600)
embed_u = matrix(finite, hm1, 1600)

for j, entry in enumerate(Bm1):
    k = key(entry)
    j0 = idx0[k]
    for r in range(hm1):
        embed_one[r, j0] += Fm1[r, j]

    a, b, d = k
    ku = (a, b, d + 1)
    assert ku in idx0
    j1 = idx0[ku]
    for r in range(hm1):
        embed_u[r, j1] += Fm1[r, j]

images = embed_one.stack(embed_u)
mul_rank = int(images.rank())
contained_rank = int(F0.stack(images).rank())
contained = int(contained_rank == F0.rank())
new_generators = int(F0.rank() - mul_rank) if contained else None

print(
    "LOWERTWIST_MULT|"
    f"prime={P}|hm1={hm1}|h0=10|"
    f"mul_rows={images.nrows()}|mul_rank={mul_rank}|"
    f"contained={contained}|"
    f"new_degree0_generators={new_generators if new_generators is not None else 'NA'}",
    flush=True,
)

# Also report the ranks of the two individual multiplication images and their
# intersection dimension.
rank_one = int(embed_one.rank())
rank_u = int(embed_u.rank())
intersection = rank_one + rank_u - mul_rank
print(
    "LOWERTWIST_STRUCTURE|"
    f"rank_1={rank_one}|rank_u={rank_u}|"
    f"intersection={intersection}|quotient_dim={new_generators if new_generators is not None else 'NA'}",
    flush=True,
)

out = GEN / "zz-h92-q8-lower-twist-multiplication-mod-43.json"
out.write_text(json.dumps({
    "schema": "elkies-k3.h92-q8-lower-twist-multiplication-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_Q8_LOWER_TWIST_MULTIPLICATION",
    "prime": P,
    "dimensions": {
        "h0_minus_1": hm1,
        "h0_0": 10,
        "multiplication_image_rank": mul_rank,
        "new_degree0_generators": new_generators,
    },
    "structure": {
        "rank_times_1": rank_one,
        "rank_times_u": rank_u,
        "intersection_dimension": intersection,
        "contained_in_n0_space": bool(contained),
    },
    "interpretation_boundary": (
        "A quotient dimension 2 identifies two minimal degree-zero generators "
        "of the currently compiled graded module. It does not by itself prove "
        "that these two generators are H^0 of the intended q8 divisor; the "
        "local-lattice normalization still has to be corrected/verified."
    ),
}, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{out}", flush=True)
