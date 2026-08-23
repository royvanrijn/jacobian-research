#!/usr/bin/env sage -python
"""
Characteristic-zero reconstruction of the corrected H92 q6-child q8 pencil.

This is the QQ counterpart of the two-prime corrected2cover modular probe.
It reconstructs the corrected q8 marked section directly from the exact q6
child points, fixes the binary-quartic 2-covering factor and the missing Dx
in the q-normalizer, derives the exact 13-column global RR intersection over
QQ, and if h0=2 eliminates the resulting pencil to a binary quartic over
QQ(U).  Its Jacobian is then minimized/classified with the repository's
elliptic-neighbor compiler.

Expected endpoint:
    exact RR kernel dimension = 2
    branch quartic degree     = 4
    child root lattice        = D13
    child MW rank (rho=19)    = 4

Run:
  sage -python ~/Downloads/derive_h92_q6_child_q8_corrected2cover_qq.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    EllipticCurve, PolynomialRing, QQ, ZZ, gcd, lcm, matrix, vector
)


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


def polynomial(ring, values):
    return ring([QQ(v) for v in values])


def rational(field, ring, data, nk, dk):
    return field(polynomial(ring, data[nk])) / field(polynomial(ring, data[dk]))


def point_data(point, ring):
    x, y = point.xy()
    return {
        "x_numerator_coefficients_low_to_high":
            [str(v) for v in ring(x.numerator()).list()],
        "x_denominator_coefficients_low_to_high":
            [str(v) for v in ring(x.denominator()).list()],
        "y_numerator_coefficients_low_to_high":
            [str(v) for v in ring(y.numerator()).list()],
        "y_denominator_coefficients_low_to_high":
            [str(v) for v in ring(y.denominator()).list()],
    }


def monic_power_root(value, exponent):
    root = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0, (factor, multiplicity, exponent)
        root *= factor.monic() ** (multiplicity // exponent)
    return root.monic()


def degree_or_minus_one(value):
    return -1 if not value else int(value.degree())


def infinity_order(value, ring):
    if not value:
        return 10**9
    return int(
        ring(value.denominator()).degree()
        - ring(value.numerator()).degree()
    )


def common_denominator(values, ring):
    result = ring.one()
    for value in values:
        if value:
            result = result.lcm(ring(value.denominator()))
    return result.monic()


def primitive_integer_row(row):
    dens = [QQ(value).denominator() for value in row]
    scale = lcm(dens) if dens else ZZ.one()
    values = [ZZ(QQ(value) * scale) for value in row]
    content = gcd(values)
    if content:
        values = [value // content for value in values]
    first = next((value for value in values if value), ZZ.zero())
    if first < 0:
        values = [-value for value in values]
    return vector(ZZ, values)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--component-target", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts" / "generated-results"
LOCAL = ROOT / "artifacts" / "local" / "elkies-k3"
CORE = ROOT / "elkies-k3" / "scripts" / "elliptic_neighbor_compiler.sage"

CHILD = GEN / "elkies-k3-h92-q6-child-jacobian.json"
ZERO = GEN / "elkies-k3-h92-q6-child-zero-section.json"
COMPONENTS = GEN / "elkies-k3-h92-q6-child-e7-infinity-sections.json"
TARGET = (
    args.component_target.resolve()
    if args.component_target
    else LOCAL / "q8-target-component-nef.json"
)
OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q8-corrected2cover-qq-child.json"
)

for path in (CORE, CHILD, ZERO, COMPONENTS, TARGET):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

# Reuse only generic, independently verified compiler routines:
# binary-quartic squarefree reduction/invariants and exact short-Weierstrass
# fibre classification.
exec(compile(CORE.read_text(), str(CORE), "exec"))

child = json.loads(CHILD.read_text())
zero = json.loads(ZERO.read_text())
components = json.loads(COMPONENTS.read_text())
target = json.loads(TARGET.read_text())

assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero["status"] == "PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"
assert components["status"] == "PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT"
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert target["normalization"]["representative"] == "component-nef"

# ---------------------------------------------------------------------------
# 1. Reconstruct the corrected marked section on the exact q6 child.
# ---------------------------------------------------------------------------

R = PolynomialRing(QQ, "T")
T = R.gen()
K = R.fraction_field()

model = child["minimal_short_weierstrass"]
A = polynomial(R, model["A_coefficients_low_to_high"])
Bcurve = polynomial(R, model["B_coefficients_low_to_high"])
Delta = polynomial(R, model["Delta_coefficients_low_to_high"])
assert A.degree() <= 8 and Bcurve.degree() <= 12 and Delta.degree() <= 24

E = EllipticCurve(K, [0, 0, 0, K(A), K(Bcurve)])

zdata = zero["section"]
P0 = E(
    rational(K, R, zdata,
             "x_numerator_coefficients_low_to_high",
             "x_denominator_coefficients_low_to_high"),
    rational(K, R, zdata,
             "y_numerator_coefficients_low_to_high",
             "y_denominator_coefficients_low_to_high"),
)
points = {}
for entry in components["sections"]:
    points[entry["sign"]] = E(
        rational(K, R, entry,
                 "x_numerator_coefficients_low_to_high",
                 "x_denominator_coefficients_low_to_high"),
        rational(K, R, entry,
                 "y_numerator_coefficients_low_to_high",
                 "y_denominator_coefficients_low_to_high"),
    )
affine = points[components["source"]["affine_E7_sign"]]
e77 = points[components["source"]["E7_7_sign"]]

Pmap = e77 - P0
Qmap = e77 - affine
S = Pmap + Qmap
assert not S.is_zero()
sx, sy = S.xy()
assert sy**2 == sx**3 + K(A)*sx + K(Bcurve)

nx, dx = R(sx.numerator()), R(sx.denominator())
ny, dy = R(sy.numerator()), R(sy.denominator())
h = monic_power_root(dx, 2)
assert h == monic_power_root(dy, 3)
assert h.degree() == 10
assert dx // h**2 in QQ and dy // h**3 in QQ
assert nx.gcd(h) in QQ and ny.gcd(h) in QQ and h.gcd(Delta) in QQ
assert (nx.degree()-dx.degree(), ny.degree()-dy.degree()) == (4, 6)

# Exact O.S=10 and height 24; both additive fibres are identity-component.
zS = -sx/sy
assert R((zS/K(h)).numerator()).gcd(h) in QQ
assert R((zS/K(h)).denominator()).gcd(h) in QQ
assert QQ(4 + 2*h.degree()) == 24

# ---------------------------------------------------------------------------
# 2. Exact component-nef finite module + globally regular q frame.
# ---------------------------------------------------------------------------

selected = target["selected_q8"]
assert int(selected["vertical_fibre_coefficient"]) == -2
assert list(map(int, selected["E8"]["vertical_cycle"])) == [-4,-5,-7,-10,-8,-6,-4,-2]
assert list(map(int, selected["E6"]["vertical_cycle"])) == [-2,-3,-4,-3,-2,-2]
assert int(selected["E8"]["affine_component_degree"]) == 0
assert int(selected["E6"]["affine_component_degree"]) == 0

ii = R(next(
    item for item in child["finite_fibres"] if item["kodaira"] == "II*"
)["factor"]).monic()
iv = R(next(
    item for item in child["finite_fibres"] if item["kodaira"] == "IV*"
)["factor"]).monic()
assert ii.degree() == iv.degree() == 1 and ii.gcd(iv) in QQ
assert h.gcd(ii*iv) in QQ
M = (ii**2 * iv**2).monic()
assert M.degree() == 4 and nx.gcd(M) in QQ

# q=(m-p)/h and the CORRECT Nx-pole cancellation:
#     Rnorm*h*Dy == Ny*Dx mod Nx.
normalizer = (ny*dx*(h*dy).inverse_mod(nx)).mod(nx)
assert (normalizer*h*dy - ny*dx) % nx == 0
assert normalizer.degree() < nx.degree()

p_fun = -sy/sx
alpha = -p_fun/K(h) - K(normalizer)/K(nx)
beta = K(T**2)/K(h)
assert R(alpha.denominator()).gcd(nx) in QQ

# q is zero in both (u^2,X,Y) quotients.  Hence qreg has residue -R/Nx,
# and C+B*qreg is finite-allowed iff C-B*R/Nx == 0 mod M.
rho = (normalizer * nx.inverse_mod(M)).mod(M)
assert rho.degree() < 4

# Convenient exact finite generators:
#   g1 = qreg + rho
#   g2 = M.
abase = alpha + K(rho)
bbase = beta
tbase = K(M)

assert infinity_order(bbase, R) == 8
assert infinity_order(tbase, R) == -4
abase_order = infinity_order(abase, R)

# ---------------------------------------------------------------------------
# 3. Complete exact infinity intersection over QQ.
# ---------------------------------------------------------------------------

required = 2
max_s = int(infinity_order(bbase, R) - required)
assert max_s == 6

most_polar_s_a = int(abase_order - max_s)
max_t = int(infinity_order(tbase, R) - most_polar_s_a)
assert max_t == 5

labels = [("s", degree) for degree in range(max_s+1)]
labels += [("t", degree) for degree in range(max_t+1)]
assert len(labels) == 13

a_columns = []
b_columns = []
for kind, degree in labels:
    if kind == "s":
        a_columns.append(K(T**degree)*abase)
        b_columns.append(K(T**degree)*bbase)
    else:
        a_columns.append(K(T**degree)*tbase)
        b_columns.append(K.zero())

def principal_rows(values):
    den = common_denominator(values, R)
    nums = [R(value*den) if value else R.zero() for value in values]
    cutoff = int(den.degree()-required)
    top = max((degree_or_minus_one(value) for value in nums), default=-1)
    rows = [
        [
            value[degree] if value and degree <= value.degree() else QQ.zero()
            for value in nums
        ]
        for degree in range(cutoff+1, top+1)
    ]
    return rows, den, cutoff, top

a_rows, a_den, a_cutoff, a_top = principal_rows(a_columns)
b_rows, b_den, b_cutoff, b_top = principal_rows(b_columns)
condition = matrix(QQ, a_rows+b_rows, ncols=len(labels))
assert condition.nrows() == 11
assert condition.rank() == 11

kernel = condition.right_kernel_matrix()
assert kernel.nrows() == 2 and kernel.ncols() == 13

# Normalize the basis to primitive integer rows for reproducibility.
kernel_Z = matrix(ZZ, [primitive_integer_row(row) for row in kernel.rows()])
assert kernel_Z.rank() == 2
assert condition * kernel_Z.transpose() == matrix(QQ, condition.nrows(), 2, 0)

def st_from_row(row):
    s = R.zero()
    t = R.zero()
    for index, (kind, degree) in enumerate(labels):
        if kind == "s":
            s += QQ(row[index])*T**degree
        else:
            t += QQ(row[index])*T**degree
    return R(s), R(t)

def old_chord_pair(row):
    s, t = st_from_row(row)
    # f=s*(qreg+rho)+t*M = Acoef+Bcoef*m.
    Bcoef = K(s)/K(h)
    Acoef = (
        -K(s)*p_fun/K(h)
        - K(s)*normalizer/K(nx)
        + K(s*rho)
        + K(t*M)
    )
    # Check infinity in the smooth frame.
    aa = K(s)*abase + K(t)*tbase
    bb = K(s)*bbase
    assert infinity_order(aa, R) >= 2
    assert not bb or infinity_order(bb, R) >= 2
    return s, t, Acoef, Bcoef

pairs = [old_chord_pair(row) for row in kernel_Z.rows()]

print(
    "Q8QQRR|ambient=13|rows={}|rank={}|kernel={}|"
    "h_deg=10|Nx_deg={}|Ny_deg={}|R_deg={}|rho_deg={}|"
    "s_max=6|t_max=5".format(
        condition.nrows(), condition.rank(), kernel_Z.nrows(),
        nx.degree(), ny.degree(), normalizer.degree(), rho.degree(),
    ),
    flush=True,
)
for index, (row, pair) in enumerate(zip(kernel_Z.rows(), pairs)):
    s, t, unused_A, unused_B = pair
    print(
        "Q8QQKERNEL|index={}|row={}|s_deg={}|t_deg={}|s={}|t={}".format(
            index,
            ",".join(map(str, row)),
            degree_or_minus_one(s), degree_or_minus_one(t),
            s, t,
        ),
        flush=True,
    )

# ---------------------------------------------------------------------------
# 4. Eliminate the exact pencil to a binary quartic over QQ(U).
# ---------------------------------------------------------------------------

U_ring = PolynomialRing(QQ, "U")
U = U_ring.gen()
U_field = U_ring.fraction_field()
T_ring = PolynomialRing(U_field, "T")
TT = T_ring.gen()
T_field = T_ring.fraction_field()

def lift_poly(value):
    value = R(value)
    return T_ring([U_field(coefficient) for coefficient in value.list()])

def lift_rat(value):
    value = K(value)
    return T_field(lift_poly(value.numerator())) / T_field(
        lift_poly(value.denominator())
    )

A0, B0 = pairs[0][2], pairs[0][3]
A1, B1 = pairs[1][2], pairs[1][3]

# New base U=f1/f0 => f1-U*f0=0.
m_value = -(
    lift_rat(A1) - U_field(U)*lift_rat(A0)
) / (
    lift_rat(B1) - U_field(U)*lift_rat(B0)
)

sxU, syU = lift_rat(sx), lift_rat(sy)
AU, BU = lift_poly(A), lift_poly(Bcurve)

X_ring = PolynomialRing(T_field, "x")
x = X_ring.gen()
y_line = X_ring(m_value)*(x-X_ring(sxU))-X_ring(syU)
relation = y_line**2 - x**3 - X_ring(AU)*x - X_ring(BU)
quadratic, remainder = relation.quo_rem(x-X_ring(sxU))
assert not remainder and quadratic.degree() == 2

disc = T_field(
    quadratic[1]**2 - 4*quadratic[2]*quadratic[0]
)
assert disc

quartic, square_factor = squarefree_binary_quartic(disc, T_ring)
assert quartic.degree() == 4
I, J = binary_quartic_invariants(quartic)
jacobian_A = U_field(-27*I)
jacobian_B = U_field(-27*J)

print(
    "Q8QQQUARTIC|degree={}|quartic={}".format(quartic.degree(), quartic),
    flush=True,
)

# ---------------------------------------------------------------------------
# 5. Exact child minimization and root-lattice classification.
# ---------------------------------------------------------------------------

classification = classify_finite_short_weierstrass_fibres(
    U_ring, jacobian_A, jacobian_B
)

finite_data = [
    {
        "factor": str(item["factor"]),
        "degree": int(item["degree"]),
        "raw_orders": list(map(int, item["raw_orders"])),
        "scaling": int(item["scaling"]),
        "minimal_orders": list(map(int, item["minimal_orders"])),
        "kodaira": item["kodaira"],
    }
    for item in classification["finite_fibres"]
]

root_rank = int(classification["finite_root_rank"])
root_euler = int(classification["finite_euler_number"])
root_det = int(classification["finite_root_determinant"])

infinity = classification["infinity_boundary"]
infinity_orders = tuple(map(int, infinity["normalized_orders"]))
infinity_kind = "smooth"
if infinity_orders[2] > 0:
    irank, ieuler, idet, infinity_kind = kodaira_data_from_short_orders(
        *infinity_orders
    )
    root_rank += int(irank)
    root_euler += int(ieuler)
    root_det *= int(idet)

minimal_A = classification["finite_minimization"]["minimal_a"]
minimal_B = classification["finite_minimization"]["minimal_b"]
minimal_Delta = classification["finite_minimization"]["minimal_discriminant"]

assert minimal_A.degree() <= 8
assert minimal_B.degree() <= 12
assert minimal_Delta.degree() <= 24
assert root_euler == 24

all_kodaira = [entry["kodaira"] for entry in finite_data]
if infinity_kind != "smooth":
    all_kodaira.append(infinity_kind)

# D13 is exactly I9* as a single reducible fibre.
assert "I9*" in all_kodaira, all_kodaira
assert (root_rank, root_det) == (13, 4)
mw_rank_if_rho19 = 19 - 2 - root_rank
assert mw_rank_if_rho19 == 4

print(
    "Q8QQCHILD|finite={}|infinity={}|root_rank={}|root_euler={}|"
    "root_det={}|MW_rank={}|status=PASS_EXACT_CORRECTED_Q8_D13_CHILD".format(
        [(entry["degree"], entry["minimal_orders"], entry["kodaira"])
         for entry in finite_data],
        (infinity_orders, infinity_kind),
        root_rank, root_euler, root_det, mw_rank_if_rho19,
    ),
    flush=True,
)

def jsonable(value):
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    try:
        parent = value.parent()
        if parent is ZZ:
            return int(value)
        if parent is QQ:
            return str(value)
    except Exception:
        pass
    return value

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-corrected2cover-qq-child.v1",
    "status": "PASS_EXACT_CORRECTED_Q8_D13_CHILD",
    "corrections": {
        "binary_quartic_2cover": (
            "selected standard-Jacobian point is "
            "(phi(E7_7)-phi(old_O))+(phi(E7_7)-phi(affine_E7)); "
            "the previous point was its double"
        ),
        "q_normalizer": "R*h*Dy == Ny*Dx mod Nx",
    },
    "marking": {
        "height": "24",
        "O_intersection": 10,
        "collision_degree": 10,
        "section": point_data(S, R),
    },
    "rr": {
        "ambient_dimension": 13,
        "condition_rows": int(condition.nrows()),
        "condition_rank": int(condition.rank()),
        "kernel_dimension": int(kernel_Z.nrows()),
        "labels": [[kind, degree] for kind, degree in labels],
        "kernel_basis_integer_rows": [
            [int(value) for value in row] for row in kernel_Z.rows()
        ],
        "kernel_polynomials": [
            {
                "s": str(pair[0]),
                "t": str(pair[1]),
                "s_coefficients_low_to_high":
                    [str(value) for value in pair[0].list()],
                "t_coefficients_low_to_high":
                    [str(value) for value in pair[1].list()],
            }
            for pair in pairs
        ],
    },
    "pencil": {
        "new_base": "U=f1/f0",
        "branch_quartic_degree": int(quartic.degree()),
        "branch_quartic": str(quartic),
        "binary_quartic_I": str(I),
        "binary_quartic_J": str(J),
    },
    "child": {
        "weierstrass": "Y^2=X^3+A(U)X+B(U)",
        "minimal_A_coefficients_low_to_high":
            [str(value) for value in minimal_A.list()],
        "minimal_B_coefficients_low_to_high":
            [str(value) for value in minimal_B.list()],
        "minimal_Delta_coefficients_low_to_high":
            [str(value) for value in minimal_Delta.list()],
        "finite_fibres": finite_data,
        "infinity": {
            "raw_orders": list(map(int, infinity["raw_orders"])),
            "scaling": int(infinity["scaling"]),
            "minimal_orders": list(infinity_orders),
            "kodaira": infinity_kind,
        },
        "root_lattice": "D13",
        "root_rank": root_rank,
        "root_determinant": root_det,
        "root_euler": root_euler,
        "MW_rank_if_rho_19": mw_rank_if_rho19,
    },
    "boundary": (
        "This certifies the corrected q8 Riemann-Roch pencil and its exact "
        "D13/MW4 Jacobian child over QQ. Later neighbour/rootless stages and "
        "new MW section coordinates are separate."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(jsonable(payload), indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}")
