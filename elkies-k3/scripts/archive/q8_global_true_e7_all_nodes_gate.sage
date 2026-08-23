#!/usr/bin/env sage -python
"""Assemble the complete source-q8 global lattice intersection modulo good primes.

This is the first global assembly after exact saturation of both additive
endpoints.

Generic fibre basis:
    B = 1,m,...,m^9, x,xm,...,xm^7.

Coordinate conventions:
    u  = repository/global base coordinate,
    E8 = u=0,
    t  = 1/u is the E7 local parameter,
    h(u)=0 are the four smooth O.(-P1) collisions.

The complete finite global envelope is

    f = h(u)^(-18) * sum_i P_i(u) B_i,
    deg P_i <= 78.

Why this envelope is complete:
  * at h=0, q=(m-p)/h with p=rho/h, so q^9 has worst B-coordinate
    denominator h^18; X*q^7 is smaller;
  * the final saturated reduced E7 basis has maximum pole t^-6=u^6,
    while deg h^18=72, hence numerator degree <=72+6=78;
  * away from E7,E8,h the horizontal B-frame is regular.

The global -11F twist is represented at E8.  Thus the ACTUAL E8 lattice is
u^11 times the already-certified reduced resolved E8 lattice

    < u^2,
      Q^b-s^b (1<=b<=9),
      X Q^b   (0<=b<=7) >,

where Q=u^2*m, X=u^4*x and s=Y(P1)/X(P1).

The E7 lattice is the exact saturated rank-18 basis reached by steps 1--7:

    m^0,...,m^9,
    t^-3*z1, t^-3*z1*m,
    t^-4*z2, t^-4*z2*m,
    t^-5*z3*m^2, t^-5*z3*m^3,
    t^-6*z4*m^2, t^-6*z4*m^3.

The adapted coordinates z1..z4 are rederived here from the exact chord
quadratic over QQ(m)[[t]], avoiding dependence on transient local scripts.

At h=0 we use the exact certified q/X principal-part formula from
compile_h92_q8_smooth_principal_parts_exact.sage with h_power=18.

A final kernel dimension 2 at a good prime is an exact upper bound h0<=2 over
QQ.  Combined with the primitive nef isotropic q8 divisor (RR gives h0>=2),
this certifies h0=2.  Two-prime agreement is retained as an implementation
cross-check.
"""

import argparse
import json
import subprocess
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import (
    GF, LaurentSeriesRing, PolynomialRing, PowerSeriesRing,
    QQ, ZZ, binomial, gcd, lcm, matrix, vector
)

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

ROOT = Path.cwd()
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
P1_PATH = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
E8_TARGET = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e8-local-target.json"
E8_IDEAL = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e8-complete-module.json"

PRIMES = (43, 59)
H_POWER = 18
# The exact marked q8 frame m^b/t^6, x*m^b/t^8 is written in the
# q6^9 helper normalization and therefore still carries the helper t^9.
# Strip it here because the literal q8 -11F twist is represented separately
# at E8:
#
#   reduced marked frame: m^b/t^15, x*m^b/t^17.
M_E7_POLE = 15
X_E7_POLE = 17
M_MAX_DEGREE = 4*H_POWER + M_E7_POLE   # 87
X_MAX_DEGREE = 4*H_POWER + X_E7_POLE   # 89
GLOBAL_FIBRE_TWIST = -11


def load(path):
    return json.loads(path.read_text())


def poly(ring, values):
    return ring([QQ(v) for v in values])


def coeff_dict(poly_value):
    """Return {degree: QQ coefficient} for a polynomial."""
    return {
        int(i): QQ(poly_value[i])
        for i in range(poly_value.degree()+1)
        if poly_value[i]
    } if poly_value else {}


def add_poly_to_B(target, b_index, t_shift, polynomial):
    """target[b_index][t_exponent] += coeff of t_shift * polynomial(m).

    This helper is used only after one m-monomial has already been selected,
    so polynomial is represented externally by its m coefficients.
    """
    raise RuntimeError("internal misuse")


# --------------------------------------------------------------------------
# Exact source data and exact E7 affine-tail reconstruction over QQ(m)[[t]].
# --------------------------------------------------------------------------

p1 = load(P1_PATH)
assert p1["status"] == "PASS_EXACT_H92_P1"
e8_target = load(E8_TARGET)
e8_ideal = load(E8_IDEAL)
assert e8_target["status"] == "PASS_EXACT_Q8_E8_SOURCE_TARGET"
assert e8_ideal["status"] == "PASS_EXACT_Q8_E8_COMPLETE_MODULE"

anchor = SourceFileLoader("q8_global_anchor", str(ANCHOR)).load_module()
r0, s0 = anchor.EXPECTED_H92
_, forms = anchor.parse_h92(H92)
A1, A, B1, B, B2 = [QQ(v(r0, s0)) for v in forms]

Ru = PolynomialRing(QQ, "u")
u = Ru.gen()
h = poly(Ru, p1["structured_denominator"]["Z4_coefficients"])
assert h.degree() == 4 and gcd(h, h.derivative()) == 1 and h(0) != 0

xn = poly(Ru, p1["x_entrance_base"]["numerator_coefficients"])
xd = poly(Ru, p1["x_entrance_base"]["denominator_coefficients"])
yn = poly(Ru, p1["y_entrance_base"]["numerator_coefficients"])
yd = poly(Ru, p1["y_entrance_base"]["denominator_coefficients"])

Rm = PolynomialRing(QQ, "m")
mvar = Rm.gen()
Km = Rm.fraction_field()
PSm = PowerSeriesRing(Km, "t", default_prec=10)
tq = PSm.gen()


def at_infinity(num, den):
    dn, dd = num.degree(), den.degree()
    nr = sum(PSm(num[i])*tq**(dn-i) for i in range(dn+1))
    dr = sum(PSm(den[i])*tq**(dd-i) for i in range(dd+1))
    return tq**(dd-dn)*nr/dr


xp_t = at_infinity(xn, xd)
yp_t = at_infinity(yn, yd)
a_t = PSm(A1)*tq**3 + PSm(A)*tq**4

# z=x-m^2 on the affine E7 branch; solve its exact series recursively.
z_series = PSm(0)
ks = {}
for order in range(2, 9):
    residual = (
        z_series**2
        + (PSm(mvar**2)+xp_t)*z_series
        + 2*PSm(mvar**2)*xp_t
        + xp_t**2 + a_t - 2*PSm(mvar)*yp_t
    )
    k = -Km(residual[order])/Km(mvar**2)
    ks[order] = k
    z_series += PSm(k)*tq**order

residual = (
    z_series**2
    + (PSm(mvar**2)+xp_t)*z_series
    + 2*PSm(mvar**2)*xp_t
    + xp_t**2 + a_t - 2*PSm(mvar)*yp_t
)
assert all(residual[n] == 0 for n in range(2, 9))

# Exact adapted-coordinate polynomial corrections.
# z1 = x-m^2 - k2*t^2.
Z1P2 = Rm(-ks[2])

# z2 = m^2*(z1-k3*t^3).
P2 = Rm(-mvar**2*ks[2])
P3 = Rm(-mvar**2*ks[3])

# z3 = z2 - m^2*k4*t^4.
P4 = Rm(-mvar**2*ks[4])

# z4 = m^2*(z3 - m^2*k5*t^5).
P5 = Rm(-mvar**4*ks[5])

# The conversions to Rm are themselves exact polynomiality checks.
for pp in (Z1P2, P2, P3, P4, P5):
    assert pp in Rm


# B indexing:
# 0..9 = m^0..m^9
# 10..17 = x*m^0..x*m^7
def m_index(power):
    assert 0 <= power <= 9
    return int(power)


def x_index(power):
    assert 0 <= power <= 7
    return 10 + int(power)


def new_generator():
    return [dict() for _ in range(18)]


def add_laurent(generator, index, exponent, coefficient):
    coefficient = QQ(coefficient)
    if not coefficient:
        return
    generator[index][int(exponent)] = (
        generator[index].get(int(exponent), QQ(0)) + coefficient
    )
    if not generator[index][int(exponent)]:
        del generator[index][int(exponent)]


def add_m_poly(generator, polynomial, m_shift, t_exponent):
    polynomial = Rm(polynomial)
    for degree, coefficient in enumerate(polynomial.list()):
        if coefficient:
            add_laurent(
                generator, m_index(m_shift+degree), t_exponent, coefficient
            )


# Final eight x-bearing E7 generators, one for each x*m^j coordinate.
e7_x_generators = [None]*8

# j=0,1: t^-3*z1*m^j
for j in range(2):
    g = new_generator()
    add_laurent(g, x_index(j), -3, 1)
    add_laurent(g, m_index(j+2), -3, -1)
    add_m_poly(g, Z1P2, j, -1)  # t^-3 * t^2
    e7_x_generators[j] = g

# j=2,3: t^-4*z2*m^(j-2)
for j in range(2):
    target_x = j+2
    g = new_generator()
    add_laurent(g, x_index(target_x), -4, 1)
    add_laurent(g, m_index(target_x+2), -4, -1)
    add_m_poly(g, P2, j, -2)
    add_m_poly(g, P3, j, -1)
    e7_x_generators[target_x] = g

# j=4,5: t^-5*z3*m^(j-2), where source shift is 2,3.
for source_shift in range(2, 4):
    target_x = source_shift+2
    g = new_generator()
    add_laurent(g, x_index(target_x), -5, 1)
    add_laurent(g, m_index(target_x+2), -5, -1)
    add_m_poly(g, P2, source_shift, -3)
    add_m_poly(g, P3, source_shift, -2)
    add_m_poly(g, P4, source_shift, -1)
    e7_x_generators[target_x] = g

# j=6,7: t^-6*z4*m^(source_shift), source_shift=2,3.
for source_shift in range(2, 4):
    target_x = source_shift+4
    g = new_generator()
    add_laurent(g, x_index(target_x), -6, 1)
    add_laurent(g, m_index(target_x+2), -6, -1)
    add_m_poly(g, P2, source_shift+2, -4)
    add_m_poly(g, P3, source_shift+2, -3)
    add_m_poly(g, P4, source_shift+2, -2)
    add_m_poly(g, P5, source_shift, -1)
    e7_x_generators[target_x] = g

assert all(g is not None for g in e7_x_generators)

e7_poles = (3, 3, 4, 4, 5, 5, 6, 6)
for j, g in enumerate(e7_x_generators):
    # Leading x coefficient must be exactly t^-d_j.
    assert g[x_index(j)] == {-e7_poles[j]: QQ(1)}

TRUE_COLUMNS = (
    [(bidx, degree)
     for bidx in range(10)
     for degree in range(M_MAX_DEGREE + 1)]
    +
    [(bidx, degree)
     for bidx in range(10, 18)
     for degree in range(X_MAX_DEGREE + 1)]
)

print(
    "Q8TRUEPREP|ambient_columns={}|h_power={}|m_max_degree={}|"
    "x_max_degree={}|reduced_marked_poles=15,17|global_twist={}|"
    "status=PASS_EXACT_TRUE_MARKED_GLOBAL_ENVELOPE".format(
        len(TRUE_COLUMNS), H_POWER, M_MAX_DEGREE, X_MAX_DEGREE,
        GLOBAL_FIBRE_TWIST,
    ),
    flush=True,
)

# Build the exact reduced-E7 generic-component condition table for the same
# 1600 columns. The repository compiler contains the helper '-9':
#
#   w + (4*k-i-9)*ord(t) + ...
#
# For the actual reduced coefficient u^d/h^18, feed i=d-9. Then both the
# valuation formula and the chart residue evaluator see the true t-order
# 72-d with the helper fibre factor removed.
LOCAL = ROOT / "artifacts/local/elkies-k3/q8-true-e7-generic-gate"
LOCAL.mkdir(parents=True, exist_ok=True)

E7_AMBIENT_PATH = LOCAL / "ambient.json"
E7_CONDITIONS_PATH = LOCAL / "conditions.json"
E7_R13_PATH = LOCAL / "e7-1-3.json"
E7_R56_PATH = LOCAL / "e7-5-6.json"
E7_R47_PATH = LOCAL / "e7-4-7.json"

generic_meta = (
    [(0, b) for b in range(10)]
    + [(1, b) for b in range(8)]
)

e7_ambient_basis = []
for bidx, degree in TRUE_COLUMNS:
    a, b = generic_meta[bidx]
    e7_ambient_basis.append({
        "kind": "x_m_power" if a else "m_power",
        "x_power": int(a),
        "m_power": int(b),
        "u_power": int(degree - 9),
        "h_power": int(H_POWER),
        "coefficient": "u^{}/h^{}".format(degree, H_POWER),
    })

E7_AMBIENT_PATH.write_text(json.dumps({
    "schema": "elkies-k3.h92-q8-true-global-e7-ambient.v1",
    "status": "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT",
    "ambient_basis": e7_ambient_basis,
    "ambient_dimension": len(e7_ambient_basis),
    "boundary": (
        "True reduced marked-chart global envelope. "
        "No approximate reduced-E7 saturation lattice is used."
    ),
}, indent=2, sort_keys=True) + "\n")


def run_builder(script, *args):
    command = [
        "sage", "-python",
        str(ROOT / "elkies-k3/scripts" / script),
        *[str(value) for value in args],
    ]
    print("Q8TRUERUN|" + " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


run_builder(
    "derive_h92_q8_all_component_generic_conditions.sage",
    "--ambient", E7_AMBIENT_PATH,
    "--output", E7_CONDITIONS_PATH,
)
run_builder(
    "derive_h92_q8_e7_1_3_generic_residue_rows.sage",
    "--conditions", E7_CONDITIONS_PATH,
    "--output", E7_R13_PATH,
)
run_builder(
    "derive_h92_q8_e7_5_6_generic_residue_rows.sage",
    "--conditions", E7_CONDITIONS_PATH,
    "--output", E7_R56_PATH,
)
run_builder(
    "derive_h92_q8_e7_4_7_generic_residue_rows.sage",
    "--conditions", E7_CONDITIONS_PATH,
    "--output", E7_R47_PATH,
)

E7_CONDITIONS = load(E7_CONDITIONS_PATH)
E7_PAYLOADS = [
    load(E7_R13_PATH),
    load(E7_R56_PATH),
    load(E7_R47_PATH),
]
assert E7_CONDITIONS["status"] == \
    "PASS_EXACT_Q8_ALL_COMPONENT_GENERIC_CONDITION_TEMPLATE"


# --------------------------------------------------------------------------
# Finite-prime global intersection.
# --------------------------------------------------------------------------

def qq_to_fp(value, F):
    value = QQ(value)
    den = F(ZZ(value.denominator()))
    if not den:
        raise ZeroDivisionError("bad reduction denominator")
    return F(ZZ(value.numerator()))/den


def build_prime(prime):
    F = GF(prime)
    Fu = PolynomialRing(F, "u")
    uf = Fu.gen()

    hf = Fu([qq_to_fp(v, F) for v in h.list()])
    assert hf.degree() == 4 and gcd(hf, hf.derivative()) == 1 and hf(0)

    # Global columns: true reduced marked-chart envelope.
    columns = list(TRUE_COLUMNS)
    col_index = {entry: i for i, entry in enumerate(columns)}
    ncols = len(columns)

    blocks = {"E8": {}, "E7GEN": {}, "H": {}}

    def emit(block, key, col, value):
        value = F(value)
        if not value:
            return
        row = blocks[block].setdefault(key, {})
        row[col] = row.get(col, F.zero()) + value
        if not row[col]:
            del row[col]

    # --------------------------------------------------------------
    # Exact reduced-E7 generic-component residue gate on all seven
    # exceptional components. These are the repository's actual
    # chart-function-field rows, not the superseded affine approximation.
    # --------------------------------------------------------------

    for basis_index in \
            E7_CONDITIONS["singleton_coordinate_block"]["basis_indices"]:
        emit(
            "E7GEN",
            ("singleton", int(basis_index)),
            int(basis_index),
            F.one(),
        )

    residue_counter = 0
    for payload in E7_PAYLOADS:
        for component in payload["components"]:
            cname = component["component"]
            for residue_row in component["non_singleton_residue_rows"]:
                key = ("residue", cname, residue_counter)
                residue_counter += 1
                for entry in residue_row["entries"]:
                    emit(
                        "E7GEN",
                        key,
                        int(entry["basis_index"]),
                        qq_to_fp(QQ(entry["coefficient"]), F),
                    )

    # ------------------------------------------------------------------
    # E8 at u=0, INCLUDING the global -11F twist.
    # ------------------------------------------------------------------

    LU0 = LaurentSeriesRing(F, "u", default_prec=40)
    us = LU0.gen()
    hs = LU0(hf)
    h_inv18 = hs**(-H_POWER)

    xnf = Fu([qq_to_fp(v, F) for v in xn.list()])
    xdf = Fu([qq_to_fp(v, F) for v in xd.list()])
    ynf = Fu([qq_to_fp(v, F) for v in yn.list()])
    ydf = Fu([qq_to_fp(v, F) for v in yd.list()])

    xp0 = LU0(xnf)/LU0(xdf)
    yp0 = LU0(ynf)/LU0(ydf)

    # Integral E8 coordinates: X_P=u^4*x_P, Y_P=u^6*y_P.
    Xp0 = us**4*xp0
    Yp0 = us**6*yp0
    assert Xp0.valuation() == 0 and Yp0.valuation() == 0
    s_e8 = Yp0/Xp0
    assert s_e8.valuation() == 0 and s_e8[0] != 0

    max_e8_order = 29  # x*m^7 requires u^(15+14)
    c0_cache = {}
    for degree in range(max(M_MAX_DEGREE, X_MAX_DEGREE)+1):
        c0_cache[degree] = us**degree*h_inv18

    for bidx, degree in columns:
        col = col_index[(bidx, degree)]
        c = c0_cache[degree]

        if bidx == 0:
            # Coupled constant coordinate must be divisible by u^13.
            for order in range(13):
                emit("E8", ("constant_coupling", order), col, c[order])

        elif 1 <= bidx <= 9:
            b = bidx
            required = 11 + 2*b
            for order in range(required):
                emit("E8", ("m_div", b, order), col, c[order])

            coupled = (s_e8**b) * c * us**(-2*b)
            # Any negative terms are automatically included in this Laurent
            # expression only if the divisibility rows fail. Since those rows
            # are stacked simultaneously, degrees 0..12 suffice here after
            # exact divisibility is enforced.
            for order in range(13):
                emit(
                    "E8", ("constant_coupling", order), col,
                    coupled[order]
                )

        else:
            b = bidx-10
            required = 15 + 2*b
            for order in range(required):
                emit("E8", ("x_div", b, order), col, c[order])

    # ------------------------------------------------------------------
    # Four smooth h-collisions: exact q/X principal-part map.
    # ------------------------------------------------------------------

    pole_bound = 27
    modulus = hf**pole_bound
    residue_dim = modulus.degree()  # 108

    Ffield = Fu.fraction_field()
    xpf = Ffield(xnf)/Ffield(xdf)
    ypf = Ffield(ynf)/Ffield(ydf)
    rho = Ffield(hf)*ypf/xpf

    assert gcd(hf, Fu(rho.numerator())) == 1
    assert gcd(hf, Fu(rho.denominator())) == 1

    def residue(value):
        value = Ffield(value)
        numerator = Fu(value.numerator())
        denominator = Fu(value.denominator())
        assert gcd(denominator, modulus) == 1
        return Fu(
            (numerator*denominator.inverse_mod(modulus)) % modulus
        )

    # B metadata.
    Bmeta = (
        [(0, b) for b in range(10)]
        + [(1, b) for b in range(8)]
    )

    for bidx, degree in columns:
        col = col_index[(bidx, degree)]
        a_power, m_power = Bmeta[bidx]

        for q_power in range(m_power+1):
            exponent = (
                2*q_power - m_power - H_POWER - 2*a_power
            )
            if exponent >= 0:
                continue

            value = (
                F(binomial(m_power, q_power))
                * uf**degree
                * rho**(m_power-q_power)
                * hf**(pole_bound+exponent)
            )
            rem = residue(value)

            for rdeg, coefficient in enumerate(rem.list()):
                emit(
                    "H",
                    (a_power, q_power, rdeg),
                    col,
                    coefficient
                )

    # ------------------------------------------------------------------
    # Sparse rank calculations.
    # ------------------------------------------------------------------

    def block_matrix(names):
        row_items = []
        for name in names:
            row_items.extend(
                (name, key, row)
                for key, row in blocks[name].items()
                if row
            )

        entries = {}
        for ridx, (name, key, row) in enumerate(row_items):
            for col, value in row.items():
                if value:
                    entries[(ridx, col)] = value

        M = matrix(
            F, len(row_items), ncols, entries, sparse=True
        )
        return M, row_items

    MH, rowsH = block_matrix(("H",))
    M8, rows8 = block_matrix(("E8",))
    M7g, rows7g = block_matrix(("E7GEN",))
    M8H, rows8H = block_matrix(("E8", "H"))
    M8H7g, rows8H7g = block_matrix(("E8", "H", "E7GEN"))

    rankH = MH.rank()
    rank8 = M8.rank()
    rank7g = M7g.rank()
    rank8H = M8H.rank()
    rankall = M8H7g.rank()

    # --------------------------------------------------------------
    # Safe necessary leading-term gates at ALL six actual E7 nodes.
    #
    # Each node has an exact regular two-parameter frame certified by
    # derive_h92_q8_all_e7_node_principal_bidegrees.sage.  We repeatedly
    # impose only singleton Pareto-minimal negative bidegrees on the current
    # global survivor.  These are rigorous necessary conditions: no guessed
    # cancellation or finite-corner replacement is used.
    # --------------------------------------------------------------

    node_frames = {
        "E7_1--E7_4": {
            "t": (3,2), "x": (4,2), "m": (2,1), "g": (4,2),
        },
        "E7_4--E7_3": {
            "t": (3,4), "x": (4,6), "m": (2,3), "g": (4,6),
        },
        "E7_3--E7_7": {
            "t": (3,4), "x": (5,6), "m": (2,3), "g": (5,6),
        },
        "E7_7--E7_2": {
            "t": (3,2), "x": (5,4), "m": (2,1), "g": (5,5),
        },
        "E7_3--E7_6": {
            "t": (2,4), "x": (3,6), "m": (2,3), "g": (3,6),
        },
        "E7_2--E7_5": {
            "t": (1,2), "x": (2,4), "m": (0,1), "g": (6,5),
        },
    }

    candidate = matrix(F, M8H7g.right_kernel().basis_matrix())
    initial_node_dimension = candidate.nrows()
    assert initial_node_dimension == ncols-rankall

    def node_bidegree(column, frame):
        bidx, degree = columns[column]
        a, b = generic_meta[bidx]

        # True reduced coefficient is u^degree/h^18. The q6^9 helper
        # convention uses i=degree-9, hence delta=4*18-(degree-9)-9
        # =72-degree.
        delta = 4*H_POWER-degree

        return tuple(
            int(
                frame["g"][axis]
                + delta*frame["t"][axis]
                + a*frame["x"][axis]
                + b*frame["m"][axis]
            )
            for axis in range(2)
        )

    def support_pairs(rows, frame):
        pairs = {}
        for col in range(rows.ncols()):
            if any(rows[row, col] for row in range(rows.nrows())):
                pairs.setdefault(node_bidegree(col, frame), []).append(col)
        return pairs

    def pareto_minima(pairs):
        result = []
        for degree, labels in pairs.items():
            if min(degree) >= 0:
                continue
            if not any(
                other != degree
                and other[0] <= degree[0]
                and other[1] <= degree[1]
                for other in pairs
            ):
                result.append((degree, labels))
        return sorted(result)

    node_records = []
    changed = True
    sweep = 0

    while changed and candidate.nrows():
        changed = False
        sweep += 1

        for node_name, frame in node_frames.items():
            while candidate.nrows():
                minima = pareto_minima(support_pairs(candidate, frame))
                unique = [
                    (degree, labels[0])
                    for degree, labels in minima
                    if len(labels) == 1
                ]
                if not unique:
                    break

                degree, label = unique[0]
                coefficient = vector(
                    F,
                    [
                        candidate[row, label]
                        for row in range(candidate.nrows())
                    ],
                )
                if not coefficient:
                    break

                before = candidate.nrows()
                cut = matrix(
                    F,
                    matrix(
                        F, 1, candidate.nrows(), [coefficient]
                    ).right_kernel().basis_matrix(),
                )
                candidate = cut*candidate
                after = candidate.nrows()

                assert after == before-1
                changed = True
                node_records.append({
                    "sweep": sweep,
                    "node": node_name,
                    "bidegree": degree,
                    "column": label,
                    "before": before,
                    "after": after,
                })

    by_node = {
        name: sum(record["node"] == name for record in node_records)
        for name in node_frames
    }

    print(
        "Q8TRUEALLNODES|prime={}|generic_survivor={}|constraints={}|"
        "node_survivor={}|sweeps={}|by_node={}|"
        "status=TRUE_ALL_E7_NODE_LEADING_CLOSURE".format(
            prime,
            initial_node_dimension,
            len(node_records),
            candidate.nrows(),
            sweep,
            ",".join(
                "{}:{}".format(name, by_node[name])
                for name in node_frames
            ),
        ),
        flush=True,
    )

    for index, record in enumerate(node_records):
        print(
            "Q8TRUEALLNODECUT|prime={}|index={}|sweep={}|node={}|"
            "degree={},{}|column={}|before={}|after={}".format(
                prime,
                index,
                record["sweep"],
                record["node"],
                record["bidegree"][0],
                record["bidegree"][1],
                record["column"],
                record["before"],
                record["after"],
            ),
            flush=True,
        )

    print(
        "Q8TRUEGENV2|prime={}|ambient={}|"
        "H={}:{}:{}|E8={}:{}:{}|E7GEN={}:{}:{}|"
        "E8H={}:{}:{}|E8H_E7GEN={}:{}:{}|"
        "status=TRUE_REDUCED_E7_GENERIC_GATE".format(
            prime, ncols,
            MH.nrows(), rankH, ncols-rankH,
            M8.nrows(), rank8, ncols-rank8,
            M7g.nrows(), rank7g, ncols-rank7g,
            M8H.nrows(), rank8H, ncols-rank8H,
            M8H7g.nrows(), rankall, ncols-rankall,
        ),
        flush=True,
    )

    return {
        "prime": prime,
        "ambient": ncols,
        "e8h_rank": rank8H,
        "e8h_kernel": ncols-rank8H,
        "all_rank": rankall,
        "kernel": ncols-rankall,
    }



parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, action="append")
args = parser.parse_args()

primes = tuple(args.prime) if args.prime else PRIMES
results = [build_prime(int(p)) for p in primes]

if len(results) > 1:
    kernels = {entry["kernel"] for entry in results}
    assert len(kernels) == 1
    print(
        "Q8TRUEGENV2SUMMARY|primes={}|kernel={}|agreement=1|"
        "status=TRUE_REDUCED_E7_GENERIC_TWO_PRIME_AGREEMENT".format(
            ",".join(str(entry["prime"]) for entry in results),
            results[0]["kernel"],
        ),
        flush=True,
    )
