from sage.all import *
import argparse
import json
from pathlib import Path


# Exact local deformation audit for the first canonical A10/MW3 target triple.
# The point was found by seed 17 of the normalized P1 search over GF(23).
ap = argparse.ArgumentParser()
ap.add_argument("--hensel-digits", type=int, default=0)
ap.add_argument("--target", type=int, choices=(1, 2), default=1)
ap.add_argument("--hensel-json", default=None)
ap.add_argument("--series-order", type=int, default=0)
ap.add_argument("--series-json", default=None)
ap.add_argument("--relation-x-degree", type=int, default=8)
ap.add_argument("--relation-y-degree", type=int, default=4)
args = ap.parse_args()

p = 23
Kp = GF(p)
base_names = [
    "rho", "r1", "s1", "lam", "x2", "x3", "y3", "y4", "y5",
    "a4", "a5", "a6",
]
p2_names = ["r2", "p20", "p21", "p22", "c20", "c21"]
p3_names = ["r3", "p30", "p31", "p32", "p33"] + [f"v3{i}" for i in range(9)]
# Build over QQ first.  Reducing these equations modulo 23 audits the target;
# clearing their numeric contents gives the genuine integral system for lift.
S = PolynomialRing(QQ, base_names + p2_names + p3_names, order="degrevlex")
d = S.gens_dict()
SF = FractionField(S)
St = PolynomialRing(SF, "t")
t = St.gen()


def V(name):
    return SF(d[name])


rho, r1, s1, lam = [V(name) for name in ("rho", "r1", "s1", "lam")]
x2, x3, y3, y4, y5 = [V(name) for name in ("x2", "x3", "y3", "y4", "y5")]
a4, a5, a6 = [V(name) for name in ("a4", "a5", "a6")]
s0 = SF(3)
sinf = 3 * rho**2

# This is the unreduced normalized P1 chart.  Keeping its twelve coordinates
# makes the extended-section Jacobian considerably smaller than substituting
# the large four-equation elimination formulas.
X_at_1 = r1**2 - 2 * s1
Y_at_1 = r1 * (r1**2 - 3 * s1)
x1 = X_at_1 - s0 - x2 - x3 - sinf
y2 = -(
    lam * Y_at_1
    + (lam**3 - lam) * y3
    + (lam**4 - lam) * y4
    + (lam**5 - lam) * y5
) / (lam**2 - lam)
y1 = Y_at_1 - y2 - y3 - y4 - y5
X1 = s0 + x1 * t + x2 * t**2 + x3 * t**3 + sinf * t**4
Y1 = y1 * t + y2 * t**2 + y3 * t**3 + y4 * t**4 + y5 * t**5
sl = X1(lam)

a0 = -3 * s0**2
a1 = 6 * y1 - 6 * s0 * x1
a7 = 6 * rho * y5 - 6 * sinf * x3
a8 = -3 * sinf**2
a2_plus_a3 = -3 * s1**2 - (a0 + a1 + a4 + a5 + a6 + a7 + a8)
a3 = -(
    a0 + a1 * lam + a2_plus_a3 * lam**2 + a4 * lam**4
    + a5 * lam**5 + a6 * lam**6 + a7 * lam**7 + a8 * lam**8
    + 3 * sl**2
) / (lam**3 - lam**2)
a2 = a2_plus_a3 - a3
A = sum(value * t**i for i, value in enumerate((a0, a1, a2, a3, a4, a5, a6, a7, a8)))
B = Y1**2 - X1**3 - A * X1
Delta = -16 * (4 * A**3 + 27 * B**2)


def numerator_equations(polynomial):
    return [S(SF(coefficient).numerator()) for coefficient in polynomial.list() if coefficient != 0]


equations = [S(d["y5"])]
equations += numerator_equations(St([B.derivative(t)(1) + s1 * A.derivative(t)(1)]))
equations += numerator_equations(sum(Delta[k] * t**k for k in range(14, 22)))

# Canonical P2: profile (6,2,1,1), P2.O=1.  Its numerator Y vanishes at all
# three finite reducible fibers and has degree at most four.
r2, p20, p21, p22, c20, c21 = [V(name) for name in p2_names]
finite_points = [SF(0), SF(1), lam]
finite_nodes = [s0, s1, sl]
F3 = t * (t - 1) * (t - lam)
C2 = St(0)
for i, (point, node) in enumerate(zip(finite_points, finite_nodes)):
    basis = St(1)
    denominator = SF(1)
    for j, other in enumerate(finite_points):
        if i == j:
            continue
        basis *= t - other
        denominator *= point - other
    C2 += node * (point - r2)**2 * basis / denominator
X2 = C2 + F3 * (p20 + p21 * t + p22 * t**2 + sinf * t**3)
Y2 = F3 * (c20 + c21 * t)
z2 = t - r2
H2 = X2**3 + A * X2 * z2**4 + B * z2**6 - Y2**2
equations += numerator_equations(H2)

# Canonical P3: profile (10,2,0,1), P3.O=1.  The finite X incidences are at
# 0 and lambda; Y is kept coefficient-wise so the full square identity is a
# polynomial deformation system rather than a post-hoc square-root test.
r3, p30, p31, p32, p33 = [V(name) for name in p3_names[:5]]
v3 = [V(f"v3{i}") for i in range(9)]
C3 = s0 * r3**2 * (t - lam) / (-lam)
C3 += sl * (lam - r3)**2 * t / lam
G3 = t * (t - lam)
X3 = C3 + G3 * (p30 + p31 * t + p32 * t**2 + p33 * t**3 + sinf * t**4)
Y3 = sum(v3[i] * t**i for i in range(9))
z3 = t - r3
H3 = X3**3 + A * X3 * z3**4 + B * z3**6 - Y3**2
equations += numerator_equations(H3)
equations += [S(SF(Y3(0)).numerator()), S(SF(Y3(lam)).numerator())]
raw_equations_all = (
    [y5, B.derivative(t)(1) + s1 * A.derivative(t)(1)]
    + [Delta[k] for k in range(21, 13, -1)]
    + [H2[k] for k in range(19)]
    + [H3[k] for k in range(19)]
    + [Y3(0), Y3(lam)]
)
raw_nonzero_indices = [i for i, equation in enumerate(raw_equations_all) if equation != 0]
raw_equations = [raw_equations_all[i] for i in raw_nonzero_indices]

# Remove exact duplicates and zero identities without changing the local ideal.
unique_equations = []
seen = set()
for equation in equations:
    if equation == 0:
        continue
    normalized = equation / equation.lc()
    key = str(normalized)
    if key in seen:
        continue
    seen.add(key)
    unique_equations.append(equation)
equations = unique_equations

target_data = {
    1: {
        "base": {
            "rho": 16, "r1": 13, "s1": 13, "lam": 9,
            "x2": 7, "x3": 17, "y3": 13, "y4": 5, "y5": 0,
            "a4": 21, "a5": 1, "a6": 17,
            "r2": 11, "p20": 16, "p21": 16, "p22": 1,
            "r3": 8, "p30": 8, "p31": 22, "p32": 6, "p33": 21,
        },
        "Y2": (0, 7, 22, 10, 7),
        "Y3": (0, 5, 18, 8, 18, 11, 5, 21, 4),
    },
    2: {
        "base": {
            "rho": 16, "r1": 18, "s1": 10, "lam": 8,
            "x2": 4, "x3": 9, "y3": 1, "y4": 18, "y5": 0,
            "a4": 3, "a5": 7, "a6": 4,
            "r2": 21, "p20": 22, "p21": 2, "p22": 11,
            "r3": 9, "p30": 11, "p31": 11, "p32": 17, "p33": 12,
        },
        "Y2": (0, 11, 21, 3, 11),
        # The raw P3 has profile (1,1); negate it to orient it as (10,2).
        "Y3": (0, 16, 4, 10, 11, 12, 2, 5, 21),
    },
}
selected = target_data[args.target]
target = dict(selected["base"])

# Recover the two P2 Y parameters from the certified numerator, avoiding a
# second hand-maintained representation in the certificate.
Kt = PolynomialRing(Kp, "u")
u = Kt.gen()
Y2_target = Kt(selected["Y2"])
quotient, remainder = Y2_target.quo_rem(
    u * (u - 1) * (u - Kp(target["lam"]))
)
if remainder or quotient.degree() > 1:
    raise RuntimeError("certified P2 numerator lost its finite-fiber factors")
target["c20"] = int(quotient[0])
target["c21"] = int(quotient[1])
for i, value in enumerate(selected["Y3"]):
    target[f"v3{i}"] = value

SZ = PolynomialRing(ZZ, S.variable_names(), order="degrevlex")
integral_equations = []
for equation in equations:
    denominator = lcm([coefficient.denominator() for coefficient in equation.coefficients()])
    integral = SZ(denominator * equation)
    content = gcd([abs(coefficient) for coefficient in integral.coefficients()])
    if content > 1:
        integral //= content
    integral_equations.append(integral)

Sp = PolynomialRing(Kp, S.variable_names(), order="degrevlex")
dp = Sp.gens_dict()
mod_equations = [Sp(equation) for equation in integral_equations]
point = {dp[name]: Kp(value) for name, value in target.items()}
values = [Kp(equation.subs(point)) for equation in mod_equations]
if any(values):
    failures = [i for i, value in enumerate(values) if value]
    raise RuntimeError(f"target fails combined equations {failures}")


def structured_residuals(coordinates, coefficient_ring):
    """Evaluate the unexpanded section system over a coefficient ring."""
    values = dict(zip(S.variable_names(), coordinates))
    CT = PolynomialRing(coefficient_ring, "tt")
    tt = CT.gen()
    rrho = values["rho"]
    rr1 = values["r1"]
    ss1 = values["s1"]
    llam = values["lam"]
    xx2, xx3 = values["x2"], values["x3"]
    yy3, yy4, yy5 = values["y3"], values["y4"], values["y5"]
    aa4, aa5, aa6 = values["a4"], values["a5"], values["a6"]
    ss0 = coefficient_ring(3)
    ssinf = 3 * rrho**2
    xx_at_1 = rr1**2 - 2 * ss1
    yy_at_1 = rr1 * (rr1**2 - 3 * ss1)
    xx1 = xx_at_1 - ss0 - xx2 - xx3 - ssinf
    yy2 = -(
        llam * yy_at_1 + (llam**3 - llam) * yy3
        + (llam**4 - llam) * yy4 + (llam**5 - llam) * yy5
    ) / (llam**2 - llam)
    yy1 = yy_at_1 - yy2 - yy3 - yy4 - yy5
    XX1 = ss0 + xx1 * tt + xx2 * tt**2 + xx3 * tt**3 + ssinf * tt**4
    YY1 = yy1 * tt + yy2 * tt**2 + yy3 * tt**3 + yy4 * tt**4 + yy5 * tt**5
    ssl = XX1(llam)
    aa0 = -3 * ss0**2
    aa1 = 6 * yy1 - 6 * ss0 * xx1
    aa7 = 6 * rrho * yy5 - 6 * ssinf * xx3
    aa8 = -3 * ssinf**2
    aa23 = -3 * ss1**2 - (aa0 + aa1 + aa4 + aa5 + aa6 + aa7 + aa8)
    aa3 = -(
        aa0 + aa1 * llam + aa23 * llam**2 + aa4 * llam**4
        + aa5 * llam**5 + aa6 * llam**6 + aa7 * llam**7 + aa8 * llam**8
        + 3 * ssl**2
    ) / (llam**3 - llam**2)
    aa2 = aa23 - aa3
    AA = sum(value * tt**i for i, value in enumerate((aa0, aa1, aa2, aa3, aa4, aa5, aa6, aa7, aa8)))
    BB = YY1**2 - XX1**3 - AA * XX1
    DD = -16 * (4 * AA**3 + 27 * BB**2)

    rr2 = values["r2"]
    points = [coefficient_ring(0), coefficient_ring(1), llam]
    nodes = [ss0, ss1, ssl]
    FF3 = tt * (tt - 1) * (tt - llam)
    CC2 = CT(0)
    for i, (fiber_point, node) in enumerate(zip(points, nodes)):
        basis = CT(1)
        denominator = coefficient_ring(1)
        for j, other in enumerate(points):
            if i == j:
                continue
            basis *= tt - other
            denominator *= fiber_point - other
        CC2 += node * (fiber_point - rr2)**2 * basis / denominator
    XX2 = CC2 + FF3 * (
        values["p20"] + values["p21"] * tt + values["p22"] * tt**2
        + ssinf * tt**3
    )
    YY2 = FF3 * (values["c20"] + values["c21"] * tt)
    zz2 = tt - rr2
    HH2 = XX2**3 + AA * XX2 * zz2**4 + BB * zz2**6 - YY2**2

    rr3 = values["r3"]
    CC3 = ss0 * rr3**2 * (tt - llam) / (-llam)
    CC3 += ssl * (llam - rr3)**2 * tt / llam
    XX3 = CC3 + tt * (tt - llam) * (
        values["p30"] + values["p31"] * tt + values["p32"] * tt**2
        + values["p33"] * tt**3 + ssinf * tt**4
    )
    YY3 = sum(values[f"v3{i}"] * tt**i for i in range(9))
    zz3 = tt - rr3
    HH3 = XX3**3 + AA * XX3 * zz3**4 + BB * zz3**6 - YY3**2
    all_residuals = (
        [yy5, BB.derivative(tt)(1) + ss1 * AA.derivative(tt)(1)]
        + [DD[k] for k in range(21, 13, -1)]
        + [HH2[k] for k in range(19)]
        + [HH3[k] for k in range(19)]
        + [YY3(0), YY3(llam)]
    )
    return [all_residuals[i] for i in raw_nonzero_indices]

jacobian = matrix(Kp, [
    [equation.derivative(variable).subs(point) for variable in Sp.gens()]
    for equation in mod_equations
])
rank = jacobian.rank()
nullity = S.ngens() - rank
kernel = jacobian.right_kernel()

print(
    f"MW3A10LIFT|p={p}|variables={S.ngens()}|equations={len(equations)}"
    f"|jacobian_rank={rank}|tangent_dimension={nullity}",
    flush=True,
)
print(
    "MW3A10LIFT|target=" + ",".join(f"{name}:{target[name]}" for name in S.variable_names()),
    flush=True,
)
for index, tangent_vector in enumerate(kernel.basis()):
    support = [
        f"{name}:{int(value)}"
        for name, value in zip(S.variable_names(), tangent_vector)
        if value
    ]
    print(f"MW3A10LIFT|tangent={index}|" + ",".join(support), flush=True)

fixed_index = S.variable_names().index("rho")
unknown_indices = [i for i in range(S.ngens()) if i != fixed_index]
unknown_jacobian = jacobian.matrix_from_columns(unknown_indices)
row_indices = list(unknown_jacobian.transpose().pivots())
if len(row_indices) != len(unknown_indices):
    raise RuntimeError("fixing rho did not give a full-rank transverse system")
square_jacobian = unknown_jacobian.matrix_from_rows(row_indices)

if args.series_order:
    if args.series_order < 3:
        raise SystemExit("--series-order must be at least 3")
    PS = PowerSeriesRing(Kp, "z", default_prec=args.series_order)
    z = PS.gen()
    series = [PS(target[name]) for name in S.variable_names()]
    series[fixed_index] += z
    raw_jacobian_rows = []
    for equation in raw_equations:
        fraction = SF(equation)
        numerator = S(fraction.numerator())
        denominator = Sp(fraction.denominator()).subs(point)
        raw_jacobian_rows.append([
            Sp(numerator.derivative(variable)).subs(point) / denominator
            for variable in S.gens()
        ])
    raw_jacobian = matrix(Kp, raw_jacobian_rows)
    raw_unknown_jacobian = raw_jacobian.matrix_from_columns(unknown_indices)
    raw_row_indices = list(raw_unknown_jacobian.transpose().pivots())
    if len(raw_row_indices) != len(unknown_indices):
        raise RuntimeError("structured residual Jacobian lost transverse rank")
    raw_square_jacobian = raw_unknown_jacobian.matrix_from_rows(raw_row_indices)
    for degree in range(1, args.series_order):
        residuals = structured_residuals(series, PS)
        rhs = vector(Kp, [
            -PS(residuals[row])[degree]
            for row in raw_row_indices
        ])
        correction = raw_square_jacobian.solve_right(rhs)
        for column, coordinate_index in enumerate(unknown_indices):
            series[coordinate_index] += PS(correction[column]) * z**degree
    failures = []
    for row, value in enumerate(structured_residuals(series, PS)):
        value = PS(value)
        if value != 0 and value.valuation() < args.series_order:
            failures.append((row, value.valuation()))
    if failures:
        raise RuntimeError(f"formal series failed full equations: {failures}")
    print(
        f"MW3A10SERIES|order={args.series_order}|verified_equations={len(raw_equations)}",
        flush=True,
    )
    for name in ("r1", "lam", "s1", "r2", "r3"):
        value = series[S.variable_names().index(name)]
        print(f"MW3A10SERIES|{name}={value}", flush=True)

    # Search smallest bidegree relations F(rho,w)=0.  A one-dimensional
    # kernel using substantially more series coefficients than monomials is a
    # Hermite--Pade candidate, not yet a characteristic-zero equation.
    relation_ring = PolynomialRing(Kp, ("RHO", "W"))
    RHO, W = relation_ring.gens()
    rho_series = series[fixed_index]
    relations = {}
    for name in ("r1", "lam", "s1", "r2", "r3"):
        value = series[S.variable_names().index(name)]
        found = None
        for y_degree in range(1, args.relation_y_degree + 1):
            for x_degree in range(0, args.relation_x_degree + 1):
                count = (x_degree + 1) * (y_degree + 1)
                if count > args.series_order - 5:
                    continue
                monomial_series = [
                    rho_series**i * value**j
                    for j in range(y_degree + 1)
                    for i in range(x_degree + 1)
                ]
                coefficient_matrix = matrix(Kp, [
                    [term[k] for term in monomial_series]
                    for k in range(args.series_order)
                ])
                relation_kernel = coefficient_matrix.right_kernel()
                if relation_kernel.dimension() != 1:
                    continue
                relation_vector = relation_kernel.basis()[0]
                relation = sum(
                    relation_vector[j * (x_degree + 1) + i] * RHO**i * W**j
                    for j in range(y_degree + 1)
                    for i in range(x_degree + 1)
                )
                if relation.degree(W) != y_degree:
                    continue
                found = relation
                break
            if found is not None:
                break
        if found is not None:
            relations[name] = str(found)
            print(f"MW3A10SERIES|relation_{name}={found}", flush=True)

    if args.series_json:
        output = {
            "p": int(p),
            "order": int(args.series_order),
            "target_index": int(args.target),
            "parameter": "rho=16+z",
            "coefficients": {
                name: [int(series[i][degree]) for degree in range(args.series_order)]
                for i, name in enumerate(S.variable_names())
            },
            "relations": relations,
        }
        output_path = Path(args.series_json)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
        print(f"MW3A10SERIES|json={output_path}", flush=True)

if args.hensel_digits:
    if args.hensel_digits < 2:
        raise SystemExit("--hensel-digits must be at least 2")
    coordinates = [ZZ(target[name]) for name in S.variable_names()]
    modulus = ZZ(p)
    for digit in range(1, args.hensel_digits):
        selected_rhs = []
        for row in row_indices:
            value = ZZ(integral_equations[row](*coordinates))
            if value % modulus:
                raise RuntimeError(f"equation {row} is not zero modulo {modulus}")
            selected_rhs.append(Kp(-(value // modulus)))
        correction = square_jacobian.solve_right(vector(Kp, selected_rhs))
        for column, coordinate_index in enumerate(unknown_indices):
            coordinates[coordinate_index] += modulus * ZZ(correction[column])
        next_modulus = modulus * p
        failures = [
            row for row, equation in enumerate(integral_equations)
            if ZZ(equation(*coordinates)) % next_modulus
        ]
        if failures:
            raise RuntimeError(
                f"selected transverse equations do not lift full system at digit {digit + 1}: {failures}"
            )
        modulus = next_modulus
        print(f"MW3A10HENSEL|digits={digit + 1}|modulus={modulus}", flush=True)

    reconstructed = []
    for coordinate in coordinates:
        residue = ZZ(coordinate % modulus)
        try:
            reconstructed.append(residue.rational_reconstruction(modulus))
        except ArithmeticError:
            reconstructed.append(None)
    reconstruction_text = ",".join(
        f"{name}:{value}" for name, value in zip(S.variable_names(), reconstructed)
    )
    print(f"MW3A10HENSEL|reconstruction={reconstruction_text}", flush=True)
    exact = all(value is not None for value in reconstructed)
    if exact:
        exact = all(equation(*reconstructed) == 0 for equation in integral_equations)
    print(f"MW3A10HENSEL|exact_rational_point={int(exact)}", flush=True)
    if args.hensel_json:
        output = {
            "p": int(p),
            "digits": int(args.hensel_digits),
            "modulus": int(modulus),
            "fixed_coordinate": "rho",
            "target_index": int(args.target),
            "residues": {
                name: int(value % modulus)
                for name, value in zip(S.variable_names(), coordinates)
            },
        }
        output_path = Path(args.hensel_json)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
        print(f"MW3A10HENSEL|json={output_path}", flush=True)
