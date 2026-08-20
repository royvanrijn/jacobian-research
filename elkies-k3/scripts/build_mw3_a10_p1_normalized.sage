from sage.all import *
from pathlib import Path
import argparse


ap = argparse.ArgumentParser(
    description=(
        "Build the normalized twelve-variable P1 chart for the semistable "
        "A10+A2+A1^2 MW3 fibration."
    )
)
ap.add_argument("--p", type=int, default=31, help="finite-field prime; use 0 for QQ")
ap.add_argument("--export", default=None, help="optional msolve-format output (requires --p)")
ap.add_argument("--show", action="store_true", help="print the nine residual equations")
ap.add_argument(
    "--verify",
    default=None,
    help="comma-separated active-coordinate assignment to reconstruct and verify",
)
ap.add_argument(
    "--stage",
    choices=("base", "triangular", "component2"),
    default="component2",
    help=(
        "stop at the 12x9 chart, apply the sparse infinity eliminations, "
        "or also impose the target I11 component-2 first jet"
    ),
)
args = ap.parse_args()

if args.p:
    if not is_prime(args.p) or args.p in (2, 3, 11, 79):
        raise SystemExit("choose a good prime not in {2,3,11,79}")
    K = GF(args.p)
else:
    K = QQ

# The I3 tangent branch and the remaining Weierstrass scaling are normalized by
#
#   s0 = 3,
#   a1 + 6*s0*x1 = 6*y1.
#
# At infinity the analogous tangent parameter is rho:
#
#   sinf = 3*rho^2,
#   a7 + 6*sinf*x3 = 6*rho*y5.
#
# The point P1(1) is parametrized on the nodal cubic.  Incidence at lambda and
# A(1)=A(lambda)=-3*s^2 then eliminate y1,y2,sl,a2,a3 and x1.
names = [
    "rho", "r1", "s1", "lam",
    "x2", "x3", "y3", "y4", "y5",
    "a4", "a5", "a6",
]
R = PolynomialRing(K, names, order="degrevlex")
d = R.gens_dict()
RF = FractionField(R)
Rt = PolynomialRing(RF, "t")
t = Rt.gen()


def V(name):
    return RF(d[name])


rho, r1, s1, lam = [V(n) for n in ("rho", "r1", "s1", "lam")]
x2, x3 = [V(n) for n in ("x2", "x3")]
y3, y4, y5 = [V(n) for n in ("y3", "y4", "y5")]
a4, a5, a6 = [V(n) for n in ("a4", "a5", "a6")]

s0 = RF(3)
sinf = 3 * rho**2

# Parametrize the nonsingular point P1(1) on
# y^2=(x-s1)^2*(x+2*s1).  The node is excluded by r1^2 != 3*s1.
X_at_1 = r1**2 - 2 * s1
Y_at_1 = r1 * (r1**2 - 3 * s1)

x1 = X_at_1 - s0 - x2 - x3 - sinf

# Solve Y(1)=Y_at_1 and Y(lambda)=0 for y1,y2.  The denominator records the
# open condition lambda*(lambda-1) != 0.
y2 = -(
    lam * Y_at_1
    + (lam**3 - lam) * y3
    + (lam**4 - lam) * y4
    + (lam**5 - lam) * y5
) / (lam**2 - lam)
y1 = Y_at_1 - y2 - y3 - y4 - y5

X = s0 + x1 * t + x2 * t**2 + x3 * t**3 + sinf * t**4
Y = y1 * t + y2 * t**2 + y3 * t**3 + y4 * t**4 + y5 * t**5
sl = X(lam)

# The two local tangent-cone parametrizations make ord_0(Delta)>=3 and the
# first three infinity cancellations identities.
a0 = -3 * s0**2
a1 = 6 * y1 - 6 * s0 * x1
a7 = 6 * rho * y5 - 6 * sinf * x3
a8 = -3 * sinf**2

# Solve A(1)=-3*s1^2 and A(lambda)=-3*sl^2 for a2,a3.
a2_plus_a3 = -3 * s1**2 - (a0 + a1 + a4 + a5 + a6 + a7 + a8)
a3 = -(
    a0
    + a1 * lam
    + a2_plus_a3 * lam**2
    + a4 * lam**4
    + a5 * lam**5
    + a6 * lam**6
    + a7 * lam**7
    + a8 * lam**8
    + 3 * sl**2
) / (lam**3 - lam**2)
a2 = a2_plus_a3 - a3

aa = [a0, a1, a2, a3, a4, a5, a6, a7, a8]
A = sum(aa[i] * t**i for i in range(9))

# Section-first representation: B is defined by P1, so all thirteen section
# coefficient equations disappear identically.
B = Y**2 - X**3 - A * X
Delta = -16 * (4 * A**3 + 27 * B**2)

# Exact identities behind the reduction.  These are deliberately checked in
# the fraction field before any residual numerator is expanded.
assert X(0) == s0 and Y(0) == 0
assert A(0) == -3 * s0**2 and B(0) == 2 * s0**3
assert Delta[0] == 0 and Delta[1] == 0 and Delta[2] == 0

assert X(1) == X_at_1 and Y(1) == Y_at_1
assert A(1) == -3 * s1**2 and B(1) == 2 * s1**3

assert X(lam) == sl and Y(lam) == 0
assert A(lam) == -3 * sl**2 and B(lam) == 2 * sl**3
assert Delta(lam) == 0 and Delta.derivative(t)(lam) == 0

assert X[4] == sinf and Y.degree() <= 5
assert A[8] == -3 * sinf**2 and B[12] == 2 * sinf**3
assert Delta[24] == 0 and Delta[23] == 0 and Delta[22] == 0

# Only nine equations remain: the I2 derivative condition at t=1 and the
# eight further cancellations giving I11 at infinity.
residuals = [("I2_1_first", B.derivative(t)(1) + s1 * A.derivative(t)(1))]
residuals += [(f"I11_inf_D{k}", Delta[k]) for k in range(21, 13, -1)]

triangular_history = []
component_history = []
open_factors = [rho, lam, lam - 1, s1, sl, r1**2 - 3 * s1]
active_names = list(names)
if args.stage in ("triangular", "component2"):
    substitutions = {}

    def settle(equation):
        equation = RF(equation)
        for _ in range(6):
            old = equation
            equation = RF(equation.subs(substitutions))
            if equation == old:
                break
        return equation

    for variable_name, source_tag in (
        ("a6", "I11_inf_D21"),
        ("a5", "I11_inf_D20"),
        ("a4", "I11_inf_D19"),
    ):
        source = next(equation for tag, equation in residuals if tag == source_tag)
        source = settle(source)
        variable = d[variable_name]
        numerator = R(source.numerator())
        if numerator.degree(variable) != 1:
            raise RuntimeError(
                f"{source_tag} is not affine-linear in {variable_name}"
            )
        coefficient = RF(numerator.derivative(variable))
        if R(coefficient.numerator()).degree(variable) != 0:
            raise RuntimeError(f"coefficient of {variable_name} is not independent")
        rhs = RF(-numerator.subs({variable: 0}) / coefficient)
        substitutions[RF(variable)] = rhs
        active_names.remove(variable_name)
        triangular_history.append((variable_name, source_tag, rhs))
        print(
            f"MW3A10P1_ELIM|var={variable_name}|from={source_tag}"
            f"|degree={numerator.total_degree()}|terms={len(numerator.monomials())}",
            flush=True,
        )

    used_tags = {source_tag for _, source_tag, _ in triangular_history}
    residuals = [
        (tag, settle(equation))
        for tag, equation in residuals
        if tag not in used_tags
    ]

if args.stage == "component2":
    # In local split-node coordinates uv=t^11, component class 2 (up to sign)
    # has both branches vanishing beyond first order.  After the tangent-cone
    # parametrization above this first deeper incidence is exactly y5=0.  The
    # D18 equation then becomes an identity.  Exact class 2 additionally needs
    # a nonzero second-jet open condition, checked after finding a point.
    y5_variable = d["y5"]
    active_names.remove("y5")
    component_residuals = []
    for tag, equation in residuals:
        specialized = RF(equation.subs({y5_variable: 0}))
        if specialized == 0:
            print(f"MW3A10P1_DEP|tag={tag}|after=y5=0", flush=True)
            continue
        component_residuals.append((tag, specialized))
    residuals = component_residuals

    # The next cancellation is y4*F=0.  The branch y4=0 continues to a deeper
    # component, while exact class 2 has y4 != 0.  Divide by that open factor
    # and use the sparse quotient to eliminate y3.
    source = next(equation for tag, equation in residuals if tag == "I11_inf_D17")
    source_numerator = R(source.numerator())
    quotient, remainder = source_numerator.quo_rem(R(d["y4"]))
    if remainder != 0:
        raise RuntimeError("D17 lost its expected exact-class-2 y4 factor")
    y3_variable = d["y3"]
    if quotient.degree(y3_variable) != 1:
        raise RuntimeError("D17/y4 is not affine-linear in y3")
    coefficient = RF(quotient.derivative(y3_variable))
    y3_rhs = RF(-quotient.subs({y3_variable: 0}) / coefficient)
    open_factors.extend((RF(d["y4"]), coefficient))
    substitutions[RF(y3_variable)] = y3_rhs
    active_names.remove("y3")
    component_history.append(("y3", "I11_inf_D17/y4", y3_rhs))
    print(
        f"MW3A10P1_ELIM|var=y3|from=I11_inf_D17/y4"
        f"|degree={quotient.total_degree()}|terms={len(quotient.monomials())}",
        flush=True,
    )
    residuals = [
        (tag, settle(equation))
        for tag, equation in residuals
        if tag != "I11_inf_D17"
    ]

# Export the small open factors separately.  Multiplying them would create a
# large polynomial for no benefit.  Delta_13 != 0 and the earlier sparse
# triangular coefficients are checked only on the tiny reconstructed survivor
# set; rho != 0 already covers their generic leading factors.
open_polynomials = []
for factor in open_factors:
    factor = RF(factor)
    if args.stage == "component2":
        factor = RF(factor.subs({d["y5"]: 0}))
    polynomial = R(factor.numerator())
    if polynomial not in (R(0), R(1)):
        open_polynomials.append(polynomial)

reduced = []
for tag, equation in residuals:
    equation = RF(equation)
    numerator = R(equation.numerator())
    if tag == "I2_1_first":
        # P1(1) is parametrized on the nodal cubic.  The derivative equation
        # includes the node branch r1^2=3*s1, which is explicitly excluded.
        node_factor = R(d["r1"]**2 - 3 * d["s1"])
        quotient, remainder = numerator.quo_rem(node_factor)
        if remainder != 0:
            raise RuntimeError("I2_1_first lost its expected node factor")
        numerator = quotient
    if numerator == 0:
        raise RuntimeError(f"unexpected dependent residual equation: {tag}")
    reduced.append((tag, numerator))

field = "QQ" if args.p == 0 else f"GF({args.p})"
print(
    f"MW3A10P1|field={field}|stage={args.stage}"
    f"|vars={len(active_names)}|eqs={len(reduced)}"
    f"|expected_dim={len(active_names)-len(reduced)}",
    flush=True,
)
print("MW3A10P1|remaining=" + ",".join(active_names), flush=True)
print(
    "MW3A10P1|open="
    "rho*lam*(lam-1)*s1*sl*(r1^2-3*s1)*Delta13 != 0",
    flush=True,
)
print(
    "MW3A10P1|identities="
    "section;I3_order3;I2_lambda_order2;I11_first3",
    flush=True,
)

for i, (tag, equation) in enumerate(reduced):
    linear = []
    support = []
    for name in active_names:
        variable = d[name]
        degree = equation.degree(variable)
        if degree:
            support.append(f"{name}:{degree}")
            if degree == 1 and equation.derivative(variable).degree(variable) == 0:
                linear.append(name)
    print(
        f"MW3A10P1_EQ|i={i}|tag={tag}|degree={equation.total_degree()}"
        f"|terms={len(equation.monomials())}|linear={','.join(linear)}"
        f"|vars={','.join(support)}",
        flush=True,
    )
    if args.show:
        print(f"MW3A10P1_EXPR|i={i}|{equation}", flush=True)

if args.verify:
    if args.p == 0:
        raise SystemExit("--verify currently requires a finite field")
    assignments = {}
    for item in args.verify.split(","):
        name, value = item.split("=", 1)
        assignments[name.strip()] = K(int(value.strip()))
    if set(assignments) != set(active_names):
        raise RuntimeError(
            f"verify assignment must give exactly {active_names}; got {sorted(assignments)}"
        )
    values = {d[name]: value for name, value in assignments.items()}
    if args.stage == "component2":
        values[d["y5"]] = K(0)

    def evaluate_fraction(expression):
        value = RF(expression).subs(values)
        return K(value)

    # The component quotient is already expressed after all triangular
    # substitutions, while some stored triangular right sides still contain
    # y3.  Evaluate the component coordinate first.
    for name, _, rhs in component_history + triangular_history:
        values[d[name]] = evaluate_fraction(rhs)

    residual_values = [K(equation.subs(values)) for _, equation in reduced]
    open_values = [K(polynomial.subs(values)) for polynomial in open_polynomials]
    print(
        "MW3A10VERIFY|residuals=" + ",".join(str(int(value)) for value in residual_values),
        flush=True,
    )
    print(
        "MW3A10VERIFY|opens=" + ",".join(str(int(value)) for value in open_values),
        flush=True,
    )
    if any(residual_values) or any(value == 0 for value in open_values):
        raise RuntimeError("candidate fails a reduced equation or chart open condition")

    Kt = PolynomialRing(K, "t")

    def specialize_polynomial(poly):
        return Kt([evaluate_fraction(coefficient) for coefficient in poly.list()])

    A_special = specialize_polynomial(A)
    X_special = specialize_polynomial(X)
    Y_special = specialize_polynomial(Y)
    B_special = Y_special**2 - X_special**3 - A_special * X_special
    Delta_special = -16 * (4 * A_special**3 + 27 * B_special**2)

    if Y_special**2 != X_special**3 + A_special * X_special + B_special:
        raise RuntimeError("reconstructed section identity failed")

    def valuation_at(poly, point):
        factor = Kt.gen() - point
        valuation = 0
        while poly != 0 and poly(point) == 0:
            poly = poly // factor
            valuation += 1
        return valuation

    lam_value = assignments["lam"]
    valuations = {
        "zero": valuation_at(Delta_special, K(0)),
        "one": valuation_at(Delta_special, K(1)),
        "lambda": valuation_at(Delta_special, lam_value),
        "infinity": 24 - Delta_special.degree(),
    }
    expected = {"zero": 3, "one": 2, "lambda": 2, "infinity": 11}
    print(
        "MW3A10VERIFY|valuations="
        + ",".join(f"{name}:{value}" for name, value in valuations.items()),
        flush=True,
    )
    if valuations != expected:
        raise RuntimeError(f"wrong exact reducible-fiber valuations: {valuations}")

    fixed_factor = (
        Kt.gen()**3
        * (Kt.gen() - 1)**2
        * (Kt.gen() - lam_value)**2
    )
    residual_delta = Delta_special // fixed_factor
    residual_squarefree = gcd(residual_delta, residual_delta.derivative()).degree() == 0
    print(
        f"MW3A10VERIFY|residual_degree={residual_delta.degree()}"
        f"|residual_squarefree={int(residual_squarefree)}",
        flush=True,
    )
    if residual_delta.degree() != 6 or not residual_squarefree:
        raise RuntimeError("six residual I1 fibers are not squarefree")

    jacobian = matrix(
        K,
        [
            [equation.derivative(d[name]).subs(values) for name in active_names]
            for _, equation in reduced
        ],
    )
    print(
        f"MW3A10VERIFY|jacobian_rank={jacobian.rank()}"
        f"|expected={len(reduced)}",
        flush=True,
    )
    if jacobian.rank() != len(reduced):
        raise RuntimeError("candidate is singular on the reduced chart")

    def coeff_string(poly):
        return ",".join(str(int(coefficient)) for coefficient in poly.list())

    print(f"MW3A10VERIFY|A_coeffs={coeff_string(A_special)}", flush=True)
    print(f"MW3A10VERIFY|B_coeffs={coeff_string(B_special)}", flush=True)
    print(f"MW3A10VERIFY|X_coeffs={coeff_string(X_special)}", flush=True)
    print(f"MW3A10VERIFY|Y_coeffs={coeff_string(Y_special)}", flush=True)

if args.export:
    if args.p == 0:
        raise SystemExit("--export requires a finite field")
    out = Path(args.export)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        handle.write(",".join(active_names) + "\n")
        handle.write(str(args.p) + "\n")
        for i, (_, equation) in enumerate(reduced):
            handle.write(str(equation).replace("**", "^"))
            handle.write(",\n" if i + 1 < len(reduced) else "\n")

    meta = out.with_suffix(".meta.txt")
    with meta.open("w") as handle:
        handle.write(f"prime={args.p}\n")
        handle.write("remaining=" + repr(active_names) + "\n")
        handle.write("expected_dimension=3\n")
        handle.write(
            "open=rho*lam*(lam-1)*s1*sl*(r1^2-3*s1)*Delta13 != 0\n"
        )
        handle.write("branch=I3 tangent sign normalized positive\n")
        handle.write("\nDERIVED\n")
        for name, value in [
            ("s0", s0), ("sinf", sinf), ("x1", x1),
            ("y1", y1), ("y2", y2), ("sl", sl),
            ("a0", a0), ("a1", a1), ("a2", a2), ("a3", a3),
            ("a7", a7), ("a8", a8),
        ]:
            handle.write(f"{name} <- {value}\n")
        for name, tag, value in triangular_history:
            handle.write(f"{name} <- {value}    # {tag}\n")
        for name, tag, value in component_history:
            handle.write(f"{name} <- {value}    # {tag}\n")
        handle.write("\nEQUATIONS\n")
        for tag, equation in reduced:
            handle.write(f"{tag}: {equation}\n")
    open_out = out.with_suffix(".open.ms")
    with open_out.open("w") as handle:
        handle.write(",".join(active_names) + "\n")
        handle.write(str(args.p) + "\n")
        for i, polynomial in enumerate(open_polynomials):
            handle.write(str(polynomial).replace("**", "^"))
            handle.write(",\n" if i + 1 < len(open_polynomials) else "\n")
    print(
        f"MW3A10P1|export={out}|meta={meta}|open_export={open_out}",
        flush=True,
    )
