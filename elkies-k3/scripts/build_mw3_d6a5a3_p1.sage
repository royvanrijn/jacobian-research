from sage.all import *
import argparse
from pathlib import Path


ap = argparse.ArgumentParser(
    description="Build the normalized P1 chart for the D6+A5+A3 MW3 fibration."
)
ap.add_argument("--p", type=int, default=23)
ap.add_argument("--show", action="store_true")
ap.add_argument("--export", default=None)
args = ap.parse_args()

if not is_prime(args.p) or args.p in (2, 3, 7, 79):
    raise SystemExit("choose a prime away from 2*3*7*79")
K = GF(args.p)
names = [
    "rho", "x1", "x2", "x3", "y2", "y3", "y4", "y5",
    "a1", "a2", "a3", "a4", "a5", "a6",
]
R = PolynomialRing(K, names, order="degrevlex")
d = R.gens_dict()
RF = FractionField(R)
Rt = PolynomialRing(RF, "t")
t = Rt.gen()


def V(name):
    return RF(d[name])


rho, x1, x2, x3 = [V(name) for name in ("rho", "x1", "x2", "x3")]
y2, y3, y4, y5 = [V(name) for name in ("y2", "y3", "y4", "y5")]
a1, a2, a3, a4, a5, a6 = [V(f"a{i}") for i in range(1, 7)]

# I6 is at t=0 with node (3,0).  On the first blowup, component 2 passes
# through the remaining singular center.  Its gradient equations are
# y1=0 and x1=-a1/18.
s0 = RF(3)
y1 = RF(0)
X = s0 + x1*t + x2*t**2 + x3*t**3 + rho**2*t**4
Y = y1*t + y2*t**2 + y3*t**3 + y4*t**4 + y5*t**5 + rho**3*t**6

# P1 is on the identity component of I2* at infinity: after x=u^-4 X,
# y=u^-6 Y its specialization is the smooth cusp point (rho^2,rho^3).
A = -27 + a1*t + a2*t**2 + a3*t**3 + a4*t**4 + a5*t**5 + a6*t**6
B = Y**2 - X**3 - A*X
Delta = -16 * (4*A**3 + 27*B**2)

# Put I4 at t=1.  P1 meets its component-1 node, so its specialization is
# (s1,0), where s1=X(1); the B node equation then follows from the section.
s1 = X(1)
tagged = [
    ("P1_I4_Y", Y(1)),
    ("I4_A", A(1) + 3*s1**2),
]
tagged += [(f"I6_D{k}", Delta[k]) for k in range(1, 6)]
tagged += [(f"I4_D{k}", Delta.derivative(t, k)(1)) for k in range(1, 4)]
tagged += [(f"I2star_B{k}", B[k]) for k in range(10, 13)]
tagged += [(f"I2star_D{k}", Delta[k]) for k in (18, 17)]

equations = []
for tag, expression in tagged:
    numerator = R(RF(expression).numerator())
    if numerator == 0:
        print(f"MW3D6P1_DEP|tag={tag}", flush=True)
        continue
    equations.append((tag, numerator))

# The sparse local equations give a deterministic triangular chain.  The
# component-2 blowup relation makes I6_D2 an identity after substitution.
active_names = list(names)
substitutions = {RF(d["x1"]): RF(-d["a1"] / 18)}
active_names.remove("x1")
open_factors = [RF(d["rho"]), RF(d["a1"])]
elimination_history = [
    ("x1", "I6 component-2 first blowup", RF(-d["a1"] / 18))
]


def settle(expression):
    expression = RF(expression)
    for _ in range(12):
        old = expression
        expression = RF(expression.subs(substitutions))
        if expression == old:
            break
    return expression


def eliminate(source_tag, variable_name, divide_by=None):
    source = next(expression for tag, expression in equations if tag == source_tag)
    numerator = R(settle(source).numerator())
    if divide_by is not None:
        factor = R(settle(divide_by).numerator())
        numerator, remainder = numerator.quo_rem(factor)
        if remainder:
            raise RuntimeError(f"{source_tag} lost expected factor {factor}")
        open_factors.append(RF(factor))
    variable = d[variable_name]
    if numerator.degree(variable) != 1:
        raise RuntimeError(
            f"{source_tag} has degree {numerator.degree(variable)} in {variable_name} "
            f"after prior eliminations: {numerator}"
        )
    coefficient = RF(numerator.derivative(variable))
    if R(coefficient.numerator()).degree(variable):
        raise RuntimeError(f"coefficient still depends on {variable_name}")
    rhs = RF(-numerator.subs({variable: 0}) / coefficient)
    substitutions[RF(variable)] = rhs
    elimination_history.append((variable_name, source_tag, rhs))
    open_factors.append(coefficient)
    active_names.remove(variable_name)
    print(
        f"MW3D6P1_ELIM|var={variable_name}|from={source_tag}"
        f"|degree={numerator.total_degree()}|terms={len(numerator.monomials())}",
        flush=True,
    )


eliminate("I2star_B11", "y5")
eliminate("P1_I4_Y", "y2")
eliminate("I6_D5", "a3")
eliminate("I2star_B10", "a6")
eliminate("I4_A", "a4")

used_tags = {
    "I2star_B11", "P1_I4_Y", "I6_D5", "I2star_B10",
    "I4_A",
}
reduced = []
for tag, expression in equations:
    if tag in used_tags:
        continue
    specialized = settle(expression)
    numerator = R(specialized.numerator())
    if numerator == 0:
        print(f"MW3D6P1_DEP|tag={tag}|after=triangular", flush=True)
        continue
    reduced.append((tag, numerator))

print(
    f"MW3D6P1|field=GF({args.p})|vars={len(names)}|eqs={len(equations)}"
    f"|naive_dimension={len(names)-len(equations)}",
    flush=True,
)
print(
    f"MW3D6P1_REDUCED|vars={len(active_names)}|eqs={len(reduced)}"
    f"|expected_dimension={len(active_names)-len(reduced)}"
    f"|remaining={','.join(active_names)}",
    flush=True,
)
for index, (tag, equation) in enumerate(reduced):
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
        f"MW3D6P1_EQ|i={index}|tag={tag}|degree={equation.total_degree()}"
        f"|terms={len(equation.monomials())}|linear={','.join(linear)}"
        f"|vars={','.join(support)}",
        flush=True,
    )
    if args.show:
        print(f"MW3D6P1_EXPR|i={index}|{equation}", flush=True)

# Report the actual tangent codimension at a generic point only after a point
# is found; equation count is not a proof of independence.
if args.export:
    output = Path(args.export)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        handle.write(",".join(active_names) + "\n")
        handle.write(str(args.p) + "\n")
        for i, (_, equation) in enumerate(reduced):
            handle.write(str(equation).replace("**", "^"))
            handle.write(",\n" if i + 1 < len(reduced) else "\n")
    open_output = output.with_suffix(".open.ms")
    open_polynomials = []
    for factor in open_factors:
        numerator = R(settle(factor).numerator())
        if numerator not in (R(0), R(1)):
            open_polynomials.append(numerator)
    with open_output.open("w") as handle:
        handle.write(",".join(active_names) + "\n")
        handle.write(str(args.p) + "\n")
        for i, equation in enumerate(open_polynomials):
            handle.write(str(equation).replace("**", "^"))
            handle.write(",\n" if i + 1 < len(open_polynomials) else "\n")
    meta_output = output.with_suffix(".meta.txt")
    with meta_output.open("w") as handle:
        handle.write(f"prime={args.p}\n")
        handle.write("remaining=" + repr(active_names) + "\n")
        handle.write("\nDERIVED\n")
        for name, tag, rhs in elimination_history:
            handle.write(f"{name} <- {rhs}    # {tag}\n")
    print(
        f"MW3D6P1|export={output}|open_export={open_output}|meta={meta_output}",
        flush=True,
    )
