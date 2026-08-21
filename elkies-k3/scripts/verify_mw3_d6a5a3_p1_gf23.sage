from sage.all import *
from pathlib import Path
import argparse


ap = argparse.ArgumentParser()
ap.add_argument(
    "--candidate",
    default="rho=1,x3=2,a1=2,a2=11,x2=13,y3=17,y4=17,a5=7",
)
args = ap.parse_args()

p = 23
K = GF(p)
artifact = Path("artifacts/local/elkies-k3/mw3-d6a5a3-p1-component2/p23-p1.ms")
meta_path = artifact.with_suffix(".meta.txt")

active_names = artifact.read_text().splitlines()[0].split(",")
Q = PolynomialRing(K, active_names, order="degrevlex")
qd = Q.gens_dict()
equations = [
    Q(line.rstrip(",").replace("^", "**"))
    for line in artifact.read_text().splitlines()[2:]
    if line.strip()
]
print(f"MW3D6P1VERIFY|phase=parsed|equations={len(equations)}", flush=True)
candidate = {
    name: int(value)
    for name, value in (item.split("=", 1) for item in args.candidate.split(","))
}
if set(candidate) != set(active_names):
    raise SystemExit(f"candidate must give exactly {active_names}")
qpoint = {qd[name]: K(value) for name, value in candidate.items()}
point_values = [int(candidate[name]) % p for name in active_names]


def evaluate_at_point(equation):
    value = 0
    for exponents, coefficient in equation.dict().items():
        term = int(coefficient) % p
        for coordinate, exponent in zip(point_values, exponents):
            if exponent:
                term = (term * pow(coordinate, exponent, p)) % p
        value = (value + term) % p
    return value


residuals = [evaluate_at_point(equation) for equation in equations]
if any(residuals):
    raise RuntimeError(f"candidate fails reduced equations: {residuals}")
print("MW3D6P1VERIFY|phase=residuals", flush=True)
def gradient_at_point(equation):
    gradient = [0] * len(active_names)
    for exponents, coefficient in equation.dict().items():
        for variable_index, exponent in enumerate(exponents):
            if not exponent:
                continue
            term = (int(coefficient) * exponent) % p
            for i, power in enumerate(exponents):
                adjusted = power - (1 if i == variable_index else 0)
                if adjusted:
                    term = (term * pow(point_values[i], adjusted, p)) % p
            gradient[variable_index] = (gradient[variable_index] + term) % p
    return gradient


jacobian = matrix(K, [gradient_at_point(equation) for equation in equations])
print("MW3D6P1VERIFY|phase=jacobian", flush=True)

base_names = [
    "rho", "x1", "x2", "x3", "y2", "y3", "y4", "y5",
    "a1", "a2", "a3", "a4", "a5", "a6",
]
R = PolynomialRing(K, base_names, order="degrevlex")
d = R.gens_dict()
RF = FractionField(R)
values = {d[name]: K(value) for name, value in candidate.items()}
inside = False
for line in meta_path.read_text().splitlines():
    if line == "DERIVED":
        inside = True
        continue
    if not inside or " <- " not in line:
        continue
    name, expression = line.split(" <- ", 1)
    expression = expression.split("    #", 1)[0]
    value = K(RF(expression).subs(values))
    values[d[name]] = value

missing = [name for name in base_names if d[name] not in values]
if missing:
    raise RuntimeError(f"unreconstructed coordinates: {missing}")
print("MW3D6P1VERIFY|phase=reconstructed", flush=True)

Kt = PolynomialRing(K, "t")
t = Kt.gen()
rho = values[d["rho"]]
X = Kt([
    K(3), values[d["x1"]], values[d["x2"]], values[d["x3"]], rho**2,
])
Y = Kt([
    K(0), K(0), values[d["y2"]], values[d["y3"]],
    values[d["y4"]], values[d["y5"]], rho**3,
])
A = Kt([K(-27)] + [values[d[f"a{i}"]] for i in range(1, 7)])
B = Y**2 - X**3 - A*X
Delta = -16 * (4*A**3 + 27*B**2)
print(
    f"MW3D6P1VERIFY|phase=model|degA={A.degree()}|degB={B.degree()}|degDelta={Delta.degree()}",
    flush=True,
)
if Delta == 0:
    raise RuntimeError("candidate lies on the identically singular discriminant boundary")
if B.degree() > 9:
    raise RuntimeError(f"B has degree {B.degree()}, expected at most 9")


def valuation_at(polynomial, point):
    factor = t - point
    valuation = 0
    while polynomial != 0 and polynomial(point) == 0:
        polynomial //= factor
        valuation += 1
    return valuation


valuations = [
    valuation_at(Delta, K(0)),
    valuation_at(Delta, K(1)),
    24 - Delta.degree(),
]
if valuations != [6, 4, 8]:
    raise RuntimeError(f"wrong reducible-fiber valuations: {valuations}")
residual_delta = Delta // (t**6 * (t - 1)**4)
if residual_delta.degree() != 6 or gcd(residual_delta, residual_delta.derivative()).degree():
    raise RuntimeError("the six residual discriminant roots are not simple")
print("MW3D6P1VERIFY|phase=fibers", flush=True)


def multiplicative_component(fiber_point, node):
    """Count singular blowup centers followed by the section in the t-chart."""
    P = PolynomialRing(K, ("u", "xx", "yy"))
    u, xx, yy = P.gens()
    shifted_A_t = Kt(A(t + fiber_point))
    shifted_B_t = Kt(B(t + fiber_point))
    shifted_A = sum(P(coefficient) * u**i for i, coefficient in enumerate(shifted_A_t.list()))
    shifted_B = sum(P(coefficient) * u**i for i, coefficient in enumerate(shifted_B_t.list()))
    surface = yy**2 - (node + xx)**3 - shifted_A * (node + xx) - shifted_B
    shifted_X = Kt(X(t + fiber_point))
    shifted_Y = Kt(Y(t + fiber_point))
    section_x = Kt((shifted_X - node) // t)
    section_y = Kt(shifted_Y // t)
    surface = P(surface(u, u*xx, u*yy) // u**2)
    steps = 1
    while True:
        center_x = K(section_x(0))
        center_y = K(section_y(0))
        point = {u: K(0), xx: center_x, yy: center_y}
        gradient = [surface.derivative(variable).subs(point) for variable in (u, xx, yy)]
        if any(gradient):
            return steps
        section_x = Kt((section_x - center_x) // t)
        section_y = Kt((section_y - center_y) // t)
        transformed = surface(u, center_x + u*xx, center_y + u*yy)
        surface = P(transformed // u**2)
        steps += 1
        if steps > 5:
            raise RuntimeError("section followed too many unresolved centers")


component_I6 = multiplicative_component(K(0), K(3))
component_I4 = multiplicative_component(K(1), K(X(1)))
if component_I6 != 2 or component_I4 != 1:
    raise RuntimeError(
        f"wrong P1 components: I6={component_I6}, I4={component_I4}"
    )

# In the I2* minimal chart at infinity, P1 specializes to the finite point
# (rho^2,rho^3) on y^2=x^3.  It is smooth because rho is nonzero, hence lies
# on the identity component.
infinity_gradient = (-3 * rho**4, 2 * rho**3)
if not any(infinity_gradient):
    raise RuntimeError("P1 specializes to the I2* cusp instead of its identity component")

print(
    f"MW3D6P1VERIFY|p={p}|valuations={','.join(map(str, valuations))}"
    f"|residual_degree={residual_delta.degree()}|residual_squarefree=1",
    flush=True,
)
print(
    f"MW3D6P1VERIFY|jacobian_rank={jacobian.rank()}|variables={len(active_names)}"
    f"|tangent_dimension={len(active_names)-jacobian.rank()}"
    f"|components=I2star:0,I6:{component_I6},I4:{component_I4}",
    flush=True,
)
print("MW3D6P1VERIFY|A=" + ",".join(map(str, map(int, A.list()))), flush=True)
print("MW3D6P1VERIFY|B=" + ",".join(map(str, map(int, B.list()))), flush=True)
print("MW3D6P1VERIFY|X=" + ",".join(map(str, map(int, X.list()))), flush=True)
print("MW3D6P1VERIFY|Y=" + ",".join(map(str, map(int, Y.list()))), flush=True)
