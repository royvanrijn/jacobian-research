from sage.all import *
import argparse
from pathlib import Path


ap = argparse.ArgumentParser(
    description="Build the normalized first-section chart for the E6+D4+2A2+A1 MW2 frame."
)
ap.add_argument("--p", type=int, default=23)
ap.add_argument("--show", action="store_true")
ap.add_argument("--export", default=None)
args = ap.parse_args()

p = args.p
if not is_prime(p) or p in (2,3,79):
    raise SystemExit("choose a good odd prime away from 3 and 79")
K = GF(p)
names = ["lam","mu","s","w","x1","x2","c1","c2","a2","a3","a4","a5"]
R = PolynomialRing(K,names,order="degrevlex")
d = R.gens_dict()
RF = FractionField(R)
Rt = PolynomialRing(RF,"t")
t = Rt.gen()

lam,mu,s,w,x1,x2,c1,c2,a2,a3,a4,a5 = [RF(d[name]) for name in names]
X = 1+x1*t+x2*t**2
# P1 meets the nonidentity components at the second I3 and at I2.  Encoding
# those two node incidences as roots keeps the model sparse and normalizes
# its smooth I0* specialization to (1,1).
Y = (t-lam)*(t-mu)*(1/(lam*mu)+c1*t+c2*t**2)
A = a2*t**2+a3*t**3+a4*t**4+a5*t**5
B = Y**2-X**3-A*X
Delta = -16*(4*A**3+27*B**2)

tagged = [
    ("I0star_B1",B[1]),
    ("I0star_B2",B[2]),
    ("I3lam_nodeA",A(lam)+3*X(lam)**2),
    ("I3lam_D1",Delta.derivative(t)(lam)),
    (
        "I3lam_order3",
        6*X(lam)**2*A.derivative(t,2)(lam)
        +6*X(lam)*B.derivative(t,2)(lam)-A.derivative(t)(lam)**2,
    ),
    ("I2mu_nodeA",A(mu)+3*X(mu)**2),
    ("I2mu_D1",Delta.derivative(t)(mu)),
    ("I3one_nodeA",A(1)+3*s**2),
    ("I3one_P1X",X(1)-(w**2-2*s)),
    ("I3one_P1Y",Y(1)-w*(w**2-3*s)),
    ("I3one_tangent",s*A.derivative(t)(1)+B.derivative(t)(1)),
    (
        "I3one_order3",
        6*s**2*A.derivative(t,2)(1)
        +6*s*B.derivative(t,2)(1)-A.derivative(t)(1)**2,
    ),
]

equations=[]
for tag,expr in tagged:
    numerator=R(RF(expr).numerator())
    if numerator==0:
        print(f"MW2P1_DEP|tag={tag}",flush=True)
    else:
        equations.append((tag,numerator))
        support=[]
        linear=[]
        for name in names:
            degree=numerator.degree(d[name])
            if degree:
                support.append(f"{name}:{degree}")
                if degree==1 and numerator.derivative(d[name]).degree(d[name])==0:
                    linear.append(name)
        print(
            f"MW2P1_RAW|tag={tag}|degree={numerator.total_degree()}"
            f"|terms={len(numerator.monomials())}|linear={','.join(linear)}"
            f"|vars={','.join(support)}",flush=True,
        )
        if args.show:
            print(f"MW2P1_EXPR|tag={tag}|{numerator}",flush=True)

print(
    f"MW2P1|p={p}|vars={len(names)}|equations={len(equations)}"
    f"|naive_dimension={len(names)-len(equations)}|expected_dimension=2",flush=True,
)

# Sparse deterministic reduction.  The I0* conditions solve c1,c2, while the
# two node-A equations are a 2-by-2 linear system for a3,a4.  Its determinant
# is supported on the open set lam*mu*(lam-mu) != 0.
substitutions = {}
elimination_order = []


def settle(expression):
    expression = RF(expression)
    for _ in range(8):
        old = expression
        expression = RF(expression.subs(substitutions))
        if expression == old:
            break
    return expression


def solve_linear(tag, variable_name):
    expression = next(expr for source,expr in tagged if source==tag)
    numerator = R(settle(expression).numerator())
    variable = d[variable_name]
    if numerator.degree(variable)!=1:
        raise RuntimeError(f"{tag} is not linear in {variable_name}")
    coefficient = RF(numerator.derivative(variable))
    rhs = RF(-numerator.subs({variable:0})/coefficient)
    substitutions[RF(variable)] = rhs
    elimination_order.append(variable_name)
    print(f"MW2P1_ELIM|var={variable_name}|from={tag}",flush=True)


solve_linear("I0star_B1","c1")
solve_linear("I0star_B2","c2")
rhs_lam = -3*X(lam)**2-a2*lam**2-a5*lam**5
rhs_mu = -3*X(mu)**2-a2*mu**2-a5*mu**5
node_matrix = matrix(RF,[[lam**3,lam**4],[mu**3,mu**4]])
node_solution = node_matrix.solve_right(vector(RF,[rhs_lam,rhs_mu]))
substitutions[RF(d["a3"])] = node_solution[0]
substitutions[RF(d["a4"])] = node_solution[1]
elimination_order.extend(["a3","a4"])
print("MW2P1_ELIM|vars=a3,a4|from=two_nodeA",flush=True)
solve_linear("I3one_P1X","x2")
solve_linear("I3one_P1Y","a2")
solve_linear("I3one_nodeA","a5")

used={
    "I0star_B1","I0star_B2","I3lam_nodeA","I2mu_nodeA",
    "I3one_nodeA","I3one_P1X","I3one_P1Y",
}
reduced=[]
active=["lam","mu","s","w","x1"]
for tag,expr in tagged:
    if tag in used:
        continue
    numerator=R(settle(expr).numerator())
    if numerator==0:
        print(f"MW2P1_DEP|tag={tag}|after=triangular",flush=True)
        continue
    reduced.append((tag,numerator))
    support=[]
    linear=[]
    for name in active:
        degree=numerator.degree(d[name])
        if degree:
            support.append(f"{name}:{degree}")
            if degree==1 and numerator.derivative(d[name]).degree(d[name])==0:
                linear.append(name)
    print(
        f"MW2P1_REDUCED_EQ|tag={tag}|degree={numerator.total_degree()}"
        f"|terms={len(numerator.monomials())}|linear={','.join(linear)}"
        f"|vars={','.join(support)}",flush=True,
    )

print(
    f"MW2P1_REDUCED|vars={len(active)}|equations={len(reduced)}"
    f"|remaining={','.join(active)}",flush=True,
)

if args.export:
    output=Path(args.export)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open("w") as handle:
        handle.write(",".join(active)+"\n")
        handle.write(str(p)+"\n")
        for index,(_,equation) in enumerate(reduced):
            handle.write(str(equation).replace("**","^"))
            handle.write(",\n" if index+1<len(reduced) else "\n")
    open_output=output.with_suffix(".open.ms")
    open_factors=[lam,mu,lam-1,mu-1,lam-mu,s,X(lam),X(mu)]
    with open_output.open("w") as handle:
        handle.write(",".join(active)+"\n")
        handle.write(str(p)+"\n")
        for index,factor in enumerate(open_factors):
            numerator=R(settle(factor).numerator())
            handle.write(str(numerator).replace("**","^"))
            handle.write(",\n" if index+1<len(open_factors) else "\n")
    meta_output=output.with_suffix(".meta.txt")
    with meta_output.open("w") as handle:
        handle.write(f"prime={p}\nremaining={active!r}\n\nDERIVED\n")
        for name in elimination_order:
            handle.write(f"{name} <- {settle(RF(d[name]))}\n")
    print(
        f"MW2P1|export={output}|open_export={open_output}|meta={meta_output}",flush=True,
    )
