from sage.all import *
from pathlib import Path
import argparse


ap = argparse.ArgumentParser(
    description="Build the four-variable fixed-surface P2 component-6 system over GF(31)."
)
ap.add_argument("--surface", type=int, choices=(1, 2, 3, 4), default=1)
ap.add_argument("--out", required=True)
ap.add_argument("--verify", default=None, help="assignment r=...,q0=...,q1=...,q2=...")
ap.add_argument("--solve", action="store_true", help="solve over GF(31) and test the residual square")
ap.add_argument("--A", default=None, help="custom comma-separated A coefficients")
ap.add_argument("--B", default=None, help="custom comma-separated B coefficients")
ap.add_argument("--lam", type=int, default=None)
ap.add_argument("--nodes", default=None, help="custom nodes at 0,1,lambda")
ap.add_argument("--sinf", type=int, default=None)
ap.add_argument("--X1", default=None, help="P1 X coefficients for the intersection gate")
ap.add_argument("--Y1", default=None, help="P1 Y coefficients for the intersection gate")
args = ap.parse_args()

K = GF(31)
Q = PolynomialRing(K, ("r", "q0", "q1", "q2"), order="degrevlex")
r, q0, q1, q2 = Q.gens()
Qt = PolynomialRing(Q, "t")
t = Qt.gen()

surfaces = {
    1: {
        "A": [4, 9, 5, 22, 0, 8, 17, 4, 1],
        "B": [23, 4, 26, 24, 17, 7, 14, 8, 25, 3, 15, 25, 30],
        "lam": 27,
        "nodes": [3, 7, 19],
        "sinf": 17,
    },
    2: {
        "A": [4, 28, 7, 18, 23, 26, 15, 7, 1],
        "B": [23, 9, 18, 12, 2, 15, 3, 0, 25, 18, 25, 5, 30],
        "lam": 27,
        "nodes": [3, 22, 4],
        "sinf": 17,
    },
    3: {
        "A": [4, 30, 20, 27, 13, 13, 5, 28, 25],
        "B": [23, 3, 27, 3, 21, 28, 24, 18, 13, 29, 6, 7, 30],
        "lam": 3,
        "nodes": [3, 21, 20],
        "sinf": 23,
    },
    4: {
        "A": [4, 23, 18, 28, 12, 18, 23, 20, 19],
        "B": [23, 24, 27, 10, 6, 8, 23, 26, 26, 16, 19, 9, 15],
        "lam": 23,
        "nodes": [3, 21, 10],
        "sinf": 29,
    },
}
if args.A is not None:
    if None in (args.B, args.lam, args.nodes, args.sinf):
        raise SystemExit("custom data requires --A,--B,--lam,--nodes,--sinf")
    data = {
        "A": [int(x) for x in args.A.split(",")],
        "B": [int(x) for x in args.B.split(",")],
        "lam": args.lam,
        "nodes": [int(x) for x in args.nodes.split(",")],
        "sinf": args.sinf,
    }
else:
    data = surfaces[args.surface]
A = Qt(data["A"])
B = Qt(data["B"])
lam = K(data["lam"])
sinf = K(data["sinf"])
fiber_points = [K(0), K(1), lam]
fiber_nodes = [K(value) for value in data["nodes"]]

F = t * (t - 1) * (t - lam)
C = Qt(0)
for i, (point, node) in enumerate(zip(fiber_points, fiber_nodes)):
    basis = Qt(1)
    denominator = K(1)
    for j, other in enumerate(fiber_points):
        if i == j:
            continue
        basis *= t - other
        denominator *= point - other
    C += Q(node) * (Q(point) - r)**2 * basis / Q(denominator)

# The three finite node incidences are built into C.  At infinity x6=sinf;
# component 6~=-5 has ybar order at least five, so deg(Y)<=4 and H has degree
# at most eight.
X = C + F * (q0 + q1 * t + q2 * t**2 + Q(sinf) * t**3)
z = t - r
H = X**3 + A * X * z**4 + B * z**6
assert H[18] == 0 and H[17] == 0

equations = [Q(H[k]) for k in range(16, 8, -1)]
if any(equation == 0 for equation in equations):
    raise RuntimeError("unexpected dependent high coefficient before solving")

out = Path(args.out)
out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w") as handle:
    handle.write("r,q0,q1,q2\n31\n")
    for i, equation in enumerate(equations):
        handle.write(str(equation).replace("**", "^"))
        handle.write(",\n" if i + 1 < len(equations) else "\n")

meta = out.with_suffix(".meta.txt")
with meta.open("w") as handle:
    handle.write(f"surface={args.surface}\n")
    handle.write("open=r*(r-1)*(r-27) != 0\n")
    handle.write("X=" + str(X) + "\n")
    handle.write("H_low=" + str(Qt([H[k] for k in range(9)])) + "\n")

print(
    f"MW3A10P2C6|surface={args.surface}|vars=4|eqs={len(equations)}"
    f"|out={out}|meta={meta}",
    flush=True,
)
for i, equation in enumerate(equations):
    print(
        f"MW3A10P2C6_EQ|i={i}|coefficient={16-i}"
        f"|degree={equation.total_degree()}|terms={len(equation.monomials())}",
        flush=True,
    )

if args.verify:
    assignment = {}
    for item in args.verify.split(","):
        name, value = item.split("=", 1)
        assignment[Q(name.strip())] = K(int(value.strip()))
    if set(assignment) != set(Q.gens()):
        raise RuntimeError("verify assignment must give r,q0,q1,q2")
    residual_values = [K(equation.subs(assignment)) for equation in equations]
    if any(residual_values):
        raise RuntimeError(f"high coefficient equations fail: {residual_values}")

    Kt = PolynomialRing(K, "t")
    X_special = Kt([K(coefficient.subs(assignment)) for coefficient in X.list()])
    H_special = Kt([K(coefficient.subs(assignment)) for coefficient in H.list()])
    factorization = H_special.factor()
    unit = K(factorization.unit())
    unit_roots = unit.sqrt(all=True)
    square = bool(unit_roots) and all(exponent % 2 == 0 for _, exponent in factorization)
    print(
        f"MW3A10P2C6_VERIFY|square={int(square)}|H_degree={H_special.degree()}"
        f"|factorization={factorization}",
        flush=True,
    )
    if square:
        Y_special = Kt(unit_roots[0])
        for factor, exponent in factorization:
            Y_special *= factor ** (exponent // 2)
        incidence = [Y_special(point) == 0 for point in fiber_points]
        print(
            f"MW3A10P2C6_VERIFY|incidence={','.join(str(int(x)) for x in incidence)}"
            f"|X={','.join(str(int(c)) for c in X_special.list())}"
            f"|Y={','.join(str(int(c)) for c in Y_special.list())}",
            flush=True,
        )
        if not all(incidence) or Y_special**2 != H_special:
            raise RuntimeError("square root fails finite component incidence")

if args.solve:
    solutions = Q.ideal(equations).variety()
    square_hits = 0
    rational_candidates = 0
    Kt = PolynomialRing(K, "t")
    for solution in solutions:
        pole = K(solution[r])
        if pole in fiber_points:
            continue
        rational_candidates += 1
        X_special = Kt([K(coefficient.subs(solution)) for coefficient in X.list()])
        H_special = Kt([K(coefficient.subs(solution)) for coefficient in H.list()])
        factorization = H_special.factor()
        unit = K(factorization.unit())
        roots = unit.sqrt(all=True)
        square = bool(roots) and all(exponent % 2 == 0 for _, exponent in factorization)
        if not square:
            continue
        Y_special = Kt(roots[0])
        for factor, exponent in factorization:
            Y_special *= factor ** (exponent // 2)
        if any(Y_special(point) != 0 for point in fiber_points):
            continue
        square_hits += 1
        intersection_text = ""
        if args.X1 is not None and args.Y1 is not None:
            X1_special = Kt([int(x) for x in args.X1.split(",")])
            Y1_special = Kt([int(x) for x in args.Y1.split(",")])
            z_special = Kt.gen() - pole
            x2_function = FractionField(Kt)(X_special) / FractionField(Kt)(z_special)**2
            y2_function = FractionField(Kt)(Y_special) / FractionField(Kt)(z_special)**3
            intersections = []
            for sign in (K(1), K(-1)):
                dx = Kt((FractionField(Kt)(X1_special) - x2_function).numerator())
                dy = Kt((FractionField(Kt)(Y1_special) - sign * y2_function).numerator())
                common = gcd(dx, dy)
                smooth = common
                for point in fiber_points:
                    while smooth(point) == 0:
                        smooth //= Kt.gen() - point
                # P1 and P2 share only the nonidentity component at lambda.
                # Equal first jets give one generic intersection on the
                # exceptional component; unequal jets separate after blow-up.
                tangent_equal = (
                    X1_special.derivative()(lam)
                    == x2_function.derivative()(lam)
                    and Y1_special.derivative()(lam)
                    == (sign * y2_function).derivative()(lam)
                )
                intersections.append(smooth.degree() + int(tangent_equal))
            intersection_text = (
                f"|P1int_plus={intersections[0]}|P1int_minus={intersections[1]}"
                f"|target_intersection={int(1 in intersections)}"
            )
        print(
            f"MW3A10P2C6_HIT|r={int(solution[r])}|q0={int(solution[q0])}"
            f"|q1={int(solution[q1])}|q2={int(solution[q2])}"
            f"|X={','.join(str(int(c)) for c in X_special.list())}"
            f"|Y={','.join(str(int(c)) for c in Y_special.list())}"
            f"{intersection_text}",
            flush=True,
        )
    print(
        f"MW3A10P2C6_SOLVE|candidates={rational_candidates}|square_hits={square_hits}",
        flush=True,
    )
