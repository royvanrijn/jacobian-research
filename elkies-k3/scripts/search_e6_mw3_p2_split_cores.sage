from sage.all import *
from pathlib import Path
import argparse


ap = argparse.ArgumentParser(
    description=(
        "Reconstruct the seven exhaustive GF(31) E6/P1 split-root cores and "
        "test the reduced 837-case canonical P2 ansatz on every resulting surface."
    )
)
ap.add_argument(
    "--meta", default="artifacts/local/elkies-k3/e6-base.meta.txt",
    help="triangular P1 elimination metadata from export_e6_p1_sliced.sage",
)
args = ap.parse_args()

K = GF(31)
T = PolynomialRing(K, "t")
t = T.gen()

# This is the coordinate slice on which the four-variable core was exhausted.
r0, s0, x1 = K(4), K(18), K(27)
x0 = r0**2 - 2*s0
y0 = r0*(r0**2 - 3*s0)

cores = [
    (2, 12, 25, 30),
    (4, 16, 6, 23),
    (14, 22, 30, 21),
    (15, 4, 17, 30),
    (17, 23, 4, 19),
    (21, 27, 0, 4),
    (23, 11, 8, 5),
]

names = [
    "a1", "a2", "a3", "a4", "b4", "b5", "lam", "mu", "s0", "s1",
    "sl", "sm", "x0", "x1", "y0", "y1", "y2",
]
R = PolynomialRing(K, names)
d = R.gens_dict()
F = FractionField(R)

metadata = {}
for line in Path(args.meta).read_text().splitlines():
    if "=" in line:
        key, value = line.split("=", 1)
        metadata[key.strip()] = value.strip()
for key in ("a3", "y2", "b5", "b4"):
    if key not in metadata:
        raise RuntimeError(f"metadata lacks {key}")


def evaluate_expression(text, values):
    expression = F(text.replace("^", "**"))
    substitutions = {d[key]: F(value) for key, value in values.items() if key in d}
    for _ in range(8):
        old = expression
        expression = F(expression.subs(substitutions))
        if expression == old:
            break
    if expression.denominator() == 0:
        raise ZeroDivisionError(text)
    try:
        return K(expression)
    except (TypeError, ValueError) as error:
        raise RuntimeError(
            f"expression did not fully specialize: {expression}"
        ) from error


def reconstruct(core):
    a1, a2, a4, s1 = map(K, core)
    values = {
        "r0": r0, "s0": s0, "x1": x1, "x0": x0, "y0": y0,
        "a1": a1, "a2": a2, "a4": a4, "s1": s1,
    }
    values["y1"] = (a1 + 3*(r0**2-s0)*x1)/(2*r0)
    for key in ("a3", "y2", "b5", "b4"):
        values[key] = evaluate_expression(metadata[key], values)

    a3, y2, b5, b4 = (values[key] for key in ("a3", "y2", "b5", "b4"))
    a0 = -3*s0**2
    b0 = 2*s0**3
    b1 = -s0*a1
    b2 = a1**2/(12*s0) - s0*a2
    b3 = (a1**3 + 36*a1*a2*s0**2 - 216*a3*s0**4)/(216*s0**3)
    a5 = -3*s1**2 - (a0+a1+a2+a3+a4)

    U = PolynomialRing(K, "B6")
    B6 = U.gen()
    Ut = PolynomialRing(U, "u")
    u = Ut.gen()
    AA = sum(U(c)*u**i for i, c in enumerate([a0,a1,a2,a3,a4,a5]))
    B7 = 2*U(s1)**3 - (
        U(b0)+U(b1)+U(b2)+U(b3)+U(b4)+U(b5)+B6
    ) - 1
    BB = (
        U(b0)+U(b1)*u+U(b2)*u**2+U(b3)*u**3+U(b4)*u**4
        + U(b5)*u**5+B6*u**6+B7*u**7+u**8
    )
    equation = U(BB.derivative(u)(1) + U(s1)*AA.derivative(u)(1))
    coefficient = equation.derivative(B6)
    if coefficient == 0 or coefficient.derivative(B6) != 0:
        raise RuntimeError("b6 equation is not affine-linear")
    b6 = K(-equation.subs({B6: 0})/coefficient)
    b7 = K(B7.subs({B6: U(b6)}))

    A = sum(T(c)*t**i for i, c in enumerate([a0,a1,a2,a3,a4,a5]))
    B = sum(T(c)*t**i for i, c in enumerate([b0,b1,b2,b3,b4,b5,b6,b7,K(1)]))
    X1 = T(x0) + T(x1)*t + T(s1-x0-x1)*t**2
    Y1 = T(y0) + T(values["y1"])*t + T(y2)*t**2
    Y1 += T(-1-y0-values["y1"]-y2)*t**3 + t**4
    if Y1**2 - X1**3 - A*X1 - B != 0:
        raise RuntimeError("reconstructed P1 identity failed")
    return A, B, X1, Y1


def valuation_at(poly, point):
    multiplicity = 0
    current = poly
    while current != 0 and current(point) == 0:
        multiplicity += 1
        current = current.derivative(t)
    return multiplicity


def singular_nodes(A, B):
    derivative_A = A.derivative(t)
    derivative_B = B.derivative(t)
    nodes = {}
    for point in K:
        matches = [
            node for node in K
            if A(point) == -3*node**2
            and B(point) == 2*node**3
            and derivative_B(point) + node*derivative_A(point) == 0
        ]
        if matches:
            nodes[point] = matches
    return nodes


def square_roots_by_coefficients(poly):
    if poly == 0:
        return [T(0)]
    degree = poly.degree()
    if degree % 2:
        return []
    root_degree = degree // 2
    roots = []
    for leading in K(poly[degree]).sqrt(all=True):
        coefficients = [K(0)]*(root_degree+1)
        coefficients[root_degree] = leading
        for target_degree in range(degree-1, root_degree-1, -1):
            index = target_degree-root_degree
            known = K(0)
            for i in range(index+1, root_degree+1):
                j = target_degree-i
                if 0 <= j <= root_degree and j != index:
                    known += coefficients[i]*coefficients[j]
            coefficients[index] = (K(poly[target_degree])-known)/(2*leading)
        candidate = T(coefficients)
        if candidate**2 == poly:
            roots.append(candidate)
    return roots


def interpolate_x(points, nodes, pole):
    result = T(0)
    for index, (point, node) in enumerate(zip(points, nodes)):
        basis = T(1)
        denominator = K(1)
        for other_index, other in enumerate(points):
            if index == other_index:
                continue
            basis *= t-other
            denominator *= point-other
        result += node*(point-pole)**2*basis/denominator
    return result


def search_p2(A, B, lam, mu, sl, sm):
    points = [K(0), K(1), lam, mu]
    nodes = [s0, None, sl, sm]
    # The I4 node at t=1 is the core parameter and is recovered from A,B.
    nodes[1] = next(
        node for node in K
        if A(1) == -3*node**2 and B(1) == 2*node**3
    )
    fiber_product = prod(t-point for point in points)
    tested = 0
    hits = []
    for pole in K:
        if pole in points:
            continue
        constant = interpolate_x(points, nodes, pole)
        z = t-pole
        for q0 in K:
            tested += 1
            X2 = constant + q0*fiber_product
            if X2(pole) == 0:
                continue
            H = X2**3 + A*X2*z**4 + B*z**6
            for Y2 in square_roots_by_coefficients(H):
                if any(Y2(point) != 0 for point in points):
                    continue
                hits.append((pole, q0, X2, Y2))
    return tested, hits


total_surfaces = 0
total_hits = 0
for core_index, core in enumerate(cores, start=1):
    try:
        A, B, X1, Y1 = reconstruct(core)
    except ZeroDivisionError:
        print(
            f"E6P2CORE|core={core_index}|parameters={','.join(map(str, core))}"
            "|status=OUTSIDE_TRIANGULAR_CHART",
            flush=True,
        )
        continue
    Delta = -16*(4*A**3 + 27*B**2)
    node_map = singular_nodes(A, B)
    repeated = sorted(
        point for point in K
        if point not in (K(0), K(1)) and valuation_at(Delta, point) >= 2
    )
    print(
        f"E6P2CORE|core={core_index}|parameters={','.join(map(str, core))}"
        f"|A={','.join(str(int(c)) for c in A.list())}"
        f"|B={','.join(str(int(c)) for c in B.list())}"
        f"|X1={','.join(str(int(c)) for c in X1.list())}"
        f"|Y1={','.join(str(int(c)) for c in Y1.list())}"
        f"|mult0={valuation_at(Delta, K(0))}|mult1={valuation_at(Delta, K(1))}"
        f"|repeated={','.join(str(int(point)) for point in repeated)}",
        flush=True,
    )
    if valuation_at(Delta, K(0)) != 4 or valuation_at(Delta, K(1)) != 4:
        print(f"E6P2CORE|core={core_index}|status=BOUNDARY", flush=True)
        continue

    for lam_index in range(len(repeated)):
        for mu_index in range(lam_index+1, len(repeated)):
            lam, mu = repeated[lam_index], repeated[mu_index]
            for sl in node_map.get(lam, []):
                for sm in node_map.get(mu, []):
                    total_surfaces += 1
                    tested, hits = search_p2(A, B, lam, mu, sl, sm)
                    total_hits += len(hits)
                    print(
                        f"E6P2SURFACE|core={core_index}|lambda={int(lam)}|mu={int(mu)}"
                        f"|sl={int(sl)}|sm={int(sm)}|tested={tested}|hits={len(hits)}",
                        flush=True,
                    )
                    for pole, q0, X2, Y2 in hits:
                        print(
                            f"E6P2_HIT|core={core_index}|lambda={int(lam)}|mu={int(mu)}"
                            f"|r={int(pole)}|q0={int(q0)}"
                            f"|X={','.join(str(int(c)) for c in X2.list())}"
                            f"|Y={','.join(str(int(c)) for c in Y2.list())}",
                            flush=True,
                        )

print(
    f"E6P2SUMMARY|cores={len(cores)}|surfaces={total_surfaces}|hits={total_hits}",
    flush=True,
)
