from sage.all import *
from pathlib import Path
import argparse
import ast
import re


ap = argparse.ArgumentParser(
    description="Validate P1 scan logs and fuse the exact P2/P3 searches."
)
ap.add_argument("--glob", default="p1-seed*.scan.log")
ap.add_argument("--max-surfaces", type=int, default=0)
ap.add_argument("--show-p2", action="store_true")
ap.add_argument("--jacobian", action="store_true")
ap.add_argument("--progress-every", type=int, default=100)
args = ap.parse_args()

p = 23
K = GF(p)
Kt = PolynomialRing(K, "t")
t = Kt.gen()
root = Path("artifacts/local/elkies-k3/mw3-d6a5a3-p1-component2")
artifact = root / "p23-p1.ms"
chart_meta = artifact.with_suffix(".meta.txt")
active_names = artifact.read_text().splitlines()[0].split(",")
Q = PolynomialRing(K, active_names, order="degrevlex")
chart_equations = [
    Q(line.rstrip(",").replace("^","**"))
    for line in artifact.read_text().splitlines()[2:] if line.strip()
]
base_names = [
    "rho", "x1", "x2", "x3", "y2", "y3", "y4", "y5",
    "a1", "a2", "a3", "a4", "a5", "a6",
]
R = PolynomialRing(K, base_names, order="degrevlex")
d = R.gens_dict()
RF = FractionField(R)

derived = []
inside = False
for line in chart_meta.read_text().splitlines():
    if line == "DERIVED":
        inside = True
        continue
    if inside and " <- " in line:
        name, expression = line.split(" <- ", 1)
        derived.append((name, RF(expression.split("    #", 1)[0])))


def valuation_at(poly, point):
    factor = t-point
    value = 0
    while poly and poly(point) == 0:
        poly //= factor
        value += 1
    return value


def multiplicative_steps_polynomial(A, B, X, Y, fiber_point, node):
    P = PolynomialRing(K, ("u", "xx", "yy"))
    u, xx, yy = P.gens()
    shifted_A_t = Kt(A(t+fiber_point))
    shifted_B_t = Kt(B(t+fiber_point))
    shifted_A = sum(P(c)*u**i for i, c in enumerate(shifted_A_t.list()))
    shifted_B = sum(P(c)*u**i for i, c in enumerate(shifted_B_t.list()))
    surface = yy**2-(node+xx)**3-shifted_A*(node+xx)-shifted_B
    shifted_X = Kt(X(t+fiber_point))
    shifted_Y = Kt(Y(t+fiber_point))
    if shifted_X(0) != node or shifted_Y(0) != 0:
        return 0
    section_x = Kt((shifted_X-node)//t)
    section_y = Kt(shifted_Y//t)
    surface = P(surface(u, u*xx, u*yy)//u**2)
    steps = 1
    while True:
        cx, cy = K(section_x(0)), K(section_y(0))
        point = {u: K(0), xx: cx, yy: cy}
        if any(surface.derivative(v).subs(point) for v in (u, xx, yy)):
            return steps
        section_x = Kt((section_x-cx)//t)
        section_y = Kt((section_y-cy)//t)
        surface = P(surface(u, cx+u*xx, cy+u*yy)//u**2)
        steps += 1
        if steps > 6:
            return 99


def fast_square_roots(H, root_degree):
    if H.degree() != 2*root_degree:
        return []
    roots = []
    for leading in K(H[2*root_degree]).sqrt(all=True):
        coeffs = [K(0)]*(root_degree+1)
        coeffs[root_degree] = leading
        for degree in range(2*root_degree-1, root_degree-1, -1):
            unknown = degree-root_degree
            known = K(0)
            for i in range(root_degree+1):
                j = degree-i
                if 0 <= j <= root_degree and i != root_degree and j != root_degree:
                    known += coeffs[i]*coeffs[j]
            coeffs[unknown] = (H[degree]-known)/(2*leading)
        root_poly = Kt(coeffs)
        if root_poly**2 == H:
            roots.append(root_poly)
    return roots


def reconstruct(candidate):
    values = {d[name]: K(value) for name, value in candidate.items()}
    try:
        for name, expression in derived:
            values[d[name]] = K(expression.subs(values))
    except (ZeroDivisionError, ValueError):
        return None
    if any(d[name] not in values for name in base_names):
        return None
    rho = values[d["rho"]]
    X = Kt([3, values[d["x1"]], values[d["x2"]], values[d["x3"]], rho**2])
    Y = Kt([0, 0, values[d["y2"]], values[d["y3"]], values[d["y4"]], values[d["y5"]], rho**3])
    A = Kt([19] + [values[d[f"a{i}"]] for i in range(1, 7)])
    B = Y**2-X**3-A*X
    if B.degree() > 9 or Y**2 != X**3+A*X+B:
        return None
    Delta = -16*(4*A**3+27*B**2)
    if not Delta or [valuation_at(Delta,K(0)),valuation_at(Delta,K(1)),24-Delta.degree()] != [6,4,8]:
        return None
    residual = Delta//(t**6*(t-1)**4)
    if residual.degree() != 6 or gcd(residual,residual.derivative()).degree():
        return None
    node1 = X(K(1))
    if multiplicative_steps_polynomial(A,B,X,Y,K(0),K(3)) != 2:
        return None
    if multiplicative_steps_polynomial(A,B,X,Y,K(1),node1) != 1:
        return None
    return A, B, X, Y


def chart_jacobian_rank(candidate):
    coordinates = [int(candidate[name]) % p for name in active_names]
    rows = []
    for equation in chart_equations:
        gradient = [0]*len(active_names)
        for exponents, coefficient in equation.dict().items():
            for variable_index, exponent in enumerate(exponents):
                if not exponent:
                    continue
                term = (int(coefficient)*exponent) % p
                for i, power in enumerate(exponents):
                    adjusted = power-(1 if i==variable_index else 0)
                    if adjusted:
                        term = term*pow(coordinates[i],adjusted,p) % p
                gradient[variable_index] = (gradient[variable_index]+term) % p
        rows.append(gradient)
    return matrix(K,rows).rank()


def find_p2(A, B, P1X, P1Y):
    node1 = P1X(K(1))
    z = polygen(K, "z")
    vector_roots = [
        K(c) for c, multiplicity in (z**3+A[6]*z+B[9]).roots()
        if 3*K(c)**2+A[6] != 0
    ]
    for c in vector_roots:
        for q1 in K:
            for q2 in K:
                X = Kt([3,q1,q2,c])
                H = X**3+A*X+B
                roots = fast_square_roots(H,4)
                for Y in roots[:1]:
                    if Y(0) != 0:
                        continue
                    for signed_Y in (Y,-Y):
                        if X(1) == node1 and signed_Y(1) == 0:
                            continue
                        if multiplicative_steps_polynomial(A,B,X,signed_Y,K(0),K(3)) != 1:
                            continue
                        if gcd(P1X-X,P1Y-signed_Y).degree() != 1:
                            continue
                        yield X, signed_Y


def local_polynomial(numerator, denominator_power, lam, point, precision=8):
    S = PowerSeriesRing(K, "u", default_prec=precision)
    u = S.gen()
    num = sum(S(numerator[i])*(S(point)+u)**i for i in range(numerator.degree()+1))
    den = (S(point)+u-lam)**denominator_power
    expansion = num/den
    return Kt([expansion[i] for i in range(precision)])


def multiplicative_steps_rational(A, B, N, M, lam, fiber_point, node):
    local_X = local_polynomial(N,2,lam,fiber_point)
    local_Y = local_polynomial(M,3,lam,fiber_point)
    local_A = Kt(A(t+fiber_point))
    local_B = Kt(B(t+fiber_point))
    return multiplicative_steps_polynomial(local_A,local_B,local_X,local_Y,K(0),node)


def rational_intersection(poly_x, poly_y, N, M, lam):
    q = t-lam
    return gcd(poly_x*q**2-N,poly_y*q**3-M).degree()


def find_p3(A, B, P1X, P1Y, P2X, P2Y):
    node1 = P1X(K(1))
    x1 = K(-A[1])/K(18)
    z = polygen(K, "z")
    vector_roots = [
        K(c) for c, multiplicity in (z**3+A[6]*z+B[9]).roots()
        if 3*K(c)**2+A[6] != 0
    ]
    for lam in K:
        if lam in (K(0),K(1)):
            continue
        q = t-lam
        for c in vector_roots:
            for n2 in K:
                for n3 in K:
                    n0 = 3*lam**2
                    n1 = lam**2*x1-2*n0/lam
                    n5 = c
                    n4 = node1*(1-lam)**2-n0-n1-n2-n3-n5
                    N = Kt([n0,n1,n2,n3,n4,n5])
                    if N(lam) == 0:
                        continue
                    H = N**3+A*N*q**4+B*q**6
                    for M in fast_square_roots(H,7):
                        if M.degree()!=7 or M(lam)==0 or M[0]!=0 or M[1]!=0 or M(1)!=0:
                            continue
                        if multiplicative_steps_rational(A,B,N,M,lam,K(0),K(3)) != 2:
                            continue
                        if multiplicative_steps_rational(A,B,N,M,lam,K(1),node1) != 1:
                            continue
                        p1dot = rational_intersection(P1X,P1Y,N,M,lam)
                        p2dot = rational_intersection(P2X,P2Y,N,M,lam)
                        if p1dot == 1 and p2dot == 1:
                            yield lam, N, M


def candidates_from_logs():
    for log_path in sorted(root.glob(args.glob)):
        meta_path = log_path.with_suffix("").with_suffix(".meta.txt")
        if not meta_path.exists():
            continue
        fixed = None
        for line in meta_path.read_text().splitlines():
            if line.startswith("values="):
                fixed = ast.literal_eval(line.split("=",1)[1])
        if fixed is None:
            continue
        for line in log_path.read_text().splitlines():
            if not line.startswith("MW3A10SCAN_HIT|"):
                continue
            hit = {}
            for item in line.split("|",1)[1].split(","):
                name, value = item.split("=",1)
                hit[name] = int(value)
            yield log_path.stem, {**fixed, **hit}


raw = valid = p2_count = p3_count = 0
seen_models = set()
for source, candidate in candidates_from_logs():
    raw += 1
    model = reconstruct(candidate)
    if model is None:
        continue
    A,B,P1X,P1Y = model
    key = (tuple(A.list()),tuple(B.list()),tuple(P1X.list()),tuple(P1Y.list()))
    if key in seen_models:
        continue
    seen_models.add(key)
    valid += 1
    tangent_suffix = ""
    if args.jacobian:
        tangent_rank = chart_jacobian_rank(candidate)
        tangent_suffix = (
            f"|jacobian_rank={tangent_rank}"
            f"|tangent_dimension={len(active_names)-tangent_rank}"
        )
    if args.jacobian or (args.progress_every and valid % args.progress_every == 0):
        print(f"MW3D6CHAIN|valid_p1={valid}|source={source}{tangent_suffix}",flush=True)
    for P2X,P2Y in find_p2(A,B,P1X,P1Y):
        p2_count += 1
        if args.show_p2:
            print(
                "MW3D6CHAIN_P2|A="+",".join(map(str,map(int,A.list())))
                +"|B="+",".join(map(str,map(int,B.list())))
                +"|P1X="+",".join(map(str,map(int,P1X.list())))
                +"|P1Y="+",".join(map(str,map(int,P1Y.list())))
                +"|P2X="+",".join(map(str,map(int,P2X.list())))
                +"|P2Y="+",".join(map(str,map(int,P2Y.list()))),flush=True,
            )
        for lam,N,M in find_p3(A,B,P1X,P1Y,P2X,P2Y):
            p3_count += 1
            print(
                "MW3D6CHAIN_TARGET|A="+",".join(map(str,map(int,A.list())))
                +"|B="+",".join(map(str,map(int,B.list())))
                +"|P1X="+",".join(map(str,map(int,P1X.list())))
                +"|P1Y="+",".join(map(str,map(int,P1Y.list())))
                +"|P2X="+",".join(map(str,map(int,P2X.list())))
                +"|P2Y="+",".join(map(str,map(int,P2Y.list())))
                +f"|lambda={int(lam)}|N="+",".join(map(str,map(int,N.list())))
                +"|M="+",".join(map(str,map(int,M.list()))),flush=True,
            )
    if args.max_surfaces and valid >= args.max_surfaces:
        break

print(
    f"MW3D6CHAINSUMMARY|raw={raw}|valid_p1={valid}|p2={p2_count}|p3={p3_count}",
    flush=True,
)
