from sage.all import *
import argparse


ap = argparse.ArgumentParser(
    description=(
        "Search the one-pole P3 on the certified GF(23) D6+A5+A3 surface "
        "with the exact MW3 component and intersection profile."
    )
)
ap.add_argument("--p", type=int, default=23)
ap.add_argument("--all", action="store_true")
args = ap.parse_args()

if args.p != 23:
    raise SystemExit("the pinned surface is currently defined only over GF(23)")

p = args.p
K = GF(p)
Kt = PolynomialRing(K, "t")
t = Kt.gen()
A = Kt([19, 17, 4, 15, 12, 15, 10])
B = Kt([8, 18, 12, 9, 9, 4, 2, 18, 3, 9])
P1X = Kt([3, 8, 14, 18, 3])
P1Y = Kt([0, 0, 2, 9, 9, 5, 21])
P2X = Kt([3, 18, 10, 5])
P2Y = Kt([0, 7, 12, 5, 13])
Delta = -16*(4*A**3 + 27*B**2)
node0 = K(3)
node1 = P1X(K(1))
x1_component2 = K(-17)/K(18)


def fast_square_roots_degree14(H):
    """Recover square roots top-down without factoring H."""
    if H.degree() != 14:
        return []
    leading_roots = K(H[14]).sqrt(all=True)
    roots = []
    for leading in leading_roots:
        coefficients = [K(0)] * 8
        coefficients[7] = leading
        for degree in range(13, 6, -1):
            unknown_index = degree - 7
            known = K(0)
            for i in range(8):
                j = degree - i
                if 0 <= j < 8 and i != 7 and j != 7:
                    known += coefficients[i]*coefficients[j]
            coefficients[unknown_index] = (H[degree] - known)/(2*leading)
        root = Kt(coefficients)
        if root**2 == H:
            roots.append(root)
    return roots


def local_polynomial(numerator, denominator_power, point, precision=8):
    S = PowerSeriesRing(K, "u", default_prec=precision)
    u = S.gen()
    num = sum(S(numerator[i])*(S(point)+u)**i for i in range(numerator.degree()+1))
    den = (S(point)+u-lam)**denominator_power
    expansion = num/den
    return Kt([expansion[i] for i in range(precision)])


def multiplicative_steps_rational(N, M, qpower_x, qpower_y, fiber_point, node):
    """Resolve the component followed by a rational section at a finite I_n."""
    local_X = local_polynomial(N, qpower_x, fiber_point)
    local_Y = local_polynomial(M, qpower_y, fiber_point)
    if local_X(0) != node or local_Y(0) != 0:
        return 0
    P = PolynomialRing(K, ("u", "xx", "yy"))
    u, xx, yy = P.gens()
    shifted_A_t = Kt(A(t + fiber_point))
    shifted_B_t = Kt(B(t + fiber_point))
    shifted_A = sum(P(c)*u**i for i, c in enumerate(shifted_A_t.list()))
    shifted_B = sum(P(c)*u**i for i, c in enumerate(shifted_B_t.list()))
    surface = yy**2 - (node + xx)**3 - shifted_A*(node + xx) - shifted_B
    section_x = Kt((local_X-node)//t)
    section_y = Kt(local_Y//t)
    surface = P(surface(u, u*xx, u*yy)//u**2)
    steps = 1
    while True:
        center_x = K(section_x(0))
        center_y = K(section_y(0))
        point = {u: K(0), xx: center_x, yy: center_y}
        gradient = [surface.derivative(v).subs(point) for v in (u, xx, yy)]
        if any(gradient):
            return steps
        section_x = Kt((section_x-center_x)//t)
        section_y = Kt((section_y-center_y)//t)
        transformed = surface(u, center_x+u*xx, center_y+u*yy)
        surface = P(transformed//u**2)
        steps += 1
        if steps > 6:
            raise RuntimeError("component resolution did not terminate")


def intersection_degree(poly_x, poly_y, N, M, sign=1):
    q = t-lam
    xeq = poly_x*q**2-N
    yeq = poly_y*q**3-sign*M
    return gcd(xeq, yeq).degree()


# The vector component at infinity makes the leading X coefficient a simple
# root of c^3+A_6*c+B_9.  Parameterizing by those roots removes one full scan
# coordinate.  I6 component 2 fixes N_0,N_1; I4 component 1 fixes N(1).
z = polygen(K, "z")
vector_roots = [K(root) for root, multiplicity in (z**3+A[6]*z+B[9]).roots()]
tested = 0
square_hits = 0
profile_hits = 0
intersection_hits = 0

for lam in K:
    if lam in (K(0), K(1)):
        continue
    q = t-lam
    for c in vector_roots:
        for n2 in K:
            for n3 in K:
                tested += 1
                n0 = node0*lam**2
                n1 = lam**2*x1_component2 - 2*n0/lam
                n5 = c
                n4 = node1*(1-lam)**2 - n0-n1-n2-n3-n5
                N = Kt([n0, n1, n2, n3, n4, n5])
                if N(lam) == 0:
                    continue
                H = N**3 + A*N*q**4 + B*q**6
                roots = fast_square_roots_degree14(H)
                if not roots:
                    continue
                square_hits += 1
                for M in roots:
                    if M.degree() != 7 or M(lam) == 0:
                        continue
                    # y1=0 at I6 and node incidence at I4.
                    if M[0] != 0 or M[1] != 0 or M(K(1)) != 0:
                        continue
                    i6steps = multiplicative_steps_rational(N, M, 2, 3, K(0), node0)
                    i4steps = multiplicative_steps_rational(N, M, 2, 3, K(1), node1)
                    if i6steps != 2 or i4steps != 1:
                        continue
                    profile_hits += 1
                    p1dot = intersection_degree(P1X, P1Y, N, M)
                    p2dot = intersection_degree(P2X, P2Y, N, M)
                    if p1dot == 1 and p2dot == 1:
                        intersection_hits += 1
                    if args.all or (p1dot == 1 and p2dot == 1):
                        print(
                            f"MW3D6P3|lambda={int(lam)}"
                            + "|N=" + ",".join(map(str, map(int, N.list())))
                            + "|M=" + ",".join(map(str, map(int, M.list())))
                            + f"|I6steps={i6steps}|I4steps={i4steps}"
                            + f"|P1dotP3={p1dot}|P2dotP3={p2dot}",
                            flush=True,
                        )

print(
    f"MW3D6P3SUMMARY|p={p}|vector_roots={','.join(map(str,map(int,vector_roots)))}"
    f"|tested={tested}|square_hits={square_hits}|profile_hits={profile_hits}"
    f"|intersection_hits={intersection_hits}",
    flush=True,
)
