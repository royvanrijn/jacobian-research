from sage.all import *
import argparse


ap = argparse.ArgumentParser(
    description=(
        "Exhaust the polynomial P2 search on the certified GF(23) "
        "D6+A5+A3 surface."
    )
)
ap.add_argument("--p", type=int, default=23)
ap.add_argument("--all", action="store_true", help="print candidates of every intersection number")
args = ap.parse_args()

if args.p != 23:
    raise SystemExit("the pinned surface is currently defined only over GF(23)")

p = args.p
K = GF(p)
Kt = PolynomialRing(K, "t")
t = Kt.gen()

# Certified by verify_mw3_d6a5a3_p1_gf23.sage (seed-3 chart point).
A = Kt([19, 17, 4, 15, 12, 15, 10])
B = Kt([8, 18, 12, 9, 9, 4, 2, 18, 3, 9])
P1X = Kt([3, 8, 14, 18, 3])
P1Y = Kt([0, 0, 2, 9, 9, 5, 21])

if P1Y**2 != P1X**3 + A*P1X + B:
    raise RuntimeError("pinned P1 data are inconsistent")

# P2 has profile (vector,5,0) at (I2*,I6,I4) and P2.O=0.
# The vector class at infinity forces deg(X)=3 and deg(Y)=4.  At I6 it
# passes through the node (3,0); component 5 is the negative orientation of
# component 1, so both signs of Y are retained until the intersection gate.
node0 = K(3)
node1 = P1X(K(1))


def polynomial_square_root(H):
    if not H.is_square():
        return None
    factors = H.factor()
    unit = K(factors.unit())
    roots = unit.sqrt(all=True)
    if not roots:
        return None
    root = Kt(roots[0])
    for factor, exponent in factors:
        if exponent % 2:
            return None
        root *= factor ** (exponent // 2)
    return root


def multiplicative_steps(X, Y, fiber_point, node):
    """Number of consecutive singular centers followed in an I_n fiber."""
    P = PolynomialRing(K, ("u", "xx", "yy"))
    u, xx, yy = P.gens()
    shifted_A_t = Kt(A(t + fiber_point))
    shifted_B_t = Kt(B(t + fiber_point))
    shifted_A = sum(P(c) * u**i for i, c in enumerate(shifted_A_t.list()))
    shifted_B = sum(P(c) * u**i for i, c in enumerate(shifted_B_t.list()))
    surface = yy**2 - (node + xx)**3 - shifted_A*(node + xx) - shifted_B
    shifted_X = Kt(X(t + fiber_point))
    shifted_Y = Kt(Y(t + fiber_point))
    if shifted_X(0) != node or shifted_Y(0) != 0:
        return 0
    section_x = Kt((shifted_X - node) // t)
    section_y = Kt(shifted_Y // t)
    surface = P(surface(u, u*xx, u*yy) // u**2)
    steps = 1
    while True:
        center_x = K(section_x(0))
        center_y = K(section_y(0))
        point = {u: K(0), xx: center_x, yy: center_y}
        gradient = [surface.derivative(v).subs(point) for v in (u, xx, yy)]
        if any(gradient):
            return steps
        section_x = Kt((section_x - center_x) // t)
        section_y = Kt((section_y - center_y) // t)
        transformed = surface(u, center_x + u*xx, center_y + u*yy)
        surface = P(transformed // u**2)
        steps += 1
        if steps > 6:
            raise RuntimeError("multiplicative component resolution did not terminate")


def i2star_vector_check(X, Y):
    """Verify that the section follows one cusp center, then exits smoothly."""
    u = Kt.gen()
    # Minimal coordinates at infinity: xbar=u^4*x(1/u), ybar=u^6*y(1/u).
    xbar = Kt(sum(X[i] * u**(4-i) for i in range(X.degree()+1)))
    ybar = Kt(sum(Y[i] * u**(6-i) for i in range(Y.degree()+1)))
    Abar = Kt(sum(A[i] * u**(8-i) for i in range(A.degree()+1)))
    Bbar = Kt(sum(B[i] * u**(12-i) for i in range(B.degree()+1)))
    if xbar.valuation() != 1 or ybar.valuation() != 2:
        return False, "wrong-valuations"
    c = K((xbar // u)(0))
    d = K((ybar // u**2)(0))
    cubic = c**3 + Abar[2]*c + Bbar[3]
    if cubic:
        return False, "misses-first-singular-center"
    # After the second blowup, x=c*u+..., y=d*u^2+....  A simple root of
    # the exceptional cubic is precisely the smooth outer vector component.
    derivative = 3*c**2 + Abar[2]
    if derivative == 0:
        return False, "non-vector-deeper-center"
    return True, f"c={int(c)},d={int(d)}"


tested = 0
degree8 = 0
squares = 0
profile_hits = 0
target_hits = 0
for q1 in K:
    for q2 in K:
        for q3 in K:
            tested += 1
            X = Kt([node0, q1, q2, q3])
            H = X**3 + A*X + B
            if H.degree() != 8:
                continue
            degree8 += 1
            Y0 = polynomial_square_root(H)
            if Y0 is None or Y0.degree() != 4 or Y0(0) != 0:
                continue
            squares += 1
            for Y in (Y0, -Y0):
                if X(K(1)) == node1 and Y(K(1)) == 0:
                    continue
                i6steps = multiplicative_steps(X, Y, K(0), node0)
                vector_ok, vector_detail = i2star_vector_check(X, Y)
                if i6steps != 1 or not vector_ok:
                    continue
                profile_hits += 1
                same = gcd(P1X-X, P1Y-Y).degree()
                opposite = gcd(P1X-X, P1Y+Y).degree()
                if same == 1 or opposite == 1:
                    target_hits += 1
                if args.all or same == 1 or opposite == 1:
                    print(
                        "MW3D6P2|X=" + ",".join(map(str, map(int, X.list())))
                        + "|Y=" + ",".join(map(str, map(int, Y.list())))
                        + f"|I6steps={i6steps}|I2star={vector_detail}"
                        + f"|P1dotP2={same}|P1dotNegP2={opposite}",
                        flush=True,
                    )

print(
    f"MW3D6P2SUMMARY|p={p}|tested={tested}|degree8={degree8}"
    f"|squares={squares}|profile_hits={profile_hits}|target_hits={target_hits}",
    flush=True,
)
