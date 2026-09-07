from sage.all import GF, PolynomialRing, is_prime
from itertools import product
import argparse


ap = argparse.ArgumentParser(
    description="Exhaust the normalized all-I3 family and its canonical P1 x-ansatz over GF(p)."
)
ap.add_argument("--p", type=int, default=5)
ap.add_argument("--max-hits", type=int, default=0,
                help="stop after this many signed-square-free x hits; zero is exhaustive")
ap.add_argument("--trivial-root", choices=("all", "0", "1", "lam"), default="all",
                help="A2 factor on which P1 has label zero")
ap.add_argument("--stage", choices=("p1", "p2", "triples"), default="p1")
args = ap.parse_args()
if not is_prime(args.p) or args.p in (2, 3):
    raise SystemExit("choose a prime other than 2 or 3")

F = GF(args.p)
Rt = PolynomialRing(F, "t")
t = Rt.gen()
values = list(F)
inv12 = F(12) ** -1
q0_representatives = (F(1), F.multiplicative_generator())

surfaces = 0
x_tests = 0
hits = []
p1_objects = []

for q0, q1, q2, a0, a1 in product(q0_representatives, values, values, values, values):
    # Weierstrass scaling changes q0 by a square, so one representative of
    # each nonzero square class is required over the finite field.
    q_at_1 = q0 + q1 + q2
    if q_at_1 == 0:
        continue

    lam_coefficient = q2 * (
        -q1 * a0**2
        - q2 * a0**2
        + 2 * q0 * a0 * a1
        + q0 * a1**2
    )
    if lam_coefficient == 0:
        continue
    constant = (
        -q1**2 * a0**2
        + q0 * q2 * a0**2
        - q1 * q2 * a0**2
        + 2 * q0 * q1 * a0 * a1
        + 2 * q0 * q2 * a0 * a1
        - q0**2 * a1**2
    )
    lam = -constant / lam_coefficient
    if lam in (0, 1):
        continue

    q = q0 + q1 * t + q2 * t**2
    aa = a0 + a1 * t
    if q(lam) == 0:
        continue
    b0 = a0**2 * inv12 / q0
    b1 = (a0 + a1) ** 2 * inv12 / q_at_1 - b0
    bb = b0 + b1 * t
    P = t * (t - 1) * (t - lam)
    relation = 12 * q * bb - aa**2
    quotient, remainder = relation.quo_rem(P)
    if remainder != 0 or quotient.degree() > 0:
        raise RuntimeError("all-I3 coefficient chart identity failed")
    c = quotient[0] if quotient else F(0)

    A = -3 * q**2 + P * aa
    B = 2 * q**3 - P * q * aa + P**2 * bb
    residual = 9 * q**2 * c + 4 * aa**3 - 54 * q * aa * bb + 27 * P * bb**2
    if A.degree() != 4 or B.degree() != 7 or residual.degree() != 5:
        continue
    if not residual.is_squarefree() or residual.gcd(P) != 1:
        continue
    surfaces += 1

    roots = (("0", F(0)), ("1", F(1)), ("lam", lam))
    profiles = roots if args.trivial_root == "all" else tuple(
        item for item in roots if item[0] == args.trivial_root
    )
    for trivial_name, trivial_root in profiles:
        nontrivial_roots = [root for name, root in roots if name != trivial_name]
        node_factor = (t - nontrivial_roots[0]) * (t - nontrivial_roots[1])
        for u0, u1, u2 in product(values, repeat=3):
            U = u0 + u1 * t + u2 * t**2
            X = q + node_factor * U
            H = X**3 + A * X + B
            x_tests += 1
            square, Y = H.is_square(root=True)
            if not square or Y.degree() > 6:
                continue
            if any(Y(root) != 0 for root in nontrivial_roots):
                continue
            if X(trivial_root) == q(trivial_root) and Y(trivial_root) == 0:
                continue
            record = (
                tuple(int(x) for x in (q0, q1, q2, a0, a1, lam, b0, b1, c)),
                trivial_name,
                tuple(int(x) for x in (u0, u1, u2)),
                tuple(int(x) for x in X.list()),
                tuple(int(x) for x in Y.list()),
            )
            hits.append(record)
            p1_objects.append({
                "family": record[0],
                "q": q,
                "A": A,
                "B": B,
                "roots": roots,
                "trivial_name": trivial_name,
                "trivial_root": trivial_root,
                "nontrivial_roots": tuple(nontrivial_roots),
                "U1": U,
                "X1": X,
                "Y1": Y,
            })
            print(
                "E8A2P1|hit={}|family={}|trivial={}|U={}|X={}|Y={}".format(
                    len(hits), record[0], record[1], record[2], record[3], record[4]
                ),
                flush=True,
            )
            if args.max_hits and len(hits) >= args.max_hits:
                print(
                    "E8A2P1|status=EARLY_STOP|p={}|surfaces={}|x_tests={}|hits={}".format(
                        args.p, surfaces, x_tests, len(hits)
                    ),
                    flush=True,
                )
                raise SystemExit(0)

print(
    "E8A2P1|status=EXHAUSTIVE|p={}|surfaces={}|x_tests={}|hits={}".format(
        args.p, surfaces, x_tests, len(hits)
    ),
    flush=True,
)

if args.stage == "p1":
    raise SystemExit(0)

# P2 has one nonzero A2 label, on either of the two factors where P1 is
# nonzero.  Its label is opposite P1's label, so choose the sign of Y2 giving
# the opposite tangent branch at the common I3 node.
p2_tests = 0
p2_objects = []
for p1 in p1_objects:
    q = p1["q"]
    A = p1["A"]
    B = p1["B"]
    X1 = p1["X1"]
    Y1 = p1["Y1"]
    roots = p1["roots"]
    for shared_root in p1["nontrivial_roots"]:
        for u0, u1, u2, u3 in product(values, repeat=4):
            U2 = u0 + u1*t + u2*t**2 + u3*t**3
            X2 = q + (t-shared_root)*U2
            H2 = X2**3 + A*X2 + B
            p2_tests += 1
            square, Y2base = H2.is_square(root=True)
            if not square or Y2base.degree() > 6 or Y2base(shared_root) != 0:
                continue
            if any(
                X2(root) == q(root) and Y2base(root) == 0
                for _, root in roots if root != shared_root
            ):
                continue

            dx1 = (X1-q).derivative()(shared_root)
            dy1 = Y1.derivative()(shared_root)
            dx2 = (X2-q).derivative()(shared_root)
            dy2 = Y2base.derivative()(shared_root)
            if (dx1,dy1) == (0,0) or (dx2,dy2) == (0,0):
                continue
            signs = [sign for sign in (F(1),F(-1))
                     if dy1*dx2 != sign*dy2*dx1]
            if len(signs) != 1:
                continue
            Y2 = signs[0]*Y2base
            record = {
                "p1": p1,
                "shared_root": shared_root,
                "U2": U2,
                "X2": X2,
                "Y2": Y2,
            }
            p2_objects.append(record)
            print(
                "E8A2P2|hit={}|family={}|P1trivial={}|shared={}|U={}|X={}|Y={}".format(
                    len(p2_objects),
                    p1["family"],
                    p1["trivial_name"],
                    int(shared_root),
                    tuple(int(x) for x in (u0,u1,u2,u3)),
                    tuple(int(x) for x in X2.list()),
                    tuple(int(x) for x in Y2.list()),
                ),
                flush=True,
            )

print(
    "E8A2P2|status=EXHAUSTIVE|p={}|p1_hits={}|tests={}|hits={}".format(
        args.p, len(p1_objects), p2_tests, len(p2_objects)
    ),
    flush=True,
)

if args.stage == "p2":
    raise SystemExit(0)

# P3 has trivial component at all three I3 fibers. Search it once per surface
# surviving P1+P2, then pair the bounded height-4 hits with every P1+P2 record.
surface_groups = {}
for p2 in p2_objects:
    surface_groups.setdefault(p2["p1"]["family"], []).append(p2)

p3_tests = 0
triple_count = 0
for family, p2_group in surface_groups.items():
    p1 = p2_group[0]["p1"]
    q, A, B, roots = p1["q"], p1["A"], p1["B"], p1["roots"]
    p3_hits = []
    for x0,x1,x2,x3,x4 in product(values, repeat=5):
        X3 = x0+x1*t+x2*t**2+x3*t**3+x4*t**4
        H3 = X3**3+A*X3+B
        p3_tests += 1
        square,Y3 = H3.is_square(root=True)
        if not square or Y3.degree()>6:
            continue
        if any(X3(root)==q(root) and Y3(root)==0 for _,root in roots):
            continue
        p3_hits.append((X3,Y3))
    print(
        "E8A2P3|family={}|p2_hits={}|p3_hits={}".format(
            family,len(p2_group),len(p3_hits)
        ),
        flush=True,
    )
    for p2 in p2_group:
        for X3,Y3 in p3_hits:
            triple_count += 1
            print(
                "E8A2TRIPLE|hit={}|family={}|P1trivial={}|P2shared={}|X3={}|Y3={}".format(
                    triple_count,
                    family,
                    p2["p1"]["trivial_name"],
                    int(p2["shared_root"]),
                    tuple(int(x) for x in X3.list()),
                    tuple(int(x) for x in Y3.list()),
                ),
                flush=True,
            )

print(
    "E8A2TRIPLE|status=EXHAUSTIVE|p={}|families={}|p2_hits={}|p3_tests={}|triples={}".format(
        args.p,len(surface_groups),len(p2_objects),p3_tests,triple_count
    ),
    flush=True,
)
