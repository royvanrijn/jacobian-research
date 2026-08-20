from sage.all import *
import argparse
import itertools
import numpy as np


ap = argparse.ArgumentParser(
    description="Exhaust the canonical P3 square-root slices on the verified GF(31) P1+P2 surface."
)
ap.add_argument("--max-slices", type=int, default=0, help="0 means all admissible slices")
ap.add_argument("--max-hits", type=int, default=20)
ap.add_argument("--A", default="4,23,18,28,12,18,23,20,19")
ap.add_argument("--B", default="23,24,27,10,6,8,23,26,26,16,19,9,15")
ap.add_argument("--lam", type=int, default=23)
ap.add_argument("--nodes", default="3,21,10")
ap.add_argument("--sinf", type=int, default=29)
args = ap.parse_args()

K = GF(31)
p = 31
Q = PolynomialRing(K, ("q0", "q1", "q2"), order="degrevlex")
q0, q1, q2 = Q.gens()
Qt = PolynomialRing(Q, "t")
t = Qt.gen()

A = Qt([int(x) for x in args.A.split(",")])
B = Qt([int(x) for x in args.B.split(",")])
lam = K(args.lam)
node_values = [K(int(x)) for x in args.nodes.split(",")]
if len(node_values) != 3:
    raise SystemExit("--nodes must give the nodes at 0,1,lambda")
s0, s1, sl = node_values
sinf = K(args.sinf)
bad = {K(0), K(1), lam}

grid = np.array(list(itertools.product(range(p), repeat=3)), dtype=np.int64).T


def evaluate(poly, points):
    values = np.zeros(points.shape[1], dtype=np.int64)
    cache = {}
    for exponents, coefficient in poly.dict().items():
        term = np.full(points.shape[1], int(coefficient), dtype=np.int64)
        for i, exponent in enumerate(exponents):
            if exponent == 0:
                continue
            key = (i, exponent)
            if key not in cache:
                cache[key] = np.power(points[i], exponent, dtype=np.int64) % p
            term = (term * cache[key]) % p
        values += term
        values %= p
    return values


hits = []
slices = 0
for pole in K:
    if pole in bad:
        continue

    C = Q(s0 * pole**2) * (t - Q(lam)) / Q(-lam)
    C += Q(sl * (lam - pole)**2) * t / Q(lam)
    G = t * (t - Q(lam))
    z = t - Q(pole)

    for q3_value in K:
        X = C + G * (
            q0 + q1 * t + q2 * t**2 + Q(q3_value) * t**3 + Q(sinf) * t**4
        )
        H = X**3 + A * X * z**4 + B * z**6
        if H[18] != 0 or H[17] != 0:
            raise RuntimeError("infinity incidence failed")
        h16 = K(H[16])
        leading_roots = h16.sqrt(all=True)
        if not leading_roots or h16 == 0:
            continue

        # +/-Y is the same section up to sign, so one leading square root is
        # enough for existence and independence discovery.
        y = [Q(0)] * 9
        y[8] = Q(leading_roots[0])
        for degree in range(15, 7, -1):
            index = degree - 8
            known = Q(0)
            for i in range(index + 1, 9):
                j = degree - i
                if 0 <= j <= 8 and j != index:
                    known += y[i] * y[j]
            y[index] = (Q(H[degree]) - known) / (2 * y[8])

        Y = sum(y[i] * t**i for i in range(9))
        residuals = [Q(Y(0)), Q(Y(lam))]
        for degree in range(8):
            square_coefficient = sum(
                y[i] * y[degree - i]
                for i in range(degree + 1)
                if i <= 8 and degree - i <= 8
            )
            residuals.append(Q(H[degree]) - square_coefficient)
        residuals = [residual for residual in residuals if residual != 0]

        slices += 1
        survivors = grid
        order = sorted(range(len(residuals)), key=lambda i: len(residuals[i].monomials()))
        for residual_index in order:
            values = evaluate(residuals[residual_index], survivors)
            survivors = survivors[:, values == 0]
            if survivors.shape[1] == 0:
                break

        if survivors.shape[1]:
            for column in range(survivors.shape[1]):
                assignment = {
                    q0: K(int(survivors[0, column])),
                    q1: K(int(survivors[1, column])),
                    q2: K(int(survivors[2, column])),
                }
                X_special = PolynomialRing(K, "t")(
                    [K(coefficient.subs(assignment)) for coefficient in X.list()]
                )
                Y_special = PolynomialRing(K, "t")(
                    [K(coefficient.subs(assignment)) for coefficient in Y.list()]
                )
                # Exact pole order and identity component at t=1 are open
                # conditions, not equations in the search ideal.
                if X_special(pole) == 0:
                    continue
                if X_special(1) == s1 * (K(1) - pole)**2 and Y_special(1) == 0:
                    continue
                hit = (
                    int(pole), int(q3_value),
                    int(survivors[0, column]), int(survivors[1, column]),
                    int(survivors[2, column]),
                    [int(c) for c in X_special.list()],
                    [int(c) for c in Y_special.list()],
                )
                hits.append(hit)
                if len(hits) <= args.max_hits:
                    print(
                        f"MW3A10P3_HIT|r={hit[0]}|q3={hit[1]}"
                        f"|q0={hit[2]}|q1={hit[3]}|q2={hit[4]}"
                        f"|X={','.join(map(str, hit[5]))}"
                        f"|Y={','.join(map(str, hit[6]))}",
                        flush=True,
                    )

        if slices % 50 == 0:
            print(f"MW3A10P3_SCAN|slices={slices}|hits={len(hits)}", flush=True)
        if args.max_slices and slices >= args.max_slices:
            break
    if args.max_slices and slices >= args.max_slices:
        break

print(f"MW3A10P3_SCAN|done=1|slices={slices}|hits={len(hits)}", flush=True)
