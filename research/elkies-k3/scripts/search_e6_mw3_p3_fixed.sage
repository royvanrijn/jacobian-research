from sage.all import *
import argparse


ap = argparse.ArgumentParser(
    description="Exhaust the polynomial canonical P3 ansatz on the GF(31) E6/P1+P2 seed."
)
ap.add_argument("--max-hits", type=int, default=20)
args = ap.parse_args()

K = GF(31)
T = PolynomialRing(K, "t")
t = T.gen()

A = T([20, 23, 11, 10, 8, 8])
B = T([8, 20, 17, 0, 19, 27, 6, 28, 1])
s0 = K(18)
nodes = {K(1): K(5), K(10): K(16), K(23): K(14)}


def square_roots_by_coefficients(poly):
    if poly == 0:
        return [T(0)]
    degree = poly.degree()
    if degree % 2:
        return []
    root_degree = degree//2
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


# The section-system infinity conditions give deg(X)<=3 and deg(Y)<=5.
# A nonzero cubic term in X would make X^3 the unique degree-9 term, which
# cannot occur in a polynomial square.  Hence deg(X)<=2 exactly.
tested = 0
hits = []
for x1 in K:
    for x2 in K:
        tested += 1
        X3 = s0+x1*t+x2*t**2
        H = X3**3+A*X3+B
        for Y3 in square_roots_by_coefficients(H):
            if Y3(0) != 0:
                continue
            if X3.degree() != 2 or Y3.degree() != 4:
                continue
            if any(X3(point) == node and Y3(point) == 0 for point, node in nodes.items()):
                continue
            hits.append((x1, x2, X3, Y3))
            if len(hits) <= args.max_hits:
                print(
                    f"E6P3_HIT|x1={int(x1)}|x2={int(x2)}"
                    f"|X={','.join(str(int(c)) for c in X3.list())}"
                    f"|Y={','.join(str(int(c)) for c in Y3.list())}",
                    flush=True,
                )

print(f"E6P3SUMMARY|tested={tested}|hits={len(hits)}", flush=True)
