#!/usr/bin/env sage -python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from sage.all import *

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from newfamily_rank11_minimal_common import build_finite_minimal_family

DEFAULT_ROOTS = (-47, -43, -31, 30, 45, 46)

def parse_roots(text):
    roots = tuple(ZZ(x.strip()) for x in text.split(","))
    if len(roots) != 6:
        raise ValueError("need exactly six comma-separated roots")
    return roots

def reduce_poly(poly, K, p):
    R = PolynomialRing(K, "t")
    t = R.gen()
    ans = R(0)
    for i, c in enumerate(poly.list()):
        c = QQ(c)
        den = ZZ(c.denominator())
        if den % p == 0:
            raise ZeroDivisionError(f"coefficient denominator divisible by p={p}")
        cc = K(ZZ(c.numerator())) / K(den)
        ans += cc * t**i
    return ans

def infinity_split(Ap, Bp):
    K = Ap.base_ring()
    R = PolynomialRing(K, "X")
    X = R.gen()
    a8 = Ap[8] if Ap.degree() >= 8 else K(0)
    b12 = Bp[12] if Bp.degree() >= 12 else K(0)
    g = X**3 + a8*X + b12
    fac = list(g.factor())
    doubles = [f for f,e in fac if e == 2 and f.degree() == 1]
    singles = [f for f,e in fac if e == 1 and f.degree() == 1]
    if len(doubles) != 1 or len(singles) != 1:
        return False, str(g.factor())
    def root_of_linear(f):
        return -f[0] / f[1]
    r = root_of_linear(doubles[0])
    s = root_of_linear(singles[0])
    gap = r - s
    return bool(gap.is_square()), str(g.factor())

def good_prime(roots, p):
    if p <= 3:
        return None
    D = build_finite_minimal_family(roots)
    A, B, Delta = D["Amin"], D["Bmin"], D["Deltamin"]
    K = GF(p)
    try:
        Ap = reduce_poly(A, K, p)
        Bp = reduce_poly(B, K, p)
        Dp = reduce_poly(Delta, K, p)
    except ZeroDivisionError:
        return None
    if Ap.degree() != 8 or Bp.degree() != 12 or Dp.degree() != 20:
        return None
    if gcd(Dp, Dp.derivative()).degree() != 0:
        return None
    if gcd(Dp, Ap).degree() != 0:
        return None
    split, inf = infinity_split(Ap, Bp)
    if not split:
        return None
    return inf

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roots", default=",".join(map(str, DEFAULT_ROOTS)))
    ap.add_argument("--max-prime", type=int, default=250)
    args = ap.parse_args()

    roots = parse_roots(args.roots)
    print("ROOTS", ",".join(map(str, roots)))
    print("GOOD_SPLIT_PRIMES")
    good = []
    for p in prime_range(5, args.max_prime + 1):
        inf = good_prime(roots, int(p))
        if inf is not None:
            good.append(int(p))
            print(f"{int(p):4d}  infinity={inf}")
    print()
    print("PRIMES=" + ",".join(map(str, good)))
    print("COUNT=" + str(len(good)))

if __name__ == "__main__":
    main()
