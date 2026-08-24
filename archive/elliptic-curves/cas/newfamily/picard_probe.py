#!/usr/bin/env sage -python
"""
Point-count / Frobenius Picard probe for the six-root K3.

For a clean prime p, count S(F_{p^n}), n=1,2,3.  The known 11 generic
sections plus the split I4 trivial lattice give 16 Q-defined divisor classes,
so only a reciprocal degree-6 H^2 factor is unknown.

Three traces reconstruct that degree-6 factor exactly.  Cyclotomic factors of
P_res(p*Z) give the reduction's geometric Picard contribution beyond 16.
If needed, an Artin-Tate square class is also computed after extending the
field enough to make all algebraic divisor classes rational.

This is intended as exact arithmetic research code; no numerical rank claim is
promoted.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

from sage.all import *

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))

from newfamily_rank11_minimal_common import build_finite_minimal_family

DEFAULT_ROOTS = (-47, -43, -31, 30, 45, 46)
KNOWN_MW_RANK = 11
TRIVIAL_RANK = 5
KNOWN_DIVISOR_RANK = 16

def parse_roots(text):
    roots = tuple(ZZ(x.strip()) for x in text.split(","))
    if len(roots) != 6:
        raise ValueError("need exactly six comma-separated roots")
    return roots

def reduce_poly(poly, K, p, name="t"):
    R = PolynomialRing(K, name)
    t = R.gen()
    ans = R(0)
    for i, c in enumerate(poly.list()):
        c = QQ(c)
        den = ZZ(c.denominator())
        if den % p == 0:
            raise ZeroDivisionError(f"coefficient denominator divisible by p={p}")
        ans += (K(ZZ(c.numerator())) / K(den)) * t**i
    return ans

def clean_model(roots, p):
    D = build_finite_minimal_family(roots)
    A, B, Delta = D["Amin"], D["Bmin"], D["Deltamin"]
    if (A.degree(), B.degree(), Delta.degree()) != (8, 12, 20):
        raise RuntimeError(
            f"expected degrees (8,12,20), got {(A.degree(),B.degree(),Delta.degree())}"
        )
    K = GF(p)
    Ap = reduce_poly(A, K, p)
    Bp = reduce_poly(B, K, p)
    Dp = reduce_poly(Delta, K, p)
    if (Ap.degree(), Bp.degree(), Dp.degree()) != (8, 12, 20):
        raise RuntimeError("leading coefficient vanished modulo p")
    if gcd(Dp, Dp.derivative()).degree() != 0:
        raise RuntimeError("finite discriminant is not squarefree modulo p")
    if gcd(Dp, Ap).degree() != 0:
        raise RuntimeError("finite discriminant meets c4=0 modulo p")

    RX = PolynomialRing(K, "X")
    X = RX.gen()
    g = X**3 + Ap[8]*X + Bp[12]
    fac = list(g.factor())
    doubles = [f for f,e in fac if e == 2 and f.degree() == 1]
    singles = [f for f,e in fac if e == 1 and f.degree() == 1]
    if len(doubles) != 1 or len(singles) != 1:
        raise RuntimeError(f"infinity is not the expected nodal cubic: {g.factor()}")
    def linroot(f):
        return -f[0] / f[1]
    r = linroot(doubles[0])
    s = linroot(singles[0])
    if not (r-s).is_square():
        raise RuntimeError("infinity I4 is nonsplit at this prime; skip it")
    return D

def field(p, n):
    if n == 1:
        return GF(p)
    return GF(p**n, name=f"z{n}", modulus="conway")

def frobenius_orbits(K, p, n):
    seen = set()
    for a in K:
        if a in seen:
            continue
        orb = []
        b = a
        while b not in orb:
            orb.append(b)
            b = b**p
        for x in orb:
            seen.add(x)
        yield a, len(orb)

def chi(v):
    if v == 0:
        return 0
    return 1 if v.is_square() else -1

def singular_character_sum(K, a, b):
    total = ZZ(0)
    for x in K:
        total += chi(x**3 + a*x + b)
    return total

def count_extension(D, p, n, progress_every):
    K = field(p, n)
    q = ZZ(K.cardinality())
    A = reduce_poly(D["Amin"], K, p)
    B = reduce_poly(D["Bmin"], K, p)

    total_character = ZZ(0)
    orbit_count = 0
    singular_orbits = 0
    started = time.time()

    for t, multiplicity in frobenius_orbits(K, p, n):
        a = A(t)
        b = B(t)
        disc_core = 4*a**3 + 27*b**2
        if disc_core == 0:
            c = singular_character_sum(K, a, b)
            singular_orbits += 1
        else:
            E = EllipticCurve(K, [0, 0, 0, a, b])
            c = ZZ(E.cardinality()) - q - 1
        total_character += multiplicity * c
        orbit_count += 1

        if progress_every and orbit_count % progress_every == 0:
            elapsed = time.time() - started
            print(
                f"COUNT|p={p}|n={n}|orbits={orbit_count}|"
                f"singular_orbits={singular_orbits}|seconds={elapsed:.3f}",
                flush=True,
            )

    # finite fibres: q(q+1)+sum chi; split I4 at infinity contributes 4q
    points = q*q + 5*q + total_character
    h2_trace = points - 1 - q*q

    return {
        "extension_degree": n,
        "field_size": int(q),
        "frobenius_orbits": orbit_count,
        "singular_orbits": singular_orbits,
        "finite_character_sum": int(total_character),
        "surface_points": int(points),
        "h2_trace": int(h2_trace),
        "seconds": time.time() - started,
    }

def residual_polynomial(p, counts):
    traces = {}
    for n in (1,2,3):
        q = ZZ(p)**n
        traces[n] = ZZ(counts[n]["h2_trace"]) - KNOWN_DIVISOR_RANK*q

    s1, s2, s3 = traces[1], traces[2], traces[3]
    c1 = -s1
    n2 = s1*s1 - s2
    if n2 % 2:
        raise ArithmeticError("Newton c2 is not integral")
    c2 = n2 // 2
    n3 = -(s1**3 - 3*s1*s2 + 2*s3)
    if n3 % 6:
        raise ArithmeticError("Newton c3 is not integral")
    c3 = n3 // 6

    R = PolynomialRing(QQ, "X")
    X = R.gen()
    P = (
        X**6 + c1*X**5 + c2*X**4 + c3*X**3
        + p**2*c2*X**2 + p**4*c1*X + p**6
    )

    # Replay first three power sums from Newton identities.
    coeff = [ZZ(P[6-k]) for k in range(7)]
    replay = power_sums_from_coefficients(coeff, 3)
    if tuple(replay[n] for n in (1,2,3)) != (s1,s2,s3):
        raise ArithmeticError("residual trace replay failed")

    return P, {str(n): int(traces[n]) for n in traces}

def power_sums_from_coefficients(c, max_n):
    # monic polynomial x^d + c1 x^(d-1)+...+cd
    d = len(c)-1
    if c[0] != 1:
        raise ValueError("polynomial must be monic")
    S = [ZZ(0)]*(max_n+1)
    for n in range(1, max_n+1):
        if n <= d:
            total = sum(c[k]*S[n-k] for k in range(1,n)) + n*c[n]
        else:
            total = sum(c[k]*S[n-k] for k in range(1,d+1))
        S[n] = -ZZ(total)
    return S

def powered_charpoly(P, exponent):
    d = P.degree()
    c = [ZZ(P[d-k]) for k in range(d+1)]
    S = power_sums_from_coefficients(c, d*exponent)
    T = [ZZ(0)] + [S[exponent*n] for n in range(1,d+1)]
    cp = [ZZ(1)]
    for n in range(1,d+1):
        numerator = T[n] + sum(cp[k]*T[n-k] for k in range(1,n))
        if numerator % n:
            raise ArithmeticError("powered Newton coefficient is not integral")
        cp.append(-numerator//n)
    R = P.parent()
    X = R.gen()
    return sum(cp[k]*X**(d-k) for k in range(d+1))

def cyclotomic_part(p, Pres):
    Rz = PolynomialRing(QQ, "Z")
    Z = Rz.gen()
    norm = Rz(Pres(p*Z) / QQ(p**6))
    rem = norm
    factors = []
    total_degree = 0
    orders = []
    for m in range(1, 31):
        if euler_phi(m) > 6:
            continue
        phi = Rz(cyclotomic_polynomial(m))
        mult = 0
        while True:
            q, r = rem.quo_rem(phi)
            if r != 0:
                break
            rem = q
            mult += 1
        if mult:
            degree = mult*phi.degree()
            total_degree += degree
            orders.extend([m]*mult)
            factors.append({
                "order": m,
                "multiplicity": mult,
                "degree": degree,
                "factor": str(phi),
            })
    return norm, rem, factors, total_degree, orders

def squareclass_rational(v):
    v = QQ(v)
    if v == 0:
        return 0
    sign = -1 if v < 0 else 1
    n = ZZ(abs(v.numerator()))
    d = ZZ(v.denominator())
    sf = ZZ(sign)
    for z in (n,d):
        for prime, exp in factor(z):
            if exp % 2:
                sf *= prime
    return int(sf)

def artin_tate_square_class(p, fullP, rho, orders):
    extension = 1
    for m in orders:
        extension = lcm(extension, m)
    q = ZZ(p)**extension
    Pq = powered_charpoly(fullP, extension)
    R = fullP.parent()
    X = R.gen()
    algebraic = (X-q)**rho
    quotient, remainder = Pq.quo_rem(algebraic)
    if remainder != 0:
        raise ArithmeticError("expected algebraic (X-q)^rho factor did not divide")
    limit = QQ(quotient(q))
    at = limit / QQ(q**(21-rho))
    signed = ((-1)**(rho-1))*at
    return {
        "extension_degree_defining_NS": int(extension),
        "q": int(q),
        "artin_tate_value_mod_Brauer_square": str(at),
        "signed_NS_discriminant_square_class": squareclass_rational(signed),
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roots", default=",".join(map(str, DEFAULT_ROOTS)))
    ap.add_argument("--prime", type=int, required=True)
    ap.add_argument("--progress-every", type=int, default=1000)
    ap.add_argument("--output")
    args = ap.parse_args()

    roots = parse_roots(args.roots)
    p = ZZ(args.prime)
    if not p.is_prime() or p <= 3:
        raise ValueError("--prime must be an odd prime > 3")

    D = clean_model(roots, int(p))
    print(
        f"MODEL|roots={','.join(map(str,roots))}|p={p}|"
        f"degA={D['Amin'].degree()}|degB={D['Bmin'].degree()}|"
        f"degDelta={D['Deltamin'].degree()}",
        flush=True,
    )

    counts = {}
    for n in (1,2,3):
        counts[n] = count_extension(D, int(p), n, args.progress_every)
        row = counts[n]
        print(
            f"POINTS|p={p}|n={n}|q={row['field_size']}|"
            f"N={row['surface_points']}|traceH2={row['h2_trace']}|"
            f"charsum={row['finite_character_sum']}|seconds={row['seconds']:.3f}",
            flush=True,
        )

    Pres, residual_traces = residual_polynomial(int(p), counts)
    norm, rem, cycfactors, cycdegree, orders = cyclotomic_part(int(p), Pres)
    rho_upper = int(KNOWN_DIVISOR_RANK + cycdegree)

    R = Pres.parent()
    X = R.gen()
    fullP = (X-p)**KNOWN_DIVISOR_RANK * Pres

    result = {
        "roots": [int(x) for x in roots],
        "prime": int(p),
        "known_generic_MW_rank": KNOWN_MW_RANK,
        "trivial_lattice_rank": TRIVIAL_RANK,
        "known_Q_divisor_rank": KNOWN_DIVISOR_RANK,
        "counts": {str(k): v for k,v in counts.items()},
        "residual_traces": residual_traces,
        "residual_degree6_characteristic_polynomial": str(Pres),
        "normalized_residual_polynomial": str(norm),
        "cyclotomic_factors": cycfactors,
        "noncyclotomic_remainder": str(rem),
        "residual_cyclotomic_degree": cycdegree,
        "reduction_geometric_Picard_rank": rho_upper,
        "reduction_geometric_MW_rank": rho_upper - TRIVIAL_RANK,
        "characteristic_zero_geometric_Picard_upper_bound_from_this_prime": rho_upper,
        "characteristic_zero_geometric_MW_upper_bound_from_this_prime": rho_upper - TRIVIAL_RANK,
    }

    # Tate is known for K3 surfaces in this setting, so the root-of-unity
    # multiplicity gives the reduction's geometric Picard rank.  Artin-Tate
    # then gives its NS discriminant square class.
    result["artin_tate"] = artin_tate_square_class(
        int(p), fullP, rho_upper, orders
    )

    print(
        f"FROB|p={p}|Pres={Pres}|cyclotomic_degree={cycdegree}|"
        f"rho_reduction={rho_upper}|mw_reduction={rho_upper-TRIVIAL_RANK}",
        flush=True,
    )
    print(
        f"ARTIN_TATE|p={p}|disc_square_class="
        f"{result['artin_tate']['signed_NS_discriminant_square_class']}|"
        f"extension={result['artin_tate']['extension_degree_defining_NS']}",
        flush=True,
    )

    if rho_upper == KNOWN_DIVISOR_RANK:
        print(
            "CONCLUSION|rho_Qbar=16|MW_Qbar(T)=11|status=exact_from_good_reduction",
            flush=True,
        )
    else:
        print(
            f"CONCLUSION|16<=rho_Qbar<={rho_upper}|"
            f"11<=MW_Qbar(T)<={rho_upper-TRIVIAL_RANK}|"
            "status=needs_more_prime_or_extra_section",
            flush=True,
        )

    if args.output:
        path = Path(args.output)
    else:
        path = (
            REPO / "artifacts/local/elliptic-curves/newfamily"
            / f"picard_probe_p{p}.json"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n")
    print(f"OUTPUT|{path}", flush=True)

if __name__ == "__main__":
    main()
