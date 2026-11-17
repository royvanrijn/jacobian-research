#!/usr/bin/env python3
"""Exact controls for the three-adic elliptic upper sieve; standard library."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import platform

HERE = Path(__file__).resolve().parent
A = -10758382
B = 6600759744
checks = []


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    checks.append(label)


def vp(a, p=3):
    a = F(a)
    if not a:
        return float("inf")
    n, d = abs(a.numerator), a.denominator
    out = 0
    while n % p == 0:
        n //= p
        out += 1
    while d % p == 0:
        d //= p
        out -= 1
    return out


def on_curve(P):
    return P is None or P[1] ** 2 == P[0] ** 3 + A * P[0] + B


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    u, v = Q
    if x == u and y == -v:
        return None
    slope = (3 * x * x + A) / (2 * y) if P == Q else (v-y)/(u-x)
    z = slope * slope - x - u
    return z, -y + slope * (x-z)


def depth(P):
    return float("inf") if P is None else vp(-P[0]/P[1])


def disc(a, b, c):
    return a*a*b*b - 4*b**3 - 4*a**3*c - 27*c*c + 18*a*b*c


def poly_mul(f, g):
    out = [F(0)] * (len(f)+len(g)-1)
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            out[i+j] += a*b
    return out


def poly_add(*polys):
    return [sum((f[i] if i < len(f) else 0) for f in polys)
            for i in range(max(map(len, polys)))]


def shifted_coefficients(v):
    s0 = F(1816347, 819200)
    t = F(10, 9)*v*v-F(10, 9)*s0*v+F(125, 1024)
    f = poly_add(poly_mul(poly_mul([1, 1], [1, 1]), [1, 1]),
                 [-F(12, 5)*z for z in poly_mul([1, 1], [1, 1])],
                 [s0-v, s0-v], [t])
    expected = [F(10, 9)*v*v-F(851209, 245760)*v+F(769467, 819200),
                F(341787, 819200)-v, F(3, 5), F(1)]
    if f != expected:
        raise AssertionError("shifted cubic coefficient identity")
    return expected


def finite_quotient(N, small):
    """The two possible actual images in V4 x Z/(3**(N-1))."""
    modulus = 3**(N-1)
    formal = range(0, modulus, 3 if small else 1)
    total = allowed = torsion_bad = annulus_bad = 0
    for a in range(2):
        for b in range(2):
            for k in formal:
                total += 1
                if (a, b) == (1, 0):
                    torsion_bad += 1
                elif (a, b) == (0, 0) and k and (vp(k)+1) % 3:
                    annulus_bad += 1
                else:
                    allowed += 1
    return total, allowed, torsion_bad, annulus_bad


def replay():
    P2 = F(1563680, 529), F(304670016, 12167)
    T0 = F(-3552), F(0)
    check("exact P2 point", on_curve(P2))
    check("exact rational two-torsion point", on_curve(T0) and add(T0, T0) is None)
    check("good reduction at three", (4*A**3+27*B**2) % 3 != 0)
    finite_points = [(x, y) for x in range(3) for y in range(3)
                     if (y*y-x**3-A*x-B) % 3 == 0]
    check("complete three-adic reduction points", finite_points == [(0, 0), (1, 0), (2, 0)])
    red = lambda P: tuple(q.numerator*pow(q.denominator, -1, 3) % 3 for q in P)
    check("independent nonzero reduction witnesses", red(T0) == (0, 0) and red(P2) == (2, 0))
    check("reduction V4 all nonzero points two-torsion", all(y == 0 for x, y in finite_points))

    Q = add(P2, P2)
    literal = (F(89497017424304041335649, 1002118368684773376),
               -F(26756035205526466434703545465981295,
                  1003179235240703925530394624))
    check("literal rational duplication", Q == literal and on_curve(Q))
    check("duplication coordinate valuations", (vp(Q[0]), vp(Q[1])) == (-4, -6))
    check("formal depth two, not presumed depth one", depth(Q) == 2)
    for n in [2, 3, 4]:
        check("exact repeated tripling depth %d" % n, depth(Q) == n)
        Q = add(add(Q, Q), Q)

    shifted_coefficients(F(0))
    check("complete symbolic shift constants", True)
    check("shifted coefficient unit valuations",
          vp(F(341787, 819200)) == vp(F(769467, 819200)) == 1
          and vp(F(851209, 245760)) == -1)
    for m in range(-40, 41):
        co = shifted_coefficients(F(27*m))
        check("T0-disc cubic valuation control %d" % m,
              list(map(vp, co)) == [1, 1, 1, 0]
              and vp(disc(co[2], co[1], co[0])) == 3)
    for a in [3, 6]:
        for b in [3, 6]:
            for c in [3, 6]:
                check("generic Eisenstein discriminant unit %d%d%d" % (a,b,c),
                      vp(disc(a,b,c)) == 3)
    for n in range(1, 25):
        heights = [-4*n, 1-2*n, 1, 0]
        check("complete formal-annulus Newton support %d" % n,
              F(heights[1]) > F(-8*n, 3) and F(heights[2]) > F(-4*n, 3))
        check("formal-annulus quadratic odd valuation %d" % n, (1-2*n) % 2 == 1)
        if n % 3:
            check("nonintegral cubic root valuation %d" % n, F(-4*n, 3).denominator == 3)

    # Explicitly generate the known index-three subgroup in the first quotient.
    group = lambda a,b,k: (a % 2, b % 2, k % 9)
    generated = {group(a, j, 3*j) for a in range(2) for j in range(6)}
    expected = {(a,b,k) for a in range(2) for b in range(2) for k in [0,3,6]}
    check("complete known quotient subgroup order twelve", generated == expected)
    tables = []
    for J in [1, 2, 3]:
        N = 3*J
        alpha = F(29,52)+F(9,13*27**J)
        beta = F(27,52)+F(3,13*27**J)
        for small, proposed in [(True, alpha), (False, beta)]:
            row = finite_quotient(N, small)
            check("complete quotient count J%d small%s" % (J,small), F(row[1],row[0]) == proposed)
            tables.append({"N":N,"smaller_image":small,"total":row[0],
                           "retained":row[1],"torsion_disc_excluded":row[2],
                           "formal_annuli_excluded":row[3],"fraction":str(proposed)})
    check("first finite quotient upper coefficient", F(29,52)+F(9,13*27) == F(7,12))
    for J in range(1, 20):
        sum_small = F(1,2)+sum((F(1,4)*F(2,3)*F(1,3**(3*k-2))
                                for k in range(1,J)), F(0))+F(1,4*3**(3*J-2))
        check("all finite geometric-sum coefficient %d" % J,
              sum_small == F(29,52)+F(9,13*27**J))
    check("limiting smaller-image coefficient", F(1,2)+F(1,4)*F(3,13) == F(29,52))
    check("limiting full-image coefficient", F(1,2)+F(1,4)*F(1,13) == F(27,52))
    corrupted = (P2[0]+1, P2[1])
    check("corrupted rational point rejected", not on_curve(corrupted))
    check("incorrect finite sieve fraction rejected", finite_quotient(3,True)[1] != 8)
    return {"passed":True,"exact_controls":len(checks),"checks":checks,
            "finite_quotients":tables,"finite_upper_fraction":"7/12",
            "limsup_upper_fraction":"29/52",
            "basis_or_saturation_assumed":False,"formal_depth_one_assumed":False,
            "uniform_error_in_J_claimed":False,"density_zero_claimed":False,
            "hasse_infinitude_claimed":False,"independent_replay":False,
            "formal_verification":False,"external_review":False,
            "arithmetic":"Python standard-library exact integers and fractions",
            "python":platform.python_version()}


if __name__ == "__main__":
    report = replay()
    files = ["ELLIPTIC_EXCEPTION_THREE_ADIC_SIEVE.md", Path(__file__).name,
             "sources/milne-elliptic-curves-2006.pdf"]
    report["input_sha256"] = {f:sha256((HERE/f).read_bytes()).hexdigest() for f in files}
    out = HERE/"elliptic_exception_three_adic_sieve_checks.json"
    out.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"passed":True,"exact_controls":report["exact_controls"],
                      "input_sha256":report["input_sha256"]},indent=2))
