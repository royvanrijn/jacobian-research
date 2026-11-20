#!/usr/bin/env python3
"""Replay the contour constants, forest determinant, and ordinary gap-six prefix."""
from argparse import ArgumentParser
from fractions import Fraction as F
from functools import cache
from hashlib import sha256
from itertools import product
from math import comb, factorial, isqrt, prod
from pathlib import Path
import json

from flint import fmpz_mat, fmpz_poly, nmod_poly

def exp_bounds(x):
    x = F(x)
    assert abs(x) < 122
    # For x>=0, bound the tail after the degree-120 Taylor sum by a geometric series.
    if x<0:
        lo,hi=exp_bounds(-x)
        return 1/hi,1/lo
    lo=sum(x**j/factorial(j) for j in range(121))
    hi=lo+(x**121/factorial(121))/(1-x/122)
    return lo,hi
def check_contour_constants():
    jmin,jmax=F(34,5),F(36,5)
    kmin,kmax=F(99,100),F(101,100)
    dj,dk=F(1,5),F(1,100)
    rj,rk=F(3,4),F(1,5)
    cj,ck=F(7,200),F(1,4)
    bj,bk=F(1,30),F(5,6)
    aj,ak=F(3,25),F(9,10)
    amaxj,amaxk=F(1,5),F(11,10)
    xmaxj,xminj=F(9),F(28,5)
    assert dj*jmin**6/(6+dj)>3180
    assert kmax**6*(6+dk)/1000<dk
    jvar=dj*(6*jmax-7)/(7*jmin**2)
    kvar=7*dk/kmin**2
    assert F(1,7)-jvar>aj and F(1,7)+jvar<amaxj
    assert 1-kvar>ak and 1+kvar<amaxk
    def log_tail(r):
        return r*r/(2*(1-r))
    assert dj/6+7*log_tail(dj/7)+log_tail(dj/6)<F(1,25)
    assert exp_bounds(F(681,350))[1]<7
    assert log_tail(dk)+ (kmax**6/1000)/(1-kmax**6/1000)<F(1,500)
    def j_env(x):
        return x*(3200+x**6)*exp_bounds(-x)[1]
    assert max(j_env(F(3)),j_env(xminj))<705
    assert 705*(1+dj*dj/18)**7<exp_bounds(F(263,40))[0]
    assert exp_bounds(F(329,50))[0]>720
    # Taylor derivative bounds for log terms at the two saddles.
    rho=3200/(jmin-rj)**6
    y=rho/(1-rho)
    assert (14+336*y+756*y*y+432*y**3)/(jmin-rj)**3<F(1,5)
    rho=(kmax+rk)**6/1000
    y=rho/(1-rho)
    assert (2+120*y+540*y*y+432*y**3)/(kmin-rk)**3<5
    assert aj/2-bj*rj>=cj and ak/2-bk*rk>=ck
    def curvature_log(terms,rho,r):
        a,b=terms; t=rho/(1-rho)
        return (1-r*r)/(1+r*r)**2-a*t-b*t*t
    assert (7*curvature_log((0,0),F(0),dj/xminj)
            - (42*3200/(xminj**6-3200)
               + 36*(3200/(xminj**6-3200))**2))>F(3,2)
    assert curvature_log((30,36),F(151,100)**6/1000,2*dk)>F(1,2)
    assert F(7,9)+F(6,9)*3200/(9**6-3200)-1<-F(1,5)
    # Envelope comparisons for the short contour and the outer/inner domains.
    def k_env(x):
        return x*(1+x**6/1000)*exp_bounds(-x)[1]
    assert max(k_env(F(3,2)),k_env(F(3)))*(1+dk*dk/F(9,2))**7<F(17,50)
    assert (F(1,2)+dk)*(1+(F(1,2)+dk)**6/1000)*exp_bounds(-F(1,2))[1]<F(17,50)
    assert F(17,50)<exp_bounds(-F(501,500)-F(1,225))[0]
    assert F(9,25)<exp_bounds(-F(501,500))[0]
    assert 1100*exp_bounds(F(-1))[1]+3**7*exp_bounds(F(-3))[1]<520
    assert 3*exp_bounds(F(-3))[1]+F(4**7,3080)*exp_bounds(F(-4))[1]<F(8,25)
    assert 4*exp_bounds(F(-4))[1]+F(7**7,3080)*exp_bounds(F(-7))[1]<F(8,25)
    def local_bound(s,r,b,c):
        return b/((s-r)*c*c)+1/(s*(s-r)*c)+1/(s*c*r)
    assert local_bound(jmin,rj,bj,cj)<11
    assert local_bound(kmin,rk,bk,ck)<43
    assert 3*factorial(4)*200**4<10000**3
    assert 5*factorial(4)*225**4<10000**3
    assert 6/amaxj/jmax**2>F(3,4)**2
    assert 6/amaxk/kmax**2>4
    assert F(13)/(F(3,4))<25 and F(44)/2<25
    assert F(1,3102)-1/(jmin**6-3100)>F(1,3300)
    assert kmax**6<2
    for s,v,a in [(jmin,6+dj,aj),(kmin,dk,ak)]:
        vprime=v*v/(6*s**7*a)
        assert vprime*(1/s+(7/s**2+14/s**3)/(2*a))<F(1,50)
    assert F(1,50)+F(1,300)<F(1,20)
    assert F(10000,3300)>F(1,10)

CUTOFF = 10000


def prime_table(limit):
    table = bytearray(b"\x01") * (limit + 1)
    table[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if table[p]:
            count = (limit - p * p) // p + 1
            table[p * p::p] = b"\x00" * count
    return table


def physical_polynomial_modulo(order, prime):
    coefficients = [1]
    rising = 1
    for k in range(1, order + 1):
        for j in range(order + 6 * (k - 1), order + 6 * k):
            rising = rising * j % prime
        coefficients.append(comb(order, k) * rising % prime)
    polynomial = nmod_poly(coefficients, prime)
    assert polynomial.degree() == order
    return polynomial


def check_finite_prefix():
    fallback_prime = 100003
    primes = prime_table(fallback_prime)
    assert primes[fallback_prime] and fallback_prime > 7 * (CUTOFF - 1)
    preceding = 2
    fallback = []
    small_gcds = []
    residual_count = 0
    for order in range(2, CUTOFF):
        if primes[order]:
            preceding = order
        offset = order - preceding
        if preceding > max(6, 7 * offset) and preceding**2 > 7 * order - 1:
            residual = physical_polynomial_modulo(offset, preceding)
            degree = residual.gcd(residual.derivative()).degree()
            if degree <= 1:
                residual_count += 1
                if degree == 1:
                    small_gcds.append([order, preceding, offset])
                continue
        polynomial = physical_polynomial_modulo(order, fallback_prime)
        degree = polynomial.gcd(polynomial.derivative()).degree()
        assert degree <= 1, (order, degree)
        fallback.append([order, degree])
    assert residual_count + len(fallback) == CUTOFF - 2
    return {
        "first_terminal": 2, "last_terminal": CUTOFF - 1,
        "preceding_prime_tests": residual_count,
        "degree_one_residual_gcds": small_gcds,
        "full_degree_prime": fallback_prime, "full_degree_fallbacks": fallback,
    }


@cache
def compositions(total, parts):
    if parts == 1:
        return ((total,),)
    return tuple((i, *tail) for i in range(total + 1)
                 for tail in compositions(total - i, parts - 1))


def forest_cumulants(gamma, gap, count):
    @cache
    def tree(level, degree):
        if not degree:
            return 1
        weight = gap * rising(gamma + level + 1, gap - 1)
        return weight * sum(
            prod(tree(level + i, amount) for i, amount in enumerate(allocation))
            for allocation in compositions(degree - 1, gap)
        )

    return [
        rising(gamma, gap) * sum(
            prod(tree(i, amount) for i, amount in enumerate(allocation))
            for allocation in compositions(degree, gap)
        )
        for degree in range(count)
    ]


def rising(shape, count):
    return prod(range(shape, shape + count))


def reciprocal(order, gamma, gap):
    return fmpz_poly([comb(order, k) * rising(gamma, gap * k) for k in range(order, -1, -1)])


def check_forest_and_resultant():
    count = 0
    for gap, terminal, offset in product((1, 2, 6), range(2, 7), (0, 2)):
        gamma = terminal + offset
        upper = reciprocal(terminal, gamma + 1, gap)
        fixed = reciprocal(terminal, gamma, gap)
        lower = reciprocal(terminal - 1, gamma, gap)
        assert fixed.derivative() == terminal * lower
        assert gamma * upper == (gamma + gap * terminal) * fixed - gap * fmpz_poly([0, 1]) * fixed.derivative()
        actual = lower.resultant(upper)
        natural = prod(rising(gamma + 1, gap * j - 1) for j in range(1, terminal))
        extra = (gamma + gap * terminal)**(terminal - 1)
        fixed_divisor = prod(rising(gamma, gap * j) for j in range(1, terminal))
        assert actual * fixed_divisor == natural * extra * lower.resultant(fixed)
        q = forest_cumulants(gamma, gap, terminal)
        for degree in (terminal - 1, terminal):
            literal = reciprocal(degree, gamma, gap)
            for z in range(degree + 1):
                matrix = fmpz_mat([
                    [z + q[0] if i == j else q[i - j] if i > j else -i - 1 if j == i + 1 else 0
                     for j in range(degree)]
                    for i in range(degree)
                ])
                assert matrix.det() == literal(z)
        count += 1
    return count


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    check_contour_constants()
    finite = check_finite_prefix()
    forest_controls = check_forest_and_resultant()
    proof = Path(__file__).with_name("GAP_SIX_RESULTANT_STRICT_POSITIVITY.md")
    report = {
        "schema_version": 1, "arithmetical_evidence": "internal exact replay",
        "analytic_terminal_cutoff": CUTOFF, "literal_coupling_and_forest_controls": forest_controls,
        "contour_constant_inequalities": "passed by exact rational Taylor bounds",
        "finite_prefix": finite,
        "source_sha256": {p.name: sha256(p.read_bytes()).hexdigest() for p in (proof, Path(__file__))},
    }
    if args.write_receipt:
        args.write_receipt.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("PASS: exact contour constants, forest determinant and all 9998 finite terminal orders")


if __name__ == "__main__":
    main()
