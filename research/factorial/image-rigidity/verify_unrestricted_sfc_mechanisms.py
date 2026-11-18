#!/usr/bin/env python3
"""Exact controls for UNRESTRICTED_SFC_MECHANISMS.tex."""
from fractions import Fraction
from math import comb, factorial, prod
from pathlib import Path
import hashlib, json

def moments(nmax):
    return [sum(comb(n,k)*factorial(2*k) for k in range(n+1))
            for n in range(nmax+1)]


def rising(q, k):
    return prod(q+j for j in range(k))


def gamma_block_moment(columns, k, n):
    """Mixed moment from independent Gamma/Beta moment formulas."""
    b = [sum(a[i]*kj for a, kj in zip(columns, k))
         for i in range(len(columns[0]))]
    return (Fraction(n)**(sum(b)-sum(k))
            * prod(rising(n+1, kj) for kj in k)
            * prod(Fraction(factorial(z), rising(n+1, z)) for z in b))


def euler_block_multiplier(columns, k, n):
    """The finite product of diagonal Euler factors, separately evaluated."""
    b = [sum(a[i]*kj for a, kj in zip(columns, k))
         for i in range(len(columns[0]))]
    value = Fraction(n)**(sum(b)-sum(k))
    for ell in range(1, n+1):
        value *= prod(Fraction(ell+kj, ell) for kj in k)
        value /= prod(Fraction(ell+z, ell) for z in b)
    return value


def check_gamma_controls():
    f = lambda k: (k+1)**2*(k+2)
    coefficients = [[f(k)*f(m-k) for k in range(m+1)]
                    for m in range(5)]
    assert coefficients[0] == [4]
    assert coefficients[1] == [24, 24]
    assert coefficients[2] == [72, 144, 72]
    values = [sum((-1)**k*term for k, term in enumerate(row))
              for row in coefficients]
    assert values == [4, 0, 0, 0, -24]
    scalar = [Fraction(factorial(k)*f(k), f(0)) for k in range(4)]
    assert scalar == [1, 6, 36, 240]
    assert scalar[2] == scalar[1]**2
    assert scalar[3] != scalar[1]**3
    bad_f = lambda k: f(k)-(1 if k == 2 else 0)
    assert sum((-1)**k*bad_f(k)*bad_f(2-k) for k in range(3)) != 0

    multiplier_controls = 0
    samples = [
        (((0,), (1,), (2,)), [(0,0,0), (1,2,1), (0,4,0)]),
        (((1,0), (0,1), (1,1)), [(1,1,1), (2,0,2), (0,0,3)]),
        (((0,0), (2,1)), [(0,0), (2,3), (1,1)]),
    ]
    for columns, ks in samples:
        for k in ks:
            for n in (1, 2, 5):
                assert (euler_block_multiplier(columns, k, n)
                        == gamma_block_moment(columns, k, n)
                           / prod(factorial(kj) for kj in k))
                multiplier_controls += 1

    columns = ((0,), (1,), (2,))
    for n in range(1, 33):
        a = gamma_block_moment(columns, (2,0,2), n)
        b = gamma_block_moment(columns, (1,2,1), n)
        c = gamma_block_moment(columns, (0,4,0), n)
        target = Fraction(144*(n+2), (n+3)*(n+4))
        assert a-2*b+c == target > 0
        assert a-b+c != target
    assert target < Fraction(144, 32)
    assert (gamma_block_moment(columns, (2,0,2), 1)
            - 2*gamma_block_moment(columns, (1,2,1), 1)
            + gamma_block_moment(columns, (0,4,0), 1)) == Fraction(108,5)
    return {
        "euler_values_at_minus_one_one": values,
        "finite_block_multiplier_controls": multiplier_controls,
        "toric_degeneration_controls": 32,
        "gamma_corrupt_coefficients_rejected": True,
    }



def valuation(a, p):
    count = 0
    while a % p == 0:
        a //= p
        count += 1
    return count


def prime_factors(a):
    factors = set()
    p = 2
    while p*p <= a:
        while a % p == 0:
            factors.add(p)
            a //= p
        p += 1
    if a > 1:
        factors.add(a)
    return factors


def determinant(rows):
    a = [[Fraction(x) for x in row] for row in rows]
    value = Fraction(1)
    for j in range(len(a)):
        pivot = next(i for i in range(j, len(a)) if a[i][j])
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            value = -value
        factor = a[j][j]
        value *= factor
        a[j] = [x/factor for x in a[j]]
        for i in range(j+1, len(a)):
            factor = a[i][j]
            a[i] = [x-factor*y for x, y in zip(a[i], a[j])]
    return value


def check_diagonal_and_saddle_controls():
    a, b = [1,4,24], [1,6,72,720]
    assert [(-2*a[k] if k < 3 else 0)
            + (30*a[k-1] if k else 0)
            + ([3,-16][k] if k < 2 else 0) for k in range(4)] == b
    assert sum(Fraction(c)*Fraction(3,16)**k for k,c in enumerate(a)) == Fraction(83,32)
    sylvester = []
    for coefficients, shifts in ((a[::-1],3), (b[::-1],2)):
        for shift in range(shifts):
            sylvester.append([0]*shift+coefficients
                             +[0]*(5-shift-len(coefficients)))
    resultant = determinant(sylvester)
    assert resultant == 382464 and resultant % 5 != 0

    table = {
        2: (Fraction(-3,2),Fraction(-3,2),Fraction(-1),Fraction(-3,2)),
        3: (Fraction(0),Fraction(-1),Fraction(-2,3),Fraction(-2,3)),
        5: (Fraction(0),Fraction(0),Fraction(0),Fraction(-1)),
        7: (Fraction(0),Fraction(0),Fraction(0),Fraction(0)),
    }
    exceptional = set().union(*(prime_factors(c) for c in a+b))
    assert exceptional == {2,3,5}
    assert a[1]**2-a[0]*a[2] == -8
    assert b[1]**2-b[0]*b[2] == -36
    assert prime_factors(8) | prime_factors(36) <= exceptional
    for p, expected in table.items():
        row = []
        for coefficients in (a,b):
            v = [valuation(c,p) for c in coefficients]
            m = len(v)-1
            row += [max(-Fraction(v[k],k) for k in range(1,m+1)),
                    min(Fraction(v[k]-v[m],m-k) for k in range(m))]
        assert tuple(row) == expected
        l2,u2,l3,u3 = row
        assert l2 >= u3 and l3 >= u2
    for coefficients, primes in ((a,(2,3)), (b,(2,3,5))):
        for p in primes:
            v = [valuation(c,p) for c in coefficients]
            assert len({v[k+1]-v[k] for k in range(len(v)-1)}) > 1

    g = lambda r: -7*r**6+12*r**5-19*r+14
    gp = lambda r: -42*r**5+60*r**4-19
    assert g(Fraction(5,4)) == Fraction(689,4096) > 0
    assert g(Fraction(4,3)) == -Fraction(70,729) < 0
    assert gp(Fraction(5,4)) == -Fraction(353,512) < 0
    assert 8-7*Fraction(5,4) < 0
    u = Fraction(1,3)
    power_sum = 3*u**5
    diagonal = -7/u**2+20*u**3/power_sum-(5*u**4/power_sum)**2
    off_diagonal = -(5*u**4/power_sum)**2
    assert diagonal == -28 and off_diagonal == -25
    for v in ((1,-1,0),(1,1,-2)):
        quadratic = sum(v[i]*v[j]*(diagonal if i==j else off_diagonal)
                        for i in range(3) for j in range(3))
        assert quadratic == -3*sum(x*x for x in v) < 0
    assert -7+Fraction(20,3) == -Fraction(1,3)
    # An intentional sign corruption loses the isolating-interval test.
    assert not (-g(Fraction(5,4)) > 0 and -g(Fraction(4,3)) < 0)
    return {
        "diagonal_normal_form_resultant": int(resultant),
        "diagonal_threshold_prime_classes": ["2","3","5","all other primes"],
        "four_maxima_rational_controls": True,
        "saddle_sign_corruption_rejected": True,
    }


def main():
    route_controls = check_diagonal_and_saddle_controls()
    gamma_controls = check_gamma_controls()
    ms=moments(32)
    for n in range(2,33):
        assert ms[n]-2*n*(2*n-1)*ms[n-1]+4*n*(n-1)*ms[n-2]==1
        assert ms[n]>=factorial(2*n)
    for n in range(3,33):
        assert (ms[n]-(4*n*n-2*n+1)*ms[n-1]
                +(n-1)*(8*n-6)*ms[n-2]
                -4*(n-1)*(n-2)*ms[n-3])==0
    assert ms[2]-(2*2*(2*2-1)+1)*ms[1]+4*2*(2-1)*ms[0]!=1
    print(json.dumps({
        "status":"PASS_exact_controls_only",
        "inhomogeneous_controls":31,
        "homogeneous_controls":30,
        "first_moments":ms[:5],
        "corrupt_coefficient_rejected":True,
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        **gamma_controls,
        **route_controls,
        "scope":("Exact finite recurrence, operator, arithmetic and saddle controls only; "
                 "all-index proofs are written. Neither finite-block nor "
                 "ordinary unrestricted SFC is established.")
    },indent=2))

if __name__=="__main__":
    main()
