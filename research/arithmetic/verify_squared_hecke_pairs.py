#!/usr/bin/env python3
"""Finite exact arithmetic checks. No asymptotic estimate is tested."""
if not __debug__:
    raise RuntimeError('Exact arithmetic checks require assertions; disable optimization.')

from fractions import Fraction as Q
from math import gcd
import json

LIMIT = 40
MAXIMUM = LIMIT * LIMIT

def tau_values(limit):
    # Derived from q*d/dq log(Delta/q)=-24 sum sigma_1(j) q^j.
    sigma = [0] * (limit + 1)
    for d in range(1, limit + 1):
        for n in range(d, limit + 1, d):
            sigma[n] += d
    tau = [0] * (limit + 1)
    tau[1] = 1
    for n in range(2, limit + 1):
        numerator = -24 * sum(sigma[j] * tau[n-j] for j in range(1, n))
        assert numerator % (n-1) == 0
        tau[n] = numerator // (n-1)
    assert [tau[n] for n in range(1, 11)] == [1,-24,252,-1472,4830,-6048,-16744,84480,-113643,-115920]
    return tau

TAU = tau_values(MAXIMUM)
assert TAU[11] == 534612
assert Q(TAU[11]**2,11**11) > 1

def c(n): return Q(TAU[n*n], n**11)
def a(n): return Q(TAU[n]**2, n**11)
def valuation(n,p):
    assert n > 0
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k

def local_a(x,k):
    values = [Q(1), x, (x-1)**2]
    for j in range(3,k+1):
        values.append((x-1)*values[-1]-(x-1)*values[-2]+values[-3])
    return values[k]

def local_total(x,z):
    return (1+z)/((1-z)*(1-(x-2)*z+z*z))

def local_pair(p,h,x,y):
    z = Q(1,p)
    v = valuation(h,p)
    av = local_a(x,v)
    bv = local_a(y,v)
    below = (1-z)*sum(local_a(x,k)*local_a(y,k)*z**k for k in range(v))
    equal = z**v*(1-2*z)*av*bv
    tailx = local_total(x,z)-sum(local_a(x,k)*z**k for k in range(v+1))
    taily = local_total(y,z)-sum(local_a(y,k)*z**k for k in range(v+1))
    return below+equal+(1-z)*(bv*tailx+av*taily)

def square_tail(z,k):
    return z**k*((k+1)**2/(1-z)+2*(k+1)*z/(1-z)**2+z*(1+z)/(1-z)**3)

classification = 0
for n in range(1,LIMIT+1):
    assert a(n) == sum(c(d) for d in range(1,n+1) if n % d == 0)
    classification += 1

crt = 0
for X in (31,64,127):
    for Z in range(1,9):
        truncated = lambda n: sum(c(d) for d in range(1,Z+1) if n % d == 0)
        observed = sum(truncated(n)*truncated(n+1) for n in range(1,X+1))
        main = sum(c(d)*c(e)/Q(d*e) for d in range(1,Z+1) for e in range(1,Z+1) if gcd(d,e)==1)
        expanded = Q(0)
        for d in range(1,Z+1):
            for e in range(1,Z+1):
                count = sum(n % d == 0 and (n+1) % e == 0 for n in range(1,X+1))
                if gcd(d,e)>1: assert count == 0
                else: assert abs(Q(count)-Q(X,d*e)) <= 1
                expanded += c(d)*c(e)*count
                crt += 1
        assert observed == expanded
        assert abs(observed-X*main) <= sum(abs(c(d)) for d in range(1,Z+1))**2

local_enclosures = normalization = prime_powers = 0
for p in (2,3,5,7):
    x = Q(TAU[p]**2,p**11)
    assert 0 <= x <= 4
    z = Q(1,p)
    power = p
    exponent = 1
    while power <= MAXIMUM:
        assert local_a(x,exponent) == a(power)
        prime_powers += 1
        exponent += 1
        power *= p
    assert a(p*p) == (a(p)-1)**2
    m = (1-z)*local_total(x,z)
    assert m == (1+z)/(1-(x-2)*z+z*z)
    assert local_pair(p,1,x,x) == 2*m-1
    assert (2*m-1)/(m*m) == 1-((x-1-z)/(p+1))**2
    assert 0 < (2*m-1)/(m*m) <= 1
    normalization += 1
    for y in (Q(0),Q(1),Q(2),Q(4),x):
        for h in range(1,9):
            v = valuation(h,p)
            for depth in (v+1,v+2):
                modulus = p**depth
                weightsx = [None if r==0 else local_a(x,valuation(r,p)) for r in range(modulus)]
                weightsy = [None if r==0 else local_a(y,valuation(r,p)) for r in range(modulus)]
                lower = sum(weightsx[r]*weightsy[(r+h)%modulus]
                            for r in range(modulus)
                            if weightsx[r] is not None and weightsy[(r+h)%modulus] is not None)/modulus
                upper = lower+(1-z)*(local_a(x,v)+local_a(y,v))*square_tail(z,depth)
                predicted = local_pair(p,h,x,y)
                assert lower <= predicted <= upper, (p,h,depth,x,y)
                local_enclosures += 1

print(json.dumps({'status':'PASS','tau_coefficients_from_q_log_derivative':MAXIMUM,
                  'normalized_hecke_convolution_checks':classification,
                  'finite_CRT_divisor_pair_checks':crt,
                  'rational_local_Haar_enclosures':local_enclosures,
                  'Delta_shift_one_prime_normalizations':normalization,
                  'local_prime_power_checks_against_q_expansion':prime_powers,
                  'literal_unit_bound_counterexample':'tau(11)=534612; a(11)>1',
                  'asymptotic_tested':False,'arithmetic':'exact integers and fractions'},indent=2))
