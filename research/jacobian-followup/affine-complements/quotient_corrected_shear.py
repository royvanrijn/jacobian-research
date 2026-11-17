"""Exact arithmetic interface for PF127; geometry is proved in its source note."""
from math import gcd,lcm


def shear(e,f,n,delta):
    if e < 2 or f < 1 or n < 1 or delta < 0 or gcd(e,n) != 1:
        raise ValueError("invalid ramification, residue, quotient or delta data")
    return (f*(2*delta-1)+pow(n,-1,e)) % e


def quotient_order(e,f,a,delta):
    if e < 2 or f < 1 or not 0 <= a < e or delta < 0:
        raise ValueError("invalid peripheral data")
    coefficient = (a-f*(2*delta-1)) % e
    if gcd(coefficient,e) != 1:
        return None
    return pow(coefficient,-1,e)


def orbit_invariant_factors(e,f,a):
    if e < 1 or f < 1 or not 0 <= a < e:
        raise ValueError("invalid normalized torus lattice")
    d = gcd(gcd(e,f),a)
    return d,e*f//d


def combine_congruences(r,m,s,n):
    """Combine T=r mod m and T=s mod n, without a parameter sweep."""
    g = gcd(m,n)
    if (s-r) % g:
        return None
    nn = n//g
    k = ((s-r)//g*pow(m//g,-1,nn)) % nn if nn > 1 else 0
    period = lcm(m,n)
    return (r+m*k) % period,period


def shared_image_class(packets):
    """Necessary odd T=2 delta-1 progression for one common branch image.

    A passing result supplies no curve or actual finite cover.
    """
    result = (1,2)
    for e,f,a,n in packets:
        if e < 2 or f < 1 or not 0 <= a < e or not 1 <= n < e or gcd(e,n) != 1:
            return None
        if n == 1 and f != 1:
            return None  # Hereditary original-Keller location, not a local-germ rule.
        rhs = a-pow(n,-1,e)
        g = gcd(e,f)
        if rhs % g:
            return None
        modulus = e//g
        residue = (rhs//g*pow(f//g,-1,modulus)) % modulus if modulus > 1 else 0
        result = combine_congruences(*result,residue,modulus)
        if result is None:
            return None
    return result


def joint_permutations(e,f,a):
    """An abstract commuting peripheral pair, not a Keller realization."""
    orbit_invariant_factors(e,f,a)
    mu,lam = [],[]
    for j in range(f):
        for i in range(e):
            mu.append(j*e+(i+1)%e)
            lam.append((j+1)*e+i if j+1 < f else (i-a)%e)
    return mu,lam
