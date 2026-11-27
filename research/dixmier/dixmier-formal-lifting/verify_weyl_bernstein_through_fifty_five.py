"""Exact corner and successor lattice checks for the Weyl cutoff through 55."""
from fractions import Fraction as F
from math import gcd

def corners(total):
    out = []
    for u in range(3, (total + 1) // 2):
        v = total - u
        div = gcd(u, v)
        for k in range(1, div):
            f1, f2 = k*u//div, k*v//div
            if f1 < 2:
                continue
            d = gcd(f1-1, f2-1)
            rho, sigma = (f2-1)//d, (1-f1)//d
            for r0 in range(1, u):
                for s0 in range(r0):
                    if rho*(r0-u)+sigma*(s0-v):
                        continue
                    off = F(v-s0,rho)-s0
                    next_x, next_y = r0-off*F(sigma,rho), s0+off
                    if d == 1 and next_y >= 2 and next_x == next_y-F(1,rho):
                        continue
                    out.append(((u,v),(r0,s0),(f1,f2),(rho,sigma)))
    return out

assert corners(18) == [((6,12),(2,0),(2,4),(3,-1))]
assert corners(27) == [((9,18),(3,0),(2,4),(3,-1))]

def successors(ell, mu, beta, previous):
    upper = 1-F(beta,ell*mu)
    single = [F(k,ell) for k in range(1,ell) if previous<F(k,ell)<upper]
    multiple = []
    for M in range(2,beta*mu+2):
        slope = 1-F(beta*M,ell*mu*(M-1))
        D = (ell*slope).denominator
        if previous < slope < upper and (M-1) % D == 0:
            for j in range(mu):
                if (mu-j) % D == 0:
                    multiple.append((M,D,j,slope))
    return single, multiple

assert successors(3,4,2,F(1,3)) == (
    [F(2,3)], [(2,1,j,F(2,3)) for j in range(4)] + [(4,3,1,F(7,9))])
assert successors(3,3,1,F(2,3)) == ([],[(3,2,1,F(5,6))])
assert successors(3,6,3,F(1,3)) == (
    [F(2,3)], [(2,1,j,F(2,3)) for j in range(6)] +
    [(4,3,0,F(7,9)),(4,3,3,F(7,9)),(6,5,1,F(4,5))])
assert successors(3,5,2,F(2,3)) == ([],[(5,2,1,F(5,6)),(5,2,3,F(5,6))])
assert successors(3,4,1,F(2,3)) == ([],[(4,3,1,F(8,9))])
assert successors(9,2,1,F(7,9)) == ([F(8,9)],[(2,1,0,F(8,9)),(2,1,1,F(8,9))])
assert successors(3,3,1,F(7,9)) == ([],[(3,2,1,F(5,6))])
