from fractions import Fraction as F
from math import gcd

def corners(total):
    out = []
    for u in range(3, (total + 1) // 2):
        v = total - u
        d = gcd(u, v)
        for k in range(1, d):
            f1, f2 = k*u//d, k*v//d
            if f1 < 2:
                continue
            cut = gcd(f1-1, f2-1)
            rho, sigma = (f2-1)//cut, (1-f1)//cut
            for r in range(1, u):
                for h in range(r):
                    if rho*(r-u) + sigma*(h-v):
                        continue
                    delta = F(v-h, rho)-h
                    next_x, next_y = r-delta*F(sigma,rho), h+delta
                    if cut == 1 and next_y >= 2 and next_x == next_y-F(1,rho):
                        continue
                    out.append(((u,v),(r,h),(f1,f2),(rho,sigma)))
    return out

assert corners(16) == [
    ((4,12),(1,0),(3,9),(4,-1)),
    ((6,10),(1,0),(3,5),(2,-1)),
]
assert corners(17) == []
assert corners(25) == [
    ((5,20),(1,0),(4,16),(5,-1)),
    ((10,15),(3,1),(2,3),(2,-1)),
]

def scalar(ell, mu, beta):
    out = []
    for N in range(5, 2*beta + 1, 2):
        slope = 1-F(beta*N,ell*(N*mu-1))
        D = (ell*slope).denominator
        b = N-2
        if F(1,ell) < slope < 1-F(beta,ell*mu) and D <= 2*mu:
            if 2*mu % D == 0 or b*mu % D == 0:
                out.append((b,slope,D))
    return out

assert scalar(5,2,3) == [(3,F(2,3),3)]
assert scalar(5,3,7) == [(3,F(1,2),2)]
assert scalar(5,4,11) == []

def commuting(ell, mu, beta):
    out = []
    for M in range(2,beta*mu+2):
        slope = 1-F(beta*M,ell*mu*(M-1))
        D = (ell*slope).denominator
        if F(1,ell) < slope < 1-F(beta,ell*mu) and (M-1) % D == 0:
            for j in range(mu):
                if (mu-j) % D == 0:
                    out.append((M,D,j))
    return out

assert commuting(5,2,3) == [(2,1,0),(2,1,1),(4,1,0),(4,1,1)]
assert commuting(5,3,7) == [(3,2,1),(15,2,1)]
assert commuting(5,4,11) == [(4,3,1),(12,1,0),(12,1,1),(12,1,2),(12,1,3)]
assert commuting(2,7,2) == [(7,3,1),(7,3,4)]
