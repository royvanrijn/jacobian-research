"""Uniform higher-gap adjoints and exact rank exporters (no all-gap closure)."""
import sympy as s
from itertools import combinations

def auxiliary(n,gamma,d,z):
    B=s.Integer(1)
    for k in range(n):
        B=s.expand((z+d*k+gamma+d-1)*B-z*s.diff(B,z))
    return s.Poly(B,z)

def tilted(n,gamma,d,t):
    return s.Poly(sum(s.binomial(n,k)*s.rf(gamma,d*k)*t**k for k in range(n+1)),t)

def reversed_tilted(n,gamma,d,z):
    return s.Poly(sum(s.binomial(n,k)*s.rf(gamma,d*k)*(-1)**((d+1)*k)*z**(d*(n-k))
                      for k in range(n+1)),z)

def Q_power(B,k,z):
    return s.expand(sum((-1)**j*s.binomial(k,j)*s.diff(B,z,j) for j in range(k+1)))

def polynomial_operator(B,n,M,d,z):
    return s.expand(z*(Q_power(B,d+1,z)-Q_power(B,1,z))
                    +(M-1)*(Q_power(B,d,z)-B)+d*n*B)

def scalar_operator(R,n,M,d,z):
    delta=lambda expr:z*s.diff(expr,z)
    v=R
    for ell in range(d):v=s.expand(delta(v)-(M+ell)*v)
    return s.cancel((delta(v)-z**d*(delta(R)-d*n*R))/z**(d+1))

def roots_of_unity(d):
    if d==2:u=-s.Integer(1)
    elif d==3:u=(-1+s.sqrt(-3))/2
    elif d==4:u=s.I
    else:u=s.exp(2*s.pi*s.I/d)
    return [s.simplify(u**j) for j in range(d)]

def adjoint_matrix(B,d,z,full=False):
    roots=roots_of_unity(d)
    return s.Matrix([[s.expand(w**k*Q_power(B,k,z).subs(z,w*z))
                       for k in range(d if full else d-1)] for w in roots])

def rank_minors(B,d,z):
    A=adjoint_matrix(B,d,z)
    return [s.expand(A[list(rows),:].det(method='domain-ge')) for rows in combinations(range(d),d-1)]

def exact_rank_gcd(B,d,z):
    minors=[s.Poly(v,z,extension=True) for v in rank_minors(B,d,z)]
    gcd=minors[0]
    for v in minors[1:]:gcd=s.gcd(gcd,v)
    return gcd.monic()

def gap_three_W(B,z):
    u=(1+s.sqrt(-3))/2;v=s.conjugate(u)
    p=B.subs(z,u*z);q=B.subs(z,v*z)
    val=s.expand(p*q-(u*s.diff(B,z).subs(z,u*z)*q-v*s.diff(B,z).subs(z,v*z)*p)/s.sqrt(-3))
    return s.Poly(s.simplify(val),z)

def cyclic_parts(poly,d,z,u):
    return [s.Poly(sum(poly.nth(k)*u**((k-r)//d) for k in range(r,poly.degree()+1,d)),u)
            for r in range(d)]
