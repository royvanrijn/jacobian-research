"""Doubled quartic-fibre contact matrix over a rational polynomial ring."""
from itertools import permutations
import sympy as S
from sympy.polys.rings import ring

R,u,v,w,z=ring('u,v,w,z',S.QQ)
t=S.Symbol('t')


def power_remainders(maximum=20):
    # t^8 plus the following lower coefficients is (t^4+u*t^3+v*t^2+w*t-z)^2.
    lower=[z*z,-2*w*z,w*w-2*v*z,2*v*w-2*u*z,
           v*v+2*u*w-2*z,2*u*v+2*w,u*u+2*v,2*u]
    result=[]
    for n in range(maximum+1):
        if n<8: result.append([R.one if i==n else R.zero for i in range(8)])
        else:
            result.append([-sum((lower[i]*result[n-8+i][j] for i in range(8)),R.zero)
                           for j in range(8)])
    return result


def remainders(polynomials):
    powers=power_remainders(max(S.degree(f,t) for f in polynomials))
    result=[]
    for f in polynomials:
        p=S.Poly(f,t)
        result.append([sum((R.ground_new(p.nth(n))*powers[n][j]
                           for n in range(p.degree()+1)),R.zero) for j in range(8)])
    return result


def contact_matrix(p,q):
    raw=remainders((q,p*p,p*q,q*q))
    pp=S.Poly(p,t)
    assert pp.degree()==6 and pp.nth(6)==1
    matrix=[[column[j]-column[6]*R.ground_new(pp.nth(j)) for column in raw]
            for j in [1,2,3,4,5,7]]
    return matrix,raw


def determinant4(matrix):
    result=R.zero
    for perm in permutations(range(4)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
        term=R.ground_new(sign)
        for i,j in enumerate(perm): term*=matrix[i][j]
        result+=term
    return result


def selected_minor(matrix,rows):
    return determinant4([[matrix[i][j] for j in range(4)] for i in rows])


def z_coefficient(polynomial,exponent):
    return R.from_dict({m[:-1]+(0,):c for m,c in polynomial.items() if m[-1]==exponent})
