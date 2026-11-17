"""Exact continuants used in the all-weight quadratic-pair proof."""
import sympy as S
from quadratic_scalar_obstruction import curved_sector

def leading_matrix(weight):
    n=weight//2+1
    matrix=S.zeros(n)
    for j in range(n):
        i=weight-2*j
        matrix[j,j]=2*(i+1)*(j+2)
        if j+1<n:
            matrix[j,j+1]=(j+2)*(j+3)
            matrix[j+1,j]=i*(i-1)
    return matrix

def bottom_continuants(weight,k):
    n=weight//2+1
    H=[S.Poly(1,k,domain=S.QQ)]
    for j in range(n):
        i=weight-2*j
        E=2*(i+1)*(i+2+k*(j+2))
        nxt=S.Poly(E,k)*H[-1]
        if j:
            D=(i+4)*(i+3)+2*k*(i+3)*(j+1)+k*k*(j+1)*(j+2)
            previous_i=i+2
            nxt-=S.Poly(previous_i*(previous_i-1)*D,k)*H[-2]
        H.append(nxt)
    return H
