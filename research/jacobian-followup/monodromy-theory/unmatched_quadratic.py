"""Canonical quadratic chart forms and independent-boundary controls."""
import sympy as S
def canonical_forms(ps):
    u,x,y,z=ps
    return {"A":S.Matrix([u*u,u*x,u*y+x*x,u*z+y*y]),
            "B":S.Matrix([u*u,u*x,u*y+x*x,u*z+x*y]),
            "C":S.Matrix([u*u,u*x,u*y,u*z+x*x]),
            "D":S.Matrix([u*u,u*x,u*y,u*z+x*y])}
def chart_inverses(xs):
    x,y,z=xs
    return {"A":S.Matrix([x,y-x*x,z-(y-x*x)**2]),
            "B":S.Matrix([x,y-x*x,z-x*(y-x*x)]),
            "C":S.Matrix([x,y,z-x*x]),
            "D":S.Matrix([x,y,z-x*y])}
def projective_direction(V,ps):
    p=S.Matrix(ps);u=ps[0]
    return (u*S.Matrix(V)-V[0]*p).applyfunc(S.expand)
def swap_output(r,i):
    A=S.eye(r)
    A.row_swap(0,i)
    return A
