"""Canonical invariant signatures for chart deduplication.

These invariants are necessary, not sufficient, for equivalence. Exact
triangular-chart equivalence is handled separately by TriangularChart.signature.
"""
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

def _support(expr,variables):
    num,den=sp.fraction(sp.cancel(expr));out=[]
    for poly,sign in ((num,1),(den,-1)):
      for monomial,_ in sp.Poly(poly,*variables).terms(): out.append(tuple(sign*e for e in monomial))
    return out

def valuation(expr,variables,weights):
    num,den=sp.fraction(sp.cancel(expr))
    def val(poly): return min(sum(e*w for e,w in zip(m,weights)) for m,_ in sp.Poly(poly,*variables).terms())
    return val(num)-val(den)

def chart_signature(expressions,variables,exceptional_weight=(1,0),sample_box=range(-2,3)):
    jac=sp.factor(sp.Matrix(expressions).jacobian(variables).det())
    supports=sum((_support(e,variables) for e in expressions),[])
    lattice=sp.Matrix(supports) if supports else sp.zeros(0,len(variables))
    snf=smith_normal_form(lattice,domain=sp.ZZ) if lattice.rows else lattice
    diagonal=tuple(abs(int(snf[i,i])) for i in range(min(snf.shape)) if snf[i,i])
    vals=tuple(valuation(e,variables,exceptional_weight) for e in expressions)
    swapped=tuple(reversed(vals))
    U,V=sp.symbols("U V")
    cone=[]
    for i in sample_box:
      for j in sample_box:
        pulled=sp.cancel((U**i*V**j).subs({U:expressions[0],V:expressions[1]},simultaneous=True))
        if sp.denom(pulled)==1: cone.append((i,j))
    return {"coordinate_valuations":min(vals,swapped),"jacobian_valuation":valuation(jac,variables,exceptional_weight),
            "lattice_snf":diagonal,"sampled_polynomializable_cone":tuple(cone)}
