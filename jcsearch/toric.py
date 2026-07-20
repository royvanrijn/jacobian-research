"""Unimodular two-divisor toric charts and exact monomial cones."""
from dataclasses import dataclass
import sympy as sp

x,y,u,v=sp.symbols("x y u v")

@dataclass(frozen=True)
class ToricChart:
    # u=x^a y^b, v=x^c y^d
    a:int;b:int;c:int;d:int
    def __post_init__(self):
        if abs(self.a*self.d-self.b*self.c)!=1: raise ValueError("exponent matrix must be unimodular")
    @property
    def matrix(self): return sp.Matrix([[self.a,self.b],[self.c,self.d]])
    @property
    def expressions(self): return (x**self.a*y**self.b,x**self.c*y**self.d)
    @property
    def jacobian(self): return sp.factor(sp.Matrix(self.expressions).jacobian((x,y)).det())
    def pullback_exponent(self,ij):
        i,j=ij;return (self.a*i+self.c*j,self.b*i+self.d*j)
    def chart_exponent(self,pq):
        vec=self.matrix.T.inv()*sp.Matrix(pq);return tuple(map(int,vec))
    def polynomializable(self,ij): return all(e>=0 for e in self.pullback_exponent(ij))
    def pullback(self,f): return sp.cancel(f.subs({u:self.expressions[0],v:self.expressions[1]},simultaneous=True))
    @property
    def reciprocal_outer_jacobian(self):
        det=self.a*self.d-self.b*self.c
        return sp.factor((1/self.jacobian).subs({x:u**(self.d//det)*v**(-self.b//det),y:u**(-self.c//det)*v**(self.a//det)}))
    @property
    def signature(self):
        columns=tuple(sorted((tuple(self.matrix.col(0)),tuple(self.matrix.col(1)))))
        return {"kind":"unimodular_toric_two_divisor","abs_det":1,"sorted_columns":columns,
                "jacobian_divisor_exponents":(self.a+self.c-1,self.b+self.d-1)}

def catalog():
    return (ToricChart(-1,1,2,-1),ToricChart(-1,1,1,-2),ToricChart(-2,1,1,-1))

def cone_support(chart,bound=8,limit=15):
    required=[chart.chart_exponent((0,0)),chart.chart_exponent((1,0)),chart.chart_exponent((0,1))]
    candidates=[]
    for i in range(-bound,bound+1):
      for j in range(-bound,bound+1):
        if chart.polynomializable((i,j)):
          pq=chart.pullback_exponent((i,j));candidates.append(((sum(pq),max(pq),abs(i)+abs(j)),(i,j)))
    selected=list(dict.fromkeys(required))
    for _,ij in sorted(candidates):
      if ij not in selected:selected.append(ij)
      if len(selected)>=limit:break
    return tuple(selected)
