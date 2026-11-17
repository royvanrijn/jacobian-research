"""Small finite-field arithmetic on explicitly monic triangular quotients."""
from fractions import Fraction
from functools import lru_cache
from itertools import product
import sympy as S

def modular(value,prime):
    q=Fraction(value)
    return q.numerator*pow(q.denominator,-1,prime)%prime

def polynomial(expression,variables,prime):
    return {m:modular(co,prime) for m,co in S.Poly(expression,*variables).terms() if modular(co,prime)}

class Algebra:
    def __init__(self,variables,degrees,relations,prime):
        self.variables=variables;self.degrees=tuple(degrees);self.prime=prime
        self.basis=list(product(*(range(n) for n in degrees)))
        self.positions={m:i for i,m in enumerate(self.basis)}
        self.tails=[]
        for i,(degree,relation) in enumerate(zip(degrees,relations)):
            exact=S.Poly(relation,*variables)
            exact_lead=tuple(degree if j==i else 0 for j in range(len(degrees)))
            assert exact.coeff_monomial(exact_lead)==1
            assert all(mon==exact_lead or (mon[i]<degree and all(mon[j]==0 for j in range(i+1,len(degrees)))) for mon in exact.monoms())
            poly=polynomial(relation,variables,prime)
            lead=tuple(degree if j==i else 0 for j in range(len(degrees)))
            assert poly.pop(lead)==1
            assert all(m[i]<degree and all(m[j]==0 for j in range(i+1,len(degrees))) for m in poly)
            self.tails.append({m:-co%prime for m,co in poly.items()})
        self.reduce_monomial=lru_cache(None)(self._reduce_monomial)

    def _reduce_monomial(self,monomial):
        for i in range(len(monomial)-1,-1,-1):
            if monomial[i]<self.degrees[i]:continue
            rest=list(monomial);rest[i]-=self.degrees[i];out={}
            for tail,co in self.tails[i].items():
                for m,value in self.reduce_monomial(tuple(a+b for a,b in zip(rest,tail))).items():
                    out[m]=(out.get(m,0)+co*value)%self.prime
            return {m:value for m,value in out.items() if value}
        return {monomial:1}

    def reduce(self,poly):
        out={}
        for mon,co in poly.items():
            for m,value in self.reduce_monomial(mon).items():out[m]=(out.get(m,0)+co*value)%self.prime
        return {m:co for m,co in out.items() if co}

    def matrix(self,expression):
        poly=self.reduce(polynomial(expression,self.variables,self.prime))
        size=len(self.basis);matrix=[[0]*size for _ in range(size)]
        for j,monomial in enumerate(self.basis):
            value=self.reduce({tuple(a+b for a,b in zip(monomial,m)):co for m,co in poly.items()})
            for mon,co in value.items():matrix[self.positions[mon]][j]=co
        return matrix

def rank(matrix,prime):
    matrix=[row[:] for row in matrix];nrows=len(matrix);ncols=len(matrix[0]);r=0;pivots=[]
    for col in range(ncols):
        pivot=next((i for i in range(r,nrows) if matrix[i][col]),None)
        if pivot is None:continue
        matrix[r],matrix[pivot]=matrix[pivot],matrix[r]
        scale=pow(matrix[r][col],-1,prime)
        matrix[r]=[(v*scale)%prime for v in matrix[r]]
        for i in range(r+1,nrows):
            coefficient=matrix[i][col]
            if coefficient:matrix[i]=[(v-coefficient*w)%prime for v,w in zip(matrix[i],matrix[r])]
        pivots.append(col);r+=1
        if r==nrows:break
    return r,pivots
