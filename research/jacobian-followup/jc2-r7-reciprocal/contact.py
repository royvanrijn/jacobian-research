"""Exact weighted minors for three conic contacts with one fixed infinity point.

The rational cubic base is t^2+v*t+w/(t-h), with w nonzero.
Coefficient filtering is exact by the nonnegative grading described in PROOF.md.
"""
from pathlib import Path
from itertools import permutations,product,combinations
import sys,json,time
import sympy as S
from sympy.polys.rings import ring

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'archive/non-elliptic/scripts'))
from verify_f2_degree_6_10_gap import construct_obstructions

R,h0,v0,w0,Z=ring('h,v,w,Z',S.QQ)
A=v0-h0;B0=-h0*v0-Z;C=h0*Z+w0
lower=[C*C,2*B0*C,B0*B0+2*A*C,2*C+2*A*B0,A*A+2*B0,2*A]
remainders=[]
for n in range(17):
    if n<6:remainders.append([R.one if j==n else R.zero for j in range(6)])
    else:remainders.append([-sum((lower[i]*remainders[n-6+i][j] for i in range(6)),R.zero) for j in range(6)])

P,*gens=ring('a,b,c,d,B,D,E,h,v,w',S.QQ)
a,b,c,d,B,D,E,h,v,w=gens
weights=(2,3,4,5,2,6,8,1,1,3)
t=S.Symbol('t')
_,rho,_,_,powers=construct_obstructions()
p=t**6+a.as_expr()*t**4+b.as_expr()*t**3+c.as_expr()*t**2+d.as_expr()*t
q=S.expand(sum(co*sum(powers[j][n].as_expr()*t**(2*j-n) for n in range(2*j+1))
               for j,co in [(5,1),(4,B.as_expr()),(2,D.as_expr()),(1,E.as_expr())]))
degrees=(6,10,12,16)
columns=[]
for f in [p,q,p*p,p*q]:
    poly=S.Poly(f,t);columns.append([P.from_expr(poly.nth(j)) for j in range(poly.degree()+1)])


def remco(n,r,k):
    return P.from_dict({(0,)*7+m[:3]:co for m,co in remainders[n][r].items() if m[-1]==k})


def coefficient(rows,power,subs=None):
    total=sum(degrees)-sum(rows)-2*power
    entries=[]
    for r in rows:
        row=[]
        for col,deg in enumerate(degrees):
            cs={}
            for k in range(9):
                if not 0<=deg-r-2*k<=total:continue
                f=sum((columns[col][n]*remco(n,r,k) for n in range(max(0,deg-total),deg+1)),P.zero)
                if subs:f=f.compose(subs)
                if f:cs[k]=f
            row.append(cs)
        entries.append(row)
    result=P.zero
    for perm in permutations(range(4)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
        for chosen in product(*(list(entries[i][j].items()) for i,j in enumerate(perm))):
            if sum(k for k,f in chosen)!=power:continue
            term=P.ground_new(sign)
            for k,f in chosen:term*=f
            result+=term
    return result


