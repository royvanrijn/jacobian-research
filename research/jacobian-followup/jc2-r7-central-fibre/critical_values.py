"""Good-prime critical-value controls for both contact34 coefficient roots."""
from pathlib import Path
import sys
import sympy as S

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'jc2-r7-elliptic-contact'))
import algebra as E
sys.path.insert(0,str(HERE.parent/'jc2-r7-quintic-incidence'))
from finite import modular

def multiply(a,b,prime):
    n=len(a)
    return [[sum(a[i][k]*b[k][j] for k in range(n))%prime for j in range(n)] for i in range(n)]

def charpoly(matrix,prime):
    """Newton traces, with every division a checked good-prime unit."""
    n=len(matrix);assert n<prime
    power=[[int(i==j) for j in range(n)] for i in range(n)]
    traces=[]
    for _ in range(n):
        power=multiply(power,matrix,prime)
        traces.append(sum(power[i][i] for i in range(n))%prime)
    coefficients=[1]
    for k in range(1,n+1):
        coefficients.append(-sum(coefficients[k-j]*traces[j-1] for j in range(1,k+1))*pow(k,-1,prime)%prime)
    # Separate Horner evaluation of the characteristic polynomial on its matrix.
    result=[[0]*n for _ in range(n)]
    for co in coefficients:
        result=multiply(result,matrix,prime)
        for i in range(n):result[i][i]=(result[i][i]+co)%prime
    return coefficients,traces,not any(any(row) for row in result)

def witnesses(name):
    row=E.family(name);par=row['parameter'];s=E.s;U,V=E.G.U,E.G.V
    prime=103 if name=='L' else 109
    f=S.Poly(row['f'].subs({par:1,E.A:0})/-9,E.M)
    assert f.LC()==1 and f.degree()==2
    fm=S.Poly.from_dict({mon:modular(co,prime) for mon,co in f.terms()},E.M,modulus=prime)
    roots=[m for m in range(prime) if int(fm.eval(m))%prime==0]
    assert len(roots)==2 and S.gcd(fm,fm.diff()).degree()==0
    h=S.Poly(S.expand(row['key'].subs(E.W,row['G'])).subs(par,1),U,V)
    out=[]
    for m in roots:
        sub={par:1,E.A:0,E.M:m}
        def univariate(expr):
            return S.Poly.from_dict({mon:modular(co,prime) for mon,co in S.Poly(expr.subs(sub),s).terms()},s,modulus=prime)
        P,Q=univariate(row['P']),univariate(row['Q'])
        def evaluate(poly):
            value=S.Poly(0,s,modulus=prime)
            for (i,j),co in poly.terms():value+=modular(co,prime)*P**i*Q**j
            return value
        hb=evaluate(h);critical=hb.diff().monic()
        gu,gv=evaluate(h.diff(U)),evaluate(h.diff(V))
        # Full multiplication matrix in the monic derivative quotient.
        columns=[]
        for j in range(15):
            reduced=(hb*S.Poly(s**j,s,modulus=prime)).rem(critical)
            columns.append([int(reduced.nth(i))%prime for i in range(15)])
        matrix=[list(col) for col in zip(*columns)]
        coefficients,traces,ch=charpoly(matrix,prime)
        value=S.Symbol('value');cp=S.Poly.from_list(coefficients,value,modulus=prime)
        checks={
            'intersection_degree16':hb.degree()==16,
            'critical_degree15':critical.degree()==15,
            'critical_points_simple':S.gcd(critical,critical.diff()).degree()==0,
            'critical_values_nonzero':S.gcd(critical,hb).degree()==0 and int(cp.nth(0))%prime!=0,
            'critical_values_pairwise_distinct':S.gcd(cp,cp.diff()).degree()==0,
            'avoids_target_critical_points':S.gcd(S.gcd(critical,gu),gv).degree()==0,
            'critical_value_Cayley_Hamilton':ch,
        }
        out.append(dict(M_mod_prime=m,hB_coefficients=[int(v)%prime for v in hb.all_coeffs()],
            critical_value_charpoly=coefficients,newton_traces=traces,checks=checks))
    return dict(family=name,parameter_value=1,A=0,coefficient_equation=str(f.as_expr()),
        prime=prime,split_roots=roots,coefficient_choices=out)
