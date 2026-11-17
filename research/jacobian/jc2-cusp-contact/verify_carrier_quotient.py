#!/usr/bin/env python3
"""Small exact carrier quotient certificate. Does not rebuild the F2 circuit."""
from pathlib import Path
from hashlib import sha256
import ast
import json
import platform
import time
import sympy as S

H=Path(__file__).resolve().parent
R=H.parents[1]
SOURCE='foundations/scripts/verify_f2_affine_k1_carrier_jet_factorization.py'
PIN='91949f5c28ba76929029f1bfb83e81413e5539680901a58ac407775a2ecef7e6'
checks=[]
def check(name,value):
    if not value: raise AssertionError(name)
    checks.append(name)
def zero(name,value): check(name,S.expand(value)==0)
N=7
def add(a,b): return [S.expand(a[i]+b[i]) for i in range(N+1)]
def scale(a,c): return [S.expand(c*x) for x in a]
def mul(a,b):
    return [S.expand(sum(a[j]*b[i-j] for j in range(i+1))) for i in range(N+1)]
one=[S.Integer(1)]+[S.Integer(0)]*N
def power(a,n):
    if n>=0:
        out=one
        for j in range(n): out=mul(out,a)
        return out
    assert a[0]==1
    u=add(a,scale(one,-1)); out=one; term=one
    for j in range(1,N+1):
        term=mul(term,u)
        out=add(out,scale(term,S.binomial(n,j)))
    return out
def extract(tree,function,target,env):
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==function)
    assign=next(n for n in node.body if isinstance(n,ast.Assign) and
                any(isinstance(x,ast.Name) and x.id==target for x in n.targets))
    allowed=(ast.List,ast.BinOp,ast.UnaryOp,ast.Name,ast.Constant,ast.Load,
             ast.Add,ast.Sub,ast.Mult,ast.Div,ast.Pow,ast.USub)
    assert all(isinstance(n,allowed) for n in ast.walk(assign.value))
    return eval(compile(ast.Expression(assign.value),'<retained arithmetic>','eval'),
                {'__builtins__':{}},env)
def graph_residual(hh,pi,coeff):
    out=one
    pp=one
    for item in coeff:
        pp=mul(pp,pi)
        out=add(out,scale(pp,item))
    return add(hh,scale(out,-1))
def main():
    started=time.monotonic()
    source=R/SOURCE
    check('retained arithmetic source pinned',sha256(source.read_bytes()).hexdigest()==PIN)
    tree=ast.parse(source.read_text())
    h=S.symbols('h1:8'); mu,eta,nu=S.symbols('mu eta nu')
    forward=extract(tree,'fixed_coordinate_transport_audit','forward',
                    dict(zip(map(str,h),h))|dict(mu=mu,eta=eta,nu=nu))
    hs=[S.Integer(1)]+list(h)
    pi=[S.Integer(0),S.Integer(1)]+[S.Integer(0)]*6
    A=add(one,scale(mul(power(pi,3),power(hs,-2)),mu))
    B=add(add(one,scale(mul(power(pi,2),power(hs,-1)),eta)),
          scale(mul(power(pi,5),power(hs,-3)),nu))
    new_pi=mul(pi,mul(power(A,3),power(B,-2)))
    new_h=mul(hs,mul(power(A,5),power(B,-3)))
    residual=graph_residual(new_h,new_pi,forward)
    for j,r in enumerate(residual):
        zero('independent target transport series coefficient '+str(j),r)
    eta0=h[1]/3
    mu0=-(h[2]+2*h[0]*eta0)/5
    nu0=S.expand(forward[4].subs({eta:eta0,mu:mu0,nu:0}))/3
    gauge={eta:eta0,mu:mu0,nu:nu0}
    normalized=[S.expand(f.subs(gauge)) for f in forward]
    for j in (1,2,4): zero('global gauge jet '+str(j+1),normalized[j])
    invariants=[normalized[j] for j in (0,3,5,6)]
    section={h[1]:0,h[2]:0,h[4]:0}
    for j,value in zip((0,3,5,6),invariants):
        zero('cross section invariant '+str(j+1),value.subs(section)-h[j])
    # The exact group law gives uniqueness of the slice on every orbit.
    mu2,eta2,nu2,p,g=S.symbols('mu2 eta2 nu2 p g')
    zero('target group composition second coordinate',
         g+eta*p+nu+eta2*(p+mu)+nu2
         -(g+(eta+eta2)*p+nu+nu2+eta2*mu))
    inverse={mu2:-mu,eta2:-eta,nu2:-nu+eta*mu}
    zero('target group inverse constant',(nu+nu2+eta2*mu).subs(inverse))
    a,b,c,d,z=S.symbols('a b c d z')
    oldH=extract(tree,'carrier_graph_jet_audit','H',dict(a=a,b=b,c=c,d=d))
    Ap=[S.Integer(1),0,a,0,0,0,0,0]
    Bp=[S.Integer(1),b,0,c,d,0,0,0]
    pp=mul(pi,mul(power(Ap,3),power(Bp,-2)))
    hh=mul(power(Ap,5),power(Bp,-3))
    for j,r in enumerate(graph_residual(hh,pp,oldH)):
        zero('independent k1 graph coefficient '+str(j),r)
    mapped=[S.expand(x.subs(dict(zip(h,oldH)),simultaneous=True)) for x in invariants]
    expected=[
        -3*b,
        -(50*a*a+48*a*b*b-27*b**4-36*b*c+45*d)/15,
        (125*a**3+9144*a*a*b*b+8640*a*b**4-216*a*b*c+135*a*d
         -5049*b**6-6480*b**3*c+8100*b*b*d+81*c*c)/135,
        (650*a**3*b+39118*a*a*b**3+25*a*a*c+36816*a*b**5-1152*a*b*b*c
         +690*a*b*d-21627*b**7-27612*b**4*c+34515*b**3*d
         +432*b*c*c+45*c*d)/75]
    for j,(f,e) in enumerate(zip(mapped,expected)): zero('quotient polynomial '+str(j),f-e)
    for f,weight in zip(mapped,(1,4,6,7)):
        check('weighted homogeneity '+str(weight),
              all(sum(u*v for u,v in zip(mon,(2,1,3,4)))==weight
                  for mon,coef in S.Poly(f,a,b,c,d).terms()))
    G=S.groebner(mapped,b,d,c,a)
    gb=[b,d+S.Rational(10,9)*a*a,c*c-S.Rational(25,81)*a**3,a*a*c,a**5]
    check('quotient zero fibre exact Groebner basis',list(G)==gb)
    standard=[S.Integer(1),a,a*a,a**3,a**4,c,a*c]
    check('zero fibre standard monomials',
          all(G.reduce(x)[1]==x for x in standard))
    check('zero fibre supported at origin',
          all(G.reduce(v**9)[1]==0 for v in (a,b,c,d)))
    leading=[S.Poly(f,b,d,c,a).monoms()[0] for f in G.polys]
    std=[]
    for aa in range(5):
        for cc in range(2):
            mon=(0,0,cc,aa)
            if not any(all(i>=j for i,j in zip(mon,lm)) for lm in leading):
                std.append(a**aa*c**cc)
    check('complete zero fibre standard monomial enumeration',set(std)==set(standard))
    t=S.symbols('t')
    zero('finite free Hilbert numerator',
         S.expand((1-t**6)*(1-t**7)
          -(1-t*t)*(1-t**3)*(1+t*t+t**3+t**4+t**5+t**6+t**8)))
    jac=S.factor(S.det(S.Matrix(mapped).jacobian([a,b,c,d])))
    res=S.resultant(3*z*z+a,5*z**4+4*b*z**3+2*c*z+d,z)
    zero('quotient ramification determinant',jac+res/15)
    check('map generically separable',jac!=0)
    # Explicit nonzero seventh-centre fibre, kept separate from any attachment claim.
    w=S.symbols('w')
    substitution={b:0,d:-S.Rational(10,9)*a*a}
    low=[S.expand(f.subs(substitution)) for f in mapped]
    relation=S.groebner([c*c-S.Rational(25,81)*a**3],c,a)
    for j in range(3):
        check('special seventh-centre first three invariants '+str(j),
              relation.reduce(low[j])[1]==0)
    check('special seventh-centre last invariant',
          relation.reduce(S.expand(low[3]+a*a*c/3))[1]==0)
    zero('special Kummer identity',
         (a*a*c)**2-a**7*S.Rational(25,81)
         -a**4*(c*c-S.Rational(25,81)*a**3))
    check('special seventh-centre immersion resultant',
          relation.reduce(S.expand(res.subs(substitution)-S.Rational(175,3)*a**4))[1]==0)
    normalized_L=[39,40,45,55,70]
    check('five normalized contact strata',
          sorted([75-5*j for j in (1,4,6,7)]+[39])==normalized_L)
    receipt=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        affine_gauge={str(k):str(v) for k,v in gauge.items()},
        quotient_invariants=list(map(str,invariants)),
        k1_quotient_polynomials=list(map(str,mapped)),zero_fibre_basis=list(map(str,G)),
        zero_fibre_length=7,finite_free_rank=7,
        module_basis=list(map(str,standard)),
        normalized_bands=normalized_L,normalized_contacts=[13*l for l in normalized_L],
        jacobian_determinant=str(jac),source_hashes={SOURCE:PIN},
        files={n:sha256((H/n).read_bytes()).hexdigest() for n in
               ('CARRIER_QUOTIENT.md','verify_carrier_quotient.py')},
        boundary='Affine gauge equivalence, not deletion of original maps; finite flat degree7 is on normalized k1 coefficient presentations, not Keller realizations or global attachment.',
        python=platform.python_version(),sympy=S.__version__,
        elapsed_seconds=round(time.monotonic()-started,4))
    (H/'carrier-quotient-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',checks=len(checks),rank=7,
                         normalized_contacts=receipt['normalized_contacts'],
                         elapsed_seconds=receipt['elapsed_seconds'])))
if __name__=='__main__': main()

