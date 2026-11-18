#!/usr/bin/env python3
"""Exact coupled-ODE reductions for mixed reciprocal paired pencils.

Run .venv/bin/python research/arithmetic/verify_mixed_reciprocal_coupled_ode.py
Optional --write-receipt and --write-degree26-system.
The written theorem supplies all-degree equivalence and open conditions.
"""
from __future__ import annotations
import argparse,hashlib,json,platform
from pathlib import Path
import sympy as sp

w,v,r=sp.symbols("w v r")

def chebyshev(maximum,variable=w):
    C=[sp.Integer(1),variable-1]
    D=[sp.Integer(1),variable+1]
    for j in range(1,maximum):
        C.append(sp.expand(variable*C[-1]-C[-2]))
        D.append(sp.expand(variable*D[-1]-D[-2]))
    return C[:maximum+1],D[:maximum+1]

def from_critical(h,s):
    """s(w) monic degree h -> reciprocal critical polynomial and F,X,Y."""
    L=4*h+1;d=2*h
    s=sp.Poly(s,w)
    assert s.degree()==h and s.LC()==1
    coeff=[sp.Integer(0)]*(d+1)
    for k in range(h+1):
        for ell in range(k+1):
            coeff[h+k-2*ell]+=s.nth(k)*sp.binomial(k,ell)
    coeff=[sp.expand(x)for x in coeff]
    F=[sp.expand(L*((coeff[j-1]if j else 0)-
            r*(coeff[j]if j<=d else 0))/sp.Integer(2*j-1))
       for j in range(d+2)]
    C,D=chebyshev(d)
    X=sp.expand(F[0]+sum(F[j]*C[j-1]for j in range(1,d+2)))
    Y=sp.expand(-F[0]+sum(F[j]*D[j-1]for j in range(1,d+2)))
    return coeff,F,X,Y

def operator_inverse(rhs,sign=1,scale=2):
    shift=2*sign
    shifted=sp.Poly(sp.expand(rhs.subs(w,w-shift)),w)
    return sp.expand(sum(shifted.nth(j)*(w+shift)**j/sp.Integer(scale*j+1)
                         for j in range(shifted.degree()+1)))

def primitive(expr,variables):
    expr=sp.cancel(expr)
    num,den=sp.fraction(expr)
    assert not(set(variables)&den.free_symbols),"variable denominator in polynomial certificate"
    _,poly=sp.Poly(num,*variables,domain=sp.QQ).primitive()
    return poly

def degree26_system():
    """Complete square-lambda sign split 2+10 on its stated open set."""
    h=6;q=5;L=25
    a,b,c=sp.symbols("a b c")
    V=w*w+a*w+b
    C,_=chebyshev(h)
    def T(Q):
        return sp.Poly(sp.expand(4*(w-2)*V*sp.diff(Q,w)
               +(V+2*(w-2)*sp.diff(V,w))*Q),w)
    target=sp.Poly(sp.expand(L*(C[h]-r*C[h-1])*(w-c)),w)
    Q=w**q;error=T(Q)-target
    denominators=[]
    for j in reversed(range(q)):
        divisor=4*j+5;denominators.append(divisor)
        correction=sp.expand(-error.nth(j+2)/sp.Integer(divisor))
        Q+=correction*w**j
        error+=T(correction*w**j)
    assert error.degree()<=1
    Q=sp.expand(Q)
    critical,F,X,Y=from_critical(h,(w-c)*Q)
    x_at_c=sp.Poly(X.subs(w,c),c)
    quotient,remainder=x_at_c.div(sp.Poly(c+2,c))
    assert remainder.is_zero
    at_zero=sp.cancel(sum(F[j]*r**j for j in range(len(F)))/(2*r))
    exprs=[error.nth(0),error.nth(1),quotient.as_expr(),at_zero]
    variables=(a,b,c,r)
    equations=[primitive(expr,variables)for expr in exprs]
    stats=[dict(total_degree=int(poly.total_degree()),
                term_count=len(poly.terms()),
                individual_degrees=[int(poly.degree(z))for z in variables])
           for poly in equations]
    return dict(schema_version=1,h=6,degree=26,sign_degrees=[2,10],
       variables=list(map(str,variables)),
       V=str(V),Q=str(Q),equations=[str(poly.as_expr())for poly in equations],
       top_down_divisors=denominators,equation_statistics=stats,
       open_conditions=[
          "r!=0; c!=2,-2; Q(c)!=0; Q is squarefree; Q(2)Q(-2)!=0",
          "s(r+1/r)!=0, where s=(w-c)Q",
          "Y(c)!=0 and gcd(Q,X)=1",
          "The h normalized squared critical values in K[w]/((w-c)Q) are distinct and nonzero."
       ],
       equivalence="All rational nondegenerate square-lambda mixed models with sign degrees2+10 correspond to rational solutions of these four equations on the written open set; rho=1/r and N(u)=rho^13 F(u/rho).",
       finiteness_claim=False,rational_solution_found=False)

E18_COEFFICIENTS=[
20480,-5529600,435456000,-11695959360,88617286656,-4418525873664,
94720761695376,-886766024513376,3726746216490384,-2445510775053888,
-38911013766339552,93176136148671072,173618172500962707,-435934454220836550,
33480118841707629,413138775247514316,111327035757469245,5813728713056250,
-2906864356528125]

def degree10_elimination():
    """Three Sylvester-resultants give the global mixed-sign obstruction."""
    c,d,S,t,delta=sp.symbols("c d S t delta")
    _,F,X,Y=from_critical(2,(w-c)*(w-d))
    plus=sp.cancel(sp.Rational(35,2)*X.subs(w,c)/(c+2))
    minus=sp.cancel(-sp.Rational(35,2)*Y.subs(w,d)/(d-2))
    zero=sp.cancel(-35*sum(F[j]*r**j for j in range(6))/(2*r))
    expected_zero=21*c*d*r*r-9*(c+d)*r**3-105*(c+d)*r+5*r**4+42*r*r-315
    assert sp.expand(zero-expected_zero)==0
    A=-5*S**3+24*t*S+9*r*S*S-60*r*t-12*S+240*r
    B=-5*S*S+14*t+9*r*S-42
    C=15*S*S-84*t+36*r*S+112
    D=15*(S-6*r)
    substitution={c:(S+delta)/2,d:(S-delta)/2}
    torus=(S*S-delta*delta)/4
    assert sp.expand((plus-minus).subs(substitution)-(A+D*delta).subs(t,torus))==0
    assert sp.expand((plus+minus).subs(substitution)-(C+B*delta).subs(t,torus))==0
    t_from_zero=(S*r*(9*r*r+105)-5*r**4-42*r*r+315)/(21*r*r)
    def numerator(expr):
        expr=sp.cancel(expr.subs(t,t_from_zero))
        num,_=sp.fraction(expr)
        return sp.Poly(num,S,r,domain=sp.QQ).primitive()[1].as_expr()
    equations=[numerator(A*B-D*C),
               numerator(A*A-D*D*(S*S-4*t)),
               numerator(C*C-B*B*(S*S-4*t))]
    resultants=[sp.Poly(sp.resultant(equations[i],equations[j],S),r,domain=sp.QQ).primitive()[1]
                for i,j in[(0,1),(0,2),(1,2)]]
    g=sp.gcd(resultants[0],sp.gcd(resultants[1],resultants[2]))
    x=sp.symbols("x")
    eliminant=sp.Poly.from_list(E18_COEFFICIENTS,x,domain=sp.ZZ)
    expected=sp.Poly(r**24*eliminant.as_expr().subs(x,r*r),r)
    assert g.monic()==expected.monic()
    values=[int(eliminant.eval(j))%11 for j in range(11)]
    assert values==[7,8,4,3,2,2,10,7,6,3,4]
    assert int(eliminant.LC())%11==9 and all(values)
    cubic=sp.Poly(4*r**3-6*r**2-9*r-15,r)
    cubic_values=[int(cubic.eval(j))%7 for j in range(7)]
    assert cubic_values==[6,2,3,5,4,3,5] and int(cubic.LC())%7!=0
    gb=sp.groebner([plus,minus,zero],c,d,r,order="lex",domain=sp.QQ)
    assert gb.is_zero_dimensional and len(gb.polys)==3
    assert sp.expand(20480*gb.polys[-1].as_expr()-eliminant.as_expr().subs(x,r*r))==0
    gb_polynomials=[str(poly.as_expr())for poly in gb.polys]
    return dict(groebner_order="lex(c,d,r)",
        groebner_basis_polynomials=gb_polynomials,
        groebner_zero_dimensional=True,
        groebner_eliminant_degree=36,
        resultant_degrees=[int(x.degree())for x in resultants],
        resultant_gcd="constant * r^24 * E18(r^2)",
        eliminant_coefficients_descending=E18_COEFFICIENTS,
        mod11_values=values,mod11_infinity_value=9,
        uniform_case_cubic_mod7_values=cubic_values,
        intermediate_equations=[str(e)for e in equations],
        covers_square_and_nonsquare_lambda=True,
        conclusion="The previously known numerator is the unique rational normalized nondegenerate reciprocal degree10 model.")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt",type=Path)
    parser.add_argument("--write-degree26-system",type=Path)
    args=parser.parse_args()
    assert sp.__version__=="1.14.0"
    controls=[]
    for h in(1,2,3,6):
        coeff=sp.symbols("s0:"+str(h))
        s=w**h+sum(coeff[j]*w**j for j in range(h))
        critical,F,X,Y=from_critical(h,s)
        L=4*h+1;d=2*h
        Sv=sum(critical[j]*v**j for j in range(d+1))
        Fv=sum(F[j]*v**j for j in range(d+2))
        C,D=chebyshev(d)
        assert sp.Poly(Fv,v).LC()==1
        assert sp.expand(2*v*sp.diff(Fv,v)-Fv-L*(v-r)*Sv)==0
        assert sp.expand(v**d*Sv.subs(v,1/v)-Sv)==0
        assert sp.expand(Fv-((v+1)*X+(v-1)*Y).subs(w,v+1/v)/2)==0
        ax=L*(D[h]-r*D[h-1])*s
        by=L*(C[h]-r*C[h-1])*s
        assert sp.expand(X+2*(w+2)*sp.diff(X,w)-ax)==0
        assert sp.expand(Y+2*(w-2)*sp.diff(Y,w)-by)==0
        assert sp.expand(X-operator_inverse(ax,1))==0
        assert sp.expand(Y-operator_inverse(by,-1))==0
        norm=(w+2)*X*X-(w-2)*Y*Y
        assert sp.expand(norm.subs(w,v+1/v)-4*Fv*Fv.subs(v,1/v))==0
        controls.append(h)

    uniform_controls=[]
    for h in(2,6):
        C,D=chebyshev(h)
        scalar=operator_inverse((4*h+1)*(D[h]-r*D[h-1]),1,scale=4)
        critical,F,X,Y=from_critical(h,scalar)
        assert sp.expand(X-scalar**2)==0
        expected=[sp.Integer(1),-((4*h+1)*r+4*h-1)/sp.Integer(4*h-3)]
        for j in range(2,2*h):
            expected.append(sp.cancel(-(4*j-4*h-3)*expected[-1]/sp.Integer(4*j-4*h-1)))
        expected.append(sp.Integer(1))
        assert all(sp.expand(a-b)==0 for a,b in zip(critical,expected))
        uniform_controls.append(h)
        if h==2:
            Fr=sp.factor(sum(F[j]*r**j for j in range(6)))
            assert sp.expand(Fr-2*r*(2*r-15)*(4*r**3-6*r**2-9*r-15)/25)==0
            rho0=sp.Rational(2,15);r0=1/rho0
            N0=[sp.cancel(rho0**(5-j)*F[j].subs(r,r0))for j in range(6)]
            assert N0==[sp.Rational(16,5625),sp.Rational(1804,5625),
                -sp.Rational(14006,5625),sp.Rational(3129,625),
                -sp.Rational(96,25),sp.Integer(1)]

    # Generic nonsquare descent in the degree26 6+6 case.
    rho,lam=sp.symbols("rho lam")
    aa=sp.symbols("a0:3");bb=sp.symbols("b0:3")
    P=w**3+sum((aa[j]+bb[j]*rho)*w**j for j in range(3))
    Q=sp.expand(-P.subs({rho:-rho,w:-w},simultaneous=True))
    _,F,X,Y=from_critical(6,P*Q)
    domain=sp.QQ.frac_field(lam,*aa,*bb)
    modulus=sp.Poly(rho*rho-lam,rho,domain=domain)
    def quadratic(expr):
        return sp.Poly(sp.cancel(expr),rho,domain=domain).rem(modulus).as_expr()
    s=sp.expand(P*Q)
    assert all(quadratic(sp.Poly(s.subs(rho,-rho)-s.subs(w,-w),w).nth(j))==0
               for j in range(7))
    rational_N=[quadratic(rho**(13-j)*F[j].subs(r,rho/lam))for j in range(14)]
    assert all(not x.has(rho)for x in rational_N)
    Fr=quadratic(sum(F[j].subs(r,rho/lam)*(rho/lam)**j for j in range(14)))
    assert sp.Poly(Fr,rho,domain=domain).nth(0)==0
    assert sp.cancel(rational_N[0]-25*lam**6)==0

    small=degree10_elimination()
    system=degree26_system()
    if args.write_degree26_system:
        args.write_degree26_system.write_text(json.dumps(system,indent=2)+"\n")
    dependencies=["PAIRED_CRITICAL_RECIPROCAL_OBSTRUCTION.md",
                  "RATIONAL_DEGREE_TEN_ALTERNATING_PENCIL.md"]
    here=Path(__file__).resolve().parent
    report=dict(schema_version=1,result="PASS",kind="exact_mixed_reciprocal_reduction",
       python=platform.python_version(),sympy=sp.__version__,
       checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       dependency_sha256={p:hashlib.sha256((here/p).read_bytes()).hexdigest()for p in dependencies},
       symbolic_uniform_h_controls=controls,
       uniform_sign_recovery_controls=uniform_controls,
       quadratic_descent="Generic h6 cubic factor over Q(rho), rho^2=lam; reconstructed N has rational coefficient functions and F(r) is anti-invariant.",
       degree10_global_classification=small,
       degree26_2_10_equation_statistics=system["equation_statistics"],
       degree26_2_10_top_down_divisors=system["top_down_divisors"],
       degree26_solution_found=False,degree26_global_nonexistence_proved=False,
       scope="Written proofs establish uniform mixed-sign equivalence, quadratic descent and four-variable2+10 reduction. Exact resultants/projective certificates establish global degree10 reciprocal uniqueness. No claim of finite Q-points, global rational solution or nonexistence is made in degree26.",
       independent_replay=False,external_review=False,formal_verification=False)
    text=json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:args.write_receipt.write_text(text)
    print(text,end="")

if __name__=="__main__":
    main()
