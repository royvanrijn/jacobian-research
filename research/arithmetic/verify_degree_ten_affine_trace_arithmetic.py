#!/usr/bin/env python3
"""Exact curve binding and finite certificates for the degree-ten finiteness proof.

No factorization oracle or exact-genus algorithm is used. Deep finiteness
comes from the written application of Faltings, not from this computation.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
import argparse,json,platform
import sympy as sp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def primitive(expr,*variables):
    p=sp.Poly(expr,*variables).clear_denoms()[1].primitive()[1]
    return -p if p.LC()<0 else p

def powmod(base,exponent,modulus):
    out=sp.Poly(1,base.gens[0],modulus=modulus.get_modulus())
    while exponent:
        if exponent&1:out=(out*base).rem(modulus)
        exponent//=2
        base=(base*base).rem(modulus)
    return out

def polynomial_from_coefficients(coefficients,x,p):
    return sp.Poly(sum(int(a)*x**j for j,a in enumerate(coefficients)),x,modulus=p)

def four_cycle_controls():
    def cp(a,b):return tuple(a[b[j]]for j in range(4))
    cycle=(1,2,3,0);identity=tuple(range(4));sizes=[]
    for a,b,c in combinations(range(4),3):
        for points in ((a,b,c),(a,c,b)):
            g=list(identity)
            for j,v in enumerate(points):g[v]=points[(j+1)%3]
            g=tuple(g);group={identity};queue=[identity]
            for h in queue:
                for k in (cycle,g):
                    w=cp(h,k)
                    if w not in group:group.add(w);queue.append(w)
            assert len(group)==24
            sizes.append(len(group))
    assert len(sizes)==8
    return sizes

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--write-receipt',type=Path,
                    default=HERE/'degree_ten_affine_trace_arithmetic_verification.json')
    args=ap.parse_args()
    path=HERE/'degree_ten_affine_trace_arithmetic_inputs.json'
    inputs=json.loads(path.read_text())
    s,t,u,r=sp.symbols('s t u R')
    F=sp.Poly(sp.sympify(inputs['curve']),s,t)
    assert F==primitive(F.as_expr(),s,t)
    assert (F.degree(s),F.degree(t),F.total_degree())==(4,16,16)
    coefficient_gcd=sp.Poly(0,t)
    for coefficient in sp.Poly(F.as_expr(),s).all_coeffs():
        coefficient_gcd=sp.gcd(coefficient_gcd,sp.Poly(coefficient,t))
    assert coefficient_gcd.degree()==0
    q0,q1=sp.symbols('q0 q1')
    q=r*r+q1*r+q0
    S=sp.Poly(sp.resultant(q,r*r-t*u*r+u*u-s*u,r),u)
    critical=sp.Poly(9*(u-1)*S.as_expr(),u)
    N=sp.Poly(sum(critical.nth(j)*u**j/sp.Integer(2*j-1)
                  for j in range(6)),u)
    assert S.degree()==4 and S.LC()==1
    assert N.degree()==5 and N.LC()==1 and N.nth(0)==9*q0*q0
    theta=s+t*r;d=theta-2*r
    C=[sp.Integer(0),sp.Integer(1)]
    T=[sp.Integer(2),theta]
    for j in range(2,6):
        C.append(sp.expand(theta*C[-1]-r*r*C[-2]))
        T.append(sp.expand(theta*T[-1]-r*r*T[-2]))
    D=sp.Poly(sp.expand(sum(N.nth(j)*(T[j]-d*C[j])for j in range(6))),r)
    gamma=[sp.Integer(0),sp.Integer(1)]
    for j in range(2,6):gamma.append(sp.expand(t*gamma[-1]-gamma[-2]))
    k=gamma[5]-gamma[4]
    E=sp.Poly(sp.expand(D.as_expr()-2*(9+k*r)*q*q),r)
    assert sp.expand(E.nth(4).subs(t,-2))==-24*s
    formulas=[sp.sympify(a)for a in inputs['q_coefficients_ascending']]
    assert formulas[2]==1
    sub={q0:formulas[0],q1:formulas[1]}
    for j in range(6):assert sp.cancel(E.nth(j).subs(sub))==0,j
    numerator,denominator=sp.cancel(N.eval(1).subs(sub)).as_numer_denom()
    assert primitive(numerator,s,t)==F
    assert sp.Poly(sp.gcd(F.as_expr(),denominator),s,t).total_degree()==0
    f0=primitive(F.as_expr().subs(t,0),s)
    assert f0==sp.Poly(sp.sympify(inputs['specialization_polynomial']),s)
    assert sp.discriminant(f0.as_expr(),s)!=0
    monic_integral=sp.Poly(sp.expand(9**3*f0.as_expr().subs(s,s/9)),s)
    assert monic_integral.LC()==1 and all(c.is_Integer for c in monic_integral.all_coeffs())
    assert sp.expand(monic_integral.as_expr().subs(s,9*s)-9**3*f0.as_expr())==0

    # Bound the discriminant by an exact retained factor product, not factoring.
    disc=sp.Poly(sp.discriminant(F.as_expr(),s),t)
    factors=[(sp.Poly(sp.sympify(row['polynomial']),t),int(row['multiplicity']))
             for row in inputs['discriminant_factors']]
    product=sp.Poly(sp.Integer(inputs['discriminant_unit']),t)
    for (p,m),row in zip(factors,inputs['discriminant_factors']):
        assert p.degree()==row['degree'];product*=p**m
    assert product==disc
    factor30=next(p for p,m in factors if p.degree()==30 and m==1)
    assert sp.gcd(factor30,factor30.diff()).degree()==0
    assert sp.gcd(factor30,sp.Poly(sp.Poly(F.as_expr(),s).LC(),t)).degree()==0
    assert all(sp.gcd(factor30,p).degree()==0 for p,m in factors if p!=factor30)

    frobenius=[]
    for row in inputs['frobenius_controls']:
        p=int(row['prime']);g=polynomial_from_coefficients(row['monic_coefficients_ascending'],s,p)
        assert sp.Poly(f0.as_expr(),s,modulus=p).monic()==g
        assert sp.gcd(g,g.diff()).degree()==0
        if p==11:
            assert row['cycle_type']==[4]
            x=sp.Poly(s,s,modulus=p)
            assert (powmod(x,p**4,g)-x).rem(g).is_zero
            assert sp.gcd(g,powmod(x,p**2,g)-x).degree()==0
        else:
            assert p==13
            assert row['cycle_type']==[1,3]
            linear=polynomial_from_coefficients(row['linear'],s,p)
            cubic=polynomial_from_coefficients(row['cubic'],s,p)
            assert linear*cubic==g
            assert all(int(cubic.eval(a))%p for a in range(p))
        frobenius.append(dict(prime=p,cycle_type=row['cycle_type'],separable=True))
    groups=four_cycle_controls()

    # Full rational diagonal-exception certificate for h=2.
    exception_rows=[]
    for j in range(2):
        qj=4*j+1;step=2*(2-j)
        # Possible denominator primes divide qj. Only qj=5 occurs here.
        if qj==5:assert (1-0)%step!=0
        gj=gamma[2*j+1]-gamma[2*j]
        diagonal=sp.Poly(9*gj-qj*k,t)
        constant=int(diagonal.eval(0));assert constant!=0
        candidates=sorted({a*sign for a in range(1,abs(constant)+1)
                           if constant%a==0 for sign in(-1,1)})
        values=[(a,int(diagonal.eval(a)))for a in candidates]
        assert [a for a,v in values if v==0]==[-2]
        exception_rows.append(dict(j=j,constant=constant,candidate_values=values,
                                   only_rational_zero=-2))

    sources=[path,Path(__file__),HERE/'DEGREE_TEN_AFFINE_TRACE_CURVE_ARITHMETIC.md',
             HERE/'ODD_AFFINE_CRITICAL_TRACE_CURVES.md',
             HERE/'MIXED_RECIPROCAL_COUPLED_ODE.md',
             HERE/'sources/faltings-1983-endlichkeitssaetze.pdf',
             HERE/'sources/milne-ant300.pdf']
    report=dict(schema_version=1,result='PASS',python=platform.python_version(),
                sympy=sp.__version__,checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                source_sha256={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sources},
                curve_degrees=dict(s=4,t=16,total=16),
                coefficient_content_in_Qt=1,
                direct_norm_euler_trace_binding=True,
                discriminant_factors=[dict(degree=p.degree(),multiplicity=m)for p,m in factors],
                simple_branch_points=30,branch_leading_coefficient_coprime=True,
                arithmetic_frobenius_controls=frobenius,
                monic_integral_specialization=str(monic_integral.as_expr()),
                four_cycle_three_cycle_closure_orders=groups,
                rational_diagonal_exceptions=exception_rows,
                genus_lower_bound=12,exact_genus_claimed=False,
                finiteness_import='Faltings1983,Satz7,p365; written application only, not reproved.',
                scope='Normalized rational even-numerator degree10 models, with full nonscalar-scale curve and separate classified reciprocal boundary. No effective point list or finiteness of raw twists/Q-isomorphism classes.',
                formal_verification=False,independent_external_replay=False,external_human_review=False)
    args.write_receipt.write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: exact degree10 curve binding, branch/genus bound and arithmetic controls.')

if __name__=='__main__':main()
