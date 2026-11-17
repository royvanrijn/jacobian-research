#!/usr/bin/env python3
"""Exact bounded controls, not an asymptotic slope or general gap solver."""
from pathlib import Path
import hashlib,json
import sympy as S

HERE=Path(__file__).resolve().parent
x,y,z,X,Y,Z,t,a,b,c,d,f,u,eps=S.symbols('x y z X Y Z t a b c d f u eps')
xs=(x,y,z);zs=(X,Y,Z)


def action(symbol,poly):
    assert not S.sympify(symbol).has(S.Float) and not S.sympify(poly).has(S.Float)
    poly=S.expand(poly);assert S.Poly(poly,*xs).total_degree()<=64
    result=0
    for powers,coefficient in S.Poly(S.expand(symbol),*zs).terms():
        term=poly
        for v,k in zip(xs,powers):term=S.diff(term,v,k)
        result+=coefficient*term
    return S.expand(result)


def iterate(symbol,poly,k):
    assert isinstance(k,int) and 0<=k<=16
    for _ in range(k):poly=action(symbol,poly)
    return poly


def height(symbol,phase,m):
    if m==0:return 0
    poly=S.expand(phase**m);assert poly!=0
    for k in range(1,17):
        poly=action(symbol,poly)
        if poly==0:return k-1
    raise AssertionError('height exceeds control bound')


def power_basis(poly,H):
    answer={};k=0
    while poly!=0:
        quotient,remainder=S.div(poly,H,t)
        for (i,),coefficient in S.Poly(remainder,t).terms():
            if coefficient!=0:answer[i,k]=coefficient
        poly=quotient;k+=1
    return answer


def dual_reduce(poly):
    poly=S.Poly(S.expand(poly),eps)
    return S.expand(poly.nth(0)+eps*poly.nth(1))


def main():
    basis=[]
    for H in [t**2+a*t+b,t**3-t+1]:
        e=S.degree(H,t)
        for N in range(9):
            terms=power_basis(t**N,H)
            assert S.expand(sum(v*t**i*H**k for (i,k),v in terms.items())-t**N)==0
            assert all(0<=i<e and e*k+i<=N for i,k in terms)
            basis.append(dict(degree=int(e),N=N,triangular=True))
    models=[('unbalanced_safe',X**2,x*y**3,t**2+t+1,4),
            ('formal_graph_jet',X-Y**2-Y**6-Y**24,z**3,t**2-t,3),
            ('nonpure_boundary',X,x,t**3-t,1)]
    heights=[]
    for tag,symbol,phase,H,D in models:
        e=int(S.degree(H,t));R=H.subs(t,phase)
        original=[height(symbol,phase,j) for j in range(10)]
        recomposed=[height(symbol,R,k) for k in range(4)]
        for m in range(1,4):
            hp=recomposed[m];g=height(S.expand(symbol**e),R,m)
            assert g==hp//e and hp<=max(original[:e*m+1])
            heights.append(dict(model=tag,m=m,original_height=original[m],
                                recomposed_T_height=hp,paired_height=g,
                                floor_equal=True,forward_height_bound=True))
        for N in range(1,7):
            assert max(original[:N+1])<=max(recomposed[:N//e+1])+D*(e-1)
    vanish=[]
    phase=x*y**3;H=t**2+t+1;R=H.subs(t,phase)
    for Q in [x,y,z]:
        assert iterate(X**2,Q*phase**2,2)==0
        assert iterate(X**4,Q*R**4,4)==0
        vanish.append(dict(marker=str(Q),old_index=2,new_index=4,zero=True))
    prefix=[]
    symbol=a*X+b*Y**2;phase=c*x+d*y**2+f
    for m in range(1,4):
        lhs=iterate(symbol,S.expand((phase+u)**m),m)-iterate(symbol,S.expand(phase**m),m)
        rhs=sum(S.binomial(m,j)*u**(m-j)*iterate(symbol,iterate(symbol,S.expand(phase**j),j),m-j)
                for j in range(1,m))
        assert S.expand(lhs-rhs)==0
        prefix.append(m)
    # Nonreduced model: original pure/cutoff0 uses eps^2=0 genuinely.
    symbol=eps*X**2;phase=eps*x**2;R=phase**2+phase+1
    for Q in [S.Integer(1),x,y**2]:
        for m in [1,2]:
            assert dual_reduce(iterate(symbol,Q*phase**m,m))==0
            assert dual_reduce(iterate(S.expand(symbol**2),Q*R**m,m))==0
    Nq,Di,delta,ee,CC=S.symbols('Nq Di delta ee CC')
    B=(1+ee*Di*CC)*(1+Di*(ee-1))
    assert S.expand(B*delta-(1+ee*Di*CC)*(delta+Di*Nq)
                    -(1+ee*Di*CC)*Di*((ee-1)*delta-Nq))==0
    cutoff=[]
    for C,D,e in [(1,4,2),(3,2,5),(0,3,2),(2,3,1)]:
        forward=int(S.ceiling(S.Rational(C*D+1,e)))
        backward=(1+e*D*C)*(1+D*(e-1))
        assert e*forward>=C*D+1
        cutoff.append(dict(C=C,D=D,e=e,forward=forward,backward=backward))
    guards=0
    for symbol,poly in [(X,S.Float(1)),(X,x**65)]:
        try:action(symbol,poly)
        except AssertionError:guards+=1
        else:raise AssertionError('guard failed')
    try:iterate(X,x,17)
    except AssertionError:guards+=1
    else:raise AssertionError('guard failed')
    deps=['verify_formal_apolar_gaps.py','verify_operator_tail_safety.py']
    report=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        supporting_checker_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        sympy_version=S.__version__,triangular_basis_controls=basis,
        paired_height_controls=heights,forward_zero_controls=vanish,
        pure_prefix_identity_controls=prefix,square_zero_cutoff_zero_control=True,
        backward_budget_factorization=True,cutoff_arithmetic_controls=cutoff,
        guard_refusals=guards,exact_bounded_calculator_implemented=True,
        universal_transfer_proved_by_finite_controls=False,
        asymptotic_slope_identity_inferred_from_fitting=False,
        general_formal_gap_decided=False,
        original_purity_reflected_from_new_pure_tests=False,
        same_index_mixed_ideals_claimed_equal=False,
        nonunit_leading_reflection_claimed=False,
        general_ideal_or_gap_solver_implemented=False,
        full_exceptional_MZ_decided=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    target=HERE/'gvc_paired_recomposition_verification.json';assert not target.exists()
    target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(completed=True,basis=len(basis),heights=len(heights),
                          prefix=prefix,dual_number=True,guards=guards)))


if __name__=='__main__':main()
