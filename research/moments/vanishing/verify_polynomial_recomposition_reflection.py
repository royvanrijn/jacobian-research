"""Small exact controls for the written universal recomposition proof."""
from pathlib import Path
import hashlib,json
import sympy as s


def divided_difference(H,T,z):
    H=s.sympify(H)
    if H.has(s.Float):raise ValueError('Exact coefficients required')
    p=s.Poly(H,T)
    if p.degree()<1 or p.degree()>8:raise ValueError('Nonconstant degree<=8 control required')
    quotient,rem=s.div(H.subs(T,z)-H,z-T,T)
    assert rem==0
    assert s.expand((z-T)*quotient-H.subs(T,z)+H)==0
    return s.Poly(quotient,T)


def haar(expr,y,z):
    p=s.Poly(s.expand(expr),y,z)
    if p.total_degree()>30:raise ValueError('Control degree cap')
    return s.expand(sum(c*s.Rational((-1)**a,a+1) for (a,b),c in p.terms() if a==b))


def run():
    T,z,t,u=s.symbols('T z t u');y=s.symbols('y');F=(1+y)*(1+z+y*z);G=-y
    laws=[]
    for n in range(7):
        pure=haar(F**n,y,z);mixed=haar(G*F**n,y,z)
        assert pure==(1 if n==0 else 0)
        assert mixed==(0 if n==0 else s.Rational((-1)**(n-1),n+1))
        laws.append(dict(index=n,pure=str(pure),marked=str(mixed)))
    reconstruction=[]
    for H in [T,T**2,T**3+2*T,T**4-3*T*T+T,T**3+1,2*T*T-2]:
        q=divided_difference(H,T,z)
        reconstruction.append(dict(polynomial=str(H),degree=int(s.degree(H,T)),coefficients=[str(q.nth(j)) for j in range(s.degree(H,T))]))
        assert s.cancel(q.as_expr()/(H.subs(T,z)-H)-1/(z-T))==0
    # A normalized finite-evaluation functional lacks the time-zero rank-one premise.
    # P has two distinct spectral centres in this relaxed model, and H(P)=0 can be safe
    # without P-c being safe for any c. This guards the essential source hypothesis.
    H=T*(T-1);q=divided_difference(H,T,z)
    C=s.Rational(1,2)/z+s.Rational(1,2)/(z-1)
    rhs=sum(q.nth(j)*((1 if j==0 else s.Rational(1,2))/H.subs(T,z)) for j in range(2))
    assert s.cancel(C-rhs)==0
    # Exact CRT projection on (D-c)^N: use spectral variable u and verify residues.
    factors=[(u-2)**3,(u+1)**2,u**2];projectors=[]
    total=s.prod(factors)
    for i,f in enumerate(factors):
        other=s.cancel(total/f);a=s.invert(other,f);proj=s.rem(s.expand(a*other),total)
        assert all(s.rem(proj,g)==(1 if j==i else 0) for j,g in enumerate(factors))
        projectors.append(str(proj))
    assert s.rem(sum(s.sympify(p) for p in projectors)-1,total)==0
    # Laplace/exponential normalization: (z-c)^-k <-> t^(k-1)e^(ct)/(k-1)!.
    partial=[]
    for c,k in [(0,1),(0,3),(2,1),(-1,2),(2,3)]:
        rational=1/(z-c)**k;egf=t**(k-1)*s.exp(c*t)/s.factorial(k-1)
        laurent=s.series(rational.subs(z,1/u)/u,u,0,7).removeO().expand()
        exponential=s.series(egf,t,0,7).removeO().expand()
        assert all(s.expand(exponential).coeff(t,n)*s.factorial(n)==laurent.coeff(u,n) for n in range(7))
        partial.append(dict(centre=c,pole_order=k,exponential=str(egf)))
    # Long control's logarithmic Cauchy derivative has nonzero simple residues.
    deriv=-1/z**2+1/(z*(z+1))
    residues=[s.residue(deriv,z,c) for c in [0,-1]];assert residues==[1,-1]
    logseries=s.series(1-s.log(1+u)/u,u,0,7).removeO().expand()
    assert [logseries.coeff(u,n) for n in range(7)]==[s.sympify(v['marked']) for v in laws]
    # Finite nilpotence model, two coefficient directions; no universal exponent fit.
    a,b=s.symbols('a b');R=a*T+b*T*T
    nil=[]
    for L in [2,3,4]:
        for m in [L,L+1]:
            p=s.Poly(s.expand(R**m),a,b,T)
            assert all(alpha[0]+alpha[1]>=L for alpha,c in p.terms())
        nil.append(dict(ideal_power=L,vanishing_indices=[L,L+1]))
    # Finite jets do not prove rationality: a delayed perturbation matches a rational
    # prefix. The receipt deliberately certifies identities, not a rationality oracle.
    delayed=1/(1-u)+u**7*s.log(1-u)
    assert s.series(delayed-1/(1-u),u,0,8).removeO()==0
    refusals=0
    for action in [lambda:divided_difference(1,T,z),lambda:divided_difference(s.Float('0.5')*T,T,z),lambda:divided_difference(T**9,T,z)]:
        try:action()
        except ValueError:refusals+=1
        else:raise AssertionError('Expected refusal')
    return dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                divided_difference_controls=reconstruction,haar_controls=laws,crt_projectors=projectors,
                partial_fraction_controls=partial,logarithmic_derivative_residues=[str(v) for v in residues],nilpotence_controls=nil,
                rank_one_hypothesis_countercontrol=str(C),guard_refusals=refusals,supporting_checker_hashes={},
                exact_divided_difference_calculator_implemented=True,finite_marker_block_rule_implemented=True,
                universal_theorem_proved_by_finite_controls=False,same_single_marker_for_every_recomposition_claimed=False,
                eventually_every_index_nonzero_claimed=False,full_exceptional_MZ_decided=False,
                generic_nonrational_marker_solver_implemented=False,independent_replay=False,
                external_review=False,formal_verification=False,novelty_claimed=False,finite_prefix_rationality_inference_used=False)


if __name__=='__main__':
    receipt=run();output=Path(__file__).with_name('polynomial_recomposition_reflection_verification.json')
    if output.exists():assert json.loads(output.read_text())==receipt,'Preserve immutable receipt'
    else:output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(completed=True,checker_sha256=receipt['checker_sha256'],divided_differences=len(receipt['divided_difference_controls']),haar_indices=len(receipt['haar_controls']),partial_fractions=len(receipt['partial_fraction_controls']),guard_refusals=receipt['guard_refusals'])))
