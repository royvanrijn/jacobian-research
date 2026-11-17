#!/usr/bin/env python3
"""Small transcription controls for scalar time equations and actual deficits."""
from pathlib import Path
import hashlib,json,math
import sympy as S

HERE=Path(__file__).resolve().parent
t,s,z,n,x,y,X,Y,A,B=S.symbols('t s z n x y X Y A B')

def apply_time(op,f):return S.expand(sum(c*t**i*S.diff(f,t,j) for (i,j),c in op.items()))

def indicial(op,source=s):
    assert op and all(isinstance(i,int) and isinstance(j,int) and i>=0 and j>=0 for i,j in op)
    assert all(not S.sympify(c).has(S.Float) for c in op.values())
    # Only supplied finite coefficient models are accepted by this checker.
    order=min(min(a[0] for a,c in S.Poly(v,source).terms() if c) for v in op.values() if v!=0)
    leading={ij:S.expand(c).coeff(source,order) for ij,c in op.items() if S.expand(c).coeff(source,order)!=0}
    d=max(i-j for i,j in leading);M=max(i-j for (i,j),c in op.items() if c!=0)
    p=S.expand(sum(c*S.prod(n-k for k in range(j)) for (i,j),c in leading.items() if i-j==d))
    assert p!=0
    roots=S.solve(p,n);integers=[int(a) for a in roots if a.is_Integer and a>=0]
    bound=max([0]+integers);C=max(bound,M-d)
    return dict(source_order=order,leading_gain=d,global_gain=M,polynomial=str(p),B=bound,C=C)

def diff_action(lam,f):
    return S.expand(sum(c*S.diff(f,x,a[0],y,a[1]) for a,c in S.Poly(lam,X,Y).terms()))

def main():
    phi=s**2+s**6+s**24
    ops=[{(1,2):1,(0,1):1,(0,0):-phi},
         {(0,1):1,(1,1):-s,(0,0):-s**2},
         {(2,2):1,(1,1):-4},
         {(1,2):s**3,(0,1):s**3,(0,0):-s**3*phi},
         {(0,1):1,(0,0):-z*x,(1,0):-2*z}]
    reports=[indicial(op,s if i!=4 else z) for i,op in enumerate(ops)]
    assert [a['C'] for a in reports]==[1,1,5,1,2]
    assert reports[3]['source_order']==3
    multi={(1,2):(A+B)*s,(0,1):(A+B)*s,(1,0):-(A*A+2*B*B)*s*s}
    mr=indicial(multi);assert mr['source_order']==1 and mr['C']==2
    # Check the leading-degree formula above the supplied integer roots.
    leadchecks=0
    for op,rep in zip(ops,reports):
        q=z if op is ops[-1] else s
        lead={ij:S.expand(c).coeff(q,rep['source_order']) for ij,c in op.items() if S.expand(c).coeff(q,rep['source_order'])!=0}
        for N in [rep['B']+1,rep['B']+3]:
            f=t**N+2*t**max(0,N-1)
            assert S.Poly(apply_time(lead,f),t).degree()==N+rep['leading_gain'];leadchecks+=1
    # Formal marker coefficients, finite time degree, Bessel eigenoperator.
    bessel=S.expand(sum(t**m*phi**m/S.factorial(m)**2 for m in range(6)))
    error=apply_time(ops[0],bessel)
    for k in range(9):assert S.expand(error).coeff(s,k)==0
    for k in range(9):
        c=S.expand(bessel).coeff(s,k)
        if c:assert S.Poly(c,t).degree()<=k
    # A common rational-time equation with no joint-source premise.
    logexp=sum(s**(j+1)*t**j/S.Integer(j) for j in range(1,7))
    F=S.series(S.exp(logexp),s,0,7).removeO()
    error=apply_time(ops[1],F)
    for k in range(7):assert S.expand(error).coeff(s,k)==0
    assert F.coeff(s,0)==1
    for k in range(1,7):
        c=S.expand(F).coeff(s,k)
        if c:assert S.Poly(c,t).degree()<=k
    assert apply_time(ops[2],1+s*t**5)==0  # Indicial resonance cannot be ignored.
    delayed=1+s**7*(S.exp(t)-1)
    assert S.diff(delayed,t,2)-S.diff(delayed,t)==0
    assert all(S.expand(delayed).coeff(s,k)==(1 if k==0 else 0) for k in range(7))
    assert S.expand(delayed).coeff(s,7).has(S.exp(t))
    # Actual reverse source for T=D_x^2,P=x: exp(zxt+zt^2).
    reverse=[]
    for m in range(1,7):
        f=x**m;terms=[];k=0
        while f!=0:
            terms.append(f*z**(m-k)/S.factorial(k));f=S.diff(f,x,2);k+=1
        R=S.expand(sum(terms));b=min(a[0] for a,c in S.Poly(R,z).terms())
        assert b==m-m//2 and m<=2*b
        expression=S.diff(S.exp(z*x*t+z*t*t),t,m).subs(t,0)
        assert S.expand(expression-R)==0
        reverse.append(dict(m=m,height=m//2,deficit=b))
    assert S.simplify(S.diff(S.exp(z*x*t+z*t*t),t)-z*(x+2*t)*S.exp(z*x*t+z*t*t))==0
    # Full raw-source top output survives using ONLY the time direction.
    u,v=S.symbols('u v');lam=X*X+Y**3;P=x*y;tops=0
    for m in range(1,4):
        shifted=S.Poly(S.expand((lam.subs({X:X+u,Y:Y+v}))**m),X,Y)
        out=S.expand(sum(c*S.diff(P**m,x,a[0],y,a[1]) for a,c in shifted.terms()))
        top=sum(c*x**a[0]*y**a[1] for a,c in S.Poly(out,x,y).terms() if sum(a)==2*m)
        assert S.expand(top-(u*u+v**3)**m*P**m)==0;tops+=1
    h=sum((t*(u*u+v**3))**m/S.factorial(m)**2 for m in range(6))
    err=S.expand(t*S.diff(h,t,2)+S.diff(h,t)-(u*u+v**3)*h)
    assert all(err.coeff(t,k)==0 for k in range(5))
    # These prefixes illustrate the infinite written countermodel; no fitting claim.
    bad=1+sum(t**(k*k)*s**k for k in range(1,4))
    assert S.expand(bad*S.diff(bad,t)-S.diff(bad,t)*bad)==0
    guards=0
    for badop in [{(-1,0):1},{(0,1):S.Float(1)}]:
        try:indicial(badop)
        except AssertionError:guards+=1
        else:raise AssertionError('guard failure')
    deps=['verify_formal_apolar_gaps.py','verify_formal_symbol_source_rigidity.py','verify_algebraic_formal_symbols.py']
    report=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        supporting_checker_hashes={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in deps},
        indicial_controls=reports,generic_marker_line_control=mr,leading_degree_controls=leadchecks,
        bessel_marker_coefficients_checked=9,parameter_power_marker_coefficients_checked=7,
        indicial_resonance_degree=5,delayed_bad_marker_degree=7,
        actual_reversed_profiles=reverse,reversed_source_cutoff=2,
        raw_top_output_controls=tops,bessel_time_eigen_controls=5,guard_refusals=guards,
        general_formal_time_equation_produced=False,general_formal_gap_decided=False,
        scalar_equation_implies_finite_marker_test=False,formal_time_coefficients_unbounded_degree_allowed=False,
        raw_rational_time_method_exceeds_algebraic_symbols=False,full_exceptional_MZ_decided=False,
        unrestricted_frontier_endpoint_changed=False,universal_proof_by_finite_samples=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    assert guards==2
    (HERE/'scalar_time_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(completed=True,indicial_models=6,leading_degree_controls=leadchecks,actual_profiles=6,raw_top_controls=tops,guards=guards)))

if __name__=='__main__':main()
