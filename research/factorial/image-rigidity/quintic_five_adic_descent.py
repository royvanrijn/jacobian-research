#!/usr/bin/env python3
"""Exact internal replay for QUINTIC_FIVE_ADIC_DESCENT.tex.

The written proof supplies unbounded scope.  No all-terminal quintic
closure or ordinary integer counterexample is inferred from these checks.
"""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
from math import factorial

HERE=Path(__file__).resolve().parent
D=range(-4,3)
def transition(delta,s,j):
    k=(2*(delta-1-j))%5
    beta=(j+3*k-delta+1)//5
    assert 5*beta==j+3*k-delta+1
    eps=int(k<s)
    return k,beta,eps,s-beta+1-5*eps

def mul_linear(poly,a,b):
    out=[0]*(len(poly)+1)
    for i,x in enumerate(poly):
        out[i]+=a*x
        out[i+1]+=b*x
    return out

def direct(delta,n,j):
    out=[1]
    for k in range(n):out=mul_linear(out,5*n+delta-1-j-3*k,5)
    return out

def split(delta,n,j,corrupt=False):
    ell,s=divmod(n,5)
    k,beta,eps,child=transition(delta,s,j)
    t=ell+eps
    out=[5**t]
    for a in range(t):
        out=mul_linear(out,5*t+child-1-3*a+int(corrupt),1)
    for a in range(n):
        if a%5!=k:out=mul_linear(out,5*n+delta-1-j-3*a,5)
    return out

transitions=[]
for delta in D:
    for s in range(5):
        for j in range(5):
            k,beta,eps,child=transition(delta,s,j)
            assert -4<=child<=2
            transitions.append(dict(delta=delta,s=s,j=j,k=k,beta=beta,
                                    epsilon=eps,child_delta=child))
assert len(transitions)==175
assert transition(-3,1,1)[-1]==-4
identities=0
for delta in D:
    for n in range(31):
        for j in range(5):
            assert direct(delta,n,j)==split(delta,n,j)
            if n>=2:assert n//5+transition(delta,n%5,j)[2]<n
            identities+=1
assert direct(-3,6,1)!=split(-3,6,1,True)

# Formal rational Artin-Hasse coefficients, with no integrality assumption.
M=30
E=[Fraction(0) for _ in range(M+1)]
E[0]=Fraction(1)
logder={5**r:1 for r in range(4) if 5**r<=M}
for n in range(1,M+1):
    E[n]=sum(E[n-k]*v for k,v in logder.items() if k<=n)/n
def vp_int(n):
    if not n:return 10**9
    v=0
    while n%5==0:n//=5;v+=1
    return v
def vp_fraction(x):
    x=Fraction(x)
    return vp_int(x.numerator)-vp_int(x.denominator)
assert all(vp_fraction(x)>=0 for x in E)
Einv=[Fraction(0) for _ in E];Einv[0]=Fraction(1)
for n in range(1,len(E)):
    Einv[n]=-sum(E[k]*Einv[n-k] for k in range(1,n+1))
assert all(vp_fraction(x)>=0 for x in Einv)
assert Einv[:5]==[Fraction((-1)**j,factorial(j)) for j in range(5)]
seed_controls=0
for h0 in [1,2,3,4,6]:
    h=[Fraction(h0)]+[Fraction(5*((j*7)%11-5)) for j in range(1,M+1)]
    gauged=[sum(Einv[k]*h[n-k] for k in range(n+1)) for n in range(M+1)]
    for j in range(5):
        c=gauged[j]
        assert vp_fraction(c)==0
        assert vp_fraction(c-h0*Fraction((-1)**j,factorial(j)))>=1
        for k in range(1,(M-j)//5+1):
            assert vp_fraction(5**k*gauged[j+5*k])>=1
        seed_controls+=1
potentials=0
for a1 in [0,5,10,25,-5]:
    b=a1-1
    assert vp_fraction(b)==0
    for r in range(1,5):
        for anext in [0,5,15,-10]:
            aprime=Fraction((anext-1)*5**(5**r-1),b**(5**r))
            assert vp_fraction(aprime)>=5**r-1
            potentials+=1



from math import factorial
import sympy as sp
m,w,y=sp.symbols("m w y")
p=5

def vp(n):
    v=0
    while n and n%p==0:n//=p;v+=1
    return v

def poly_mul(a,b,mod):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,z in enumerate(b):out[i+j]=(out[i+j]+x*z)%mod
    while len(out)>1 and not out[-1]:out.pop()
    return out

def coeff_poly(j,K):
    if j==0:return [1]
    t=(3*j+4)//5
    e=j-t-vp(factorial(j))
    if e>=K:return [0]
    assert e>=0
    mod=p**K
    fac=factorial(j)//p**vp(factorial(j))
    out=[pow(p,e,mod)*pow(fac,-1,mod)%mod]
    for a in range(j-t):out=poly_mul(out,[-t-a,1],mod)
    for r in range(1,3*j):
        if r%p:
            out=poly_mul(out,[-pow(r,-1,mod)*pow(5*pow(r,-1,mod),a,mod)%mod for a in range(K-e)],mod)
    return out

def Qpoly(K):
    mod=p**K
    terms=[]
    # Bound j-ceil(3j/5)-v5(j!)>=3j/20-1, so j>=ceil(20(K+1)/3) suffices.
    limit=(20*(K+1)+2)//3
    for j in range(limit):
        poly=coeff_poly(j,K)
        terms.extend(c*m**a*w**j for a,c in enumerate(poly) if c)
    return sp.Poly(sum(terms),m,w, modulus=mod).as_expr()


def Qeval(mm,ww,K):
    """Exact mod5^K evaluation using a proved all-tail bound."""
    mod=5**K
    total=1
    limit=(20*(K+1)+2)//3
    for j in range(1,limit):
        t=(3*j+4)//5
        vfac=vp(factorial(j))
        e=j-t-vfac
        if e>=K:continue
        assert e>=0
        term=pow(5,e,mod)*pow(factorial(j)//5**vfac,-1,mod)%mod
        for a in range(j-t):term=term*(mm-t-a)%mod
        denominator=1
        for r in range(1,3*j):
            if r%5:denominator=denominator*(5*mm-r)%mod
        term=term*pow(denominator,-1,mod)%mod
        total=(total+term*pow(ww,j,mod))%mod
    return total

def eval_FG(c,mm,yy,K):
    mod=5**K
    ww=c+25*yy
    q=Qeval(mm,ww,K+2)
    qnext=Qeval(mm+1,ww,K+2)
    assert q%25==0 and qnext%25==0
    return (q//25)%mod,((qnext-q)%(25*mod)//25)%mod

def lift_branch(c,digits):
    mm,yy=(3,1) if c==2 else (1,4)
    trajectory=[]
    for k in range(1,digits):
        f,g=eval_FG(c,mm,yy,k+1)
        assert f%(5**k)==g%(5**k)==0
        rf,rg=(f//5**k)%5,(g//5**k)%5
        if c==2:
            dm=(-3*rg)%5
            dy=(-rf-4*dm)%5
        else:
            dm=rg
            dy=(rf+3*dm)%5
        mm+=5**k*dm
        yy+=5**k*dy
        assert eval_FG(c,mm,yy,k+1)==(0,0)
        trajectory.append(dict(precision=k+1,m=mm,y=yy,w=(c+25*yy)%(5**(k+3))))
    return mm,yy,trajectory

def exact_physical(mm,ww,K):
    mod=5**K
    total=1
    for j in range(1,mm+1):
        t=(3*j+4)//5
        e=j-t-vp(factorial(j))
        if e>=K:continue
        term=pow(5,e,mod)*pow(factorial(j)//5**vp(factorial(j)),-1,mod)%mod
        for a in range(j-t):term=term*(mm-t-a)%mod
        den=1
        for r in range(1,3*j):
            if r%5:den=den*(5*mm-r)%mod
        total=(total+term*pow(den,-1,mod)*pow(ww,j,mod))%mod
    return total

def seed_blind_counterexample():
    # Literal coefficients: ascending powers. All arithmetic is integral.
    den=1349729058338
    seed=[den,270777950790,-19101814440,-542394950285]
    quotient=[-132651371853458640,-62036179766657136,
              -11633463676444856,-1088171938849410,
              -54296800471820,-1627184850855]
    expected_sections=[-296770817069940,-112236517464324,-13462409209854]
    assert den%5!=0 and all(a%5==0 for a in seed[1:])
    def Xplus(poly,c):
        out=[0]*(len(poly)+1)
        for j,a in enumerate(poly):
            out[j]+=(c-j)*a
            out[j+1]+=a
        return out
    B=Xplus(Xplus(seed,7),10)
    C=Xplus(B,13)
    D=Xplus(C,16)
    H=[70*a for a in D]+[0]
    for j,a in enumerate(B):H[j+3]+=3*a
    lhs=H[:]
    for j,a in enumerate(C):lhs[j]-=40*a
    rhs=[0]*(len(quotient)+3)
    for j,a in enumerate(quotient):
        rhs[j]-=10*a
        rhs[j+3]+=a
    assert lhs==rhs
    corrupted=rhs[:]
    corrupted[0]+=5
    assert lhs!=corrupted
    def section_at10(poly):
        out=[0,0,0]
        for j,a in enumerate(poly):out[j%3]+=a*10**(j//3)
        return out
    cs=section_at10(C)
    hs=section_at10(H)
    assert cs==expected_sections and any(cs)
    assert hs==[40*a for a in cs]
    assert all(cs[i]*hs[j]-cs[j]*hs[i]==0 for i in range(3) for j in range(i+1,3))
    return {
        "actual_rank":1,"cyclic_argument":10,"argument_valuation":1,
        "seed_in_one_plus_5Z5_z":True,
        "seed_is_actual_seed_one":False,
        "denominator":den,"seed_numerator":seed,
        "literal_quotient_numerator":quotient,
        "nonzero_C_sections_numerator":cs,
        "corrupted_identity_rejected":True,
        "stronger_seed_blind_assertion_refuted":True,
        "ordinary_SFC_counterexample":False,
    }

def formal_checks_and_receipt():
    from pathlib import Path
    import hashlib,json
    q3=Qpoly(3)
    literal_rows=[
     ("1","1"),
     ("-50*m**2-55*m+1","25*m**2-60*m-58"),
     ("-52","-53"),
     ("-50*m**2+40*m-5","-25*m**2+20*m+60"),
     ("25*m**2-45*m+35","25*m**2+55*m-15"),
     ("35*m**2+5*m+45","-5*m**2+35*m-60"),
     ("-50*m**3-40*m**2+35*m-50","50*m**3-10*m**2+40*m+50"),
     ("-45*m**2-5*m+25","-35*m**2+10*m-50"),
     ("25*m**3+50*m**2+50*m","25*m**3+50*m**2+50*m"),
     ("50*m**3-50*m**2+50*m-50","-25*m**3+25*m**2-25*m+25"),
     ("50*m**4-50","50-50*m**4"),
     ("50*m**4+50*m**3+50*m**2+50*m","25*m**4+25*m**3+25*m**2+25*m"),
     ("25*m**4+50*m**3-25*m**2-50*m","25*m**4+50*m**3-25*m**2-50*m"),
    ]
    for j,row in enumerate(literal_rows):
        for c,literal in zip([2,9],row):
            actual=sum(a*m**i for i,a in enumerate(coeff_poly(j,3)))*c**j
            assert sp.Poly(actual-sp.sympify(literal),m,modulus=125).is_zero
    expected={
        2:25*(m-1)**2,
        9:25*(2*m*m-m-2),
    }
    assert sp.Poly(q3-(1-2*w+2*w*w),m,w,modulus=5).is_zero
    for c,v in expected.items():
        assert sp.Poly(q3.subs(w,c)-v,m,modulus=125).is_zero
    assert not sp.Poly(q3.subs(w,2)-expected[2]-25,m,modulus=125).is_zero
    q4=Qpoly(4)
    residuals={2:((m-1)**2+y,2*m-1),9:(2*m*m-m-2-y,1-m)}
    for c,(ef,eg) in residuals.items():
        loc=sp.Poly(q4.subs(w,c+25*y),m,y,modulus=625)
        assert all(int(co)%25==0 for co in loc.coeffs())
        ff=sp.Poly(loc.as_expr()/25,m,y,modulus=25)
        diff=sp.Poly((q4.subs(m,m+1)-q4).subs(w,c+25*y),m,y,modulus=625)
        assert all(int(co)%25==0 for co in diff.coeffs())
        gg=sp.Poly(diff.as_expr()/25,m,y,modulus=25)
        assert sp.Poly(ff.as_expr()-ef,m,y,modulus=5).is_zero
        assert sp.Poly(gg.as_expr()-eg,m,y,modulus=5).is_zero
    # Independent literal factorial quotients at actual integer indices.
    coefficient_controls=0
    from fractions import Fraction
    for mm in [8,9,10,11,12,13,17,24,25,31,40]:
        for j in range(mm+1):
            actual=Fraction(factorial(mm),factorial(j)*factorial(mm-j))*Fraction(factorial(5*mm-3*j),factorial(5*mm))*5**j
            t=(3*j+4)//5
            denominator=factorial(j)
            for r in range(1,3*j):
                if r%5:denominator*=5*mm-r
            numerator=5**(j-t)
            for a in range(j-t):numerator*=mm-t-a
            canceled=Fraction(numerator,denominator)
            assert actual==canceled
            coefficient_controls+=1
    # Demonstrate the essential nonzero post-resonance tail, not just truncation.
    mm=8
    resumed_j=(5*mm)//3+1
    t=(3*resumed_j+4)//5
    num=1
    for a in range(resumed_j-t):num*=mm-t-a
    assert num!=0
    tail_controls=0
    for mm in range(8,81):
        L=max(3,mm//4)
        for j in range((5*mm)//3+1,(5*mm)//3+41):
            assert j-(3*j+4)//5-vp(factorial(j))>=L
            tail_controls+=1
    comparison_controls=0
    for mm in [8,9,12,16,24,25,32,40,60]:
        L=max(3,mm//4)
        for ww in [1,2,3,4,9,27,109]:
            assert (Qeval(mm,ww,L+1)-exact_physical(mm,ww,L+1))%(5**L)==0
            comparison_controls+=1
    branches={}
    for c in [2,9]:
        mm,yy,tr=lift_branch(c,14)
        assert mm=={2:3644508588,9:2716918281}[c]
        assert eval_FG(c,mm+5**13,yy,14)!=(0,0)
        branches[str(c)]={"m_mod_5_14":mm,"y_mod_5_14":yy,
                         "w_mod_5_16":(c+25*yy)%5**16,
                         "Hensel_trajectory":tr}
    here=Path(__file__).resolve().parent
    files=[here/"QUINTIC_FIVE_ADIC_DESCENT.tex",Path(__file__).resolve()]
    out={
      "seven_offset_transfer": {
       "complete_offset_residue_transitions":len(transitions),
       "exact_integer_polynomial_factorization_controls":identities,
       "six_offset_nonclosure_detected":True,
       "corrupted_factorization_rejected":True,
       "Artin_Hasse_and_inverse_integral_coefficients_checked":len(E),
       "constant_unit_seed_controls":seed_controls,
       "closed_potential_coefficient_controls":potentials,
       "joint_three_cyclic_rank_descent_proved":False,
       "transition_table":transitions,
      },
      "status":"passed",
      "seed_blind_rank_counterexample":seed_blind_counterexample(),
      "symbolic_residual_branches":2,
      "literal_coefficientwise_residual_table_entries":26,
      "corrupted_residual_and_Hensel_controls_rejected":True,
      "complete_nonsingular_residual_common_pairs":[[3,1],[1,4]],
      "exact_literal_factorial_coefficient_controls":coefficient_controls,
      "tail_valuation_controls":tail_controls,
      "finite_polynomial_comparison_controls":comparison_controls,
      "Hensel_digits_each_branch":14,
      "analytic_common_pairs_exactly_two_by_written_proof":True,
      "analytic_family_equals_actual_finite_polynomial":False,
      "post_resonance_tail_explicitly_nonzero_at_m8_j14":True,
      "ordinary_integer_common_pair_identified":False,
      "all_terminal_quintic_SFC":False,
      "independent_replay":False,"formal_verification":False,"external_review":False,
      "unbounded_scope":"Proved Tate convergence, exactly two common analytic indices and all-integer necessary closeness; finite checks support the written proof.",
      "branches":branches,
      "input_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    }
    (here/"quintic_five_adic_descent_verification.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:v for k,v in out.items() if k not in {"branches","seven_offset_transfer"}},indent=2))
    print("complete_seven_offset_transitions",len(transitions))
    print("branch_m_residues_mod5^14",{c:v["m_mod_5_14"] for c,v in branches.items()})

if __name__=="__main__":
    formal_checks_and_receipt()
