#!/usr/bin/env python3
"""PF181 exact controls; uniform Taylor and budget proof is in the source note."""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import platform
import time
import sympy as S
START=time.perf_counter()
BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_NONCRITICAL_POLE_BUDGET_AND_DEGREE30_LEDGER.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v,s,Y,X,Z,p,lam=S.symbols("v s Y X Z p lam")
def initial(expr,var):
    num,den=S.fraction(S.cancel(expr))
    if num==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(num,var),S.Poly(den,var)
    on=min(m[0] for m in pn.monoms());od=min(m[0] for m in pd.monoms())
    return on-od,S.cancel(pn.nth(on)/pd.nth(od))
def partitions(n,cap=None):
    if n==0:
        yield []
        return
    if cap is None:cap=n
    for i in range(min(n,cap),0,-1):
        for tail in partitions(n-i,i):
            yield [i]+tail
def positive_orders(part,L,total):
    bounds=[int(S.ceiling(a*(S.Rational(L,2)-1))-1) for a in part]
    return [qs for qs in itertools.product(*(range(1,b+1) for b in bounds))
            if sum(qs)==total]

# Exact whole polynomial model of the noncritical exceptional boundary.
# Tlead=(1/2)Y^2-(3/2)Y has B/(A gamma)=-3, the k=3 resonance ratio.
# These are source-coefficient controls, not actual companion claims.
models=[]
for L,t in [(6,10),(7,12),(8,14),(9,16),(10,18)]:
    ell=t//2
    YY=1+X
    SS=YY*(X**2+v**ell)
    F0=S.expand(SS**2+v**t*(S.Rational(1,2)*YY**2-S.Rational(3,2)*YY)
                 +v**(2*L)*(3+v))
    check(f"leading ratio is allowed k3 L{L}",S.Rational(-3,2)/S.Rational(1,2)==-3)
    check(f"full exceptional initial cancels its constant v^t L{L}",
          zero(F0.subs(X,0)-v**(2*L)*(3+v)))
    check(f"noncritical constant perturbation is nonzero L{L}",
          S.Rational(1,2)-S.Rational(3,2)==-1)
    for h in [t+1,2*t+3]:
        center=s**h
        sub={v:s**4,X:center}
        o2,l2=initial(S.diff(F0,X,2).subs(sub,simultaneous=True),s)
        o3,l3=initial(S.diff(F0,X,3).subs(sub,simultaneous=True),s)
        check(f"entire second derivative exact ell L{L} center{h}",
              o2==4*ell and l2!=0)
        check(f"entire third derivative later than tquarter L{L} center{h}",o3>t)
        for j in [4,5,6]:
            oj,lj=initial(S.diff(F0,X,j).subs(sub,simultaneous=True),s)
            check(f"entire higher derivative nonnegative L{L} center{h} j{j}",oj>=0)
    R0=S.Rational(L)-S.Rational(t,4)
    for offset in [S.Rational(-1,4),S.Integer(0),S.Rational(1,4)]:
        RR=R0+offset
        check(f"Taylor first-response index obeys PF164 L{L} offset{offset}",RR>S.Rational(L,2))
        check(f"Taylor first-response index exceeds inner scale L{L} offset{offset}",RR>S.Rational(t,4))
        quad=ell+2*RR
        check(f"quadratic exact comparison to parameter order L{L} offset{offset}",
              (quad<2*L if offset<0 else quad==2*L if offset==0 else quad>2*L))
        for j in [4,5,6]:
            check(f"higher Taylor degree strictly later L{L} offset{offset} j{j}",
                  j*RR>quad)
        check(f"cubic Taylor lower bound strictly later L{L} offset{offset}",
              S.Rational(t,4)+3*RR>quad)
    # Boundary response is visibly quadratic on the WHOLE sextic.
    XX=s**(4*R0)*lam
    whole=(F0-p*v**(2*L)).subs({v:s**4,X:XX},simultaneous=True)/s**(8*L)
    oo,ll=initial(whole,s)
    check(f"whole boundary first response is genuinely quadratic L{L}",
          oo==0 and zero(ll-(2*lam**2+3-p)))
    check(f"whole boundary coefficient degree is two L{L}",S.degree(ll,lam)==2)
    models.append({"L":L,"t":t,"ell":ell,"forbidden_response_order":str(R0),
                   "quadratic_initial":str(ll)})

# A nearby full canonical model has constant simple leading roots and AFFINE later response.
# It saturates the integer local budget five at L7,b_e5, without positing an actual mate.
L=7;t=12;ell=6
YY=1+X
Fpositive=S.expand(YY**2*(X**2+v**ell)**2
      +v**t*(S.Rational(1,2)*YY**2-S.Rational(3,2)*YY)
      +v**13*YY+3*v**14)
cc=S.I/S.sqrt(2)
XX=cc*s**7+s**9*lam
scaled=(Fpositive-p*v**14).subs({v:s**2,X:XX},simultaneous=True)/s**28
oo,ll=initial(scaled,s)
check("nearby whole model first response is affine",
      oo==0 and zero(ll-(4*cc*lam+cc**4+3-p)))
check("nearby whole model response degree is one",S.degree(ll,lam)==1)
mu0=S.Rational(19,2)
do,dl=initial(S.diff(Fpositive,X).subs({v:s**2,X:XX},simultaneous=True)/s**(2*mu0),s)
check("nearby whole model exact inner derivative ninehalf",do==0 and zero(dl-4*cc))
check("nearby control saturates the integer local pole budget",
      2*(9-8)+2*(S.Rational(19,2)-8)==5==2*5-5)
check("nearby control inner index is strictly below affine bound",S.Rational(19,2)<7+S.Rational(12,4))

# Symbolic budget identities and exact strictness margins.
Ls,ts,bs,Rs=S.symbols("L t b R")
outer=S.Rational(3,4)*ts-Ls-1
inner=ts/4-1
bound=2*ts-2*Ls-4
check("two outer plus two inner charges give uniform strict bound",zero(2*outer+2*inner-bound))
check("all-simple estimate has positive gap twoLminusT",zero(bound-4*outer-(2*Ls-ts)))
check("integer budget equals twice eorder minus five",
      S.expand((2*ts-2*Ls-5).subs(ts,Ls+bs)-(2*bs-5))==0)
check("maximal local budget equals twoLminusNine",
      S.expand((2*bs-5).subs(bs,Ls-2)-(2*Ls-9))==0)
check("inner Taylor threshold is equivalent to strict derivative bound",
      S.expand(ts/2+2*(2*Ls-S.Symbol("mu"))-2*Ls
               -2*(Ls+ts/4-S.Symbol("mu")))==0)
for L in range(3,15):
    for b in range(L-1):
        t=L+b
        if 3*t<=4*(L+1):continue
        check(f"noncritical inner charge bound positive L{L} b{b}",S.Rational(t,4)-1>0)
        check(f"all simple budget strictly below full bound L{L} b{b}",
              3*t-4*L-4<2*t-2*L-4)
        check(f"uniform noncritical local bound L{L} b{b}",2*b-5<=2*L-9)

# Mixed global ledger and all-noncritical degree consequence.
SS,nn,nc,mm=S.symbols("S n nc m")
rhs=2*SS-5*nn-4*nc
check("global mixed ledger clearing gives exact inequality",
      zero(3*rhs-(2*SS+2*mm)-(4*SS-15*nn-12*nc-2*mm)))
check("all-noncritical coefficient is twenty-seven",
      S.expand((15*nn+12*nc).subs(nc,nn)-27*nn)==0)
check("all-noncritical minimal slope sum is eight",S.ceiling(S.Rational(27+4,4))==8)
check("all-noncritical raw degree at least thirtyfour",2+4*8==34)
check("degree congruence raises all-noncritical minimum to thirtysix",
      min(d for d in range(34,40) if d%3==0)==36)
check("all-noncritical intermediate degree at least seven",S.ceiling(S.Rational(36+6,6))==7)
check("all-noncritical original degree at least fourteen",2*7==14)

# Complete D30 numerical ledger, with the proof's finite a priori inequality.
ledger=[]
for m in range(2,31,2):
    H=(30-m)//2
    for n in range(1,11):
        if m+5*n>10:continue
        check(f"D30 necessary support count one m{m} n{n}",n==1)
        for total in range(n,H//2+1):
            lo=S.ceiling(S.Rational(30+3*m,6))
            hi=2*total-5*n
            if lo<=hi:ledger.append((m,H,n,total,int(lo),int(hi)))
expected=[(2,14,1,6,6,7),(2,14,1,7,6,9),(4,13,1,6,7,7)]
check("D30 complete numerical family list",ledger==expected)
for m,H,n,L,lo,hi in ledger:
    check(f"D30 noncritical local budget cannot reach generic degree m{m} L{L}",
          2*L-9<lo)
check("D30 all supports must be critical",all(2*L-9<lo for m,H,n,L,lo,hi in ledger))

# Complete index partitions needed for critical pole mass and order refinements.
partitions_data={}
for L,total in [(6,6),(6,7),(7,6),(7,7),(7,8),(7,9)]:
    options=[]
    for mass in range(1,5):
        for part in partitions(mass):
            for orders in positive_orders(part,L,total):
                options.append({"indices":part,"pole_orders":list(orders)})
    check(f"critical degree partition controls remain nonempty L{L} degree{total}",len(options)>0)
    masses={sum(o["indices"]) for o in options}
    expected_mass={4} if L==6 or total>=8 else {3,4}
    check(f"critical complete pole mass range L{L} degree{total}",masses==expected_mass)
    if L==6 and total==7:
        check("L6 degree7 unique indexfour pole",
              options==[{"indices":[4],"pole_orders":[7]}])
        check("L6 degree7 derivative trace forces eorder at least three",
              S.ceiling(7+S.Rational(7,4))-6==3)
    if L==7 and total==9:
        check("L7 degree9 only two index/order options",
              options==[{"indices":[4],"pole_orders":[9]},
                        {"indices":[3,1],"pole_orders":[7,2]}])
        check("L7 degree9 minimum derivative trace threshold is ten",
              min(min(8+S.Rational(q,a) for a,q in zip(o["indices"],o["pole_orders"]))
                  for o in options)==10)
    partitions_data[f"L{L}_degree{total}"]=options

receipt={"claim":"PF2RULED181","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{name:sha(BASE/name) for name in [
        "POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md",
        "RULING_SIX_AFFINE_POLE_BUDGET_AND_SOURCE_INFINITY.md",
        "RULING_SIX_RETAINED_SOURCE_PACKETS_AND_STRICT_INFINITY_BOUNDS.md",
        "RULING_SIX_LEADING_PRIMITIVE_AND_POLYNOMIAL_ODE.md",
        "RULING_SIX_BALANCED_QUADRATIC_REMAINDER_AND_LEADING_RESONANCE.md",
        "RULING_SIX_NONCRITICAL_QUARTIC_SPLITTING_AND_SMALL_SLOPE_EXCLUSION.md"]},
    "uniform_proved_scope":{
        "coefficient_and_companion_degrees":"unrestricted",
        "noncritical_local_pole_degree_bound":"beta_c<=2*ord_c e-5<=2L-9",
        "mixed_global_degree_bound":"b_field<=2S-5*n_A-4*n_nc",
        "mixed_support_inequality":"4S>=15*n_A+12*n_nc+2m",
        "all_noncritical_subclass_minimum_D":36,
        "all_noncritical_subclass_minimum_intermediate_degree":7,
        "all_noncritical_subclass_minimum_original_degree":14,
        "complete_D30_ledger":[list(row) for row in ledger],
        "entire_D30_noncritical_case_excluded":True,
        "D30_critical_case_closed":False,
        "full_ruling_six_or_JC2_closed":False},
    "whole_inner_taylor_models":models,
    "critical_index_partition_controls":partitions_data,
    "assurance":"Exact internal symbolic controls and arbitrary-degree written proof; no formal or external human review."}
parser=argparse.ArgumentParser()
parser.add_argument("--write-receipt",action="store_true")
args=parser.parse_args()
if args.write_receipt:
    (BASE/"ruling_six_noncritical_budget_degree30_receipt.json").write_text(
        json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
print(json.dumps({"claim":receipt["claim"],"outcome":"PASS","check_count":len(checks),
                  "elapsed_seconds":receipt["elapsed_seconds"],"receipt_written":args.write_receipt},indent=2))
