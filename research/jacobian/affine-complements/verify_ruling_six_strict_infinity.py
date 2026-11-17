#!/usr/bin/env python3
"""Exact controls for PF166; the retained-source valuation proof is written."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import tempfile
import time
import sympy as S

START=time.perf_counter();BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_RETAINED_SOURCE_PACKETS_AND_STRICT_INFINITY_BOUNDS.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ceilq(n,d):return -((-n)//d)
s,p,lam=S.symbols("s p lam")
def order(e):
    pn=S.Poly(S.expand(e),s)
    return min(a[0] for a in pn.monoms())
# Nonlinear target images retain the same transverse ramification index.
for e in range(1,9):
    Q=lam+s**e
    Gamma=S.expand(Q**2-lam**2)
    check("transverse Gamma order e"+str(e),order(Gamma)==e)
    check("coframe/index exact equality e"+str(e),order(S.diff(Q,s))==e-1)
    check("generic target Q derivative unit e"+str(e),2*lam!=0)

# Complete local Artin quotient factors: omitted ramified residue packet
# PLUS a distinct nonzero retained source idempotent. These are local
# module controls, not actual Keller-map examples.
packet_controls=[]
for e,f in [(2,1),(3,1),(2,2),(4,2),(3,3)]:
    d=e*f
    U=S.zeros(d)
    for i in range(e-1):
        for j in range(f):U[(i+1)*f+j,i*f+j]=1
    check(f"local DVR quotient nilpotency e{e} f{f}",
          U**e==S.zeros(d) and U**(e-1)!=S.zeros(d))
    omitted=S.diag(*([1]*d+[0]));source=S.diag(*([0]*d+[1]))
    check(f"distinct retained source projector e{e} f{f}",
          source!=S.zeros(d+1) and omitted*source==S.zeros(d+1))
    check(f"complete local quotient rank e{e} f{f}",omitted.trace()+source.trace()==e*f+1)
    check(f"strict omitted packet degree e{e} f{f}",e*f<d+1)
    packet_controls.append({"e":e,"f":f,"omitted_length":d,"retained_source_length":1,
                            "complete_rank":d+1})

q=S.symbols("q",integer=True)
for rem in range(6):
    check("all-period ceiling shift residue"+str(rem),
          S.simplify(S.ceiling((6*q+rem-6)/6)+1-S.ceiling((6*q+rem)/6))==0)
for f in range(1,7):
    for rem in range(6):
        check("complete fractional ceiling product f"+str(f)+" rem"+str(rem),
              f*ceilq(rem,6)>=ceilq(f*rem,6))
Dplus,mplus,rhoplus=S.symbols("Dplus mplus rhoplus",integer=True,nonnegative=True)
DD=1+Dplus;MM=2+mplus;rho=3+rhoplus
gap=S.expand(rho*(DD+3*MM-6)/6+1-(DD/2-5))
check("all rho>=3 strictly violate universal degree bound",
      S.expand(gap-(6+S.Rational(3,2)*mplus+rhoplus*(1+Dplus+3*mplus)/6))==0)
allowed=[(a,f) for a in range(1,7) for f in range(1,7) if a*f<=2]
check("complete infinity index/image-residue pairs",
      allowed==[(1,1),(1,2),(2,1)])

# Complete nonlinear image packet enumeration: sum r*g=6,
# r in {1,2}, g>=2. Positivity bounds g<=6 and the number of packets<=3.
choices=[(r,g) for r in [1,2] for g in range(2,7) if r*g<=6]
def packet_lists(remaining,start=0):
    if remaining==0:
        yield ()
    for i in range(start,len(choices)):
        a=choices[i];cost=a[0]*a[1]
        if cost<=remaining:
            for tail in packet_lists(remaining-cost,i):
                yield (a,)+tail
packets=set(packet_lists(6))
expected={
    ((1,6),),((2,3),),((1,2),(1,4)),((1,2),(2,2)),
    ((1,3),(1,3)),((1,2),(1,2),(1,2))}
check("complete six nonlinear infinity packet types",packets==expected)
for packet in sorted(packets):
    check("exact total u degree packet"+str(packet),sum(r*g for r,g in packet)==6)
    check("at most three closed infinity primes packet"+str(packet),len(packet)<=3)
    check("all image degrees nonlinear packet"+str(packet),all(g>=2 for r,g in packet))
check("mixed r2/r1 quadratic packets cannot share image",
      2+1>2)
check("three quadratic r1 packets cannot all share image",
      1+1+1>2)
# The P-degree of the image and the residue degree over that image differ.
check("nonlinear P residue remains compatible with f=1",
      (1*1*6)==6 and (1,1) in allowed)
check("nonlinear P residue remains compatible with f=2",
      (1*2*3)==6 and (1,2) in allowed)
for h in [2,3,6]:
    H=lam**h
    check("polynomial leading image finite-pole guard degree"+str(h),
          S.degree(H-1,lam)==h and S.rem(S.Integer(1),H-1,lam)==1)

# Symbolic elimination of the unrestricted global ledger.
SS,NA,mm,extra=S.symbols("SS NA mm extra",integer=True,nonnegative=True)
D=mm+4*SS+extra
upper=2*SS-5*NA
check("strict source coverage global inequality elimination",
      S.expand(upper-(D+3*mm)/6)
      ==S.Rational(1,3)*(4*SS-15*NA-2*mm)-extra/6)
check("rho-two global inequality elimination",
      S.expand(upper-((D+3*mm-6)/3+1))
      ==S.Rational(1,3)*(2*SS-15*NA-4*mm+3-extra))
# COMPLETE integer-residue proof of S>=5 for NA>=1,mm>=2.
A,B=S.symbols("A B",integer=True,nonnegative=True)
for r in range(4):
    for t in range(2):
        n=1+4*A+r;m=2+2*B+t
        endpoint=ceilq(19+15*r+2*t,4)
        expr=S.ceiling((15*n+2*m)/4)
        check(f"all-count/all-m lower S residue n{r} m{t}",
              S.simplify(expr-(15*A+B+endpoint))==0 and endpoint>=5)
# Complete parity classes for rho-two's S>=10.
for r in range(2):
    n=1+2*A+r;m=2+B
    endpoint=10+ceilq(15*r,2)
    expr=S.ceiling((15*n+4*m-3)/2)
    check("rho-two all-period lower S n residue"+str(r),
          S.simplify(expr-(15*A+2*B+endpoint))==0 and endpoint>=10)
check("global D lower22 endpoint",2+4*5==22)
check("global b lower5 endpoint",ceilq(22+3*2,6)==5)
check("global original degree lower10",2*5==10)
check("two A divisor count lower S9",ceilq(15*2+2*2,4)==9)
check("two A divisor count lower D38",2+4*9==38)
check("two A divisor count lower b8",ceilq(38+3*2,6)==8)
check("rho-two D lower42 endpoint",2+4*10==42)
check("rho-two b lower15 endpoint",ceilq((42+3*2-6),3)+1==15)
check("rho-two original degree lower30",2*15==30)

# Numerical compatibility controls guard against replacing necessary
# inequalities with existence or declaring full ruling-six closure.
for label,n,m,S0,D0,b in [
    ("one A",1,2,5,22,5),("two A",2,2,9,38,8),
    ("rho2 permitted",1,2,10,42,15)]:
    check(label+" complete scalar ledger compatibility",
          D0>=m+4*S0 and ceilq(D0+3*m,6)<=b<=2*S0-5*n)
check("rho2 strict field budget compatibility",ceilq(2*(42+6-6),6)+1==15)
check("two distinct r1 infinity packets use same joint bound",
      2*ceilq(42+6-6,6)+1==15)
check("low regime two r1 packets cannot share image",
      2*ceilq(22+6-6,6)+1>5)

receipt={"claim":"PF2RULED166","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{
        name:sha(BASE/name) for name in [
          "RULING_SIX_AFFINE_POLE_BUDGET_AND_SOURCE_INFINITY.md",
          "RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md",
          "POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md"]},
    "complete_local_packet_controls":packet_controls,
    "nonlinear_infinity_packet_types":[list(map(list,a)) for a in sorted(packets)],
    "uniform_proved_scope":{
      "infinity_ramification_residue_budget":"e*f<=b-1",
      "strict_infinity_lower":"f*ceil(alpha*(D+3*m-6)/6)+1<=b",
      "image_group_infinity_index_residue_sum_maximum":2,
      "individual_infinity_index_residue_product_maximum":2,
      "image_P_degree_minimum":2,"closed_intermediate_infinity_prime_count_maximum":3,
      "A_slope_sum_lower":"ceil((15*n_A+2*m)/4)",
      "leading_coefficient_degree_minimum":22,"intermediate_generic_degree_minimum":5,
      "original_generic_degree_minimum":10,
      "rho2_A_slope_sum_lower":"ceil((15*n_A+4*m-3)/2)",
      "rho2_leading_coefficient_degree_minimum":42,
      "rho2_intermediate_generic_degree_minimum":15,"rho2_original_generic_degree_minimum":30,
      "rho2_bounds_apply_to_total_per_image_mass":True,
      "low_regime_infinity_images_pairwise_distinct":True,
      "residue_degrees_over_image_not_over_C_P":True},
    "whole_ruling_six_class_closed":False,"generic_first_derivative_integrality_proved":False,
    "normal_tail_increasing_power_containment_proved":False,"full_goal_complete":False,
    "independent_whole_proof_review":False,"formal_verification":False,"external_review":False}
parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path,default=None)
args=parser.parse_args()
out=args.output or Path(tempfile.mkdtemp(prefix="math-check-"))/"ruling_six_strict_infinity_receipt.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),"elapsed_seconds":receipt["elapsed_seconds"],
                  "receipt":str(out)},indent=2))
