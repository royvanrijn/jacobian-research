"""PF178 exact controls; the complete unbounded proof is in the source note."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import tempfile
import time
import sympy as S

START=time.perf_counter();BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_BALANCED_QUADRATIC_REMAINDER_AND_LEADING_RESONANCE.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(value):return S.cancel(S.expand(value))==0
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
Y,x=S.symbols("Y x")
A,B,C,gamma=S.symbols("A B C gamma",nonzero=True)
r=S.symbols("r",integer=True,positive=True)
check("all-r derivative bracket identity",zero((r-1)*(A*gamma+B)+r*A*gamma-((2*r-1)*A*gamma+(r-1)*B)))
check("global odd-index ratio identity",zero((2*r-1)/(r-1)-2*(2*r-1)/((2*r-1)-1)))
F=C*Y*(Y-gamma)**2
for rr in range(1,8):
    numerator=Y**(rr-1)*(A*Y+B)**rr
    direct=S.diff(numerator,Y).subs(Y,gamma)/C
    formula=gamma**(rr-2)*(A*gamma+B)**(rr-1)*((2*rr-1)*A*gamma+(rr-1)*B)/C
    check("full balanced coefficient r"+str(rr),zero(direct-formula))
    # Independent Laurent expansion of one actual rational initial.
    numeric=(Y*(2*Y+3))**rr/(7*Y*(Y-5)**2)
    actual=S.expand(S.series(numeric.subs(Y,1/x),x,0,2).removeO()).coeff(x,1)
    predicted=formula.subs({A:2,B:3,C:7,gamma:5})
    check("whole Laurent initial r"+str(rr),zero(actual-predicted))
    pure_quad=S.diff(Y**(2*rr-1),Y).subs(Y,gamma)*A**rr/C
    check("pure quadratic coefficient r"+str(rr),zero(pure_quad-(2*rr-1)*A**rr*gamma**(2*rr-2)/C))
    check("pure quadratic nonzero r"+str(rr),pure_quad!=0)
    pure_lin=S.diff(Y**(rr-1),Y).subs(Y,gamma)*B**rr/C
    check("pure linear coefficient r"+str(rr),zero(pure_lin-(rr-1)*B**rr*gamma**(rr-2)/C))
    if rr>=2:
        check("pure linear nonzero r"+str(rr),pure_lin!=0)
        check("critical-root cancellation r"+str(rr),zero(formula.subs(B,-A*gamma)))
        ratio=-S.Rational(2*rr-1,rr-1)
        check("top companion-ratio cancellation r"+str(rr),zero(formula.subs(B,ratio*A*gamma)))
        previous=gamma**(rr-3)*(A*gamma+B)**(rr-2)*((2*rr-3)*A*gamma+(rr-2)*B)/C
        check("ratio cancellation does not kill previous coefficient r"+str(rr),S.factor(previous.subs(B,ratio*A*gamma))!=0)
        following=gamma**(rr-1)*(A*gamma+B)**rr*((2*rr+1)*A*gamma+rr*B)/C
        check("ratio cancellation does not kill following coefficient r"+str(rr),S.factor(following.subs(B,ratio*A*gamma))!=0)

# Full binomial finite-range and strict valuation bounds, without degree cutoffs
# in the proof. The loops are exact boundary controls, not that proof.
for L in range(3,9):
    for q in range(1,L+1):
        for rr in range(1,6):
            n=2*rr-1
            first=-L*(n+1)+(2*L-q)*rr
            check("first binomial order L"+str(L)+" q"+str(q)+" r"+str(rr),first==-rr*q)
            if rr>=2:
                check("top pole strictly dominates L"+str(L)+" q"+str(q)+" r"+str(rr),-rr*q<-(rr-1)*q)
            for jj in range(rr+1,(3*n+1)//4+1):
                check("all later binomial margins L"+str(L)+" q"+str(q)+" r"+str(rr)+" j"+str(jj),
                      -L*(n+1)+(2*L-q)*jj>-rr*q)
for L in range(3,10):
    allowed=[(a,b) for a in range(L,2*L+1) for b in range(0,L+1)
             if a==b+L and (a<2*L or b<L) and a!=2*L-1]
    check("balanced valuation interval L"+str(L),allowed==[(b+L,b) for b in range(L-1)])
check("slope-five exact valuation pairs",[(b+5,b) for b in range(4)]==[(5,0),(6,1),(7,2),(8,3)])
# Coefficient residue scaling is exact even on a full finite Laurent sum.
v,w=S.symbols("v w",nonzero=True)
for L in [3,5,7]:
    series=2*w**3+3*w+5+7/w+11/w**2
    transformed=S.expand(series.subs(w,v**(-L)*Y))
    coefficient=S.expand(transformed).coeff(Y,-1)
    check("coefficient-residue scaling L"+str(L),zero(v**(-L)*coefficient-7))

receipt={
 "claim":"PF2RULED178","outcome":"PASS","check_count":len(checks),"checks":checks,
 "elapsed_seconds":time.perf_counter()-START,
 "environment":{"python":platform.python_version(),"sympy":S.__version__},
 "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
 "retained_proof_sources":{name:sha(BASE/name) for name in [
   "RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md",
   "RULING_SIX_LEADING_PRIMITIVE_AND_POLYNOMIAL_ODE.md",
   "RULING_SIX_QUADRATIC_REMAINDER_LOCAL_ORDERS.md",
   "RULING_SIX_ODD_HALF_POWER_RESIDUE_HIERARCHY.md",
   "RULING_SIX_COMPANION_RESONANCE_AND_LINEAR_REMAINDER_CLOSURE.md",
 ]},
 "uniform_proved_scope":{
   "both_remainder_coefficients_nonzero":True,
   "every_actual_A_support":"ord d=ord e+L; 0<=ord e<=L-2",
   "leading_ratio":"A*gamma+B=0 or B/(A*gamma)=-2k/(k-1), k>=3 odd",
   "companion_index_bounded":False,
 },
 "remaining_balanced_quadratic_class_closed":False,
 "automatic_target_integrality_proved":False,"full_goal_complete":False,
 "independent_whole_proof_review":False,"formal_verification":False,"external_review":False,
}
parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path)
args=parser.parse_args();out=args.output or Path(tempfile.mkdtemp(prefix="math-check-"))/"ruling_six_balanced_quadratic_resonance_receipt.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),"elapsed_seconds":receipt["elapsed_seconds"],"receipt":str(out)},indent=2))
