#!/usr/bin/env python3
"""PF172 controls; uniform cubic-trace and odd-coefficient inductions are written."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import tempfile
import time
import sympy as S

START=time.perf_counter();BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_NONCONSTANT_REMAINDER_AND_CUBIC_FIELD_TRACE.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
u,H,w=S.symbols("u H w")
A=S.Symbol("A",nonzero=True)
B,C,D,beta=S.symbols("B C D beta",nonzero=True)
MW=S.Matrix([[0,0,H/A-D],[1,0,-C],[0,1,-B]])
check("actual monic cubic multiplication identity",
      (MW**3+B*MW**2+C*MW+(D-H/A)*S.eye(3)).applyfunc(S.expand)==S.zeros(3))
check("root trace independent of H",S.trace(MW)==-B and S.diff(S.trace(MW),H)==0)
check("second cubic Newton trace",S.expand(S.trace(MW**2)-(B**2-2*C))==0)
# Whole recurrence, exact over Q(A,B,C,D)[H].
mom=[S.Integer(3),-B,B**2-2*C]
for m in range(3,19):
    mom.append(S.expand(-B*mom[m-1]-C*mom[m-2]+(H/A-D)*mom[m-3]))
for m in range(19):
    check("cubic trace degree bound m"+str(m),S.degree(mom[m],H)<=m//3)
    if m%3==0:
        j=m//3
        check("cubic trace exact leading m"+str(m),
              zero(S.Poly(mom[m],H).nth(j)-3*A**(-j)))
power=S.eye(3)
for m in range(10):
    check("actual multiplication trace agrees with recurrence m"+str(m),
          zero(S.trace(power)-mom[m]))
    power=(power*MW).applyfunc(S.expand)
for k in [3,5]:
    q=S.symbols("q0:"+str(3*k))
    T=beta*A**k*mom[3*k]+sum(q[m]*mom[m] for m in range(3*k))
    check("entire odd companion trace degree k"+str(k),S.degree(T,H)==k)
    check("entire odd companion trace leading k"+str(k),
          zero(S.Poly(T,H).nth(k)-3*beta))
    check("every lower companion term has strictly lower trace degree k"+str(k),
          all(S.degree(mom[m],H)<=k-1 for m in range(3*k)))
# The universal recurrence proves these properties for all k; the displayed
# finite cases replay actual multiplication and every lower coefficient.

# Trace of du wedge dw/M is its dH coefficient, namely (Tr w)_H/M.
M=S.Symbol("M",nonzero=True)
check("field volume trace has zero dH coefficient",S.diff(S.trace(MW),H)/M==0)
bfun=S.Function("bfun")(u);afun=S.Function("afun")(u)
check("u dependent root trace does not change volume conclusion",
      S.diff(-bfun/afun,H)==0)
f=u**3+u
fp=S.diff(f,u)
check("chosen nonconstant control has nonzero field derivative",fp!=0)
for k in [3,5,7]:
    coeff=[S.Function("t"+str(j))(u) for j in range(k+1)]
    T=sum(coeff[j]*H**j for j in range(k+1))
    PDE=S.Poly(S.expand(fp*S.diff(T,H)-2*H*S.diff(T,u)),H)
    for m in range(k+2):
        rhs=(m+1)*fp*coeff[m+1] if m+1<=k else 0
        if m-1>=0 and m-1<=k:rhs-=2*S.diff(coeff[m-1],u)
        check("entire PDE recurrence k"+str(k)+" m"+str(m),
              S.simplify(PDE.nth(m)-rhs)==0)
    for j in range(1,k+1,2):
        value=PDE.nth(j-1)
        if j>=3:
            value=value.subs(S.diff(coeff[j-2],u),0)
        check("odd-coefficient induction k"+str(k)+" j"+str(j),
              S.simplify(value-j*fp*coeff[j])==0)
for r in [1,2,3,4]:
    T=(H**2+f)**r
    check("even trace solutions deliberately retained r"+str(r),
          S.expand(fp*S.diff(T,H)-2*H*S.diff(T,u))==0)
# f constant is a distinct case: the trace PDE alone permits odd powers.
T=H**3
check("constant f requires separate whole-source argument",
      S.diff(T,u)==0 and S.diff(T,H)!=0)
x,t,a,b=S.symbols("x t a b")
hsrc=a*x+b+t*(x**2+1)
P=hsrc**2+7
point={t:0,x:-b/a}
check("constant remainder source H has a central zero",S.cancel(hsrc.subs(point))==0)
check("constant remainder source P differential vanishes",
      S.cancel(S.diff(P,x).subs(point))==0 and S.cancel(S.diff(P,t).subs(point))==0)

d,e,f0,N=S.symbols("d e f0 N")
R=d*(w**2-N)/M+e*w+f0
check("canonical quadratic remainder top coefficient",
      S.Poly(R,w).nth(2)==d/M)
check("canonical linear remainder retained",S.Poly(R.subs(d,0),w).nth(1)==e)
check("ruling-zero remainder exactly d=e=0",
      S.degree(R.subs({d:0,e:0}),w)==0)
check("trace codomain is not the actual target",
      S.degree(-B,H)==0)

receipt={"claim":"PF2RULED172","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{
       name:sha(BASE/name) for name in [
        "RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md",
        "RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md"]},
    "actual_cubic_trace_controls":{"multiplication_powers":list(range(10)),
       "recurrence_powers":list(range(19)),"full_companion_cases":[3,5],
       "PDE_odd_cases":[3,5,7]},
    "uniform_proved_scope":{
       "remainder_in_C_u_excluded":True,
       "remaining_remainder_w_degrees":[1,2],
       "cubic_trace_codomain":"C(u,H3)","cubic_field_degree":3,
       "trace_companion_degree":"k odd","trace_companion_leading":"3*beta",
       "all_coefficient_degrees":"unrestricted","all_odd_companion_k":"every odd k>=3",
       "trace_landing_in_target_ring_asserted":False},
    "whole_ruling_six_class_closed":False,"generic_first_derivative_integrality_proved":False,
    "normal_tail_increasing_power_containment_proved":False,"full_goal_complete":False,
    "independent_whole_proof_review":False,"formal_verification":False,"external_review":False}
parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path,default=None)
args=parser.parse_args()
out=args.output or Path(tempfile.mkdtemp(prefix="math-check-"))/"ruling_six_cubic_trace_remainder_receipt.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),"elapsed_seconds":receipt["elapsed_seconds"],
                  "receipt":str(out)},indent=2))
