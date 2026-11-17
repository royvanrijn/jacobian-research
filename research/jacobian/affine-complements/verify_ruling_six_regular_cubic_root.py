"""PF171 exact nested-escape controls; unrestricted proofs are in the note."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import tempfile
import time
import sympy as S

START=time.perf_counter()
BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md"
checks=[]
models=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(value):return S.cancel(S.expand(value))==0
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def initial(value,parameter):
    terms=S.Add.make_args(S.expand(value))
    orders=[t.as_powers_dict().get(parameter,0) for t in terms]
    order=min(orders)
    lead=S.expand(sum(t/parameter**a for t,a in zip(terms,orders) if a==order))
    return order,lead

v,s,w,Y,p=S.symbols("v s w Y p")

# The exact regular cubic has TWO separate Laurent roots; no truncated
# approximate root is substituted into the purported full equation.
for k,L,a,d,dc in [
    (2,3,2,None,0),       # cubic remainder dominant; half-integral refinement
    (3,5,2,None,0),       # near the essential outer-face inequality
    (3,5,7,1,1),          # quadratic remainder dominant
    (3,5,6,1,2),          # equal orders, leading constants do not cancel
    (2,4,5,1,2),
]:
    j=k+L
    W=w*(v**k*w+1)*(v**L*w+1)
    D=S.Integer(0) if d is None else dc*v**d
    P=S.expand(W**2+v**a*w**3+D*w**2)
    wo=-v**(-L)
    nu=min(a-3*L,d-2*L if d is not None else S.oo)
    beta=(S.Rational(nu)-2*k+2*L)/2
    derivative_order=S.Rational(nu,2)+k-L
    label=f"k{k} L{L} a{a} d{d} c{dc}"
    check("outer degree-six face bound "+label,a+L>2*k)
    check("outer root is exact "+label,zero(W.subs(w,wo)))
    oo,ol=initial(S.diff(W,w).subs(w,wo),v)
    check("outer cubic derivative exact order "+label,oo==k-L and ol==1)
    check("outer displacement finer than centre "+label,beta+L>0)
    check("outer coframe strictly holomorphic "+label,derivative_order<0)
    ww=-s**(-2*L)+s**(2*beta)*Y
    scaled=S.expand(s**(-2*nu)*(P.subs({v:s**2,w:ww})-p))
    order,lead=initial(scaled,s)
    check("outer entire equation has regular scaled coefficients "+label,order==0)
    check("outer entire initial is quadratic plus nonzero constant "+label,
          S.degree(lead,Y)==2 and lead.coeff(Y,2)==1 and lead.coeff(Y,1)==0
          and lead.coeff(Y,0)!=0)
    root=S.sqrt(-lead.coeff(Y,0))
    check("outer nonzero constant root simple "+label,
          zero(lead.subs(Y,root)) and S.diff(lead,Y).subs(Y,root)!=0)
    do,dl=initial(S.diff(P,w).subs({v:s**2,w:ww}).subs(Y,root),s)
    check("outer entire derivative order "+label,do==2*derivative_order)
    check("outer derivative leading unit "+label,zero(dl-2*root))
    check("outer du keeps holomorphic coframe "+label,-do+1>=0)
    models.append({"kind":"outer","k":k,"L":L,"a":a,"d":d,
                   "nu":nu,"displacement_v_order":str(beta),
                   "wp_v_order":str(-derivative_order),
                   "formal_model_not_Keller_realization":True})

# Exact outer cancellation makes the SECOND root reveal a simple graph.
for k,L,d in [(2,3,1),(3,5,1),(3,5,2),(4,6,3)]:
    a=L+d;j=L+k
    W=w*(v**k*w+1)*(v**L*w+1)
    P=S.expand(W**2+v**a*w**3+v**d*w**2)
    wo=-v**(-L);wi=-v**(-k)
    label=f"k{k} L{L} d{d}"
    check("outer remainder cancellation exact "+label,zero(P.subs(w,wo)))
    check("inner root of regular cubic exact "+label,zero(W.subs(w,wi)))
    oo,ol=initial(S.diff(W,w).subs(w,wi),v)
    check("inner cubic derivative unit "+label,oo==0 and ol==-1)
    check("inner cubic remainder higher "+label,a-d-k==L-k>0)
    check("inner linear remainder margin "+label,k-d>0)
    beta=S.Rational(d,2)-k
    ww=-s**(-2*k)+s**(2*beta)*Y
    scaled=S.expand(s**(4*k-2*d)*(P.subs({v:s**2,w:ww})-p))
    order,lead=initial(scaled,s)
    check("inner entire scaled equation "+label,order==0 and zero(lead-Y**2-1))
    check("inner roots simple nonzero "+label,S.diff(lead,Y).subs(Y,S.I)==2*S.I)
    do,dl=initial(S.diff(P,w).subs({v:s**2,w:ww}).subs(Y,S.I),s)
    check("inner entire derivative order "+label,do==d-2*k)
    check("inner derivative leading unit "+label,zero(dl-2*S.I))
    check("inner coframe holomorphic after Puiseux "+label,-do+1>0)
    check("inner displacement smaller than centre "+label,beta+k>0)
    models.append({"kind":"inner_after_outer_cancellation","k":k,"L":L,"d":d,
                   "a":a,"wp_v_order":str(-beta),
                   "formal_model_not_Keller_realization":True})

# The exceptional boundary d=0 is excluded by the ENTIRE lower initial.
for k,L in [(1,2),(2,3),(3,5)]:
    W=w*(v**k*w+1)*(v**L*w+1)
    P=S.expand(W**2+v**L*w**3+w**2)
    label=f"k{k} L{L}"
    whole=S.expand(v**(2*k)*(P.subs(w,v**(-k)*Y)-p))
    order,lead=initial(whole,v)
    check("d0 entire lower initial "+label,
          order==0 and zero(lead-Y**2*((Y+1)**2+1)))
    root=-1+S.I
    check("d0 nonzero simple lower root "+label,
          zero(lead.subs(Y,root)) and S.diff(lead,Y).subs(Y,root)!=0)
    do,dl=initial(S.diff(P,w).subs(w,v**(-k)*root),v)
    check("d0 exact coframe order "+label,do==-k and dl!=0)

# Exact regular full polynomial-part square root in the surface algebra.
h,K,M,N,B,e,z,D,E,F=S.symbols("h K M N B e z D E F")
W=(h*w+K)*z+B*w
P=W**2+2*h*e*w*z+D*z+E*w+F
HH=W+e
RR=(D-2*e*K)*z+(E-2*e*B)*w+F-e**2
check("full regular cubic square identity",zero(P-HH**2-RR))
check("quadratic remainder ruling bound",
      S.degree(S.expand(RR.subs(z,(w**2-N)/M)),w)<=2)
check("full square root ruling degree",
      S.degree(S.expand(HH.subs(z,(w**2-N)/M)),w)==3)
x,sigma=S.symbols("x sigma")
check("central affine restriction of full root",
      zero(HH.subs({w:0,z:2*sigma*x})-(2*sigma*K*x+e)))

receipt={
    "claim":"PF2RULED171","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{name:sha(BASE/name) for name in [
        "RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md",
        "RULING_SIX_OUTER_SQUARES_AND_EVEN_SPECIALIZATIONS.md",
        "RULING_SIX_GLOBAL_QUARTIC_DIVISIBILITY_AND_CUBIC_REMAINDER.md",
        "RULING_SIX_QUARTIC_TRACE_AND_CONSTANT_SQUARE_ROOT_CORRECTION.md",
    ]},
    "checked_local_models":models,
    "uniform_proved_scope":{
        "next_coefficient_divisibility":"h divides C=a1-2B0*K",
        "regular_full_cubic_approximate_root":True,
        "regular_remainder_ruling_upper_bound":2,
        "full_root_in_target_field":False,
        "coefficientwise_regularity_in_C_u_w_claimed":False,
    },
    "whole_ruling_six_class_closed":False,
    "automatic_derivative_or_trace_integrality_proved":False,
    "full_goal_complete":False,"independent_whole_proof_review":False,
    "formal_verification":False,"external_review":False,
}
parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path)
args=parser.parse_args()
out=args.output or Path(tempfile.mkdtemp(prefix="math-check-"))/"ruling_six_regular_cubic_root_receipt.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),
                 "elapsed_seconds":receipt["elapsed_seconds"],"receipt":str(out)},indent=2))
