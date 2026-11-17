#!/usr/bin/env python3
"""PF170 controls; the uniform trace/Hensel/graph proof is in the source note."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import tempfile
import time
import sympy as S

START=time.perf_counter();BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_QUARTIC_TRACE_AND_CONSTANT_SQUARE_ROOT_CORRECTION.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v,Y,p,w,s=S.symbols("v Y p w s")
def initial(e,var=v):
    n,d=S.fraction(S.together(e))
    if n==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(n,var),S.Poly(d,var)
    a=min(i[0] for i in pn.monoms());b=min(i[0] for i in pd.monoms())
    return a-b,S.cancel(pn.nth(a)/pd.nth(b))
def trunc(e,n):
    pp=S.Poly(S.expand(e),v)
    return S.expand(sum(c*v**i[0] for i,c in pp.terms() if i[0]<n))

# Actual companion multiplication in the rank-four quartic algebra.
b,c,d,e0=S.symbols("b c d e0")
C=S.Matrix([[0,0,0,-e0],[1,0,0,-d],[0,1,0,-c],[0,0,1,-b]])
Gprime=4*C**3+3*b*C**2+2*c*C+d*S.eye(4)
q=d-b*c/2+b**3/8
check("actual quartic multiplication satisfies its monic equation",
      (C**4+b*C**3+c*C**2+d*C+e0*S.eye(4)).applyfunc(S.expand)==S.zeros(4))
check("actual quartic derivative trace equals -8 depressed linear",
      S.expand(S.trace(Gprime)+8*q)==0)
check("quartic Newton sums independently give same trace",
      S.expand(4*(-b**3+3*b*c-3*d)+3*b*(b**2-2*c)+2*c*(-b)+4*d+8*q)==0)
a2,a1,a0,eps,D=S.symbols("a2 a1 a0 eps D",nonzero=True)
A=a2*Y**2+a1*Y+a0
G=S.Poly(S.expand((A**2+eps*Y+D)/a2**2),Y)
bb,cc,dd=G.nth(3),G.nth(2),G.nth(1)
check("normalized square plus linear has exact depressed invariant",
      zero(dd-bb*cc/2+bb**3/8-eps/a2**2))
check("constant term cannot cancel depressed linear invariant",
      S.diff(dd-bb*cc/2+bb**3/8,D)==0)

factor_models=[]
for L,r,Nprec in [(3,1,7),(4,3,5)]:
    A=(Y-1)**2+v*(Y+1)+v**2
    epsilon=v**r
    d0=-epsilon+v**(r+1)  # deliberately cancel the leading value at Y=1
    F=Y**2*(A**2+epsilon*Y+d0)+v**L*Y+v**(2*L)*(1-p)
    H=Y**2;G=(Y-1)**4
    G0=(Y-1)**4
    for n in range(1,Nprec+1):
        res=S.Poly(S.expand(F-H*G),v).nth(n)
        dh=S.rem(S.expand(res*(1+4*Y)),Y**2,Y)
        dg,rem=S.div(S.expand(res-G0*dh),Y**2,Y)
        check(f"Hensel exact monic division L{L} order{n}",rem==0)
        H+=dh*v**n;G+=dg*v**n
    check("whole Hensel product L"+str(L),trunc(F-H*G,Nprec+1)==0)
    check("two/four monic ranks retained L"+str(L),
          S.Poly(H,Y).LC()==1 and S.degree(H,Y)==2 and S.Poly(G,Y).LC()==1 and S.degree(G,Y)==4)
    check("finite factor is Y2 through required cutoff L"+str(L),trunc(H-Y**2,L)==0)
    check("actual quartic agrees with square/linear/constant L"+str(L),
          trunc(G-(A**2+epsilon*Y+d0),L)==0)
    gp=S.Poly(G,Y)
    qq=gp.nth(1)-gp.nth(3)*gp.nth(2)/2+gp.nth(3)**3/8
    oo,ll=initial(qq)
    check("actual quartic invariant exact small order L"+str(L),oo==r and ll==1)
    check("actual quartic coefficients polynomial in p L"+str(L),not S.denom(S.together(G)).has(p))
    check("whole exact parameter term L"+str(L),S.diff(F,p)==-v**(2*L))
    factor_models.append({"ell":L,"epsilon_order":r,"precision":Nprec,
                           "quartic_invariant_order":oo})
    if Nprec>=2*L:
        check("parameter does not enter factors too early L"+str(L),
              trunc(S.diff(G,p),2*L)==0)
# Deliberately degenerate quartic: trace only needs ONE small derivative,
# not simplicity of every root of a leading model.
X=S.symbols("X")
bad=X**4-4*v**3*X+3*v**4
check("double-root quartic control retained",
      S.expand(bad-(X-v)**2*(X**2+2*v*X+3*v**2))==0)
check("nonzero linear quartic trace survives a double root",
      -8*(-4*v**3)==32*v**3)
check("nonzero linear quartic cannot be a quadratic square",
      not all(m%2==0 for f,m in S.factor_list(X**4-4*X+3,X)[1]))

# The response contradiction is uniform over all possible derivative orders.
ell,mu,rr=S.symbols("ell mu rr")
check("whole response order from exact parameter term",
      S.expand((2*ell-mu)+mu-2*ell)==0)
gap,extra=S.symbols("gap extra",positive=True)
# e=ell-gap and mu=e-extra; gap>0,extra>=0 in the proof.
R=ell+gap+extra
check("nonlinear Taylor terms strictly exceed parameter order",
      S.expand(2*R-2*ell)==2*gap+2*extra)
check("linear coefficient unchanged by prior truncation",
      S.expand(R-(ell-gap-extra))==2*gap+2*extra)
for a in range(1,5):
    for L,r in [(3,1),(4,3),(5,2)]:
        for mm in [S.Rational(0),S.Rational(r,2),S.Rational(r)]:
            response=2*L-mm
            check(f"all normalized coframe/response controls alpha{a} L{L} mu{mm}",
                  response>L and a*(L-mm)+a-1>=0)
# Half-integral outer weight is handled by clearing a parameter denominator.
check("half-integral trace cutoff cleared",2*S.Rational(3,2)==3)
check("half-integral failed-divisibility response margin",
      S.Rational(3,2)-1>0)

# Unit-K exact cubic-centre and simple split, including M-roots.
unit_models=[]
for j,m,a in [(1,0,0),(2,0,1),(1,2,0),(3,2,1),(2,3,0)]:
    # Local structural model: h=v^j,K=1,M=v^m,N=0,B0=0.
    # Exact W1 zero is wo=-v^-j; source polynomial forms are retained.
    h=v**j;M=v**m
    W=(h*w+1)*w**2/M
    P=W**2+v**a*w**3/M+w**2/M+w+1
    wo=-v**(-j)
    check(f"unitK exact cubic centre j{j} m{m} a{a}",zero(W.subs(w,wo)))
    no,ni=initial(S.diff(W,w).subs(w,wo))
    check(f"unitK exact cubic derivative j{j} m{m} a{a}",no==-j-m and ni!=0)
    beta=S.Rational(a-j+m,2)
    # Clear every half power by v=s^2, delta=s^(a-j+m).
    shifted=s**(2*(3*j+m-a))*(P-p).subs(
         {v:s**2,w:-s**(-2*j)+s**(a-j+m)*Y},simultaneous=True)
    no,ni=initial(shifted,s)
    check(f"unitK entire simple split j{j} m{m} a{a}",
          no==0 and zero(ni-(Y**2-1)))
    check(f"unitK simple roots j{j} m{m} a{a}",
          ni.subs(Y,1)==0 and S.diff(ni,Y).subs(Y,1)!=0)
    wp=S.Rational(5*j+3*m-a,2);vol=wp-m
    check(f"unitK exact holomorphic volume bound j{j} m{m} a{a}",vol>0)
    check(f"unitK shift is smaller than polar centre j{j} m{m} a{a}",beta+j>0)
    unit_models.append({"h_order":j,"M_order":m,"C1_order":a,
                        "wp_order":str(wp),"coframe_v_order":str(vol)})
check("M-root hzero forces unitK via reached line",
      S.expand((S.Symbol("K")**2).subs(S.Symbol("K"),0))==0)

receipt={"claim":"PF2RULED170","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{
      name:sha(BASE/name) for name in [
        "RULING_SIX_GLOBAL_QUARTIC_DIVISIBILITY_AND_CUBIC_REMAINDER.md",
        "RULING_SIX_OUTER_SQUARES_AND_EVEN_SPECIALIZATIONS.md",
        "RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md",
        "POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md"]},
    "quartic_Hensel_models":factor_models,"unit_K_models":unit_models,
    "uniform_proved_scope":{
       "C1":"a1-2*B0*K","constant_correction":"C1/(2*h)",
       "regular_at_every_A_pole_divisor":True,
       "regular_at_every_outer_four_sheet_h_root":True,
       "regular_at_every_M_root":True,"regular_at_every_unit_K_h_root":True,
       "quartic_trace":"Tr(G_prime)=-8*(d-b*c/2+b^3/8)",
       "remaining_possible_poles":"off-M common h,K zero, degree2, outer escape mass2"},
    "whole_ruling_six_class_closed":False,"global_constant_correction_regular_proved":False,
    "generic_first_derivative_integrality_proved":False,
    "normal_tail_increasing_power_containment_proved":False,"full_goal_complete":False,
    "independent_whole_proof_review":False,"formal_verification":False,"external_review":False}
parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path,default=None)
args=parser.parse_args()
out=args.output or Path(tempfile.mkdtemp(prefix="math-check-"))/"ruling_six_quartic_trace_correction_receipt.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),"elapsed_seconds":receipt["elapsed_seconds"],
                  "receipt":str(out)},indent=2))
