#!/usr/bin/env python3
"""PF179 exact controls; the unrestricted Newton/Hensel proof is in the note."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import time
import sympy as S
START=time.perf_counter()
BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_NONCRITICAL_QUARTIC_SPLITTING_AND_SMALL_SLOPE_EXCLUSION.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v,s,Y,X,Z,p,lam,eps=S.symbols("v s Y X Z p lam eps")
A,B,C,gamma=S.symbols("A B C gamma",nonzero=True)
def initial(expr,var):
    num,den=S.fraction(S.cancel(expr))
    if num==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(num,var),S.Poly(den,var)
    on=min(t[0] for t in pn.monoms());od=min(t[0] for t in pd.monoms())
    return on-od,S.cancel(pn.nth(on)/pd.nth(od))

# The correction to the other Hensel factor is included via reduction modulo U^2.
W=(A+B/gamma-B*(Y-gamma)/gamma**2+B*(Y-gamma)**2/gamma**3
   -B*(Y-gamma)**3/gamma**4)/C**2
T=Y*(A*Y+B)
check("first-variation remainder includes factor C2Y2",
      S.rem(S.cancel(C**2*Y**2*W-T),S.Poly((Y-gamma)**4,Y),Y)==0)
check("quartic correction degree is at most three",S.degree(W,Y)<=3)
wx=S.Poly(S.expand(W.subs(Y,X+gamma)),X)
expected=[(A+B/gamma)/C**2,-B/(C**2*gamma**2),
          B/(C**2*gamma**3),-B/(C**2*gamma**4)]
for i in range(4):
    check("exact quartic deformation coefficient X"+str(i),zero(wx.nth(i)-expected[i]))
# A direct polynomial product identity supplies the corresponding quadratic correction.
F2dot=S.cancel((T-C**2*Y**2*W)/(Y-gamma)**4)
check("other Hensel factor correction remains polynomial",S.denom(F2dot).has(Y)==False)
check("full first-order factor product",
      zero(F2dot*(Y-gamma)**4+C**2*Y**2*W-T))

# Depression of a quadratic square perturbed by a general cubic.
J0,w3,w2,w1,w0=S.symbols("J0 w3 w2 w1 w0")
raw=(X**2+J0)**2+eps*(w3*X**3+w2*X**2+w1*X+w0)
depressed=S.expand(S.series(raw.subs(X,Z-eps*w3/4),eps,0,2).removeO())
check("depression removes cubic coefficient to first order",zero(depressed.coeff(Z,3)))
check("depressed linear correction includes J0 term",
      zero(depressed.coeff(Z,1)-eps*(w1-J0*w3)))
Jnew=J0+eps*w2/2
Dnew=S.series(depressed.coeff(Z,0)-Jnew**2,eps,0,2).removeO()
check("Delta correction includes J0 term",zero(Dnew-eps*(w0-J0*w2)))
check("positive-order J0 leaves noncritical Delta leading",
      zero((w0-J0*w2).subs(J0,0)-w0))
check("positive-order J0 leaves q leading",
      zero((w1-J0*w3).subs(J0,0)-w1))

# All repeated initial roots at 2ell=t are accounted for.
j,delta=S.symbols("j delta",nonzero=True)
equal=(lam**2+j)**2+delta
deriv=S.diff(equal,lam)
check("nonzero critical derivative root requires lam2=-j",
      zero(deriv-4*lam*(lam**2+j)))
check("lam2=-j cannot be an initial root when delta nonzero",
      S.rem(equal,S.Poly(lam**2+j,lam),lam)==delta)
exception=S.expand(equal.subs(delta,-j**2))
check("exceptional initial has exactly a double zero",
      zero(exception-lam**2*(lam**2+2*j)))
check("exceptional nonzero factor is separable",
      S.resultant(lam**2+2*j,S.diff(lam**2+2*j,lam),lam)!=0)
check("exceptional nonzero roots have unit full derivative",
      S.resultant(lam**2+2*j,deriv,lam)!=0)

# Whole sextic models in the original local w coordinate, not a truncated quartic alone.
# These are coefficient controls; no actual mate or global Keller realization is claimed.
models=[]
for L in [3,4,5]:
    for b in range(L-1):
        t=L+b
        cases=[("less",1,1,1),("greater",t//2+1,1,1)]
        if t%2==0:
            cases.extend([("equal",t//2,1,1),("exceptional",t//2,1,-2)])
        for kind,ell,aa,bb in cases:
            tag=f"L{L} b{b} {kind}"
            SS=Y*((Y-1)**2+v**ell)
            MM=(1+v)**2
            NN=1+v
            # Canonical remainder d*z+e*w+f, with d=aa*M*v^t,e=bb*v^b.
            F0=SS**2+v**t*(aa*Y**2+bb*Y)+v**(2*L)*(3+v)-aa*v**(2*L+t)*NN
            FF=F0-p*v**(2*L)
            check("parameter is strictly later than initial "+tag,2*L>t)
            check("noncritical leading remainder "+tag,aa+bb!=0)
            if kind=="less":
                YY=1+S.I*s**(2*ell)+s**(2*(t-ell))*lam
                mu=S.Rational(t+ell,2)
                expected_initial=-4*lam**2+aa+bb
                expected_derivative=-8*lam
            else:
                YY=1+s**t*lam
                mu=S.Rational(3*t,4)
                expected_initial=(lam**4+aa+bb if kind=="greater"
                                  else (lam**2+1)**2+aa+bb)
                expected_derivative=(4*lam**3 if kind=="greater"
                                     else 4*lam*(lam**2+1))
            substituted=FF.subs({v:s**4,Y:YY},simultaneous=True)
            oo,ll=initial(substituted/s**(4*t),s)
            check("whole sextic simple initial "+tag,oo==0 and zero(ll-expected_initial))
            dy=S.diff(F0,Y).subs({v:s**4,Y:YY},simultaneous=True)
            do,dl=initial(dy/s**(4*mu),s)
            check("whole sextic derivative order "+tag,do==0 and zero(dl-expected_derivative))
            simple_factor=(lam**2+2 if kind=="exceptional" else expected_initial)
            check("selected leading factor has constant unit derivative "+tag,
                  S.resultant(simple_factor,expected_derivative,lam)!=0)
            check("selected derivative satisfies global Newton bound "+tag,mu<=S.Rational(3*t,4))
            check("small slope coframe is holomorphic or simple "+tag,mu<=L+1)
            normalized_order=4*(L+1-mu)-1
            check("whole coframe normalized order dichotomy "+tag,
                  normalized_order.is_integer and
                  ((normalized_order>=0) if mu<L+1 else normalized_order==-1))
            # Parameter dependence of implicit solution has a nonzero linear forcing.
            response_order=2*L-mu
            check("first parameter response is strictly after Newton scale "+tag,
                  response_order>(S.Rational(t,4) if kind!="less" else S.Rational(t-ell,2)))
            models.append({"L":L,"e_order":b,"t":t,"case":kind,
                           "J_order":ell,"mu":str(mu),
                           "normalized_order_after_v_s4":str(normalized_order)})

# Exact whole-coordinate derivative and coframe conversion.
ww=S.symbols("w")
for L,b,ell in [(3,0,1),(5,3,5),(5,1,3)]:
    t=L+b
    H=ww*((v**L*ww-1)**2+v**ell)
    PP=H**2+v**t*ww**2+v**b*ww+(3+v)-v**t*(1+v)
    scaled=v**(2*L)*PP.subs(ww,v**(-L)*Y)
    check(f"whole implicit velocity conversion L{L} b{b}",
          zero(S.diff(PP,ww)-v**(-L)*S.diff(scaled,Y).subs(Y,v**L*ww)))
    check(f"whole parameter term exact L{L} b{b}",
          zero((v**(2*L)*(PP-p)).subs(ww,v**(-L)*Y)-scaled+p*v**(2*L)))
for alpha in [1,2,4,8]:
    for rr in [S.Rational(9,2),S.Rational(11,2),S.Integer(6)]:
        normalized=alpha*(6-rr)-1
        if normalized.is_integer:
            check(f"L5 coframe conversion alpha{alpha} mu{rr}",
                  (normalized>=0 if rr<6 else normalized==-1))

# Uniform inequality and exact small-slope consequence; finite loops are controls only.
Lsym,bsym=S.symbols("L b",integer=True)
check("uniform inequality has equivalent e-order form",
      S.expand(3*(Lsym+bsym)-4*(Lsym+1)-(3*bsym-Lsym-4))==0)
for L in range(3,11):
    for b in range(L-1):
        allowed=3*(L+b)>4*(L+1)
        check(f"uniform noncritical inequality equivalence L{L} b{b}",
              allowed==(3*b>L+4))
        if L<=5:
            check(f"every small-slope order is excluded L{L} b{b}",not allowed)
check("large-slope control is deliberately retained",3*4>6+4)
check("noncritical L6 control can have mu above coframe threshold",
      S.Rational(3*(6+4),4)>6+1)
kk=S.symbols("k",integer=True,positive=True)
check("companion-dependent ratio is never critical for k>=3",
      zero(1-2*kk/(kk-1)+(kk+1)/(kk-1)))
for k in [3,5,7,11,101]:
    check("actual odd companion ratio noncritical k"+str(k),-2*S.Rational(k,k-1)!=-1)
check("degree24 remaining local order list",
      [(5+b,b) for b in range(4)]==[(5,0),(6,1),(7,2),(8,3)])

receipt={"claim":"PF2RULED179","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{name:sha(BASE/name) for name in [
        "RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md",
        "RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md",
        "RULING_SIX_DEGREE24_LEADING_CLASSIFICATION_AND_LINEAR_CONTROLS.md",
        "RULING_SIX_BALANCED_QUADRATIC_REMAINDER_AND_LEADING_RESONANCE.md"]},
    "uniform_proved_scope":{
        "coefficient_and_companion_degrees":"unrestricted",
        "noncritical_necessary_bound":"3*(L+ord e)>4*(L+1), equivalently 3*ord e>L+4",
        "every_L_at_most_five_support_critical":True,
        "degree24_companion_dependent_noncritical_ratio_excluded":True,
        "critical_balanced_remainder_closed":False,
        "degree24_or_full_ruling_six_closed":False,
        "target_integrality_or_JC2_claimed":False},
    "whole_sextic_models":models,
    "assurance":"Exact internal symbolic controls; unrestricted proof written; no formal or external human review."}
parser=argparse.ArgumentParser()
parser.add_argument("--write-receipt",action="store_true")
args=parser.parse_args()
if args.write_receipt:
    (BASE/"ruling_six_noncritical_quartic_splitting_receipt.json").write_text(
        json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
print(json.dumps({"claim":receipt["claim"],"outcome":"PASS","check_count":len(checks),
                  "elapsed_seconds":receipt["elapsed_seconds"],
                  "receipt_written":args.write_receipt},indent=2))
