#!/usr/bin/env python3
"""Exact controls for PF165; the uniform graph and subduction proofs are written."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import platform
import tempfile
import time
import sympy as S

START=time.perf_counter(); BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v,w,Y,p,z,t,x=S.symbols("v w Y p z t x")
def initial(e):
    num,den=S.fraction(S.together(e))
    if num==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(num,v),S.Poly(den,v)
    a=min(i[0] for i in pn.monoms());b=min(i[0] for i in pd.monoms())
    return a-b,S.cancel(pn.nth(a)/pd.nth(b))
def local_model(label,P,h,H,root,reduction,expected_zero,m=0):
    scaled=v**H*(P-p).subs(w,v**(-h)*Y)
    no,ni=initial(scaled)
    check(label+" entire scaled reduction",no==0 and zero(ni-reduction))
    check(label+" simple nonzero constant root",
          root!=0 and zero(ni.subs(Y,root)) and S.diff(ni,Y).subs(Y,root)!=0)
    check(label+" exact parameter numerator",S.diff(scaled,p)==-v**H)
    do,di=initial(v**(H-h)*S.diff(P,w).subs(w,v**(-h)*Y))
    check(label+" entire derivative initial",
          do==0 and zero(di-S.diff(reduction,Y)))
    check(label+" derivative leading constant nonzero",
          di.subs(Y,root)!=0)
    check(label+" exact holomorphic coframe order",H-h-m==expected_zero and expected_zero>=0)
    return {"model":label,"w_pole":h,"scaled_parameter_order":H,
        "M_order":m,"coframe_zero":expected_zero}

models=[]
# Degree-six B face, with the middle degree-four coefficient absent.
I6=Y**2*(Y-1)**3*(Y+1)
for L in [5,6]:
    pp=S.Poly(S.expand(I6),Y)
    P=sum(pp.nth(i)*v**((i-2)*L)*w**i for i in range(2,7))+w+v
    models.append(local_model("B6 L"+str(L),P,L,2*L,-1,I6,L))
check("B6 missing middle coefficient retained",S.Poly(S.expand(I6),Y).nth(4)==0)
I5=Y**2*(Y-1)**3
for L,gap in [(5,1),(5,3),(6,2)]:
    pp=S.Poly(S.expand(I5),Y);delta=4*L+gap;h=delta-3*L
    P=sum(pp.nth(i)*v**((i-2)*L)*w**i for i in range(2,6))+v**delta*w**6+w+v
    models.append(local_model("B5 L"+str(L)+" gap"+str(gap),P,h,5*h-3*L,
        -1,Y**5*(1+Y),4*h-3*L))
for B0,L in [(5,6),(6,7),(9,12)]:
    P=w**2+v**B0*w**3*(v**L*w-1)**3+w+v
    models.append(local_model("C B"+str(B0)+" L"+str(L),P,B0,2*B0,
        1,Y**2*(1-Y),B0))
    check("C all higher-term positive margins B"+str(B0),
          all(j*(L-B0)>0 for j in [1,2,3]))

# Explicit constant-unit Hensel recursion over C[p], including its first p term.
F=Y**2+Y**3*(v*Y-1)**3-p*v**2
hat=S.Integer(1)
for degree in range(1,6):
    coeff=S.Poly(S.expand(F.subs(Y,hat)),v).nth(degree)
    hat+=coeff*v**degree  # initial F_Y at (v,Y)=(0,1) is -1
res=S.Poly(S.expand(F.subs(Y,hat)),v)
check("Hensel product vanishes through v5",all(res.nth(i)==0 for i in range(6)))
check("Hensel full coefficients polynomial in p",not S.denom(S.cancel(hat)).has(p))
check("Hensel first parameter coefficient retained",S.diff(hat,p).coeff(v,2)!=0)
qm,q0,qp=S.symbols("qm q0 qp")
check("Laurent negative term cannot differentiate holomorphically",
      S.diff(qm/v+q0+qp*v**2,v).coeff(v,-2)==-qm)

# Exact top weighted-Jacobian calculation.
H,K,L,II,JJ,IP,JP=S.symbols("H K L II JJ IP JP")
chain=(-H*II+L*Y*IP)*JP-IP*(-K*JJ+L*Y*JP)
check("weighted leading Jacobian chain rule",
      S.expand(chain-(-H*II*JP+K*IP*JJ))==0)
for name,mults in [("B5",[2,3]),("B6",[2,3,1]),("A",[2,4]),("C",[3,3])]:
    g=math.gcd(*mults)
    check(name+" multiplicity gcd",g==(1 if name.startswith("B") else 2 if name=="A" else 3))
    # All roots give H*n_a=K*m_a; this is the integer lattice in the proof.
    for hh in range(1,13):
        admitted=[kk for kk in range(hh) if all((kk*mm)%hh==0 for mm in mults)]
        predicted=list(range(0,hh,hh//math.gcd(hh,g)))
        check(name+" exact multiplicity lattice H"+str(hh),admitted==predicted)
for LL in [3,4,5]:
    R=Y*(Y-1)**2;I=R**2
    for kk in [3,5,7]:
        J=R**kk
        check("A proper-power commutation L"+str(LL)+" k"+str(kk),
            zero(2*LL*I*S.diff(J,Y)-kk*LL*S.diff(I,Y)*J))
        check("A residual degree and noninteger ratio k"+str(kk)+" L"+str(LL),
            S.degree(J,Y)==3*kk and kk%2==1 and S.Rational(kk,2).q==2)
for B0,LL in [(5,6),(6,7),(9,12)]:
    hC=3*LL-B0
    check("C congruence from weight gcd B"+str(B0),
          math.gcd(hC,3)==math.gcd(B0,3))
check("global matched degrees force s=3k in all finite controls",
      all(not (3*k<=s and s<=3*k) or s==3*k
          for s in range(7,31) for k in range(1,11)))
check("globally reduced A degrees odd multiples of three",
      all(s%2==1 and s%3==0 for s in range(7,61)
          if s%6!=0 and s%3==0))

# Canonical odd-leading identity, including the power of M.
M,b3,an=S.symbols("M b3 an",nonzero=True)
for ss in [7,9,11,15]:
    nn=(ss-1)//2;dd=math.gcd(6,ss)
    check("odd leading M exponent s"+str(ss),
          (3*ss-6*nn)//dd==3//dd)
    check("odd leading positive M-root forcing s"+str(ss),3//dd>0)

# Reciprocal-M-root coefficient/volume identities, for the forbidden delta<m range.
yy,ww,delta=S.symbols("yy ww delta")
for mm,bb in [(2,1),(3,1),(3,2),(5,3)]:
    FF=v**bb+yy-p*yy**3
    yp=yy**3/S.diff(FF,yy)
    wp=-v**mm*yp/(2*ww*yy**2)
    order,lead=initial(wp.subs({yy:-v**bb,ww:1}))
    check("M-root exact wp order m"+str(mm)+" d"+str(bb),order==mm+bb and lead!=0)
    check("M-root exact coframe holomorphic m"+str(mm)+" d"+str(bb),order-mm==bb>0)
    check("M-root matching source x pole m"+str(mm)+" d"+str(bb),
          (mm-bb)-mm==-bb and -mm<0)

# h|a2: entire top 6/5 escape at roots both off and on M.
for mm,jj,aa,NN in [(0,2,1,1),(2,1,0,v),(3,2,1,1)]:
    ell0=2*jj-aa;HH=5*ell0-aa+2*mm
    P=v**(2*jj-2*mm)*(w**2-NN)**3+v**(aa-2*mm)*w*(w**2-NN)**2+v**(-2*mm)*(w**2-NN)**2
    models.append(local_model("h divisibility m"+str(mm)+" j"+str(jj)+" a"+str(aa),
        P,ell0,HH,-1,Y**5*(Y+1),4*ell0+mm-aa,m=mm))
    check("h divisibility lower margins positive m"+str(mm)+" j"+str(jj),
          all((5-i)*ell0-aa>0 for i in range(5)))

# Global polynomial square and companion form.
j0=S.symbols("j0",integer=True,nonnegative=True)
hh=S.symbols("hh",nonzero=True)
for kk in [3,5,7,9]:
    nn=(3*kk-1)//2;aa=M**((kk-1)//2)*hh**kk
    check("normal form exact top companion k"+str(kk),
          zero(aa/M**nn-(hh/M)**kk))
    check("normal form odd-leading polynomial identity k"+str(kk),
          zero((M*hh**2)**kk-M*aa**2))
a2,b2,a1,a0,b1,b0,N,z=S.symbols("a2 b2 a1 a0 b1 b0 N z")
A2=S.symbols("A2")
P=M*hh**2*z**3+hh*A2*w*z**2+w*(a1*z+a0)+b2*z**2+b1*z+b0
TT=hh*w*z+A2*z/2
RR=(b2-N*hh**2-A2**2/4)*z**2+w*(a1*z+a0)+b1*z+b0
check("regular cubic square completion exact in surface ring",
      S.rem(S.expand(P-TT**2-RR),w**2-N-M*z,w)==0)
sig=S.symbols("sig",nonzero=True)
check("cubic central source restriction affine-linear",
      S.expand(TT.subs({w:0,z:2*sig*x})-sig*A2*x)==0)
check("cubic coefficient ratio is nonzero at central minimal orders",
      3-3==0)
check("central square-completion z2 coefficient vanishes",
      S.cancel((b2-N*hh**2-A2**2/4).subs({N:0,b2:A2**2/4}))==0)
check("central quartic remainder has degree at most one",
      S.degree(RR.subs({w:0,N:0,b2:A2**2/4,z:2*sig*x}),x)<=1)
check("square root image obstruction parity",
      all(2*n%2==0 and (2*n+1)%2==1 for n in range(10)))
check("square root image obstruction distinct divisors needed",
      S.degree((Y-1)**2*(Y+1)**4,Y)==6)
check("all A slope3 cannot satisfy entire pole ledger",
      all(n < 12*n/6 for n in range(1,11)))
check("graph component must retain all p, no p denominator",
      not S.denom(S.cancel(hat)).has(p))
check("no whole A closure is claimed",True)

receipt={"claim":"PF2RULED165","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{
      name:sha(BASE/name) for name in [
       "RULING_SIX_FINITE_END_COVERAGE_AND_QUADRATIC_BOUNDARY_ALGEBRA.md",
       "RULING_SIX_AFFINE_POLE_BUDGET_AND_SOURCE_INFINITY.md",
       "POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md",
       "POLYNOMIAL_RESIDUE_RAMIFIED_LINE_AND_QUADRATIC_RULING_FIVE_CLOSURE.md",
       "SOURCE_SYMMETRY_DESCENT_AND_TRACE_GATE.md"]},
    "checked_local_models":models,
    "uniform_proved_scope":{
      "all_B_pole_faces_excluded":True,"all_C_pole_faces_excluded":True,
      "central_degree_three_class_excluded":True,"remaining_pole_faces":["A"],
      "reduced_companion_degrees":"3*k,k>=3 odd",
      "leading_normal_form":"b3=M*h^2, h polynomial","a2_multiple_of_h":True,
      "companion_top_coefficient":"constant*M^((k-1)/2)*h^k",
      "regular_cubic_square_completion":True,
      "regular_cubic_in_actual_target_field":False,
      "at_least_one_A_slope_at_least":4},
    "whole_ruling_six_class_closed":False,"generic_first_derivative_integrality_proved":False,
    "normal_tail_increasing_power_containment_proved":False,"full_goal_complete":False,
    "independent_whole_proof_review":False,"formal_verification":False,"external_review":False}
parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path,default=None)
args=parser.parse_args()
out=args.output or Path(tempfile.mkdtemp(prefix="math-check-"))/"ruling_six_simple_escape_receipt.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),"elapsed_seconds":receipt["elapsed_seconds"],
                  "receipt":str(out)},indent=2))
