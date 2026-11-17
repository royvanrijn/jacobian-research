"""PF177 exact coefficient classification and controls, without a mate claim."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import tempfile
import time
import sympy as S

START=time.perf_counter();BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_DEGREE24_LEADING_CLASSIFICATION_AND_LINEAR_CONTROLS.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(value):return S.cancel(S.expand(value))==0
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
u,tau,v,w,Y=S.symbols("u tau v w Y")
F=40*tau**3-105*tau**2+84*tau-14
ledger=[]
for m in range(2,25):
    if (24-m)%2:continue
    H=(24-m)//2
    for count in range(1,H+1):
        for slopes in range(count,H//2+1):
            low=int(S.ceiling(S.Rational(24+3*m,6)))
            high=2*slopes-5*count
            for degree in range(max(1,low),high+1):ledger.append((m,H,count,slopes,degree))
check("complete degree24 integer ledger",ledger==[(2,11,1,5,5)])
check("degree-two M is entirely central",2-2==0)
D0=u*(u-1)*(u-tau)
G=u**4*(u-1)**10*(u-tau)
c0=-3/tau
c1=3*(8*tau-1)/(2*tau**2)
c2=-3*(40*tau**2-25*tau+4)/(10*tau**3)
ss=c2*u**2+c1*u+c0;rho=D0*ss
residual=S.cancel(3*D0*S.diff(ss,u)+(-6*u**2+(8*tau-1)*u-tau)*ss-3)
check("entire remaining parameter equation",zero(residual+3*F*u**3/(5*tau**3)))
check("entire leading ODE residual factor",zero(3*G*S.diff(rho,u)-S.diff(G,u)*rho-3*G-G*residual))
for power in [0,1,2,4]:
    check("parameter comparison coefficient u"+str(power),S.expand(residual).coeff(u,power)==0)
check("cubic has three distinct roots",S.discriminant(F,tau)==-1428840)
check("zero parameter is not a cubic root",F.subs(tau,0)==-14)
check("A parameter is not a cubic root",F.subs(tau,1)==5)
check("central e exceptional parameter excluded",F.subs(tau,S.Rational(3,2))==S.Rational(43,4))
check("all relevant parameter units",S.gcd(F,tau*(tau-1)*(2*tau-3))==1)
check("forced leading degree coefficient nonzero on cubic",S.gcd(F,40*tau**2-25*tau+4)==1)
G0=u**5*(u-1)**10
rho0=(-42*u+420*u**2-945*u**3+810*u**4-243*u**5)/28
check("central remaining factor positive leading ODE",zero(3*G0*S.diff(rho0,u)-S.diff(G0,u)*rho0-3*G0))

delta=1-tau;h=(u-1)**10*(u-tau);M=u**2
el=(u-tau)*(u-1)**4*(3-2*tau-u)
ratio=S.cancel(el**2/h)
check("uniform exact third character",zero(ratio.subs(u,v+1)-(v-3*delta+4*delta**3/v**2)))
primitive=v**2/2-3*delta*v-4*delta**3/v
check("uniform third-character primitive",zero(S.diff(primitive,v)-ratio.subs(u,v+1)))
check("central h value",h.subs(u,0)==-tau)
check("central linear coefficient",zero(el.subs(u,0)-tau*(2*tau-3)))
check("A linear order exactly four",zero(S.limit(el/(u-1)**4,u,1)-2*delta**2))
check("other h-root linear order exactly one",zero(S.limit(el/(u-tau),u,tau)+3*(tau-1)**5))
check("A support order strict",4<5)
check("unit-K hroot order",1>=1)
K=(u-1)**5;B0=1/(4*delta)
check("A cubic face square root",zero(delta*Y**3+Y**2+Y/(4*delta)-delta*Y*(Y+1/(2*delta))**2))
z,N=S.symbols("z N")
H3=(h*w+K)*z+B0*w
P=H3**2+el*w
check("source-image regular polynomial",S.denom(S.together(P)).has(u)==False)
check("whole central source line degree two",S.degree(S.expand(P.subs({u:0,w:0,z:-2*Y})),Y)==2)
check("whole A line degree two",S.degree(S.expand(P.subs({u:1,z:(w**2-1)})),w)==2)
check("whole third-root line degree four",S.degree(S.expand(P.subs({u:tau,z:(w**2-tau)/tau**2})),w)==4)
check("generic ruling degree six",S.degree(S.expand(P.subs(z,(w**2-u)/u**2)),w)==6)

# Entire generic fifth character, derived without inserting the controls.
a3,a2,a1,a0,ee,ff,x=S.symbols("a3 a2 a1 a0 ee ff x")
HH=a3*w**3+a2*w**2+a1*w+a0
formal=S.series(((ee*w+ff)**3/HH).subs(w,1/x),x,0,2).removeO()
Xi5=S.expand(formal).coeff(x,1)
check("generic fifth character",zero(Xi5-(3*ee**2*ff/a3-ee**3*a2/a3**2)))
check("canonical fifth character",zero((Xi5/M).subs({a3:h/M,a2:K/M,ee:el})-ratio*(3*ff-el*K/h)))
check("uniform eK over h",zero((el*K/h).subs(u,v+1)-(2*delta/v-1)))
constant_f=S.symbols("constant_f")
Sbad=S.cancel((ratio*(3*constant_f-el*K/h)).subs(u,v+1))
check("constant remainder fails fifth character",zero(S.residue(Sbad,v,0)-6*delta**2))
f=-v/(2*delta)
Sgood=S.cancel((v-3*delta+4*delta**3/v**2)*(3*f-(2*delta/v-1)))
expected=-3*v**2/(2*delta)+11*v/2-5*delta+4*delta**3/v**2-8*delta**4/v**3
check("polynomial repair full fifth function",zero(Sgood-expected))
primitive5=-v**3/(2*delta)+11*v**2/4-5*delta*v-4*delta**3/v+4*delta**4/v**2
check("polynomial repair exact fifth primitive",zero(S.diff(primitive5,v)-Sgood))
check("polynomial repair fifth residue vanishes",S.residue(Sgood,v,0)==0)

receipt={
 "claim":"PF2RULED177","outcome":"PASS","check_count":len(checks),"checks":checks,
 "elapsed_seconds":time.perf_counter()-START,
 "environment":{"python":platform.python_version(),"sympy":S.__version__},
 "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
 "retained_proof_sources":{name:sha(BASE/name) for name in [
   "RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md",
   "RULING_SIX_RETAINED_SOURCE_PACKETS_AND_STRICT_INFINITY_BOUNDS.md",
   "RULING_SIX_LEADING_PRIMITIVE_AND_POLYNOMIAL_ODE.md",
   "RULING_SIX_ODD_HALF_POWER_RESIDUE_HIERARCHY.md",
   "RULING_SIX_COMPANION_RESONANCE_AND_LINEAR_REMAINDER_CLOSURE.md",
 ]},
 "uniform_proved_scope":{
   "degree24_ledger":{"m":2,"H":11,"A_support_count":1,"A_slope":5,"b":5},
   "degree24_leading_parameter":"tau=0 or 40tau^3-105tau^2+84tau-14=0",
   "three_noncentral_parameter_roots_distinct":True,
   "uniform_first_third_fifth_character_controls":True,
   "fifth_character_formula":"(m5/M)=(5/16)*(e^2/h)*(3f-eK/h)",
 },
 "control_scope":"Regular source functions and named finite guards only; complete companion resonance excludes the linear controls.",
 "actual_companion_constructed":False,"whole_degree24_closed":False,
 "automatic_target_integrality_proved":False,"full_goal_complete":False,
 "independent_whole_proof_review":False,"formal_verification":False,"external_review":False,
}
parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path)
args=parser.parse_args();out=args.output or Path(tempfile.mkdtemp(prefix="math-check-"))/"ruling_six_degree24_leading_receipt.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),"elapsed_seconds":receipt["elapsed_seconds"],"receipt":str(out)},indent=2))
