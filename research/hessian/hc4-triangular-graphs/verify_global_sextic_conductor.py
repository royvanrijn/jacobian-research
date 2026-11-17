"""Exact controls for the global sextic conductor pole-repair proof."""
from pathlib import Path
import argparse,hashlib,json,platform,tempfile,time
import sympy as s
START=time.perf_counter()
BASE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    if not bool(value):raise AssertionError(name)
    checks.append(name)
def zero(name,value):check(name,s.cancel(value)==0)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def initial(expr,var):
    num,den=s.fraction(s.cancel(expr))
    if num==0:return s.oo,s.Integer(0)
    pn,pd=s.Poly(num,var),s.Poly(den,var)
    n=min(m[0] for m in pn.monoms());d=min(m[0] for m in pd.monoms())
    return n-d,s.cancel(pn.nth(n)/pd.nth(d))

u,y,z,A,B,C,U0,S0,T0=s.symbols("u y z A B C U0 S0 T0")
p=u**3+A*u
q=u**4+B*u*u+C*u
g=z**3-(B-2*A)*z-C
K=3*z*z+4*A
R=4*A*(2*A-3*B)**2+27*C*C
kap=u**6+(A+B)*u**4+2*C*u**3+(A*A-A*B+B*B)*u*u+C*(A+B)*u+A*(A-B)**2+C*C
zero("exact monic conductor from all collision sums",
     s.resultant(g,u*u-z*u+z*z+A,z)-kap)
zero("exact source immersion resultant",
     s.resultant(s.diff(p,u),s.diff(q,u),u)-R)
zero("same immersion factor from collision gaps",s.resultant(g,K,z)-R)
zero("integral coefficient descent product",
     ((A+B)**2-(A*A-A*B+B*B))/3-A*B)
zero("A satisfies the monic coefficient polynomial",
     A*A-(A+B)*A+((A+B)**2-(A*A-A*B+B*B))/3)
zero("B satisfies the same monic coefficient polynomial",
     B*B-(A+B)*B+((A+B)**2-(A*A-A*B+B*B))/3)

# An exact arbitrary-coefficient inclusion matrix over C[p].
W=s.symbols("W")
columns=[s.Poly(s.rem(q**j,p-W,u),u) for j in range(3)]
matrix=s.Matrix(3,3,lambda i,j:columns[j].nth(i))
index=s.Poly(s.expand(matrix.det()),W)
check("arbitrary-coefficient monogenic index has degree three",index.degree()==3)
check("index leading coefficient is a constant unit",index.LC()==1)
check("three leading residue classes are distinct",
      [4*j%3 for j in range(3)]==[0,1,2])
check("the exact small missing source degrees are one two five",
      [n for n in range(12) if not any(n==3*i+4*j for i in range(5) for j in range(3))]==[1,2,5])

# Cable physical period differentiation with BOTH A and C varying.
Ap,Cp,L,Lp,J=s.symbols("Ap Cp L Lp J")
gc=z**3+4*A*z/3-C
zprime=(3*Cp-4*Ap*z)/(9*z*z+4*A)
gap2=-3*z*z-4*A
gap2prime=-6*z*zprime-4*Ap
RC=C*C*(8*A*z+9*C)
RCprime=2*C*Cp*(8*A*z+9*C)+C*C*(8*Ap*z+8*A*zprime+9*Cp)
expression=s.cancel(Lp*RC+L*(RCprime+gap2prime*RC/(2*gap2))+J)
numer,denom=s.fraction(expression)
rem=s.Poly(s.rem(numer,gc,z),z)
E1=9*C**3*Lp+30*C*C*Cp*L+J
E2=4*A*C*Lp+10*A*Cp*L+5*Ap*C*L
zero("cable constant-coefficient remainder identity",
     rem.nth(0)-8*A*(2*A*E1+3*C*C*E2))
zero("cable linear-coefficient remainder gives exact physical ODE",
     rem.nth(1)-27*C*E1)
zero("cable quadratic-coefficient remainder identity",
     rem.nth(2)-6*(2*A*E1+9*C*C*E2))
zero("cable finite-pair coefficient relation",
     s.rem(C-K*z/3,gc,z))
zero("cable weighted factor at a finite separated pair",
     (RC-K*K*z**3*(5*K-6*z*z)/9).subs(C,K*z/3))
zero("constant C physical equation has nonzero affine scale slope",
     E1.subs(Cp,0)-(9*C**3*Lp+J))

# The quartic pole identity keeps every possible W case.
Wh=s.symbols("Wh")
H=15*z*z+6*A-7*B
zero("quartic H in terms of finite gap and first integral",
     H-(s.Rational(35,2)*z*z-s.Rational(5,6)*K+s.Rational(7,3)*(4*A-3*B)))
zero("three distinct fine-limit collision sums",
     z**3-2*z+4-(z+2)*(z*z-2*z+2))
check("the other fine-limit gaps are nonzero",
      s.expand((-12+6*s.I)*(-12-6*s.I))==180)
n=s.symbols("n",integer=True,positive=True)
zero("quartic fine pole coefficient never cancels",
     45*(2*n)+105*(-3*n)+225*n)
zero("cable pole coefficient never cancels",
     9*(5*n)+30*(-n)-15*n)
check("positive integer cable pole cannot have half order",
      s.solve(s.Eq(2*n-1,0),n)==[])
check("quartic compared orders are strictly separated",
      s.simplify(4*n-(2*n-1)).is_positive is True)

v=s.symbols("v")
for order in [1,2,3,7]:
    hh=v**(-order);sigma=-2*hh
    AA=(1-3*sigma*sigma)/4
    for mode in ["zero","fine"]:
        if mode=="fine" and order!=7:continue
        WW=0 if mode=="zero" else v**-4
        BB=(4*AA-WW)/3
        CC=sigma**3-(BB-2*AA)*sigma
        RR=4*AA*(2*AA-3*BB)**2+27*CC*CC
        lam=v**(2*order)
        ro,rl=initial(RR,v)
        check(f"quartic exact resultant pole order n{order} W{mode}",ro==-4*order)
        zero(f"quartic exact resultant leading coefficient n{order} W{mode}",rl-180)
        lhs=45*s.diff(lam,v)+105*s.diff(CC,v)*lam/CC
        oo,ll=initial(lhs,v)
        check(f"quartic period left exact order n{order} W{mode}",oo==2*order-1)
        zero(f"quartic period forbidden coefficient n{order} W{mode}",ll+225*order)
        check(f"quartic left and right cannot balance n{order} W{mode}",
              initial(1/RR,v)[0]>oo)
    # A cable finite pair can have no positive integral centroid pole.
    AA=(1-3*sigma*sigma)/4
    CC=sigma/3
    lam=v**(5*order)
    weighted=CC*CC*(8*AA*sigma+9*CC)
    check(f"cable canonical action has exact pole five n n{order}",
          initial(weighted,v)[0]==-5*order)
    lhs=9*CC**3*s.diff(lam,v)+30*CC*CC*s.diff(CC,v)*lam
    oo,ll=initial(lhs,v)
    check(f"cable left exact positive order n{order}",oo==2*order-1>0)
    zero(f"cable left forbidden coefficient n{order}",ll+s.Rational(40,9)*order)

for order in [1,2,4]:
    AA=-s.Rational(1,2);BB=v**(-order)
    WW=4*AA-3*BB;HH=15+6*AA-7*BB;lam=v**order
    check(f"finite unit action forces scale order k at B pole k{order}",
          initial(lam*HH,v)[0]==0)
    check(f"quartic first integral contradicts that B pole k{order}",
          initial(lam*lam*WW**7,v)[0]==-5*order!=0)

# A genuine non-Keller immersion family passes the entire constant-degree
# conductor architecture and has a UNIT immersion resultant.
AA=3*(y*y-1);BB=2*(y*y-1)-s.I;CC=2*y
P=p.subs(A,AA);Q=q.subs({B:BB,C:CC})
zero("non-Keller geometric control has unit resultant",
     R.subs({A:AA,B:BB,C:CC})-108)
Jac=s.expand(s.diff(P,u)*s.diff(Q,y)-s.diff(P,y)*s.diff(Q,u))
check("geometric constant-degree control is not a Keller pair",
      Jac.has(u,y) and Jac!=0)
zero("geometric control exhibits vanishing source Jacobian on u zero",
     Jac.subs(u,0))
check("geometric control retains a degree-six monic conductor",
      s.Poly(kap.subs({A:AA,B:BB,C:CC}),u).degree()==6)

# Source inputs are pinned; finite calculations do not establish their proofs.
audit_path=BASE/"sources/neumann-genus-three-source-audit.json"
audit=json.loads(audit_path.read_text())
pdfs={}
for source in audit["sources"]:
    path=BASE/"sources"/source["local_pdf"]
    actual=sha(path)
    check("primary PDF identity "+source["local_pdf"],actual==source["sha256"])
    pdfs[source["local_pdf"]]=actual
proof=BASE/"GLOBAL_SEXTIC_CONDUCTOR_EXCLUSION_AND_PERIOD_POLE_REPAIR.md"
parents=[
"GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md",
"UNRESTRICTED_NODAL_ACTION_AND_DEGREE_TWO_CONDUCTOR_EXCLUSION.md",
"CONSTANT_DEGREE_FOUR_CONDUCTORS_AND_PERIOD_ODE.md",
"GLOBAL_DEGREE_FOUR_CONDUCTOR_EXCLUSION.md",
"HYPERELLIPTIC_SEXTIC_CONDUCTORS_AND_THREE_PERIOD_RIGIDITY.md",
"QUARTIC_GENUS_THREE_CHORD_ACTIONS_AND_PERIOD_EQUATIONS.md",
"CABLE_GENUS_THREE_WEIGHTED_PERIODS_AND_NONVANISHING.md"]
record={
"outcome":"PASS","check_count":len(checks),"checks":checks,
"elapsed_seconds":round(time.perf_counter()-START,3),
"environment":{"python":platform.python_version(),"sympy":s.__version__},
"proof_sha256":sha(proof),"checker_sha256":sha(Path(__file__)),
"source_audit_sha256":sha(audit_path),"source_sha256":pdfs,
"retained_proof_sources":{p:sha(BASE/p) for p in parents},
"scope":{
"all_global_directional_degree_six_conductors_excluded_by_written_proof":True,
"special_degree_six_to_four_or_two_drops_retained":True,
"generic_family_target_automorphism_assumed":False,
"pointwise_classification_coverage_proved_by_filtered_rank_and_closed_period_conditions":True,
"noninvertible_global_directional_conductor_degree_at_least_eight":True,
"every_special_fibre_conductor_degree_at_least_eight":False,
"primary_classification_proved_by_symbolic_checker":False,
"automatic_trace_integrality_or_JC2_proved":False,
"independent_replay":False,"formal_verification":False,"external_review":False}}
parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path)
args=parser.parse_args()
out=args.output or Path(tempfile.mkdtemp(prefix="global-sextic-conductor-check-"))/"global_sextic_conductor_verification.json"
out.write_text(json.dumps(record,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),"receipt":str(out)}))
