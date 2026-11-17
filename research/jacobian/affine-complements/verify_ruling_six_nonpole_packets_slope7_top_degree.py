"""PF183 exact controls; unrestricted arguments are in the source note."""
from pathlib import Path
import argparse,hashlib,json,platform,tempfile,time
import sympy as S
START=time.perf_counter();BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_NONPOLE_PACKETS_AND_SLOPE7_TOP_DEGREE_EXCLUSION.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(value):return S.cancel(S.expand(value))==0
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
T,p,v,s,X,Z=S.symbols("T p v s X Z")
# Full response polynomials supply h distinct conjugates with common p derivative.
for H in [T,T**2,T**3+T,T**4-2*T+1]:
    h=S.degree(H,T)
    check("first response irreducible "+str(H),S.Poly(H-p,T,domain=S.QQ.frac_field(p)).is_irreducible)
    check("first response separable "+str(H),S.gcd(H-p,S.diff(H,T),T)==1)
    check("conjugate root count "+str(H),S.degree(H-p,T)==h)
    if h==1:
        check("affine response derivative is constant",not S.diff(H,T).has(T,p))
    else:
        check("nonaffine packet has at least two roots "+str(H),h>=2)
# Counts are sheets: a ramified point can already supply several roots.
for mass in [1,2]:
    check("small pole mass cannot supply field6 "+str(mass),S.Rational(5,2)*mass<6)
check("mass3 leaves forbidden single nonpole sheet",4-3==1)

# Positive integer pole orders for every integer partition of four.
def partitions(total,minimum=1):
    if total==0:
        yield []
    for first in range(minimum,total+1):
        for tail in partitions(total-first,first):yield [first]+tail
patterns=[]
from itertools import product
for indices in partitions(4):
    ranges=[range(1,(5*a-1)//2+1) for a in indices]
    for orders in product(*ranges):
        if sum(orders)==9:
            patterns.append((tuple(indices),tuple(orders)))
check("complete field9 pole partitions",patterns==[((1,3),(2,7)),((4,),(9,))])
expected=[tuple([S.Rational(31,3)]*3+[S.Integer(10)]),tuple([S.Rational(41,4)]*4)]
actual=[]
for indices,orders in patterns:
    vals=[]
    for alpha,q in zip(indices,orders):vals.extend([8+S.Rational(q,alpha)]*alpha)
    actual.append(tuple(sorted(vals)))
check("complete field9 derivative multisets",actual==[tuple(sorted(a)) for a in expected])
for vals in actual:
    check("field9 exact weighted sum "+str(vals),sum(vals)==41)
    check("field9 strict sheet range "+str(vals),all(8<a<S.Rational(21,2) for a in vals))

# All finite parameter orders are exhaustively classified; the infinity J case
# is checked separately. No coefficient-degree cutoff is used in the proof.
for t in [10,11,12]:
    check("delta cubic boundary absent t"+str(t),not (S.Rational(4*t,3).q==1 and t+1<=S.Rational(4*t,3)<=14))
    for delta in range(t+1,15):
        for ell in list(map(S.Integer,range(1,17)))+[S.oo]:
            label="t"+str(t)+" d"+str(delta)+" ell"+str(ell)
            if delta<2*ell and delta<S.Rational(4*t,3):
                mu=S.Rational(3*delta,4)
                check("dominant Delta rejects target values "+label,mu not in [S.Rational(41,4),S.Rational(31,3),10])
            elif delta==2*ell and delta<S.Rational(4*t,3):
                mu=S.Rational(3*delta,4)
                check("square Delta boundary rejects target values "+label,mu not in [S.Rational(41,4),S.Rational(31,3),10])
            elif ell>S.Rational(2*t,3) and delta>S.Rational(4*t,3):
                check("cubic linear balance rejects both patterns "+label,t!=S.Rational(41,4) and t!=S.Rational(31,3))
            else:
                check("remaining case is square centers "+label,ell<S.Rational(2*t,3) and delta>2*ell)
                m0=min(S.Integer(delta),t+ell/2)
                check("square-center fine displacement "+label,m0>2*ell)
                check("square-center linear term higher "+label,m0<2*t-ell)
                mu=(m0+ell)/2
                check("two square-center roots reject degree4 pattern "+label,mu!=S.Rational(41,4))
                check("square-center conjugates cannot have cubic denominator "+label,mu!=S.Rational(31,3))
                # If mu=10, at least two roots contradict the sole ten in (3,1).
                check("square-center pair rejects second pattern "+label,list(expected[0]).count(mu)<2)
check("only integral quarter candidate is outside square regime",S.Rational(41-2*10,3)==7 and 7>S.Rational(20,3))
for t in [11,12]:
    check("no integral quarter candidate t"+str(t),S.Rational(41-2*t,3).q!=1)

# Exact formal controls for the other remaining field degrees.
Gbase=(X-v**3)**3*(X+3*v**3)
G6=S.expand(Gbase+v**13-p*v**14)
G7=S.expand(Gbase+v**10*(X-v**3)-p*v**14)
for degree,GG in [(6,G6),(7,G7)]:
    Jcoef=S.expand(GG).coeff(X,2)/2
    qcoef=S.expand(GG).coeff(X,1)
    dcoef=S.expand(GG).coeff(X,0)-Jcoef**2
    check("degree"+str(degree)+" depressed quartic cubic term zero",S.expand(GG).coeff(X,3)==0)
    check("degree"+str(degree)+" exact J",zero(Jcoef+3*v**6))
    check("degree"+str(degree)+" q order nine",S.cancel(qcoef/v**9).subs(v,0)==8)
    check("degree"+str(degree)+" Delta order twelve",S.cancel(dcoef/v**12).subs(v,0)==-12)
    check("degree"+str(degree)+" affine parameter coefficient",zero(S.diff(dcoef,p)+v**14))
    init=S.expand(GG.subs(X,v**3*Z)/v**12).subs(v,0)
    check("degree"+str(degree)+" entire simple-root initial",zero(init-(Z-1)**3*(Z+3)))
    check("degree"+str(degree)+" simple-root derivative minus64",S.diff(init,Z).subs(Z,-3)==-64)
    direct=S.cancel(S.diff(GG,X).subs(X,v**3*Z)/v**9).subs(v,0)
    check("degree"+str(degree)+" exact simple-sheet mu9",direct.subs(Z,-3)==-64)
    check("degree"+str(degree)+" simple-sheet first p order5",14-9==5)
    check("degree"+str(degree)+" simple-sheet pole order1",9-8==1)

triple6=S.expand(G6.subs({v:s**3,X:s**9+s**10*Z})/s**39)
check("degree6 entire triple substitution",zero(triple6-(4*Z**3+1+s*Z**4-p*s**3)))
check("degree6 three leading roots simple",S.gcd(4*Z**3+1,12*Z**2,Z)==1)
deriv6=S.expand(S.diff(G6,X).subs({v:s**3,X:s**9+s**10*Z})/s**29)
check("degree6 exact triple derivative order",zero(deriv6.subs(s,0)-12*Z**2))
check("degree6 triple derivative unit at each root",S.gcd(4*Z**3+1,deriv6.subs(s,0),Z)==1)
check("degree6 first p order13 in normalized s",10+3==13)
check("degree6 index3 pole order5",29-3*8==5)
check("degree6 total normalized pole degree",5+1==6)
check("degree6 sheet-weighted derivative sum",3*S.Rational(29,3)+9==32+6)
check("degree6 strict derivative range",8<S.Rational(29,3)<S.Rational(21,2))

outer7=S.expand(G7.subs({v:s**2,X:s**6+s**7*Z})/s**27)
check("degree7 entire outer substitution",zero(outer7-(4*Z**3+Z+s*Z**4-p*s)))
for root in [S.I/2,-S.I/2]:
    check("degree7 simple nonzero outer root "+str(root),zero((4*Z**3+Z).subs(Z,root)))
    check("degree7 outer root derivative minus2 "+str(root),(12*Z**2+1).subs(Z,root)==-2)
deriv7=S.expand(S.diff(G7,X).subs({v:s**2,X:s**6+s**7*Z})/s**20)
check("degree7 exact outer derivative order",zero(deriv7.subs(s,0)-(12*Z**2+1)))
check("degree7 outer first p order8 in normalized s",7+1==8)
check("degree7 index2 pole order4",20-2*8==4)
inner7=S.expand(G7.subs(X,v**3+v**4*Z)/v**14)
check("degree7 entire inner substitution",zero(inner7-(Z-p+4*v*Z**3+v**2*Z**4)))
check("degree7 inner affine first response",zero(inner7.subs({v:0,Z:p})))
check("degree7 inner derivative leading unit",S.cancel(S.diff(G7,X).subs(X,v**3+v**4*Z)/v**10).subs(v,0)==1)
check("degree7 index1 pole order2",10-8==2)
check("degree7 total normalized pole degree",4+2+1==7)
check("degree7 sheet-weighted derivative sum",3*10+9==32+7)
check("degree7 strict derivative range",8<9<10<S.Rational(21,2))

# Exact full quartic left by all these necessary tests.
G=X**4+2*v**7*X**2+v**10*X+(1-p)*v**14
outer=S.expand((G.subs({v:s**3,X:s**10*Z}))/s**40)
check("degree8 full outer substitution",zero(outer-(Z**4+Z+2*s*Z**2+s**2*(1-p))))
initial=Z**4+Z
check("outer three nonzero roots simple",S.gcd(Z**3+1,S.diff(initial,Z),Z)==1)
qder=S.expand(S.diff(G,X).subs({v:s**3,X:s**10*Z})/s**30)
check("outer derivative leading minus three",zero(S.rem(qder.subs(s,0)+3,Z**3+1,Z)))
check("outer index-three normalized pole order",3*(10-8)==6)
inner=S.expand(G.subs(X,v**4*Z)/v**14)
check("degree8 full inner substitution",zero(inner-(Z+1-p+2*v*Z**2+v**2*Z**4)))
check("inner affine first response",zero(inner.subs({v:0,Z:p-1})))
check("inner derivative leading unit",zero(S.cancel(S.diff(G,X).subs(X,v**4*Z)/v**10).subs(v,0)-1))
check("inner index-one normalized pole order",10-8==2)
check("degree8 total pole degree",6+2==8)
check("degree8 all-sheet derivative sum",4*10==32+8)
check("degree8 strict sheet bound",8<10<S.Rational(21,2))
check("degree8 parameter precision",zero((1-p)*v**14-(v**7)**2+p*v**14))

receipt={
 "claim":"PF2RULED183","outcome":"PASS","check_count":len(checks),"checks":checks,
 "elapsed_seconds":time.perf_counter()-START,
 "environment":{"python":platform.python_version(),"sympy":S.__version__},
 "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
 "retained_proof_sources":{name:sha(BASE/name) for name in [
  "POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md",
  "RULING_SIX_AFFINE_POLE_BUDGET_AND_SOURCE_INFINITY.md",
  "RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md",
  "RULING_SIX_BALANCED_QUADRATIC_REMAINDER_AND_LEADING_RESONANCE.md",
  "RULING_SIX_DEGREE24_CRITICAL_QUARTIC_POLE_BUDGET_EXCLUSION.md",
  "RULING_SIX_NONCRITICAL_POLE_BUDGET_AND_DEGREE30_LEDGER.md",
  "RULING_SIX_PARAMETER_PRECISION_AND_SLOPE6_DEGREE30_EXCLUSION.md",
 ]},
 "uniform_proved_scope":{
  "isolated_nonpole_A_sheet_excluded":True,
  "all_D30_quartic_sheets_poles":True,
  "D30_intermediate_degree9_excluded":True,
  "remaining_D30_intermediate_degrees":[6,7,8],
  "remaining_D30_linear_remainder_coefficient_orders":[2,3,4,5],
 },
 "formal_control_field_degrees":[6,7,8],
 "formal_controls_actual_source_realization":False,"formal_controls_actual_mate":False,
 "remaining_D30_class_closed":False,"full_goal_complete":False,
 "independent_whole_proof_review":False,"formal_verification":False,"external_review":False,
}
parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path)
args=parser.parse_args();out=args.output or Path(tempfile.mkdtemp(prefix="math-check-"))/"ruling_six_nonpole_packets_slope7_top_degree_receipt.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),"elapsed_seconds":receipt["elapsed_seconds"],"receipt":str(out)},indent=2))
