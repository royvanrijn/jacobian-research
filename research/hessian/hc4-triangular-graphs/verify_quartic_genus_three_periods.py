"""Exact intrinsic-quartic genus-three chord and period controls."""
from pathlib import Path
import hashlib,json,tempfile,time
import sympy as sp

HERE=Path(__file__).resolve().parent
started=time.perf_counter()
checks=[]
def check(name,test):
    if not test:
        raise AssertionError(name)
    checks.append(name)
def zero(name,expr):
    check(name,sp.cancel(expr)==0)

t,s,A,B,C,Ap,Bp,Cp,L,Lp,J,U,V=sp.symbols(
    "t s A B C Ap Bp Cp L Lp J U V")
p=t**3+A*t
q=t**4+B*t*t+C*t
g=s**3-(B-2*A)*s-C
K=3*s*s+4*A
H=15*s*s+6*A-7*B
uu=2*A-3*B
vv=2*B-A
R=4*A*uu*uu+27*C*C
zero("normalization immersion resultant",
     sp.resultant(sp.diff(p,t),sp.diff(q,t),t)-R)
zero("collision-gap immersion resultant",sp.resultant(g,K,s)-R)
F=sp.resultant(p-U,q-V,t)
kap=-sp.cancel(sp.diff(F,V).subs({U:p,V:q},simultaneous=True)/sp.diff(p,t))
expected=t**6+(A+B)*t**4+2*C*t**3+(A*A-A*B+B*B)*t*t+C*(A+B)*t+A*(A-B)**2+C*C
zero("universal monic quartic-class conductor",kap-expected)
zero("ordinary triple-point conductor power",kap.subs(B,A)-(t**3+A*t+C)**2)
zero("order-three contact conductor power",kap.subs({B:2*A,C:0})-(t*t+A)**3)

h={0:sp.Integer(0),1:sp.Integer(1),2:s}
for i in range(3,8):
    h[i]=sp.expand(s*h[i-1]-(s*s+A)*h[i-2])
primitive=sp.integrate(sp.expand(p*sp.diff(q,t)),t)
raw=sum(coef*h[mon[0]] for mon,coef in sp.Poly(primitive,t).terms())
zero("compact chord-action formula modulo collision cubic",
     sp.rem(raw-K*K*H/420,g,s))
zero("repeated collision root has a nonzero action factor",
     H+2*K-7*sp.diff(g,s))

sy=((Bp-2*Ap)*s+Cp)/sp.diff(g,s)
total=lambda f:sp.diff(f,s)*sy+sp.diff(f,A)*Ap+sp.diff(f,B)*Bp+sp.diff(f,C)*Cp
gaplog=(-3*s*sy-2*Ap)/(-3*s*s-4*A)
inner=14*(6*A-5*B-3*s*s)*Ap+28*(3*s*s-A)*Bp+105*s*Cp
zero("full action derivative cancellation",
     total(K*K*H)+gaplog*K*K*H-K*inner)
invnum=3*uu*s*s+9*C*s+uu*uu
zero("explicit inverse of the gap square modulo cubic",
     sp.rem(K*invnum-R,g,s))
eq=sp.rem(Lp*K*H+L*inner+J*invnum/R,g,s)
E2=12*vv*Lp+42*(2*Bp-Ap)*L+3*J*uu/R
E1=45*C*Lp+105*Cp*L+9*J*C/R
E0=4*A*(6*A-7*B)*Lp+(14*(6*A-5*B)*Ap-28*A*Bp)*L+J*uu*uu/R
zero("exact quadratic coefficient equation",sp.Poly(eq,s).nth(2)-E2)
zero("exact linear coefficient equation",sp.Poly(eq,s).nth(1)-E1)
zero("exact constant coefficient equation",sp.Poly(eq,s).nth(0)-E0)
W=4*A-3*B
Wp=4*Ap-3*Bp
zero("retained triple-point factor and first integral",
     E0-uu*E2/3-4*(A-B)*(2*W*Lp+7*Wp*L))

# Divided differences reconstruct the collision equations and tangent determinant.
alpha,beta=sp.symbols("alpha beta")
ss=alpha+beta
vpair=alpha*beta
zero("cubic equality gives the stated pair product",
     (p.subs(t,beta)-p.subs(t,alpha))/(beta-alpha)
       -(ss*ss-vpair+A))
zero("quartic equality reduces to the collision cubic",
     sp.rem(sp.cancel((q.subs(t,beta)-q.subs(t,alpha))/(beta-alpha)
            +g.subs(s,ss)),beta**2+alpha*beta+alpha**2+A,beta))
tangentdet=3*s*(-3*s**3-(B+2*A)*s)-2*(B-2*A)*(-3*s*s-2*A)
zero("exact tangent determinant factor",
     tangentdet+K*sp.diff(g,s))
zero("canonical derivatives at distinct endpoints cannot be equal for nonzero sum",
     (3*beta*beta+A)-(3*alpha*alpha+A)-3*(alpha+beta)*(beta-alpha))

# Exact examples, including ordinary nodes, ordinary triple point and A5 contact.
examples=[]
for aa,bb,cc in [(1,0,1),(1,1,1),(1,2,0)]:
    pg=p.subs({A:aa,B:bb,C:cc})
    qg=q.subs({A:aa,B:bb,C:cc})
    rg=R.subs({A:aa,B:bb,C:cc})
    check("representative quartic normalization remains immersed",rg!=0)
    examples.append({"A":aa,"B":bb,"C":cc,
                     "immersion_resultant":str(rg),
                     "collision_cubic":str(g.subs({A:aa,B:bb,C:cc})),
                     "conductor":str(sp.factor(kap.subs({A:aa,B:bb,C:cc})))})
check("ordinary triple-point cubic is separable",
      sp.discriminant(g.subs({A:1,B:1,C:1}),s)!=0)
check("A5 case has exactly a triple zero collision sum",
      g.subs({A:1,B:2,C:0})==s**3)

# Persistent low-genus coefficient constancy via the coprime derivative degree bounds.
m,n=sp.symbols("m n",integer=True,positive=True)
check("triple-point polynomial derivative degree inequalities are incompatible",
      sp.simplify((m-1)-((2*m)+1)).is_negative is True)
zero("the discriminant derivative has the stated factors",
     sp.diff(4*A**3+27*C*C,A)*Ap+sp.diff(4*A**3+27*C*C,C)*Cp
       -12*A*A*Ap-54*C*Cp)

# The cable physical action lemma used in the fixed-order closure is weighted.
pp=t**3+A*t+3*C/4
qq=t**4+sp.Rational(2,3)*A*t*t+C*t
ant=sp.integrate(sp.expand(pp**2*sp.diff(qq,t)),t)
hh={0:sp.Integer(0),1:sp.Integer(1),2:s}
for i in range(3,12):
    hh[i]=sp.expand(s*hh[i-1]-(s*s+A)*hh[i-2])
weighted=sum(coef*hh[mon[0]] for mon,coef in sp.Poly(ant,t).terms())
gcable=s**3+sp.Rational(4,3)*A*s-C
zero("retained cable weighted action formula",
     sp.rem(weighted-C*C*(8*A*s+9*C)/40,gcable,s))
check("nonzero cable C prevents a triple-zero collision cubic",
      sp.Poly(gcable,s).nth(0)==-C)

audit_path=HERE/"sources/neumann-genus-three-source-audit.json"
audit=json.loads(audit_path.read_text())
source_hashes={}
for source in audit["sources"]:
    path=HERE/"sources"/source["local_pdf"]
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    check("pinned primary source "+source["local_pdf"],actual==source["sha256"])
    source_hashes[source["local_pdf"]]=actual
cable_proof=HERE/"CABLE_GENUS_THREE_WEIGHTED_PERIODS_AND_NONVANISHING.md"
cable_pin="9cadd3f5aed9c34845c0b143ebd74a50a8d99a658bdecf7e562b1f508eb96355"
check("retained cable physical-action proof pin",
      hashlib.sha256(cable_proof.read_bytes()).hexdigest()==cable_pin)
proof=HERE/"QUARTIC_GENUS_THREE_CHORD_ACTIONS_AND_PERIOD_EQUATIONS.md"
record={
    "outcome":"PASS","check_count":len(checks),"checks":checks,"examples":examples,
    "source_sha256":source_hashes,"cable_action_proof_sha256":cable_pin,
    "source_audit_sha256":hashlib.sha256(audit_path.read_bytes()).hexdigest(),
    "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
    "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "sympy_version":sp.__version__,
    "elapsed_seconds":round(time.perf_counter()-started,3),
    "scope":{
        "quartic_chord_action_and_three_period_equations_checked":True,
        "canonical_coefficients_proved_rational_not_everywhere_regular":True,
        "persistent_nonnodal_quartic_branches_excluded_by_written_argument":True,
        "ordinary_triple_point_and_A5_fixed_order_closures_use_cable_action_lemma":True,
        "entire_generic_ordinary_three_node_quartic_branch_excluded_here":False,
        "all_conductor_degree_six_cases_excluded":False,
        "uniform_family_target_automorphism_assumed":False,
        "automatic_trace_integrality_or_JC2_proved":False,
        "independent_replay":False,"independent_mathematical_review":False,
        "formal_verification":False,"external_review":False
    }
}
outdir=Path(tempfile.mkdtemp(prefix="quartic-g3-period-check-"))
receipt=outdir/"quartic_genus_three_period_verification.json"
receipt.write_text(json.dumps(record,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","checks":len(checks),"receipt":str(receipt)}))
