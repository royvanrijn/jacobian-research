"""Exact cable weighted-period identities; classification remains imported."""
from pathlib import Path
import argparse,hashlib,json,platform,tempfile,time
import sympy as S
START=time.perf_counter();HERE=Path(__file__).resolve().parent
PROOF=HERE/"CABLE_GENUS_THREE_WEIGHTED_PERIODS_AND_NONVANISHING.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(value):return S.cancel(S.expand(value))==0
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
X,Y,H,a,b,c,D,level,beta,alpha=S.symbols("X Y H a b c D level beta alpha")
HH=X**3+Y**2+a*X/2+b/2
f1=(X**3+Y**2)**2+a*X**4+b*X**3+a**2*X**2/4+c*X+a*X*Y**2+b*Y**2+D*Y
check("primary cable square completion",zero(f1-(HH**2+(c-a*b/2)*X+D*Y-b*b/4)))
Psi=(beta-H**2-alpha*X)**2+D**2*X**3+a*D**2*X/2+b*D**2/2-D**2*H
Yexpr=(beta-H**2-alpha*X)/D
check("actual quotient-ring elimination",zero(Psi-D**2*(X**3+Yexpr**2+a*X/2+b/2-H)))
check("intrinsic exact H degree four",S.degree(Psi,H)==4)
check("intrinsic exact X degree three",S.degree(Psi,X)==3)
for term in [H*X**2,H**3,H*X]:
    check("forbidden intrinsic coefficient "+str(term),S.Poly(Psi,H,X).coeff_monomial(term)==0)
check("intrinsic nonzero linear H coefficient",S.Poly(Psi,H,X).coeff_monomial(H)==-D**2)

# Independent resultant and all physical normalization coefficients.
t,U,V,A,B,C,L,eta,E=S.symbols("t U V A B C L eta E")
p=t**3+A*t;q=t**4+B*t**2+C*t
R=S.resultant(p-U,q-V,t)
expected=U**4+3*C*U**3+(3*B-4*A)*U**2*V+(B*(A-B)**2+3*C**2)*U**2+C*(3*B-5*A)*U*V+C*(A*(A-B)**2+C**2)*U-V**3+2*A*(A-B)*V**2-A*(A*(A-B)**2+C**2)*V
check("complete intrinsic quartic resultant",zero(R-expected))
physical=S.Poly(S.expand(R.subs({U:U-eta,V:V-L*U-E},simultaneous=True)),U,V)
check("physical mixed UV2 coefficient",zero(physical.coeff_monomial(U*V**2)-3*L))
physical0=S.Poly(physical.as_expr().subs(L,0),U,V)
check("physical cubic U coefficient",zero(physical0.coeff_monomial(U**3)-(3*C-4*eta)))
centered=S.Poly(S.expand(physical0.as_expr().subs(eta,3*C/4)),U,V)
check("physical mixed UV coefficient",zero(centered.coeff_monomial(U*V)-C*(2*A-3*B)/2))
check("C zero kills every odd U term",zero(centered.as_expr().subs(C,0)-centered.as_expr().subs({C:0,U:-U})))
cable=S.Poly(S.expand(centered.as_expr().subs(B,S.Rational(2,3)*A)),U,V)
check("restricted normalized linear U coefficient",zero(cable.coeff_monomial(U)+C**3/8))

# Exact polynomial one-form modulo the defining fibre.
FF=H**2+alpha*X+D*Y-beta
primitive=-alpha*X**2/2-D*X*Y+beta*X
check("weighted form dX coefficient",zero(H**2-FF-S.diff(primitive,X)))
check("weighted form dY coefficient",zero(-D*X-S.diff(primitive,Y)))

# Universal matching-pair weighted period, using symmetric recurrences.
sigma,pi,gap=S.symbols("sigma pi gap")
pp=t**3+A*t+3*C/4
qq=t**4+S.Rational(2,3)*A*t*t+C*t
integ=S.Rational(2,5)*t**10+S.Rational(7,6)*A*t**8+C*t**7+S.Rational(10,9)*A*A*t**6+2*A*C*t**5+(A**3/3+S.Rational(15,16)*C*C)*t**4+A*A*C*t**3+S.Rational(9,8)*A*C*C*t*t+S.Rational(9,16)*C**3*t
check("entire weighted primitive",zero(S.diff(integ,t)-pp**2*S.diff(qq,t)))
seq=[S.Integer(0),S.Integer(1)]
for n in range(2,11):seq.append(S.expand(sigma*seq[-1]-(sigma**2+A)*seq[-2]))
quotient=S.expand(sum(coef*seq[exponent[0]] for exponent,coef in S.Poly(integ,t).terms()))
g=sigma**3+S.Rational(4,3)*A*sigma-C
reduced=S.rem(quotient,g,sigma)
formula=C**2*(8*A*sigma+9*C)/40
check("universal weighted action quotient",zero(reduced-formula))
check("matching-pair cubic",zero(g-(sigma*(sigma**2+S.Rational(4,3)*A)-C)))
check("zero gap would force C zero",zero(g-sigma*(3*sigma**2+4*A)/3+C))
check("A zero all actions nonzero",zero(formula.subs(A,0)-S.Rational(9,40)*C**3))
check("nonzero A quotient is genuinely linear",S.degree(formula,sigma)==1)
check("triple cubic root would force C zero",g.subs({sigma:0,A:0})==-C)
check("exact weighted-action resultant",zero(S.resultant(g,formula,sigma)-C**7*(1280*A**3+729*C**2)/64000))

# Actual rational cable curves, with independent original source identities.
zz=t*t+A/3
hc=t**3+A*t-S.Rational(3,2)
xc=zz**2-2*t
yc=S.I*(zz**3-3*t*zz+S.Rational(3,2))
check("rational cable normalization H identity",zero(xc**3+yc**2+A*xc+S.Rational(3,4)-hc))
check("rational cable normalization fibre identity",zero(hc**2-A*xc+S.I*yc-(S.Rational(3,4)-4*A**3/27)))
check("rational cable primary f1 value",zero(f1.subs({X:xc,Y:yc,a:2*A,b:S.Rational(3,2),c:A/2,D:S.I},simultaneous=True)-(S.Rational(3,16)-4*A**3/27)))
check("normalization degrees coprime",S.degree(hc,t)==3 and S.degree(xc,t)==4)
check("normalization immersion at H critical points",zero(S.rem(S.diff(xc,t)+2,3*t*t+A,t)))
# General monogenic conductor ratio before the harmless translations.
Rcab=R.subs(B,S.Rational(2,3)*A)
pc=t**3+A*t;qc=t**4+S.Rational(2,3)*A*t*t+C*t
kappa=S.cancel(-S.diff(Rcab,V).subs({U:pc,V:qc},simultaneous=True)/S.diff(pc,t))
kapformula=t**6+S.Rational(5,3)*A*t**4+2*C*t**3+S.Rational(7,9)*A*A*t*t+S.Rational(5,3)*A*C*t+A**3/9+C*C
check("exact intrinsic conductor ratio",zero(kappa-kapformula))
check("degree-six conductor",S.degree(kappa,t)==6)
check("triple-point conductor control",zero(kappa.subs({A:0,C:-2})-(t**3-2)**2))
check("Cminus2 weighted quotient",zero(formula.subs(C,-2)-(4*A*sigma-9)/5))
# A*sigma=9/4 and sigma^3=-5 give a genuine zero individual action.
subA=S.Rational(9,4)/sigma
check("zero-one-action matching cubic",zero(g.subs({A:subA,C:-2})-(sigma**3+5)))
check("zero-one-action weighted quotient",zero(formula.subs({A:subA,C:-2})))
check("zero-one-action nonzero gap certificate",zero((-3*sigma**2-4*subA)*sigma+3*(sigma**3+5)-6))
check("zero-one-action parameter cube",S.rem(S.cancel((subA**3+S.Rational(729,320))*sigma**3),sigma**3+5,sigma)==0)
disc=S.discriminant(g,sigma)
check("zero-one-action other pairs remain distinct",zero(disc.subs(C,-2).subs(A**3,-S.Rational(729,320))+S.Rational(432,5)))

audit_path=HERE/"sources/neumann-genus-three-source-audit.json"
audit=json.loads(audit_path.read_text());source_hashes={}
for source in audit["sources"]:
    path=HERE/"sources"/source["local_pdf"]
    actual=sha(path)
    check("pinned primary source "+source["local_pdf"],actual==source["sha256"])
    source_hashes[source["local_pdf"]]=actual
check("classification remains attributed written input",audit["source_classification_proved_by_symbolic_checker"] is False)
receipt={
 "outcome":"PASS","check_count":len(checks),"checks":checks,
 "elapsed_seconds":time.perf_counter()-START,
 "environment":{"python":platform.python_version(),"sympy":S.__version__},
 "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
 "source_audit_sha256":sha(audit_path),"source_sha256":source_hashes,
 "retained_proof_sources":{"GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md":sha(HERE/"GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md")},
 "scope":{
  "primary_cable_ambient_class_only":True,
  "normalization_A1_required":True,
  "cable_weighted_action_identity":True,
  "entire_matching_pair_action_vector_nonzero":True,
  "every_individual_action_nonzero":False,
  "intrinsic_pair_is_target_automorphism":False,
  "generic_cable_family_closed":False,
  "all_degree_six_conductors_closed":False,
  "actual_Keller_map_constructed":False,
  "automatic_target_integrality_or_JC2_proved":False,
  "independent_replay":False,"independent_mathematical_review":False,
  "formal_verification":False,"external_review":False,
 },
}
parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path)
args=parser.parse_args();out=args.output or Path(tempfile.mkdtemp(prefix="cable-weighted-period-check-"))/"cable_genus_three_weighted_periods_verification.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),"elapsed_seconds":receipt["elapsed_seconds"],"receipt":str(out)},indent=2))
