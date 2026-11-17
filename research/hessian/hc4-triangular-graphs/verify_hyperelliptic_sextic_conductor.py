"""Exact three-period rigidity controls for the generic hyperelliptic sextic class."""
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

T,u,r,y,a,b,c,d,ap,bp,cp,dp,L,Lp,J=sp.symbols(
    "T u r y a b c d ap bp cp dp L Lp J")
g=a*T**3+b*T*T+c*T+d
C=5*a*T**4-9*c*T*T-30*d*T
p=u*u
q=u*g.subs(T,u*u)
period=sp.integrate(sp.expand(p*sp.diff(q,u)),(u,-r,r))
difference=sp.together(period/r-sp.Rational(8,315)*C.subs(T,r*r))
zero("canonical three-node period modulo the root equation",
     sp.rem(difference,g.subs(T,r*r),r))
Ty=-(ap*T**3+bp*T*T+cp*T+dp)/sp.diff(g,T)
Cy=sp.diff(C,T)*Ty+sp.diff(C,a)*ap+sp.diff(C,c)*cp+sp.diff(C,d)*dp
expected=-(35*ap*T**4+45*bp*T**3+63*cp*T*T+105*dp*T)/2
numer=sp.together(Cy+Ty*C/(2*T)-expected).as_numer_denom()[0]
zero("full variable-coefficient period derivative modulo cubic",
     sp.rem(numer,g,T))

M=2*Lp*C-L*(35*ap*T**4+45*bp*T**3+63*cp*T*T)+4*J
E1=70*d*d*a*Lp-35*d*d*ap*L+4*J*a*c
E2=10*a*b*Lp+(45*a*bp-35*b*ap)*L+4*J*a*a/d
E3=28*a*c*Lp+(63*a*cp-35*c*ap)*L+4*J*a*b/d
quot=(10*Lp-35*ap*L/a)*T+4*J/d
zero("all three differential equations from polynomial divisibility",
     M-g*quot+E2/a*T**3+E3/a*T*T+E1/(a*d)*T)
drop2=105*d*d*(3*a*cp-c*ap)*L-4*J*a*(2*c*c-5*b*d)
zero("degree-six to degree-two local elimination",
     5*d*d*E3-2*c*E1-drop2)

# Complete integer order comparisons; formal leading exponents remain symbolic.
m,n=sp.symbols("m n",positive=True,integer=True)
zero("drop-four exceptional derivative exponent",10*sp.Rational(7,2)*m-35*m)
check("drop-four resonance is above the matching exponent for all positive m",
      sp.simplify(sp.Rational(7,2)*m-(m+1)).is_positive is True)
zero("drop-two forbidden coefficient",28-35*m-7*(4-5*m))
mm=sp.symbols("mm")
check("drop-two coefficient only vanishes at noninteger four fifths",
      sp.solve(sp.Eq(28-35*mm,0),mm)==[sp.Rational(4,5)])

# Concrete exact local models guard both leading-coefficient degeneration types.
v=sp.symbols("v")
for order in [1,2,4,7]:
    aa=v**order+v**(order+1)
    bb=1+v
    cc=v+v*v
    ll=v**(order+1)
    second=sp.expand(10*aa*bb*sp.diff(ll,v)
                     +(45*aa*sp.diff(bb,v)-35*bb*sp.diff(aa,v))*ll)
    zero(f"drop-four forced derivative coefficient m={order}",
         sp.expand(second).coeff(v,2*order)-(10-25*order))
    third=sp.expand(28*aa*cc*sp.diff(ll,v)
                    +(63*aa*sp.diff(cc,v)-35*cc*sp.diff(aa,v))*ll)
    check(f"drop-four third equation has too high an order m={order}",
          min(k[0] for k,coef in sp.Poly(third,v).terms() if coef!=0)>order)
    bb=v+v*v
    cc=1+v
    denom=105*(3*aa*sp.diff(cc,v)-cc*sp.diff(aa,v))
    scale=sp.cancel(4*aa*(2*cc*cc-5*bb)/denom)
    check(f"drop-two scale has exact order one m={order}",
          sp.limit(scale/v,v,0)!=0)
    third=sp.cancel(28*aa*cc*sp.diff(scale,v)
             +(63*aa*sp.diff(cc,v)-35*cc*sp.diff(aa,v))*scale)
    lead=sp.limit(third/v**order,v,0)
    zero(f"drop-two forbidden derivative coefficient m={order}",
         lead-(28-35*order)*sp.limit(scale/v,v,0))
    check(f"drop-two leading coefficient is nonzero m={order}",lead!=0)

# Constant-a polynomial Lambda degree calculation.
ell,l1,l2,d0,J0=sp.symbols("ell l1 l2 d0 J0")
bexpr=35*d0**3*(28*l1*l1+63*ell*l2)/(8*J0**2)
leadb=sp.factor(bexpr.subs({ell:1,l1:n,l2:n*(n-1)}))
zero("all-degree b leading coefficient",
     leadb-245*d0**3*n*(13*n-9)/(8*J0**2))
zero("all-degree final polynomial equation leading coefficient",
     leadb*(10*n+45*(2*n-2))
     -1225*d0**3*n*(10*n-9)*(13*n-9)/(4*J0**2))

# A persistent reflected contact forces opposite source Jacobians.
fcontact=(T-1)**2*(T-2)
P=u*u+y
Q=u*fcontact.subs(T,u*u)+y*y
Jac=sp.expand(sp.diff(P,u)*sp.diff(Q,y)-sp.diff(P,y)*sp.diff(Q,u))
zero("contact first endpoint values match",P.subs(u,1)-P.subs(u,-1))
zero("contact second endpoint values match",Q.subs(u,1)-Q.subs(u,-1))
zero("contact source-u derivatives are opposite",
     sp.diff(Q,u).subs(u,1)+sp.diff(Q,u).subs(u,-1))
zero("contact source-y derivatives agree",
     sp.diff(Q,y).subs(u,1)-sp.diff(Q,y).subs(u,-1))
zero("persistent contact Jacobians are opposite",Jac.subs(u,1)+Jac.subs(u,-1))
check("persistent contact control can have nonzero opposite Jacobians",
      Jac.subs({u:1,y:1})!=0 and Jac.subs({u:-1,y:1})!=0)

# Canonical two-sided conductor and parity of all functions below odd degree seven.
fsep=(T-1)*(T-4)*(T-9)
for degree in range(7):
    poly=u**degree
    odd=sp.expand((poly-poly.subs(u,-u))/2)
    in_order=(sp.rem(odd/u if odd!=0 else 0,fsep.subs(T,u*u),u)==0)
    check(f"source degree below seven belongs only if even degree={degree}",
          in_order==(degree%2==0))
Cquotient=sp.rem(C/T,g,T)
zero("universal nonzero canonical period-vector remainder",
     Cquotient+5*b*T*T+14*c*T+35*d)
check("canonical period-vector remainder has nonzero constant for d unit",
      sp.Poly(Cquotient,T).nth(0)==-35*d)

# Sharp non-Keller family passing the geometric sextic constraints.
gf=(y-1)*u**6+u**4+u*u+1
pf=u*u
qf=u*gf
F=sp.symbols("V")**2-sp.symbols("U")*(
    (y-1)*sp.symbols("U")**3+sp.symbols("U")**2+sp.symbols("U")+1)**2
U,V=sp.symbols("U V")
zero("unmonic sextic family exact plane equation",
     F.subs({U:pf,V:qf},simultaneous=True))
zero("unmonic sextic family exact conductor",
     sp.diff(F,V).subs({U:pf,V:qf},simultaneous=True)/sp.diff(pf,u)-gf)
zero("unmonic sextic family Jacobian is non-Keller",
     sp.diff(pf,u)*sp.diff(qf,y)-sp.diff(pf,y)*sp.diff(qf,u)-2*u**8)
check("sextic immersion family has no fibrewise normalization ramification",
      sp.gcd(sp.diff(pf,u),sp.diff(qf,u))==1)
check("sextic family retains the special degree-four drop",
      sp.degree(gf,u)==6 and sp.degree(gf.subs(y,1),u)==4)

# Primary cable form transcription has a useful exact square rewrite.
X,Y,aa,bb,cc,dd=sp.symbols("X Y aa bb cc dd")
H=X**3+Y*Y+aa*X/2+bb/2
f1=(X**3+Y*Y)**2+aa*X**4+bb*X**3+aa*aa*X*X/4+cc*X+aa*X*Y*Y+bb*Y*Y+dd*Y
zero("primary cable polynomial exact transcription",
     f1-H**2-(cc-aa*bb/2)*X-dd*Y+bb*bb/4)

audit_path=HERE/"sources/neumann-genus-three-source-audit.json"
audit=json.loads(audit_path.read_text())
source_hashes={}
for source in audit["sources"]:
    path=HERE/"sources"/source["local_pdf"]
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    check("pinned primary source "+source["local_pdf"],actual==source["sha256"])
    source_hashes[source["local_pdf"]]=actual
proof=HERE/"HYPERELLIPTIC_SEXTIC_CONDUCTORS_AND_THREE_PERIOD_RIGIDITY.md"
record={
    "outcome":"PASS","check_count":len(checks),"checks":checks,
    "source_sha256":source_hashes,
    "source_audit_sha256":hashlib.sha256(audit_path.read_bytes()).hexdigest(),
    "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
    "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "sympy_version":sp.__version__,
    "elapsed_seconds":round(time.perf_counter()-started,3),
    "scope":{
        "generic_hyperelliptic_genus_three_ambient_branch_excluded_by_written_argument":True,
        "persistent_hyperelliptic_contacts_excluded":True,
        "special_degree_drops_retained_in_theorem":True,
        "all_degree_six_conductors_excluded":False,
        "intrinsic_quartic_or_cable_genus_three_branches_excluded":False,
        "uniform_family_target_automorphism_assumed":False,
        "primary_classification_proved_by_symbolic_checker":False,
        "automatic_trace_integrality_or_JC2_proved":False,
        "independent_replay":False,"independent_mathematical_review":False,
        "formal_verification":False,"external_review":False
    }
}
outdir=Path(tempfile.mkdtemp(prefix="hyperelliptic-sextic-check-"))
receipt=outdir/"hyperelliptic_sextic_conductor_verification.json"
receipt.write_text(json.dumps(record,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","checks":len(checks),"receipt":str(receipt)}))
