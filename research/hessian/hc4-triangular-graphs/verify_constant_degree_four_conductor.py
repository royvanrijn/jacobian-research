"""Exact constant-degree-four conductor and two-period ODE controls."""
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

u,y,U,V,r,s,a,K,S,S1,S2,Lam,Lam1,J=sp.symbols(
    "u y U V r s a K S S1 S2 Lam Lam1 J")
g=u**4-S*u*u+K
p=u*u
q=u*g
f=V*V-U*(U*U-S*U+K)**2
zero("canonical full-family equation",f.subs({U:p,V:q},simultaneous=True))
zero("canonical full-family exact conductor ratio",
     sp.diff(f,V).subs({U:p,V:q},simultaneous=True)/sp.diff(p,u)-g)
index=U*U-S*U+K
zero("canonical index degree two",sp.degree(index,U)-2)
zero("canonical quartic discriminant",
     sp.discriminant(g,u)-16*K*(S*S-4*K)**2)
zero("two-node conductor when K=r^2*s^2,S=r^2+s^2",
     g.subs({S:r*r+s*s,K:r*r*s*s})-(u*u-r*r)*(u*u-s*s))
zero("tacnode conductor when S=2r^2,K=r^4",
     g.subs({S:2*r*r,K:r**4})-(u*u-r*r)**2)
zero("tacnode quotient equation",
     f.subs({S:2*r*r,K:r**4})-(V*V-U*(U-r*r)**4))

def action(pp,qq,v):
    return sp.integrate(sp.expand(pp*sp.diff(qq,u)),(u,-v,v))
qnodes=q.subs({S:r*r+s*s,K:r*r*s*s})
Ar=action(p,qnodes,r)
As=action(p,qnodes,s)
zero("first canonical moving-node period",
     Ar-sp.Rational(8,105)*r**5*(3*r*r-7*s*s))
zero("second canonical moving-node period",
     As-sp.Rational(8,105)*s**5*(3*s*s-7*r*r))
zero("fixed tacnode period",
     action(p,u*(u*u-r*r)**2,r)+sp.Rational(32,105)*r**7)
zero("canonical action normalized by r",
     Ar.subs(s*s,K/(r*r))-sp.Rational(8,105)*r*(3*r**6-7*K*r*r))

# Derive the equations from the two separate period derivatives.
b=K/a
St=a+b
d=a-b
N=3*St**2-10*K
Ca=3*a**3-7*K*a
Cb=3*b**3-7*K*b
aprime=a*S1/d
zero("paired root logarithmic derivative identity",
     sp.diff(b,a)*aprime+b*S1/d)
zero("canonical odd parts difference",Ca-Cb-d*N)
zero("canonical odd parts sum",Ca+Cb-St*(3*St**2-16*K))
E1=Lam1*Ca+Lam*(sp.diff(Ca,a)*aprime+S1*Ca/(2*d))+2*J
E2=Lam1*Cb+Lam*(sp.diff(Cb,a)*aprime-S1*Cb/(2*d))+2*J
eqdiff=2*N*Lam1+21*St*S1*Lam
zero("period-equation difference gives first Lambda equation",
     E1-E2-d*eqdiff/2)
eqsum=St*(3*St**2-16*K)*Lam1+sp.Rational(21,2)*(St**2-2*K)*S1*Lam+4*J
zero("period-equation sum",E1+E2-eqsum)
zero("exact constant coefficient after elimination",
     N*(St**2-2*K)-St**2*(3*St**2-16*K)-20*K*K)
zero("period-equation Lambda elimination",
     2*N*eqsum-St*(3*St**2-16*K)*eqdiff
     -2*(210*K*K*Lam*S1+4*J*N))

# Independently verify differentiation of Lambda*S' and final ODE.
Nfree=3*S*S-10*K
diff_equation=210*K*K*(Lam1*S1+Lam*S2)+24*J*S*S1
substituted=diff_equation.subs({
    Lam1:-21*S*S1*Lam/(2*Nfree),
    J:-210*K*K*Lam*S1/(4*Nfree)},simultaneous=True)
zero("differentiated scale equation gives final polynomial ODE",
     substituted-105*K*K*Lam/Nfree*(2*Nfree*S2-33*S*S1*S1))
n=sp.symbols("n",integer=True,positive=True)
zero("all-degree polynomial leading coefficient factor",
     6*n*(n-1)-33*n*n+3*n*(9*n+2))
for degree in [1,2,3,6]:
    Sm=y**degree+2*y+3
    if degree==1:
        Sm=y+3
    ode=sp.Poly(sp.expand(2*(3*Sm**2-10*K)*sp.diff(Sm,y,2)
                          -33*Sm*sp.diff(Sm,y)**2),y)
    check(f"nonconstant polynomial ODE model has expected degree {degree}",
          ode.degree()==3*degree-2)
    zero(f"nonconstant polynomial ODE leading coefficient {degree}",
         ode.LC()+3*degree*(9*degree+2))

# Exact moving-endpoint identity in a genuine canonical family.
# Set K=1, S=y on the algebraic root cover y=r^2+r^-2.
qvar=u*(u**4-y*u*u+1)
Jvar=sp.expand(sp.diff(p,u)*sp.diff(qvar,y))
zero("non-Keller immersion-family Jacobian",Jvar+2*u**4)
check("entire canonical family is immersed on every line",
      sp.gcd(sp.diff(p,u),sp.diff(qvar,u))==1)
ys=r*r+r**-2
rprime=1/sp.diff(ys,r)
Ivar=action(p,qvar,r).subs(y,ys)
zero("moving endpoint action derivative keeps boundary cancellation",
     sp.diff(Ivar,r)*rprime
     +sp.integrate(Jvar,(u,-r,r)).subs(y,ys))
for yc in [2,-2]:
    gc=sp.expand(g.subs({S:yc,K:1}))
    check(f"constant-degree control has a repeated-root fibre at y={yc}",
          sp.discriminant(gc,u)==0)
    check(f"constant-degree control does not lose conductor degree at y={yc}",
          sp.degree(gc,u)==4)
for yc in [0,1,3]:
    check(f"constant-degree control has four simple roots at y={yc}",
          sp.discriminant(g.subs({S:yc,K:1}),u)!=0)

# Immersion is essential for the nonzero constant-term gate.
qcusp=u*(u**4-y*u*u)
zero("zero constant term creates a ramified normalization point",
     sp.diff(qcusp,u).subs(u,0))
zero("even coordinate derivative also vanishes at the cusp",
     sp.diff(p,u).subs(u,0))
check("zero-product control still has conductor degree four",
      sp.degree(g.subs({S:y,K:0}),u)==4)

audit_path=HERE/"sources/neumann-two-node-genus-two-source-audit.json"
audit=json.loads(audit_path.read_text())
source_hashes={}
for source in audit["sources"]:
    path=HERE/"sources"/source["local_pdf"]
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    check("pinned primary source "+source["local_pdf"],actual==source["sha256"])
    source_hashes[source["local_pdf"]]=actual
proof=HERE/"CONSTANT_DEGREE_FOUR_CONDUCTORS_AND_PERIOD_ODE.md"
record={
    "outcome":"PASS","check_count":len(checks),"checks":checks,
    "source_sha256":source_hashes,
    "source_audit_sha256":hashlib.sha256(audit_path.read_bytes()).hexdigest(),
    "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
    "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "sympy_version":sp.__version__,
    "elapsed_seconds":round(time.perf_counter()-started,3),
    "scope":{
        "constant_fibre_degree_four_conductors_excluded_by_written_argument":True,
        "repeated_roots_and_tacnodes_allowed_in_theorem":True,
        "all_global_degree_four_conductors_excluded":False,
        "degree_four_to_degree_two_specialization_excluded":False,
        "period_ODE_polynomial_identity_checked":True,
        "genus_two_classification_uses_audited_primary_source":True,
        "uniform_target_automorphism_over_parameter_assumed":False,
        "classification_proved_by_symbolic_checker":False,
        "automatic_trace_integrality_or_JC2_proved":False,
        "independent_replay":False,"independent_mathematical_review":False,
        "formal_verification":False,"external_review":False
    }
}
outdir=Path(tempfile.mkdtemp(prefix="constant-degree-four-check-"))
receipt=outdir/"constant_degree_four_conductor_verification.json"
receipt.write_text(json.dumps(record,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","checks":len(checks),"receipt":str(receipt)}))
