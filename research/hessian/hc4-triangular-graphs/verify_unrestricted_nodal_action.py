"""Exact controls for unrestricted nodal action and conductor-degree-two exclusion."""
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
def period(p,q):
    return sp.integrate(sp.expand(p*sp.diff(q,t)),(t,-1,1))

t,u,y,X,Y,r,a,U=sp.symbols("t u y X Y r a U")
s=t*t-1
z=t*(t*t-1)
zero("standard nodal normalization relation",z*z-s*s*(s+1))
zero("standard nodal action",period(s,z)-sp.Rational(8,15))
check("standard normalization is immersed",sp.gcd(sp.diff(s,t),sp.diff(z,t))==1)
for epsilon in [-1,1]:
    zero(f"normalization orientation changes action by sign {epsilon}",
         period(s,epsilon*z)-epsilon*sp.Rational(8,15))
    check(f"node preimages remain exactly plus and minus one {epsilon}",
          sp.gcd(s,epsilon*z)==t*t-1)

# Closed polynomial one-forms used for arbitrary target action scaling.
def primitive_closed(A,B):
    first=sp.integrate(A,X)
    return sp.expand(first+sp.integrate(B-sp.diff(first,Y),Y))

maps=[
    (X+Y**2,Y),
    (X,Y+X**4),
    (2*X+3*Y+1,5*X+7*Y-2),
    (Y,X),
    (X+Y**2,Y+(X+Y**2)**3),
]
map_controls=[]
for i,(F,G) in enumerate(maps):
    lam=sp.expand(sp.diff(F,X)*sp.diff(G,Y)-sp.diff(F,Y)*sp.diff(G,X))
    check(f"target control has constant nonzero Jacobian {i}",
          lam!=0 and not lam.has(X,Y))
    A=sp.expand(F*sp.diff(G,X))
    B=sp.expand(F*sp.diff(G,Y)-lam*X)
    zero(f"target action correction is closed {i}",sp.diff(A,Y)-sp.diff(B,X))
    H=primitive_closed(A,B)
    zero(f"target action primitive X coefficient {i}",sp.diff(H,X)-A)
    zero(f"target action primitive Y coefficient {i}",sp.diff(H,Y)-B)
    pp=sp.expand(F.subs({X:s,Y:z},simultaneous=True))
    qq=sp.expand(G.subs({X:s,Y:z},simultaneous=True))
    zero(f"exact target action scaling {i}",period(pp,qq)-lam*sp.Rational(8,15))
    map_controls.append({"jacobian":str(lam),"p_degree":int(sp.degree(pp,t)),
                         "q_degree":int(sp.degree(qq,t))})

# An actual nodal generator of degree SIX with partner degree unbounded.
# q-H(p)=z, p-z^2=s gives the polynomial inverse of the target shears.
six_controls=[]
for N in [1,2,5,9]:
    pp=sp.expand(z*z+s)
    qq=sp.expand(z+pp**N)
    zz=qq-pp**N
    ss=pp-zz**2
    zero(f"degree-six pair recovers cubic generator N={N}",zz-z)
    zero(f"degree-six pair recovers quadratic generator N={N}",ss-s)
    zero(f"degree-six unbounded-partner action N={N}",
         period(pp,qq)-sp.Rational(8,15))
    check(f"degree-six partner exact degree N={N}",
          sp.degree(pp,t)==6 and sp.degree(qq,t)==6*N)
    six_controls.append({"p_degree":6,"q_degree":6*N,"action":"8/15"})

# Exact endpoint variation, without assuming the Jacobian is constant.
P=(u*u-1)*(u**4+y*u*u+y*y*u+2)
Q=(u*u-1)*(u**3+y*y*u*u+y*u+1)+y
J=sp.expand(sp.diff(P,u)*sp.diff(Q,y)-sp.diff(P,y)*sp.diff(Q,u))
I=sp.integrate(sp.expand(P*sp.diff(Q,u)),(u,-1,1))
zero("generic endpoint P values match",P.subs(u,1)-P.subs(u,-1))
zero("generic endpoint Q values match",Q.subs(u,1)-Q.subs(u,-1))
boundary=sp.expand((P*sp.diff(Q,y)).subs(u,1)
                   -(P*sp.diff(Q,y)).subs(u,-1))
zero("endpoint variation boundary vanishes",boundary)
zero("universal integration-by-parts variation identity",
     sp.diff(I,y)-boundary+sp.integrate(J,(u,-1,1)))
J0=sp.symbols("J0",nonzero=True)
zero("constant nonzero Jacobian gives affine action slope",
     sp.integrate(-J0,(u,-1,1))+2*J0)

# Boundary matching is necessary: the affine Keller pair omits it.
P0,Q0=u,y
I0=sp.integrate(P0*sp.diff(Q0,u),(u,-1,1))
boundary0=(P0*sp.diff(Q0,y)).subs(u,1)-(P0*sp.diff(Q0,y)).subs(u,-1)
zero("affine Keller control has nonzero endpoint boundary",boundary0-2)
zero("affine Keller control obeys full variation identity",
     sp.diff(I0,y)-boundary0+2)
check("affine Keller control fails exact nodal endpoint hypothesis",
      P0.subs(u,1)!=P0.subs(u,-1))

# Repair the source's final cubic normalization using a double root.
v=sp.symbols("v")
cubic=Y*Y-(X-r)**2*(X-(r-a*a))
zero("universal double-root cubic normalizes by affine scaling",
     cubic.subs({X:r+a*a*u,Y:a**3*v},simultaneous=True)
     -a**6*(v*v-u*u*(u+1)))
f=X**3-Y**2+X**2
for xc,value in [(0,0),(-sp.Rational(2,3),sp.Rational(4,27))]:
    zero(f"source normalized cubic critical X derivative at {xc}",
         sp.diff(f,X).subs({X:xc,Y:0}))
    zero(f"source normalized cubic critical Y derivative at {xc}",
         sp.diff(f,Y).subs({X:xc,Y:0}))
    zero(f"source normalized cubic distinct critical value at {xc}",
         f.subs({X:xc,Y:0})-value)
    check(f"source normalized cubic critical point is a node at {xc}",
          sp.hessian(f,(X,Y)).det().subs({X:xc,Y:0})!=0)

# Zero action alone does not certify the exact one-node generator algebra.
pbad=s
qbad=t*s*(t*t-sp.Rational(3,7))
zero("two-node polynomial immersion has zero action",period(pbad,qbad))
check("two-node negative control is immersed",
      sp.gcd(sp.diff(pbad,t),sp.diff(qbad,t))==1)
W=sp.Matrix([[0,U+1],[1,0]])
Qbad=W*U*(W**2-sp.Rational(3,7)*sp.eye(2))
inc=sp.Matrix.hstack(sp.Matrix([1,0]),Qbad*sp.Matrix([1,0]))
zero("zero-action control has index degree two instead of one",
     inc.det()-U*(U+sp.Rational(4,7)))
check("zero-action control has an extra node",
      sp.degree(inc.det(),U)==2)

audit_path=HERE/"sources/neumann-nodal-classification-source-audit.json"
audit=json.loads(audit_path.read_text())
source_hashes={}
for source in audit["sources"]:
    source_path=HERE/"sources"/source["local_pdf"]
    actual=hashlib.sha256(source_path.read_bytes()).hexdigest()
    check("pinned primary PDF hash "+source["local_pdf"],actual==source["sha256"])
    source_hashes[source["local_pdf"]]=actual
check("source audit does not claim symbolic proof of classification",
      audit["source_classification_proved_by_symbolic_checker"] is False)

proof=HERE/"UNRESTRICTED_NODAL_ACTION_AND_DEGREE_TWO_CONDUCTOR_EXCLUSION.md"
record={
    "outcome":"PASS","check_count":len(checks),"checks":checks,
    "target_automorphism_controls":map_controls,
    "degree_six_unbounded_partner_controls":six_controls,
    "source_sha256":source_hashes,
    "source_audit_sha256":hashlib.sha256(audit_path.read_bytes()).hexdigest(),
    "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
    "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "sympy_version":sp.__version__,
    "elapsed_seconds":round(time.perf_counter()-started,3),
    "scope":{
        "unrestricted_nodal_action_uses_audited_primary_classification":True,
        "degree_two_keller_conductor_excluded_by_written_argument":True,
        "unrestricted_classification_proved_by_finite_controls":False,
        "primary_source_critical_value_slip_repaired":True,
        "automatic_trace_integrality_or_JC2_proved":False,
        "independent_replay":False,"independent_mathematical_review":False,
        "formal_verification":False,"external_review":False
    }
}
outdir=Path(tempfile.mkdtemp(prefix="unrestricted-nodal-action-check-"))
receipt=outdir/"unrestricted_nodal_action_verification.json"
receipt.write_text(json.dumps(record,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","checks":len(checks),"receipt":str(receipt)}))
