"""Exact nodal-action, cubic-index and Keller variation controls."""
from pathlib import Path
import hashlib,json,tempfile,time
import sympy as s

HERE=Path(__file__).resolve().parent
started=time.perf_counter()
checks=[]
def check(name,test):
    if not test: raise AssertionError(name)
    checks.append(name)
def zero(name,expr):
    check(name,s.cancel(expr)==0)
def action(p,q,t):
    return s.integrate(s.expand(p*s.diff(q,t)),(t,-1,1))

t,y,U,V,b,B,C,H,a,beta=s.symbols("t y U V b B C H a beta")
# Universal cubic inclusion determinant and norm, independently formed.
W=s.Matrix([[0,0,b+U],[1,0,1],[0,1,-b]])
Q=H*s.eye(3)+B*W+C*W**2
e0=s.eye(3)[:,0]
inclusion=s.Matrix.hstack(e0,Q*e0,Q**2*e0)
Delta=s.factor(inclusion.det())
norm=s.factor(((B-b*C)*s.eye(3)-C*W).det())
zero("universal cubic index equals polynomial norm",Delta-norm)
zero("universal cubic explicit index form",
     Delta-(B**3-2*b*B**2*C+(b*b-1)*B*C**2-U*C**3))
zero("nodal cubic remainder gives exact index degree one",
     Delta.subs(B,0)+C**3*U)

# Exact all-parameter actions, with no bound on the omitted H(p) term.
p2=a*(t*t-1)
q2=beta*t*(t*t-1)
zero("quadratic generator exact action",action(p2,q2,t)-8*a*beta/s.Integer(15))
p3=(t*t-1)*(a*t+b)
q3=C*t*t
zero("cubic generator exact action",action(p3,q3,t)+8*a*C/s.Integer(15))
for p in [t*t-1,(t*t-1)*(t+2)]:
    for exponent in [2,7]:
        hpoly=p**exponent+2*p*p-3
        zero(f"arbitrarily higher generator term has zero action d={s.degree(p,t)},H={exponent}",
             action(p,hpoly,t))
    check("common node values retained",p.subs(t,1)==p.subs(t,-1))

# Independent norm degree controls include potential leading cancellations
# between B and b*C. Different residue classes prevent cancellation with t*C.
for m,n in [(0,1),(1,0),(1,1),(2,2),(3,1)]:
    pp=(t*t-1)*(t+2)
    Bu=U**m
    Cu=U**n
    r=(Bu-2*Cu).subs(U,pp)-t*Cu.subs(U,pp)
    N=s.factor(norm.subs({b:2,B:Bu,C:Cu}))
    check(f"cubic norm degree formula m={m},n={n}",
          s.degree(N,U)==s.degree(r,t))
    check(f"nonconstant cubic remainder cannot have index one m={m},n={n}",
          s.degree(N,U)>1)

# Polynomial target automorphism covariance, including both degrees high.
p=t*t-1
q=t*(t*t-1)
I0=action(p,q,t)
zero("standard nodal action",I0-s.Rational(8,15))
F=U+V*V
G=V+F**3
jac=s.det(s.Matrix([F,G]).jacobian([U,V]))
zero("nonlinear target automorphism Jacobian",jac-1)
pnew=s.expand(F.subs({U:p,V:q},simultaneous=True))
qnew=s.expand(G.subs({U:p,V:q},simultaneous=True))
zero("nonlinear target action covariance",action(pnew,qnew,t)-I0)
check("covariance control has both generator degrees above three",
      s.degree(pnew,t)>3 and s.degree(qnew,t)>3)
zero("nonunit constant target determinant covariance",
     action(2*pnew,3*qnew,t)-6*I0)
# Verify its polynomial inverse before using it as an ambient automorphism.
invV=V-U**3
invU=U-invV**2
zero("target automorphism inverse first coordinate",
     invU.subs({U:F,V:G},simultaneous=True)-U)
zero("target automorphism inverse second coordinate",
     invV.subs({U:F,V:G},simultaneous=True)-V)

# Universal action variation is just the exact product identity plus the
# common endpoint values; this test uses actual two-variable polynomials.
PP=y+(t*t-1)*(t+y*y)
QQ=y*y+(t*t-1)*(t*t+y*t+1)
J=s.det(s.Matrix([PP,QQ]).jacobian([t,y]))
II=action(PP,QQ,t)
zero("action variation with identical node endpoint functions",
     s.diff(II,y)+s.integrate(J,(t,-1,1)))
zero("variation boundary P values agree",PP.subs(t,1)-PP.subs(t,-1))
zero("variation boundary Q_y values agree",
     s.diff(QQ,y).subs(t,1)-s.diff(QQ,y).subs(t,-1))
zero("pointwise universal variation transcription",
     s.diff(PP*s.diff(QQ,t),y)+J-s.diff(PP*s.diff(QQ,y),t))

# Deliberately false generator hypothesis: zero action but extra node.
badp=t*t-1
badq=t*(t*t-1)*(t*t-s.Rational(3,7))
zero("zero-action nongenerating immersion control",action(badp,badq,t))
check("nongenerating control is everywhere an immersion",
      s.gcd(s.diff(badp,t),s.diff(badq,t))==1)
badindex=U*(U+s.Rational(4,7))
check("nongenerating control has index dimension two",
      s.degree(badindex,U)==2)
extra=s.sqrt(s.Rational(3,7))
zero("additional collision P values",badp.subs(t,extra)-badp.subs(t,-extra))
zero("additional collision Q values",badq.subs(t,extra)-badq.subs(t,-extra))
check("zero-action pair cannot generate one-node algebra",
      extra!=1 and extra!=-1 and s.degree(badindex,U)!=1)
# Direct basis inclusion over C[p] certifies the index used above.
W2=s.Matrix([[0,U+1],[1,0]])
badQ=W2*U*(U+s.Rational(4,7))
badC=s.Matrix.hstack(s.eye(2)[:,0],badQ*s.eye(2)[:,0])
zero("independent bad-pair inclusion determinant",badC.det()-badindex)

proof=HERE/"LOW_DEGREE_NODAL_ACTION_AND_CONDUCTOR_EXCLUSION.md"
record={
    "outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":round(time.perf_counter()-started,3),
    "sympy_version":s.__version__,
    "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
    "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "scope":{
        "all_degree_quadratic_generator_action_proved":True,
        "all_degree_cubic_generator_action_proved":True,
        "other_generator_degree_is_unbounded":True,
        "polynomial_target_automorphism_action_covariance_proved":True,
        "actual_degree_two_conductor_low_directional_degree_exclusion_proved":True,
        "zero_action_nongenerating_immersion_control_supplied":True,
        "all_nodal_generator_pairs_ambient_equivalent":False,
        "unrestricted_nodal_action_nonvanishing_proved":False,
        "entire_degree_two_conductor_case_excluded":False,
        "automatic_trace_integrality_or_JC2_proved":False,
        "formal_verification":False,"external_review":False
    }
}
outdir=Path(tempfile.mkdtemp(prefix="low-degree-nodal-action-check-"))
receipt=outdir/"low_degree_nodal_action_verification.json"
receipt.write_text(json.dumps(record,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","checks":len(checks),"receipt":str(receipt)}))
