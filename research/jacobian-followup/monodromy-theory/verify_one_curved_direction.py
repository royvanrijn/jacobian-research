"""Exact whole identities and sharp controls for KMFI127."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import sympy as S
import json,time,hashlib,platform
from one_curved_direction import architecture,chart_identities
from primitive_multiplier import feedback_form,logarithmic_residues
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def zero(x):return S.cancel(x)==0
def mz(M):return all(zero(e) for e in M)
generic=[]
for r,d in [(2,3),(3,2),(3,4),(4,3),(5,2)]:
    ps=list(S.symbols("p0:"+str(r)));u=ps[0];xs=ps[1:-1];t=ps[-1]
    aa=list(S.symbols("a0:"+str(r-2)));c=S.Symbol("c")
    q0,q1,q2,q3=S.symbols("q0 q1 q2 q3")
    Q=q0+sum(q1*x+q2*x*x for x in xs)
    if len(xs)>1:Q+=q3*xs[0]*xs[-1]
    if d>=3 and xs:Q+=q3*xs[0]**3
    data=architecture(ps,d,aa,c,Q);ids=chart_identities(data)
    A=data["A"];vv=data["v"];T=ids["T"]
    tag="generic rank"+str(r)+" degree"+str(d)
    check(tag+" whole arrow resolvent",mz((T*S.eye(r)-A)*ids["resolvent"]-vv))
    check(tag+" polynomial adjugate numerator",all(S.fraction(e)[1]==1 for e in ids["adjugate_vector"]))
    check(tag+" whole characteristic polynomial",zero(A.charpoly(T).as_expr()-ids["chi"]))
    C=S.Matrix.hstack(*[A**j*vv for j in range(r)])
    check(tag+" entire Krylov determinant",zero(C.det(method="domain-ge")-ids["krylov_determinant"]))
    cadj=ids["adjugate_vector"].subs(T,c)
    check(tag+" c-adjugate including repeated d",mz(cadj-S.Matrix([0]*(r-1)+[ids["K"]*t+ids["R"]])))
    for i,(a,x) in enumerate(zip(aa,xs)):
        coef=a*(a-1)*x*(a-c)*S.prod(a-b for j,b in enumerate(aa) if j!=i)
        expected=S.zeros(r,1);expected[i+1]=coef
        expected[-1]=S.cancel(coef*S.diff(Q,x)/(a-c))
        check(tag+" a-adjugate"+str(i)+" including resonance d",mz(ids["adjugate_vector"].subs(T,a)-expected))
    detfull=S.factor(data["V"].jacobian(ps).det(method="domain-ge"))
    check(tag+" whole homogeneous determinant power",zero(detfull-d*c*S.prod(aa)*u**(r*(d-1))))
    generic.append({"r":r,"d":d,"Q":str(Q),"K":str(ids["K"]),"R":str(ids["R"]),"chart_Krylov_determinant":str(ids["krylov_determinant"])})
# The coefficient contradiction is independent of rank and homogeneous degree.
a,b,c,k=S.symbols("a b c k")
check("every univariate nonlinear monomial has exact pole coefficient",
 zero(1+(k*(a-1)-(c-1))/(c-a)-(k-1)*(a-1)/(c-a)))
check("distinct weights remove every mixed second derivative",
 zero(1/(c-a)-1/(c-b)-(a-b)/((c-a)*(c-b))))
# Whole feedback and exact lambda consistency, including d resonances.
controls=[]
for r,d,aa,cval in [
 (3,3,[3],7),(4,3,[3,4],8),(5,2,[3,4,5],9),
 (4,3,[4,5],3),(4,3,[3,4],4),(4,3,[3,4],1),
 (3,4,[3],7)]:
    ps=list(S.symbols("u0:"+str(r)));u=ps[0];xs=ps[1:-1];t=ps[-1];lam=S.Symbol("lam")
    Q=1+sum(x+x*x for x in xs)
    if len(xs)>1:Q+=xs[0]*xs[-1]
    if d>=3:Q+=xs[0]**3
    data=architecture(ps,d,list(map(S.Integer,aa)),S.Integer(cval),Q)
    f=feedback_form(data["V"],ps,lam*u**(d-1))
    logs=logarithmic_residues(f,ps,lam)
    tag="numeric "+str((r,d,aa,cval))
    check(tag+" feedback solve whole source",mz(f["eta"]*f["C"]-f["moments"]))
    T=f["T"]
    check(tag+" whole feedback characteristic identity",
      zero((T-lam*u**(d-1))**r-f["chi"]+
        (f["eta"]*(T*S.eye(r)-f["A"]).adjugate()*data["V"])[0]))
    condition=logs["parameter_condition"].as_expr()
    check(tag+" isolated lambda condition kept",zero(condition-(lam-cval)**r))
    check(tag+" only candidate forces t-independence",zero(f["eta"][-1].subs(lam,cval)))
    check(tag+" consistency and residue degree bound",
       logs["parameter_condition"].degree()==r and all(S.degree(v,lam)<=r for v in logs["residues"]))
    controls.append({"r":r,"d":d,"weights":aa,"c":cval,"Q":str(Q),"condition":str(S.factor(condition))})
# Degenerate early eigenvalues must not escape via repeated/resonant parameters.
u,v,w,t,lam=S.symbols("u v w t lam");ps=[u,v,w,t]
degenerate=[]
for aa,cval in [([1,4],6),([3,3],6),([0,4],6)]:
    data=architecture(ps,3,list(map(S.Integer,aa)),S.Integer(cval),1+v*v+v*w+w**3)
    C=S.Matrix.hstack(*[data["A"]**j*data["v"] for j in range(4)])
    check("degenerate base weights rejected "+str(aa),zero(C.det(method="domain-ge")))
    degenerate.append({"weights":aa,"c":cval,"Krylov_singular":True})
# An affine projective Q really can survive for arbitrarily high displayed degree.
# Work over Q(c), c root of 4c^2-4c+d^2, using exact polynomial remainders.
positive=[]
for d in [2,3,5,10]:
    u,t,c,gam=S.symbols("u t c gam");ps=[u,t]
    mod=S.Poly(4*c*c-4*c+d*d,c)
    def reduce(e):
        num,den=S.fraction(S.cancel(e))
        return S.cancel(S.rem(S.Poly(num,c),mod).as_expr()/S.rem(S.Poly(den,c),mod).as_expr())
    def zmod(e):return reduce(e)==0
    data=architecture(ps,d,[],c,c-1);lam=S.Rational(d,2)+c;F=t+u
    h=lam*u**(d-1)
    X=S.Matrix([-F*h,*list(F*(data["V"]-h*S.Matrix(ps)))])
    P=(X[1:,:].jacobian(ps)-S.Matrix(ps)*X[:1,:].jacobian(ps)).applyfunc(S.expand)
    tag="accepted affine Q degree"+str(d)
    check(tag+" entire normalized nilpotency",all(zmod(e) for e in P**2))
    sweep=X+gam*S.Matrix([1,u,t])
    check(tag+" whole sweep Jacobian",zmod(sweep.jacobian([gam,u,t]).det(method="domain-ge")-gam**2))
    minor=X.jacobian(ps).extract([0,1],[0,1]).det()
    check(tag+" genuine full rank",not zmod(minor.subs({u:1,t:1})))
    gcd=S.gcd_list(list(data["V"]))
    check(tag+" displayed degree reduces to primitive1",zero(gcd-u**(d-1)))
    check(tag+" exact factorization H and f",mz(F*data["V"]-(F*u**(d-1))*S.Matrix([u,c*t+(c-1)*u])))
    positive.append({"d":d,"coefficient_minimal_polynomial":str(mod.as_expr()),"Q":"c-1","F":"t+u",
      "lambda":str(lam),"coordinate_gcd":str(gcd),"primitive_degree":1,"rankDX":2,
      "whole_Jacobian":"gam^2","new_Keller_action_claimed":False})
# Rank-deficient normalized multipliers remain possible for nonlinear Q.
u,v,w,t,gam=S.symbols("u v w t gam");ps=[u,v,w,t]
data=architecture(ps,2,[1+S.I,1-S.I],S.Integer(2),v*w)
F=v*w;h=2*u
X=S.Matrix([-F*h,*list(F*(data["V"]-h*S.Matrix(ps)))])
sweep=X+gam*S.Matrix([1,*ps])
check("sharp nonlinear Q admits normalized rank3 multiplier",
 zero(sweep.jacobian([gam,*ps]).det(method="domain-ge")-gam**4)
 and X.jacobian(ps).subs({p:1 for p in ps}).rank()==3
 and all(zero(S.diff(e,t)) for e in X))
model={"claims":["KMFI127"],"previous_goal_turn":"progress: KMFI124--126 established universal Euler cyclicity, a complete fixed-direction multiplier criterion and the all-parameter quadratic binary-shear obstruction.",
 "architecture":"V=(u^d,ai*u^(d-1)*xi,c*u^(d-1)*t+u^d*Q(x/u)), r>=2,d>=1,degQ<=d",
 "unbounded_result":"full rank forces Q affine; every survivor has primitive degree1 and is classified by KMFI121",
 "all_parameters_including_resonances":True,"generic_identity_controls":generic,
 "nonlinear_Q_feedback_controls":controls,"degenerate_weights":degenerate,
 "accepted_affine_Q_full_rank_controls":positive,
 "accepted_nonlinear_Q_rank_deficient_control":{"d":2,"r":4,"weights":["1+i","1-i"],"c":2,"Q":"v*w","F":"v*w","rankDX":3},
 "degree11_groups_after_full_rank":["S11","A11","regular C11"],
 "all_primitive_directions_classified":False,"multiple_curved_coordinates_classified":False,
 "new_abstract_action_claimed":False,"actual_PSL11_input_supplied":False,
 "absolute_exceptional_action_decided":False,"full_action_spectrum_complete":False,"goal_complete":False}
(ROOT/"one_curved_direction_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["ONE_CURVED_DIRECTION.md","one_curved_direction.py","verify_one_curved_direction.py",
 "one_curved_direction_model.json","PRIMITIVE_MULTIPLIER_CRITERION.md","primitive_multiplier.py",
 "FACTORIZED_HOMOGENEOUS_FIELDS.md","SWEEP_DIFFERENTIAL_RECONSTRUCTION.md"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,
 "runtime_seconds":round(time.monotonic()-started,3),"python":platform.python_version(),"sympy":S.__version__,
 "sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 "assurance":"Unbounded written all-rank/all-degree proof. Whole symbolic arrow, adjugate and Krylov identities; exact resonant, positive and rank-deficient controls. No formal verification, independent whole-proof replay or external review claimed.",
 "goal_complete":False}
(ROOT/"one_curved_direction_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],"all_rank_all_degree_architecture_closed":True,"goal_complete":False}))
