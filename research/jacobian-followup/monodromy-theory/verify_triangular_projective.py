"""Whole identities and sharp exact controls for KMFI128--129."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import sympy as S
import json,time,platform,hashlib
from triangular_projective import architecture,resolvent,last_coordinate,rigidity_residuals
from primitive_multiplier import feedback_form,logarithmic_residues
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def zero(e):return S.cancel(e)==0
def mz(M):return all(zero(e) for e in M)
generic=[]
for r,d in [(2,3),(3,2),(4,3),(5,2)]:
    ps=list(S.symbols("p0:"+str(r)));u=ps[0];xs=ps[1:];n=len(xs)
    aa=list(S.symbols("a0:"+str(n)));qq,rr=S.symbols("q r")
    Q=[qq]
    for i in range(1,n):
        e=rr*xs[i-1]**2+qq*xs[0]
        if i>1:e+=qq*xs[0]*xs[i-1]
        if d>=3:e+=rr*xs[0]**3
        Q.append(e)
    data=architecture(ps,d,aa,Q);A=data["A"];V=data["chart_V"];B=data["B"];W=data["W"]
    T=S.Symbol("T");z=resolvent(data,T);tag="generic rank"+str(r)+" degree"+str(d)
    check(tag+" entire projective resolvent",mz((T*S.eye(r)-A)*z-V))
    check(tag+" whole chart derivative block",mz(A-S.BlockMatrix([[S.Matrix([[d]]),S.zeros(1,n)], [d*data["v"]-B*S.Matrix(xs),B]]).as_explicit()))
    chi=(T-d)*S.prod(T-a for a in aa)
    check(tag+" exact constant spectrum",zero(A.charpoly(T).as_expr()-chi))
    C=S.Matrix.hstack(*[A**j*V for j in range(r)])
    U=S.eye(r)
    for j in range(1,r):U[j-1,j]=-d
    D=C*U
    expected=S.Matrix.hstack(V,*[S.Matrix([0,*list(B**j*W)]) for j in range(1,r)])
    check(tag+" unit Krylov column transformation",U.det()==1 and mz(D-expected))
    if r<=4:
        K=S.Matrix.hstack(*[B**j*W for j in range(n)])
        check(tag+" whole cyclic determinant identity",zero(C.det(method="domain-ge")-B.det()*K.det(method="domain-ge")))
    top=last_coordinate(data,n-1)
    adj=(chi*z).applyfunc(S.cancel)
    check(tag+" polynomial adjugate numerator",all(S.fraction(e)[1]==1 for e in adj))
    check(tag+" whole last eigenvalue adjugate through d resonance",
      mz(adj.subs(T,aa[-1])-S.Matrix([0]*n+[top["L"]])))
    ell=top["ell"];q=top["q"];c=top["c"]
    check(tag+" entire left eigenrow",mz(ell*B-c*ell))
    check(tag+" entire spectral graph equation",zero((ell*W)[0]-(c-1)*q))
    check(tag+" last projection coefficient exact",zero(S.diff(top["L"],xs[-1])-top["K"]))
    generic.append({"r":r,"d":d,"Q":[str(e) for e in Q],"last_graph":str(q),
                    "last_adjugate":str(top["L"]),"Krylov_unit_column_identity":True})
# Universal differentiated identities behind the affine eigen-coordinate lemma.
x,y,t,c,b=S.symbols("x y t c b")
v=S.Matrix([x*x+y,x*y+t,y*y+x*t]);q=x*x*y+y*t+t**3
rr=rigidity_residuals(v,[x,y,t],q,c)
check("universal Hessian self-adjoint residual identity",mz(rr["symmetric_identity"]))
check("universal differentiated Darboux residual identity",mz(rr["differentiated_identity"]))
sharp=rigidity_residuals(S.Matrix([x,2*t+x*x]),[x,t],t+x*x,S.Integer(2))
check("cyclicity hypothesis essential nonlinear eigen-coordinate",
 mz(sharp["R"]) and zero(sharp["Z"]) and not mz(sharp["H"]) and
 zero(S.Matrix.hstack(sharp["W"],sharp["B"]*sharp["W"]).det()))
# Coupled nonlinear feedback: isolated parameter conditions are kept exactly.
u,v,w,t,lam=S.symbols("u v w t lam");ps=[u,v,w,t]
specs=[
 ("coupled_quadratic",2,[3,4,5],[0,v*v,v*w+w*w],S.Integer(1)),
 ("early_curvature_with_diagonal_spectator",2,[3,4,5],[0,v*v,0],(lam-4)**4),
 ("repeated_last",2,[3,4,4],[0,v*v,v*w],(lam-4)**3),
 ("repeated_earlier",2,[3,3,5],[0,v*v,v*w],S.Integer(1)),
 ("coupled_cubic_resonance_a_equals_d",3,[3,4,5],[0,v**3,v*w*w],S.Integer(1)),
 ("coupled_lower_terms",3,[3,4,6],[1,v*v+v+1,v*w+w*w+v+1],S.Integer(1))]
negative=[]
for name,d,aa,Q,condition in specs:
    data=architecture(ps,d,list(map(S.Integer,aa)),Q);f=feedback_form(data["V"],ps,lam*u**(d-1))
    logs=logarithmic_residues(f,ps,lam);T=f["T"]
    check(name+" exact entire feedback solve",mz(f["eta"]*f["C"]-f["moments"]))
    check(name+" entire feedback characteristic",
      zero((T-lam*u**(d-1))**4-f["chi"]+(f["eta"]*(T*S.eye(4)-f["A"]).adjugate()*data["V"])[0]))
    check(name+" exact lambda consistency",zero(logs["parameter_condition"].as_expr()-condition))
    if condition!=1:
        at={u:1,v:1,w:1,t:1,lam:4}
        rho=f["rho"].subs(at);N=f["N"].subs(at)
        check(name+" only compatible candidate has rank collapse",zero(S.Matrix.vstack(*[rho*N**j for j in range(4)]).det()))
    negative.append({"name":name,"r":4,"d":d,"weights":aa,"Q":[str(e) for e in Q],
                    "condition":str(S.factor(logs["parameter_condition"].as_expr()))})
# Rank5 coupled curvature, beyond one curved coordinate and rank4.
pp=list(S.symbols("y0:5"));u0=pp[0];xx=pp[1:]
data=architecture(pp,2,list(map(S.Integer,[3,4,5,6])),[0,xx[0]**2,xx[0]*xx[1],xx[1]*xx[2]])
f=feedback_form(data["V"],pp,lam*u0);logs=logarithmic_residues(f,pp,lam)
check("rank5 several curved coordinates reject every lambda",logs["parameter_condition"].degree()==0)
negative.append({"name":"rank5_coupled","r":5,"condition":str(logs["parameter_condition"].as_expr())})
# Affine projective couplings survive honestly, after a constant linear change.
u,v,w,gam=S.symbols("u v w gam");ps=[u,v,w];L=S.Matrix([[1,0,0],[1,1,0],[3,2,1]])
M=(L.inv()*S.diag(1,1+S.I,1-S.I)*L).applyfunc(S.expand)
V=u*M*S.Matrix(ps);F=(u+v)*(3*u+2*v+w);h=2*u
X=S.Matrix([-F*h,*list(F*(V-h*S.Matrix(ps)))])
P=(X[1:,:].jacobian(ps)-S.Matrix(ps)*X[:1,:].jacobian(ps)).applyfunc(S.expand)
sweep=X+gam*S.Matrix([1,*ps])
check("affine coupled survivor whole normalized nilpotency",mz((P**3).applyfunc(S.expand)))
check("affine coupled survivor whole sweep Jacobian",zero(sweep.jacobian([gam,*ps]).det(method="domain-ge")-gam**3))
check("affine coupled survivor full rank",X.jacobian(ps).subs({u:1,v:0,w:0}).rank()==3)
check("affine coupled survivor reduces to primitive1",zero(S.gcd_list(list(V))-u))
data=architecture(ps,2,[M[1,1],M[2,2]],[M[1,0],M[2,0]+M[2,1]*v])
check("affine survivor belongs to exact triangular architecture",mz(data["V"]-V))
top=last_coordinate(data,1)
check("surviving eigen-coordinate really affine",mz(S.hessian(top["q"],[v,w])) and mz(S.Matrix([[S.diff(top["q"],v),S.diff(top["q"],w)]])-top["ell"]))
# Nonlinear coupled inputs can remain normalized and have rank r-1.
u,v,w,t,gam=S.symbols("u v w t gam");ps=[u,v,w,t]
data=architecture(ps,2,[S.Integer(2)]*3,[0,v*v,v*w]);F=S.Integer(1);h=2*u
X=S.Matrix([-h,*list(data["V"]-h*S.Matrix(ps))]);sweep=X+gam*S.Matrix([1,*ps])
check("nonlinear coupled rank3 control whole normalized Jacobian",
 zero(sweep.jacobian([gam,*ps]).det(method="domain-ge")-gam**4))
check("nonlinear coupled rank3 control preserves full-rank boundary",
 X.jacobian(ps).subs({p:1 for p in ps}).rank()==3 and all(zero(S.diff(e,t)) for e in X))
# a_i=d and repeated spectra in the whole adjugate formula, with no division by eigenvalue differences.
for aa in [[2,3,4],[3,4,3],[3,3,3],[3,4,1]]:
    data=architecture(ps,3,list(map(S.Integer,aa)),[0,v*v,v*w+w*w]);top=last_coordinate(data,2);T=S.Symbol("T")
    adj=(T*S.eye(4)-data["A"]).adjugate()*data["chart_V"]
    check("whole repeated/resonant adjugate "+str(aa),
      mz(adj.subs(T,aa[-1])-S.Matrix([0,0,0,top["L"]])))
model={"claims":["KMFI128","KMFI129"],
 "previous_goal_turn":"progress: KMFI127 closed all homogeneous one-curved-coordinate architectures in every rank and degree.",
 "affine_eigen_coordinate_lemma":{"hypotheses":"B=Dv,W=v-x cyclic; dq*B=c*dq; W(q)=(c-1)q, constantc","conclusion":"q affine","triangularity_required":False,"cyclicity_countercontrol_retained":True},
 "homogeneous_architecture":"V0=u^d,Vi=ai*u^(d-1)*xi+u^d*Qi(x1/u,...,x(i-1)/u); eachdegQi<=d",
 "unbounded_result":"all Qi affine under normalization and full rank; primitive degree1",
 "all_ranks_and_all_scalar_and_direction_degrees":True,"repeated_and_resonant_parameters_covered":True,
 "generic_whole_identity_controls":generic,"nonlinear_feedback_controls":negative,
 "affine_coupled_full_rank_survivor":{"r":3,"d":2,"F":"(u+v)*(3*u+2*v+w)","h":"2*u","linear_frame":[list(L.row(i)) for i in range(3)],"rankDX":3,"whole_Jacobian":"gam^3"},
 "nonlinear_coupled_rank_deficient_control":{"r":4,"d":2,"weights":[2,2,2],"Q":["0","v^2","v*w"],"F":"1","rankDX":3,"whole_Jacobian":"gam^4"},
 "degree11_groups":["S11","A11","regular C11"],
 "nonconstant_projective_diagonal_classified":False,"arbitrary_nontriangular_directions_classified":False,
 "new_abstract_action_claimed":False,"actual_PSL11_input_supplied":False,
 "absolute_exceptional_action_decided":False,"full_action_spectrum_complete":False,"goal_complete":False}
# All JSON coefficients are integers or displayed exact strings.
model["affine_coupled_full_rank_survivor"]["linear_frame"]=[[int(e) for e in L.row(i)] for i in range(3)]
(ROOT/"triangular_projective_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["TRIANGULAR_PROJECTIVE_RIGIDITY.md","triangular_projective.py","verify_triangular_projective.py",
 "triangular_projective_model.json","PRIMITIVE_MULTIPLIER_CRITERION.md","primitive_multiplier.py",
 "FACTORIZED_HOMOGENEOUS_FIELDS.md","ONE_CURVED_DIRECTION.md"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,"runtime_seconds":round(time.monotonic()-started,3),
 "python":platform.python_version(),"sympy":S.__version__,"sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 "assurance":"Unbounded written affine eigen-coordinate lemma and all-rank/all-degree triangular projective rigidity proof. Whole identities and sharp controls supplement the proof; no formal verification, independent whole-proof replay or external review claimed.","goal_complete":False}
(ROOT/"triangular_projective_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],"entire_constant_diagonal_triangular_architecture_closed":True,"goal_complete":False}))
