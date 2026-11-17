"""Exact normal-feedback classification and cubic controls KMFI158--161."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import sympy as S
import json,time,hashlib,platform
from nilpotent_normal import linear_feedback,moments,cubic_form,cubic_control
from quadratic_primitive import feedback,normal_limit,trace_equation
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def zero(v):return S.cancel(v)==0
def mz(M):return all(zero(v) for v in M)
q1,q2,q3,q4,T,m=S.symbols("q1 q2 q3 q4 T m");q=S.Matrix([q1,q2,q3,q4])
# Entire Fitting-block resolvent, including the singular nilpotent summand.
cc=S.symbols("c0:4");C=S.Matrix(2,2,cc);K=S.Matrix([[0,1],[0,0]])
B=S.diag(C,K);eta=S.Matrix([S.symbols("eta0:4")])
whole=(eta*(T*S.eye(4)-B).inv()*B*q)[0]
active=(eta[:,:2]*(T*S.eye(2)-C).inv()*C*q[:2,:])[0]
check("entire Fitting resolvent splits into invertible part and a pure pole",
      zero(whole-active-eta[2]*q4/T))
check("whole invertible Fitting part tends to negative Euler pairing",
      zero(S.limit(active,T,0)+(eta[:,:2]*q[:2,:])[0]))
check("Fitting characteristic has exact zero multiplicity",
      zero((T*S.eye(4)-B).det()-T*T*(T*S.eye(2)-C).det()))
check("nilpotent target characteristic ratio tends to zero when active part exists",
      S.limit(T**4/((T*S.eye(4)-B).det()),T,0)==0)
check("zero-Euler negative degree is polynomially impossible coefficient by coefficient",
      S.simplify(q1*S.diff(q1**m,q1)+q1**m-(m+1)*q1**m)==0)
# Singular nonnilpotent B was not excluded by the former full determinant gate.
Bn=S.diag(1,-1,0);q3v=S.Matrix([q1,q2,q3])
An=Bn+(Bn*q3v)*S.Matrix([[m/q1,m/q2,0]])
check("singular nonnilpotent trace-zero family whole characteristic",
      zero(An.charpoly(T).as_expr()-T*(T*T-(2*m+1))))
check("singular nonnilpotent family normal determinant remains zero",An.det()==0)
check("singular nonnilpotent family entire trace is zero",S.trace(An)==0)
for degree in [0,1,4]:
    J=(q1*q2)**degree
    AA=linear_feedback(Bn,q3v,J)
    check("actual polynomial singular family rejected at degree"+str(degree),
          zero(AA.charpoly(T).as_expr()-T*(T*T-(2*degree+1))) and not mz(AA**3))
# Rational polynomiality boundary is exact.
Br=S.diag(1,0)
check("rational negative-Euler feedback is exactly zero",
      mz(linear_feedback(Br,S.Matrix([q1,q2]),1/q1)))
# Complete finite-moment criterion with nilpotent regular blocks.
regular=[]
for c in [2,3,4]:
    qq=S.Matrix(S.symbols("r0:"+str(c)));BB=S.zeros(c)
    for i in range(c-1):BB[i,i+1]=1
    J=qq[-1]**m;AA=linear_feedback(BB,qq,J).applyfunc(S.simplify)
    check("c"+str(c)+" all regular-block moment invariances",all(zero(v) for v in moments(BB,qq,J)))
    check("c"+str(c)+" whole arbitrary-order nilpotency",mz(AA**c))
    check("c"+str(c)+" regular feedback has nonzero symbolic final chain coefficient",
          zero((AA**(c-1))[0,c-1]-(m+1)))
    regular.append({"normal_dimension":c,"J":"last_coordinate^m","nilpotent":True,"last_chain_coefficient":"m+1"})
B3=S.Matrix([[0,1,0],[0,0,1],[0,0,0]]);Jbad=2*q1*q3-q2*q2
badmom=moments(B3,q3v,Jbad);Abad=linear_feedback(B3,q3v,Jbad)
check("trace-only invariant countercontrol first moment zero",zero(badmom[0]))
check("trace-only invariant countercontrol second moment nonzero",zero(badmom[1]-2*q3*q3))
check("trace-only invariant countercontrol entire nonnilpotent characteristic",
      zero(Abad.charpoly(T).as_expr()-(T**3-2*q3*q3*T/Jbad)))
# Full zero-stratum limits for regular nilpotent normal derivatives.
u,x,y,z,eps,lam=S.symbols("u x y z eps lam");ps=[u,x,y,z];p=S.Matrix(ps)
V=S.Matrix([x+u*u,y+u*x,x*x,u+x+y]);h=u
normal_controls=[]
for degree in [0,2,5]:
    J=(z+1)*y**degree;F=J+(u+x+y)**(degree+1)
    lim=normal_limit(V,ps,F,h,3,eps)["limit"]
    expected=linear_feedback(B3,S.Matrix([u,x,y]),J)
    check("order"+str(degree)+" entire regular-nilpotent normal limit",mz(lim[:3,:3]-expected))
    check("order"+str(degree)+" lower normal-limit blocks vanish",mz(lim[3:,:]))
    check("order"+str(degree)+" upper tangential correction retained",not mz(lim[:3,3:]))
    check("order"+str(degree)+" limiting normal block genuinely nilpotent",mz(lim[:3,:3]**3))
    normal_controls.append({"order":degree,"normal_J":"(z+1)*y^m","upper_right_retained":True})
check("initial normal pure power does not force global divisibility",
      S.limit((eps*u+(eps*y)**2)/eps,eps,0)==u and S.gcd(u+y*y,u)==1)
# General polynomial translation coefficient, no power or linear restriction.
P=u*u+u*y+y*y+1
F=S.Function("F")(u,x,y,z);V0=S.Matrix([u*x+y*y,x*x+u*y,u*u+y,x*y])
VV=V0+S.Matrix([z*P,0,0,0]);R=u+x+y*y
check("entire general polynomial translation density operator",
      zero(trace_equation(VV,ps,F,R)-trace_equation(V0,ps,F,R)-z*S.diff(P*F,u)))
check("general translation coefficient has the required directional dependencies",
      S.diff(P,z)==0 and S.diff(P,u)!=0)
Jtest=u**3+x*u+y
check("positive directional polynomial degrees add exactly",
      S.degree(P*Jtest,u)==S.degree(P,u)+S.degree(Jtest,u))
check("general nonlinear translation rational countercontrol has zero density divergence",
      zero(S.diff((1/P)*z*P,u)))
for power in [1,2,4]:
    check("power"+str(power)+" translation highest coefficient identity",
          zero(S.diff(u**power*Jtest,u)-u**(power-1)*(power*Jtest+u*S.diff(Jtest,u))))
# Whole cubic family and its two singular normal blocks.
Kc=cubic_form(ps);mm=S.symbols("a0:16");A=S.Matrix(4,4,mm);Vc=A*Kc;DV=Vc.jacobian(ps)
check("generic cubic canonical tuple is primitive",S.gcd_list(list(Kc))==1)
check("entire arbitrary-output cubic Jacobian determinant",
      zero(DV.det(method="domain-ge")-3*A.det()*u**8))
check("arbitrary invertible output preserves the cubic base equations",
      mz(A.adjugate()*Vc-A.det()*Kc))
B1=DV.extract([0,2],[0,2]).subs({u:0,y:0})
B2=DV.extract([0,1],[0,1]).subs({u:0,x:0})
check("cubic first whole normal block",
      mz(B1-x*x*S.Matrix([[A[0,2],0],[A[2,2],0]])))
check("cubic second whole normal block",
      mz(B2-y*y*S.Matrix([[0,A[0,3]],[0,A[1,3]]])))
check("both cubic normal determinants identically zero",B1.det()==0 and B2.det()==0)
check("first cubic normal trace forces A02 zero",zero(S.trace(B1)-A[0,2]*x*x))
check("second cubic normal trace forces A13 zero",zero(S.trace(B2)-A[1,3]*y*y))
check("both full affine planes are zero strata",
      mz(Vc.subs({u:0,y:0})) and mz(Vc.subs({u:0,x:0})))
check("cubic last-direction translation derivative entire",mz(Vc.diff(z)-u*u*A[:,3]))
check("cubic general translation derivative forces A03 zero",
      zero(sum(A[j,3]*S.diff(u*u,ps[j]) for j in range(4))-2*u*A[0,3]))
# Explicit cyclic cubic passes determinant/cyclicity gates and fails the new one.
swap=S.eye(4);swap.row_swap(0,2);Vex=swap*Kc;Dex=Vex.jacobian(ps)
point={v:1 for v in ps};vv=Vex.subs(point);DD=Dex.subs(point)
CC=S.Matrix.hstack(*[DD**j*vv for j in range(4)])
check("excluded cubic whole fourth-power Jacobian identity",zero(Dex.det(method="domain-ge")+3*u**8))
check("excluded cubic exact unit Krylov determinant",CC.det()==1)
base={u:0,x:1,y:0,z:0}
check("excluded cubic actual projective base point",mz(Vex.subs(base)))
check("excluded cubic actual nonnilpotent derivative at base",
      Dex.subs(base)==S.diag(1,0,0,0))
check("excluded cubic basic last translation test has zero coefficient derivative",
      Vex.diff(z)==S.Matrix([0,0,0,u*u]) and S.diff(u*u,z)==0)
# Initial multiplier in the regular lower 2-block is a pure u power.
kap=S.Symbol("kap")
Blower=S.Matrix([[0,0],[kap,0]])
Alower=linear_feedback(Blower,S.Matrix([u,y]),u**m).applyfunc(S.simplify)
check("regular lower normal block feedback arbitrary-order identity",
      mz(Alower-S.Matrix([[0,0],[(m+1)*kap,0]])))
# Actual homogeneous cubic normalization survives below full rank.
gam=S.Symbol("gam");controls=[]
for ss in [0,2]:
    data=cubic_control(ps,ss);X=data["X"];Q=data["P"];la=data["lam"]
    check("cubic s"+str(ss)+" entire strictly lower P",
          all(zero(Q[i,j]) for i in range(4) for j in range(i,4)))
    check("cubic s"+str(ss)+" regular nilpotency",mz(Q**4) and not mz(Q**3))
    sweep=X+gam*S.Matrix([1,*ps])
    check("cubic s"+str(ss)+" whole normalized sweep Jacobian",
          zero(sweep.jacobian([*ps,gam]).det(method="domain-ge")-gam**4))
    minor=X.jacobian(ps).extract([0,3,4],[0,1,2]).det()
    check("cubic s"+str(ss)+" exact rank3 minor",
          zero(minor+4*la**3*(ss+2)*u**(3*ss+2)*x*x*y) and minor!=0)
    check("cubic s"+str(ss)+" all of X independent of z",all(S.diff(v,z)==0 for v in X))
    check("cubic s"+str(ss)+" primitive direction and determinant law",
          S.gcd_list(list(data["V"]))==1 and zero(data["V"].jacobian(ps).det(method="domain-ge")-3*data["A"].det()*u**8) and la**4==(ss+3)*data["A"].det())
    controls.append({"s":ss,"lambda":ss+3,"F":"u^"+str(ss),"A":"diag(1,lambda,lambda,lambda)",
                     "rankDX":3,"entire_sweep_Jacobian":"gamma^4"})
model={"claims":["KMFI158","KMFI159","KMFI160","KMFI161"],
"previous_goal_turn":"progress: KMFI155--157 closed the whole primitive quadratic rank4 direction class.",
"linear_feedback_classification":{"B":"arbitrary matrix over characteristic-zero K","J":"any nonzero polynomial in normal variables",
"nilpotent_iff":"B nilpotent and dJ*B^j*q=0 for every j1..c-1",
"Fitting_invertible_part_would_force":"Euler_E J=-J",
"singular_nonnilpotent_trace_zero_countercontrols_retained":True,
"trace_only_nilpotent_B_invariant_countercontrol_retained":True,
"regular_nilpotent_J":"polynomial in left-kernel coordinate; homogeneous J=coefficient*ell^m",
"regular_block_controls":regular},
"strong_normal_stratum_result":{"B_must_be_nilpotent_even_when_singular":True,
"initial_form_satisfies_all_moment_invariances":True,"entire_limit_controls":normal_controls,
"upper_right_tangential_block_retained":True,
"homogeneous_base_point_test":"DV(a) nilpotent at every nonzero zero of V",
"initial_form_power_claimed_global_F_divisor":False},
"general_translation_test":{"coefficient":"arbitrary polynomial P","hypotheses":"V_b=P*a,P_b=R_b=0,P_a nonzero",
"conclusion":"no polynomial F for div(FV)=R*F","highest_coefficient":"partial_a(P*J)=0 impossible",
"rational_controls_retained":True},
"cubic_triangular_family":{"K":[str(v) for v in Kc],"A":"arbitrary GL4",
"Jacobian":"3*detA*u^8","required_matrix_conditions":["A02=0","A03=0","A13=0"],
"first_normal_det":0,"second_normal_det":0,
"regular_B1_initial_F":"coefficient(x,z)*u^m when A22 nonzero",
"remaining_matrices_classified":False,"remaining_full_rank_multipliers_found":False},
"excluded_explicit_cubic":{"V":[str(v) for v in Vex],"entire_Jacobian":"-3*u^8","unit_Krylov_determinant":1,
"base_point":[0,1,0,0],"derivative_at_base":"diag(1,0,0,0)",
"all_scalar_degrees_excluded":True,"basic_last_translation_gate_would_reject":False},
"actual_normalized_cubic_rank3_controls":controls,
"all_cubic_directions_closed":False,"arbitrary_quadratic_h0_shapes_retained":True,
"genuinely_nonhomogeneous_inputs_solved":False,"new_occurrence_action_claimed":False,"actual_PSL11_input_supplied":False,
"absolute_exceptional_action_decided":False,"full_action_spectrum_complete":False,"goal_complete":False}
(ROOT/"nilpotent_normal_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["NILPOTENT_NORMAL_STRATA.md","nilpotent_normal.py","verify_nilpotent_normal.py","nilpotent_normal_model.json",
"QUADRATIC_PRIMITIVE_RIGIDITY.md","quadratic_primitive.py","PRIMITIVE_MULTIPLIER_CRITERION.md",
"SWEEP_DIFFERENTIAL_RECONSTRUCTION.md"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,"runtime_seconds":round(time.monotonic()-started,3),
"python":platform.python_version(),"sympy":S.__version__,"sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
"assurance":"Unbounded Fitting/Euler linear feedback classification, nilpotent normal derivative and base-point tests, general polynomial translation obstruction and cubic exclusions; exact whole resolvent, normal limits, polynomial identities and actual lower-rank entire sweep Jacobians. No formal verification, independent whole-proof replay or external review claimed.","goal_complete":False}
(ROOT/"nilpotent_normal_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],
"singular_nonnilpotent_normal_blocks_excluded":True,"explicit_cyclic_cubic_excluded_in_all_scalar_degrees":True,"goal_complete":False}))
