"""Exact nonlinear-head, rational-flag and multi-fiber identities KMFI166--168."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import sympy as S
import json,time,hashlib,platform
from rational_flag import nonlinear_fiber,multi_direction,split_cubic
from quadratic_primitive import feedback
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def zero(v):return S.cancel(v)==0
def mz(M):return all(zero(v) for v in M)
x,y,w,z,mu,lam=S.symbols("x y w z mu lam");q=S.Matrix([x,y,w]);ps=[x,y,w,z];p=S.Matrix(ps)
f=x*x/y;v=S.Matrix([x,y,w+f]);T=w-f
data=nonlinear_fiber(v,q,S.Integer(2),T);B=data["B"]
check("nonlinear head derivative determinant is nonzero constant",zero(B.det()-1))
check("nonlinear head and rational fiber are homogeneous degree1",
      mz(B*q-v) and zero(sum(q[i]*S.diff(T,q[i]) for i in range(3))-T))
check("nonlinear T can have a linear graph shift",zero(data["S"]-w))
check("linear graph shift has zero full gradient defect",mz(data["defect"]))
check("fiber term itself remains nonlinear rational",not mz(S.hessian(T,list(q))))
W=S.Matrix([*list(v),2*z+T]);C=S.eye(4);C[3,2]=1
Wnew=C*W.subs(z,z-w)
check("entire allowed constant shear decouples nonlinear rational fiber",
      mz(Wnew-S.Matrix([*list(v),2*z])))
check("decoupling shear is genuinely invertible constant linear",C.det()==1)
# Nonconstant-head Hessian terms must not be omitted.
HS=S.hessian(data["S"],list(q));HT=S.hessian(T,list(q))
extra=sum((data["dS"][i]*S.hessian(v[i],list(q)) for i in range(3)),S.zeros(3))
check("entire nonlinear-head Hessian identity with all extra terms",
      mz(HT-2*HS+B.T*HS+extra))
check("nonlinear-head extra symmetric Hessian term is genuinely nonzero",not mz(extra))
check("omitting the nonlinear-head term would give a false identity",
      not mz(HT-2*HS+B.T*HS))
# An arbitrary nonzero defect has exactly the derivative predicted by the formula.
Sgen=x**3/(x*x+y*y+w*w)+w;Tgen=w-f
ds=S.Matrix([[S.diff(Sgen,t) for t in q]]);dt=S.Matrix([[S.diff(Tgen,t) for t in q]])
defect=ds*B+dt-2*ds
extra=sum((ds[i]*S.hessian(v[i],list(q)) for i in range(3)),S.zeros(3))
check("whole Hessian identity also retains derivative of a nonzero gradient defect",
      mz(S.hessian(Tgen,list(q))-2*S.hessian(Sgen,list(q))+B.T*S.hessian(Sgen,list(q))+extra-defect.jacobian(list(q))))
# Cyclicity remains essential: nonconstant head with a nonlinear graph satisfying residues.
nc=nonlinear_fiber(v,q,S.Integer(2),f)
check("noncyclic nonlinear head permits nonlinear residue shift",
      zero(nc["S"]-f) and mz(nc["defect"]) and not mz(S.hessian(nc["S"],list(q))))
Khead=S.Matrix.hstack(q,B*q,B**2*q)
check("noncyclic head really has zero Krylov determinant",zero(Khead.det()))
Wnc=S.Matrix([*list(v),2*z+f]);Anc=Wnc.jacobian(ps);ell=S.Matrix([[*list(nc["dS"]),1]])
for j in range(4):
    residue=(ell*Anc**j*Wnc)[0].subs(z,-nc["S"])
    check("noncyclic boundary entire graph moment"+str(j),zero(residue))
# Euler and determinant identities apply to the full rational map.
A=W.jacobian(ps);G=x*y*w*z
check("full rational degree1 direction satisfies entire Euler identity",mz(A*p-W))
check("full rational degree1 derivative is invertible",zero(A.det()-2))
left=(G*W).jacobian(ps);right=A*(G*S.eye(4)+p*S.Matrix([[S.diff(G,t) for t in ps]]))
check("whole homogeneous scalar derivative factorization",mz(left-right))
check("whole scalar Euler determinant factor",
      zero((G*S.eye(4)+p*S.Matrix([[S.diff(G,t) for t in ps]])).det(method="domain-ge")-5*G**4))
# Flag spectrum and descending-shear parameter identities retain repeated branches.
a1,a2,a3=S.symbols("a1 a2 a3");entries=S.symbols("b0:3")
Bh=S.Matrix([[a1,0,0],[entries[0],a2,0],[entries[1],entries[2],a3]])
check("entire flag last-leading determinant",
      zero((mu*S.eye(3)-Bh).det()-(mu-a1)*(mu-a2)*(mu-a3)))
for ai in [a1,a2,a3]:
    check("flag repeated-eigenvalue branch explicitly has zero leading coefficient "+str(ai),
          zero((mu*S.eye(3)-Bh).det().subs(mu,ai)))
L0,L1,L2,L3,d1,d2=S.symbols("L0 L1 L2 L3 d1 d2")
s1,s2=S.symbols("s1 s2");xx,yy=S.symbols("xx yy")
head=S.Matrix([xx,2*yy,d1*s1,d2*s2]);row=S.Matrix([[L0,L1,L2,L3]])
shift=(row*S.Matrix([xx,yy,s1,s2]))[0];cob=mu*shift-(row*head)[0]
check("first already-decoupled coefficient retained in descending shear",
      zero(S.diff(cob,s1)-(mu-d1)*L2))
check("second already-decoupled coefficient retained in descending shear",
      zero(S.diff(cob,s2)-(mu-d2)*L3))
# Multi-fiber determinant: all perturbations and kernel directions retained.
det_controls=[]
for r,c in [(4,2),(5,2),(5,3)]:
    pp=S.Matrix(S.symbols("p0:"+str(r)));Q0=S.Symbol("Q0")
    B0=S.eye(r)[:,r-c:r];dq=S.Matrix([[*S.symbols("q0:"+str(r-c)),*([0]*c)]])
    dr=S.Matrix(c,r,lambda i,j:S.Symbol("r"+str(i)+"_"+str(j)) if j<r-c else 0)
    D=Q0*S.eye(r)+pp*dq+B0*dr
    check("rank"+str(r)+" fibers"+str(c)+" whole determinant independent of perturbations",
          zero(D.det(method="domain-ge")-Q0**r-Q0**(r-1)*(dq*pp)[0]))
    check("rank"+str(r)+" fibers"+str(c)+" all kernel contractions zero",
          mz(dq*B0) and mz(dr*B0))
    det_controls.append({"rank":r,"fibers":c,"identity":"Q^r+Q^(r-1)*dQ*p"})
# Arbitrary M and both independent perturbations in the split-quadric family.
u=S.Symbol("u");ps4=[u,x,y,z];p4=S.Matrix(ps4);Q=u*x
B0=S.eye(4)[:,2:4];R=S.Matrix([u**3,x**3]);mm=S.symbols("m0:16");M=S.Matrix(4,4,mm)
V=multi_direction(Q,R,M,p4,B0);dq=S.Matrix([[S.diff(Q,t) for t in ps4]])
dr=R.jacobian(ps4)
check("whole arbitrary-M multi-fiber derivative factorization",
      mz(V.jacobian(ps4)-M*(Q*S.eye(4)+p4*dq+B0*dr)))
check("split critical-factor constant translation kernel is exact",
      S.hessian(Q,ps4).nullspace()==[S.eye(4)[:,2],S.eye(4)[:,3]])
for i in [0,1]:
    bvec=B0[:,i]
    check("kernel direction"+str(i)+" entire translation derivative",mz(V.jacobian(ps4)*bvec-Q*M*bvec))
    check("kernel direction"+str(i)+" target derivative detects invariant subspace",
          zero((dq*M*bvec)[0]-(M[0,i+2]*x+M[1,i+2]*u)))
check("both perturbations invariant under every kernel direction",mz(dr*B0))
# Divisible perturbations leave only the scalar-times-linear map, without singular update.
ls=S.symbols("l0:4");L=S.Matrix([[ls[0],ls[1],0,0],[ls[2],ls[3],0,0]])
U=S.eye(4)+B0*L
check("whole kernel-linear update squares to zero",mz((B0*L)**2))
check("whole divisible update has determinant1 and exact polynomial inverse",
      U.det()==1 and mz(U*(S.eye(4)-B0*L)-S.eye(4)))
check("entire divisible multi-perturbation factorization",
      mz(multi_direction(Q,Q*L*p4,M,p4,B0)-Q*M*U*p4))
# Explicit primitive cubic over two distinct critical hyperplanes.
ex=split_cubic();pp=ex["ps"];VV=ex["V"];QQ=ex["Q"];AA=VV.jacobian(pp)
check("explicit split-quadric cubic is primitive",S.gcd_list(list(VV))==1)
check("explicit split-quadric entire fourth-power Jacobian",
      zero(AA.det(method="domain-ge")-72*QQ**4))
unit={t:1 for t in pp};av=AA.subs(unit);vv=VV.subs(unit)
check("explicit split-quadric cubic exact cyclic unit certificate",
      S.Matrix.hstack(*[av**j*vv for j in range(4)]).det()==1056)
check("explicit split-quadric all base derivatives vanish",mz(AA.subs({u:0,x:0})))
check("explicit split-quadric curvature cannot be divisible by Q",
      S.rem(S.Poly(u**3,x),S.Poly(u*x,x)).as_expr()!=0 and S.rem(S.Poly(x**3,u),S.Poly(u*x,u)).as_expr()!=0)
head=VV[:3,:]/QQ;Tlast=S.cancel(VV[3]/QQ-4*z)
last=nonlinear_fiber(head,S.Matrix(pp[:3]),S.Integer(4),Tlast)
check("explicit nonlinear-head fiber shift exact",zero(last["S"]-8*x*x/(3*u)))
check("explicit nonlinear-head gradient defect exact at unit",
      last["defect"].subs({u:1,x:1,y:1})==S.Matrix([[4,-S.Rational(8,3),0]]))
check("explicit nonlinear head itself passes cyclicity at unit",
      S.Matrix.hstack(S.Matrix(pp[:3]),last["B"]*S.Matrix(pp[:3]),last["B"]**2*S.Matrix(pp[:3])).subs({u:1,x:1,y:1}).det()==2)
# Actual accepted nonprimitive fully split boundary in rank3.
beta=S.Symbol("beta");mod=beta*beta+1
def red(v):
    num,den=S.cancel(v).as_numer_denom()
    assert not den.has(beta)
    return S.Poly(num,beta).rem(S.Poly(mod,beta)).as_expr()/den
def rz(v):return red(v)==0
def rmz(M):return all(rz(v) for v in M)
ps3=[x,y,z];p3=S.Matrix(ps3);weights=[S.Rational(1,2),(1+beta)/2,(1-beta)/2]
Mt=S.diag(*weights);G=x*y*z;Qsplit=x*y;Fs=z;Vsplit=Qsplit*Mt*p3
X=S.Matrix([-G,*list(G*(Mt-S.eye(3))*p3)]);N3=feedback(Mt*p3,ps3,G,S.Integer(1))
check("accepted split boundary whole nilpotent feedback",rmz(N3**3))
check("accepted split boundary actual nonzero full-rank minor",
      rz(X[1:,:].jacobian(ps3).det().subs({x:1,y:1,z:1})+1))
check("accepted split boundary exactly matches multi-fiber presentation",
      mz(Fs*Vsplit-G*Mt*p3))
check("accepted split boundary scalar divisibility exact",zero(G/Qsplit-Fs))
gam=S.Symbol("gam");sweep=X+gam*S.Matrix([1,x,y,z])
check("accepted split boundary entire normalized sweep Jacobian",
      rz(sweep.jacobian([gam,*ps3]).det(method="domain-ge")-gam**3))
check("accepted split boundary is nonprimitive as required",
      zero(S.gcd_list(list(Vsplit))-Qsplit))
model={"claims":["KMFI166","KMFI167","KMFI168"],
"previous_goal_turn":"progress: KMFI164--165 gave linear-head rational fiber rigidity and an irreducible one-kernel obstruction.",
"nonlinear_head_decoupling":{"W":"(v(q),mu*z+T(q))","v_and_T":"rational homogeneous degree1",
"polynomial_normalized_X_full_rank_required":True,"det_DW":"lambda^r/(degG+1)",
"head_B":"Dv, not assumed constant","whole_Hessian_identity_retains_sum_Si_Hvi":True,
"graph_shift_S":"constant linear Lq","fiber_term_T":"mu*Lq-Lv, not necessarily itself linear",
"allowed_equivalence":"simultaneous constant shear z_new=z+Lq",
"noncyclic_nonlinear_graph_control_is_actual_full_rank_input":False},
"rational_flag_rigidity":{"flag":"Wi=ai*pi+Ti(earlier variables)","diagonal":"constant, initially arbitrary",
"all_ranks_and_scalar_degrees":True,"repeated_eigenvalue_zero_branches_retained":True,
"descending_shifts_do_not_recouple_separated_fibers":True,
"conclusion":"every Ti linear, ai distinct/nonzero and not lambda, exactly KMFI121"},
"multi_fiber_classification":{"Q":"any nonconstant homogeneous polynomial","U":"full constant translation kernel of Q, positive dimension",
"R":"one homogeneous degree degQ+1 component per basis vector, all invariant under all U",
"V":"M*(Qp+B0R)","M":"arbitrary GLr",
"whole_Jacobian":"(degQ+1)*detM*Q^r","translation_forces_MU_equals_U":True,
"constant_block_triangularization_gives_rational_flag":True,
"full_rank_forces":"Q divides every Ri",
"remaining_data":"V=Q*Mprime*p, G=FQ among KMFI121 data with Q dividing G",
"any_nonlinear_irreducible_factor_in_Q":"no full-rank input",
"determinant_controls":det_controls,
"invariance_under_one_kernel_direction_alone_claimed_sufficient":False},
"explicit_split_quadric_cubic":{"V":[str(v) for v in VV],"Q":"u*x","whole_Jacobian":"72*u^4*x^4",
"unit_Krylov_determinant":1056,"all_basepoint_derivatives_zero":True,"entire_primitive_family_all_scalar_degrees_excluded":True,
"last_graph_gradient_defect_at_unit":["4","-8/3","0"]},
"accepted_fully_split_nonprimitive_boundary":{"rank":3,"beta":"beta^2+1=0","Q":"x*y","R":0,"F":"z","G":"x*y*z",
"whole_sweep_Jacobian":"gamma^3 with source order gamma,x,y,z","unit_full_rank_minor":-1,"new_occurrence_action":False},
"arbitrary_rational_degree1_directions_solved":False,"variable_diagonal_flags_solved":False,
"arbitrary_cubic_and_nondegenerate_critical_quadric_cases_solved":False,
"genuinely_nonhomogeneous_inputs_solved":False,"actual_PSL11_input_supplied":False,"absolute_exceptional_action_decided":False,
"full_action_spectrum_complete":False,"goal_complete":False}
(ROOT/"rational_flag_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["RATIONAL_FLAG_RIGIDITY.md","rational_flag.py","verify_rational_flag.py","rational_flag_model.json",
"RATIONAL_LINEAR_FIBER_RIGIDITY.md","rational_fiber.py","CUBIC_TRIANGULAR_CHAIN_RIGIDITY.md",
"PRIMITIVE_MULTIPLIER_CRITERION.md","NILPOTENT_NORMAL_STRATA.md","FACTORIZED_HOMOGENEOUS_FIELDS.md"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,"runtime_seconds":round(time.monotonic()-started,3),
"python":platform.python_version(),"sympy":S.__version__,"sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
"assurance":"Unbounded nonlinear-head fiber decoupling, rational constant-diagonal flag classification and multi-fiber perturbation reduction. Entire nonlinear Hessian, graph, flag and determinant identities, explicit primitive split-quadric cubic and actual normalized full-rank split nonprimitive boundary. No formal verification, independent whole-proof replay or external review claimed.","goal_complete":False}
(ROOT/"rational_flag_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],
"rational_constant_diagonal_flags_reduce_to_linear":True,"multi_fiber_primitive_curvature_excluded":True,"goal_complete":False}))
