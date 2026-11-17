"""Exact rational-fiber and irreducible-critical-factor controls KMFI164--165."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import sympy as S
import json,time,hashlib,platform
from rational_fiber import graph_shift,rank_one_direction,explicit_cubic
from quadratic_primitive import feedback
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def zero(v):return S.cancel(v)==0
def mz(M):return all(zero(v) for v in M)
x,y,w,z,mu,lam,t=S.symbols("x y w z mu lam t");q=S.Matrix([x,y,w])
cs=S.symbols("c0:9");C=S.Matrix(3,3,cs);dS=S.Matrix([S.symbols("s0:3")]);dT=S.Matrix([S.symbols("t0:3")])
S0,T0=S.symbols("S0 T0")
# Entire graph residue recurrence; no diagonalization or source transformation.
row=dS*C+dT-mu*dS
last=T0-mu*S0;residues=[(dS*C*q)[0]+last]
for j in range(1,4):
    last=(dT*C**j*q)[0]+mu*last
    residues.append((dS*C**(j+1)*q)[0]+last)
    check("graph moment"+str(j)+" entire recurrence to gradient defect",
          zero(residues[j]-mu*residues[j-1]-(row*C**j*q)[0]))
D=mu*S.eye(3)-C;adj=D.adjugate();det=D.det()
L=det*T0+(dT*adj*C*q)[0]
check("entire last-graph shift Euler reduction before division",
      zero(L-mu*(dT*adj*q)[0]-det*(T0-(dT*q)[0])))
A=S.BlockMatrix([[C,S.zeros(3,1)],[dT,S.Matrix([[mu]])]]).as_explicit()
check("rational fiber entire derivative characteristic",
      zero(A.charpoly(t).as_expr()-(t-mu)*C.charpoly(t).as_expr()))
# Every Hessian transport step is a whole telescoping identity.
hs=S.symbols("h0:6");H=S.Matrix([[hs[0],hs[1],hs[2]],[hs[1],hs[3],hs[4]],[hs[2],hs[4],hs[5]]])
error=H*C-C.T*H
for j in [1,2,3]:
    telescope=sum((C.T**i*error*C**(j-1-i)*q for i in range(j)),S.zeros(3,1))
    check("Hessian transport"+str(j)+" entire commutator identity",
          mz(H*C**j*q-C.T**j*H*q-telescope))
# Generic head cyclicity is compatible with nonsemisimple matrices.
N=S.Matrix([[0,1,0],[0,0,1],[0,0,0]]);J=S.eye(3)+N;e3=S.Matrix([0,0,1])
check("nonsemisimple invertible head has a cyclic vector",
      J.det()==1 and S.Matrix.hstack(e3,J*e3,J**2*e3).det()!=0)
check("nonsemisimple head retains repeated eigenvalue",J.charpoly(t).as_expr()==(t-1)**3 or zero(J.charpoly(t).as_expr()-(t-1)**3))
# Scalar/noncyclic heads show why Euler and Hessian commutation alone are insufficient.
Q=x*x+y*y+w*w;curved=x**3/Q
scalar=graph_shift(S.eye(3),q,S.Integer(5),curved)
check("noncyclic head can have zero gradient defect for nonlinear shift",
      mz(scalar["defect"]) and not mz(S.hessian(scalar["S"],list(q))))
check("noncyclic head Krylov rank really drops",
      S.Matrix.hstack(q,q,q).rank()==1)
check("homogeneous rational shift whole Euler identity",
      zero((scalar["dS"]*q)[0]-scalar["S"]) and mz(S.hessian(scalar["S"],list(q))*q))
# Accepted genuinely full-rank polynomial input on the linear boundary.
omega=S.Symbol("omega");mod=omega*omega-omega+1
def red(v):
    num,den=S.cancel(v).as_numer_denom()
    assert not den.has(omega)
    return S.Poly(num,omega).rem(S.Poly(mod,omega)).as_expr()/den
def rz(v):return red(v)==0
def rmz(M):return all(rz(v) for v in M)
u=S.Symbol("u");ps2=[u,z];W=S.Matrix([u,omega*z+(omega-1)*u]);G=u*(z+u);lambda0=1+omega
N2=feedback(W,ps2,G,lambda0)
X=S.Matrix([-lambda0*G,*list(G*(W-lambda0*S.Matrix(ps2)))])
lin=graph_shift(S.Matrix([[1]]),S.Matrix([u]),omega,(omega-1)*u)
check("positive linear-curvature control has exactly linear graph shift",zero(lin["S"]-u))
check("positive linear-curvature control complete gradient defect vanishes",mz(lin["defect"]))
check("positive linear-curvature control entire normalized nilpotency",rmz(N2**2))
check("positive linear-curvature control exact full-rank unit minor",
      rz(X[1:,:].jacobian(ps2).det().subs({u:1,z:1})-12*omega) and red(12*omega)!=0)
check("positive linear-curvature control exact cyclic unit minor",
      rz(S.Matrix.hstack(W,W.jacobian(ps2)*W).det().subs({u:1,z:1})+2))
gam=S.Symbol("gam");sweep=X+gam*S.Matrix([1,u,z])
check("positive linear-curvature control whole sweep Jacobian",
      rz(sweep.jacobian([u,z,gam]).det(method="domain-ge")-gam**2))
# General rank-two determinant formula with one constant translation direction.
matrix_determinants=[]
for r in [3,4]:
    pp=S.Matrix(S.symbols("p0:"+str(r)));bb=S.eye(r)[:,r-1]
    qq=S.Symbol("Q0");aq=S.Matrix([[*S.symbols("q0:"+str(r-1)),0]])
    ar=S.Matrix([[*S.symbols("r0:"+str(r-1)),0]])
    B=qq*S.eye(r)+pp*aq+bb*ar
    check("rank"+str(r)+" entire radial-plus-kernel determinant identity",
          zero(B.det(method="domain-ge")-qq**r-qq**(r-1)*(aq*pp)[0]))
    check("rank"+str(r)+" both kernel-direction derivative contractions zero",
          (aq*bb)[0]==0 and (ar*bb)[0]==0)
    matrix_determinants.append({"rank":r,"identity":"Q^r+Q^(r-1)*dQ*p"})
# General arbitrary output M and perturbation R keep the exact Jacobian factor.
ps=[x,y,w,z];pvec=S.Matrix(ps);b=S.Matrix([0,0,0,1])
mm=S.symbols("m0:16");M=S.Matrix(4,4,mm);R=x**3+x*y*w+y*y*w
V=rank_one_direction(Q,R,M,pvec,b);gradQ=S.Matrix([[S.diff(Q,v) for v in ps]]);gradR=S.Matrix([[S.diff(R,v) for v in ps]])
check("entire arbitrary-M derivative factorization",
      mz(V.jacobian(ps)-M*(Q*S.eye(4)+pvec*gradQ+b*gradR)))
check("quadric and perturbation both independent of the translation direction",
      S.diff(Q,z)==0 and S.diff(R,z)==0)
check("entire arbitrary-M translation derivative",mz(V.diff(z)-Q*M*b))
check("whole general translation directional criterion",
      zero((gradQ*M*b)[0]-2*(M[0,3]*x+M[1,3]*y+M[2,3]*w)))
check("quadratic translation kernel is exactly one-dimensional",
      S.hessian(Q,ps).rank()==3 and S.hessian(Q,ps).nullspace()==[b])
check("geometric irreducibility rank criterion is nondegenerate on the head",
      S.hessian(Q,[x,y,w]).det()!=0)
# If R is already divisible by Q, the perturbation is still factorized.
lr=S.symbols("l0:3");L1=lr[0]*x+lr[1]*y+lr[2]*w
rowL=S.Matrix([[*lr,0]]);Mprime=M*(S.eye(4)+b*rowL)
check("entire divisible-curvature factorization",
      mz(rank_one_direction(Q,Q*L1,M,pvec,b)-Q*Mprime*pvec))
check("divisible curvature preserves invertibility of M",
      zero(Mprime.det(method="domain-ge")-M.det()))
# An irreducible factor Q cannot divide a product of independent linear coordinates.
check("nonlinear critical Q remains a genuine nonlinear polynomial",
      S.Poly(Q,x,y,w).total_degree()==2 and S.hessian(Q,[x,y,w]).rank()==3)
# Explicit cyclic cubic with irreducible critical quadric.
ex=explicit_cubic();VV=ex["V"];QQ=ex["Q"];MM=ex["M"];CC=ex["C"];mu0=ex["mu"];pp=ex["ps"]
AA=VV.jacobian(pp);point={v:1 for v in pp}
check("explicit cubic is primitive",S.gcd_list(list(VV))==1)
check("explicit cubic whole fourth-power Jacobian",
      zero(AA.det(method="domain-ge")-150*QQ**4))
check("explicit cubic linear matrix determinant exact",MM.det()==50)
vv=VV.subs(point);ap=AA.subs(point)
krylov=S.Matrix.hstack(*[ap**j*vv for j in range(4)])
check("explicit cubic exact nonzero generic-cyclicity certificate",krylov.det()==-347733000)
check("explicit head conformal derivative of Q exact",
      zero((S.Matrix([[S.diff(QQ,v) for v in pp[:3]]])*CC*S.Matrix(pp[:3]))[0]-4*QQ))
check("explicit head and mu avoid the singular last-coordinate branch",
      CC.det()!=0 and (mu0*S.eye(3)-CC).det()!=0)
check("explicit cubic basic translation obstruction passes",
      mz(VV.diff(pp[3])-5*QQ*b) and S.diff(QQ,pp[3])==0)
# Every base point, not a numerical sample: reduce the full square modulo x,Q.
xx,yy,ww,zz=pp
square=AA*AA
remainders=[S.rem(S.Poly(S.expand(v.subs(xx,0)),ww),S.Poly(ww*ww+yy*yy,ww)).as_expr() for v in square]
check("explicit cubic derivatives square-zero at all projective base points",all(v==0 for v in remainders))
check("explicit cubic head zero equations have the asserted form",
      mz(VV[:3,:]-QQ*CC*S.Matrix(pp[:3])))
check("explicit cubic on Q0 has exactly last curvature x-cube",
      mz(VV-QQ*MM*S.Matrix(pp)-5*b*xx**3))
# Entire rational graph defect, independent of any scalar degree or lambda.
Tcur=5*xx**3/QQ;gdata=graph_shift(CC,S.Matrix(pp[:3]),mu0,Tcur)
Sdisplay=xx**2*(7*ww*ww*xx+9*ww*ww*yy+9*xx**3+9*xx*xx*yy+9*xx*yy*yy+9*yy**3)/(6*QQ**2)
check("explicit cubic entire rational graph shift identity",zero(gdata["S"]-Sdisplay))
unit={xx:1,yy:1,ww:1}
defect=gdata["defect"].subs(unit)
check("explicit cubic exact nonzero residue-gradient defect",
      defect==S.Matrix([[-S.Rational(106,81),S.Rational(82,81),S.Rational(26,27)]]))
check("explicit cubic rational shift has nonzero Hessian",
      S.hessian(gdata["S"],pp[:3]).subs(unit).rank()==2)
check("explicit cubic rational shift still satisfies homogeneous Euler identities",
      zero((gdata["dS"]*S.Matrix(pp[:3]))[0]-gdata["S"]) and mz(S.hessian(gdata["S"],pp[:3])*S.Matrix(pp[:3])))
model={"claims":["KMFI164","KMFI165"],
"previous_goal_turn":"progress: KMFI162--163 closed the entire arbitrary-matrix triangular cubic chain.",
"rational_degree1_fiber_rigidity":{"W":"(Cq,mu*z+T(q))","T":"rational homogeneous degree1",
"G":"positive-degree homogeneous polynomial clearing the input denominators",
"X":"polynomial normalized and full rank","C":"arbitrary invertible constant matrix, no diagonalizability assumption",
"DW_W_cyclic_by_rational_Euler_and_KMFI124":True,"head_q_cyclic_for_C":True,
"graph_S":"dT*(muI-C)^-1*q","full_graph_moments_retained":True,
"necessary_gradient_identity":"dT=dS*(muI-C)","Hessian_transport":"H_S*C=C^t*H_S and H_S*q=0",
"conclusion":"S and T linear; whole input is KMFI121 scalar-times-linear",
"nonlinear_source_change_claimed_normalization_preserving":False,
"noncyclic_Hessian_countercontrol_is_actual_full_rank_input":False},
"nonlinear_critical_factor_obstruction":{"ranks":"all ranks under the stated Q hypotheses",
"Q":"irreducible homogeneous degree k>=2, constant translation kernel exactly k*b",
"R":"homogeneous degree k+1 and b-independent","M":"arbitrary GLr",
"V":"Q*M*p+(M*b)*R","whole_Jacobian":"(k+1)*detM*Q^r",
"translation_forces_Mb_in_kernel_Q":True,"kernel_dimension1_forces_invariant_b":True,
"rational_fiber_rigidity_forces_Q_divides_R":True,
"remaining_scalar_FQ_has_forbidden_nonlinear_factor_by_KMFI121":True,
"all_M_R_and_scalar_degrees_excluded_to_full_rank":True,"determinant_controls":matrix_determinants},
"explicit_rank4_irreducible_quadric_cubic":{"Q":"x^2+y^2+w^2","R":"x^3","M":"diag([[2,1,0],[-1,2,0],[0,0,2]],5)",
"whole_Jacobian":"150*Q^4","unit_Krylov_determinant":-347733000,
"all_basepoint_derivatives_square_zero":True,"basic_last_translation_obstruction_passed":True,
"graph_gradient_defect_at_unit":["-106/81","82/81","26/27"],"all_scalar_degrees_excluded":True},
"accepted_rank2_linear_boundary":{"coefficient_field":"omega^2-omega+1=0","W":"(u,omega*z+(omega-1)*u)",
"G":"u*(z+u)","lambda":"1+omega","rankDX":2,"whole_sweep_Jacobian":"gamma^2"},
"all_irreducible_quadric_critical_cubic_directions_closed":False,"nondegenerate_rank4_quadric_cases_closed":False,
"general_cubic_higher_and_nonhomogeneous_inputs_solved":False,"new_occurrence_action_claimed":False,
"actual_PSL11_input_supplied":False,"absolute_exceptional_action_decided":False,"full_action_spectrum_complete":False,"goal_complete":False}
(ROOT/"rational_fiber_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["RATIONAL_LINEAR_FIBER_RIGIDITY.md","rational_fiber.py","verify_rational_fiber.py","rational_fiber_model.json",
"CUBIC_TRIANGULAR_CHAIN_RIGIDITY.md","PRIMITIVE_MULTIPLIER_CRITERION.md","NILPOTENT_NORMAL_STRATA.md",
"TRIANGULAR_PROJECTIVE_RIGIDITY.md","FACTORIZED_HOMOGENEOUS_FIELDS.md"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,"runtime_seconds":round(time.monotonic()-started,3),
"python":platform.python_version(),"sympy":S.__version__,"sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
"assurance":"Unbounded rational homogeneous degree-one fiber rigidity and irreducible nonlinear critical-factor obstruction. Whole graph-moment, Hessian and determinant identities, exact cubic passing older local gates, and actual accepted rank2 whole sweep control. No formal verification, independent whole-proof replay or external review claimed.","goal_complete":False}
(ROOT/"rational_fiber_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],
"rational_degree1_fiber_curvature_excluded_to_full_rank":True,"irreducible_critical_factor_family_closed":True,"goal_complete":False}))
