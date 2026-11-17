"""Exact universal certificate and sharp controls for KMFI130--132."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import sympy as S
import json,time,hashlib,platform
from rank_one_projective import homogenize,companion,companion_frame
from primitive_multiplier import direction_shape
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
def zero(e):return S.cancel(e)==0
def mz(M):return all(zero(e) for e in M)
d=companion();u,x,y,z=d["ps"];lam,tau,sig,delta=d["parameters"];f=d["feedback"];D=d["D"];eta=d["eta"]
T=f["T"];A=f["A"].subs(u,1);vchart=d["V"].subs(u,1)
check("canonical whole chart matrix and characteristic",
 zero(A.charpoly(T).as_expr()-(T-2)*(T**3-tau*T**2+(sig-2*x)*T-delta)))
check("canonical whole projective determinant constant",zero(d["v"].jacobian([x,y,z]).det()-delta))
check("canonical full homogeneous determinant power",zero(d["V"].jacobian([u,x,y,z]).det(method="domain-ge")-2*delta*u**4))
check("canonical exact whole feedback solve",mz(eta*f["C"].subs(u,1)-f["moments"].subs(u,1)))
check("canonical whole feedback characteristic",
 zero((T-lam)**4-A.charpoly(T).as_expr()+(eta*(T*S.eye(4)-A).adjugate()*vchart)[0]))
check("D is a polynomial in source and every constant parameter",S.fraction(D)[1]==1)
check("universal denominator cannot specialize to zero",S.Poly(D,x,y,z).coeff_monomial(x**5*z)==-2)
curl=d["curl_polynomial"]
check("D squared times entire xy curl is polynomial",S.fraction(curl)[1]==1)
check("parameter independent curl coefficient is exactly8",S.Poly(curl,x,y,z).coeff_monomial(x*x*y*z**4)==8)
check("entire yz curl vanishes but does not suffice for closedness",zero(S.diff(eta[3],y)-S.diff(eta[2],z)))
# Special parameter checks supplement the universal coefficient, not a sample proof.
special=[]
for vals in [(0,0,1,1),(0,0,1,0),(2,-3,1,2),(3,4,-2,3)]:
    sub=dict(zip([tau,sig,delta,lam],vals));cp=S.Poly(curl.subs(sub),x,y,z)
    check("universal curl survives specialization "+str(vals),cp.coeff_monomial(x*x*y*z**4)==8)
    special.append({"tau":vals[0],"sigma":vals[1],"delta":vals[2],"lambda":vals[3],"constant_curl_coefficient":8})
# General linear mixing is reduced by the explicit companion frame.
L=S.Matrix([[1,1,0],[0,1,1],[1,0,0]])
M0=S.Matrix([[0,1,0],[-sig,tau,1],[delta,0,0]])
M=(L.inv()*M0*L).applyfunc(S.expand);cc=L.inv()*S.Matrix([0,2,0]);ell=S.Matrix([[1,0,0]])*L
frame=companion_frame(M,cc,ell)
check("mixed nontriangular determinant condition",zero((ell*M.inv()*cc)[0]))
check("mixed nontriangular transverse trace condition",zero((ell*cc)[0]) and frame["a"]==2)
check("three canonical forms are independent",frame["frame"].det()!=0)
check("whole canonical linear conjugacy",mz(frame["transformed_M"]-M0))
check("whole canonical nonlinear conjugacy",mz(frame["transformed_c"]-S.Matrix([0,2,0])))
check("common scaling normalizes nonzero quadratic coefficient",zero(2*(2*x*x)-(2*x)**2))
s,t0,r0=S.symbols("s t0 r0");pstar=S.Matrix([s,t0,r0]);c0=S.Matrix([0,1,0])
bvec=pstar-M0*pstar-c0*s*s
xx=S.Matrix([x,y,z]);affine_v=M0*xx+bvec+c0*x*x
translated=(affine_v.subs({x:x+s,y:y+t0,z:z+r0},simultaneous=True)-pstar).applyfunc(S.expand)
check("affine fixed point translation removes every constant",
 mz(translated-S.Matrix([y,z-(sig-2*s)*x+tau*y+x*x,delta*x])))
# Quadratic determinant gate: no diagonal constancy is assumed initially.
aa,bb,kk=S.symbols("a b k")
tri=S.Matrix([aa*x+x*x,bb*y+x*x,kk*z+x*y])
V=homogenize(tri,u,[x,y,z],2)
check("triangular nonconstant diagonal obeys whole determinant identity",
 zero(V.jacobian([u,x,y,z]).det(method="domain-ge")-2*u**4*tri.jacobian([x,y,z]).det().subs({x:x/u,y:y/u,z:z/u},simultaneous=True)))
check("a genuine nonconstant diagonal fails fourth-power gate",not direction_shape(V,[u,x,y,z])["power_compatible"])
tri2=S.Matrix([aa*x+1,bb*y+x*x,kk*z+x*y+y*y])
check("constant determinant triangular curvature survives shape only",
 zero(tri2.jacobian([x,y,z]).det()-aa*bb*kk))
# Derive the negative-cofactor leading equations from entire generic fields.
alpha,kappa=S.symbols("alpha kappa")
mm=list(S.symbols("m0:9"));Mgen=S.Matrix(3,3,mm);bv=S.Matrix(S.symbols("b0:3"))
coords=S.Matrix([x,y,z]);Alead=y*y+y*z+1;degree=2
poly=x**degree*Alead+x*(y+z)+1
Q=alpha*x*x+x*(y+z)+y*y
W=Mgen*coords+bv+S.Matrix([Q,0,0])-coords
residual=(S.Matrix([[S.diff(poly,v) for v in coords]])*W)[0]-(kappa-S.diff(Q,x))*poly
actual=S.Poly(S.expand(residual),x).coeff_monomial(x**(degree+1))
expected=alpha*(degree+2)*Alead+Mgen[1,0]*S.diff(Alead,y)+Mgen[2,0]*S.diff(Alead,z)
check("whole longitudinal quadratic field gives negative constant-derivative eigenvalue",zero(actual-expected))
Alead=x*x+x*z+1;poly=y**degree*Alead+y*(x+z)+1;Q=x*y+z*z
W=Mgen*coords+bv+S.Matrix([Q,0,0])-coords
residual=(S.Matrix([[S.diff(poly,v) for v in coords]])*W)[0]-(kappa-y)*poly
actual=S.Poly(S.expand(residual),y).coeff_monomial(y**(degree+1))
expected=(x+Mgen[0,1])*S.diff(Alead,x)+Mgen[2,1]*S.diff(Alead,z)+Alead
check("whole isotropic quadratic field gives negative Euler eigenvalue",zero(actual-expected))
for qdeg in [3,5]:
    Alead=y*y+y*z+1;poly=x**degree*Alead+x*(y+z)+1;Q=alpha*x**qdeg+x+1
    W=Mgen*coords+bv+S.Matrix([Q,0,0])-coords
    residual=(S.Matrix([[S.diff(poly,v) for v in coords]])*W)[0]-(kappa-S.diff(Q,x))*poly
    actual=S.Poly(S.expand(residual),x).coeff_monomial(x**(degree+qdeg-1))
    check("whole univariate degree"+str(qdeg)+" forces negative polynomial degree",zero(actual-alpha*(degree+qdeg)*Alead))
# Concrete constant-Jacobian rank-one projective maps cover every quadratic-form case.
cases=[]
for name,Q,cvec,Mcase,expectedrank in [
 ("rank2_nonisotropic",x*x+y*y,S.Matrix([1,0,1]),S.Matrix([[1,0,1],[0,1,0],[0,0,1]]),2),
 ("rank2_isotropic",x*y,S.Matrix([1,0,0]),S.Matrix([[0,0,1],[1,0,0],[0,1,0]]),2),
 ("rank2_radical",x*y,S.Matrix([0,0,1]),S.diag(1+S.I,1-S.I,2),2),
 ("rank1_invariant_plane",x*x,S.Matrix([0,0,1]),S.Matrix([[2,0,0],[0,3,1],[0,1,4]]),1),
 ("rank1_invariant_line",x*x,S.Matrix([0,0,1]),S.Matrix([[2,1,0],[1,3,0],[0,0,4]]),1),
 ("rank1_companion",x*x,S.Matrix([0,1,0]),S.Matrix([[0,1,0],[-2,3,1],[5,0,0]]),1)]:
    grad=S.Matrix([[S.diff(Q,v) for v in [x,y,z]]]);B=Mcase+cvec*grad
    check(name+" whole rank-one determinant gate",
      zero(B.det()-Mcase.det()) and mz(S.hessian(Q,[x,y,z])*Mcase.inv()*cvec))
    check(name+" quadratic rank exhaustive row",S.hessian(Q,[x,y,z]).rank()==expectedrank)
    cases.append({"name":name,"Q":str(Q),"rankQ":expectedrank,"detM":str(Mcase.det()),
       "DcQ":str((grad*cvec)[0]),"quadratic_on_c":str(Q.subs(dict(zip([x,y,z],cvec)),simultaneous=True))})
for name,Mcase in [
 ("invariant_plane",S.Matrix([[2,0,0],[0,3,1],[0,1,4]])),
 ("invariant_line",S.Matrix([[2,1,0],[1,3,0],[0,0,4]]))]:
    cvec=S.Matrix([0,0,1]);ell=S.Matrix([[1,0,0]])
    K=S.Matrix.hstack(cvec,Mcase*cvec,Mcase**2*cvec)
    check(name+" cyclic space lies in radical hyperplane",mz(ell*K) and K.rank() in [1,2])
# The full-rank qualification is sharp in the radical case.
ps=[u,x,y,z];gam=S.Symbol("gam");Msharp=S.diag(1+S.I,1-S.I,2)
V=homogenize(Msharp*S.Matrix([x,y,z])+S.Matrix([0,0,x*y]),u,[x,y,z],2)
F=x*y;h=2*u;X=S.Matrix([-F*h,*list(F*(V-h*S.Matrix(ps)))])
sweep=X+gam*S.Matrix([1,*ps])
check("rank-one nonlinear normalized rank3 boundary control",
 zero(sweep.jacobian([gam,*ps]).det(method="domain-ge")-gam**4) and
 X.jacobian(ps).subs({p:1 for p in ps}).rank()==3 and all(zero(S.diff(e,z)) for e in X))
model={"claims":["KMFI130","KMFI131","KMFI132"],
 "previous_goal_turn":"progress: KMFI128--129 proved affine eigen-coordinate rigidity and closed all constant-diagonal triangular homogeneous projective directions.",
 "universal_companion_certificate":{"V":[str(e) for e in d["V"]],"D":str(D),
 "angular_eta":[str(e) for e in eta[1:]],"D_squared_xy_curl":str(curl),
 "D_x5_z_coefficient":-2,"curl_x2_y_z4_coefficient":8,
 "polynomial_parameter_ring":"Q[lambda,tau,sigma,delta,x,y,z]","delta_nonzero_required":True},
 "parameter_specializations":special,"all_quadratic_form_cases":cases,
 "triangular_quadratic_result":"no constant diagonal assumed; every full-rank survivor affine and primitive1",
 "negative_cofactor_result":"no polynomial for nonradical quadratic direction; also all-degree univariate longitudinal curvature",
 "rank4_image_dimension1_quadratic_class":"fully excluded to full rank for every homogeneous scalar degree",
 "no_normalized_multiplier_claimed_for_every_case":False,
 "rank_deficient_normalized_radical_control_retained":True,
 "quadratic_nonimage1_directions_classified":False,
 "all_V0_without_pure_square_classified":False,
 "new_abstract_action_claimed":False,"actual_PSL11_input_supplied":False,
 "absolute_exceptional_action_decided":False,"full_action_spectrum_complete":False,"goal_complete":False}
(ROOT/"rank_one_projective_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["RANK_ONE_PROJECTIVE_CURVATURE.md","rank_one_projective.py","verify_rank_one_projective.py","rank_one_projective_model.json",
 "primitive_multiplier.py","PRIMITIVE_MULTIPLIER_CRITERION.md","TRIANGULAR_PROJECTIVE_RIGIDITY.md","FACTORIZED_HOMOGENEOUS_FIELDS.md"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,"runtime_seconds":round(time.monotonic()-started,3),
 "python":platform.python_version(),"sympy":S.__version__,"sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 "assurance":"Unbounded triangular quadratic and negative-cofactor proofs, exhaustive rank4 rank-one quadratic normal-form cases and universal nonzero curl coefficient. Exact whole identities supplement written proofs; no formal verification, independent whole-proof replay or external review claimed.","goal_complete":False}
(ROOT/"rank_one_projective_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],"universal_curl_coefficient":8,"rank4_quadratic_projective_image1_closed":True,"goal_complete":False}))
