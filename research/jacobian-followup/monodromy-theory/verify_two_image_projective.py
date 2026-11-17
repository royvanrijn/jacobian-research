"""Exact source-pinned reduction and coupled two-image curl gates, KMFI133--134."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import sympy as S
import json,time,hashlib,platform
from two_image_projective import nilpotent_normalizers,conjugate,companion_certificate
from rank_one_projective import homogenize
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def zero(e):return S.cancel(e)==0
def mz(M):return all(zero(e) for e in M)
x,y,z,u=S.symbols("x y z u");xs=[x,y,z];a,b,c,d=S.symbols("a b c d",nonzero=True)
H=S.Matrix([0,a*x*x,b*x*x+c*x*y+d*y*y])
SA,SB=nilpotent_normalizers(x,y,z,a,b,c,d)
check("generic nonzero y-square branch normalizes exactly",mz(conjugate(H,SA,xs)-S.Matrix([0,x*x,y*y])))
check("generic zero y-square mixed branch normalizes exactly",mz(conjugate(H.subs(d,0),SB,xs)-S.Matrix([0,x*x,x*y])))
check("square normalizer has exact nonzero determinant",zero(SA.det()-1/(d*a**3)))
check("mixed normalizer has exact nonzero determinant",zero(SB.det()-1/(a*a*c)))
# All M and all b are retained after normalizing the quadratic part.
mm=S.symbols("m0:9");M=S.Matrix(3,3,mm);bb=S.Matrix(S.symbols("b0:3"));xx=S.Matrix(xs)
normal=[]
for name,H0 in [("square",S.Matrix([0,x*x,y*y])),("mixed",S.Matrix([0,x*x,x*y]))]:
    J=H0.jacobian(xs)
    check(name+" whole nilpotent quadratic Jacobian",mz(J**3))
    check(name+" exact homogeneous scaling determinant",zero((S.eye(3)+S.Symbol("t")*J).det()-1))
    vec=M*(xx+H0)+bb;B=vec.jacobian(xs)
    check(name+" arbitrary linear mixing retained in whole derivative",mz(B-M*(S.eye(3)+J)))
    check(name+" entire projective determinant constant",zero(B.det(method="domain-ge")-M.det()))
    V=homogenize(vec,u,xs,2)
    check(name+" entire homogeneous Jacobian-power shape",zero(V.jacobian([u,*xs]).det(method="domain-ge")-2*M.det()*u**4))
    coeff=S.Matrix([[S.Poly(e,*xs).coeff_monomial(mon) for e in H0] for mon in [(2,0,0),(1,1,0),(0,2,0)]])
    check(name+" quadratic image exactly two-dimensional",coeff.rank()==2)
    L=S.Matrix([[1,1,0],[0,1,1],[0,0,1]])
    changed=conjugate(vec,L,xs)
    expected=L*M*L.inv()*(xx+conjugate(H0,L,xs))+L*bb
    check(name+" whole conjugacy preserves arbitrary M and constants",mz(changed-expected))
    cyc=S.Matrix([[0,1,0],[0,0,1],[1,0,0]])*(xx+H0)
    check(name+" triangular H need not give constant-diagonal triangular v",any(S.diff(S.trace(cyc.jacobian(xs)),v)!=0 for v in xs))
    normal.append({"type":name,"H":[str(e) for e in H0],"determinant_compatible_arbitrary_M_and_b":True})
manifest=json.loads((ROOT/"sources/two-image-projective-manifest.json").read_text())
primary=manifest["sources"][0];raw=(ROOT/primary["path"]).read_bytes()
check("primary theorem PDF source bytes and hash pinned",raw.startswith(b"%PDF") and len(raw)==primary["bytes"] and hashlib.sha256(raw).hexdigest()==primary["sha256"])
# Universal all-parameter companion coefficient certificate.
data=companion_certificate();u,x,y,z=data["ps"];lam,beta,gamma,tau,sig,delta=data["parameters"]
f=data["feedback"];D=data["D"];eta=data["eta"];R=data["curl_numerator"];T=f["T"]
check("coupled projective Jacobian is exactly delta",zero(data["v"].jacobian([x,y,z]).det()-delta))
check("coupled full homogeneous determinant is exactly 2delta u4",zero(data["V"].jacobian([u,x,y,z]).det(method="domain-ge")-2*delta*u**4))
check("whole feedback moment solve at u1",mz(eta*f["C"].subs(u,1)-f["moments"].subs(u,1)))
check("coupled universal beta obstruction coefficient",zero(R.coeff_monomial(x**4*z**5)+8*beta**2))
check("coupled universal gamma obstruction coefficient",zero(R.coeff_monomial(x*x*y**6*z**3)-32*gamma**4))
# Denominators in the numerator certificate may only use the nonzero delta.
allgens=[x,y,z,lam,beta,gamma,tau,sig,delta]
for name,expr in [("Nx",delta*data["Nx"]),("Ny",delta*data["Ny"]),("curl",delta**2*R.as_expr())]:
    try:S.Poly(S.cancel(expr),*allgens)
    except S.PolynomialError:raise AssertionError("unexpected denominator "+name)
    checks.append("no beta gamma or spectral denominator in "+name+" certificate")
Mcomp=S.Matrix([[0,1,0],[-sig,tau,1],[delta,0,0]])
quad=S.Matrix([x*x,beta*x*y+gamma*y*y,0])
check("companion nonlinear shapes meet both normal forms exactly",
 mz(Mcomp.inv()*quad-S.Matrix([0,x*x,-tau*x*x+beta*x*y+gamma*y*y])))
special=[]
for vals in [(1,0,0,0,1,1),(0,1,0,0,1,2),(-3,2,4,-2,3,0),(0,0,1,2,2,3),(0,0,0,0,1,1)]:
    # beta,gamma,tau,sigma,delta,lambda
    sub=dict(zip([beta,gamma,tau,sig,delta,lam],vals))
    dp=D.subs(sub)
    if dp==0:
        check("singular Krylov specialization takes the full-rank rejection branch "+str(vals),
              f["C"].subs(u,1).subs(sub).det()==0)
        special.append({"parameters":list(vals),"Krylov_singular":True})
        continue
    ep=eta.subs(sub);rp=R.as_expr().subs(sub)
    point={x:1,y:2,z:3}
    if dp.subs(point)==0:point={x:2,y:3,z:5}
    check("specialization keeps nonzero denominator "+str(vals),dp.subs(point)!=0)
    actual=(S.diff(ep[2],x)-S.diff(ep[1],y)).subs(point)
    check("polynomial numerator equals actual curl "+str(vals),zero(dp.subs(point)**2*actual-rp.subs(point)))
    special.append({"parameters":list(vals),"source_point":[point[x],point[y],point[z]],"numerator_identity_checked":True})
# At beta=gamma=0 a negative constant-derivative cofactor excludes every degree.
kappa=S.Symbol("kappa");J=y*y+y*z+1;F=x*x*J+x*(y+z)+1
vv=data["v"].subs({beta:0,gamma:0});W=vv-S.Matrix([x,y,z])
res=(S.Matrix([[S.diff(F,v) for v in [x,y,z]]])*W)[0]-(kappa-2*x)*F
coef=S.Poly(S.expand(res),x).coeff_monomial(x**3)
expected=-sig*S.diff(J,y)+delta*S.diff(J,z)+4*J
check("entire beta gamma zero field gives negative cofactor",zero(coef-expected))
Dj=-sig*S.diff(J,y)+delta*S.diff(J,z)
D2=-sig*S.diff(Dj,y)+delta*S.diff(Dj,z);D3=-sig*S.diff(D2,y)+delta*S.diff(D2,z)
check("constant derivative in cofactor argument is locally nilpotent",zero(D3))
# Explicit automorphy of both nonlinear normal forms is used in KMFI135.
xx=S.Matrix([x,y,z])
for name,Phi,Inv in [
 ("square",S.Matrix([x,y+x*x,z+y*y]),S.Matrix([x,y-x*x,z-(y-x*x)**2])),
 ("mixed",S.Matrix([x,y+x*x,z+x*y]),S.Matrix([x,y-x*x,z-x*(y-x*x)]))]:
    check(name+" explicit inverse in both directions",
      mz(Phi.subs(dict(zip([x,y,z],Inv)),simultaneous=True)-xx) and
      mz(Inv.subs(dict(zip([x,y,z],Phi)),simultaneous=True)-xx))
# A preserved hyperplane may still have a nonlinear map in the remaining variables.
q=x+1;vexample=S.Matrix([2*x+1,3*y,4*z+y*y]);Wexample=vexample-xx
check("affine divisor includes the projective constant correctly",
 zero((S.Matrix([[S.diff(q,w) for w in [x,y,z]]])*Wexample)[0]-q))
L=S.Matrix([[1,0,0],[0,1,0],[0,0,1]])
u0=S.Symbol("u0");pp=[u0,x,y,z]
Vexample=homogenize(vexample,u0,[x,y,z],2);li=u0+x
gradli=S.Matrix([[S.diff(li,w) for w in pp]])
check("preserved affine divisor homogenizes to eigenrow",
 zero((gradli*Vexample)[0]-2*u0*li))
# Tangent E to du=dl=0 is invariant; df of a split scalar in u,l kills it.
E=S.Matrix([[0,0],[0,0],[1,0],[0,1]])
check("homogeneous derivative preserves tangent complement",
 mz((gradli*Vexample.jacobian(pp)*E)) and
 mz(S.Matrix([[1,0,0,0]])*Vexample.jacobian(pp)*E))
Fsplit=u0**2*li**3;fscalar=-lam*u0*Fsplit
check("split scalar differential annihilates invariant complement",
 mz(S.Matrix([[S.diff(fscalar,w) for w in pp]])*E))
# The imported theorem must not be confused with the claimed new conclusion.
model={"claims":["KMFI133","KMFI134","KMFI135"],
 "previous_goal_turn":"progress: KMFI130--132 closed quadratic triangular inputs and rank4 projective quadratic curvature of one-dimensional nonlinear image.",
 "imported_primary_result":{"source":primary["path"],"url":primary["url"],"theorem":"3.2","applied_to":"H=M^-1*(homogeneous quadratic part)","applied_to_whole_v":False},
 "exhaustive_remaining_shapes":normal,"arbitrary_M_and_constants_b_still_open":True,
 "normalization_basis_matrices":{"square":[[str(e) for e in SA.row(i)] for i in range(3)],"mixed":[[str(e) for e in SB.row(i)] for i in range(3)]},
 "coupled_companion":{"v":[str(e) for e in data["v"]],"D":str(D),"Nx":str(data["Nx"]),"Ny":str(data["Ny"]),
 "curl_numerator":str(R.as_expr()),"beta_coefficient":"-8*beta^2","gamma_coefficient":"32*gamma^4",
 "source_polynomial_term_count":len(R.terms()),"parameters_other_than_delta_may_vanish":True,
 "delta_nonzero_required":True,"singular_Krylov_branch_excluded_by_full_rank":True,"all_homogeneous_multiplier_degrees_closed":True,
 "zero_beta_gamma_negative_cofactor_proved":True,"specializations":special},
 "split_scalar_result":{"all_ranks_and_degrees":True,"requires_projective_polynomial_automorphism":True,
 "completely_split_scalar_forces_affine_v":True,"both_explicit_normal_form_inverses_checked":True,
 "remaining_primitive_quadratic_scalar_requires_nonlinear_irreducible_factor":True},
 "all_arbitrary_M_b_families_excluded":False,"generic_groups_of_all_two_image_directions_computed":False,
 "new_abstract_action_claimed":False,"actual_PSL11_input_supplied":False,
 "absolute_exceptional_action_decided":False,"full_action_spectrum_complete":False,"goal_complete":False}
(ROOT/"two_image_projective_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["TWO_IMAGE_PROJECTIVE_REDUCTION.md","two_image_projective.py","verify_two_image_projective.py","two_image_projective_model.json",
 "sources/two-image-projective-manifest.json",primary["path"],"primitive_multiplier.py","rank_one_projective.py",
 "RANK_ONE_PROJECTIVE_CURVATURE.md","TRIANGULAR_PROJECTIVE_RIGIDITY.md"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,"runtime_seconds":round(time.monotonic()-started,3),
 "python":platform.python_version(),"sympy":S.__version__,"sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 "assurance":"Source-attributed constant-field triangularization plus unbounded two-shape reduction, exact normalizers, whole determinant/conjugacy controls and universal coupled companion curl coefficients. No formal verification, independent whole-proof replay or external review claimed.","goal_complete":False}
(ROOT/"two_image_projective_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],"remaining_nonlinear_shapes":2,"companion_beta_coefficient":"-8 beta^2","companion_gamma_coefficient":"32 gamma^4","goal_complete":False}))
