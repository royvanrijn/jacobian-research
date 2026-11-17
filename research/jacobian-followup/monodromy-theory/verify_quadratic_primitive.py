"""Exact translation, zero-stratum and complete quadratic gates KMFI155--157."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import sympy as S
import json,time,hashlib,platform
from quadratic_primitive import feedback,trace_equation,normal_limit,lower_rank_family
from unmatched_quadratic import canonical_forms
from divisor_homotopy import affine_survivor
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def zero(v):return S.cancel(v)==0
def mz(M):return all(zero(v) for v in M)
u,x,y,z,lam,eps=S.symbols("u x y z lam eps");ps=[u,x,y,z];p=S.Matrix(ps)
# General constant frame: preserve u, send a to the u axis and b to e_z.
aa=S.symbols("a0:4");a=S.Matrix(aa);b=S.Matrix([0,0,0,1])
C=S.Matrix([[1,0,0,0],[-aa[1],aa[0],0,0],[-aa[2],0,aa[0],0],[-aa[3]/aa[0],0,0,1]])
check("translation frame retains the critical linear form",C[0,:]==S.Matrix([[1,0,0,0]]))
check("translation frame sends a to its nonzero u component",mz(C*a-S.Matrix([aa[0],0,0,0])))
check("translation frame retains b",C*b==b)
check("translation frame is invertible at u(a) nonzero",zero(C.det()-aa[0]**2))
F=S.Function("F")(u,x,y,z);J=S.Function("J")(u,x,y)
V0=S.Matrix([u*u+x*y,u*x+y*y,x*x+u*y,x*y+u*u]);cc=S.Symbol("cc")
V=V0+S.Matrix([cc*u*z,0,0,0]);R=u+x*x+y
whole=trace_equation(V,ps,F,R)
rest=trace_equation(V0,ps,F,R)
check("entire arbitrary-polynomial density operator split",
      zero(whole-rest-cc*z*(u*S.diff(F,u)+F)))
k=S.Symbol("k",integer=True,nonnegative=True)
check("all nonnegative Euler degrees have eigenvalue k+1",
      S.simplify(u*S.diff(u**k,u)+u**k-(k+1)*u**k)==0)
rho=S.Matrix([u*z,0])
check("translation rational control has exact zero divergence",
      zero(S.diff((1/u)*rho[0],u)+S.diff((1/u)*rho[1],z)))
# Every canonical family has the translation derivative, with arbitrary A.
mm=S.symbols("m0:16");A=S.Matrix(4,4,mm);forms=canonical_forms(ps)
for name,K in forms.items():
    V=A*K
    check(name+" whole last-direction translation derivative",mz(V.diff(z)-u*A[:,3]))
    check(name+" whole divergence last-direction coefficient",
          zero(S.diff(S.trace(V.jacobian(ps)),z)-A[0,3]))
# Matrix determinant lemma at arbitrary normal codimension; entire identities.
normaldet=[]
for codim in [1,2,3]:
    qq=S.Matrix(S.symbols("q0:"+str(codim)));hh=S.Matrix([S.symbols("eta0:"+str(codim))])
    bb=S.symbols("b0:"+str(codim*codim));B=S.Matrix(codim,codim,bb)
    N=B+(B*qq)*hh
    check("c"+str(codim)+" entire normal rank-one determinant",
          zero(N.det(method="domain-ge")-B.det()*(1+(hh*qq)[0])))
    normaldet.append({"codimension":codim,"identity":"detB*(1+eta*q)"})
# Full zero-stratum limits retain the upper-right block and remove lower blocks.
B=S.Matrix([[y+1,z],[z+1,y+2]]);q=S.Matrix([u,x]);Ctan=S.Matrix([[y,z+1],[z,y+1]])
V=S.Matrix([*list(B*q+S.Matrix([u*u+u*x,x*x+u*u])),
            *list(Ctan*q+S.Matrix([u*x,x*x+u*u]))])
h=u+x+u*u
initials=[y+z+1,u*(y+1)+x*z,u*u+u*x*y+x*x*(z+1),u**3+x**3*y+u*x*x*z]
limits=[]
for m,Jm in enumerate(initials):
    FF=Jm+(u+x)**(m+1)
    data=normal_limit(V,ps,FF,h,2,eps);lim=data["limit"]
    rowq=S.Matrix([[S.diff(Jm,v)/Jm for v in [u,x]]])
    rowt=S.Matrix([[S.diff(Jm,v)/Jm for v in [y,z]]])
    Bq=B*q;normal=B+Bq*rowq
    top=(Bq.jacobian([y,z])+Bq*rowt)
    expected=S.BlockMatrix([[normal,top],[S.zeros(2,2),S.zeros(2,2)]]).as_explicit()
    check("order"+str(m)+" entire zero-stratum gauge limit",mz(lim-expected))
    check("order"+str(m)+" leading Euler degree exact",zero((rowq*q)[0]-m))
    check("order"+str(m)+" whole nonzero normal determinant",
          zero(normal.det(method="domain-ge")-(m+1)*B.det()) and B.det()!=0)
    check("order"+str(m)+" tangential correction retained",not mz(top))
    limits.append({"initial_order":m,"whole_limit":"checked","determinant":"(m+1)*detB","upper_right_retained":True})
# Sharp boundary controls for the zero-stratum lemma.
Vsing=S.Matrix([x,0,0]);p3=[u,x,y]
N=feedback(Vsing,p3,S.Integer(1),S.Integer(0))
check("singular normal derivative admits genuine polynomial nilpotency",mz(N**2) and not mz(N))
Vrat=S.Matrix([u,0]);Nr=feedback(Vrat,[u,x],1/u,S.Integer(0))
check("negative rational normal order admits zero feedback",mz(Nr))
# Exact mixed-chain normal derivative along its full affine zero plane.
V=A*forms["B"];Bmix=V[:2,:].jacobian([u,x]).subs({u:0,x:0})
Bexpected=S.Matrix([[A[0,2]*y+A[0,3]*z,A[0,3]*y],
                    [A[1,2]*y+A[1,3]*z,A[1,3]*y]])
check("mixed chain whole normal derivative along u=x0",mz(Bmix-Bexpected))
check("mixed chain entire normal determinant with all parameters",
      zero(Bmix.det()-(A[0,2]*A[1,3]-A[0,3]*A[1,2])*y*y))
check("mixed chain full zero-stratum condition",mz(V.subs({u:0,x:0})))
check("mixed chain normal determinant after translation gate",
      zero(Bmix.det().subs(A[0,3],0)-A[0,2]*A[1,3]*y*y))
# Both single shears, and the d0 chain branch, reduce to a matched row.
eta=S.Matrix([S.symbols("eta0:4")])
for name in ["B","C","D"]:
    AA=A.copy();AA[0,3]=0
    if name=="B":AA[0,2]=0
    VV=AA*forms[name];L=S.cancel(VV[0]/u)
    NN=VV.jacobian(ps)+VV*eta-lam*u*S.eye(4)
    check(name+" entire regular invariant-critical-hyperplane row",
          mz(NN[:1,:].subs(u,0)-S.Matrix([[L.subs(u,0),0,0,0]])))
# Whole remaining mixed-chain trace and highest (y,z)-degree operator.
AA=A.copy();AA[0,3]=AA[1,3]=0
VV=AA*forms["B"];d=AA[0,2];f=AA[1,2]
div=S.trace(VV.jacobian(ps))
expected=(2*AA[0,0]+AA[1,1]+AA[2,2]+AA[3,3])*u+(AA[0,1]+2*AA[1,2]+AA[2,3])*x+d*y
check("entire remaining mixed-chain divergence",zero(div-expected))
scale=S.Symbol("scale")
for m in [0,2,4]:
    top=sum((u**(j+1)+x*x+j+1)*y**j*z**(m-j) for j in range(m+1))
    FF=top+(u+x) if m else top
    T=trace_equation(VV,ps,FF,4*lam*u).subs({y:scale*y,z:scale*z},simultaneous=True)
    coef=S.Poly(T,scale).coeff_monomial(scale**(m+1))
    raising=u*y*(d*S.diff(top,u)+f*S.diff(top,x))+d*y*top
    check("plane-degree"+str(m)+" complete highest trace coefficient",zero(coef-raising))
D=lambda expr:d*S.diff(expr,u)+f*S.diff(expr,x)
check("last Euler change preserves its transverse invariant",zero(D(x-f*u/d)))
W=S.Matrix([d*u*y,f*u*y,0,0])
check("last leading field satisfies the translation lemma exactly",
      mz(W.diff(y)-u*S.Matrix([d,f,0,0])) and zero(S.trace(W.jacobian(ps))-d*y))
# Exact normalized lower-rank controls in both newly closed families.
gam=S.Symbol("gam");controls=[]
for name in ["B","D"]:
    for ss in [1,3]:
        data=lower_rank_family(ps,ss,name);X=data["X"];P=data["P"]
        check(name+"s"+str(ss)+" whole P strictly lower triangular",
              all(zero(P[r,co]) for r in range(4) for co in range(r,4)))
        check(name+"s"+str(ss)+" entire nilpotency",mz(P**4))
        sweep=X+gam*S.Matrix([1,*ps])
        check(name+"s"+str(ss)+" entire normalized sweep Jacobian",
              zero(sweep.jacobian([*ps,gam]).det(method="domain-ge")-gam**4))
        rank=3 if name=="B" else 2
        check(name+"s"+str(ss)+" exact input rank and missing coordinate",X.jacobian(ps).rank()==rank and all(S.diff(v,z)==0 for v in X))
        check(name+"s"+str(ss)+" primitive quadratic direction",S.gcd_list(list(data["V"]))==1)
        controls.append({"family":name,"s":ss,"lambda":ss+2,"rankDX":rank,"entire_sweep_Jacobian":"gamma^4"})
# A genuine full-rank quadratic presentation remains allowed on the nonprimitive boundary.
pos=affine_survivor();pp=pos["ps"];uu,xx,yy,zz=pp;ll=pos["l"];T=pos["T"];red=pos["red"]
chi=(T-1)*pos["C"].charpoly(T).as_expr()
check("nonprimitive positive control coefficient polynomial irreducible",pos["modulus"].is_irreducible)
check("nonprimitive positive control entire spectral feedback identity",
      red(5*chi-T*S.diff(chi,T)-(T-ll)**4)==0)
grad=S.Matrix([[S.diff(pos["G"],v) for v in pp]])
Mp=S.eye(4);Cp=S.eye(3)
for j in range(1,5):
    Mp=(Mp*pos["M"]).applyfunc(red);Cp=(Cp*pos["C"]).applyfunc(red)
    check("nonprimitive positive control whole norm moment"+str(j),
          red((grad*Mp*S.Matrix(pp))[0]-(1+S.trace(Cp))*pos["G"])==0)
point={uu:1,xx:1,yy:0,zz:0}
check("nonprimitive positive control exact nonzero full-rank minor",
      red(pos["X"][1:,:].jacobian(pp).subs(point).det(method="domain-ge")-ll**4)==0 and red(ll**4)!=0)
check("nonprimitive positive control has coordinate gcd u",zero(S.gcd_list(list(pos["V"]))-uu))
check("nonprimitive positive control H scalar-linear factorization exact",
      mz(pos["F"]*pos["V"]-pos["G"]*pos["M"]*S.Matrix(pp)))
model={"claims":["KMFI155","KMFI156","KMFI157"],
"previous_goal_turn":"progress: KMFI153--154 closed the entire arbitrary-matrix single-square family.",
"translation_density_lemma":{"V_b":"u*a","u(b)":0,"u(a)":"nonzero","R_b":0,
"F":"nonzero polynomial, no homogeneity or degree bound","conclusion":"div(FV)=R*F impossible",
"highest_coefficient":"u*J_u=-J","rational_negative_degree_control_retained":True,
"all_four_canonical_families_gate":"A03=0"},
"zero_stratum_lemma":{"V_and_h_vanish_on_linear_stratum":True,"normal_derivative":"invertible over k(t)",
"F":"any nonzero polynomial","leading_normal_order":"any m>=0",
"normal_block":"B+(Bq)*d_q logJ","whole_normal_determinant":"(m+1)*detB",
"arbitrary_codimension_determinant_controls":normaldet,"entire_gauge_limits":limits,
"upper_right_tangential_block_retained":True,"singular_normal_and_negative_rational_controls_retained":True},
"single_shear_short_closure":{"families":["K_C","K_D"],"all_A_and_scalar_degrees":True,
"A03=0_then_invariant_critical_hyperplane_forces_matched_square":True},
"mixed_chain_closure":{"A03":0,"A02":"nonzero unless matched","normal_determinant":"A02*A13*y^2",
"normal_lemma_forces_A13":0,"highest_y_z_degree":"negative u Euler eigenvalue","all_A_and_scalar_degrees":True},
"complete_rank4_quadratic_result":{"all_four_primitive_families_closed":True,"primitive_degree2_full_rank":"impossible",
"nonprimitive_survivors":"V=u*M*p, G=u*F, exactly KMFI121 with G divisible by u",
"degree11_rank4_group":"A11 by KMFI122","all_degree_field":"exact KMFI121 one-variable polynomial",
"new_unbounded_group_classification_claimed":False,"homogeneous_rank4_exceptional11_minimum_primitive_direction_degree":3,
"general_cubic_and_higher_inputs_solved":False,"genuinely_nonhomogeneous_inputs_solved":False},
"actual_normalized_lower_rank_controls":controls,
"actual_nonprimitive_full_rank_control":{"source":"divisor_homotopy.affine_survivor","normalization":"entire spectral identity and four norm moments by KMFI121",
"full_rank_minor":"l^4","coordinate_gcd":"u","expanded_full_sweep_Jacobian_claimed":False},
"new_occurrence_action_claimed":False,"actual_PSL11_input_supplied":False,"absolute_exceptional_action_decided":False,
"full_action_spectrum_complete":False,"goal_complete":False}
(ROOT/"quadratic_primitive_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["QUADRATIC_PRIMITIVE_RIGIDITY.md","quadratic_primitive.py","verify_quadratic_primitive.py","quadratic_primitive_model.json",
"MIXED_BOUNDARY_MULTIPLIERS.md","UNMATCHED_QUADRATIC_CHARTS.md","unmatched_quadratic.py",
"QUADRATIC_PURE_SQUARE_RIGIDITY.md","PRIMITIVE_MULTIPLIER_CRITERION.md","CHAINED_SQUARE_QUADRATIC_RIGIDITY.md",
"FACTORIZED_HOMOGENEOUS_FIELDS.md","divisor_homotopy.py"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,"runtime_seconds":round(time.monotonic()-started,3),
"python":platform.python_version(),"sympy":S.__version__,"sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
"assurance":"Unbounded written translation-density and normal-zero-stratum obstructions, complete primitive quadratic rank4 closure, exact entire gauge and polynomial trace identities, actual lower-rank whole sweep Jacobians and retained nonprimitive full-rank norm control. No formal verification, independent whole-proof replay or external review claimed.","goal_complete":False}
(ROOT/"quadratic_primitive_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],
"entire_primitive_quadratic_rank4_class_closed_to_full_rank":True,"goal_complete":False}))
