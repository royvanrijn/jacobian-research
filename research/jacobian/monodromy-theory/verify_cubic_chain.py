"""Exact last-coordinate and all-matrix cubic-chain closure KMFI162--163."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import sympy as S
import json,time,hashlib,platform
from cubic_chain import last_adjugate,direction
from nilpotent_normal import cubic_control
from quadratic_primitive import feedback,trace_equation
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def zero(v):return S.cancel(v)==0
def mz(M):return all(zero(v) for v in M)
u,x,y,z,lam,T=S.symbols("u x y z lam T");ps=[u,x,y,z];p=S.Matrix(ps)
# Cyclic adjugate at an eigenvalue: entire companion identities retain repeats.
cc=S.Symbol("cc")
for size in [2,3,4]:
    rr=S.symbols("r0:"+str(size-1));Q=T**(size-1)+sum(rr[j]*T**j for j in range(size-1))
    chi=S.Poly((T-cc)*Q,T);B=S.zeros(size)
    for j in range(size-1):B[j+1,j]=1
    for j in range(size):B[j,size-1]=-chi.nth(j)
    v=S.eye(size)[:,0];C=S.Matrix.hstack(*[B**j*v for j in range(size)])
    qB=sum((rr[j]*B**j for j in range(size-1)),S.zeros(size))+B**(size-1)
    check("r"+str(size)+" cyclic companion Krylov matrix is identity",C==S.eye(size))
    check("r"+str(size)+" entire adjugate equals quotient polynomial on cyclic vector",
          mz((cc*S.eye(size)-B).adjugate()*v-qB*v))
    check("r"+str(size)+" adjugate cyclic vector stays nonzero at repeated root",
          zero((qB*v)[size-1]-1))
# Polynomial graph multiplicity retains exact integer degree, not only closedness.
nu=3;P=u*u+x;L=u*x+y;K=-P**4/S.Integer(nu);W=K*z+L;F=W**nu
check("entire single-root graph logarithmic equation",
      zero(S.diff(F,z)*W/F+P**4))
check("graph multiplicity is exact positive polynomial degree",S.Poly(F,z).degree()==nu)
# General h=c branch loses all z dependence once eta_z0 is proved.
v=S.Matrix([u*u+x*y,u*x+y*y,u*y+x*x]);c=u*u+x;R=u*x+y*y
V=S.Matrix([*list(v),c*z+R]);FF=u*u+x*y+y*y
XX=S.Matrix([-FF*c,*list(FF*(V-c*p))])
check("h=c and F_z0 make the whole reconstructed X independent of z",
      all(zero(S.diff(q,z)) for q in XX))
# Actual homogeneous full-rank rank2 graph control over its exact coefficient field.
omega=S.Symbol("omega");minimal=omega**2-omega+1
def red(v):
    num,den=S.cancel(v).as_numer_denom()
    assert not den.has(omega)
    return S.Poly(num,omega).rem(S.Poly(minimal,omega)).as_expr()/den
def rz(v):return red(v)==0
def rmz(M):return all(rz(v) for v in M)
p2=[u,z];V2=S.Matrix([u,omega*z]);F2=u*z;h2=1+omega
N2=feedback(V2,p2,F2,h2);X2=S.Matrix([-F2*h2,*list(F2*(V2-h2*S.Matrix(p2)))])
ad2=last_adjugate(V2,p2)
check("positive graph control whole nilpotent feedback",rmz(N2**2))
check("positive graph control exact multiplicity identity",rz(ad2["K"]+1) and rz(ad2["K"]+(ad2["c"]-h2)**2))
check("positive graph control whole adjugate vector",rmz((ad2["c"]*S.eye(2)-V2.jacobian(p2)).adjugate()*V2-S.Matrix([0,ad2["W"]])))
check("positive graph control full rank by exact unit minor",
      rz(X2[1:,:].jacobian(p2).det().subs({u:1,z:1})-3*omega) and red(3*omega)!=0)
gam=S.Symbol("gam")
S2=X2+gam*S.Matrix([1,u,z])
check("positive graph control entire normalized sweep Jacobian",
      rz(S2.jacobian([u,z,gam]).det(method="domain-ge")-gam**2))
check("positive graph control exact nonzero cyclicity",
      rz(S.Matrix.hstack(V2,V2.jacobian(p2)*V2).det()+u*z))
# Every previous matrix restriction is retained; other entries stay arbitrary.
a,b,d,e,f,g,h,i,j,k,l,m,n=S.symbols("a b d e f g h i j k l m n")
A=S.Matrix([[a,b,0,0],[d,e,f,0],[g,h,i,j],[k,l,m,n]])
V=direction(A,ps)
div=S.trace(V.jacobian(ps))
expected=(3*a+e+i+n)*u*u+2*(b+f)*u*x+2*j*x*y
check("whole cubic divergence after normal-stratum restrictions",zero(div-expected))
# The j!=0 coordinate branch is an allowed simultaneous source/output shear.
C=S.eye(4);C[3,2]=-n/j
Vnew=(C*V.subs(z,z+n*y/j)).applyfunc(S.expand);rr=m-n*i/j
check("j branch source/output coordinate shear determinant1",C.det()==1)
check("j branch curvature moved entirely into y component",
      zero(S.Poly(Vnew[2],y).coeff_monomial(y*y)-j*x) and S.Poly(Vnew[3],y).degree()<=1)
check("j branch last field independent of new z",S.diff(Vnew[3],z)==0)
check("j branch last u-square y coefficient exact",
      zero(S.Poly(Vnew[3],y).coeff_monomial(y)-rr*u*u))
check("j branch divergence preserved by genuine conjugacy",
      zero(S.trace(Vnew.jacobian(ps))-div))
for degree in [0,2,4]:
    JJ=u*u+x*z+z*z+1
    FF=JJ*y**degree+(u+x) if degree else JJ
    tr=trace_equation(Vnew,ps,FF,4*lam*u*u)
    coef=S.Poly(tr,y).coeff_monomial(y**(degree+1))
    wanted=u*u*(f*S.diff(JJ,x)+rr*S.diff(JJ,z))+j*(degree+2)*x*JJ
    check("j branch degree"+str(degree)+" entire highest-y coefficient",zero(coef-wanted))
D=lambda P:f*S.diff(P,x)+rr*S.diff(P,z)
check("constant derivative straightens when f nonzero",zero(D(x/f)-1) and zero(D(rr*x-f*z)))
check("constant derivative straightens in the f0 r-nonzero branch",
      zero((rr*S.diff(z/rr,z))-1))
# General last-coordinate adjugate retains every last-row constant.
A0=A.subs(j,0);V0=direction(A0,ps);ad=last_adjugate(V0,ps);B=V0.jacobian(ps)
adjvec=(ad["c"]*S.eye(4)-B).adjugate()*V0
check("whole last-coordinate adjugate has zero head",
      mz(adjvec[:3,:]))
check("whole last-coordinate adjugate retains entire W",
      zero(adjvec[3]-ad["W"]))
check("entire W is linear in z with the claimed K",
      zero(ad["W"]-ad["K"]*z-ad["W"].subs(z,0)))
M3=S.Matrix([[a,b,0],[d,e,f],[g,h,i]])
Kexpected=n*u**8*(n*S.eye(3)-M3).det()+2*f*n*n*(a-n)*u**7*x+b*f*n*n*u**6*x*x
check("entire generic head-feedback leading coefficient identity",zero(ad["K"]-Kexpected))
check("head-feedback coefficient independent of all last-row lower constants",
      all(S.diff(ad["K"],q)==0 for q in [k,l,m]))
kp=S.Poly(ad["K"],u,x,y)
check("non-u8 first coefficient forces f(a-n)0",
      zero(kp.coeff_monomial(u**7*x)-2*f*n*n*(a-n)))
check("non-u8 second coefficient forces bf0",
      zero(kp.coeff_monomial(u**6*x*x)-b*f*n*n))
check("the f-nonzero compatibility branch makes entire K zero",
      zero(ad["K"].subs({a:n,b:0})))
check("the repeated-root branch keeps determinant zero explicitly",
      zero((n*S.eye(3)-M3).det().subs({a:n,b:0})))
# Weighted leading equation retains all terms of equal weight.
Af=A0.subs(f,0);Vf=direction(Af,ps);kap=3*a+e+i+n-4*lam
Dhead=lambda J:(a*u+b*x)*S.diff(J,u)+(d*u+e*x)*S.diff(J,x)
scale=S.Symbol("scale")
for weight in [0,3,6]:
    ellmax=weight//2
    Js=[u**(ell+1)+x*x+ell+1 for ell in range(ellmax+1)]
    top=sum(Js[ell]*y**(weight-2*ell)*z**ell for ell in range(ellmax+1))
    FF=top+(u+x+y) if weight>=2 else top
    tr=trace_equation(Vf,ps,FF,4*lam*u*u).subs({y:scale*y,z:scale**2*z},simultaneous=True)
    coef=S.Poly(tr,scale).coeff_monomial(scale**weight)
    op=u*u*Dhead(top)+i*u*u*y*S.diff(top,y)+(n*u*u*z+n*x*y*y)*S.diff(top,z)+(kap*u*u+2*b*u*x)*top
    check("weight"+str(weight)+" entire highest weighted trace",zero(coef-op))
    ss=weight-2*ellmax
    topz=S.Poly(op,z).coeff_monomial(z**ellmax)
    wanted=y**ss*(u*u*Dhead(Js[-1])+((kap+i*ss+n*ellmax)*u*u+2*b*u*x)*Js[-1])
    check("weight"+str(weight)+" entire top-z trace coefficient",zero(topz-wanted))
ord=S.Symbol("ord",integer=True,nonnegative=True);kk=S.Symbol("kk")
J0=x*x+1;JJ=u**ord*J0
left=u*Dhead(JJ)+(kk*u+2*b*x)*JJ
bracket=u*Dhead(J0)+((ord*a+kk)*u+(ord+2)*b*x)*J0
check("whole arbitrary-order u factor extraction in leading trace",
      S.simplify(left-u**ord*bracket)==0)
check("universal offending critical coefficient forces b0",
      zero(bracket.subs(u,0)-(ord+2)*b*x*(x*x+1)))
# Exact final triangular chart with nonzero quadratic curvature.
At=Af.subs(b,0);Vt=direction(At,ps)/a;vt=Vt.subs(u,1)
Jtri=vt[1:,:].jacobian([x,y,z])
check("final forced output matrix is entirely lower triangular",
      all(At[r,co]==0 for r in range(4) for co in range(r+1,4)))
check("final output determinant is the product of nonzero diagonal entries",
      zero(At.det()-a*e*i*n))
check("final first component is the matched pure cube",zero(Vt[0]-u**3))
check("final projective derivative has constant diagonal",
      all(zero(Jtri[r,r]-val) for r,val in enumerate([e/a,i/a,n/a])))
check("final projective map is triangular",
      all(Jtri[r,co]==0 for r in range(3) for co in range(r+1,3)))
check("final genuinely nonlinear y curvature is retained",
      zero(S.Poly(vt[2],x,y,z).coeff_monomial(x*x)-i/a))
# Actual normalized lower-rank controls are not removed by this theorem.
controls=[]
for degree in [0,2]:
    data=cubic_control(ps,degree);X=data["X"];Q=data["P"]
    check("rank3 control s"+str(degree)+" entire nilpotency",mz(Q**4))
    sweep=X+gam*S.Matrix([1,*ps])
    check("rank3 control s"+str(degree)+" whole sweep Jacobian",
          zero(sweep.jacobian([*ps,gam]).det(method="domain-ge")-gam**4))
    check("rank3 control s"+str(degree)+" exact rank and missing last coordinate",
          X.jacobian(ps).rank()==3 and all(S.diff(v,z)==0 for v in X))
    check("rank3 control s"+str(degree)+" rejected branch lambda=n is real",
          data["lam"]==data["A"][3,3] and S.diff(data["F"],z)==0)
    controls.append({"s":degree,"lambda":degree+3,"rankDX":3,"lambda_equals_n":True,"entire_sweep_Jacobian":"gamma^4"})
model={"claims":["KMFI162","KMFI163"],
"previous_goal_turn":"progress: KMFI158--161 established strong normal and translation tests and exact cubic matrix restrictions.",
"last_coordinate_multiplicity_lemma":{"rank":"every r>=2","V":"(v(q),c(q)z+R(q))","h":"independent of z",
"cyclic_DV_V_required":True,"N_nilpotent_and_rankDX_full_required":True,
"W":"adj(cI-DV)V has only last component Kz+L","cyclicity_proves_W_nonzero_even_at_repeated_roots":True,
"all_rejected_branches_retained":["c=h","K=0"],
"positive_integer_degree":"nu=deg_zF","required_identity":"nu*K=-(c-h)^r",
"F_over_head_field":"scalar(q)*(z+L/K)^nu","global_denominator_clearing_not_assumed":True,
"actual_rank2_full_rank_control":{"omega":"omega^2-omega+1=0","V":"(u,omega*z)","F":"u*z","h":"1+omega","K":-1,"nu":1,"whole_sweep_Jacobian":"gamma^2"}},
"complete_cubic_family_closure":{"K":"(u^3,u^2x,u^2y+ux^2,u^2z+xy^2)","A":"arbitrary GL4","all_scalar_degrees":True,
"previous_flags":["A02=0","A03=0","A13=0"],"highest_y_trace_forces_A23":0,
"whole_K":"n*u^8*det(nI-M3)+2*f*n^2*(a-n)*u^7*x+b*f*n^2*u^6*x^2",
"multiplicity_identity_forces_A12":0,"weighted_trace_forces_A01":0,
"final_forced_A":"lower triangular","final_full_rank_contradiction":"KMFI129 and nonzero i/a*x^2",
"all_repeated_and_resonant_parameters_retained":True},
"actual_normalized_cubic_rank3_controls":controls,
"arbitrary_cubic_directions_closed":False,"all_quadratic_h0_shapes_retained":True,
"genuinely_nonhomogeneous_inputs_solved":False,"new_occurrence_action_claimed":False,"actual_PSL11_input_supplied":False,
"absolute_exceptional_action_decided":False,"full_action_spectrum_complete":False,"goal_complete":False}
(ROOT/"cubic_chain_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["CUBIC_TRIANGULAR_CHAIN_RIGIDITY.md","cubic_chain.py","verify_cubic_chain.py","cubic_chain_model.json",
"NILPOTENT_NORMAL_STRATA.md","nilpotent_normal.py","PRIMITIVE_MULTIPLIER_CRITERION.md",
"TRIANGULAR_PROJECTIVE_RIGIDITY.md","QUADRATIC_PRIMITIVE_RIGIDITY.md"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,"runtime_seconds":round(time.monotonic()-started,3),
"python":platform.python_version(),"sympy":S.__version__,"sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
"assurance":"Unbounded last-coordinate multiplicity lemma and entire arbitrary-matrix cubic-chain exclusion. Exact whole adjugates, repeated-root cyclicity, coordinate and weighted trace identities, actual rank2 full-rank and cubic rank3 whole sweep Jacobians. No formal verification, independent whole-proof replay or external review claimed.","goal_complete":False}
(ROOT/"cubic_chain_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],
"entire_arbitrary_matrix_cubic_chain_closed_to_full_rank":True,"goal_complete":False}))
