"""Exact fixed-locus, ruling and trace identities for KMFI153--154."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
from itertools import combinations
import sympy as S
import json,time,hashlib,platform
from single_square import direction,wedge,triple,fixed_quadric_form,rank3_control
from unmatched_quadratic import canonical_forms
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def zero(v):return S.cancel(v)==0
def mz(M):return all(zero(v) for v in M)
u,x,y,z,lam=S.symbols("u x y z lam");ps=[u,x,y,z];p=S.Matrix(ps)
# Complete linear coefficient kernels for the elementary exterior argument.
exterior=[]
for dim in [3,4]:
    zz=S.symbols("p0:"+str(dim));pairs=list(combinations(range(dim),2))
    coeff=S.symbols("B0:"+str(len(pairs)*dim))
    Bij={ij:sum(coeff[j*dim+k]*zz[k] for k in range(dim)) for j,ij in enumerate(pairs)}
    eq=[]
    for i,j,k in combinations(range(dim),3):
        q=zz[i]*Bij[j,k]-zz[j]*Bij[i,k]+zz[k]*Bij[i,j]
        eq+=S.Poly(q,*zz).coeffs()
    mat,_=S.linear_eq_to_matrix(eq,coeff)
    check("n"+str(dim)+" entire linear two-form kernel dimension",len(coeff)-mat.rank()==dim)
    cv=S.symbols("c0:"+str(dim));model={ij:zz[ij[0]]*cv[ij[1]]-zz[ij[1]]*cv[ij[0]] for ij in pairs}
    check("n"+str(dim)+" all p-wedge-c forms satisfy exterior identity",
          all(zero(zz[i]*model[j,k]-zz[j]*model[i,k]+zz[k]*model[i,j]) for i,j,k in combinations(range(dim),3)))
    cn=S.symbols("C0:"+str(len(pairs)));const=dict(zip(pairs,cn));eq=[]
    for i,j,k in combinations(range(dim),3):
        eq+=S.Poly(zz[i]*const[j,k]-zz[j]*const[i,k]+zz[k]*const[i,j],*zz).coeffs()
    cm,_=S.linear_eq_to_matrix(eq,cn)
    check("n"+str(dim)+" constant two-form kernel is zero",cm.rank()==len(cn))
    # Entire fixed quadratic map normal form and its rank-one matrix.
    aa=S.symbols("a0:"+str(dim));rr=S.symbols("r0:"+str(dim));kap=S.Symbol("kap")
    data=fixed_quadric_form(zz,kap,aa,rr);pp=S.Matrix(zz)
    check("n"+str(dim)+" whole fixed-quadric wedge form",
          mz(wedge(pp,data["V"])-data["D"]*wedge(pp,S.Matrix(aa))))
    bb=S.Matrix([0,*S.symbols("b1:"+str(dim))])
    rb=(S.Matrix([rr])*bb)[0];res=data["M"]*bb-S.Matrix(aa)
    check("n"+str(dim)+" transverse Mb relation forces rb1",
          zero(res[0]-aa[0]*(rb-1)))
    check("n"+str(dim)+" whole remaining Mb relation forces kappa*b0",
          mz(res-S.Matrix(aa)*(rb-1)-kap*bb))
    check("n"+str(dim)+" rank-one fixed matrix determinant identity",
          zero(data["M"].det(method="domain-ge")-kap**(dim-1)*(kap+sum(aa[j]*rr[j] for j in range(dim)))))
    exterior.append({"dimension":dim,"linear_two_form_kernel":dim,"constant_two_form_kernel":0})
# A precise square-curvature countercontrol when transversality is dropped.
W=S.Matrix([u*u,u*x,u*y,2*u*z+2*x*x]);D=u*z+2*x*x
check("nontransverse square control is in canonical family",
      mz(W-direction(S.diag(1,1,1,2),ps)))
check("nontransverse control fixes a nonlinear quadric pointwise",
      mz(W-u*p-S.Matrix([0,0,0,D])))
check("nontransverse control whole weighted divisor identity",
      zero(D.subs(dict(zip(ps,W)),simultaneous=True)-2*u*u*D))
check("nontransverse control full line homotopy identity",
      zero(D.subs(dict(zip(ps,p+S.Symbol("t")*W)),simultaneous=True)-D*(1+3*u*S.Symbol("t")+2*u*u*S.Symbol("t")**2)))
check("nontransverse square control entire determinant",
      zero(W.jacobian(ps).det(method="domain-ge")-4*u**4))
# Arbitrary output matrix remains throughout.
mm=S.symbols("m0:16");M=S.Matrix(4,4,mm);V=direction(M,ps);a=M[:,3];b=S.Matrix([0,0,0,1])
check("arbitrary square direction equals canonical form",mz(V-M*canonical_forms(ps)["C"]))
check("entire square direction Jacobian for arbitrary M",
      zero(V.jacobian(ps).det(method="domain-ge")-2*M.det()*u**4))
check("canonical curvature relation Mb=a",M*b==a)
check("entire base line lies in canonical base locus",mz(V.subs({u:0,x:0})))
check("outside base line boundary image is curvature direction",
      mz(V.subs(u,0)-a*x*x))
# The critical-hyperplane row needs only regularity of dlogF there.
R=S.Matrix(V);MM=M.copy();MM[0,3]=0;R=direction(MM,ps)
eta=S.Matrix([S.symbols("eta0:4")])
N=R.jacobian(ps)+R*eta-lam*u*S.eye(4)
L=S.cancel(R[0]/u)
check("entire invariant critical-boundary feedback row",
      mz(N[:1,:].subs(u,0)-S.Matrix([[L.subs(u,0),0,0,0]])))
# Generic base point: the degree-one cancelled line map sends q to Mq.
ss,tt,eps=S.symbols("s t eps");q=S.Matrix([0,0,ss,tt]);mq=M*q
line=q+eps*mq
vl=V.subs(dict(zip(ps,line)),simultaneous=True)
cancelled=(vl/line[0]).applyfunc(S.cancel)
check("whole cancelled line map value at its base point",
      mz(cancelled.subs(eps,0)-mq))
# Variable-base ruling gives a Segre quadric or a plane; the skew case has
# an explicit limiting line violating u preservation.
Delta=M[0,2]*M[1,3]-M[0,3]*M[1,2]
q0=S.Matrix([0,0,M[0,3],-M[0,2]]);r0=M*q0
basis=S.Matrix.hstack(S.Matrix([0,0,1,0]),b,M[:,2],M[:,3])
check("directrix skewness determinant is exact Delta",zero(basis.det()-Delta))
check("limiting base point has zero critical coordinate",q0[0]==0)
check("limiting image base point also has zero critical coordinate",zero(r0[0]))
check("limiting image point has nonzero x precisely in skew case",zero(r0[1]+Delta))
vr=V.subs(dict(zip(ps,r0)),simultaneous=True)
check("whole critical-coordinate limiting-line contradiction",
      zero(vr[0]-M[0,3]*Delta**2))
al,be=S.symbols("al be")
segre=S.Matrix([al*ss,al*tt,be*ss,be*tt])
check("corresponding-point family whole Segre quadric equation",
      zero(segre[0]*segre[3]-segre[1]*segre[2]))
check("two skew directrices parameterize the entire corresponding line",
      mz(basis*segre-(al*q+be*M*q)))
# Constant-base case: quotient by a genuine invariant eigenline in B.
mu=S.Symbol("mu");Q=M.copy()
for j in range(4):Q[j,2]=0
Q[2,2]=mu
VV=direction(Q,ps);indices=[0,1,3];Qbar=Q.extract(indices,indices)
check("whole invariant vertex eigenline",Q[:,2]==S.Matrix([0,0,mu,0]))
check("quotient square direction has exactly the same form",
      mz(VV.extract(indices,[0])-(u*Qbar*S.Matrix([u,x,z])+Qbar[:,2]*x*x)))
check("quotient matrix invertible whenever original is",
      zero(Q.det()-mu*Qbar.det()))
check("quotient keeps curvature outside critical hyperplane",Qbar[0,2]==M[0,3])
check("quotient direction independent of vertex coordinate",
      all(S.diff(v,y)==0 for v in VV.extract(indices,[0])))
# Exact linear-factor identities imply all F annihilate a and b.
dd=S.symbols("d0:4");row=S.Matrix([dd]);c=S.Symbol("c")
Dp=(row*p)[0];res=(row*V)[0]-c*u*Dp
check("whole weighted linear-factor coefficient identity",
      zero(res-u*((row*M-c*row)*p)[0]-(row*a)[0]*x*x))
check("linear-factor x-square coefficient is D(a)",
      zero(S.Poly(res,*ps).coeff_monomial(x*x)-(row*a)[0]))
check("entire Mb relation forces D(b)0 for every eigenfactor",
      zero((row*a)[0]-c*(row*b)[0]-((row*M-c*row)*b)[0]))
# Differentiated trace identity is whole, with an arbitrary function F.
div=S.trace(V.jacobian(ps))
expected=u*S.trace(M)+(M*p)[0]+2*x*a[1]
check("entire direction divergence identity",zero(div-expected))
check("whole constant-base-direction derivative of V",mz(V.diff(z)-u*a))
check("whole divergence derivative is transverse curvature",zero(S.diff(div,z)-a[0]))
F=S.Function("F")(u,x,y)
DaF=sum(a[j]*S.diff(F,ps[j]) for j in range(4))
trace=sum(V[j]*S.diff(F,ps[j]) for j in range(4))+F*div-4*lam*u*F
check("entire differentiated trace for arbitrary z-independent F",
      zero(S.diff(trace,z)-u*DaF-a[0]*F))
I1=a[1]*u-a[0]*x;I2=a[2]*u-a[0]*y;FF=I1**2*I2**3
Da=lambda f:sum(a[j]*S.diff(f,ps[j]) for j in range(4))
check("all-factor directional-annihilation model exact",
      zero(S.diff(FF,z)) and zero(Da(FF)))
T=sum(V[j]*S.diff(FF,ps[j]) for j in range(4))+FF*div-4*lam*u*FF
check("entire final trace contradiction under both annihilations",
      zero(S.diff(T,z)-a[0]*FF))
# Actual normalized rank3 control over the exact quadratic coefficient field.
beta=S.Symbol("beta");minimal=beta*beta-beta+1
def reduced(v):
    num,den=S.cancel(v).as_numer_denom()
    assert not den.has(beta)
    return S.Poly(num,beta).rem(S.Poly(minimal,beta)).as_expr()
def rz(v):return reduced(v)==0
def rmz(A):return all(rz(v) for v in A)
data=rank3_control(ps,beta);X=data["X"];P=data["P"];la=data["lam"];gam=S.Symbol("gam")
check("control coefficient field irreducible and parameters nonzero",
      S.discriminant(minimal,beta)==-3 and S.gcd(minimal,beta*(1+beta))==1)
check("control lambda-square relation exact",rz(la**2-3*beta))
check("control whole nilpotent P",rmz(P**4))
regular=not rmz(P**3)
check("control actual normalized nonconstant P",not rmz(P) and all(S.diff(v,z)==0 for v in X))
sweep=X+gam*S.Matrix([1,*ps])
check("control entire sweep Jacobian in exact coefficient field",
      rz(sweep.jacobian([*ps,gam]).det(method="domain-ge")-gam**4))
minor=X.jacobian(ps).extract([0,3,4],[0,2,1]).det()
check("control exact nonzero rank3 minor",
      rz(minor-2*la**2*u*x*y**3) and reduced(minor)!=0)
check("control primitive canonical single-square direction",
      S.gcd_list(list(data["V"]))==1 and mz(data["V"]-data["M"]*canonical_forms(ps)["C"]))
model={"claims":["KMFI153","KMFI154"],
"previous_goal_turn":"progress: KMFI151--152 closed the entire arbitrary-matrix chained-square family.",
"fixed_hypersurface_lemma":{"dimensions":"every n>=3","map":"u*M*p+a*ell^2","M":"invertible","Mb":"a","u(b)":0,"u(a)":"nonzero",
"exterior_coefficient_kernels_checked":exterior,"degree2_fixed_form":"L*p+D*c","higher_degree_fixed_form":"L*p",
"rank_one_M_contradiction":True,"nontransverse_square_fixed_quadric_control_retained":True},
"single_square_family":{"direction":"A*(u^2,ux,uy,uz+x^2)","A":"arbitrary GL4",
"matched_critical_hyperplane_branch_closed":True,
"variable_base_ruling":"skew directrices produce a limiting line in u0 sent outside it",
"whole_offending_critical_value":"M03*Delta^2",
"constant_base_ruling":"invariant cone vertex gives fixed nonlinear quotient curve, excluded by the n3 lemma",
"all_nonlinear_scalar_factors_excluded":True,"individual_factor_degree_bound_used":False,
"linear_factors":"D(a)=D(b)=0 and D*M=cD",
"whole_scalar":"partial_a F=partial_b F=0",
"differentiated_trace":"u(a)*F=0, contradiction"},
"actual_normalized_control":{"coefficient_field":"Q(beta), beta^2-beta+1=0","lambda":"1+beta","F":"y","A":"diag(1,lambda,beta,lambda)",
"rankDX":3,"entire_sweep_Jacobian":"gamma^4","P_generically_regular_nilpotent":regular,"minor":"2*lambda^2*u*x*y^3"},
"remaining_primitive_quadratic_families":["K_B mixed chain","K_D product shear"],
"remaining_scalar_form":"nonconstant, coprime to u, irreducible degrees1..3 and unbounded multiplicities",
"all_quadratic_directions_closed":False,"new_occurrence_action_claimed":False,"actual_PSL11_input_supplied":False,
"absolute_exceptional_action_decided":False,"full_action_spectrum_complete":False,"goal_complete":False}
(ROOT/"single_square_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["SINGLE_SQUARE_QUADRATIC_RIGIDITY.md","single_square.py","verify_single_square.py","single_square_model.json",
"MIXED_BOUNDARY_MULTIPLIERS.md","UNMATCHED_QUADRATIC_CHARTS.md","unmatched_quadratic.py",
"QUADRATIC_PURE_SQUARE_RIGIDITY.md","PRIMITIVE_MULTIPLIER_CRITERION.md","CHAINED_SQUARE_QUADRATIC_RIGIDITY.md"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,"runtime_seconds":round(time.monotonic()-started,3),
"python":platform.python_version(),"sympy":S.__version__,"sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
"assurance":"Unbounded written transverse fixed-hypersurface lemma and arbitrary-matrix square-shear obstruction; exact exterior kernels, entire ruling-limit coordinates and linear-factor/trace identities, with a whole normalized rank3 sweep Jacobian over its exact coefficient field. No formal verification, independent whole-proof replay or external review claimed.","goal_complete":False}
(ROOT/"single_square_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],
"arbitrary_matrix_single_square_family_closed_to_full_rank":True,"goal_complete":False}))
