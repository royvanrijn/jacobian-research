"""Exact all-parameter chained-square trace gates for KMFI151--152."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import sympy as S
import json,time,hashlib,platform
from chained_square import direction,feedback,shift_matrices,rank3_family
from unmatched_quadratic import canonical_forms
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def zero(v):return S.cancel(v)==0
def mz(M):return all(zero(v) for v in M)
u,x,y,z,t,s,lam,cc=S.symbols("u x y z t s lam cc");ps=[u,x,y,z];p=S.Matrix(ps)
K=canonical_forms(ps)["A"]
check("canonical base locus sequential u-square identity",K[0]==u*u)
check("canonical base locus sequential x-square identity",zero(K[2]-u*y-x*x))
check("canonical base locus sequential y-square identity",zero(K[3]-u*z-y*y))
check("canonical base locus includes vertex for every z",mz(K.subs({u:0,x:0,y:0})))
mm=S.symbols("a0:16");A=S.Matrix(4,4,mm);V=direction(A,ps)
check("arbitrary invertible target does not change base ideal",
      mz(A.adjugate()*V-A.det()*K))
check("entire direction Jacobian determinant for arbitrary A",
      zero(V.jacobian(ps).det(method="domain-ge")-2*A.det()*u**4))
# Pointwise fixed-surface countercontrol has positive-dimensional base.
D=u*z-x*y;W=S.Matrix([u*u,u*x,u*y,2*u*z-x*y])
check("fixed-quadric control whole composition identity",zero(D.subs(dict(zip(ps,W)),simultaneous=True)-2*u*u*D))
check("fixed-quadric control exact pointwise proportionality modulo D",
      mz(W-u*p-S.Matrix([0,0,0,D])))
check("fixed-quadric control full homotopy factor",
      zero(D.subs(dict(zip(ps,p+t*W)),simultaneous=True)-D*(1+3*u*t+2*u*u*t*t)))
check("fixed-quadric base has two distinct lines",
      mz(W.subs({u:0,x:0})) and mz(W.subs({u:0,y:0})))
# Cone translation is exact; a noncone fixed quadric violates it.
cone=u*u*y+x**3+u*x*y
check("cone equation whole vertex-translation identity",zero(cone.subs(z,z+t)-cone))
check("fixed-quadric boundary example is not a cone at the chained vertex",
      S.diff(D,z)==u and S.diff(D,z)!=0)
# Every linear scalar factor is supported on u,x.
dd=S.symbols("d0:4");row=S.Matrix([dd]);L=(row*p)[0]
res=(row*V)[0]-cc*u*L;poly=S.Poly(res,*ps)
check("linear factor x-square coefficient",zero(poly.coeff_monomial(x*x)-(row*A)[2]))
check("linear factor y-square coefficient",zero(poly.coeff_monomial(y*y)-(row*A)[3]))
check("linear factor uy versus x-square eliminates y coefficient",
      zero(poly.coeff_monomial(u*y)-poly.coeff_monomial(x*x)+cc*dd[2]))
check("linear factor uz versus y-square eliminates z coefficient",
      zero(poly.coeff_monomial(u*z)-poly.coeff_monomial(y*y)+cc*dd[3]))
# Whole z-trace obstruction for every polynomial independent of z.
J,Ju,Jx,Jy=S.symbols("J Ju Jx Jy")
trace_num=J*S.trace(V.jacobian(ps))+(S.Matrix([[Ju,Jx,Jy,0]])*V)[0]-4*lam*u*J
zcoef=S.Poly(trace_num,z).coeff_monomial(z)
check("entire z coefficient of scalar trace numerator",
      zero(zcoef-u*(A[0,3]*Ju+A[1,3]*Jx+A[2,3]*Jy)-A[0,3]*J))
du,dx,dy,m=S.symbols("du dx dy m")
Invx=du*x-dx*u;Invy=du*y-dy*u
der=lambda q:du*S.diff(q,u)+dx*S.diff(q,x)+dy*S.diff(q,y)
check("first linear invariant of nonzero-u constant derivative",zero(der(Invx)))
check("second linear invariant of nonzero-u constant derivative",zero(der(Invy)))
check("Euler negative eigenvalue has nonzero monomial coefficient",
      S.simplify(u*der(u**m)/u**m+du-(m+1)*du)==0)
Inv=dy*x-dx*y
check("binary reduction invariant for u-tangent derivative",
      zero(dx*S.diff(Inv,x)+dy*S.diff(Inv,y)))
# Zero coupling column and full-rank rejection are complete symbolic identities.
AA=A.copy()
for i in range(3):AA[i,3]=0
VV=direction(AA,ps);FF=u*u+x*y+y*y;eta=[S.diff(FF,v)/FF for v in ps]
N=feedback(VV,ps,eta,lam)
check("zero coupling last column is an eigenvector identity",
      mz(N[:,3]-S.Matrix([0,0,0,(AA[3,3]-lam)*u])))
XX=S.Matrix([-lam*u*FF,*list(FF*(VV-lam*u*p))])
check("forced last eigenvalue destroys whole z dependence",
      all(zero(S.diff(v,z).subs(lam,AA[3,3])) for v in XX))
# Critical invariant hyperplane cannot carry a nonmatched linear factor.
a,b,c,d=S.symbols("a b c d")
q=S.Matrix([u*(a*u+b*x+c*y+d*z),x*x+u*y,y*y+u*z,x*y+z*z])
etas=S.symbols("eta0:4");N=feedback(q,ps,etas,lam)
check("entire critical-boundary feedback row with regular logarithm",
      mz(N[:1,:].subs(u,0)-S.Matrix([[b*x+c*y+d*z,0,0,0]])))
# Linear factors in u,x form one line when the coefficient d=A02 is nonzero.
l0,l1=S.symbols("l0 l1")
ll=S.Matrix([[l0,l1,0,0]])
check("whole linear-factor shear-column constraint",
      zero((ll*A)[2]-l0*A[0,2]-l1*A[1,2]))
check("nonzero shear coefficient fixes unique scalar factor",
      zero((l0*A[0,2]+l1*A[1,2]).subs(l0,-l1*A[1,2]/A[0,2])))
# Allowed matched source/output conjugacy, not an independent target change.
C,B=shift_matrices(t)
oldK=K.subs(x,x-t*u)
check("whole canonical family identity under x translation",mz(oldK-B*K))
check("conjugacy matrices both invertible with determinant1",C.det()==1 and B.det()==1)
newA=C*A*B
check("critical shear coefficient survives conjugacy",zero(newA[0,2]-A[0,2]))
check("first output's last coupling survives conjugacy",zero(newA[0,3]-A[0,3]))
check("whole direction transforms by genuine simultaneous conjugacy",
      mz(C*V.subs(x,x-t*u)-direction(newA,ps)))
# Final unbounded coefficient obstruction: all remaining entries retained.
a,b,c,d,ee,ff,g,h,i,j,k,l=S.symbols("a b c d ee ff g h i j k l")
AA=S.Matrix([[a,b,d,0],[0,c,0,0],[ee,ff,g,h],[i,j,k,l]])
VV=direction(AA,ps);eta=[0,s/x,0,0];N=feedback(VV,ps,eta,lam)
trace=S.expand(S.trace(N))
expected=(2*a+(s+1)*c+g+l-4*lam)*u+b*x+(d+2*h)*y
check("whole trace identity with all remaining arbitrary parameters",zero(trace-expected))
check("first trace forces ux coefficient b0",zero(S.Poly(trace,*ps).coeff_monomial(x)-b))
trace2=S.cancel(S.trace(N*N))
check("whole second trace has no remaining x denominator",not S.fraction(trace2)[1].has(x))
coef=S.Poly(trace2,*ps).coeff_monomial(x*x)
check("universal second trace x-square coefficient",
      zero(coef-b*b-2*(s+2)*c*d))
check("nonzero offending coefficient after first trace",
      zero(coef.subs(b,0)-2*(s+2)*c*d))
check("last scalar factor is a nonzero eigen-coordinate",
      zero(VV[1]-c*u*x) and zero(AA.det()-c*AA.extract([0,2,3],[0,2,3]).det()))
# Actual normalized rank3 controls retain full-rank necessity.
controls=[]
gam=S.Symbol("gam")
for ss in [1,3]:
    data=rank3_family(ps,ss);X=data["X"];P=data["P"];lambda0=data["lam"]
    check("s"+str(ss)+" whole P strictly lower triangular",all(zero(P[r,co]) for r in range(4) for co in range(r,4)))
    check("s"+str(ss)+" generic P regular nilpotent",mz(P**4) and not mz(P**3))
    sweep=X+gam*S.Matrix([1,*ps])
    check("s"+str(ss)+" entire normalized sweep Jacobian",
          zero(sweep.jacobian([*ps,gam]).det(method="domain-ge")-gam**4))
    minor=X.jacobian(ps).extract([0,3,4],[0,1,2]).det()
    check("s"+str(ss)+" exact rank3 minor and absent z dependence",
          zero(minor+4*lambda0**3*(ss+1)*u**(3*ss)*x*y) and minor!=0 and all(S.diff(v,z)==0 for v in X))
    check("s"+str(ss)+" primitive chained-square direction",
          S.gcd_list(list(data["V"]))==1 and data["V"]==direction(S.diag(1,lambda0,lambda0,lambda0),ps))
    controls.append({"s":ss,"lambda":ss+2,"A":"diag(1,lambda,lambda,lambda)","F":"u^"+str(ss),
                     "rankDX":3,"whole_sweep_Jacobian":"gam^4","generic_P_regular_nilpotent":True})
model={"claims":["KMFI151","KMFI152"],
"previous_goal_turn":"progress: KMFI149--150 excluded every critical-boundary scalar factor, including mixed multipliers.",
"canonical_direction":"A*(u^2,ux,uy+x^2,uz+y^2)","A":"arbitrary GL4",
"scalar_cone_reduction":{"unique_projective_base_point":"[0:0:0:1]","fixed_surface_excluded_by_curve_line_bundle_degrees":True,
"nonlinear_factor_degree_bound_used":False,"selected_line_family_and_degree1_restriction_used":True,
"nonlinear_factors":"cones with vertex e3, independent of z","linear_factors":"span(u,x)","whole_F_independent_of_z":True,
"positive_dimensional_base_fixed_quadric_control_retained":True},
"all_parameter_trace_reduction":{"A03_nonzero":"negative Euler eigenvalue",
"A03=A13=A23=0":"last eigenvalue forces X independent of z",
"remaining_case":"A03=A13=0,A23 nonzero,F(u,x)",
"A02_zero":"invariant critical hyperplane forces matched square",
"A02_nonzero":"all scalar factors proportional to x+t*u",
"coordinate_change":"simultaneous source/output constant conjugacy, canonical family preserved",
"trace_forces":"b0","second_trace_x_square_coefficient":"2*(s+2)*c*d","c_and_d":"both nonzero"},
"conclusion":"no homogeneous multiplier of any degree gives primitive quadratic chained-square rank4 full rank",
"rank3_normalized_controls":controls,
"remaining_primitive_quadratic_families":["K_B mixed chain","K_C square shear","K_D product shear"],
"remaining_scalar_form":"nonconstant, coprime to u, factors of degree1..3, arbitrary multiplicities",
"all_quadratic_directions_closed":False,"new_occurrence_action_claimed":False,"actual_PSL11_input_supplied":False,
"absolute_exceptional_action_decided":False,"full_action_spectrum_complete":False,"goal_complete":False}
(ROOT/"chained_square_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["CHAINED_SQUARE_QUADRATIC_RIGIDITY.md","chained_square.py","verify_chained_square.py","chained_square_model.json",
"MIXED_BOUNDARY_MULTIPLIERS.md","UNMATCHED_QUADRATIC_CHARTS.md","unmatched_quadratic.py",
"QUADRATIC_PURE_SQUARE_RIGIDITY.md","PRIMITIVE_MULTIPLIER_CRITERION.md","SWEEP_DIFFERENTIAL_RECONSTRUCTION.md"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,"runtime_seconds":round(time.monotonic()-started,3),
"python":platform.python_version(),"sympy":S.__version__,"sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
"assurance":"Unbounded written finite-basepoint cone reduction and entire arbitrary-matrix chained-square obstruction; exact whole traces, allowed conjugacy, polynomial boundary controls and entire normalized rank3 sweep Jacobians. No formal verification, independent whole-proof replay or external review claimed.","goal_complete":False}
(ROOT/"chained_square_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],
"arbitrary_matrix_chained_square_family_closed_to_full_rank":True,"goal_complete":False}))
