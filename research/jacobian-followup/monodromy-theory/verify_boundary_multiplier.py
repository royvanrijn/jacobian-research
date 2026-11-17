"""Exact weighted-plane and boundary controls for KMFI146--148."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import sympy as S
import json,time,hashlib,platform
from boundary_multiplier import weighted_matrix,boundary_limit,rank3_family
from unmatched_quadratic import canonical_forms
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def zero(e):return S.cancel(e)==0
def mz(M):return all(zero(e) for e in M)
u,x,y,z,e,lam,a=S.symbols("u x y z e lam a");ps=[u,x,y,z];xs=[x,y,z]
ell=S.Matrix([[S.Symbol("l0"),S.Symbol("l1"),S.Symbol("l2")]]);L=a*u+(ell*S.Matrix(xs))[0]
Q=S.Matrix([x*x+y*y,x*y+z*z,x*z+y*y])
R=S.Matrix([x+y,y+z,z+x]);b=S.Matrix([1,2,3])
V=S.Matrix([u*L,*list(Q+u*R+u*u*b)])
pure=boundary_limit(V,ps,e,lam,S.Integer(1))
C=weighted_matrix(Q,ell,xs,e)
check("entire pure-power rescaled boundary matrix",mz(pure["limit"]-C))
check("pure-power gauge removes every source u denominator",
 all(not S.fraction(v)[1].has(u) for v in pure["rescaled"]))
G=u*u+u*(x+y+z)+x*x+y*z
mixed=boundary_limit(V,ps,e,lam,G);G0=G.subs(u,0)
xi=S.Matrix([[S.diff(G0,v)/G0 for v in xs]])
outer=S.Matrix([(ell*S.Matrix(xs))[0],*list(Q)])*S.Matrix([[0,*list(xi)]])
check("entire mixed-factor logarithmic boundary correction",mz(mixed["limit"]-C-outer))
check("mixed correction is genuinely present and cannot be discarded",not mz(outer))
# Explicit row restrictions and independent target boundary for every A*K form.
forms=canonical_forms(ps);mm=S.symbols("m0:16")
restrictions=[]
for name,K in forms.items():
    A=S.Matrix(4,4,mm)
    entries=[(0,2),(0,3)] if name in ["A","B"] else [(0,3)]
    ufree=(A*K)[0].subs(u,0)
    monoms=[x*x,(y*y if name=="A" else x*y)] if name in ["A","B"] else [x*x if name=="C" else x*y]
    check(name+" exact boundary-factor row constraints",
       all(S.Poly(ufree,x,y,z).coeff_monomial(mon)==A[i,j] for (i,j),mon in zip(entries,monoms)))
    for i,j in entries:A[i,j]=0
    VV=A*K;Qfull=VV.subs(u,0);alpha=S.Matrix([[1,0,0,0]])*A.adjugate()
    check(name+" full boundary vector annihilated by source coordinate",Qfull[0]==0)
    check(name+" full boundary vector annihilated by independent target covector",
       zero((alpha*Qfull)[0]) and any(v!=0 for v in alpha[1:]))
    restrictions.append({"type":name,"zero_entries":[list(v) for v in entries]})
# All coefficients are independent in the binary determinant identity.
s=S.Symbol("s");aa=S.symbols("a0:3");bb=S.symbols("b0:3");t=S.Symbol("t")
R1=aa[0]*x*x+aa[1]*x*y+aa[2]*y*y
R2=bb[0]*x*x+bb[1]*x*y+bb[2]*y*y
J=S.Matrix([R1,R2]).jacobian([x,y]).det()
C3=S.Matrix([[(s+1)*x,1,0],[s*R1,S.diff(R1,x),S.diff(R1,y)],
             [s*R2,S.diff(R2,x),S.diff(R2,y)]])
check("whole generic binary boundary determinant identity",zero(C3.det()-(s+2)*x*J/2))
r1=aa[0]+aa[1]*t+aa[2]*t*t;r2=bb[0]+bb[1]*t+bb[2]*t*t
check("whole binary Jacobian Wronskian identity",
 zero(J-2*x*x*(r1*S.diff(r2,t)-S.diff(r1,t)*r2).subs(t,y/x)))
# Coordinate conjugacy removing a proportional second component preserves the matrix form.
k=S.Symbol("k");Qp=S.Matrix([R1,k*R1,0]);row=S.Matrix([[1,0,0]])
B=S.Matrix([[1,0,0],[-k,1,0],[0,0,1]]);D=S.diag(1,B)
Cp=weighted_matrix(Qp,row,xs,s)
changed=(D*Cp.subs(dict(zip(xs,B.inv()*S.Matrix(xs))),simultaneous=True)*D.inv()).applyfunc(S.cancel)
Qnew=(B*Qp.subs(dict(zip(xs,B.inv()*S.Matrix(xs))),simultaneous=True)).applyfunc(S.cancel)
check("whole proportional-image conjugacy keeps weighted matrix structure",
 mz(changed-weighted_matrix(Qnew,row,xs,s)) and zero(Qnew[1]))
# Trace and determinant of the final two-dimensional block.
c=S.Symbol("c");RR=-(s+1)*x*x/2+c*y*y
C2=S.Matrix([[(s+1)*x,1],[s*RR,S.diff(RR,x)]])
check("trace-normalized rank-one block trace exact",zero(S.trace(C2)))
det2=S.expand(C2.det())
check("universal offending x-square coefficient",zero(S.Poly(det2,x,y).coeff_monomial(x*x)+(s+1)*(s+2)/2))
check("remaining y-square coefficient exact",zero(S.Poly(det2,x,y).coeff_monomial(y*y)+s*c))
# If ell annihilates the image plane, its nonzero scalar is already an eigenvalue.
lz,T=S.symbols("lz T");Qplane=S.Matrix([R1,R2,0])
Cann=weighted_matrix(Qplane,S.Matrix([[0,0,lz]]),xs,s)
D2=S.Matrix([R1,R2]).jacobian([x,y])
pred=T*(T-(s+1)*lz*z)*(T*S.eye(2)-D2).det()
check("annihilating-image-plane case has the claimed eigenvalue",zero(Cann.charpoly(T).as_expr()-pred))
# The two exceptional negative weights are not polynomial-scalar cases.
exceptional=[]
for weight,quad in [(-1,S.zeros(3,1)),(-2,S.Matrix([x*x/2,0,0]))]:
    E=weighted_matrix(quad,row,xs,S.Integer(weight))
    check("exceptional negative weight really permits nilpotency "+str(weight),mz(E*E))
    exceptional.append({"s":weight,"Q":[str(v) for v in quad],"nilpotency":"C^2=0","polynomial_scalar_case":False})
# Actual normalized nonlinear inputs show why full rank is indispensable.
controls=[]
gam=S.Symbol("gam")
for ss in [1,2,5]:
    data=rank3_family(ps,ss);P=data["P"];X=data["X"]
    check("s"+str(ss)+" whole P is strictly lower triangular",all(zero(P[i,j]) for i in range(4) for j in range(i,4)))
    check("s"+str(ss)+" generic nilpotent block is regular",not mz(P**3) and mz(P**4))
    sweep=X+gam*S.Matrix([1,*ps])
    check("s"+str(ss)+" entire normalized sweep Jacobian",
      zero(sweep.jacobian([gam,*ps]).det(method="domain-ge")-gam**4))
    at={p:1 for p in ps};minor=X.jacobian(ps).extract([0,3,4],[0,1,2]).det()
    check("s"+str(ss)+" rank exactly3",minor.subs(at)!=0 and all(zero(S.diff(v,z)) for v in X))
    check("s"+str(ss)+" direction is primitive quadratic and target chart matched",
      S.gcd_list(list(data["V"]))==1 and data["V"][0]==u*u)
    check("s"+str(ss)+" determinant degree relation exact",
      zero(data["V"].jacobian(ps).det(method="domain-ge")-2*data["lam"]**3*u**4) and data["lam"]**4==(ss+2)*data["lam"]**3)
    controls.append({"s":ss,"lambda":ss+2,"F":"u^"+str(ss),"V":[str(v) for v in data["V"]],
                     "X":[str(v) for v in X],"rankDX":3,"whole_Jacobian":"gam^4","regular_P":True})
model={"claims":["KMFI146","KMFI147","KMFI148"],
 "previous_goal_turn":"progress: KMFI142--145 gave four exhaustive quadratic chart forms, global weighted scalar preservation and the factor degree bound3.",
 "boundary_factor":{"e_positive_required":True,"conclusion":"u divides V_u","all_four_output_row_restrictions":restrictions,
 "mixed_limit":"C_e + column(ell*x,Q) row(0,dlogG0)","mixed_rank_one_term_retained":True},
 "plane_image_matrix_lemma":{"matrix":"[(s+1)*ell*x,ell;s*Q,DQ]","Q_image_dimension_at_most":2,
 "excluded_weights":[-1,-2],"all_nonnegative_s_covered":True,
 "binary_determinant":"(s+2)*x*detD(R1,R2)/2","rank_one_bad_coefficient":"-(s+1)*(s+2)/2",
 "exceptional_negative_weight_controls":exceptional},
 "pure_boundary_scalar":{"all_quadratic_rank4_directions":True,"all_s_ge0":True,
 "normalization_for_s_ge1_forces_matched_boundaries_without_full_rank":True,
 "full_rank_primitive_case":"excluded by KMFI141","s0":"full rank excluded by KMFI117"},
 "actual_nonlinear_rank3_boundary_family":controls,
 "pure_u_with_other_factors_claimed_excluded":False,
 "all_unmatched_low_degree_factor_arrangements_solved":False,
 "new_abstract_action_claimed":False,"actual_PSL11_input_supplied":False,
 "absolute_exceptional_action_decided":False,"full_action_spectrum_complete":False,"goal_complete":False}
(ROOT/"boundary_multiplier_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["PURE_BOUNDARY_MULTIPLIER.md","boundary_multiplier.py","verify_boundary_multiplier.py","boundary_multiplier_model.json",
 "UNMATCHED_QUADRATIC_CHARTS.md","unmatched_quadratic.py","QUADRATIC_PURE_SQUARE_RIGIDITY.md",
 "PRIMITIVE_MULTIPLIER_CRITERION.md","SWEEP_DIFFERENTIAL_RECONSTRUCTION.md"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,"runtime_seconds":round(time.monotonic()-started,3),
 "python":platform.python_version(),"sympy":S.__version__,"sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 "assurance":"Unbounded boundary-row and mixed-limit identities, plane-image nilpotency obstruction and pure-boundary scalar exclusion. Exact whole polynomial identities, negative-weight countercontrols and normalized rank3 examples supplement written proofs. No formal verification, independent whole-proof replay or external review claimed.","goal_complete":False}
(ROOT/"boundary_multiplier_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],
 "all_pure_boundary_scalar_exponents_closed_to_primitive_full_rank":True,"mixed_factors_excluded":False,"goal_complete":False}))
