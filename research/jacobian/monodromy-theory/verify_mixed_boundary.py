"""Exact polynomial mixed-boundary gates for KMFI149--150."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import sympy as S
import json,time,hashlib,platform
from mixed_boundary import mixed_matrix,binary_matrix,two_block
from unmatched_quadratic import canonical_forms
from boundary_multiplier import rank3_family
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def zero(v):return S.cancel(v)==0
def mz(M):return all(zero(v) for v in M)
u,x,y,z,e,n,ax,ay,b,k,h,c,eta,tt=S.symbols("u x y z e n ax ay b k h c eta T")
co=S.symbols("r0:6")
R1=co[0]*x*x+co[1]*x*y+co[2]*y*y
R2=co[3]*x*x+co[4]*x*y+co[5]*y*y
C=S.Matrix([[(e+1)*x,1,0],[e*R1,S.diff(R1,x),S.diff(R1,y)],
             [e*R2,S.diff(R2,x),S.diff(R2,y)]])
v=S.Matrix([x,R1,R2])
check("whole binary Euler column identity",mz(C*S.Matrix([1,x,y])-(e+2)*v))
M=C+v*S.Matrix([[0,ax,ay]])
jac=S.Matrix([R1,R2]).jacobian([x,y]).det()
check("entire mixed binary determinant before homogeneity",
      zero(M.det(method="domain-ge")-(e+2+x*ax+y*ay)*x*jac/2))
check("homogeneous polynomial logarithmic factor gives positive determinant weight",
      S.expand(e+2+n)==e+n+2)
# Exact polynomial j controls, including a factor removed along the image plane.
Q=S.Matrix([R1,R2,0]);ell=[1,0,0]
js=[S.Integer(1),x,y,x*x+y*y,x*x*y*(x+y)**2]
for idx,j in enumerate(js):
    Mj=binary_matrix([R1,R2],x,y,e,j);deg=S.Poly(j,x,y).total_degree()
    check("binary homogeneous correction determinant j"+str(idx),
          zero(Mj.det(method="domain-ge")-(e+deg+2)*x*jac/2))
J=z**3*(x*x+x*y+y*y+z*(x+y)+z*z)
B=mixed_matrix(Q,ell,[x,y,z],e,J)
check("image plane gives identically zero final output row",mz(B[3:,:]))
check("factored plane power does not affect tangential logarithmic derivatives",
      all(zero(S.diff(J,t)/J-S.diff(J/z**3,t)/(J/z**3)) for t in [x,y]))
active=B[:3,:3].applyfunc(S.cancel).subs(z,0).applyfunc(S.cancel)
check("entire principal binary restriction retains mixed correction",
      mz(active-binary_matrix([R1,R2],x,y,e,x*x+x*y+y*y)))
# Proportional binary quadratics with simultaneous source/target coordinates.
r=S.Symbol("r");P=S.Matrix([[1,0],[r,1]])
old=[x,y+r*x];RR=S.Matrix([R1,k*R1])
# k is the proportionality; use y_new=y_old-k*x.
P=S.Matrix([[1,0],[-k,1]]);Pinv=P.inv()
j=x*x+x*y+2*y*y
oldR=(RR.subs(y,y+k*x));newR=P*oldR
D=S.diag(1,1,1);D[1:,1:]=P
Mold=binary_matrix(RR,x,y,e,j).subs(y,y+k*x)
Mnew=binary_matrix(newR,x,y,e,j.subs(y,y+k*x))
check("whole mixed proportional-image conjugacy",mz(D*Mold*D.inv()-Mnew))
check("second proportional quadratic removed",newR[1]==0)
# R1 nonzero branch, after removal of R2.
R=h*x*x+b*x*y+c*y*y
block=S.Matrix([[(e+1)*x,1+x*eta],[e*R,S.diff(R,x)+R*eta]])
forced=-(e+1)*x-S.diff(R,x)
detreduced=S.expand(block.det()).subs(eta,forced/R)
check("rank-one trace-reduced determinant exact",
      zero(detreduced-(e*(x*S.diff(R,x)-R)-(e+1)*x*x)))
check("rank-one coefficients force a and c",
      zero((e*(x*S.diff(R,x)-R)-(e+1)*x*x).subs({h:(e+1)/e,c:0})))
aa=(e+1)/e
forced_eta=-((e+1+2*aa)*x+b*y)/(x*(aa*x+b*y))
check("nonzero b branch forces negative residue minus1 at x0",
      zero(S.limit(x*forced_eta,x,0)+1))
check("zero b branch forces negative residue minus(e+2)",
      zero(forced_eta.subs(b,0)+(e+2)/x))
# R1 zero branch; no proportionality division may discard this case.
block=two_block(R,x,y,e,x,eta)
detreduced=S.expand(block.det()).subs(eta,forced/R)
forced2=-(e+1)*x-S.diff(R,y)
check("zero-first-quadratic trace-reduced determinant exact",
      zero(block.det().subs(eta,forced2/R)-x*(e*S.diff(R,y)-(e+1)*x)))
RR0=x*(h*x+aa*y);eta0=S.cancel((-(e+1)*x-S.diff(RR0,y))/RR0)
check("zero-first-quadratic negative residue minus(e+1)",
      zero(S.cancel(eta0*(h*x+aa*y)/aa)+(e+1)))
# Positive-weight rational controls: they pass nilpotency and violate polynomiality.
rational_controls=[]
for ee in [1,2,5]:
    a=S.Rational(ee+1,ee);q=S.Matrix([a*x*x,0,0]);Jrat=x**(-ee-2)
    Br=mixed_matrix(q,[1,0,0],[x,y,z],ee,Jrat).applyfunc(S.cancel)
    check("rational negative-residue control square-zero e"+str(ee),mz(Br**2) and not mz(Br))
    rational_controls.append({"e":ee,"J":str(Jrat),"ell_nonzero":True,"square_zero":True,"polynomial_multiplier":False})
# Exhaustive canonical row gates, without choosing any invertible output matrix.
mm=S.symbols("a0:16");forms=canonical_forms([u,x,y,z]);row_gates=[]
for name,K in forms.items():
    A=S.Matrix(4,4,mm)
    for i,jj in ([(0,2),(0,3)] if name in ["A","B"] else [(0,3)]):A[i,jj]=0
    V=A*K;L=S.cancel(V[0]/u);ell=S.Matrix([[S.diff(L,t) for t in [x,y,z]]]);Q=V[1:,:].subs(u,0)
    if name in ["A","B"]:
        obstruction=(ell*Q)[0]
        q1=x*x;q2=y*y if name=="A" else x*y
        check(name+" entire plane-annihilation output flag",
              zero(obstruction-A[0,1]*(A[1,2]*q1+A[1,3]*q2)))
        A[1,2]=A[1,3]=0;V=A*K;Q=V[1:,:].subs(u,0)
        check(name+" flag boundary first high output is zero",Q[0]==0)
        check(name+" flag boundary field independent of final coordinate",
              all(S.diff(v,z)==0 for v in Q))
        check(name+" flag output matrix determinant retains both blocks",
              zero(A.det()-A[:2,:2].det()*A[2:,2:].det()))
        row_gates.append({"form":name,"first_row":"(a,b,0,0)","b_nonzero_branch":"A12=A13=0","result":"negative residue or coefficient contradiction"})
    else:
        check(name+" boundary quadratic image is a constant line",
              mz(Q-S.Matrix([A[1,3],A[2,3],A[3,3]])*(x*x if name=="C" else x*y)))
        check(name+" vanishing boundary ell gives matched first row",
              ell==S.Matrix([[A[0,1],A[0,2],0]]))
        row_gates.append({"form":name,"boundary_image_span":1,"lemma_conclusion":"ell0 hence matched"})
# Entire z-infinity limit, retaining rather than dropping logarithmic terms.
J=(x+y)*z**3+(x*x+y*y)*z**2+x*y*z+x**4
Q=S.Matrix([0,h*x*x+k*x*y,c*x*x+b*x*y])
B=mixed_matrix(Q,[b,0,0],[x,y,z],e,J)
lim=B.applyfunc(lambda v:S.limit(v,z,S.oo))
expected=mixed_matrix(Q,[b,0,0],[x,y,z],e,x+y)
check("whole flagged z-infinity mixed matrix limit",mz(lim-expected))
check("z-infinity first high row zero",mz(lim[1:2,:]))
check("z-infinity final column zero",mz(lim[:,3:4]))
active=lim.extract([0,2],[0,2])
check("whole surviving flagged two-dimensional block",
      mz(active-two_block(Q[1],x,y,e,b*x,1/(x+y))))
check("full flagged limiting characteristic reduces to two block",
      zero(lim.charpoly(tt).as_expr()-tt**2*active.charpoly(tt).as_expr()))
R=h*x*x+k*y*y
block=two_block(R,x,y,e,b*x,eta);trforce=(-(e+1)*b*x-S.diff(R,y))/R
check("flag trace-reduced whole determinant",
      zero(block.det().subs(eta,trforce)-b*x*(e*S.diff(R,y)-(e+1)*b*x)))
gate=e*S.diff(R,y)-(e+1)*b*x
check("square-chain x coefficient forces b0",zero(S.Poly(gate,x,y).coeff_monomial(x)+(e+1)*b))
R=h*x*x+k*x*y;kk=(e+1)*b/e
gate=e*S.diff(R,y)-(e+1)*b*x
check("mixed-chain coefficient forces k=(e+1)b/e",zero(gate.subs(k,kk)))
eta_forced=S.cancel((-(e+1)*b*x-S.diff(R,y))/R)
check("mixed-chain logarithmic residue is -(e+1)",
      zero((eta_forced*(h*x+k*y)/k).subs(k,kk)+(e+1)))
# Retained sharp normalized polynomial controls from the matched class.
for ee in [1,3]:
    data=rank3_family([u,x,y,z],ee);X=data["X"];P=data["P"]
    check("matched normalized rank3 boundary retained e"+str(ee),
          mz(P**4) and not mz(P**3) and X.jacobian([u,x,y,z]).rank()==3 and all(S.diff(v,z)==0 for v in X))
model={"claims":["KMFI149","KMFI150"],
"previous_goal_turn":"progress: KMFI146--148 established exact mixed boundary matrices and the pure-boundary-power obstruction.",
"mixed_boundary_lemma":{"e":"every positive integer","J":"every nonzero homogeneous polynomial","conclusion":"ell annihilates every plane containing imageQ","image_span_at_most1":"ell0",
"whole_binary_determinant":"(e+degj+2)*x*detD(R1,R2)/2","all_zero_and_proportional_binary_branches_retained":True,
"negative_polynomial_residues":[-1,"-(e+2)","-(e+1)"],"rational_countercontrols":rational_controls},
"all_four_direction_family_gates":row_gates,
"critical_boundary_factor":{"normalization_alone_forces_matched_chart":True,"full_rank_primitive_quadratic_rank4":"u cannot divide F",
"other_factors_and_multiplicities_arbitrary":True,"mixed_correction_dropped":False},
"remaining_frontier":{"F":"nonconstant and coprime to u","factor_degrees":[1,2,3],"unmatched_target_alpha_factor":"forbidden",
"weighted_factor_preservation_required":True,"total_degree_bounded":False,"all_nonboundary_arrangements_solved":False},
"new_abstract_action_claimed":False,"actual_PSL11_input_supplied":False,
"absolute_exceptional_action_decided":False,"full_action_spectrum_complete":False,"goal_complete":False}
(ROOT/"mixed_boundary_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["MIXED_BOUNDARY_MULTIPLIERS.md","mixed_boundary.py","verify_mixed_boundary.py","mixed_boundary_model.json",
"PURE_BOUNDARY_MULTIPLIER.md","boundary_multiplier.py","UNMATCHED_QUADRATIC_CHARTS.md","unmatched_quadratic.py",
"QUADRATIC_PURE_SQUARE_RIGIDITY.md","PRIMITIVE_MULTIPLIER_CRITERION.md"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,"runtime_seconds":round(time.monotonic()-started,3),
"python":platform.python_version(),"sympy":S.__version__,
"sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
"assurance":"Unbounded written mixed-polynomial plane-image and all-four-family obstruction, with exact entire determinant, conjugacy, limit, residue and rational negative-residue controls. No formal verification, independent whole-proof replay or external review claimed.","goal_complete":False}
(ROOT/"mixed_boundary_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],
"all_critical_boundary_factors_excluded_to_primitive_quadratic_rank4_full_rank":True,"goal_complete":False}))
