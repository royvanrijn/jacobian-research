"""Exact divisor, rank-case and boundary controls for KMFI136--138."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import sympy as S
import json,time,platform,hashlib
from low_degree_multiplier import taylor_data,hyperplane_forms
from rank_one_projective import homogenize
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def zero(e):return S.cancel(e)==0
def mz(M):return all(zero(e) for e in M)
x,y,z,U,Vv,Zv=S.symbols("x y z U Vv Zv");xs=[x,y,z]
co=S.symbols("q0:10")
q=co[0]*x*x+co[1]*x*y+co[2]*x*z+co[3]*y*y+co[4]*y*z+co[5]*z*z+co[6]*x+co[7]*y+co[8]*z+co[9]
formal=taylor_data(q,S.Matrix([U,Vv,Zv]),xs)
check("entire quadratic Taylor identity for arbitrary target coordinates",zero(S.expand(formal["taylor_identity"])))
beta,gamma,delta=S.symbols("beta gamma delta")
v=S.Matrix([y+x*x,z+beta*x*y+gamma*y*y,delta*x])
data=taylor_data(q,v,xs)
check("entire second directional derivative identity",zero(S.expand(data["chain_identity"])))
check("entire quadratic Taylor identity for coupled nonlinear map",zero(S.expand(data["taylor_identity"])))
# The discrete preservation and Darboux hypotheses must not be conflated.
vbad=S.Matrix([y,x+y*y,z+x*x-y*y]);qbad=x*x+z
bad=taylor_data(qbad,vbad,xs)
check("automorphism preserves quadratic without Darboux condition",
 zero(bad["qv"]-qbad) and zero(bad["Wq"]+(x-y)**2))
invbad=S.Matrix([y-x*x,x,z-(y-x*x)**2+x*x])
check("quadratic preservation countercontrol has an actual polynomial inverse",
 mz(vbad.subs(dict(zip(xs,invbad)),simultaneous=True)-S.Matrix(xs)) and
 mz(invbad.subs(dict(zip(xs,vbad)),simultaneous=True)-S.Matrix(xs)))
check("discrete preservation alone does not supply divisor residue",
 S.rem(S.Poly(bad["Wq"],z),S.Poly(qbad,z)).as_expr()!=0)
# Rank-one leading part: simultaneous identities can hold but destroy cyclicity.
c=S.Symbol("c",nonzero=True);q1=x*x+y
v1=S.Matrix([x,c*y+(c-1)*x*x,z+x*x])
one=taylor_data(q1,v1,xs)
check("rank-one quadratic preservation and Darboux are exact",
 zero(one["qv"]-c*q1) and zero(one["Wq"]-(c-1)*q1))
check("rank-one quadratic control has fixed affine coordinate",zero(one["W"][0]))
K=S.Matrix.hstack(one["W"],one["B"]*one["W"],one["B"]**2*one["W"])
check("rank-one fixed coordinate makes projective Krylov singular",zero(K.det()))
# Rank-two leading part: the product residue forces a fixed coordinate.
a,b,c=S.symbols("a b c",nonzero=True);q2=x*y+z
v2=S.Matrix([a*x,b*y,c*z+(c-a*b)*x*y])
two=taylor_data(q2,v2,xs)
check("rank-two quadratic preservation exact",zero(two["qv"]-c*q2))
check("rank-two extra Darboux condition forces affine fixed direction",
 zero(two["Wq"]-(c-1)*q2+(a-1)*(b-1)*x*y))
check("rank-two Hessian residue is product of affine displacements",
 zero((two["W"].T*two["H"]*two["W"])[0]/2-(a-1)*(b-1)*x*y))
# Rank-three form cannot contain a plane; a plane-isotropic Gram matrix is singular.
s0,s1,s2=S.symbols("s0 s1 s2")
gram=S.Matrix([[s0,s1,s2],[s1,0,0],[s2,0,0]])
check("any plane-isotropic three-dimensional symmetric Gram matrix is singular",gram.det()==0)
HA=S.Matrix([0,x*x,y*y]);HB=S.Matrix([0,x*x,x*y])
for name,H in [("square",HA),("mixed",HB)]:
    check(name+" quadratic image has generic differential rank2",H.jacobian(xs).rank()==2)
check("rank2 form can contain an isotropic plane so rank3 hypothesis matters",zero((x*y).subs(x,0)))
# Hyperplane eigen-coordinate forces first-row matrix shape in both normal forms.
ell=S.Matrix(S.symbols("l0:3")).T
M=S.Matrix(3,3,S.symbols("m0:9"))
for name,H in [("square",HA),("mixed",HB)]:
    polynomial=(ell*M*H)[0]
    coeffs=S.Poly(polynomial,x,y,z)
    check(name+" independent quadratic shapes identify two output columns",
      coeffs.coeff_monomial(x*x)==(ell*M)[1] and
      coeffs.coeff_monomial(y*y if name=="square" else x*y)==(ell*M)[2])
# Every matrix and every affine constant in the hyperplane flag is retained.
forms=hyperplane_forms();X,Y,Z=forms["xs"]
a,b,c,d,e,f,g,k1,k2,k3=forms["parameters"];r=forms["r"];kk=S.Symbol("kappa")
expectedA=S.Matrix([(a-1)*X+k1,
 b*X+(c+g-1)*Y+d*Z+c*X*X+d*Y*Y+k2,
 (e-g*b/d)*X+r*Y-Z+r*X*X+k3-g*k2/d])
expectedB=expectedA.copy();expectedB[1]=expectedB[1]-d*Y*Y+d*X*Y
check("entire generic hyperplane square alignment",mz(forms["WA"]-expectedA))
check("entire generic hyperplane mixed alignment",mz(forms["WB"]-expectedB))
check("invertibility makes aligned constant derivative coefficient nonzero",
 zero(forms["M"].det()+a*d*r))
check("square trace cofactor has only the indicated varying term",
 zero(S.trace(forms["vA"].jacobian([X,Y,Z]))-(a+c+g+2*d*Y)))
check("mixed trace cofactor has only the indicated varying term",
 zero(S.trace(forms["vB"].jacobian([X,Y,Z]))-(a+c+g+d*X)))
for name,old in [("square",forms["oldA"]),("mixed",forms["oldB"])]:
    atzero=old.subs(d,0);B0=atzero.jacobian([X,Y,Z])
    check(name+" d0 branch is exactly constant-diagonal triangular",
       all(zero(B0[i,j]) for i in range(3) for j in range(i+1,3)) and
       all(not B0[i,i].has(X,Y,Z) for i in range(3)))
# Derive the negative cofactor from whole generic fields, not a sampled search.
J=X*X+X*Z+Z*Z+1;F=Y**2*J+Y*X**7+1
res=(S.Matrix([[S.diff(F,v) for v in [X,Y,Z]]])*forms["WA"])[0]-(kk-2*d*Y)*F
check("square whole-field highest Y equation",
 zero(S.Poly(S.expand(res),Y).coeff_monomial(Y**3)-r*S.diff(J,Z)-4*d*J))
# For the mixed case transverse plane degree must be taken before X degree.
J=Y*Y+Y*Z+Z*Z
F=X**2*J+X*(2*Y*Y+Z*Z)+X**8*Y+1
res=S.expand((S.Matrix([[S.diff(F,v) for v in [X,Y,Z]]])*forms["WB"])[0]-(kk-d*X)*F)
top=S.Add(*[coef*Y**mon[0]*Z**mon[1] for mon,coef in S.Poly(res,Y,Z).terms() if sum(mon)==2])
check("mixed whole-field transverse-top then X-leading equation",
 zero(S.Poly(S.expand(top),X).coeff_monomial(X**3)-d*(Y*S.diff(J,Y)+J)))
# An actual normalized nonlinear two-image input with a linear multiplier has rank3.
u,x,y,z,tgam,alpha=S.symbols("u x y z tgam alpha");ps=[u,x,y,z]
mod=S.Poly(alpha**2-alpha+1,alpha);lam=1+alpha
def reduce(e):
    num,den=S.fraction(S.cancel(e))
    return S.cancel(S.rem(S.Poly(num,alpha),mod).as_expr()/S.rem(S.Poly(den,alpha),mod).as_expr())
vv=S.Matrix([alpha*x,lam*y+x*x,lam*z+x*y])
V=homogenize(vv,u,[x,y,z],2);F=x;h=lam*u
surface=S.Matrix([-F*h,*list(F*(V-h*S.Matrix(ps)))])
P=surface[1:,:].jacobian(ps)-S.Matrix(ps)*surface[:1,:].jacobian(ps)
check("rank3 boundary input is normalized by whole nilpotency",all(reduce(e)==0 for e in P**4))
sweep=surface+tgam*S.Matrix([1,*ps])
check("rank3 boundary whole sweep Jacobian",reduce(sweep.jacobian([tgam,*ps]).det(method="domain-ge")-tgam**4)==0)
minor=surface.jacobian(ps).extract([0,3,4],[0,1,2]).det()
check("rank3 boundary is genuine and excludes full rank",
 reduce(minor.subs({u:1,x:1,y:1,z:1}))!=0 and all(zero(S.diff(e,z)) for e in surface))
quad=S.Matrix([sum(co*S.prod(p**e for p,e in zip([x,y,z],mon)) for mon,co in S.Poly(entry,x,y,z).terms() if sum(mon)==2) for entry in vv])
coefficients=S.Matrix([[S.Poly(entry,x,y,z).coeff_monomial(mon) for entry in quad] for mon in [(2,0,0),(1,1,0)]])
check("rank3 boundary has two independent quadratic outputs",
 reduce(vv.jacobian([x,y,z]).det())!=0 and coefficients.rank()==2)
# u factors disappear on the chart and must not be declared excluded.
cubic=x**3+u*u*y;allowed_form=u**2*cubic
check("u-factor boundary is invisible to dehomogenization",
 zero(allowed_form.subs(u,1)-(x**3+y)) and allowed_form.subs(u,0)==0)
check("necessary cubic factor is irreducible but not a realization",
 S.Poly(cubic,y).degree()==1 and S.gcd(u*u,x**3)==1 and S.total_degree(cubic)==3)
# Higher-degree factors have a genuine extra Taylor term.
qc=x**3+y+z;vt=S.Matrix([y+x*x,z,x]);ht=taylor_data(qc,vt,[x,y,z])
check("second-order Taylor argument does not automatically cover cubics",
 not zero(ht["taylor_identity"]) and zero(ht["taylor_identity"]-ht["W"][0]**3))
model={"claims":["KMFI136","KMFI137","KMFI138"],
 "previous_goal_turn":"progress: KMFI133--135 reduced the remaining quadratic direction class to two shapes, closed a coupled companion family and excluded completely split multipliers.",
 "quadratic_factor_invariance":{"ranks_and_direction_degrees":"all","projective_v_polynomial_automorphism":True,
 "full_rank_needed_for_invariance_lemma":False,"uses_feedback_residues":[0,1],
 "conclusion":"every irreducible factor of Fbar of degree1 or2 satisfies q(v)=constant*q"},
 "quadratic_factor_obstruction":{"architecture":"rank4 quadratic,V0=u^2,nonaffine projective direction","q2_ranks_exhausted":[1,2,3],
 "rank3":"dense quadratic image plane cannot be isotropic for nondegenerate q2",
 "rank1":"square displacement divisible by q forces a fixed affine direction",
 "rank2":"degree of product forces two coordinate images affine; residue forces a fixed affine direction",
 "all_total_scalar_degrees":True},
 "linear_factor_obstruction":{"non_u_factors_only":True,"all_other_factors_may_be_nonlinear":True,
 "arbitrary_M_and_affine_constants_retained":True,"d0_branch":"constant-diagonal triangular",
 "square_nonzero_d":"negative constant-derivative eigenvalue","mixed_nonzero_d":"negative Euler eigenvalue",
 "u_factor_excluded":False},
 "necessary_remaining_scalar_form":"F=constant*u^e*product D_i^bi, e>=0,bi>0,eachirreducibledegDi>=3,atleastoneDi",
 "actual_rank_deficient_boundary":{"coefficient_relation":"alpha^2-alpha+1=0","F":"x","v":[str(e) for e in vv],
 "rankDX":3,"whole_Jacobian":"tgam^4","new_Keller_action_claimed":False},
 "necessary_factor_form_control_claimed_realization":False,
 "arbitrary_higher_irreducible_factors_classified":False,"new_abstract_action_claimed":False,
 "actual_PSL11_input_supplied":False,"absolute_exceptional_action_decided":False,
 "full_action_spectrum_complete":False,"goal_complete":False}
(ROOT/"low_degree_multiplier_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["LOW_DEGREE_MULTIPLIER_DIVISORS.md","low_degree_multiplier.py","verify_low_degree_multiplier.py",
 "low_degree_multiplier_model.json","TWO_IMAGE_PROJECTIVE_REDUCTION.md","two_image_projective.py",
 "RANK_ONE_PROJECTIVE_CURVATURE.md","TRIANGULAR_PROJECTIVE_RIGIDITY.md","PRIMITIVE_MULTIPLIER_CRITERION.md","rank_one_projective.py"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,"runtime_seconds":round(time.monotonic()-started,3),
 "python":platform.python_version(),"sympy":S.__version__,"sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 "assurance":"Unbounded low-degree divisor invariance, exhaustive quadratic rank cases and generic hyperplane-flag negative cofactor proofs. Whole polynomial identities and an actual normalized rank3 boundary supplement written arguments. No formal verification, independent whole-proof replay or external review claimed.","goal_complete":False}
(ROOT/"low_degree_multiplier_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],
 "all_non_u_multiplier_factors_degree_at_least":3,"u_factor_not_excluded":True,"goal_complete":False}))
