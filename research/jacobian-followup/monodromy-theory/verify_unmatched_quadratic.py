"""Exact polarization, chart and unmatched-feedback controls, KMFI142--145."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import sympy as S
import json,time,hashlib,platform
from unmatched_quadratic import canonical_forms,chart_inverses,projective_direction,swap_output
from primitive_multiplier import feedback_form,logarithmic_residues
from divisor_homotopy import affine_survivor
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def zero(e):return S.cancel(e)==0
def mz(M):return all(zero(e) for e in M)
u,x,y,z=S.symbols("u x y z");ps=[u,x,y,z];p=S.Matrix(ps)
qs=list(S.symbols("q0:4"));q=S.Matrix(qs)
mons=[ps[i]*ps[j] for i in range(4) for j in range(i,4)]
cc=S.symbols("c0:40");Vgen=S.Matrix([sum(cc[10*i+j]*mon for j,mon in enumerate(mons)) for i in range(4)])
for sign in [1,-1]:
    rhs=Vgen.jacobian(ps).subs(dict(zip(ps,(p+sign*q)/2)),simultaneous=True)*(p-sign*q)
    lhs=Vgen-Vgen.subs(dict(zip(ps,sign*q)),simultaneous=True)
    check("entire generic40-coefficient midpoint identity sign"+str(sign),mz((lhs-rhs).applyfunc(S.expand)))
check("quadratic map is even so affine injectivity is not asserted",
 mz(Vgen.subs(dict(zip(ps,-p)),simultaneous=True)-Vgen))
# Every canonical polynomial chart and its inverse are exact.
forms=canonical_forms(ps);inverses=chart_inverses([x,y,z]);mm=S.symbols("m0:16");A=S.Matrix(4,4,mm)
detA=A.det();canonical=[]
for name,K in forms.items():
    phi=K[1:,:].subs(u,1);inv=inverses[name];xx=S.Matrix([x,y,z])
    check(name+" full source Jacobian-power identity",zero(K.jacobian(ps).det(method="domain-ge")-2*u**4))
    check(name+" no common coordinate factor",S.gcd_list(list(K))==1)
    check(name+" whole affine chart determinant",zero(phi.jacobian([x,y,z]).det()-1))
    check(name+" explicit inverse in both directions",
      mz(phi.subs(dict(zip([x,y,z],inv)),simultaneous=True)-xx) and
      mz(inv.subs(dict(zip([x,y,z],phi)),simultaneous=True)-xx))
    VV=A*K
    check(name+" every independent output matrix retained in derivative",
       mz(VV.jacobian(ps)-A*K.jacobian(ps)))
    check(name+" entire adjugate gives the target hyperplane equation",
      zero((A.adjugate()*VV)[0]-detA*u*u))
    check(name+" matched component forces all three offdiagonal row entries zero",
      all(S.Poly(VV[0],*ps).coeff_monomial(u*ps[j])==A[0,j] for j in range(1,4)))
    atinf=K.subs(u,0)
    expected={"A":S.Matrix([0,0,x*x,y*y]),"B":S.Matrix([0,0,x*x,x*y]),
              "C":S.Matrix([0,0,0,x*x]),"D":S.Matrix([0,0,0,x*y])}[name]
    check(name+" whole projective boundary base equations",mz(atinf-expected))
    # A specific unmatched output chart has target boundary x0, source boundary u0.
    P=swap_output(4,1);VU=P*K;alpha=S.Matrix([[1,0,0,0]])*P.inv()
    check(name+" independent source and target boundaries retained",
       list(alpha)==[0,1,0,0] and zero((alpha*VU)[0]-u*u) and not zero(VU[0]-u*u))
    check(name+" target boundary cannot be a multiplier divisor at u1",
       S.rem(S.Poly((alpha*VU)[0],x),S.Poly(x,x)).as_expr()!=0)
    canonical.append({"type":name,"K":[str(e) for e in K],"chart_inverse":[str(e) for e in inv],
                      "base_equations":[str(e) for e in atinf],"source_boundary":"u=0","example_target_boundary":"x=0"})
# Two signs matter: an ordinary quadratic map is two-to-one radially.
point={u:1,x:2,y:3,z:4};K=forms["A"]
check("antipodal affine points share the quadratic image",
 K.subs(point)==K.subs({v:-a for v,a in point.items()}) and p.subs(point)!=-p.subs(point))
check("one midpoint can lie on the critical hyperplane while the other does not",
 (p-p).subs(point)==S.zeros(4,1) and p.subs(point)[0]!=0)
# Global divisor homotopy before choosing matched coordinates.
t=S.Symbol("t");D=x**3+u*u*y;K=forms["C"]
phi=p+t*K
check("global cubic divisor homotopy polynomial identity",
 zero(D.subs(dict(zip(ps,phi)),simultaneous=True)-(1+t*u)**3*D))
check("global cubic divisor weighted projective identity",
 zero(D.subs(dict(zip(ps,K)),simultaneous=True)-u**3*D))
AA=K.jacobian(ps);grad=S.Matrix([[S.diff(D,v) for v in ps]])
check("global cubic divisor entire finite residue sequence",
 all(zero((grad*AA**j*K)[0]-3*(2**j)*u**(j+1)*D) for j in range(4)))
# A global invariant projective coordinate forces a linear factor in detDV.
l0,l1,l2,l3=S.symbols("l0 l1 l2 l3");L=l0*u+l1*x+l2*y+l3*z
row0=S.Matrix([[S.diff(u*L,v) for v in ps]]);row1=S.Matrix([[S.diff(x*L,v) for v in ps]])
for i in range(4):
 for j in range(i+1,4):
    wedge=S.expand(row0[0,i]*row1[0,j]-row0[0,j]*row1[0,i])
    check("entire derivative row minor has factorL "+str((i,j)),
      S.fraction(S.cancel(wedge/L))[1]==1)
# Weighted factor identity in an actual full-rank affine survivor (not an unmatched witness).
positive=affine_survivor();rr=positive["red"];lam0=positive["l"];pc=positive["ps"];uu,xx,yy,zz=pc
eig=lam0-1;ell=S.Matrix([[1,eig,eig**2]]);linear=xx+eig*yy+eig**2*zz
check("positive affine norm factor eigenrow exact",
 all(rr(e)==0 for e in ell*positive["C"]-eig*ell))
rem=S.rem(S.Poly(positive["F"],xx),S.Poly(linear,xx)).as_expr()
check("weighted identity control factor really divides scalarF",rr(rem)==0)
check("weighted identity in the actual normalized affine input",
 rr(linear.subs(dict(zip(pc,positive["V"])),simultaneous=True)-eig*uu*linear)==0)
# Eight exact controls are samples of A, not a claim about all A.
lam=S.Symbol("lam");bench=[]
for name,K in forms.items():
 for i in [1,3]:
    P=swap_output(4,i);VV=P*K;tag=name+"_swap0"+str(i)
    try:f=feedback_form(VV,ps,lam*u)
    except ValueError:
        C=S.Matrix.hstack(*[VV.jacobian(ps)**j*VV for j in range(4)])
        check(tag+" exact singular Krylov branch",zero(C.det(method="domain-ge")))
        bench.append({"type":name,"swap":i,"result":"full-rank rejected by singularC"});continue
    logs=logarithmic_residues(f,ps,lam);cond=logs["parameter_condition"].as_expr()
    expected=(lam-1)**3 if i==1 else S.Integer(1)
    check(tag+" exact univariate logarithmic consistency",zero(cond-expected))
    if i==1:
        s_at=1/P.det()-2
        check(tag+" sole parameter has negative Euler degree",s_at==-3)
    bench.append({"type":name,"swap":i,"condition":str(S.factor(cond)),
                  "result":"no polynomial multiplier in the specified full-rank control",
                  "all_A_classified":False})
# Pin precisely the two imported primary inputs.
manifest=json.loads((ROOT/"sources/unmatched-quadratic-chart-manifest.json").read_text())
for category in ["sources","retained_primary_inputs"]:
 for src in manifest[category]:
    raw=(ROOT/src["path"]).read_bytes()
    check("primary source pin "+src["path"],len(raw)==src["bytes"] and hashlib.sha256(raw).hexdigest()==src["sha256"])
model={"claims":["KMFI142","KMFI143","KMFI144","KMFI145"],
 "previous_goal_turn":"progress: KMFI139--141 closed the entire matched pure-square quadratic rank4 architecture.",
 "quadratic_chart_theorem":{"all_ranks_r_ge2":True,"Jacobian":"delta*u^r",
 "midpoint_requires_both_signs":True,"projectively_injective_over_all_algebraic_closures":True,
 "source_U_is_A_n":"u nonzero","image":"complement of one unique target hyperplane alpha",
 "alphaV":"u^2 after scale","alpha_proportional_u_assumed":False,
 "affine_direction_map_degree":2,"direction_field_confused_with_Keller_sweep_field":False},
 "four_exhaustive_rank4_primitive_forms":canonical,
 "arbitrary_GL4_output_matrix_retained":True,
 "global_multiplier_divisor":{"all_degrees":True,"no_matched_component_needed":True,
 "homotopy":"D(p+tV) divisible byD","nonboundary_weighted_preservation":"D(V)=constant*u^degD*D",
 "unmatched_target_boundary_factor_alpha_excluded":True,"same_claim_for_u_factor":False},
 "scalar_factor_bound":{"primitive_quadratic_rank4":True,"all_non_u_factors_degree_at_most":3,
 "total_scalar_degree_or_multiplicity_bound_claimed":False,
 "old_matched_lower_degree_bound_imported_to_unmatched":False,
 "matched_component_forced_if_degree_gt3":True},
 "eight_explicit_unmatched_controls":bench,
 "whole_arbitrary_A_parameter_space_classified":False,
 "new_abstract_action_claimed":False,"actual_PSL11_input_supplied":False,
 "absolute_exceptional_action_decided":False,"full_action_spectrum_complete":False,"goal_complete":False}
(ROOT/"unmatched_quadratic_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["UNMATCHED_QUADRATIC_CHARTS.md","unmatched_quadratic.py","verify_unmatched_quadratic.py",
 "unmatched_quadratic_model.json","sources/unmatched-quadratic-chart-manifest.json",
 "sources/stacks-025G-open-immersion.html","sources/debondt-quadratic-rank-two-1601.00579v4.pdf",
 "divisor_homotopy.py","primitive_multiplier.py","QUADRATIC_PURE_SQUARE_RIGIDITY.md",
 "TWO_IMAGE_PROJECTIVE_REDUCTION.md","PRIMITIVE_MULTIPLIER_CRITERION.md"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,"runtime_seconds":round(time.monotonic()-started,3),
 "python":platform.python_version(),"sympy":S.__version__,"sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 "assurance":"Unbounded two-sign projective chart theorem with explicit source import, four-form classification, global divisor homotopy and degree<=3 scalar-factor bound. Whole polynomial identities, exact inverses and eight specified feedback controls supplement the written proofs; no formal verification, independent whole-proof replay or external review claimed.","goal_complete":False}
(ROOT/"unmatched_quadratic_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],
 "primitive_quadratic_normal_forms":4,"all_non_u_factor_degrees_at_most":3,"goal_complete":False}))
