"""Exact formal, line-family and whole-moment controls for KMFI139--141."""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import sympy as S
import json,time,hashlib,platform
from divisor_homotopy import truncate,formal_inverse,residue_moments,affine_survivor
from rank_one_projective import homogenize
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def zero(e):return S.cancel(e)==0
def mz(M):return all(zero(e) for e in M)
x,y,z,t,s,a=S.symbols("x y z t s a");xs=[x,y,z];xx=S.Matrix(xs);N=6
# A formal inverse is used in R[[t]], including when t1 is not an automorphism.
W=S.Matrix([x*x,0,0]);v=xx+W;G=formal_inverse(W,xs,t,N)
Phi=xx+t*W
left=Phi.subs(dict(zip(xs,G)),simultaneous=True)-xx
right=G.subs(dict(zip(xs,Phi)),simultaneous=True)-xx
check("formal polynomial-coefficient inverse in both directions through order6",
 all(truncate(e,t,N)==0 for e in left) and all(truncate(e,t,N)==0 for e in right))
check("formal inverse is not asserted to be a global polynomial inverse",not mz(left))
invjac=S.invert(S.Poly(1+2*t*x,t),S.Poly(t**N,t)).as_expr()
Y=S.Matrix([truncate(invjac*x*x,t,N),0,0])
velocity=G.diff(t)+Y.subs(dict(zip(xs,G)),simultaneous=True)
check("whole formal inverse velocity identity through order5",
 all(truncate(e,t,N-1)==0 for e in velocity))
q=y**3+z
mom=residue_moments(q,v,xs,5)
check("nonautomorphic formal example has every checked tangent coefficient",all(zero(e) for e in mom))
check("formal ideal remains preserved for the nonlinear divisor",
 zero(q.subs(dict(zip(xs,Phi)),simultaneous=True)-q) and
 truncate(q.subs(dict(zip(xs,G)),simultaneous=True)-q,t,N)==0)
# A nontrivial unit, not equality q(Phi)=q, occurs in a scaling example.
vr=a*xx;Wr=vr-xx;qr=x**3+y**3+z**3
check("scaling example entire homotopy divisor factor",
 zero(qr.subs(dict(zip(xs,xx+t*Wr)),simultaneous=True)-(1+(a-1)*t)**3*qr))
check("scaling example whole residue sequence",all(zero(e-3*(a-1)*a**j*qr) for j,e in enumerate(residue_moments(qr,vr,xs,5))))
# A cubic passes the first two moments but fails the third.
qc=z-x*y*y;vc=S.Matrix([2*x,-y,-2*z]);wc=vc-xx;mc=residue_moments(qc,vc,xs,3)
rem=lambda e:S.rem(S.Poly(S.expand(e),z),S.Poly(qc,z)).as_expr()
check("cubic countercontrol first two residues vanish",zero(rem(mc[0])) and zero(rem(mc[1])))
check("cubic countercontrol third residue is essential",not zero(rem(mc[2])))
check("cubic countercontrol homotopy has the exact surviving Taylor remainder",
 zero(rem(qc.subs(dict(zip(xs,xx+t*wc)),simultaneous=True))+4*t**3*x*y*y))
# A genuine ruled cubic and quadratic automorphism illustrate affine line restriction.
vl=S.Matrix([x+1,y,z+y*y]);wl=vl-xx
check("ruled cubic entire homotopy containment",
 zero(qc.subs(dict(zip(xs,xx+t*wl)),simultaneous=True)-qc))
check("ruled cubic whole moment tangency",all(zero(e) for e in residue_moments(qc,vl,xs,5)))
check("ambient map has a whole polynomial inverse",
 mz(vl.subs({x:x-1,y:y,z:z-y*y},simultaneous=True)-xx) and
 mz(S.Matrix([x-1,y,z-y*y]).subs(dict(zip(xs,vl)),simultaneous=True)-xx))
check("ambient map restricts affinely to every selected generator",
 mz(vl.subs(dict(zip(xs,xx+s*wl)),simultaneous=True)-(xx+(s+1)*wl)))
check("leading quadratic part vanishes on the generator direction",zero(wl[1]**2))
check("large divisor example fails cyclicity as required",zero(S.Matrix.hstack(wl,vl.jacobian(xs)*wl,vl.jacobian(xs)**2*wl).det()))
# A plane has a two-dimensional line family; homotopy does not force line-affine action.
vp=S.Matrix([x,z,y+z*z]);wp=vp-xx
check("plane boundary homotopy remains in the divisor",zero((xx+t*wp)[0]-x))
check("plane boundary assigned lines can have nonlinear images",not zero(wp[2]**2))
invp=S.Matrix([x,z-y*y,y])
check("plane boundary is still an ambient automorphism",
 mz(vp.subs(dict(zip(xs,invp)),simultaneous=True)-xx) and
 mz(invp.subs(dict(zip(xs,vp)),simultaneous=True)-xx))
# The pairwise-incidence lemma has both possible elementary configurations.
P0=S.Matrix([0,0,0,1])
concurrent=[S.Matrix.hstack(P0,S.Matrix([1,i,i*i,0])) for i in [0,1,2]]
check("concurrent projective lines span a nonplane family and share their vertex",
 all(L.rank()==2 for L in concurrent) and S.Matrix.hstack(*concurrent).rank()==4 and
 all(S.Matrix.hstack(concurrent[0],L).rank()==3 for L in concurrent[1:]))
# Both remaining quadratic leading parts have the same hyperplane base constraint.
mm=S.symbols("m0:9");M=S.Matrix(3,3,mm);detM=M.det();w0,w1,w2=S.symbols("w0 w1 w2")
for name,H in [("square",S.Matrix([0,w0*w0,w1*w1])),("mixed",S.Matrix([0,w0*w0,w0*w1]))]:
    Q=M*H
    check(name+" entire adjugate proves hyperplane base-locus constraint",
       mz((M.adjugate()*Q-detM*H).applyfunc(S.expand)) and zero((M.adjugate()*Q)[1]-detM*w0*w0))
# A fixed affine direction kills every projective Krylov determinant.
vf=S.Matrix([x,y+x*x,z+x*y]);wf=vf-xx;Bf=vf.jacobian(xs)
check("fixed affine covector kills the complete checked span",
 all(zero((Bf**j*wf)[0]) for j in range(4)) and zero(S.Matrix.hstack(wf,Bf*wf,Bf**2*wf).det()))
# An actual affine full-rank rank4 survivor lies on the primitive-degree1 boundary.
pos=affine_survivor();ps=pos["ps"];u,x,y,z=ps;l=pos["l"];T=pos["T"];red=pos["red"]
C=pos["C"];Mt=pos["M"];F=pos["F"];GG=pos["G"];point={u:1,x:1,y:0,z:0}
chi=(T-1)*C.charpoly(T).as_expr()
check("positive control quartic coefficient polynomial irreducible",pos["modulus"].is_irreducible)
check("positive control norm is exactly cubic",S.Poly(F,x,y,z).total_degree()==3)
check("positive control entire constant spectral feedback identity",red(5*chi-T*S.diff(chi,T)-(T-l)**4)==0)
grad=S.Matrix([[S.diff(GG,p) for p in ps]])
Ppower=S.eye(4);Cpower=S.eye(3)
for j in range(1,5):
    Ppower=(Ppower*Mt).applyfunc(red);Cpower=(Cpower*C).applyfunc(red)
    identity=(grad*Ppower*S.Matrix(ps))[0]-(1+S.trace(Cpower))*GG
    check("positive control whole norm moment"+str(j),red(identity)==0)
check("positive control full rank by exact unit minor",
 red(pos["X"][1:,:].jacobian(ps).subs(point).det(method="domain-ge")-l**4)==0 and red(l**4)!=0)
check("positive control displayed quadratic tuple has common u",
 zero(S.gcd_list(list(pos["V"]))-u))
# Nonlinear rank-deficient normalization remains possible.
alpha,tgam=S.symbols("alpha tgam");mod=S.Poly(alpha**2-alpha+1,alpha);lam=1+alpha
def ra(e):
    num,den=S.fraction(S.cancel(e))
    return S.cancel(S.rem(S.Poly(num,alpha),mod).as_expr()/S.rem(S.Poly(den,alpha),mod).as_expr())
vv=S.Matrix([alpha*x,lam*y+x*x,lam*z+x*y]);V=homogenize(vv,u,[x,y,z],2)
Fs=x;h=lam*u;X=S.Matrix([-Fs*h,*list(Fs*(V-h*S.Matrix(ps)))])
sweep=X+tgam*S.Matrix([1,*ps])
check("actual nonlinear rank3 boundary whole sweep Jacobian",
 ra(sweep.jacobian([tgam,*ps]).det(method="domain-ge")-tgam**4)==0)
check("actual nonlinear rank3 boundary qualification",
 all(zero(S.diff(e,z)) for e in X) and
 ra(X.jacobian(ps).extract([0,3,4],[0,1,2]).det().subs({p:1 for p in ps}))!=0)
model={"claims":["KMFI139","KMFI140","KMFI141"],
 "previous_goal_turn":"progress: KMFI136--138 excluded individual quadratic and non-u linear scalar factors, requiring a higher-degree divisor.",
 "all_degree_divisor_homotopy":{"all_ranks_and_direction_degrees":True,"full_rank_needed":False,"automorphy_needed_for_line_containment":False,"automorphy_needed_for_eigen_polynomial":True,
 "moment_hypotheses":"dq*(Dv-I)^j*(v-x) divisible byq, finite initial Krylov sequence suffices",
 "conclusion":"q(x+t*(v-x)) divisible byq in k[x,t] for every factor",
 "formal_inverse_not_global_polynomial_inverse":True,"formal_inverse_checkpoint":[str(e) for e in G],
 "cubic_first_two_residues_insufficient_countercontrol":True},
 "large_divisor_degree_bound":{"rank":4,"projective_variables":3,"v_degree":"D>=2",
 "leading_base_locus":"contained in a finite union of fixed linear hyperplanes",
 "full_rank_required":True,"conclusion":"every irreducible Fbar factor degree<=D",
 "selected_line_family_argument":"irreducible covering family on a nonplane surface is a curve and generically has one selected line",
 "ambient_automorphism_line_restriction":"closed affine line image equals itself, hence affine action",
 "plane_boundary_retained":True,"ruled_cubic_without_cyclicity_boundary_retained":True},
 "complete_quadratic_architecture":{"V0":"u^2 matched to sourceu","all_M_b_and_all_scalar_degrees":True,
 "conclusion":"v affine, V=u*M*p; primitive degree1; all nonlinear primitive quadratic members excluded",
 "degree11_groups":["S11","A11","regular C11"],"other_unmatched_quadratic_directions_classified":False},
 "positive_affine_full_rank_control":{"coefficient_minimal_polynomial":str(pos["modulus"].as_expr()),
 "C":[[str(e) for e in C.row(i)] for i in range(3)],"F":str(F),"G":"u*F","lambda":"l",
 "whole_normalization":"exact spectral characteristic and four whole norm moments, by KMFI121",
 "full_rank_unit_minor":"l^4","expanded_full_sweep_Jacobian_claimed":False,"primitive_coordinate_degree":1,
 "new_Keller_action_claimed":False},
 "actual_nonlinear_rank3_control_retained":True,
 "new_abstract_action_claimed":False,"actual_PSL11_input_supplied":False,
 "absolute_exceptional_action_decided":False,"full_action_spectrum_complete":False,"goal_complete":False}
(ROOT/"divisor_homotopy_model.json").write_text(json.dumps(model,indent=2)+"\n")
paths=["QUADRATIC_PURE_SQUARE_RIGIDITY.md","divisor_homotopy.py","verify_divisor_homotopy.py","divisor_homotopy_model.json",
 "LOW_DEGREE_MULTIPLIER_DIVISORS.md","TWO_IMAGE_PROJECTIVE_REDUCTION.md","TRIANGULAR_PROJECTIVE_RIGIDITY.md",
 "FACTORIZED_HOMOGENEOUS_FIELDS.md","PRIMITIVE_MULTIPLIER_CRITERION.md","rank_one_projective.py"]
receipt={"result":"pass","check_count":len(checks),"checks":checks,"runtime_seconds":round(time.monotonic()-started,3),
 "python":platform.python_version(),"sympy":S.__version__,"sha256":{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 "assurance":"Unbounded formal divisor-homotopy proof, elementary selected-line-family argument, leading-base-locus degree bound and full matched pure-square rank4 quadratic rigidity. Exact identities and affine-full-rank/nonlinear-rank3 controls supplement written proofs. No formal verification, independent whole-proof replay or external review claimed.","goal_complete":False}
(ROOT/"divisor_homotopy_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"runtime_seconds":receipt["runtime_seconds"],
 "entire_matched_pure_square_quadratic_rank4_architecture_closed":True,"goal_complete":False}))
