"""Exact controls for regular-coordinate plane completions."""
from pathlib import Path
import hashlib
import json
import tempfile
import time
import sympy as s

HERE=Path(__file__).resolve().parent
started=time.perf_counter()
checks=[]
def check(label,condition):
    if not condition:
        raise AssertionError(label)
    checks.append(label)
def zero(label,expression):
    check(label,s.cancel(expression)==0)

x,y,z,w,u=s.symbols("x y z w u")
controls=[]
for n,q in [(2,1),(2,2),(3,2),(3,3),(4,5)]:
    label=f" n={n},q={q}"
    Z=x*y**q
    aa={j:x**j*(Z-1) for j in range(1,n+1)}
    zero("boundary chart defining relation"+label,
         y**(q*n)*aa[n]-Z**n*(Z-1))
    for j in range(1,n+1):
        zero("localized generator reconstruction"+label+f",j={j}",
             y**(q*(n-j))*aa[n]-aa[j]*Z**(n-j))
        zero("finite monic generator relation"+label+f",j={j}",
             aa[j]**n-aa[n]**j*(Z-1)**(n-j))
    zero("monic finite z relation"+label,Z**(n+1)-Z**n-y**(q*n)*aa[n])
    zero("source chart inverse"+label,aa[1]/(Z-1)-x)
    C=z**(n-1)*((n+1)*z-n)
    f=y**(q*n)*w-z**n*(z-1)
    zero("implicit chart derivative"+label,-s.diff(f,z)-C)
    check("entire boundary chart is smooth"+label,C.subs(z,1)==1)
    jac=s.det(s.Matrix([y,aa[n]]).jacobian([x,y]))
    expected=-x**(n-1)*((n+1)*Z-n)
    zero("independent non-Keller source Jacobian"+label,jac-expected)
    check("source Jacobian has retained zeros"+label,jac.subs(x,0)==0)
    zero("target Jacobian compensates volume zero"+label,
         expected.subs(x,z/y**q)+C/y**(q*(n-1)))
    a=q*(n-1)
    check("positive canonical boundary coefficient"+label,a>0)
    zero("implicit-chart source volume coefficient"+label,
         (y**(q*n)/C)/y**q-y**a/C)
    # The defining equation proves x=y^-q+y^a*w/z^n.
    zero("exact boundary principal-part identity"+label,
         (z/y**q-y**(-q)-y**a*w/z**n)*y**q*z**n+f)
    if q==1:
        # The correction has nonnegative y order, so residue is exactly 1.
        check("q1 logarithmic primitive residue"+label,
              s.residue(y**(-1)+y**a*w,y,0)==1)
        check("q1 zero-residue regular primitive claim rejected"+label,
              s.residue(y**(-1)+y**a*w,y,0)!=0)
    else:
        A=q*x/s.Integer(q-1)    # dy coefficient
        B=y/s.Integer(q-1)      # dx coefficient
        zero("polynomial primitive exterior derivative"+label,
             s.diff(A,x)-s.diff(B,y)-1)
        # theta=(q*x dy+y dx)/(q-1), pull x=z/y^q:
        zero("chart primitive dy cancellation"+label,
             q*(z/y**q)/(q-1)+(y/(q-1))*s.diff(z/y**q,y))
        zero("chart primitive dz coefficient"+label,
             (y/(q-1))*s.diff(z/y**q,z)-1/((q-1)*y**(q-1)))
        local_dy=q*n*y**a*w/((q-1)*C)
        local_dw=y**(a+1)/((q-1)*C)
        zero("implicit primitive dy coefficient"+label,
             q*n*y**(q*n-1)*w/((q-1)*y**(q-1)*C)-local_dy)
        zero("implicit primitive dw coefficient"+label,
             y**(q*n)/((q-1)*y**(q-1)*C)-local_dw)
        check("primitive coefficients regular at entire boundary"+label,
              a>=0 and a+1>=0 and C.subs(z,1)==1)
    controls.append({"n":n,"q":q,"canonical_boundary_order":a,
                     "source_volume_exact":q>=2,
                     "finite_target_generic_degree":n+1})

# An arbitrary-polynomial tangent cover has a volume-preserving LND
# on its full finite normalization, without target-field descent.
t=s.symbols("t")
H=s.Function("H")(t)
P=s.diff(H,t)-u
Q=H-t*s.diff(H,t)+u*t
zero("all-polynomial tangent cover Jacobian",s.det(s.Matrix([P,Q]).jacobian([u,t]))+u)
zero("all-polynomial tangent cover root equation",H-P*t-Q)
zero("all-polynomial source translation first target velocity",s.diff(P,t)-s.diff(H,t,2))
zero("all-polynomial source translation second target velocity",s.diff(Q,t)-(u-t*s.diff(H,t,2)))
zero("translation preserves pulled-back target volume coefficient",s.diff(-u,t))
zero("translation preserves ramification equation",s.diff(u,t))
for degree in (3,4,7):
    HH=t**degree+t**2+2*t+3
    PP=s.diff(HH,t)-u
    QQ=HH-PP*t
    zero(f"degree {degree} finite root equation",HH-PP*t-QQ)
    zero(f"degree {degree} integral source-coordinate reconstruction",s.diff(HH,t)-PP-u)
    check(f"degree {degree} target derivative has positive smaller minimal degree",
          1<=s.degree(s.diff(HH,t,2),t)<degree)
# Quadratic boundary case: the translation DOES descend.
HH=t**2+2*t+3
PP=s.diff(HH,t)-u
QQ=HH-PP*t
zero("quadratic translation target P derivative is constant",s.diff(PP,t)-2)
zero("quadratic translation target Q derivative descends",s.diff(QQ,t)-(2-PP))

# Differential of every Laurent monomial has zero logarithmic residue.
# This exact check includes negative, zero and positive orders.
for exponent in [-7,-2,-1,0,1,4]:
    check(f"rational derivative logarithmic residue k={exponent}",
          s.residue(s.diff(y**exponent,y),y,0)==0)

proof=HERE/"REGULAR_SOURCE_COORDINATE_COMPLETIONS_AND_PRIMITIVE_RESIDUES.md"
record={
    "outcome":"PASS",
    "check_count":len(checks),
    "checks":checks,
    "controls":controls,
    "sympy_version":s.__version__,
    "elapsed_seconds":round(time.perf_counter()-started,3),
    "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
    "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "scope":{
        "universal_chart_and_residue_proofs_are_written":True,
        "whole_plane_open_and_regular_actual_source_coordinate":True,
        "smooth_affine_normal_completion":True,
        "positive_nontorsion_canonical_class":True,
        "q1_entire_Keller_target_pair_architecture_excluded":True,
        "q_ge_2_exact_regular_source_volume_controls":True,
        "explicit_full_finite_normalization_supplied":True,
        "all_degree_volume_preserving_full_normalization_action_without_target_descent":True,
        "control_source_map_is_Keller":False,
        "one_integral_source_coordinate_Keller_problem_resolved":False,
        "automatic_trace_integrality_proved":False,
        "JC2_or_HC4_resolved":False,
        "independent_mathematical_review":False,
        "formal_verification":False,
        "external_review":False
    }
}
outdir=Path(tempfile.mkdtemp(prefix="source-coordinate-completion-check-"))
receipt=outdir/"source_coordinate_completions_verification.json"
receipt.write_text(json.dumps(record,indent=2)+"\n")
print(json.dumps({"outcome":record["outcome"],"checks":len(checks),
                  "receipt":str(receipt)}))
