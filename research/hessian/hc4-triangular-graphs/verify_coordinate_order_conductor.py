"""Exact controls for generic coordinate orders and collision conductors."""
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
def matrix_zero(label,matrix):
    check(label,all(s.cancel(v)==0 for v in matrix))

t,y,u,U,V,T,z,h=s.symbols("t y u U V T z h")
controls=[]
for k in [1,2,3,5]:
    label=f" k={k}"
    p=t*t-1
    q=t*p**k+y
    F=(V-T)**2-(U+1)*U**(2*k)
    zero("coordinate-order defining relation"+label,F.subs({U:p,V:q,T:y}))
    jac=s.det(s.Matrix([p,q]).jacobian([t,y]))
    zero("non-Keller target Jacobian"+label,jac-2*t)
    check("unramified normalization derivative gcd"+label,
          s.gcd(s.diff(p,t),s.diff(q,t))==1)
    ft=s.diff(F,T).subs({U:p,V:q,T:y})
    fv=s.diff(F,V).subs({U:p,V:q,T:y})
    px=s.diff(p,t)
    zero("finite trace-dual conductor ratio"+label,fv/px-p**k)
    zero("nonconstant-Jacobian corrected F_T ratio"+label,-ft/jac-p**k)
    check("naive F_T collision locus rejected"+label,
          ft.subs(t,0)==0 and (p**k).subs(t,0)!=0)

    # Independent trace matrices over D=C[u,y], where t^2=u+1.
    W=s.Matrix([[0,u+1],[1,0]])
    Bq=u**k*W+y*s.eye(2)
    GR=s.Matrix([[s.trace(W**(i+j)) for j in range(2)] for i in range(2)])
    GA=s.Matrix([[s.trace(Bq**(i+j)) for j in range(2)] for i in range(2)])
    C=s.Matrix([[1,y],[0,u**k]])
    matrix_zero("order Gram change of lattice"+label,GA-C.T*GR*C)
    detR=s.factor(GR.det())
    detA=s.factor(GA.det())
    zero("trace discriminant index square"+label,detA-detR*u**(2*k))
    # R-dual is (2t)^-1 R; A-dual is (2t*u^k)^-1 A.
    DR=(2*W).inv()
    DA=(2*u**k*W).inv()*C
    matrix_zero("normalization monogenic trace dual"+label,
                GR*DR-s.Matrix([[0,1],[1,0]]))
    # Pair basis (1,q) with columns DA; it need only be unimodular.
    pairing=C.T*GR*DA
    check("order monogenic trace dual determinant unit"+label,
          s.cancel(pairing.det()) in [-1,1])

    # Direct coefficient check: A=D + D*u^k*t, conductor=u^k*R.
    for B in [s.eye(2),W]:
        prod=u**k*B
        check("conductor lattice lies in A"+label,
              all(s.cancel((C.inv()*prod)[i,j]).is_polynomial(u,y)
                  for i in range(2) for j in range(2)))
    rejected=s.cancel((C.inv()*(u**(k-1)*W))[1,0])
    check("smaller conductor power rejected"+label,not rejected.is_polynomial(u,y))
    check("conductor collision support exactly plus/minus one"+label,
          s.factor(p**k)==(t-1)**k*(t+1)**k)
    check("collision fibre is reduced"+label,s.gcd(p,s.diff(p,t))==1)
    check("distinct collision pair same target and coordinate"+label,
          q.subs(t,1)==q.subs(t,-1) and p.subs(t,1)==p.subs(t,-1))
    check("no vertical conductor component"+label,
          s.Poly(p**k,t).LC()==1)
    controls.append({"k":k,"normalization_unramified":True,
                     "conductor":str(p**k),
                     "conductor_radical":k==1,
                     "source_target_Jacobian":str(jac)})

# Actual constant-Jacobian automorphism family: nontrivial D-ranks.
x=s.symbols("x")
H=x+y*y+y
for m in [2,3,5]:
    label=f" m={m}"
    P=H**m+y
    Q=-H
    F=V**m+(-1)**m*(T-U)
    jac=s.det(s.Matrix([P,Q]).jacobian([x,y]))
    zero("Keller automorphism Jacobian"+label,jac-1)
    sub={U:P,V:Q,T:y}
    zero("Keller coordinate-order relation"+label,F.subs(sub))
    fu=s.diff(F,U).subs(sub)
    fv=s.diff(F,V).subs(sub)
    ft=s.diff(F,T).subs(sub)
    a=s.diff(P,x)
    b=s.diff(Q,x)
    zero("Keller conductor gradient U identity"+label,fu-b*ft)
    zero("Keller conductor gradient V identity"+label,fv+a*ft)
    zero("Keller conductor trace ratio"+label,fv/a+ft)
    check("Keller inverse conductor is unit"+label,ft in [-1,1])

# Necessary conductor restrictions distinguish two forbidden patterns.
bad_missing=x*y+1
bad_vertical=(y-1)*(x*x+1)
check("constant nonzero fibre candidate rejected",bad_missing.subs(y,0)==1)
check("identically zero fibre candidate rejected",bad_vertical.subs(y,1)==0)
good=x*x-1
check("nonempty finite fibre pattern retained",
      s.Poly(good,x).degree()==2 and s.gcd_list(s.Poly(good,x).all_coeffs())==1)

proof=HERE/"GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md"
record={
    "outcome":"PASS","check_count":len(checks),"checks":checks,
    "controls":controls,"sympy_version":s.__version__,
    "elapsed_seconds":round(time.perf_counter()-started,3),
    "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
    "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "scope":{
        "universal_arguments_are_written_not_finitely_certified":True,
        "generic_raw_single_coordinate_reconstruction_proved":True,
        "exact_principal_coordinate_order_conductor_proved":True,
        "entire_closed_same_line_collision_support_proved":True,
        "all_source_line_conductor_restrictions_proved_using_imported_one_line_theorem":True,
        "arbitrary_tangency_finite_unramified_controls_supplied":True,
        "coordinate_order_conductor_is_always_radical":False,
        "nonnormal_coordinate_order_identified_with_full_target_normalization":False,
        "first_velocities_automatically_integral_over_target":False,
        "automatic_trace_integrality_proved":False,
        "JC2_or_HC4_resolved":False,
        "independent_mathematical_review":False,
        "formal_verification":False,"external_review":False
    }
}
outdir=Path(tempfile.mkdtemp(prefix="coordinate-order-conductor-check-"))
receipt=outdir/"coordinate_order_conductor_verification.json"
receipt.write_text(json.dumps(record,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","checks":len(checks),"receipt":str(receipt)}))
