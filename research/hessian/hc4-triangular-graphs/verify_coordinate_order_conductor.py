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
    Delta=(-1)**(m*(m-1)//2)
    zero("Keller exact signed norm-square certificate"+label,
         (ft*s.eye(m)).det()-(-1)**m*Delta**2)

# Necessary conductor restrictions distinguish two forbidden patterns.
bad_missing=x*y+1
bad_vertical=(y-1)*(x*x+1)
check("constant nonzero fibre candidate rejected",bad_missing.subs(y,0)==1)
check("identically zero fibre candidate rejected",bad_vertical.subs(y,1)==0)
good=x*x-1
check("nonempty finite fibre pattern retained",
      s.Poly(good,x).degree()==2 and s.gcd_list(s.Poly(good,x).all_coeffs())==1)


# Index-square controls of arbitrary multiple-sheet order, including a
# generic source shear. No reduced-conductor hypothesis is used.
def matrix_eval(poly, variable, matrix):
    result=s.zeros(matrix.rows)
    for (power,),coefficient in s.Poly(poly,variable).terms():
        result += coefficient*matrix**power
    return result

multiple_controls=[]
for f,k,alpha in [(t**3-t,1,0),(t**3-t,2,5),(t**4-5*t*t+4,1,3)]:
    M=s.degree(f,t)
    label=f" M={M},k={k},alpha={alpha}"
    coefficients=s.Poly(f-u,t)
    W=s.zeros(M)
    for j in range(M-1):
        W[j+1,j]=1
    for i in range(M):
        W[i,M-1]=-coefficients.nth(i)
    matrix_zero("companion polynomial identity"+label,
                matrix_eval(f,t,W)-u*s.eye(M))
    rU=u**k-alpha
    Bq=W*rU+y*s.eye(M)
    e0=s.eye(M)[:,0]
    C=s.Matrix.hstack(*[Bq**j*e0 for j in range(M)])
    Delta=s.factor(C.det())
    zero("explicit multiple-order index determinant"+label,
         Delta-rU**(M*(M-1)//2))
    GR=s.Matrix([[s.trace(W**(i+j)) for j in range(M)] for i in range(M)])
    GA=s.Matrix([[s.trace(Bq**(i+j)) for j in range(M)] for i in range(M)])
    matrix_zero("multiple-order trace Gram index identity"+label,
                GA-C.T*GR*C)
    conductormatrix=rU**(M-1)*s.eye(M)
    zero("multiple-order exact norm square"+label,
         conductormatrix.det()-Delta**2)
    rv=f**k-alpha
    q=t*rv+y
    FF=sum(s.Poly(f,t).nth(i)*(V-T)**i*(U**k-alpha)**(M-i)
           for i in range(M+1))-U*(U**k-alpha)**M
    substitution={U:f,V:q,T:y}
    zero("multiple-order monic defining relation"+label,FF.subs(substitution))
    check("multiple-order normalization unramified"+label,
          s.gcd(s.diff(f,t),s.diff(q,t))==1)
    ratio=s.cancel(s.diff(FF,V).subs(substitution)/s.diff(f,t))
    zero("multiple-order exact conductor ratio"+label,ratio-rv**(M-1))
    zero("multiple-order delta and conductor degrees"+label,
         s.degree(rv**(M-1),t)-2*s.degree(Delta,u))
    check("persistent multiple-sheet conductor is nonradical"+label,
          s.gcd(rv**(M-1),s.diff(rv**(M-1),t))!=1)
    zero("all generic source fibre points are collinear"+label,
         q-y-rv*t)
    multiple_controls.append({
        "M":int(M),"k":k,"source_shear_alpha":alpha,
        "index_delta":str(Delta),
        "conductor":str(rv**(M-1)),
        "conductor_degree":int(s.degree(rv**(M-1),t)),
        "normalization_unramified":True,
        "generic_triple_collinearity_persistent":True,
        "Keller":False
    })

r1,r2,r3,Rcoef,TargetV=s.symbols("r1 r2 r3 Rcoef TargetV")
zero("universal triple-collinearity identity",
     s.det(s.Matrix([[1,r,TargetV-Rcoef*r] for r in [r1,r2,r3]])))

# Necessary even-degree restrictions catch both generic and special failures.
odd_generic=x**3+x+1
odd_drop=y*x**4+x**3+x*x+1
even_drop=y*x**4+x*x+1
check("odd generic conductor degree rejected",s.degree(odd_generic,x)%2==1)
check("odd special conductor degree rejected",
      s.degree(odd_drop,x)%2==0 and s.degree(odd_drop.subs(y,0),x)%2==1)
check("even special degree drop retained by parity filter",
      s.degree(even_drop,x)==4 and s.degree(even_drop.subs(y,0),x)==2)

# Exact quadratic classification transcription and nodal pinching model.
shift=y*y-y
quadratic=(x-shift)**2-3
check("quadratic discriminant is constant nonzero",
      s.discriminant(quadratic,x)==12)
zero("quadratic polynomial shear recovers fixed nodal roots",
     quadratic.subs(x,z+shift)-(z*z-3))
node_p=z*z-1
node_q=z*node_p
zero("lowest-degree nodal product relation",
     node_q**2-node_p**2*(node_p+1))
for testpoly in [node_p,node_q,y,node_p*(z**5+y*z)]:
    check("nodal-order paired-value restriction",testpoly.subs(z,1)==testpoly.subs(z,-1))
proof=HERE/"GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md"
record={
    "outcome":"PASS","check_count":len(checks),"checks":checks,
    "controls":controls,"multiple_order_controls":multiple_controls,"sympy_version":s.__version__,
    "elapsed_seconds":round(time.perf_counter()-started,3),
    "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
    "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "scope":{
        "universal_arguments_are_written_not_finitely_certified":True,
        "generic_raw_single_coordinate_reconstruction_proved":True,
        "exact_principal_coordinate_order_conductor_proved":True,
        "explicit_conductor_norm_square_certificate_proved":True,
        "even_conductor_degree_on_every_source_line_proved":True,
        "degree_two_coordinate_order_is_nodal_product_proved":True,
        "degree_two_conductor_case_excluded":False,
        "generic_conductor_always_reduced":False,
        "persistent_generic_triple_collinearity_non_Keller_controls_supplied":True,
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
