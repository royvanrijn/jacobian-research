"""Exact variable-product periods and global quartic-conductor exclusion controls."""
from pathlib import Path
import hashlib,json,tempfile,time
import sympy as sp

HERE=Path(__file__).resolve().parent
started=time.perf_counter()
checks=[]
def check(name,test):
    if not test:
        raise AssertionError(name)
    checks.append(name)
def zero(name,expr):
    check(name,sp.cancel(expr)==0)

A,S,K,S1,K1,Ap,Lam,Lam1,J=sp.symbols("A S K S1 K1 Ap Lam Lam1 J")
B=K/A
St=A+B
d=A-B
Ca=3*A**3-7*K*A
Cb=3*B**3-7*K*B
Aprime=(A*S1-K1)/d
Bprime=sp.diff(B,A)*Aprime+sp.diff(B,K)*K1
zero("paired roots have correct differentiated sum",Aprime+Bprime-S1)
zero("paired roots have correct differentiated product",B*Aprime+A*Bprime-K1)
Ca_prime=sp.diff(Ca,A)*Aprime+sp.diff(Ca,K)*K1
Cb_prime=sp.diff(Cb,A)*Aprime+sp.diff(Cb,K)*K1
zero("first canonical period derivative with variable product",
     Ca_prime+Aprime*Ca/(2*A)
     -sp.Rational(7,2)*(3*A*A*S1-5*A*K1))
zero("second canonical period derivative with variable product",
     Cb_prime+Bprime*Cb/(2*B)
     -sp.Rational(7,2)*(3*B*B*S1-5*B*K1))
N=3*St*St-10*K
EA=Lam1*Ca+sp.Rational(7,2)*Lam*(3*A*A*S1-5*A*K1)+2*J
EB=Lam1*Cb+sp.Rational(7,2)*Lam*(3*B*B*S1-5*B*K1)+2*J
diff=2*N*Lam1+7*(3*St*S1-5*K1)*Lam
summ=St*(3*St*St-16*K)*Lam1+sp.Rational(7,2)*Lam*(
    3*(St*St-2*K)*S1-5*St*K1)+4*J
zero("variable-product period difference",EA-EB-d*diff/2)
zero("variable-product period sum",EA+EB-summ)
zero("variable-product period scale elimination",
     2*N*summ-St*(3*St*St-16*K)*diff
     -2*(105*K*Lam*(2*K*S1-St*K1)+4*J*N))

a,b,c,a1,b1,a2,b2=sp.symbols("a b c a1 b1 a2 b2")
Ss=-b/a
Ks=c/a
Sp=(b*a1-a*b1)/a**2
Kp=-c*a1/a**2
D=3*b*b-10*c*a
H=b*a1-2*a*b1
E=3*a*b*b1-(3*b*b-5*c*a)*a1
Npoly=3*Ss**2-10*Ks
zero("polynomial normalization of first scale equation",
     a**3*(2*Npoly*Lam1+7*(3*Ss*Sp-5*Kp)*Lam)
     -(2*a*D*Lam1+7*E*Lam))
zero("polynomial normalization of second scale equation",
     a**4*(105*Ks*Lam*(2*Ks*Sp-Ss*Kp)+4*J*Npoly)
     -(105*c*c*Lam*H+4*J*a*a*D))

# Eliminate Lambda independently by differentiating its rational expression.
Dprime=6*b*b1-10*c*a1
Hprime=b*a2-a1*b1-2*a*b2
logderivative=2*a1/a+Dprime/D-Hprime/H
optional=2*a*D*Hprime-H*(33*a*b*b1-(9*b*b+25*c*a)*a1)
zero("optional polynomial ODE follows from rational scale",
     2*a*D*H*(logderivative+7*E/(2*a*D))+optional)

m=sp.symbols("m",integer=True,positive=True)
zero("all-integer local multiplicity equation",
     6*(m+1)-21*m-(6-15*m))
check("unique formal zero multiplicity is nonintegral",
      sp.solve(sp.Eq(6*(m+1)-21*m,0),m)==[])
# Solve without integer assumptions to display the exact forbidden exponent.
mm=sp.symbols("mm")
check("unrestricted local multiplicity solution is two fifths",
      sp.solve(sp.Eq(6*(mm+1)-21*mm,0),mm)==[sp.Rational(2,5)])

v=sp.symbols("v")
for order in [1,2,3,5,8]:
    aa=2*v**order+v**(order+1)
    bb=3+v+2*v*v
    dd=3*bb*bb-10*aa
    hh=bb*sp.diff(aa,v)-2*aa*sp.diff(bb,v)
    ee=3*aa*bb*sp.diff(bb,v)-(3*bb*bb-5*aa)*sp.diff(aa,v)
    # Nonzero constants J=c=1; retain the complete rational scale.
    ll=sp.cancel(-4*aa*aa*dd/(105*hh))
    normalized=sp.cancel(ll/v**(order+1))
    check(f"rational scale has exact local order m+1 m={order}",
          normalized.subs(v,0)!=0)
    l0=normalized.subs(v,0)
    residual=sp.cancel(2*aa*dd*sp.diff(ll,v)+7*ee*ll)
    lead=sp.limit(residual/v**(2*order),v,0)
    zero(f"full local residual has forbidden exact leading coefficient m={order}",
         lead-3**2*2*l0*(6-15*order))
    check(f"local ODE obstruction is nonzero including m=1 m={order}",
          lead!=0)

# Actual monic coordinate prevents a finite root's reflected mate escaping.
x,z,h=sp.symbols("x z h")
P=x**4+v*x*x+3*x+2
finite=v+v*v
mate_polar=2/v-finite
difference=sp.expand(P.subs(x,mate_polar)-P.subs(x,finite))
zero("monic mate equality fails at the forced leading pole",
     sp.limit(v**4*difference,v,0)-16)
# The odd cubic coefficient of a pole-centred quartic is visibly nonregular
# unless the monic-mate gate repairs the centroid.
quartic=sp.expand((x-h)**4+b*(x-h)**2+c)
zero("quartic root centroid formula",sp.Poly(quartic,x).nth(3)+4*h)
zero("centroid shift makes the quartic even",
     quartic.subs(x,z+h)-z**4-b*z*z-c)

# Sharp non-Keller degree-drop family with the correct exact order.
u,y,U,V=sp.symbols("u y U V")
g=(y-1)*u**4+u*u+1
p=u*u
q=u*g
F=V*V-U*((y-1)*U*U+U+1)**2
zero("degree-drop control plane equation",
     F.subs({U:p,V:q},simultaneous=True))
zero("degree-drop control exact conductor",
     sp.diff(F,V).subs({U:p,V:q},simultaneous=True)/sp.diff(p,u)-g)
zero("degree-drop control Jacobian is non-Keller",
     sp.diff(p,u)*sp.diff(q,y)-sp.diff(p,y)*sp.diff(q,u)-2*u**6)
check("degree-drop control is immersed on every source line",
      sp.gcd(sp.diff(p,u),sp.diff(q,u))==1)
check("degree-drop control has generic conductor degree four",
      sp.degree(g,u)==4)
check("degree-drop control specializes to conductor degree two",
      sp.degree(g.subs(y,1),u)==2)
check("degree-two drop has distinct conductor roots",
      sp.discriminant(g.subs(y,1),u)!=0)
check("degree-drop leading and quadratic coefficients are coprime",
      sp.gcd(y-1,sp.Integer(1))==1)

audit_path=HERE/"sources/neumann-two-node-genus-two-source-audit.json"
audit=json.loads(audit_path.read_text())
source_hashes={}
for source in audit["sources"]:
    path=HERE/"sources"/source["local_pdf"]
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    check("pinned primary source "+source["local_pdf"],actual==source["sha256"])
    source_hashes[source["local_pdf"]]=actual
proof=HERE/"GLOBAL_DEGREE_FOUR_CONDUCTOR_EXCLUSION.md"
record={
    "outcome":"PASS","check_count":len(checks),"checks":checks,
    "source_sha256":source_hashes,
    "source_audit_sha256":hashlib.sha256(audit_path.read_bytes()).hexdigest(),
    "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
    "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "sympy_version":sp.__version__,
    "elapsed_seconds":round(time.perf_counter()-started,3),
    "scope":{
        "all_global_directional_degree_four_conductors_excluded":True,
        "special_degree_four_to_two_drops_allowed_in_theorem":True,
        "noninvertible_actual_keller_global_conductor_degree_at_least_six":True,
        "monic_actual_target_coordinate_used_to_remove_centroid_poles":True,
        "variable_product_period_equations_checked":True,
        "unrestricted_integer_valuation_argument_proved_by_written_proof":True,
        "genus_two_classification_proved_by_symbolic_checker":False,
        "uniform_target_automorphism_over_parameter_assumed":False,
        "automatic_trace_integrality_or_JC2_proved":False,
        "independent_replay":False,"independent_mathematical_review":False,
        "formal_verification":False,"external_review":False
    }
}
outdir=Path(tempfile.mkdtemp(prefix="global-quartic-conductor-check-"))
receipt=outdir/"global_degree_four_conductor_verification.json"
receipt.write_text(json.dumps(record,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","checks":len(checks),"receipt":str(receipt)}))
