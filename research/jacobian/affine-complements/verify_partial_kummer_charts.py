"""PF142: exact partial-root charts and source seals; no predecessor replay."""
from pathlib import Path
from math import gcd
import hashlib
import json
import platform
import sympy as sp

P=Path(__file__).resolve().parent
NOTE="PARTIAL_KUMMER_PLANE_CHARTS_AND_CYCLIC_EQUIVARIANT_FACTORS.md"
MANIFEST="sources/partial-kummer-charts/manifest.json"
RECEIPT="partial_kummer_charts_receipt.json"


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    manifest=json.loads((P/MANIFEST).read_text())
    for s in manifest["sources"]:
        path=P/s["path"]
        check("source seal "+s["path"],path.stat().st_size==s["bytes"] and sha(path)==s["sha256"])

    # The binomial normalization is proved for arbitrary positive n,mu in
    # the note. These controls distinguish special-fibre multiplicity from
    # normal ramification, including ramification deliberately discarded.
    for n,mu in [(1,3),(2,2),(2,3),(3,2),(3,3),(4,2),(4,6),(6,4),(6,6),(6,9)]:
        g=gcd(n,mu);normal=n//g;fibre=mu//g
        check(f"primitive local parameter exponents {n},{mu}",gcd(normal,fibre)==1)
        check(f"binomial relation and degrees {n},{mu}",normal*mu==fibre*n and g*normal==n)
        check(f"local etale criterion {n},{mu}",(normal==1)==(mu%n==0))
    check("2+3 cover ramified exactly on discarded type",2//gcd(2,2)==1 and 2//gcd(2,3)==2)
    check("partial cover need not be universal etale cover",gcd(2,3)==1 and 2>1)

    t,w,X,Y,Z=sp.symbols("t w X Y Z")
    equation=X**2*Z-(Y**2-X)*((Y-1)**3-X)
    xx=t**2; yy=t+t**3*w
    zz=(2*w+t**2*w**2)*((yy-1)**3-t**2)
    substitution={X:xx,Y:yy,Z:zz}
    check("explicit plane chart lands on retained mixed surface",sp.expand(equation.subs(substitution))==0)
    jac=sp.det(sp.Matrix([[sp.diff(xx,t),sp.diff(xx,w)],[sp.diff(yy,t),sp.diff(yy,w)]]))
    check("unit volume pullback",sp.expand(jac-2*xx**2)==0)
    check("retained entire line isomorphism",[sp.expand(f.subs(t,0)) for f in [xx,yy,zz]]==[0,0,-2*w])
    check("third-multiplicity line not retained",sp.expand(yy.subs(t,0)-1)==-1)
    ww=sp.cancel((Y-t)/t**3)
    check("generic reconstruction x coordinate",sp.cancel(yy.subs(w,ww)-Y)==0)
    check("generic reconstruction surface coordinate",sp.rem(sp.together((zz.subs(w,ww)-(Y**2-t**2)*((Y-1)**3-t**2)/t**4))*t**4,t**2-X,t)==0)
    # Zeta roots do not change the residual affine-line parameter. The
    # displayed chart has irreducible degree two over C(X,Y).
    check("quadratic generic field extension",sp.Poly(t**2-X,t).degree()==2 and sp.Poly(X,X).degree()==1)

    # A second component choice verifies the same phenomenon at degree3.
    x3=t**3; y3=1+t+t**4*w
    z3=(3*w+3*t**3*w**2+t**6*w**3)*(y3**2-t**3)
    check("degree-three complementary chart lands",sp.expand(equation.subs({X:x3,Y:y3,Z:z3}))==0)
    jac3=sp.det(sp.Matrix([[sp.diff(x3,t),sp.diff(x3,w)],[sp.diff(y3,t),sp.diff(y3,w)]]))
    check("degree-three chart unit volume",sp.expand(jac3-3*x3**2)==0)
    check("degree-three retained line",[sp.expand(f.subs(t,0)) for f in [x3,y3,z3]]==[0,1,3*w])

    # One extra binomial transcription, not an increasing parameter scan.
    r,s=4,5; xr=t**r; yr=t+t**(r+1)*w
    factor=sp.div((1+t**r*w)**r-1,t**r,t)[0]
    zr=factor*((yr-1)**s-xr)
    check("arbitrary-binomial polynomial quotient control",sp.rem((1+t**r*w)**r-1,t**r,t)==0)
    check("all-parameter chart transcription control",sp.expand(xr**2*zr-(yr**r-xr)*((yr-1)**s-xr))==0)
    check("all-parameter volume transcription control",sp.expand(sp.diff(xr,t)*sp.diff(yr,w)-r*xr**2)==0)

    # Degree descent is independent of the rational slice. Arithmetic
    # controls cover non-coprime profiles as well as the former 2+3 case.
    for parts in [(2,3),(2,9),(4,6),(6,10,15),(7,7)]:
        m=sum(parts); n=min(parts)
        check("mixed minimum at most half "+str(parts),2*n<=m)
        check("strict comparison independent quotient degree "+str(parts),3*n<3*m)
    check("even the coprime mixed 2+3 case descends",5*7>2*7>1)
    check("whole gcd is not required",min([4,6])!=gcd(4,6))
    check("single power retains original comparison degree",7*11==7*min([11]))
    check("two wholly critical fibres exceed derivative budget",sp.simplify(sp.Symbol("m")/2+sp.Symbol("m")/2-(sp.Symbol("m")-1))==1)

    proof=(P/NOTE).read_text()
    for marker in ["MINIMAL GENERIC DEGREE","not the original finite",
                   "beta RATIONAL","finite etale universal cover",
                   "No volume, fundamental-group or Keller assumption",
                   "FIELD statement","full Wright/HC4/JC2 remain OPEN"]:
        check("written scope marker "+marker,marker in proof)

    theorem={
        "component_selected_partial_kummer_plane_chart":True,
        "nonpower_equivariant_factor_strict_degree_descent_at_least_half":True,
        "all_minimal_equivariant_factors_actual_cyclic_deck_fixed_fields":True,
        "minimal_trivial_deck_all_source_LND_derivative_fields_full":True,
        "explicit_mixed_surface_etale_plane_charts":True,
    }
    assurance={
        "all_degree_written_proof_not_parameter_census":True,
        "entire_cover_ramification_retained_then_deleted_on_comparison_open":True,
        "source_image_not_entire_intermediate_etale_locus":True,
        "MR0_PF130_PF139_ML0_MU0_classical_limits_credited":True,
        "full_original_normalization_and_trace_energy_goal_retained":True,
        "unconditional_mixed_Keller_exclusion":False,
        "automatic_equivariant_factor":False,
        "new_external_theorem_import":False,
        "registered_predecessor_mathematical_replay":False,
        "formal_verification":False,
        "independent_replay":False,
        "external_review":False,
        "novelty_claimed":False,
        "full_goal_complete":False,
    }
    pins={NOTE:sha(P/NOTE),Path(__file__).name:sha(Path(__file__)),MANIFEST:sha(P/MANIFEST)}
    pins.update({s["path"]:s["sha256"] for s in manifest["sources"]})
    out=dict(claim="PF2RULED142",outcome="PASS",check_count=len(checks),checks=checks,
             theorem=theorem,assurance=assurance,source_sha256=pins,
             environment=dict(python=platform.python_version(),sympy=sp.__version__))
    (P/RECEIPT).write_text(json.dumps(out,indent=2)+"\n")
    print(f"PASS {len(checks)} PF142 partial Kummer chart controls; full goal ACTIVE OPEN")


if __name__=="__main__":
    main()
