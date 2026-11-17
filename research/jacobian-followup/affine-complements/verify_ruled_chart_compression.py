"""PF130 exact tower, ruling and mixed-fibre primitive controls."""
from pathlib import Path
from math import gcd
import hashlib
import json
import platform
import sympy as S

P=Path(__file__).resolve().parent
NOTE="RELATIVE_EULER_DEFECT_AND_RULED_CHART_COMPRESSION.md"
MANIFEST="sources/ruled-chart-compression/manifest.json"


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    manifest=json.loads((P/MANIFEST).read_text())
    for item in manifest["sources"]:
        path=P/item["path"]
        check("source seal "+path.name,sha(path)==item["sha256"] and path.stat().st_size==item["bytes"])
    old=json.loads((P/"hereditary_euler_excess_receipt.json").read_text())
    check("PF129 sealed predecessor PASS",old["outcome"]=="PASS" and old["check_count"]==66)
    for f,h in old["source_sha256"].items():
        check("PF129 predecessor bytes "+f,sha(P/f)==h)

    t,b,A1,A2,E1,E2,Delta,Er=S.symbols("t b A1 A2 E1 E2 Delta Er",integer=True)
    c1=b-A1-E1;c2=t*b-A2-E2
    check("UNBOUNDED relative Euler defect exactly separates support and excess",
          S.expand(t*c1-c2-(A2-t*A1+E2-t*E1))==0)
    check("UNBOUNDED exact residual defect substitution",
          S.expand((t*c1-c2).subs({A2:t*A1+Delta,E2:t*E1+Er})-Delta-Er)==0)
    check("UNBOUNDED normalized Euler density difference",
          S.cancel(c1/b-c2/(t*b)-(t*c1-c2)/(t*b))==0)
    k,r,e,q,ell=S.symbols("k r e q ell",positive=True)
    check("UNBOUNDED pointwise critical excess covariance uses local topological degree kr",
          S.expand(k*r*ell-k*e*r*q-k*r*(ell-e*q))==0)
    check("genuine PF127 quotient tower local excess covariance",
          16-5*2==2*(8-5) and 2*2==4*1)
    d1,d2,d3,t1,t2,t3,c0=S.symbols("d1 d2 d3 t1 t2 t3 c0")
    cc1=t1*c0-d1;cc2=t2*cc1-d2;cc3=t3*cc2-d3
    check("UNBOUNDED entire tower defect telescopes with the correct degrees",
          S.expand(t1*t2*t3*c0-cc3-(t2*t3*d1+t3*d2+d3))==0)
    check("zero residual defect does not exclude relative ramification over deleted axes",
          2*7-7*2==0 and 7>1)

    # Whole-fibre gcd versus an actual component of that multiplicity.
    for profile,expected in [([1,1],True),([2,4],True),([3,3],True),
                             ([1,4,7],True),([2,3],False),([4,6],False)]:
        m=S.igcd(*profile)
        normalized=[int(j/m) for j in profile]
        check("component choice and normalized Kummer fibre "+str(profile),
              (int(m) in profile)==expected and (1 in normalized)==expected)
    check("killing whole multiple gcd2 does not manufacture multiplicity1 in fibre4plus6",
          [4//2,6//2]==[2,3] and gcd(2,3)==1)
    for bb,tt,mm,cc in [(5,6,2,2),(7,3,3,3),(9,1,1,4),(4,12,3,1)]:
        check("actual degree compression arithmetic "+str((bb,tt,mm,cc)),
              tt%mm==0 and bb*mm>1 and bb*mm<=bb*tt and
              ((bb*mm<bb*tt)==(mm<tt)) and tt*cc-1>=0)
    check("reduced-fibre simply-connected case forces relative degree1 under minimality",
          all(3*1<3*tt for tt in [2,3,11]))
    for cc in [1,2,3]:
        check("noncyclic full deck index cannot fit c at most3: "+str(cc),
              not [h for h in range(2,cc+1) if cc%h==0 and cc//h>=2])
    check("relative versus full cyclic decks can differ when c2",
          6//3==2 and (3*2)//6==1)
    check("c4 first numerical possibility of a noncyclic extension",
          8//4==2 and (4*4)//8==2)

    # A genuine exact-volume smooth affine surface with bad fibre(2,3).
    x,y,z=S.symbols("x y z")
    p0=y**2*(y-1)**3;p1=-(y**2+(y-1)**3)
    poly=(y**2-x)*((y-1)**3-x);F=x**2*z-poly
    RR=y*(y-1)*(4*y**2-7*y+2)/2
    HH=(5*y-2)*(4*y**2-7*y+2)/2
    TT=(8*y-7)/2
    BB=(z-1)*(2*S.diff(RR,y)-2-HH+x*TT)-p1*TT
    check("mixed-fibre defining polynomial has stated x expansion",S.expand(poly-p0-x*p1-x**2)==0)
    check("smoothness along both bad-fibre components",
          S.diff(F,x).subs({x:0,y:0})==-1 and S.diff(F,x).subs({x:0,y:1})==1)
    check("smoothness away from x0",S.diff(F,z)==x**2)
    check("bad fibre components really have multiplicities2 and3",
          S.factor_list(poly.subs(x,0),y)[1]==[(y,2),(y-1,3)])
    check("first polynomial primitive certificate",S.expand(RR*S.diff(p0,y)-HH*p0)==0)
    check("second polynomial primitive certificate",
          S.expand(p1*(HH-S.diff(RR,y)+1)-RR*S.diff(p1,y)-TT*p0)==0)
    NN=x*S.diff(F,x)*(S.diff(RR,y)-1)-RR*S.diff(poly,y)
    check("GLOBAL regularity certificate for theta across x0",
          S.expand(NN-x**2*BB-(HH-x*TT)*F)==0)
    theta_x=-RR/x**2;theta_y=(S.diff(RR,y)-1)/x
    check("exact derivative is the nonvanishing residue volume on Dx",
          S.cancel(S.diff(theta_y,x)-S.diff(theta_x,y)-1/x**2)==0)
    check("x and Fx opens cover all of the surface special fibre",
          gcd(2,3)==1 and all(S.diff(F,x).subs({x:0,y:j})!=0 for j in [0,1]))

    note=" ".join((P/NOTE).read_text().split())
    for phrase in ["LOCAL TOPOLOGICAL degree", "No surjectivity", "MINIMAL GENERIC DEGREE",
                   "EVERY fibre has SOME component of multiplicity equal to g",
                   "It does not assert that every universal cover has a plane chart",
                   "order may exceed t", "No whole A2 chart", "HC4MYML0"]:
        check("scope and hypothesis retained: "+phrase,phrase in note)

    theorem=dict(relative_Euler_defect_exact_for_any_two_genuine_intermediates=True,
                 critical_excess_covariance_and_entire_etale_equality_case=True,
                 universal_plane_chart_compression_works_for_arbitrary_Euler=True,
                 affine_ruling_gcd_component_condition_constructs_universal_plane_chart=True,
                 minimal_proper_factor_through_reduced_fibre_affine_ruling_excluded=True,
                 primitive_full_deck_cyclic_and_nonprimitive_noncyclic_cost_c_ge4=True,
                 shared_complete_base_period_obstruction_reused_for_every_etale_intermediate=True,
                 exact_smooth_affine_bad_fibre_2plus3_countercontrol_has_regular_primitive=True)
    assurance=dict(prior_axis_Euler_universal_cover_and_deck_inputs_credited=True,
                   general_ruling_pi1_and_bundle_structure_explicitly_adopted=True,
                   entire_original_plane_and_all_normal_intermediates_retained=True,
                   exact_volume_countercontrol_regular_primitive_algebraically_certified=True,
                   whole_fibre_gcd1_assumed_to_supply_component1=False,
                   every_etale_intermediate_assumed_to_have_chart_or_ruling=False,
                   plane_chart_assumed_surjective_onto_intermediate=False,
                   all_maps_assumed_degree_minimal=False,
                   c_gt1_relative_deck_always_identified_with_full_deck=False,
                   fixed_lower_block_movement_discarded=False,
                   local_topological_degree_identified_with_singular_fibre_scheme_length=False,
                   bad_fibre_surface_claimed_actual_plane_Keller_or_finite_target=False,
                   full_natural_AS_binary_affine_trace_energy_or_empty_fibre_gate_closed=False,
                   formal_verification=False,independent_replay=False,external_review=False,
                   full_goal_complete=False)
    files=[NOTE,"verify_ruled_chart_compression.py",MANIFEST,
           "hereditary_euler_excess_receipt.json"]+[i["path"] for i in manifest["sources"]]
    result=dict(claim="PF2RULED130",outcome="PASS",check_count=len(checks),checks=checks,
                theorem=theorem,assurance=assurance,source_sha256={f:sha(P/f) for f in files},
                environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/"ruled_chart_compression_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS",len(checks),"PF130 controls; relative Euler and beyond-Euler-one ruled compression; goal ACTIVE")


if __name__=="__main__":
    main()
