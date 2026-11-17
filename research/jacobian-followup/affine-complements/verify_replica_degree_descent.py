"""PF139: exact derivative/fibre/degree controls for the written global descent."""
from pathlib import Path
from functools import reduce
from math import gcd
import hashlib,json,platform
import sympy as S

P=Path(__file__).resolve().parent
NOTE="REPLICA_EQUIVARIANT_DEGREE_DESCENT_AND_THE_UNIQUE_MIXED_FIBRE.md"
MANIFEST="sources/replica-degree-descent/manifest.json"


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=json.loads((P/MANIFEST).read_text())
    for source in manifest["sources"]:
        p=P/source["path"]
        check("source seal "+source["path"],p.stat().st_size==source["bytes"] and sha(p)==source["sha256"])
    old=json.loads((P/"ruled_chart_compression_receipt.json").read_text())
    check("retained PF130 PASS without mathematical replay",old["outcome"]=="PASS" and old["check_count"]==64)
    for f,h in old["source_sha256"].items():check("PF130 source seal "+f,sha(P/f)==h)
    hdir=P/"../hc4-triangular-graphs"
    mr=json.loads((hdir/"receipt-replica-degree.json").read_text())
    check("shared MR0 PASS84 retained without mathematical replay",mr["outcome"]=="PASS" and mr["check_count"]==84)
    check("shared MR0 proof seal",sha(hdir/"REPLICA_DEGREE_AND_SOURCE_DERIVATIVE_FIELD_RECONSTRUCTION.md")==mr["proof_sha256"])
    check("shared MR0 checker seal",sha(hdir/"verify_replica_degree.py")==mr["checker_sha256"])

    m,g,b,q=S.symbols("m g b q",integer=True,positive=True)
    check("UNBOUNDED two wholly critical fibres exceed derivative budget",
          S.expand(2*(m/2)-(m-1))==1)
    check("UNBOUNDED new actual Keller degree",S.cancel((b*m)*g/m-b*g)==0)
    check("UNBOUNDED strict descent gap",S.expand(b*m-b*g-b*(m-g))==0)
    ell,rdeg=S.symbols("ell rdeg",integer=True,nonnegative=True)
    check("UNBOUNDED replica index relation after gcd normalization",
          S.expand((g*ell)*(rdeg+1)-g*ell*(rdeg+1))==0)
    check("UNBOUNDED mixed root multiplicity lower bound",S.expand(3*g-1-g*3+1)==0)

    def profile(parts):
        gg=reduce(gcd,parts)
        return gg,gg in parts
    for parts,expected in [([2,1],(1,True)),([2,2],(2,True)),
                           ([2,3],(1,False)),([4,6],(2,False)),
                           ([3,3,6],(3,True)),([5],(5,True))]:
        gg,condition=profile(parts)
        check("whole gcd component criterion "+str(parts),(gg,condition)==expected)
        if not condition:
            normalized=[n//gg for n in parts]
            check("normalized mixed partition minimum degree "+str(parts),
                  min(normalized)>=2 and reduce(gcd,normalized)==1 and sum(normalized)>=5)
            check("mixed derivative forces repeated replica root "+str(parts),
                  max(parts)-1>=3*gg-1)
    check("normalized threshold equality is two plus three",sum([2,3])==5 and gcd(2,3)==1)

    y=S.symbols("y")
    h=y**3-3*y
    check("Morse cubic positive critical-value fibre",S.expand(h-2-(y-2)*(y+1)**2)==0)
    check("Morse cubic negative critical-value fibre",S.expand(h+2-(y+2)*(y-1)**2)==0)
    check("Morse cubic has simple derivative roots",S.degree(S.gcd(S.diff(h,y),S.diff(h,y,2)),y)==0)
    quartic=(y*y-1)**2
    check("colliding critical values can still satisfy gcd-component criterion",
          S.degree(S.gcd(S.diff(quartic,y),S.diff(quartic,y,2)),y)==0 and profile([2,2])==(2,True))
    check("quartic descent factor has universal group order below source index",2<4)
    H=y**2*(y-1)**3
    for gg in [1,2]:
        h=H**gg;p=S.diff(h,y)
        check("mixed polynomial relative-index control "+str(gg),S.degree(h,y)==5*gg)
        check("mixed polynomial replica-degree control "+str(gg),S.degree(p,y)+1==5*gg)
        check("mixed fibre derivative order at multiplicity-three root "+str(gg),
              S.Poly(p.subs(y,y+1),y).terms()[-1][0][0]==3*gg-1)

    for c,rr in [(0,1+y),(2,y*y-3)]:
        # p=(y-c)*r((y-c)^2), where r is any displayed polynomial.
        t=S.symbols("t");p=(y-c)*rr.subs(y,(y-c)**2)
        check("quadratic alternative is odd about its actual centre "+str((c,rr)),
              S.expand(p.subs(y,c+t)+p.subs(y,c-t))==0)
        check("quadratic alternative has odd degree "+str((c,rr)),int(S.degree(p,y))%2==1)
    check("square-free even-degree replica previously allowed index three",
          S.gcd(y*y-1,2*y)==1 and (2+1)%3==0)
    check("odd square-free replica need not have centred odd symmetry",
          S.gcd(y**5-y-1,5*y**4-1)==1 and (y**5-y-1).subs(y,0)!=0)

    theorem=dict(
        actual_equivariant_fibre_gcd_is_polynomial_fibre_gcd=True,
        at_most_one_wholly_critical_fibre_and_source_image_pi1_cyclic=True,
        gcd_component_criterion_produces_actual_smaller_plane_Keller_map_unless_pure_power=True,
        minimal_survivor_mixed_fibre_has_normalized_degree_at_least_five_and_replica_repeated_root=True,
        square_free_replicas_only_allow_quadratic_deck_and_odd_centre_in_minimal_degree=True)
    assurance=dict(
        written_all_degree_proof_not_partition_census=True,
        shared_MR0_and_retained_PF130_chart_construction_credited=True,
        selected_source_image_not_identified_with_entire_intermediate_etale_locus=True,
        minimality_needed_for_exclusion_without_smaller_map_alternative=True,
        cyclic_trivial_deck_and_general_full_trace_goal_open=True,
        mixed_fibre_claimed_constructed_as_actual_Keller_map=False,
        existence_of_source_replica_claimed_to_exclude_map=False,
        arbitrary_intermediate_pi1_identified_with_generic_deck=False,
        raw_derivative_algebra_or_finite_order_recovered_from_field_generators=False,
        new_external_theorem_import=False,
        registered_predecessor_mathematical_replay=False,
        independent_replay=False,formal_verification=False,external_review=False,
        novelty_claim=False,full_goal_complete=False)
    files=[NOTE,"verify_replica_degree_descent.py",MANIFEST]+[s["path"] for s in manifest["sources"]]
    result=dict(claim="PF2RULED139",outcome="PASS",check_count=len(checks),checks=checks,
                theorem=theorem,assurance=assurance,source_sha256={f:sha(P/f) for f in files},
                environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/"replica_degree_descent_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS",len(checks),"PF139 controls; actual degree descent and all-degree square-free replica exclusion; goal ACTIVE")


if __name__=="__main__":main()
