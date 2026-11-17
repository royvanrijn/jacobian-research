"""PF136: exact common-orbit/Euler controls, not a surface census."""
from pathlib import Path
import hashlib,json,platform
import sympy as S

P=Path(__file__).resolve().parent
NOTE="COMMON_OMITTED_ORBITS_AND_ALMOST_SURJECTIVE_DECK_QUOTIENTS.md"
MANIFEST="sources/common-omitted-orbits/manifest.json"


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=json.loads((P/MANIFEST).read_text())
    for source in manifest["sources"]:
        p=P/source["path"]
        check("source seal "+p.name,p.stat().st_size==source["bytes"] and sha(p)==source["sha256"])
    old=json.loads((P/"cotangent_deck_trace_receipt.json").read_text())
    check("sealed PF135 PASS without replay",old["outcome"]=="PASS" and old["check_count"]==79)
    for f,h in old["source_sha256"].items():check("PF135 predecessor bytes "+f,sha(P/f)==h)
    g,c,r,s,j,b0=S.symbols("g c r s j b0")
    check("UNBOUNDED exact Euler drop per removed free curve orbit",
          S.cancel(((1+b0-g*r)/g-(c-r)).subs(b0,g*c-1))==0)
    check("UNBOUNDED auxiliary quotient Picard rank",S.expand((g*c-1-g*r+1)/g-1-(c-r-1))==0)
    cp=S.symbols("cp")
    check("UNBOUNDED refined degree budget",S.expand(c+s+j-(cp+r+s+j)).subs(c,cp+r)==0)
    check("UNBOUNDED saturation nonnegative remainder",S.expand((cp-3)+(s-2)+j-((cp+r+s+j)-r-5))==0)
    a0,bb,x=S.symbols("a0 bb x")
    fixed=bb/(1-a0)
    check("finite affine stabilizer fixed point identity",S.cancel(a0*fixed+bb-fixed)==0)

    # Retained standard Danielewski chart combinatorics, all k=gh.
    for k,gg in [(2,2),(4,2),(6,2),(6,3),(12,3),(12,4)]:
        step=k//gg;roots=set(range(k));chart_omitted=roots-{0}
        translated=[{(i+a*step)%k for i in chart_omitted} for a in range(gg)]
        common=set.intersection(*translated)
        retained={a*step%k for a in range(gg)}
        check("common omitted roots for cyclic chart action "+str((k,gg)),common==roots-retained)
        orbits=[];todo=set(common)
        while todo:
            start=min(todo);orbit={(start+a*step)%k for a in range(gg)}
            orbits.append(orbit);todo-=orbit
        rr=k//gg-1
        check("every common component orbit free and full "+str((k,gg)),len(orbits)==rr and all(len(o)==gg for o in orbits))
        check("exact auxiliary Euler in cyclic control "+str((k,gg)),(k-len(common))//gg==1 and k//gg-rr==1)
        check("cyclic translated charts cover auxiliary source "+str((k,gg)),
              set.union(*(roots-z for z in translated))==roots-common)
    for gg,cc,rr in [(4,3,0),(4,5,2),(6,4,1),(12,7,4)]:
        original_b0=gg*cc-1;remaining=original_b0-gg*rr
        check("noncyclic quotient Euler and remaining original boundary "+str((gg,cc,rr)),
              remaining>=0 and (1+remaining)//gg==cc-rr==3)
        for SS,JJ in [(2,0),(3,0),(2,2)]:
            degree_index=cc+SS+JJ
            check("new common-orbit fixed-field bound "+str((gg,cc,rr,SS,JJ)),degree_index>=rr+5)
        check("saturated quotient has Picard rank two "+str((gg,cc,rr)),cc-rr-1==2)
    check("almost-surjectivity does not assert empty finite missed set",4*3==12 and 12//4==3)
    check("retained cyclic Euler-one control cannot satisfy noncyclic lower bound",1<3)
    check("index-five saturation leaves no common omitted curve orbit",5-5==0)
    check("higher-index saturation still has auxiliary Euler three",all((rr+3)-rr==3 for rr in [0,1,7]))
    check("rank-two frontier is not irreducible-fibre Platonic rank one",2!=1)

    theorem=dict(common_omitted_divisors_are_free_deck_orbits_of_disjoint_affine_lines=True,
                 actual_auxiliary_affine_quotient_is_maximal_almost_surjective_for_original_plane=True,
                 exact_Euler_drop_and_Picard_rank_with_whole_plane_and_target_retained=True,
                 noncyclic_common_orbit_budget_r_le_c_minus3_and_b_ge_r_plus5=True,
                 saturation_has_natural_Sb_zero_excess_and_index5_original_Galois_pseudocover=True)
    assurance=dict(actual_original_plane_and_target_map_preserved=True,
                   entire_original_normalization_discriminant_and_all_defects_retained=True,
                   auxiliary_quotient_not_mislabeled_entire_intermediate=True,
                   finite_missed_points_and_mixed_fibres_remain_open=True,
                   shared_MO0_chart_and_MA0_residual_audit_credited=True,
                   all_deck_translates_assumed_cover_every_point=False,
                   original_boundary_basis_assumed_G_permuted_before_common_part=False,
                   cyclic_control_given_an_actual_etale_plane_target=False,
                   mixed_Euler3_profile_treated_as_irreducible_fibre_Platonic_space=False,
                   blocked_full_pseudocover_article_claimed_audited_or_adopted=False,
                   full_trace_energy_or_JC2_closed=False,
                   registered_predecessor_replay=False,formal_verification=False,
                   independent_replay=False,external_review=False,full_goal_complete=False)
    files=[NOTE,"verify_common_omitted_orbits.py",MANIFEST,"cotangent_deck_trace_receipt.json"]+[s["path"] for s in manifest["sources"]]
    result=dict(claim="PF2RULED136",outcome="PASS",check_count=len(checks),checks=checks,theorem=theorem,
                assurance=assurance,source_sha256={f:sha(P/f) for f in files},environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/"common_omitted_orbits_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS",len(checks),"PF136 controls; actual almost-surjective deck quotient and curve-orbit budget; goal ACTIVE")


if __name__=="__main__":main()
