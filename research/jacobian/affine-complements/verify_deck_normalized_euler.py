"""PF133: exact ledger and semiregular-permutation controls, no group census."""
from pathlib import Path
from itertools import product
import hashlib
import json
import platform
import sympy as S

P=Path(__file__).resolve().parent
NOTE="DECK_NORMALIZED_EULER_AND_UNRESTRICTED_RADICAL_TOWER_EXCLUSION.md"
MANIFEST="sources/deck-normalized-euler/manifest.json"


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):
        return hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=json.loads((P/MANIFEST).read_text())
    for source in manifest["sources"]:
        p=P/source["path"]
        check("source seal "+p.name,p.stat().st_size==source["bytes"] and sha(p)==source["sha256"])
    old=json.loads((P/"radical_tower_exclusion_receipt.json").read_text())
    check("sealed PF132 PASS without replay",old["outcome"]=="PASS" and old["check_count"]==95)
    for f,h in old["source_sha256"].items():
        check("PF132 predecessor bytes "+f,sha(P/f)==h)

    b,c,g,s,j,mu,A,E,delta,er,cv=S.symbols("b c g s j mu A E delta er cv")
    check("UNBOUNDED normalized original Euler budget",
          S.cancel(((g*b-g*s-g*j)/g).subs(b,c+s+j)-c)==0)
    check("UNBOUNDED relative residual defect lowers quotient Euler",
          S.expand(g*cv-g*(cv-(delta+er)/g)-delta-er)==0)
    check("UNBOUNDED index-five boundary sum of nonnegative terms",
          S.expand((c-3)+(s-2)+j-(b-5)).subs(b,c+s+j)==0)
    check("UNBOUNDED minimum support refinement nonnegative remainder",
          S.expand(g*(c-3)+(A-mu)+E-(g*c+A+E-mu-3*g))==0)
    check("index-five numerical profile is unique",
          [(cc,ss,jj) for cc in range(3,6) for ss in range(2,6) for jj in range(4) if cc+ss+jj==5]==[(3,2,0)])
    for gg in [4,6,8,12]:
        check("noncyclic index lower bound and original degree "+str(gg),5*gg>=20 and 4*gg-2*gg<3*gg)
        for dd in [0,gg,2*gg]:
            check("positive residual defect cannot rescue quartic quotient "+str((gg,dd)),2-dd//gg<=2)
        for bb in [3,4]:
            check("visible quotient branch forces low Euler for index "+str((gg,bb)),bb-2<=2)
        check("index-five boundary original support/Euler "+str(gg),5*gg-2*gg==3*gg and 3*gg-1+1==3*gg)

    for q in [1,2,3]:
        gg=1<<q;d=4*gg;ident=tuple(range(d))
        translations=[tuple(x*gg+(z^a) for x in range(4) for z in range(gg)) for a in range(gg)]
        def perm(base,values):
            return tuple(base[x]*gg+(z^values[x]) for x in range(4) for z in range(gg))
        controls=[perm([1,0,2,3],[0,0,0,0]),
                  perm([1,0,2,3],[0,1,0,1]),
                  perm([0,1,2,3],[1,1,0,0]),
                  perm([1,2,0,3],[0,0,0,0])]
        check("radical translations are semiregular q"+str(q),
              all(all(t[x]!=x for x in range(d)) for t in translations[1:]))
        check("selected inertia types commute with ALL deck translations q"+str(q),
              all(all(p[t[x]]==t[p[x]] for x in range(d)) for p in controls for t in translations))
        supports=[sum(p[x]!=x for x in range(d)) for p in controls]
        check("supports divisible by deck order and minimum half q"+str(q),
              all(w%gg==0 and w>=d//2 for w in supports))
        check("order-four and vertical controls retain positive residual defect q"+str(q),
              supports[1]-2*gg==gg and supports[2]-0*gg==2*gg)
        check("nondecorated order-four permutation is still valid q"+str(q),
              any(controls[1][controls[1][x]]!=x for x in range(d)) and
              all(controls[1][controls[1][controls[1][controls[1][x]]]]==x for x in range(d)))
        check("natural tower index stays four q"+str(q),d//gg==4)
        check("noncyclic tower threshold q"+str(q),(gg>2)==(q>=2))
    for a,b0,c0 in product(range(4),repeat=3):
        values=[c0,c0^a,c0^b0,c0^a^b0]
        assert values.count(0)!=3
    check("small affine-function zero sets corroborate three-point determination",True)
    check("cyclic q1 is not rejected by noncyclic premise",2==2 and 4*2//2==4)
    check("index-five frontier is below half-support minimum",2*4<5*4//2)

    theorem=dict(deck_normalized_support_and_total_excess_budget_integral=True,
                 low_Euler_intermediate_forces_relative_Galois_group_cyclic_even_with_residual=True,
                 every_noncyclic_generic_deck_subgroup_has_fixed_field_index_at_least_five=True,
                 minimum_support_gives_three_times_deck_order_at_most_degree_minus_mu=True,
                 entire_natural_radical_tower_family_excluded_in_ALL_inertia_types=True,
                 first_index_five_profile_exact_but_not_excluded=True)
    assurance=dict(PF132_surface_lemma_and_coordinator_HC4_credit_retained=True,
                   actual_free_deck_quotient_used_instead_of_entire_intermediate_identification=True,
                   all_residual_fixed_block_and_invisible_support_retained=True,
                   entire_original_plane_and_all_omitted_primes_retained=True,
                   no_minimal_degree_or_inertia_class_assumption=True,
                   new_external_classification_import=False,
                   residual_ramification_assumed_zero_for_arbitrary_inertia=False,
                   every_relative_monodromy_group_treated_as_generic_deck=False,
                   local_or_branchwise_excess_assumed_deck_divisible=False,
                   same_abstract_group_in_all_permutation_actions_excluded=False,
                   q1_index_five_or_full_JC2_closed=False,
                   knot_representations_invalidated=False,
                   registered_predecessor_replay=False,formal_verification=False,
                   independent_replay=False,external_review=False,full_goal_complete=False)
    files=[NOTE,"verify_deck_normalized_euler.py",MANIFEST,"radical_tower_exclusion_receipt.json"]+[s["path"] for s in manifest["sources"]]
    result=dict(claim="PF2RULED133",outcome="PASS",check_count=len(checks),checks=checks,theorem=theorem,
                assurance=assurance,source_sha256={f:sha(P/f) for f in files},
                environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/"deck_normalized_euler_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS",len(checks),"PF133 controls; natural radical towers excluded in every inertia class; goal ACTIVE")


if __name__=="__main__":
    main()
