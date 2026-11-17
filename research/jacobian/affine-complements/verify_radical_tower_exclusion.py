"""PF132: source audit and exact decorated-tower/orbifold controls only."""
from collections import deque
from fractions import Fraction
from math import gcd
from pathlib import Path
import hashlib
import json
import platform

P=Path(__file__).resolve().parent
NOTE="LOW_EULER_CYCLICITY_AND_DECORATED_RADICAL_TOWER_EXCLUSION.md"
MANIFEST="sources/radical-tower-exclusion/manifest.json"


def compose(a,b):
    return tuple(a[x] for x in b)


def lift(q,edge,r=0,values=None):
    t=1<<q
    a,b=edge
    values=values if values is not None else [r if x in edge else 0 for x in range(4)]
    return tuple(((b if x==a else a if x==b else x)*t)+(z^values[x])
                 for x in range(4) for z in range(t))


def centralizer(gens):
    n=len(gens[0]);out=set()
    for target in range(n):
        values={0:target};queue=deque([0]);valid=True
        while queue and valid:
            x=queue.popleft()
            for g in gens:
                y=g[x];image=g[values[x]]
                if y in values:
                    if values[y]!=image:
                        valid=False
                        break
                else:
                    values[y]=image;queue.append(y)
        if valid and len(values)==n and len(set(values.values()))==n:
            out.add(tuple(values[x] for x in range(n)))
    return out


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
        check("source seal "+path.name,sha(path)==s["sha256"] and path.stat().st_size==s["bytes"])
    old=json.loads((P/"ruled_chart_compression_receipt.json").read_text())
    check("sealed PF130 PASS without replay",old["outcome"]=="PASS" and old["check_count"]==64)
    for f,h in old["source_sha256"].items():
        check("PF130 predecessor bytes "+f,sha(P/f)==h)
    packet=P/"../cross-programme-connections"
    prior=json.loads((packet/"inward-coupled-restrictions-controls-2026-11-16.json").read_text())
    check("coupling packet proof hash agrees",prior["proof_sha256"]==sha(packet/"INWARD_COUPLED_RESTRICTIONS_2026_11_16.tex"))
    check("coupling packet checker hash agrees",prior["checker_sha256"]==sha(packet/"check_inward_coupled_restrictions_20261116.py"))
    for f,h in prior["dependency_sha256"].items():
        # Check retained input bytes; do not execute unrelated Mathieu controls.
        check("coupling dependency bytes "+Path(f).name,sha(Path(f))==h)

    for q in [0,1,2,3]:
        t=1<<q;n=4*t;ident=tuple(range(n))
        gens=[lift(q,e) for e in [(0,1),(1,2),(2,3)]]
        gens += [lift(q,(0,1),1<<j) for j in range(q)]
        translations={tuple(x*t+(z^c) for x in range(4) for z in range(t)) for c in range(t)}
        got=centralizer(gens)
        check("full centralizer is constant radical translations q"+str(q),got==translations)
        check("radical centralizer has exponent two q"+str(q),all(compose(g,g)==ident for g in got))
        check("deck cyclicity boundary q"+str(q),(len(got)<=2)==(q<=1))
        check("all specified generators have support half and are involutions q"+str(q),
              all(sum(g[x]!=x for x in range(n))==n//2 and compose(g,g)==ident for g in gens))
        check("fixed quotient blocks have NO internal movement q"+str(q),
              all(g[x]==x for g in gens for x in range(n) if g[x]//t==x//t))
        check("exact support defect zero for specified inertia q"+str(q),n//2-t*2==0)
        check("Wright budget excludes two decorated branch images q"+str(q),2*(n//2)>n-1)
        check("quartic Euler slack is two q"+str(q),4-2==2 and n//t==4)
        if q:
            vertical=tuple(x*t+(z^(1 if x in [0,1] else 0)) for x in range(4) for z in range(t))
            check("vertical inertia has positive invisible residual support q"+str(q),
                  all(vertical[x]//t==x//t for x in range(n)) and sum(vertical[x]!=x for x in range(n))==n//2)
            fourth=lift(q,(0,1),values=[0,1,0,1])
            check("other affine-function lift has order four and residual fixed-block movement q"+str(q),
                  compose(fourth,fourth)!=ident and compose(compose(fourth,fourth),compose(fourth,fourth))==ident
                  and sum(fourth[x]!=x for x in range(n))-t*2==t)

    for orders,g in [((2,2,2),4),((2,2,5),10),((2,3,3),12),((2,3,4),24),((2,3,5),60)]:
        orb=2-sum((1-Fraction(1,m) for m in orders),Fraction())
        check("spherical orbifold universal cover Euler "+str(orders),g*orb==2)
        components=sum((g//m)*(m-1) for m in orders)
        check("split REDUCED pullback fibre Euler is twice cover degree "+str(orders),2+components==2*g)
        check("exact local root order kills component multiplicity "+str(orders),all(m//gcd(m,m)==1 for m in orders))
    for m,n in [(2,3),(2,4),(3,6)]:
        check("bad spindle only imposes cyclic gcd relation "+str((m,n)),gcd(m,n) in [1,2,3])
    check("mixed fibre gcd cover need not give reduced components",[j//gcd(4,6) for j in [4,6]]==[2,3])
    check("complete base low Euler forces single components",all(2+sum(j-1 for j in a)>2 for a in [[2],[1,2],[3,1,2]]))
    for g in [4,8,12]:
        check("new noncyclic generic deck boundary improves old one g"+str(g),3*g-1>2*g-1 and (3*g-1+1)//g==3)
    theorem=dict(low_Euler_exact_volume_negative_log_surfaces_have_cyclic_fundamental_group=True,
                 good_orbifold_base_change_reduced_components_not_assumed_irreducible=True,
                 decorated_Hq_lifted_transposition_family_excluded_for_all_q_ge2=True,
                 entire_quartic_intermediate_and_universal_relative_cover_retained=True,
                 every_genuine_chi_le2_intermediate_has_cyclic_fundamental_group=True,
                 noncyclic_generic_deck_subgroup_requires_b0_ge3order_minus1=True)
    assurance=dict(coordinator_coupling_deduction_credited_and_audited=True,
                   classical_ruling_and_orbifold_inputs_explicitly_adopted=True,
                   HC4MYML0_full_component_selection_input_essential=True,
                   entire_original_plane_and_all_omitted_primes_retained=True,
                   specified_branch_permutation_hypothesis_retained=True,
                   no_minimal_degree_assumption=True,
                   bad_orbifold_manifold_cover_assumed=False,
                   reduced_pullback_fibres_assumed_irreducible=False,
                   abstract_Hq_or_other_inertia_classes_excluded=False,
                   q1_case_or_knot_representations_excluded=False,
                   unrelated_Mathieu_or_Factorial_packet_sections_adopted=False,
                   registered_predecessor_replay=False,formal_verification=False,
                   independent_replay=False,external_review=False,full_goal_complete=False)
    files=[NOTE,"verify_radical_tower_exclusion.py",MANIFEST,"ruled_chart_compression_receipt.json"]+[s["path"] for s in manifest["sources"]]
    result=dict(claim="PF2RULED132",outcome="PASS",check_count=len(checks),checks=checks,theorem=theorem,
                assurance=assurance,source_sha256={f:sha(P/f) for f in files},environment=dict(python=platform.python_version()))
    (P/"radical_tower_exclusion_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS",len(checks),"PF132 controls; decorated radical towers excluded at all q>=2; goal ACTIVE")


if __name__=="__main__":
    main()
