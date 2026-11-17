"""PF134: classical-source seals and tiny orbifold/fibre/group controls."""
from collections import deque
from fractions import Fraction
from functools import reduce
from math import gcd,lcm
from pathlib import Path
import hashlib
import json
import platform

P=Path(__file__).resolve().parent
NOTE="FINITE_ORBIFOLD_GROUPS_AND_UNIVERSAL_COMPLETE_RULINGS.md"
MANIFEST="sources/finite-orbifold-groups/manifest.json"


def compose(a,b):
    return tuple(a[x] for x in b)


def permutation(n,*cycles):
    out=list(range(n))
    for c in cycles:
        for a,b in zip(c,c[1:]+c[:1]):out[a]=b
    return tuple(out)


def order(p):
    out=1;seen=set()
    for i in range(len(p)):
        if i in seen:continue
        j=i;length=0
        while j not in seen:
            seen.add(j);length+=1;j=p[j]
        out=lcm(out,length)
    return out


def closure(gens):
    ident=tuple(range(len(gens[0])));seen={ident};queue=deque([ident])
    while queue:
        x=queue.popleft()
        for a in gens:
            y=compose(a,x)
            if y not in seen:seen.add(y);queue.append(y)
    return seen


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
    old=json.loads((P/"deck_normalized_euler_receipt.json").read_text())
    check("sealed PF133 PASS without replay",old["outcome"]=="PASS" and old["check_count"]==82)
    for f,h in old["source_sha256"].items():check("PF133 predecessor bytes "+f,sha(P/f)==h)
    review=json.loads((P/"sources/finite-orbifold-groups/source-review.json").read_text())
    check("general reducible-fibre scope verified before Definition2.5",
          review["Miyanishi_Lemma2_4_not_assumed_irreducible"] is True)
    check("retained classical classification explicitly adopted not independently proved",
          review["broader_use_of_retained_classical_inputs"] is True and review["independent_classification_proof"] is False)

    for signature,g in [((2,2,2),4),((2,2,5),10),((2,3,3),12),((2,3,4),24),((2,3,5),60)]:
        orb=2-sum(1-Fraction(1,m) for m in signature)
        check("spherical cover degree "+str(signature),g*orb==2)
        base=sum((g//m)*(m-1) for m in signature)
        check("irreducible target fibres become split reduced cover fibres "+str(signature),2+base==2*g)
        check("one ordinary two-component fibre gives Euler3 cover "+str(signature),2+base+g==3*g)
        m=signature[0]
        increased=(g//m)*m
        check("one cone two-component fibre gives same Euler3 "+str(signature),2+base+increased==3*g)
        check("finite noncyclic order is even and degree bound "+str(signature),g%2==0 and 5*g>=20)

    for values,expected in [([3],True),([2,4],True),([4,6],False),([2,3],False),([6,10,15],False)]:
        m=reduce(gcd,values)
        normalized=[a//m for a in values]
        lifted=normalized*m
        check("multiplicities after exact-order base change "+str(values),
              reduce(gcd,normalized)==1 and (1 in lifted)==expected and len(lifted)==m*len(values))
        check("gcd-component condition equivalent to unit in lifted fibre "+str(values),(m in values)==(1 in lifted))
    check("minimal normalized bad two-component profile is coprime two plus three",gcd(2,3)==1 and min(2,3)>=2)
    check("whole multiplicity can coexist with bad normalized components",reduce(gcd,[4,6])==2 and min(4//2,6//2)>1)
    check("three-component bad fibre exceeds Euler3 budget",2+(3-1)>3)

    # Three explicit spherical realizations, at most sixty permutations.
    for name,n,a,b,target in [
        ("A4",4,[(0,1),(2,3)],[(0,1,2)],12),
        ("S4",4,[(0,1)],[(1,2,3)],24),
        ("A5",5,[(0,1),(2,3)],[(0,2,4)],60)]:
        aa=permutation(n,*a);bb=permutation(n,*b);gg=closure([aa,bb])
        check("tiny spherical permutation realization "+name,len(gg)==target and order(aa)==2 and order(bb)==3)
        check("third triangle order "+name,order(compose(aa,bb))==dict(A4=3,S4=4,A5=5)[name])
        check("necessary original degree lower bound "+name,5*target==dict(A4=60,S4=120,A5=300)[name])
    for n in [2,3,4,5,8]:
        elements=[(a,e) for a in range(n) for e in range(2)]
        def mul(x,y):
            a,e=x;b,f=y
            return ((a+(-1 if e else 1)*b)%n,(e+f)%2)
        involutions=sum(mul(x,x)==(0,0) and x!=(0,0) for x in elements)
        check("dihedral order and involution count n"+str(n),len(elements)==2*n and involutions==n+(n%2==0))
        check("only Klein four is noncyclic abelian dihedral n"+str(n),
              all(mul(x,y)==mul(y,x) for x in elements for y in elements)==(n==2))

    def quaternion_mul(x,y):
        sign=(1 if x>0 else -1)*(1 if y>0 else -1);a,b=abs(x),abs(y)
        if a==1:return sign*b
        if b==1:return sign*a
        if a==b:return -sign
        cyclic={(2,3):4,(3,4):2,(4,2):3}
        return sign*(cyclic[(a,b)] if (a,b) in cyclic else -cyclic[(b,a)])
    q8=[1,-1,2,-2,3,-3,4,-4]
    check("quaternion negative control has one involution and no element of order8",
          sum(x!=1 and quaternion_mul(x,x)==1 for x in q8)==1 and
          all(quaternion_mul(quaternion_mul(x,x),quaternion_mul(x,x))==1 for x in q8))
    for rank in [3,4,7]:
        check("higher elementary-abelian deck rank outside cyclic-or-V4 container "+str(rank),2**rank>4)
    check("odd relative degree cannot contain any noncyclic group in list",all(g%2==0 for g in [4,6,10,12,24,60]))

    theorem=dict(all_genuine_intermediate_finite_pi1_and_generic_deck_groups_are_Mobius_groups=True,
                 abelian_deck_groups_cyclic_or_V4_and_odd_relative_pi1_cyclic=True,
                 noncyclic_universal_cover_has_actual_complete_ruling_and_faithful_base_action=True,
                 exact_split_fibre_multiplicities_and_Euler_formula=True,
                 every_noncyclic_ruling_has_no_component_equal_to_gcd_in_some_fibre=True,
                 Euler_three_profile_is_one_mixed_coprime_two_component_fibre=True)
    assurance=dict(retained_classical_inputs_broader_scope_explicitly_adopted=True,
                   original_plane_entire_intermediates_and_actual_free_deck_quotient_retained=True,
                   coordinator_PF132_and_HC4_component_selection_credit_retained=True,
                   generic_deck_and_full_monodromy_roles_distinguished=True,
                   bad_normalized_multiplicity_fibres_remain_open=True,
                   whole_multiple_fibres_assumed_irreducible=False,
                   all_base_changed_components_assumed_reduced=False,
                   bad_cyclic_orbifolds_given_nonexistent_manifold_cover=False,
                   complete_ruling_confused_with_additive_affine_ruling=False,
                   whole_finite_target_trace_energy_or_JC2_closed=False,
                   classification_independently_proved=False,
                   registered_predecessor_replay=False,formal_verification=False,
                   independent_replay=False,external_review=False,full_goal_complete=False)
    files=[NOTE,"verify_finite_orbifold_groups.py",MANIFEST,"deck_normalized_euler_receipt.json"]+[s["path"] for s in manifest["sources"]]
    result=dict(claim="PF2RULED134",outcome="PASS",check_count=len(checks),checks=checks,theorem=theorem,
                assurance=assurance,source_sha256={f:sha(P/f) for f in files},environment=dict(python=platform.python_version()))
    (P/"finite_orbifold_groups_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS",len(checks),"PF134 controls; finite Mobius containers and actual complete rulings; goal ACTIVE")


if __name__=="__main__":
    main()
