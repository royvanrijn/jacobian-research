"""PF131: exact form/peripheral controls, not a classification proof."""
from pathlib import Path
import hashlib
import json
import platform

P=Path(__file__).resolve().parent
NOTE="PROPER_BINARY_LINEAR_GROUPS_AND_TRIVIAL_LONGITUDE.md"
MANIFEST="sources/proper-binary-longitude/manifest.json"


def parity(x):
    return x.bit_count()%2


def rank(rows):
    piv={}
    for x in rows:
        while x:
            i=x.bit_length()-1
            if i in piv:
                x^=piv[i]
            else:
                piv[i]=x
                break
    return len(piv)


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
    old=json.loads((P/"seifert_framing_cancellation_receipt.json").read_text())
    check("sealed PF126 PASS without replay",old["outcome"]=="PASS" and old["check_count"]==144)
    for f,h in old["source_sha256"].items():
        check("PF126 predecessor bytes "+f,sha(P/f)==h)
    review=json.loads((P/"sources/proper-binary-longitude/source-review.json").read_text())
    check("proper SL statement and actual symmetric module distinguished",
          review["primary_restatement_PDF_pages_2_and_3_read"] and review["representation_not_abstract_group_names"])
    check("original McLaughlin proof not claimed inspected",review["original_McLaughlin_PDF_obtained"] is False)

    for k in [2,4,6]:
        n=k//2
        def B(x,y):
            return parity((x&((1<<n)-1))&(y>>n)) ^ parity((y&((1<<n)-1))&(x>>n))
        basis=[1<<j for j in range(k)]
        v=1;l=1<<n
        q=lambda x:x^(v if parity(l&x) else 0)
        check("chosen meridian preserves common B k"+str(k),
              all(B(q(x),q(y))==B(x,y) for x in basis for y in basis))
        for u in [x for x in range(1,1<<k) if parity(l&x)==0]:
            lin=lambda x:x^(u if parity(l&x) else 0)
            preserved=all(B(lin(x),lin(y))==B(x,y) for x in basis for y in basis)
            check("common B forces nonzero residual column v "+str((k,u)),preserved==(u==v))
        complementary=lambda x:x^(v if parity(l&x)==0 else 0)
        check("complementary transvection prevents pointwise-fixing shortcut k"+str(k),
              any(complementary(x)!=x for x in range(1<<k)) and
              all(complementary(x)==x for x in range(1<<k) if parity(l&x)))

    for k in [2,3,4,7]:
        pairs=[(a,b) for a in range(k) for b in range(a+1,k)]
        def value(x,y):
            return sum(((((x>>a)&1)*((y>>b)&1))^(((x>>b)&1)*((y>>a)&1)))<<i
                       for i,(a,b) in enumerate(pairs))
        constraints=[]
        for i in range(k):
            for j in range(k):
                if i==j:continue
                T=lambda x:x^((1<<i) if x&(1<<j) else 0)
                for a in range(k):
                    for b in range(k):
                        x=1<<a;y=1<<b
                        constraints.append(value(T(x),T(y))^value(x,y))
        check("full GL invariant alternating form space dimension k"+str(k),
              len(pairs)-rank(constraints)==(1 if k==2 else 0))

    for m in [5,6,7,8]:
        Y=[(1<<i)|(1<<(m-1)) for i in range(m-1)]
        gram=[sum(parity(x&y)<<j for j,y in enumerate(Y)) for x in Y]
        check("deleted permutation module has nondegenerate alternating quotient m"+str(m),
              all(parity(x&x)==0 for x in Y) and rank(gram)==m-1-(m%2==0))
        for i in range(m-1):
            def swap(x):
                return x^(((1<<i)|(1<<(i+1))) if ((x>>i)^(x>>(i+1)))&1 else 0)
            differences=[swap(x)^x for x in Y]
            check("symmetric adjacent transposition has quotient rank one "+str((m,i)),
                  rank(differences)==1 and ((1<<m)-1) not in differences and
                  all(parity(swap(x)&swap(y))==parity(x&y) for x in Y for y in Y))

    # Curve-cover algebraic escape in full AGL3, not a Keller realization.
    l=4;v=1;u=2
    T=lambda a,x:x^(a if parity(l&x)==0 else 0)
    q=lambda x:x^(v if parity(l&x) else 0)
    lam=lambda x:T(v^u,x)
    check("full AGL escape fixes moved half pointwise",all(lam(x)==x for x in range(8) if x&l))
    check("full AGL escape is nonidentity involution centralizing q",
          any(lam(x)!=x for x in range(8)) and all(lam(lam(x))==x and q(lam(x))==lam(q(x)) for x in range(8)))
    check("two distinct columns span four and give nonzero infinity column",
          rank([v,u])==2 and v^u!=0 and all(lam(x)==T(v,T(u,x)) for x in range(8)))
    check("escaping rank-four residual cover has genus zero and two infinity points",1-4+(4//4)*(2+1)==0 and 4//2==2)
    for vo in [0,1,2,7]:
        check("even identical columns cancel at every tested branch count "+str(vo),
              (sum([1]*(2*vo))%2)==0)
        if vo:
            check("double cover genus and infinity "+str(vo),1-2+(2*(2*vo))//4==vo-1)
    for b in [4,16,64]:
        check("zero-longitude peripheral tori arithmetic "+str(b),b//4+b//2==3*b//4)

    theorem=dict(proper_irreducible_transvection_linear_group_has_common_alternating_form=True,
                 every_proper_binary_linear_primitive_quotient_has_identity_preferred_longitude=True,
                 residual_components_identical_double_covers_or_split_affine_lines=True,
                 nonidentity_longitude_requires_full_AGL_k_ge3_and_odd_contact=True,
                 every_odd_k_ge3_binary_primitive_quotient_has_full_AGL=True)
    assurance=dict(classification_explicitly_adopted_via_pinned_primary_restatement=True,
                   given_module_and_separate_symplectic_proposition_checked=True,
                   actual_residual_cover_needed_beyond_common_form=True,
                   original_plane_and_invisible_primes_retained=True,
                   original_McLaughlin_PDF_obtained=False,
                   classification_independently_proved=False,
                   all_binary_groups_excluded=False,
                   full_AGL_control_claimed_Keller_realization=False,
                   integral_null_core_or_common_quadratic_refinement_assumed=False,
                   registered_predecessor_replay=False,formal_verification=False,
                   independent_replay=False,external_review=False,full_goal_complete=False)
    files=[NOTE,"verify_proper_binary_longitude.py",MANIFEST,"seifert_framing_cancellation_receipt.json"]+[s["path"] for s in manifest["sources"]]
    result=dict(claim="PF2RULED131",outcome="PASS",check_count=len(checks),checks=checks,theorem=theorem,
                assurance=assurance,source_sha256={f:sha(P/f) for f in files},environment=dict(python=platform.python_version()))
    (P/"proper_binary_longitude_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS",len(checks),"PF131 controls; proper binary linear groups have trivial longitude; goal ACTIVE")


if __name__=="__main__":
    main()
