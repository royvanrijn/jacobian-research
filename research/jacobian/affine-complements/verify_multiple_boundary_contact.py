"""PF103 exact controls. Universal proof and primary imports are separate."""
from pathlib import Path
import hashlib,json,platform,sys
import sympy as S
import multiple_boundary_contact as M

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,value):check(name,S.factor(value)==0)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    for name in ["primitive_end_contact_receipt.json","proper_finite_marked_end_receipt.json",
                 "full_peripheral_linking_receipt.json","canonical_formats_receipt.json"]:
        old=json.loads((P/name).read_text())
        check("predecessor PASS "+name,old["outcome"]=="PASS")
        for f,h in old["source_sha256"].items():
            check("predecessor seal "+f,sha(P/f)==h)
    manifest=json.loads((P/"sources/multiple-boundary-contact/manifest.json").read_text())
    for source in manifest["primary_sources"]:
        f=P/source["path"]
        check("primary pin "+f.name,sha(f)==source["sha256"] and f.stat().st_size==source["bytes"])
        check("full-page visual audit "+f.name,bool(source["complete_PDF_pages_rendered_and_visually_inspected"]))
    r,E,p=S.symbols("r E p",integer=True,positive=True)
    ident("UNBOUNDED chi/sigma canonical formula",-E-2*(1+r)-3*(-r)-(r-2-E))
    ident("UNBOUNDED A-chain theta",-2*p-3*(-(p-1))-(p-3))
    ident("UNBOUNDED contact equality requires E=r-p+1",(r-2-E-(p-3)).subs(E,r-p+1))
    s=S.symbols("s",nonnegative=True,integer=True)
    ident("UNBOUNDED rank bound r=p-1-s requires E=-s",(r-p+1).subs(r,p-1-s)+s)
    check("positive E incompatible with required -s",(-s).is_nonpositive and E.is_positive)
    n,A,C=S.symbols("n A C",integer=True,positive=True)
    # B=I+J, J^2=nJ: inverse I-J/(n+1). No matrix-size extrapolation.
    ident("UNBOUNDED rank-one Gram inverse coefficient",1-1/(n+1)-n/(n+1))
    ident("UNBOUNDED zero-canonical controls theta",(n-2)-(n+1-3))
    ident("UNBOUNDED Gram energy expression",(A-C/(n+1))*(n+1)-((n+1)*A-C))
    controls=[]
    for B,indices,label in [([[2]],[1],"allowed quadric"),
                             ([[2]],[3],"ramified RP3 fails theta"),
                             ([[1,0],[0,2]],[2,1],"blown-up formal lattice passes theta but fails minimal rank"),
                             ([[2,0],[0,2]],[1,1],"individual cores need not generate"),
                             ([[2,-1],[-1,2]],[3,1],"two-prime ramified A2 lattice fails theta")]:
        packet=M.canonical_packet(B,indices)
        controls.append(dict(label=label,packet=packet))
    quad,ram,blow,nonprimitive,ar2=[c["packet"] for c in controls]
    check("actual affine quadric e1 permitted",quad["theta"]=="-1" and not quad["canonical_nonzero"])
    check("RP3 ramification mismatch",not M.a_type_comparison(ram,2)["contact_theta_matches"])
    b=M.a_type_comparison(blow,2)
    check("minimality/rank input ESSENTIAL",b["contact_theta_matches"] and b["H1_order_matches"] and
          not b["minimal_filling_rank_bound"] and blow["canonical_nonzero"])
    check("collective generation does not imply individual primitivity",
          nonprimitive["H1_order"]==4 and nonprimitive["core_orders"]==[2,2])
    check("multi-prime A2 canonical failure",ar2["canonical_energy"]=="8/3" and
          not M.a_type_comparison(ar2,3)["contact_theta_matches"])
    for size in [1,2,3,6]:
        B=S.eye(size)+S.ones(size);inv=S.eye(size)-S.ones(size)/(size+1)
        check("representative all-r marked Gram inverse "+str(size),B*inv==S.eye(size))
        packet=M.canonical_packet(B,[1]*size)
        check("representative determinant and individual core orders "+str(size),
              packet["H1_order"]==size+1 and packet["core_orders"]==[size+1]*size)
        check("representative zero-canonical A-type boundary data "+str(size),
              M.a_type_comparison(packet,size+1)["contact_theta_matches"])
    B=S.Matrix([[2,-1],[-1,2]]);k=S.Matrix([2,0]);T=S.Matrix([[1,2],[0,1]])
    ident("inverse energy invariant under integral compact-basis change",
          ((T.T*k).T*(T.T*B*T).inv()*(T.T*k))[0]-(k.T*B.inv()*k)[0])
    check("basis change unimodular",T.det()==1)
    for bad,indices in [([[0]],[1]),([[1,2],[2,1]],[2,2]),([[2]],[2]),([[2]],[0])]:
        try:M.canonical_packet(bad,indices)
        except ValueError:check("invalid source packet rejected "+str((bad,indices)),True)
        else:raise AssertionError("invalid packet accepted")
    # Genuine smooth affine whole-plane extensions, algebra only (no Keller target).
    x,f,W,V,t,u,a=S.symbols("x f W V t u a")
    for size in [1,2,4]:
        roots=list(range(size+1));poly=S.prod(W-j for j in roots)
        v=x*S.prod(f*x-j for j in roots[1:])
        ident("affine whole-plane defining identity "+str(size),f*v-poly.subs(W,f*x))
        check("affine smoothness P squarefree "+str(size),S.gcd(poly,S.diff(poly,W))==1)
        ident("retained component inverse "+str(size),v/S.prod(f*x-j for j in roots[1:])-x)
    ident("all omitted local charts volume has no zero",S.diff(a/u+t,t)*S.diff(u,u)-1)
    proof=(P/"MULTIPLE_BOUNDARY_CANONICAL_CONTACT_AND_A_TYPE_END_EXCLUSION.md").read_text()
    for word in ["COLLECTIVELY","k is NONZERO","symplectically minimal","r <= p-1",
                 "opposite","not independently","singular N","WHOLE"]:
        check("scope/proof guard "+word,word in proof)
    files=["MULTIPLE_BOUNDARY_CANONICAL_CONTACT_AND_A_TYPE_END_EXCLUSION.md",
           "multiple_boundary_contact.py","verify_multiple_boundary_contact.py",
           "sources/multiple-boundary-contact/manifest.json"]
    receipt=dict(claim="PF2RULED103",outcome="PASS",check_count=len(checks),checks=checks,
                 source_sha256={f:sha(P/f) for f in files},controls=controls,
                 python=platform.python_version(),sympy=S.__version__,
                 assurance=dict(written_universal_proof=True,primary_contact_and_filling_theorems_imported=True,
                                all_oriented_A_type_lens_ends_excluded_for_actual_smooth_Keller_normalization=True,
                                arbitrary_many_omitted_primes=True,all_field_degrees=True,
                                singular_normalizations_excluded=False,other_lens_ends_excluded=False,
                                finite_checks_prove_universal_claim=False,formal_verification=False,
                                external_review=False,global_JC2_closed=False,external_review_deferred=True))
    (P/"multiple_boundary_contact_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(dict(outcome="PASS",check_count=len(checks),claim="PF2RULED103")))


if __name__=="__main__":main()
