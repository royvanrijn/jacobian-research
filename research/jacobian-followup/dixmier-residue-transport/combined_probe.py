#!/usr/bin/env python3
"""Build the shear-invariant monomial envelope before seeking a strong dual."""
import json
import time
from pathlib import Path
import verify_transport as a
from explore_degree_five_quantum_residue import split_pair, ThirdOrderFamily, fifth_order_coefficients

HERE=Path(__file__).resolve().parent

def orbit_support(bound,zcap):
    native=a.laurent_monomials(bound,zcap,0,3)
    support=set(native)
    for i,j,k in native:
        support.update((i+4*r,j,k-r) for r in range(k+1))
    return sorted(support)

def main():
    start=time.monotonic()
    S,T=a.degree_five_family(a.QQ,-a.QQ(1)/2,a.QQ(1))
    source_family=a.third_order_family(S,T,a.QQ)
    sm,tm=orbit_support(25,3),orbit_support(21,2)
    u,v=orbit_support(21,1),orbit_support(17,0)
    # Retained three-variable polynomial functions, with no lambda sampling.
    from explore_degree_five_a2_subprincipal import poisson,pi_power,scale
    cols=[poisson({m:a.QQ.one},T) for m in sm]+[poisson(S,{m:a.QQ.one}) for m in tm]
    rhs=scale(pi_power(S,T,3),-a.QQ.one/24)
    base,kernel,rank=a.solve_affine(cols,rhs,a.QQ)
    packet={"status":"LOWER_SPACE_CONSTRUCTED","scope":"Full polynomial monomial envelope of both native and transported supports; stronger than their two-space sum.",
            "lower_support_counts":[len(sm),len(tm)],"current_support_counts":[len(u),len(v)],
            "lower_rank":rank,"lower_dimension":len(kernel),"elapsed_seconds":time.monotonic()-start}
    print(json.dumps(packet),flush=True)
    (HERE/"combined_probe.json").write_text(json.dumps(packet,indent=2)+"\n")
    fam=ThirdOrderFamily(split_pair(base,sm,tm),[split_pair(k,sm,tm) for k in kernel],rank,len(cols))
    O,variations=fifth_order_coefficients(S,T,fam,a.QQ)
    current=[poisson({m:a.QQ.one},T) for m in u]+[poisson(S,{m:a.QQ.one}) for m in v]
    strong=current+variations
    packet.update({"variation_columns":len(variations),"strong_columns":len(strong),"output_monomials":len(set(O).union(*(set(p) for p in strong))),"elapsed_seconds":time.monotonic()-start})
    print("ASSEMBLED",json.dumps(packet),flush=True)
    (HERE/"combined_probe.json").write_text(json.dumps(packet,indent=2)+"\n")
    ranks=[a.column_rank(strong),a.column_rank(strong+[O])]
    packet["strong_ranks"]=ranks
    print("STRONG_RANKS",ranks,flush=True)
    if ranks[1]>ranks[0]:
        mons=sorted(set(O).union(*(set(p) for p in strong)))
        equations=strong+[O]
        transpose=[{j:p[m] for j,p in enumerate(equations) if m in p} for m in mons]
        dual,_,_=a.solve_affine(transpose,{len(strong):a.QQ.one},a.QQ)
        functional={mons[i]:c for i,c in dual.items()}
        for p in strong: assert sum(c*functional.get(m,a.QQ.zero) for m,c in p.items())==0
        assert sum(c*functional.get(m,a.QQ.zero) for m,c in O.items())==1
        packet["status"]="STRONG_COCYCLE_FOUND"
        packet["functional"]=[[list(m),str(c)] for m,c in sorted(functional.items())]
        print("DUAL_TERMS",len(functional),flush=True)
    else:
        packet["status"]="NO_NONZERO_STRONG_PERIOD_ON_ENVELOPE"
    packet["elapsed_seconds"]=time.monotonic()-start
    (HERE/"combined_probe.json").write_text(json.dumps(packet,indent=2)+"\n")

if __name__=="__main__": main()
