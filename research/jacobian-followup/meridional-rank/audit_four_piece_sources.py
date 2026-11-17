#!/usr/bin/env python3
"""Seal precise classical Seifert/JSJ and rank2 inputs, without novelty claims."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    additions=[{
        'id':'BW05','title':'The structure of 3-manifolds with two-generated fundamental group',
        'authors':'Michel Boileau and Richard Weidmann','year':2005,'doi':'10.1016/j.top.2004.10.008',
        'url':'https://doi.org/10.1016/j.top.2004.10.008','snapshot':'sources/boileau-weidmann-structure-author.pdf',
        'location':'Author manuscript42pp, Theorem2 printedp3 and all10necessary cases; Lemma2 printedp5; Section3.1 separating-torus argument printedpp16-20. Published Topology44(2)(2005),283-320; Crossref DOI metadata read.',
        'audit':'Read every necessary type in Theorem2 and the relevant separating/no-Klein-bundle case proof. The necessary list has at most3JSJpieces; original exceptional irregular-hyperbolic-cover case10 has only2 and does not affect this consequence. Millichap-Worden2016 Cor5.4 eliminates that original type, separately recorded as current literature. This classical theorem is imported, not independently reproved or formally certified. Full42page theorem proof and full published-version comparison not claimed; publication metadata checked via Crossref. A bounded erratum/gap search found no correction to the theorem, not an exhaustive historical audit.',
        'scope':'A compact orientable irreducible3manifold of rank2 with nontrivial JSJ cannot have4canonical JSJpieces. Used as a necessary condition only, without rank=Heegaard-genus or sufficiency/priority assumptions.'},
        {'id':'HAT00','title':'Notes on Basic 3-Manifold Topology','authors':'Allen Hatcher','year':2000,
        'url':'https://pi.math.cornell.edu/~hatcher/3M/3M.pdf','snapshot':'sources/hatcher-basic-3manifold.pdf',
        'location':'Online notes61pp, introduction dates small revisions1999/2000; Section1.2 Theorem1.9, Propositions1.11-1.13, Lemmas1.14-1.15 and proof of uniqueness pp14,16-23; Section2.1 explicit Seifert construction pp24-26',
        'audit':'Read Seifert circle-bundle/filling construction, irreducibility argument, boundary-fiber uniqueness and its exceptional list, and relevant JSJ uniqueness proof. Classical Jaco-Shalen/Johannson theory credited through this exposition; original monographs and all61pages not independently audited. Boundary-group injectivity in the actual four-piece construction is separately proved by explicit amalgam normal forms.',
        'scope':'Realize D(2,3) and A(2) pieces, glue incompressible tori, and certify canonical JSJ using nonparallel unique fiber slopes. Negative orbifold Euler characteristics exclude solid-torus/product/Klein-bundle exceptions. Does not identify source-cover homeomorphism from a group isomorphism; only a constructed manifold model with that group is needed.'},
        {'id':'MW16','title':'Hidden symmetries and commensurability of 2-bridge link complements',
        'authors':'Christian Millichap and William Worden','year':2016,'doi':'10.2140/pjm.2016.285.453',
        'url':'https://doi.org/10.2140/pjm.2016.285.453','snapshot':'sources/millichap-worden-2bridge-2016.pdf',
        'location':'Published PacificJ.Math285(2)(2016), Corollaries1.3 and5.3-5.4; printedpp455 and479',
        'audit':'Read revised rank2 characterization and its short deduction from the irregular-cover obstruction, including arithmetic/nonarithmetic cases and stated dependencies. Full canonical-cusp/hidden-symmetry proof not independently audited; treated as a published source for current literature context. The eight-case rank proof already uses the necessary JSJpiece bound from the older Boileau-Weidmann theorem and does not require this refinement.',
        'scope':'Eliminates the original exceptional irregular hyperbolic2bridge-cover case in the Boileau-Weidmann rank2 list. Records the current refinement instead of leaving the original conjectural type open; not a new result of this programme.'}]
    catalogue=json.loads((HERE/'sources.json').read_text());by_id={s['id']:s for s in catalogue['sources']}
    for record in additions:
        assert (HERE/record['snapshot']).is_file()
        if record['id'] in by_id:
            old=by_id[record['id']]
            changed={k for k in set(old)|set(record) if old.get(k)!=record.get(k)}
            assert not changed or (record['id']=='BW05' and changed=={'audit'})
            old.update(record)
        else:catalogue['sources'].append(record)
    count=len(catalogue['sources']);assert count==len({s['id'] for s in catalogue['sources']}) and count>=45
    (HERE/'sources.json').write_text(json.dumps(catalogue,indent=2,ensure_ascii=False)+'\n')
    files=['audit_four_piece_sources.py','sources.json','sources/boileau-weidmann-structure-author.pdf','sources/hatcher-basic-3manifold.pdf','sources/millichap-worden-2bridge-2016.pdf']
    receipt={'status':'PASS scoped four-piece source audit','catalogue_distinct_ids':count,'catalogue_at_original_checkpoint':45,'new_records':additions,
        'rank2_necessary_type_piece_counts':[2,2,2,2,2,3,2,1,1,2],
        'classical_rank2_maximum_JSJ_pieces':3,
        'assurance':'Source-derived necessary condition, not a new classification. New algebra/model inputs and matching uppers checked separately. No exhaustive literature, previously-unresolved worldwide status, or global priority conclusion.',
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'four_piece_source_audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'catalogue_distinct_ids':count,'new_ids':['BW05','HAT00','MW16']}))


if __name__=='__main__':main()
