#!/usr/bin/env python3
"""Record exact classical rigidity and published graph-cover family scopes."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    path=HERE/'sources.json';catalogue=json.loads(path.read_text());known={r['id']:r for r in catalogue['sources']}
    afw={'id':'AFW15','title':'3-Manifold Groups','authors':'Matthias Aschenbrenner, Stefan Friedl and Henry Wilton',
        'year':2015,'doi':'10.4171/154','url':'https://ems.press/books/elm/217',
        'snapshot':'sources/afw-3manifold-groups-final-author-2015.pdf',
        'snapshot_url':'https://www.maths.gla.ac.uk/~mpowell/Aschenbrenner-Friedl-Wilton-3-manifold-groups-final-version-031115.pdf',
        'comparison_snapshot':'sources/afw-3manifold-groups-author-v3.pdf',
        'comparison_url':'https://friedl.app.uni-regensburg.de/papers/3-manifold-groups-v3.pdf',
        'location':'Final author book184pp, printedpp29-31/PDF33-35, Theorems2.1.1-2.1.3 and Section2.2 boundary warning; compared earlier149pp author-v3 Theorems2.1-2.3, printedpp21-22. EMS publisher metadata gives20Aug2015, DOI10.4171/154,230printed pages.',
        'audit':'Read the prime closed rigidity theorem and full arbitrary-closed classification statement with lens/orientation exceptions and nearby attribution to geometrization, Mostow, Waldhausen and Scott. Read the Kneser converse and boundary-peripheral warning. Compared only these scopes across the two retained author versions. Full184/149page manuscripts, constituent rigidity proofs and full published-version comparison not audited; imported classical results.',
        'scope':'Closed orientable3manifolds with isomorphic groups differ only by orientations of prime factors and lens-space ambiguity. A closed orientable prime irreducible infinite-group graph-manifold model therefore determines any closed orientable manifold with its group up to homeomorphism and orientation; no boundary-manifold or finite lens-group rigidity shortcut.'}
    if 'AFW15' in known:assert known['AFW15']==afw
    else:catalogue['sources'].append(afw)
    bjw=known['BJW18'];bjw.update({'doi':'10.2140/pjm.2018.292.61',
        'snapshot':'sources/bjw-graph-cover-published-2018.pdf',
        'snapshot_url':'https://msp.org/pjm/2018/292-1/pjm-v292-n1-p03-s.pdf',
        'location':'Published PacificJ.Math292(1)(2018), printedpp61-64, especially p62 Theorem1.2, Cor1.3 and prime graph-manifold definition; proof architecture read, full proof imported.',
        'audit':'Published introduction and all relevant hypotheses re-read, including the graph-manifold definition as prime, and explicit treatment of both arborescent and non-arborescent branches. Cor1.3 assumes bridge exactly4; it does not by itself turn an upper4 into equality. Theorem1.2 plus the known meridional ranks1/2 equivalences gives MRC throughout this cover class with bridge at most4. Full rank3 inversion and branch-involution classification proof not independently audited.',
        'scope':'Graph-manifold double branched cover and meridional rank3 imply bridge3; exact bridge4 in that class gives meridional rank4. Hence this existing theorem covers the sixteen programme cases once actual-cover graph-manifold membership is proved by closed rigidity. This is distinct from the2017 theorem about graph-manifold knot exteriors. Prime nontrivial-JSJ knot exteriors have meridional rank at least4.'})
    count=len(catalogue['sources']);assert count==len({r['id'] for r in catalogue['sources']})==46
    path.write_text(json.dumps(catalogue,indent=2,ensure_ascii=False)+'\n')
    files=['audit_graph_cover_family_sources.py','graph_cover_family_scope_plan.json','sources.json',
        afw['snapshot'],afw['comparison_snapshot'],bjw['snapshot']]
    receipt={'status':'PASS scoped graph-cover family and closed-rigidity source audit','catalogue_distinct_ids':count,
        'new_ids':['AFW15'],'updated_ids':['BJW18'],'records':[afw,bjw],
        'rigidity_theorems':['AFW15 final-author Theorem2.1.2','AFW15 final-author Theorem2.1.3'],
        'known_family_theorem':'BJW18 published Corollary1.3; actual prime graph-manifold double branched cover and b=4',
        'unproved_by_sources_alone':'Actual named source-to-model group identities and individual upper/lower4 require the retained certificate proofs. Named global computational priority, all other known-family memberships and fullMRC are not settled by these source reads.',
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'graph_cover_family_source_audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'catalogue_distinct_ids':count,'new_ids':['AFW15']}))


if __name__=='__main__':main()
