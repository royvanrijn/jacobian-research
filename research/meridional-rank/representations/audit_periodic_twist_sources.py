#!/usr/bin/env python3
"""Record two scoped primary-source reads without modifying frozen manifests."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    records=[{
        'id':'STOIM07','title':'Graphs, determinants of knots and hyperbolic volume',
        'authors':'Alexander Stoimenow','year':2007,'doi':'10.2140/pjm.2007.232.423',
        'url':'https://msp.org/pjm/2007/232-2/pjm-v232-n2-p11-p.pdf',
        'snapshot':'sources/stoimenow-graphs-2007.pdf',
        'location':'Published PDF pp431,433: Lemmas3.14,3.21; pp446-447: Proposition6.3 and Remark6.4',
        'audit':'Primary published statements and their elementary spanning-tree arguments read. Exact region graphs and determinant counts independently checked here. Full volume estimates and the cited Thistlethwaite arborescent classification are not audited.',
        'scope':'Alternating-link determinant equals checkerboard spanning-tree count; deletion-contraction. The explicit6star exception prevents inferring knot nonarborescence from failure of series-parallel diagram reduction alone. Supplies no MRC or global novelty theorem.'
    },{
        'id':'JP25','title':'Meridional rank and bridge number of knotted 2-spheres',
        'authors':'Jason Joseph and Puttipong Pongtanapaisan','year':2025,
        'doi':'10.4153/S0008414X23000883','url':'https://doi.org/10.4153/S0008414X23000883',
        'published':'Canadian Journal of Mathematics77(1),282-299; online27Dec2023',
        'snapshot':'sources/joseph-pongtanapaisan-v2.html',
        'version_dates':'arXiv2111.09233v2,6Feb2023, under the title Bridge number and meridional rank of knotted surfaces; publication metadata checked separately',
        'location':'Preprint Section3, especially Construction3.5, Definition3.6 and Corollary3.12',
        'audit':'Classical construction, periodic twist-box rules and sharpness argument read. The finite-cycle generation-capacity theorem remains an explicit source import. Entire published proof and four-dimensional/welded extensions not independently audited.',
        'scope':'Specified classical families K_(p,n) have mr=bridge=n via symmetric/alternating quotients with meridians p-cycles. Periodic finite-quotient twist labellings have this precedent. No proof that our four seeds belong to or lie outside these families; knotted-sphere nonadditivity is not a classical MRC counterexample.'
    }]
    catalogue=json.loads((HERE/'sources.json').read_text());known={r['id']:r for r in catalogue['sources']}
    for r in records:
        if r['id'] in known:assert known[r['id']]==r
        else:catalogue['sources'].append(r)
    count=len(catalogue['sources'])
    assert len({r['id'] for r in catalogue['sources']})==count and count>=37
    (HERE/'sources.json').write_text(json.dumps(catalogue,indent=2,ensure_ascii=False)+'\n')
    paths=['audit_periodic_twist_sources.py','sources.json']+[r['snapshot'] for r in records]
    report={'status':'PASS scoped primary source records','records':records,'catalogue_distinct_ids':count,
            'supplementary_snapshots':2,'frozen_main_manifest_snapshots':18,
            'source_read_boundary':'Scoped source read, not a full-paper correctness audit or global priority review. Frozen earlier source manifests remain unchanged; these2 supplementary snapshots have their own seals.',
            'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'periodic_twist_source_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'new_source_ids':[r['id'] for r in records],'distinct_ids':count},indent=2))


if __name__=='__main__':main()
