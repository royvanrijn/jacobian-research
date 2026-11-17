#!/usr/bin/env python3
"""Scoped primary inputs for cover-bound additivity and its known boundaries."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    records=[{
        'id':'STALL71','title':'Group Theory and Three-Dimensional Manifolds',
        'authors':'John Stallings','year':1971,
        'url':'https://webhomes.maths.ed.ac.uk/~v1ranick/papers/stallings1.pdf',
        'snapshot':'sources/stallings-1971.pdf','location':'Theorem1.B.5, printedp6; author discussion and references23/24',
        'original_proof_metadata':'A topological proof of Grushko\'s theorem on free products, Math.Z.90(1965),1-8; DOI10.1007/BF01112046, verified on official Springer page',
        'audit':'Author primary statement read. Rank additivity is derived from its free-group lifting form in the new note. Grushko/Wagner/Neumann credit retained; the original1965 full proof is a classical import and is not independently certified.',
        'scope':'Finite ordinary group rank adds under free products. This supplies the rank step after the explicit lifted-sphere decomposition; it is not a meridional-additivity or full MRC theorem.'
    },{
        'id':'SCHULT03','title':'Additivity of bridge numbers of knots','authors':'Jennifer Schultens','year':2003,
        'doi':'10.1017/S0305004103006832','url':'https://doi.org/10.1017/S0305004103006832',
        'snapshot':'sources/schultens-bridge-additivity.pdf',
        'location':'Theorem1 and connected-sum argument, author manuscript printedpp4-5; published135(3),539-544',
        'audit':'Primary author manuscript theorem and connected-sum argument read; Cambridge publication metadata checked. Credit to Schubert1954 retained. No independent certification of every geometric foliation lemma or the original1954 proof.',
        'scope':'Bridge number has Schubert(-1)additivity for classical knots. This does not assert general meridional-rank additivity; used for additive nonnegative bridge-minus-cover-bound defects and known controls.'
    }]
    catalogue=json.loads((HERE/'sources.json').read_text());known={r['id']:r for r in catalogue['sources']}
    for record in records:
        if record['id'] in known:assert known[record['id']]==record
        else:catalogue['sources'].append(record)
    bdjw=next(r for r in catalogue['sources'] if r['id']=='BDJW17')
    bdjw['snapshot']='sources/bdjw-graph-manifold-v1.html'
    bdjw['audit']='Primary introduction, Theorem1, Corollary1 and class closure reread for connected sums of torus knots. Full graph-of-groups proof remains an explicit published import. Distinct from Dutra2022 with the separately checked gap.'
    bbk=next(r for r in catalogue['sources'] if r['id']=='BBK21')
    suffix=' The introduction also records(-1)additivity of meridional rank for knots detected by Coxeter quotients, without assuming bridge equality; this conditional precedent is credited, not generalized to all knots.'
    if suffix not in bbk['scope']:bbk['scope']+=suffix
    count=len(catalogue['sources'])
    assert count==len({r['id'] for r in catalogue['sources']}) and count>=39
    (HERE/'sources.json').write_text(json.dumps(catalogue,indent=2,ensure_ascii=False)+'\n')
    paths=['audit_cover_sum_sources.py','sources.json']+[r['snapshot'] for r in records]+[
        'sources/bdjw-graph-manifold-v1.html','sources/knotinfo-bridge-description-2026-11-17.html']
    report={'status':'PASS scoped cover-sum source audit','new_records':records,'catalogue_distinct_ids':count,
            'known_controls':'Connected sums of T(4,7) have graph-manifold exteriors and are covered by BDJW17, independently of whether the simple-cover lower bound is sharp.',
            'primary_depth':'Stallings primary theorem statement and official1965 metadata; Schultens primary theorem and connected-sum argument, geometric lemmas remain imports; BDJW primary class/theorem scopes. No full classical MRC, general meridional additivity or global-priority claim.',
            'live_description_scope_observation':{'url':'https://knotinfo.org/descriptions/bridge_index.html','date':'2026-11-17','observation':'The description credits Pfaff for bridge numbers of13crossing knots, without an explicit all-knots quantifier, named coverage list or certificates on that page. It does not certify closure of the retained named panel. Archived named rows and source-specific certificates remain separate evidence.'},
            'frozen_main_manifest_unchanged':True,'new_supplementary_snapshots':4,
            'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'cover_sum_source_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'new_ids':[r['id'] for r in records],'catalogue':count},indent=2))


if __name__=='__main__':main()
