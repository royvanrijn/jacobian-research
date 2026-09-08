#!/usr/bin/env python3
"""Exact specialized section-span ranks and known-subgroup index summary."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
import certify_endpoint_section_relations as relations
import verify_endpoint_section_lattice_indices as lattice
from research_runtime.store import checkpoint
ROOT=relations.ROOT;ART=relations.ART;OUT=ART/'endpoint_section_spans_summary_v1.json'

def expected():
    r=relations.expected();lattice.main();indices=cert.read(ART/'endpoint_section_lattice_indices_v1.json')
    if cert.read(relations.OUT)!=json.loads(json.dumps(r)):raise ArithmeticError('exact relation proof differs')
    rows=[]
    for a,b in zip(r['rows'],indices['rows']):
        rows.append({'family':a['family'],'endpoint':a['endpoint'],'specialized_generic_section_span_rank':a['generic_section_span_rank'],'original_generic_section_count':len(a['generic_points']),'known_section_index':b['index_of_certified_subset_in_known_section_lattice']})
    paths=[Path(__file__).resolve(),Path(relations.__file__).resolve(),Path(lattice.__file__).resolve(),relations.OUT,ART/'endpoint_section_lattice_indices_v1.json']
    return {'schema':'elliptic-curves.endpoint-section-spans-summary.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'rows':rows,'nonsingular_endpoints':len(rows),'strict_specialization_rank_losses':sum(r['specialized_generic_section_span_rank']<r['original_generic_section_count'] for r in rows),'proper_known_section_indices':sum(r['known_section_index']>1 for r in rows),'claim_boundary':'Ranks are exact for the rational spans of the transported generic sections, using independently certified subsets and exact integer group relations for every section value. Indices concern two known subgroups after quotienting torsion, and do not saturate the full Mordell-Weil group. All five proper indices are odd, so the inclusion of free lattices induces an isomorphism modulo two. This does not prove whole-curve rank, exclude exceptional directions, produce a new near-record curve or change the frozen endpoint point trial.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('endpoint section-span aggregate differs')
    else:
        if OUT.exists():raise FileExistsError('preserve section-span aggregate')
        checkpoint(OUT,d)
    print('EXACT21 ENDPOINT SECTION SPANS;',d['strict_specialization_rank_losses'],'STRICT LOSSES;',d['proper_known_section_indices'],'PROPER KNOWN SUBSETS',flush=True)
