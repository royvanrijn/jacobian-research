#!/usr/bin/env python3
"""Exact rank19 seed replay and post-score comparison with the pinned620 catalogue."""
import argparse
from pathlib import Path
import score_retained_native19 as scores
import audit_inventory200_current_catalogue as catalogue
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint,digest
ROOT,ART=scores.ROOT,scores.ART
OUT=ART/'retained_native19_intake_v1.json'

def expected():
    p=scores.protocol();score=scores.read(scores.OUT);source=scores.read(scores.SOURCE)
    metadata=catalogue.cert.read(catalogue.D/'metadata.json');db=catalogue.cert.read(catalogue.D/'database.json')
    assert scores.hashed(catalogue.D/'database.json')==metadata['sha256'] and len(db['curves'])==db['count']==620
    index=ART/'new_high_rank_curve_index_v22.json';inventory=scores.read(index)
    byj={};local={}
    for r in db['curves']:byj.setdefault(catalogue.j(r['ainvs']),[]).append(r)
    for r in inventory['curves']:local.setdefault(catalogue.j(r['curve']),[]).append(r)
    rows=[];seen=set()
    assert len(source['rows'])==len(score['rows'])==12
    for seed,scalar in zip(source['rows'],score['rows']):
        assert seed['curve']==scalar['model'] and seed['parameter']==scalar['parameter'] and seed['word']==scalar['word']
        proof=seed['rank_certificate'];model=tuple(map(catalogue.cert.F,seed['curve']))
        points=[tuple(map(catalogue.cert.F,P)) for P in seed['independent_points']]
        actual=checked_rank(model,points,[r['prime'] for r in proof['signatures']],proof['no_rational_2_torsion_prime'])
        assert digest(actual)==digest(proof) and actual['rank_lower_bound']==19
        j=catalogue.j(model);assert j not in seen;seen.add(j)
        matches=[r['id'] for r in byj.get(j,[]) if catalogue.cert.isomorphic(model,r['ainvs'])]
        previous=[r['id'] for r in local.get(j,[]) if catalogue.cert.isomorphic(model,r['curve'])]
        rows.append({'id':scalar['id'],'word':seed['word'],'parameter':seed['parameter'],'rank_lower_bound':19,
            'model_coefficient_bits':scalar['model_coefficient_bits'],'score_units':scalar['score_units'],
            'catalogue_matches':matches,'inventory_matches':previous})
    paths=[Path(__file__).resolve(),scores.SOURCE,scores.OUT,index,catalogue.D/'database.json',catalogue.D/'metadata.json',
           Path(catalogue.__file__),ROOT/'elliptic-curves/cas/memory_rank_certificate.py']
    return {'schema':'elliptic-curves.retained-native19-intake.v1','status':'PASS','rows':rows,
        'sources':{str(q.relative_to(ROOT)):scores.hashed(q) for q in paths},'catalogue_count':620,
        'inventory_count':201,'distinct_j_count':12,
        'scope':'All twelve finite19-point seed proofs replay. Exact j indexing and Q-isomorphism comparison with620 pinned catalogue equations and201 current inventory curves happen after scoring. This is relative novelty, not universal novelty, a rank increase, or point-search exposure. No catalogue field enters the score calculation.'}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');args=p.parse_args();r=expected()
    if args.check:assert r==scores.read(OUT)
    else:
        if OUT.exists():raise FileExistsError('preserve intake')
        checkpoint(OUT,r)
    print('PASS12 rank19 seeds; catalogue matches',sum(bool(x['catalogue_matches']) for x in r['rows']),'; inventory matches',sum(bool(x['inventory_matches']) for x in r['rows']))
