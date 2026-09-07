#!/usr/bin/env python3
"""Bind six exact incidence cohorts to the current98 equations and duplicate map."""
import argparse
from pathlib import Path
from collections import Counter
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
NAMES=['compact','latest7','latest8','latest23','retention','extended20_mw16']
PROOFS=[ART/(n+'_cross_family_j_incidence_v1.json') for n in NAMES]
TRANSPORT=ART/'compact_published_r17_generic_transport_v1.json'
INDEX=ART/'new_high_rank_curve_index_v9.json'

def result():
    index=cert.read(INDEX);targets={r['id']:r for r in index['curves']};seen=set();counts=Counter();duplicates=[]
    if len(targets)!=98:raise ArithmeticError('fixed98 inventory differs')
    maps=None
    for path in PROOFS:
        data=cert.read(path)
        if data['status']!='COMPLETE_DECLARED_INCIDENCE_AUDIT':raise ArithmeticError('incomplete incidence cohort')
        current={r['family']:(r['A'],r['B']) for r in data['maps']}
        if maps is None:maps=current
        if maps!=current or len(maps)!=12:raise ArithmeticError('presentation equations changed across cohorts')
        for r in data['targets']:
            if r['id'] in seen or any(r[k]!=targets[r['id']][k] for k in r if k!='rank_lower_bound') or r['rank_lower_bound']>targets[r['id']]['rank_lower_bound']:raise ArithmeticError('cohort overlap, target model mismatch or rank regression')
            seen.add(r['id'])
        if len(data['pairs'])!=12*len(data['targets']) or len({(r['target'],r['family']) for r in data['pairs']})!=len(data['pairs']):raise ArithmeticError('cohort product differs')
        for pair in data['pairs']:
            counts[pair['status']]+=1
            if pair['status']=='NO_RATIONAL_J_PREIMAGE':continue
            if pair['status']!='RATIONAL_J_PREIMAGES_CERTIFIED' or pair['infinity'] or len(pair['rational_roots'])!=1:raise ArithmeticError('unresolved or extra incidence')
            target=targets[pair['target']];parameter=cert.F(pair['rational_roots'][0]['parameter'])
            if pair['family']==target['family']:
                if parameter!=cert.F(target['parameter']):raise ArithmeticError('own parameter changed')
            else:
                if {pair['family'],target['family']}!={'08234','published-R17'}:raise ArithmeticError('new presentation needs a separate transport proof')
                compact=parameter if pair['family']=='08234' else cert.F(target['parameter'])
                published=parameter if pair['family']=='published-R17' else cert.F(target['parameter'])
                if compact!=-26*published-50:raise ArithmeticError('duplicate base map differs')
                duplicates.append({'id':target['id'],'compact_parameter':str(compact),'published_parameter':str(published)})
    if seen!=set(targets) or counts!={'NO_RATIONAL_J_PREIMAGE':1057,'RATIONAL_J_PREIMAGES_CERTIFIED':119} or len(duplicates)!=21:raise ArithmeticError('fixed1176 coverage differs')
    transport=cert.read(TRANSPORT)
    if transport['status']!='PASS_EXACT_GENERIC_TRANSPORT' or transport['matrix_determinant']!=-1 or transport['compact_model_scale_u']!='26':raise ArithmeticError('exact generic transport missing')
    sources={str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),INDEX,TRANSPORT,*PROOFS]}
    return {'schema':'elliptic-curves.inventory98-incidence.v2','status':'PASS','sources':sources,'targets_checked':98,'pairs_checked':1176,'status_counts':dict(counts),'same_generic_subgroup_duplicates':duplicates,'historical_rank_metadata':'Incidence targets retain their historical lower bounds; only nondecreasing current bounds are accepted, and all curve, j, family and parameter fields must remain identical. Incidence itself is independent of rank.','claim_boundary':'Combines the six separately exact-replayed rational j-incidence cohorts and the exact Q(t) unimodular17-section transport. Every additional recorded presentation is the same generic subgroup. The separate incidence/transport replays remain required; this summary does not replay ranks or exclude other families, nongeneric points or higher ranks.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path);p.add_argument('--check',type=Path);a=p.parse_args();r=result()
    if a.check:
        if cert.read(a.check)!=r:raise ArithmeticError('inventory incidence summary differs')
    else:
        if a.output is None or a.output.exists():raise FileExistsError('new output required')
        checkpoint(a.output,r)
    print('EXACT COHORT BINDING98 CURVES1176 PAIRS21 DUPLICATE PRESENTATIONS',flush=True)
