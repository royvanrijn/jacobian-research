#!/usr/bin/env python3
"""Source-provenance rollup separating surfaces, fibrations and fibres."""
import argparse
from collections import Counter
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';OUT=ART/'parent_portfolio_coverage_v1.json'

def expected():
    index=ART/'new_high_rank_curve_index_v22.json';r17=ART/'compact_six_r17_atlas_v1.json';mw16=ART/'compact_five_mw16_atlas_v1.json';parents=ART/'mestre_parent_portfolio_intake_v1.json';proof=ART/'mestre_parent_and_label_independent_v1.json'
    data=cert.read(index);a=cert.read(r17);b=cert.read(mw16);new=cert.read(parents);verified=cert.read(proof)
    labels={r['family'] for r in a['families']}|{r['fibration_id'] for r in b['families']}|{'published-R17'}
    counts=Counter(r['family'] for r in data['curves']);assert set(counts)==labels and sum(counts.values())==201 and len(labels)==12
    paths=[Path(__file__).resolve(),index,r17,mw16,parents,proof,ROOT/'elkies-k3/AGENTS.md',ROOT/'elliptic-curves/notes/COMPACT_FIVE_MW16_ATLAS_2026-09-05.md']
    for r in a['families']:
        path=ROOT/r['source'];assert cert.hashed(path)==r['source_sha256'];paths.append(path)
    assert verified['status']==new['status']=='PASS' and len(new['rows'])==6 and len(new['pairwise_separation'])==21
    assert all(p['status']=='PROVED_NOT_Q_ISOMORPHIC_K3' and p['separating_good_primes'] for p in new['pairwise_separation'])
    assert [r['rank_lower_bound'] for r in new['rows']]==[r['rank_lower_bound'] for r in verified['parents']]==[11]*6
    return {'schema':'elliptic-curves.parent-portfolio-coverage.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'current_high_rank_inventory':{'curve_count':201,'family_label_counts':dict(sorted(counts.items())),'source_parent_surface':'X948','source_parent_count':1,'provenance_basis':'Existing direct fibration-hop and atlas source certificates; not inferred from equal NS determinant and not a new replay of every birational map.'},'new_parent_intake':{'count':6,'parent_ids':[r['id'] for r in new['rows']],'rank_lower_bounds':[11]*6,'pairwise_and_X948_separations':21,'count_including_current_X948':7,'point_search_boxes':0},'boundary':'Six additional verified Q-parent inputs, not six new production campaigns, geometric isomorphism classes, different-NS milestones, literature novelties or record curves. Family aliases, carrier lifts and coordinate changes do not increase the parent count.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();r=expected()
    if a.check:assert r==cert.read(OUT)
    else:
        if OUT.exists():raise FileExistsError('preserve parent coverage')
        checkpoint(OUT,r)
    print('PASS201 fibres from1 source parent;6 additional Q-distinct parent inputs')
