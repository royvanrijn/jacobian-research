#!/usr/bin/env python3
"""Verify the retained generic-only exposure of two recent outer26 discoveries."""
import argparse
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import json
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint,digest

ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
INDEX=ART/'new_high_rank_curve_index_v22.json';OUT=ART/'recent_outer26_exposure_audit_v1.json'
ROSTER=('new-20260906-189','new-20260906-192')


def expected():
    paths=set()
    def read(path):paths.add(path);return json.loads(path.read_text())
    inventory=read(INDEX);rows=[]
    for ident in ROSTER:
        entry=next(r for r in inventory['curves'] if r['id']==ident)
        source=read(ART/entry['source_certificate'])['curves'][entry['source_curve_index']]
        assert entry['rank_lower_bound']==source['rank_lower_bound']==26
        assert entry['parameter']==source['parameter']
        assert not source['icarm_matches'] and not source['previous_matches']
        assert not entry['current_catalogue_matches']
        raw_path=ROOT/source['discovery_witness']['path'];raw=read(raw_path)
        assert sha256(raw_path.read_bytes()).hexdigest()==source['discovery_witness']['sha256']
        assert raw['curve']==source['curve'] and raw['rank_lower_bound']==26
        assert raw['status']=='COMPLETE_DECLARED_POINT_ATTEMPT'
        points=source['points'];generic=source['generic_points'];m=len(generic)
        assert len(points)==26 and points[:m]==generic and m in (16,17)
        proof=source['rank_certificate']
        replay=checked_rank(tuple(map(Q,source['curve'])),[tuple(map(Q,p)) for p in points],
            [r['prime'] for r in proof['signatures']],proof['no_rational_2_torsion_prime'])
        assert digest(replay)==digest(proof)
        for c in raw['charts']:
            rep=c['centre']['representative'];search=c['search'];word=search['input']['centre']['coefficients']
            assert len(rep)==m and word[:m]==rep and not any(word[m:])
            assert search['status']=='bounded_search_complete' and search['height_bound']==125000 and search['timeout_seconds']==10
        expected_count=49 if m==17 else 43
        assert len(raw['charts'])==expected_count
        rows.append({'id':ident,'family':source['family'],'parameter':source['parameter'],
            'parameter_height':max(abs(Q(source['parameter']).numerator),Q(source['parameter']).denominator),
            'rank_lower_bound':26,'generic_prefix_dimension':m,'known_directions_beyond_generic':26-m,
            'completed_generic_center_boxes':expected_count,'source_search':str(raw_path.relative_to(ROOT)),
            'source_certificate':entry['source_certificate'],'source_curve_index':entry['source_curve_index'],
            'own26_center_exposure_in_this_source_run':0})
    paths.add(Path(__file__).resolve())
    return {'schema':'elliptic-curves.recent-outer26-exposure.v1','status':'PASS',
            'sources':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
            'rows':rows,'following_search':None,
            'scope':'The two retained source runs use only their16/17 generic directions as chart centres, despite ending with26 independent points. This proves an exposure difference available to an own26 policy; it is not an exhaustive repository-wide absence certificate, a rank-gain prediction, or an automatic following campaign.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');args=p.parse_args();r=expected()
    if args.check:assert r==json.loads(OUT.read_text())
    else:
        if OUT.exists():raise FileExistsError('preserve retained-exposure audit')
        checkpoint(OUT,r)
    print('PASS',r['rows'])
