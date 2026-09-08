#!/usr/bin/env python3
"""Post-discovery exact equation comparison with a separately downloaded catalogue."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/inventory200-current-catalogue-v1'
INDEX=ART/'new_high_rank_curve_index_v20.json'
OUT=ART/'inventory200_current_catalogue_comparison_v1.json'


def j(model):
    inv=cert.weierstrass_invariants(tuple(map(cert.F,model)))
    if not inv['discriminant']:raise ArithmeticError('singular equation')
    return inv['c4']**3/inv['discriminant']


def result():
    metadata=cert.read(D/'metadata.json');database=cert.read(D/'database.json')
    inventory=cert.read(INDEX);prior=inventory['catalogue']['equations']
    if cert.hashed(D/'database.json')!=metadata['sha256'] or database['count']!=620 or len(database['curves'])!=620:
        raise ArithmeticError('frozen620-equation download differs')
    if len(inventory['curves'])!=200 or len(prior)!=593:
        raise ArithmeticError('fixed inventory and historical comparison required')
    equations=[{'id':r['id'],'ainvs':list(map(str,r['ainvs']))} for r in database['curves']]
    if len({r['id'] for r in equations})!=620:raise ArithmeticError('duplicate catalogue IDs')
    byj={}
    for r in equations:byj.setdefault(j(r['ainvs']),[]).append(r)
    rows=[]
    for r in inventory['curves']:
        value=j(r['curve']);samej=byj.get(value,[])
        rows.append({'id':r['id'],'curve':r['curve'],'j_invariant':str(value),
                     'same_j_ids':[q['id'] for q in samej],
                     'q_isomorphism_matches':[q['id'] for q in samej if cert.isomorphic(r['curve'],q['ainvs'])]})
    old_ids={r['id'] for r in prior}
    historical=[{'id':r['id'],'identical_equation_ids':[q['id'] for q in byj.get(j(r['ainvs']),[])
                 if tuple(map(cert.F,r['ainvs']))==tuple(map(cert.F,q['ainvs']))]} for r in prior]
    paths=[Path(__file__).resolve(),Path(cert.__file__),ROOT/'elliptic-curves/cas/elliptic_candidate_record.py',INDEX,D/'database.json',D/'metadata.json']
    return {'schema':'elliptic-curves.inventory200-current-catalogue.v1','status':'PASS',
            'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
            'catalogue_intake':metadata,'catalogue_count':620,'equations':equations,
            'inventory_count':200,'inventory_comparisons':rows,'historical_equation_comparisons':historical,
            'added_catalogue_ids':sorted(r['id'] for r in equations if r['id'] not in old_ids),
            'matched_inventory_ids':[r['id'] for r in rows if r['q_isomorphism_matches']],
            'claim_boundary':'Post-discovery exact Q-isomorphism comparisons of the complete200-curve V20 inventory with620 equations in this hash-pinned download. No search, candidate selection, point import, rank recertification, conductor calculation or universal-novelty claim. Existing593-equation proofs remain unchanged; any newly found match must be reported rather than silently removed.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=result()
    if a.check:
        if cert.read(OUT)!=d:raise ArithmeticError('catalogue replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve fresh comparison')
        checkpoint(OUT,d)
    print('INVENTORY200 VS620 PASS; MATCHED',d['matched_inventory_ids'])
