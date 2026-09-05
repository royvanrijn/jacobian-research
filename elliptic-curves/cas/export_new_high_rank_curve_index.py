#!/usr/bin/env python3
"""Consolidate the32 new rank>=22 curves and replay their exact point proofs."""
import argparse
import csv
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
NAMES=('compact_r17_new_curves_v1.json','compact_r17_wide_new_curves_v1.json','compact_r17_top64_interim_curves_v1.json',
       'compact_r17_largest_gain_curve_v1.json','compact_atlas_new_curves_v1.json','prospective_mw16_results_v1.json','prospective_mw16_wide_results_v1.json')

def sources():
    paths=(Path(__file__).resolve(),Path(cert.__file__).resolve(),ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',ROOT/'elliptic-curves/cas/elliptic_candidate_record.py')
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def build(output):
    if output.exists() or output.with_suffix('.csv').exists():raise FileExistsError('preserve curve index')
    rows=[]
    for name in NAMES:
        for i,r in enumerate(cert.read(ART/name)['curves']):
            rank=r['rank_certificate']['rank_lower_bound']
            if rank<22 or r.get('icarm_matches') or r.get('previous_matches'):continue
            model=tuple(map(cert.F,r['curve']));inv=cert.weierstrass_invariants(model)
            rows.append({'family':r.get('family','published-R17'),'parameter':r['parameter'],'rank_lower_bound':rank,
                'curve':r['curve'],'points':r['points'],'rank_certificate':r['rank_certificate'],
                'j_invariant':str(inv['c4']**3/inv['discriminant']),'source_certificate':name,'source_curve_index':i})
    if len(rows)!=32:raise ArithmeticError('frozen32-curve roster changed')
    rows.sort(key=lambda r:(-r['rank_lower_bound'],r['family'],r['parameter']))
    for i,r in enumerate(rows):r['id']=f'new-20260905-{i+1:02}'
    projection=[{'id':r['id'],'ainvs':r['ainvs']} for r in cert.read(cert.DATABASE)['curves']]
    cert.write(output,{'schema':'elliptic-curves.new-high-rank-index.v1','sources':sources(),'curves':rows,
        'source_certificate_hashes':{name:cert.hashed(ART/name) for name in NAMES},
        'catalogue':{'date':'2026-09-05','raw_sha256':cert.hashed(cert.DATABASE),'equations':projection,'acknowledgement':'ICARM, supported by NSF Grant DMS2425401'},
        'claim_boundary':'32 distinct exact j-invariants,32 independent point certificates of ranks at least22–25, and no Q-isomorphism in the pinned584-equation catalogue. Labels are local inventory IDs. No universal novelty, exact rank, conductor record, rank>=28 or rank>=32 claim.'})
    with output.with_suffix('.csv').open('x',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['id','rank_lower_bound','family','parameter','a1','a2','a3','a4','a6','source_certificate','source_curve_index'])
        for r in rows:writer.writerow([r['id'],r['rank_lower_bound'],r['family'],r['parameter'],*r['curve'],r['source_certificate'],r['source_curve_index']])

def check(path):
    data=cert.read(path)
    if data['sources']!=sources():raise ArithmeticError('index checker sources changed')
    inputs={}
    for name,h in data['source_certificate_hashes'].items():
        if cert.hashed(ART/name)!=h:raise ArithmeticError('source certificate changed')
        inputs[name]=cert.read(ART/name)['curves']
    seen=set()
    for r in data['curves']:
        source=inputs[r['source_certificate']][r['source_curve_index']]
        if any(r[k]!=source[k] for k in ('parameter','curve','points','rank_certificate')) or r['family']!=source.get('family','published-R17'):raise ArithmeticError('index extraction changed')
        model=tuple(map(cert.F,r['curve']));points=tuple(tuple(map(cert.F,p)) for p in r['points']);proof=r['rank_certificate']
        actual=cert.checked_rank(model,points,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
        if json.dumps(actual,sort_keys=True)!=json.dumps(proof,sort_keys=True) or r['rank_lower_bound']!=len(points):raise ArithmeticError('rank certificate changed')
        inv=cert.weierstrass_invariants(model);j=inv['c4']**3/inv['discriminant']
        if str(j)!=r['j_invariant'] or j in seen:raise ArithmeticError('j-invariant equality or duplicate')
        if any(cert.isomorphic(model,q['ainvs']) for q in data['catalogue']['equations']):raise ArithmeticError('catalogue match in new index')
        seen.add(j);print('REPLAYED',r['id'],'rank >=',len(points),flush=True)
    if len(seen)!=32:raise ArithmeticError('index roster changed')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ART/'new_high_rank_curve_index_v1.json');p.add_argument('--check',type=Path);a=p.parse_args()
    check(a.check) if a.check else build(a.output)
