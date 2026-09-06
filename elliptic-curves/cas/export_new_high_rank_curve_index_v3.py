#!/usr/bin/env python3
"""Extend the exact curve inventory while preserving every existing local curve ID."""
import argparse
import csv
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OLD=ART/'new_high_rank_curve_index_v2.json'
ADDITIONAL=ART/'compact_r17_wide_results_v1.json'

def sources():
    paths=(Path(__file__).resolve(),Path(cert.__file__).resolve(),ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',ROOT/'elliptic-curves/cas/elliptic_candidate_record.py')
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def expected_rows():
    previous=cert.read(OLD)
    rows=[dict(r) for r in previous['curves']]
    additions=[]
    for i,r in enumerate(cert.read(ADDITIONAL)['curves']):
        rank=r['rank_certificate']['rank_lower_bound']
        if rank<22 or r['icarm_matches']:
            continue
        model=tuple(map(cert.F,r['curve']))
        match=next((q for q in rows if cert.isomorphic(model,q['curve'])),None)
        if match is not None and rank<=match['rank_lower_bound']:
            continue
        inv=cert.weierstrass_invariants(model)
        row={'family':r['family'],'parameter':r['parameter'],'rank_lower_bound':rank,
            'curve':r['curve'],'points':r['points'],'rank_certificate':r['rank_certificate'],
            'j_invariant':str(inv['c4']**3/inv['discriminant']),
            'source_certificate':ADDITIONAL.name,'source_curve_index':i}
        if match is not None:
            row['id']=match['id'];rows[rows.index(match)]=row
        else:
            additions.append(row)
    additions.sort(key=lambda r:(-r['rank_lower_bound'],r['family'],r['parameter']))
    for i,row in enumerate(additions):
        row['id']=f"new-20260905-{len(previous['curves'])+i+1:02}"
    rows.extend(additions)
    return sorted(rows,key=lambda r:(-r['rank_lower_bound'],r['family'],r['parameter'])),previous

def build(output):
    if output.exists() or output.with_suffix('.csv').exists():
        raise FileExistsError('preserve curve index')
    rows,previous=expected_rows()
    names=set(previous['source_certificate_hashes'])|{ADDITIONAL.name}
    cert.write(output,{'schema':'elliptic-curves.new-high-rank-index.v3','sources':sources(),
        'previous_index_sha256':cert.hashed(OLD),'curves':rows,
        'source_certificate_hashes':{name:cert.hashed(ART/name) for name in sorted(names)},
        'catalogue':cert.read(ADDITIONAL)['catalogue'],
        'claim_boundary':'Exact independent-point lower bounds and absence from the pinned 586-equation catalogue. Existing local curve IDs are preserved. These are neither exact ranks nor a proof of universal novelty; no new rank-28/32 curve is claimed.'})
    with output.with_suffix('.csv').open('x',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['id','rank_lower_bound','family','parameter','a1','a2','a3','a4','a6','source_certificate','source_curve_index'])
        for r in rows:
            writer.writerow([r['id'],r['rank_lower_bound'],r['family'],r['parameter'],*r['curve'],r['source_certificate'],r['source_curve_index']])

def check(path):
    data=cert.read(path)
    if data['sources']!=sources() or data['previous_index_sha256']!=cert.hashed(OLD):
        raise ArithmeticError('index sources changed')
    expected,previous=expected_rows()
    if data['curves']!=expected or data['catalogue']!=cert.read(ADDITIONAL)['catalogue']:
        raise ArithmeticError('index extraction or stable IDs differ')
    inputs={}
    for name,h in data['source_certificate_hashes'].items():
        if cert.hashed(ART/name)!=h:
            raise ArithmeticError('source certificate changed')
        inputs[name]=cert.read(ART/name)['curves']
    seen=set()
    for r in data['curves']:
        source=inputs[r['source_certificate']][r['source_curve_index']]
        if any(r[k]!=source[k] for k in ('parameter','curve','points','rank_certificate')):
            raise ArithmeticError('source point proof differs')
        model=tuple(map(cert.F,r['curve']));points=[tuple(map(cert.F,p)) for p in r['points']]
        old=r['rank_certificate']
        actual=cert.checked_rank(model,points,[s['prime'] for s in old['signatures']],old['no_rational_2_torsion_prime'])
        if json.dumps(old,sort_keys=True)!=json.dumps(actual,sort_keys=True) or r['rank_lower_bound']!=len(points):
            raise ArithmeticError('exact rank proof differs')
        inv=cert.weierstrass_invariants(model);j=inv['c4']**3/inv['discriminant']
        if str(j)!=r['j_invariant'] or j in seen:
            raise ArithmeticError('j equality or duplicate needs separate treatment')
        if any(cert.isomorphic(model,q['ainvs']) for q in data['catalogue']['equations']):
            raise ArithmeticError('catalogue match in new inventory')
        seen.add(j)
        print('REPLAYED V3 INVENTORY',r['id'],'rank >=',len(points),flush=True)
    print('EXACT CURVE INVENTORY',len(seen),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ART/'new_high_rank_curve_index_v3.json');p.add_argument('--check',type=Path)
    a=p.parse_args();check(a.check) if a.check else build(a.output)
