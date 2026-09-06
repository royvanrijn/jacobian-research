#!/usr/bin/env python3
"""Extend the exact curve inventory while preserving every existing local curve ID."""
import argparse
import csv
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OLD=ART/'new_high_rank_curve_index_v14.json'
ADDITIONAL=ART/'full11952_late64_r17_results_v1.json'

def sources():
    paths=(Path(__file__).resolve(),Path(cert.__file__).resolve(),ROOT/'elliptic-curves/cas/memory_rank_certificate.py',ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',ROOT/'elliptic-curves/cas/elliptic_candidate_record.py',ART/'full11952_late64_point_portable_replay_v1.json')
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def promotion_gate():
    base=ROOT/'artifacts/local/elliptic-curves/full11952-late64-r17-pari-v1'
    portable=cert.read(ART/'full11952_late64_point_portable_replay_v1.json')
    if portable['status']!='PASS' or portable['logical_stages']!=194:
        raise ArithmeticError('all194 isolated point-proof stages required before inventory promotion')
    if any(cert.hashed(ROOT/n)!=h for n,h in portable['sources'].items()):
        raise ArithmeticError('isolated proof source or ledger binding changed')
    protocol=cert.read(base/'post-batch/protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in protocol['sources'].items()):
        raise ArithmeticError('fixed full11952_late64 finalization sources differ')
    ledger=cert.read(base/'post-batch/ledger.json')
    stages=['certify','proof-replay','aggregate-build','aggregate-check']
    if ledger['status']!='PASS' or ledger.get('stronger_odd_prime_bounds')!=[] or [r['name'] for r in ledger['rows']]!=stages:
        raise ArithmeticError('complete full11952_late64 proof and review of any stronger odd bounds required')
    for row in ledger['rows']:
        supervision=row['supervision']
        log=base/'post-batch'/(row['name']+'.log')
        if row['status']!='PASS' or supervision['outcome']!='completed' or supervision['returncode']!=0 or cert.hashed(log)!=supervision['log_sha256']:
            raise ArithmeticError('successful finalization transcript differs')
    aggregate=cert.read(ART/'full11952_late64_experiment_v1.json')
    if aggregate['status']!='PASS' or len(aggregate['rows'])!=64:
        raise ArithmeticError('complete64 aggregate required')
    if any(cert.hashed(ROOT/n)!=h for n,h in aggregate['sources'].items()):
        raise ArithmeticError('full11952_late64 aggregate inputs changed')
    incoming=cert.read(ADDITIONAL)
    if incoming['within_batch_isomorphic_pairs'] or len(incoming['curves'])!=64:
        raise ArithmeticError('complete distinct64 point proof required')
    for row, proof in zip(aggregate['rows'],incoming['curves']):
        if row['id']!=proof['id'] or row['rank_lower_bound']!=proof['rank_lower_bound'] or max(row['odd_modulus_lower_bounds'].values())>row['rank_lower_bound']:
            raise ArithmeticError('unreviewed stronger odd bound or cohort mismatch')
        if row['previous_matches'] or proof['previous_matches']:
            raise ArithmeticError('unexpected previously measured equation requires review')

def expected_rows():
    promotion_gate()
    previous=cert.read(OLD)
    rows=[dict(r) for r in previous['curves']]
    additions=[]
    incoming=list(enumerate(cert.read(ADDITIONAL)['curves']))
    incoming.sort(key=lambda t:(-t[1]['rank_certificate']['rank_lower_bound'],t[1]['family'],t[1]['parameter']))
    for i,r in incoming:
        rank=r['rank_certificate']['rank_lower_bound']
        if rank<22 or r['icarm_matches']:continue
        model=tuple(map(cert.F,r['curve']))
        match=next((q for q in rows+additions if cert.isomorphic(model,q['curve'])),None)
        if match is not None and rank<=match['rank_lower_bound']:continue
        inv=cert.weierstrass_invariants(model)
        row={'family':r['family'],'parameter':r['parameter'],'rank_lower_bound':rank,'curve':r['curve'],'points':r['points'],'rank_certificate':r['rank_certificate'],'j_invariant':str(inv['c4']**3/inv['discriminant']),'source_certificate':ADDITIONAL.name,'source_curve_index':i}
        if match is not None:
            row['id']=match.get('id');destination=rows if match in rows else additions;destination[destination.index(match)]=row
        else:additions.append(row)

    additions.sort(key=lambda r:(-r['rank_lower_bound'],r['family'],r['parameter']))
    for i,row in enumerate(additions):
        row['id']=f"new-20260906-{len(previous['curves'])+i+1:02}"
    rows.extend(additions)
    return sorted(rows,key=lambda r:(-r['rank_lower_bound'],r['family'],r['parameter'])),previous

def build(output):
    if output.exists() or output.with_suffix('.csv').exists():
        raise FileExistsError('preserve curve index')
    rows,previous=expected_rows()
    names=set(previous['source_certificate_hashes'])|{ADDITIONAL.name}
    cert.write(output,{'schema':'elliptic-curves.new-high-rank-index.v15','sources':sources(),
        'previous_index_sha256':cert.hashed(OLD),'curves':rows,
        'source_certificate_hashes':{name:cert.hashed(ART/name) for name in sorted(names)},
        'catalogue':cert.read(ADDITIONAL)['catalogue'],
        'claim_boundary':'Exact independent-point lower bounds and absence from the pinned 593-equation catalogue. Existing local curve IDs are preserved. These are neither exact ranks nor a proof of universal novelty; the largest certified lower bound is reported by the actual point certificates.'})
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
        actual=checked_rank(model,points,[s['prime'] for s in old['signatures']],old['no_rational_2_torsion_prime'])
        if json.dumps(old,sort_keys=True)!=json.dumps(actual,sort_keys=True) or r['rank_lower_bound']!=len(points):
            raise ArithmeticError('exact rank proof differs')
        inv=cert.weierstrass_invariants(model);j=inv['c4']**3/inv['discriminant']
        if str(j)!=r['j_invariant'] or j in seen:
            raise ArithmeticError('j equality or duplicate needs separate treatment')
        if any(cert.isomorphic(model,q['ainvs']) for q in data['catalogue']['equations']):
            raise ArithmeticError('catalogue match in new inventory')
        seen.add(j)
        print('REPLAYED V15 INVENTORY',r['id'],'rank >=',len(points),flush=True)
    print('EXACT CURVE INVENTORY',len(seen),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ART/'new_high_rank_curve_index_v15.json');p.add_argument('--check',type=Path)
    a=p.parse_args();check(a.check) if a.check else build(a.output)
