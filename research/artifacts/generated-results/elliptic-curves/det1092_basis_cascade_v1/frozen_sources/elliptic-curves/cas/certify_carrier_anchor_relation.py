#!/usr/bin/env python3
"""Sage-free exact eighteen-dimensional span proof for the fixed carrier anchor."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank
from half_lattice_pointed_sieve import linear_combination_python
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/carrier-anchor-relation-v1';INPUT=ART/'soluble_pair_carrier_height_v1.json';OUT=ART/'soluble_pair_carrier_anchor_relation_v1.json'
def expected():
    p=cert.read(D/'protocol.json');q=cert.read(D/'proposal.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()) or q['protocol_sha256']!=cert.hashed(D/'protocol.json'):raise ArithmeticError('frozen relation input differs')
    r=next(r for r in cert.read(INPUT)['rows'] if r['word']==[1,1]);model=tuple(map(cert.F,r['curve']))
    if q['status']!='EXACT_GROUP_RELATION_PENDING_INDEPENDENT_REPLAY' or q['curve']!=r['curve'] or q['basis']!=r['generic_points']+[r['supplied_points'][0]] or q['target']!=r['supplied_points'][1] or q['rank_certificate']!=r['rank_certificate']:raise ArithmeticError('exact declared relation required')
    d=q['denominator'];word=q['coefficients']
    if type(d)!=int or not 1<=d<=p['maximum_denominator'] or len(word)!=18 or any(type(c)!=int or abs(c)>p['maximum_coefficient'] for c in word):raise ArithmeticError('bounded relation differs')
    basis=[tuple(map(cert.F,P)) for P in q['basis']];target=tuple(map(cert.F,q['target']));proof=q['rank_certificate'];primes=[s['prime'] for s in proof['signatures']]
    actual=checked_rank(model,basis,primes,proof['no_rational_2_torsion_prime'])
    if json.loads(json.dumps(actual))!=proof:raise ArithmeticError('independent18-point proof differs')
    generic=checked_rank(model,basis[:17],primes,proof['no_rational_2_torsion_prime'])
    if linear_combination_python(model,[*basis,target],[*word,-d]) is not None:raise ArithmeticError('rational group identity differs')
    paths=[Path(__file__).resolve(),CAS/'memory_rank_certificate.py',CAS/'half_lattice_pointed_sieve.py',INPUT,D/'protocol.json',D/'proposal.json']
    return {'schema':'elliptic-curves.soluble-pair-carrier-anchor-relation.v1','status':'PASS',
        'sources':{str(a.relative_to(ROOT)):cert.hashed(a) for a in paths},'parameter':r['compact_parameter'],
        'curve':r['curve'],'generic_points':r['generic_points'],'supplied_points':r['supplied_points'],
        'independent_points':q['basis'],'independence_certificate':proof,'generic_independence_certificate':generic,
        'denominator':d,'coefficients':word,'generic_span_rank':17,'combined_span_rank':18,'supplied_quotient_rank':1,
        'claim_boundary':'The19 supplied points at the already known anchor have rational span exactly18:17 independent generic points, one extra independent direction, and an exact integer relation for the other supplied point. Thus this pair contributes exactly one direction modulo the generic span on this fibre. The full curve rank, its remaining quotient and other carrier fibres are unchanged.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('anchor relation replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve anchor span proof')
        checkpoint(OUT,d)
    print('EXACT CARRIER ANCHOR SPAN',d['combined_span_rank'],'QUOTIENT',d['supplied_quotient_rank'],flush=True)
