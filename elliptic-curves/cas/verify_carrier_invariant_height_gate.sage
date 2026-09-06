#!/usr/bin/env sage-python
"""Independent Sage invariant and integer inequality replay of the carrier gate."""
import argparse,json,hashlib
from pathlib import Path
from sage.all import QQ,ZZ,EllipticCurve,gcd
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'soluble_pair_carrier_height_v1.json';PROOF=ART/'soluble_pair_carrier_invariant_height_v1.json';OUT=ART/'soluble_pair_carrier_invariant_height_replay_v1.json'
def hashed(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def expected():
    source=json.loads(INPUT.read_text());proof=json.loads(PROOF.read_text())
    if proof['status']!='PASS' or len(source['rows'])!=len(proof['rows']):raise ArithmeticError('complete fixed proof required')
    for name,h in proof['sources'].items():
        if hashed(ROOT/name)!=h:raise ArithmeticError('proof input binding differs')
    checked=0;excluded=[]
    for original,row in zip(source['rows'],proof['rows']):
        if original['word']!=row['word']:raise ArithmeticError('word differs')
        if 'curve' not in original:
            if row['status']!='NO_SPECIALIZED_MODEL_IN_INPUT':raise ArithmeticError('censoring differs')
            continue
        E=EllipticCurve(QQ,original['curve']);c4,c6=map(ZZ,E.c_invariants());G=gcd(c4,c6)
        if str(c4)!=row['c4'] or str(c6)!=row['c6'] or str(G)!=row['invariant_gcd']:raise ArithmeticError('independent invariants differ')
        B=QQ(c6*c6)/G**3;L=ZZ(row['normalized_coefficient_size_lower_bound'])
        if not (1224*(L-1)+521)**2<B<=(1224*L+521)**2 or L.nbits()!=row['normalized_coefficient_bits_lower_bound']:
            raise ArithmeticError('exact integer lower bound differs')
        exceeds=bool(B>(1224*(2**360-1)+521)**2)
        if exceeds!=row['all_normalized_360_bit_integral_models_excluded']:raise ArithmeticError('height comparison differs')
        if exceeds:excluded.append(original['compact_parameter'])
        checked+=1
    if checked!=11 or len(excluded)!=10 or len(set(excluded))!=7:raise ArithmeticError('fixed counts differ')
    return {'schema':'elliptic-curves.soluble-pair-invariant-height-replay.v1','status':'PASS',
        'sources':{str(p.relative_to(ROOT)):hashed(p) for p in (Path(__file__).resolve(),INPUT,PROOF)},
        'exact_invariant_rows':checked,'excluded_rows':len(excluded),'excluded_distinct_parameters':len(set(excluded)),
        'claim_boundary':'Independent Sage invariants and exact rational inequalities replay the factor-free normalized integral-model height bound. No global minimal-model algorithm, factorization or point search is used.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if json.loads(OUT.read_text())!=d:raise ArithmeticError('replay record differs')
    else:
        if OUT.exists():raise FileExistsError('preserve independent height replay')
        OUT.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print('CARRIER INVARIANT HEIGHT REPLAY PASS',d['exact_invariant_rows'],d['excluded_distinct_parameters'],flush=True)
