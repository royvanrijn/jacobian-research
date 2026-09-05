#!/usr/bin/env sage-python
"""Checkpointed, bounded PARI 2-descent probe for the new small-conductor curve.

No upper bound is accepted from incomplete BNF data. Even completed CAS
output remains unpromoted until a separate replay and source review.
"""
import argparse
from pathlib import Path
import sys
from sage.all import QQ, pari
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
INPUT = ROOT/'artifacts/generated-results/elliptic-curves/small_conductor_rank22_proof_v1.json'

def run(directory):
    output = directory/'result.json'
    if output.exists():
        raise FileExistsError('preserve descent probe')
    row = cert.read(INPUT)
    pari.allocatemem(64000000, 512000000)
    pari.set_real_precision_bits(256)
    E = pari.ellinit([QQ(q) for q in row['integral_model']])
    points = [[QQ(q) for q in P] for P in row['integral_points']]
    if len(points) != 22 or any(not E.ellisoncurve(P) for P in points):
        raise ArithmeticError('fixed22 input point gate differs')
    data = {'schema': 'elliptic-curves.small-conductor-descent-probe.v1',
            'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in (INPUT, Path(__file__).resolve())},
            'pari_version': str(pari.version()), 'known_rank_lower_bound': 22,
            'status': 'STARTING_RANK_INITIALIZATION', 'unconditional_upper_bound': None,
            'claim_boundary': 'Incomplete stages prove no upper bound. Completed CAS output requires separate replay before mathematical promotion.'}
    checkpoint(output, data)
    print('START RANK INITIALIZATION', flush=True)
    rank_data = E.ellrankinit(precision=256)
    raw = directory/'rankinit.pari';raw.write_text(str(rank_data)+'\n')
    data.update(status='STARTING_BNF_CERTIFICATION', rankinit_sha256=cert.hashed(raw))
    checkpoint(output, data)
    print('RANK INITIALIZATION COMPLETE', flush=True)
    # PARI2.17 rankinit stores the nonsplit 2-division BNF in its third
    # component. Reject a changed structure instead of guessing a location.
    if len(rank_data) != 3 or len(rank_data[2]) != 1:
        raise ArithmeticError('fixed PARI2.17 irreducible cubic structure differs')
    bnf = rank_data[2][0]
    bnf_raw = directory/'bnf.pari';bnf_raw.write_text(str(bnf)+'\n')
    certified = int(bnf.bnfcertify())
    if certified != 1:
        raise ArithmeticError('BNF unconditional certification failed')
    data.update(status='STARTING_TWO_DESCENT', bnfcertify=1, bnf_sha256=cert.hashed(bnf_raw))
    checkpoint(output, data)
    print('BNF CERTIFIED; START TWO DESCENT', flush=True)
    result = rank_data.ellrank(0, points, precision=256)
    raw = directory/'ellrank.pari';raw.write_text(str(result)+'\n')
    data.update(status='COMPLETED_CAS_DATA_UNREPLAYED', ellrank_sha256=cert.hashed(raw),
                reported_lower=int(result[0]), reported_upper=int(result[1]), reported_sha_quotient_rank=int(result[2]),
                reported_points=[[str(q) for q in P] for P in result[3]])
    checkpoint(output, data)
    print('TWO DESCENT COMPLETED; NOT YET PROMOTED', data['reported_lower'], data['reported_upper'], flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, required=True)
    run(p.parse_args().directory.resolve())
