#!/usr/bin/env python3
"""Sage-free rank, global minimality, conductor and pinned-catalogue proof.

Lucas's criterion proves each prime: a has order p-1 modulo every prime
divisor of p, hence every such divisor is at least p. All factors of p-1
are recursively certified. No probable-prime test or factorization is used
in this checker. Local conductor criteria are stated in the canonical note.
"""
import argparse
from fractions import Fraction as F
import json
from math import gcd, prod
from pathlib import Path
import certify_compact_r17_candidates as cert

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
DATABASE = ROOT/'artifacts/local/elliptic-curves/next12-current-catalogue-v1/database.json'
CONDUCTOR = ART/'next12_rank22_exact_conductor_v1.json'
PRIMES = ART/'next12_conductor_lucas_primes_v1.json'
POINTS = ART/'prospective_mw16_next12_results_v1.json'

def sources():
    return {str(p.relative_to(ROOT)): cert.hashed(p) for p in (
        Path(__file__).resolve(), Path(cert.__file__).resolve(),
        ROOT/'elliptic-curves/cas/elliptic_candidate_record.py',
        ROOT/'elliptic-curves/cas/mod2_reduction_independence.py')}

def require(ok, message):
    if not ok:
        raise ArithmeticError(message)

def verify_primes(nodes, targets):
    verified = set()
    def prove(p):
        if p in verified:
            return
        require(p >= 2, 'invalid prime')
        node = nodes[str(p)]
        require(node['prime'] == str(p), 'prime node label differs')
        if p == 2:
            require(node == {'prime': '2', 'base_case': True}, 'base case differs')
        else:
            factors = [(int(q), int(e)) for q, e in node['p_minus_one_factorization']]
            require(len({q for q, e in factors}) == len(factors), 'repeated p-1 factor')
            require(all(2 <= q < p and e > 0 for q, e in factors), 'invalid recursive factors')
            require(prod(q**e for q, e in factors) == p-1, 'incomplete p-1 factorization')
            for q, e in factors:
                prove(q)
            a = int(node['witness'])
            require(1 < a < p and pow(a, p-1, p) == 1, 'Lucas full power failed')
            require(all(gcd(pow(a, (p-1)//q, p)-1, p) == 1 for q, e in factors), 'Lucas order failed')
        verified.add(p)
    for p in targets:
        prove(p)
    require(set(map(int, nodes)) == verified, 'extraneous unverified prime nodes')
    return len(verified)

def verify(data):
    require(data['sources'] == sources(), 'checker source hash differs')
    model = tuple(map(F, data['integral_model']))
    require(all(q.denominator == 1 for q in model), 'nonintegral displayed model')
    require(model[:3] == (1, 0, 0), 'transport outside fixed model shape')
    inv = cert.weierstrass_invariants(model)
    c4, c6, delta = (inv[k] for k in ('c4', 'c6', 'discriminant'))
    require(delta != 0 and delta.denominator == 1, 'invalid discriminant')
    require(str(delta) == data['discriminant'], 'discriminant mismatch')
    factors = [(int(p), int(e)) for p, e in data['discriminant_factorization']]
    require(len({p for p, e in factors}) == len(factors) and all(e > 0 for p, e in factors), 'invalid discriminant factors')
    prime_count = verify_primes(data['prime_proof_nodes'], [p for p, e in factors])
    require(prod(p**e for p, e in factors) == abs(delta), 'incomplete discriminant factorization')
    local = []
    for p, e in factors:
        if c4 % p:
            f, reason = 1, 'c4 unit: minimal multiplicative reduction'
        else:
            # v(Delta)<12 proves minimality of this integral equation.
            # For a minimal equation at p>=5, p|c4 and p|Delta mean additive
            # reduction and conductor exponent 2. Reject wild cases.
            require(p >= 5 and e < 12, 'local criterion unresolved')
            f, reason = 2, 'v(Delta)<12: minimal; p>=5 and p|c4: additive reduction'
        local.append({'prime': str(p), 'discriminant_valuation': e,
                      'conductor_valuation': f, 'reason': reason})
    conductor = prod(int(r['prime'])**r['conductor_valuation'] for r in local)
    require(local == data['local_proof'] and str(conductor) == data['conductor'], 'conductor proof differs')
    short = (F(0), F(0), F(0), -c4/48, -c6/864)
    require(short == tuple(map(F, data['short_model'])), 'short equation transport differs')
    points = [tuple(map(F, p)) for p in data['short_points']]
    integral_points = [(X-F(1, 12), Y-(X-F(1, 12))/2) for X, Y in points]
    require(integral_points == [tuple(map(F, p)) for p in data['integral_points']], 'point transport differs')
    require(all(cert.is_on_weierstrass_curve(model, p) for p in integral_points), 'integral model point check failed')
    proof = data['rank_certificate']
    actual = cert.checked_rank(short, points, [s['prime'] for s in proof['signatures']], proof['no_rational_2_torsion_prime'])
    require(json.dumps(actual, sort_keys=True) == json.dumps(proof, sort_keys=True), 'rank certificate differs')
    require(len(points) == data['rank_lower_bound'] == 22, 'fixed rank22 proof differs')
    catalogue = data['catalogue']
    raw = ROOT/catalogue['path']
    require(cert.hashed(raw) == catalogue['sha256'], 'pinned catalogue changed')
    projection = [{k: r[k] for k in ('id', 'ainvs', 'rank_lower_bound', 'conductor')} for r in cert.read(raw)['curves']]
    require(projection == catalogue['projection'] and len(projection) == 586, 'catalogue projection differs')
    require(len({r['id'] for r in projection}) == 586, 'duplicate catalogue IDs')
    matches = [r['id'] for r in projection if cert.isomorphic(model, tuple(map(F, r['ainvs'])))]
    smaller = sorted([r for r in projection if r['rank_lower_bound'] >= 22 and r['conductor'] and int(r['conductor']) < conductor], key=lambda r: int(r['conductor']))
    missing = sorted(r['id'] for r in projection if r['rank_lower_bound'] >= 22 and not r['conductor'])
    minimum = min((r for r in projection if r['rank_lower_bound'] >= 23 and r['conductor']), key=lambda r: int(r['conductor']))
    comparison = {'q_isomorphism_matches': matches, 'smaller_recorded_rank22_ids': [r['id'] for r in smaller],
                  'rank22_missing_conductor_ids': missing, 'minimum_recorded_rank23_id': minimum['id'],
                  'minimum_recorded_rank23_conductor': minimum['conductor']}
    require(comparison == data['comparison'] and not matches, 'catalogue comparison differs')
    require([r['id'] for r in smaller] == [376, 575] and minimum['id'] == 539 and conductor < int(minimum['conductor']), 'fixed conductor placement differs')
    print('SAGE-FREE PROOF PASS: rank >=22; global minimal; exact conductor', conductor,
          '; Lucas prime nodes', prime_count, '; third among recorded rank>=22 conductors in pinned586 catalogue', flush=True)

def build(output):
    if output.exists():
        raise FileExistsError('preserve small-conductor proof')
    conductor = cert.read(CONDUCTOR)
    primes = cert.read(PRIMES)
    row = next(r for r in cert.read(POINTS)['curves'] if r['parameter'] == '3/17')
    require(primes['status'] == 'COMPLETE_LUCAS_CERTIFICATES', 'prime tree unfinished')
    require(conductor['status'] == 'PASS_EXACT_CONDUCTOR', 'Sage conductor unfinished')
    model = tuple(map(F, conductor['integral_model']))
    c4 = cert.weierstrass_invariants(model)['c4']
    local = []
    for p, e in conductor['discriminant_factorization']:
        f, reason = (1, 'c4 unit: minimal multiplicative reduction') if c4 % int(p) else (
            2, 'v(Delta)<12: minimal; p>=5 and p|c4: additive reduction')
        local.append({'prime': p, 'discriminant_valuation': e, 'conductor_valuation': f, 'reason': reason})
    projection = [{k: r[k] for k in ('id', 'ainvs', 'rank_lower_bound', 'conductor')} for r in cert.read(DATABASE)['curves']]
    minimum = min((r for r in projection if r['rank_lower_bound'] >= 23 and r['conductor']), key=lambda r: int(r['conductor']))
    data = {'schema': 'elliptic-curves.small-conductor-rank22-proof.v1', 'sources': sources(),
            'provenance': {str(p.relative_to(ROOT)): cert.hashed(p) for p in (CONDUCTOR, PRIMES, POINTS)},
            'curve_id': 'new-20260905-36', 'family': row['family'], 'parameter': '3/17',
            'integral_model': conductor['integral_model'], 'short_model': row['curve'],
            'discriminant': conductor['discriminant'], 'discriminant_factorization': conductor['discriminant_factorization'],
            'prime_proof_nodes': primes['nodes'], 'local_proof': local, 'conductor': conductor['conductor'],
            'short_points': row['points'], 'rank_certificate': row['rank_certificate'], 'rank_lower_bound': 22,
            'integral_points': [[str(F(X)-F(1, 12)), str(F(Y)-(F(X)-F(1, 12))/2)] for X, Y in row['points']],
            'catalogue': {'path': str(DATABASE.relative_to(ROOT)), 'sha256': cert.hashed(DATABASE),
                          'url': 'https://elliptic-rank.icarm.cloud/database.json', 'projection': projection},
            'comparison': {'q_isomorphism_matches': [], 'smaller_recorded_rank22_ids': [376, 575],
                           'rank22_missing_conductor_ids': [537, 543, 545, 581],
                           'minimum_recorded_rank23_id': minimum['id'], 'minimum_recorded_rank23_conductor': minimum['conductor']},
            'claim_boundary': 'Unconditional rank lower bound22 and exact conductor of this global minimal equation. Third among rank>=22 entries with recorded conductors in the pinned586 catalogue; four rank>=22 entries have no recorded conductor. No exact rank, universal novelty or absolute world placement is proved.'}
    verify(data)
    cert.write(output, data)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--output', type=Path)
    g.add_argument('--check', type=Path)
    a = p.parse_args()
    verify(cert.read(a.check)) if a.check else build(a.output)
