#!/usr/bin/env python3
"""Bind complete endpoint transport, equation comparisons and odd-prime section proofs."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'compact_endpoint_summary_v1.json'
def expected():
    paths={Path(__file__).resolve()}
    def read(p):paths.add(p);return cert.read(p)
    endpoint=read(ART/'compact_atlas_endpoints_v2.json');odd=read(LOCAL/'compact-endpoint-odd-primes-v1/ledger.json');old=read(LOCAL/'compact-atlas-endpoints-v1/run.supervisor.json');terminal=read(LOCAL/'compact-atlas-endpoints-v2/check.supervisor.json')
    if old['returncode']==0 or terminal['outcome']!='completed' or terminal['returncode']!=0 or endpoint['status']!='COMPLETE_DECLARED_ENDPOINT_AUDIT' or odd['status']!='PASS':raise ArithmeticError('preserved failure and successful endpoint replays required')
    if len(endpoint['rows'])!=22 or endpoint['within_roster_isomorphic_pairs']:raise ArithmeticError('fixed22 distinct endpoint roster differs')
    rows=[];keyed={(r['family'],r['endpoint']):r for r in odd['rows']};singular=[]
    for r in endpoint['rows']:
        if r['status']=='SINGULAR_FIBRE':singular.append([r['family'],r['endpoint']]);continue
        q=keyed[(r['family'],r['endpoint'])];ip=ROOT/q['input'];op=ROOT/q['output'];data=read(ip);proof=read(op)
        if r['status']!='CERTIFIED_SPECIALIZED_SUBGROUP' or q['status']!='PASS' or cert.hashed(ip)!=q['input_sha256'] or cert.hashed(op)!=q['output_sha256'] or proof['input_sha256']!=cert.hashed(ip):raise ArithmeticError('endpoint point proof binding differs')
        if data['points']!=r['generic_points'] or proof['points']!=r['generic_points'] or data['curve']!=r['curve'] or r['catalogue_matches'] or r['previous_matches']:raise ArithmeticError('generic transport or equation exclusion differs')
        if q['mod2_lower_bound']!=r['rank_lower_bound'] or q['odd_modulus_lower_bounds']!={str(a['modulus']):a['finite_column_rank'] for a in proof['audits']} or any(v!=r['rank_lower_bound'] for v in q['odd_modulus_lower_bounds'].values()):raise ArithmeticError('full generic-section finite bounds differ')
        rows.append({'family':r['family'],'endpoint':r['endpoint'],'finite_generic_sections':len(r['generic_points']),'rank_lower_bound':r['rank_lower_bound'],'odd_modulus_lower_bounds':q['odd_modulus_lower_bounds']})
    if len(rows)!=21 or singular!=[['a1-fibration-03','infinity']] or endpoint['previous_address_equations']!=528:raise ArithmeticError('fixed endpoint outcomes differ')
    return {'schema':'elliptic-curves.compact-endpoint-summary.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in sorted(paths)},'endpoint_count':22,'nonsingular_distinct_curves':21,'singular_endpoints':singular,'rows':rows,'catalogue_snapshot_curves':593,'previous_address_equations':528,'catalogue_or_previous_matches':0,'rational_point_searches_run':0,'claim_boundary':'Exact omitted endpoint equations and transported section values, with independent finite point lower bounds11..17 and matching modulo3/5 bounds. Lower finite ranks do not prove rational dependence or exact specialization rank. One singular displayed endpoint is not an elliptic curve in this fibre chart. The21 unmatched nonsingular curves have not received an endpoint point-search campaign. No high-rank inventory addition, rank upper bound, saturation, conductor or universal novelty.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=d:raise ArithmeticError('endpoint summary differs')
    else:
        if OUT.exists():raise FileExistsError('preserve endpoint summary')
        checkpoint(OUT,d)
    print('ENDPOINT SUMMARY PASS:21 nonsingular unmatched;one singular;no point search',flush=True)
