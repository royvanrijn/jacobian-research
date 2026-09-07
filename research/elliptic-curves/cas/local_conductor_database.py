#!/usr/bin/env python3
"""Read certified conductor information alongside the unchanged research inventory."""
import hashlib
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
from refresh_icarm_local_database import load_catalogue, load_inventory

ROOT = Path(__file__).resolve().parents[2]
CURRENT = ROOT/'elliptic-curves/data/conductor_screen_current.json'


def checked(path, digest):
    path = ROOT/path
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise ArithmeticError('conductor database hash mismatch: '+str(path))
    return json.loads(path.read_text())


def load_conductor_inventory():
    """Return all201 rows with conductor status, proof paths and unresolved bounds."""
    manifest = json.loads(CURRENT.read_text())
    screen = checked(manifest['summary'],manifest['summary_sha256'])
    submitted = checked(manifest['submitted_summary'],manifest['submitted_summary_sha256'])
    small = checked(manifest['small_rank22_certificate'],manifest['small_rank22_sha256'])
    proofs = {}
    for identifier, source in screen['certificates'].items():
        d = checked(source['path'],source['sha256'])
        if d['id'] != identifier:
            raise ArithmeticError('conductor certificate ID differs')
        proofs[identifier] = (d['integral_curve'],{
            'status':d['status'],'conductor':d['exact_conductor'],'bad_primes':d['bad_primes'],
            'conductor_divisor':d['conductor_divisor'],'conductor_upper_bound':d['conductor_upper_bound'],
            'certificate':source['path'],'recorded_benchmark':d['recorded_benchmark']})
    for source in submitted['rows']:
        d = checked(source['certificate'],source['certificate_sha256'])
        proofs[d['inventory_id']] = (d['ainvs'],{
            'status':'EXACT','conductor':d['conductor'],'bad_primes':d['bad_primes'],
            'certificate':source['certificate']})
    proofs[small['curve_id']] = (small['short_model'],{
        'status':'EXACT','conductor':small['conductor'],
        'bad_primes':[r['prime'] for r in small['local_proof'] if r['conductor_valuation']],
        'certificate':manifest['small_rank22_certificate']})
    public = {r['id']:r for r in load_catalogue()['curves']}
    rows = load_inventory()
    for r in rows:
        if r['id'] in proofs:
            model, information = proofs[r['id']]
            if not cert.isomorphic(r['curve'],model):
                raise ArithmeticError('conductor proof equation differs from research curve')
            r['conductor_information'] = information
        else:
            reports = [{'icarm_id':i,'conductor':public[i]['conductor']}
                       for i in r['current_catalogue_matches'] if public[i]['conductor']]
            r['conductor_information'] = {'status':'REPORTED' if reports else 'UNKNOWN',
                'reported_public_values':reports,'conductor':None,'bad_primes':None}
    if set(proofs)-{r['id'] for r in rows}:
        raise ArithmeticError('orphan conductor proof')
    return rows


def index_current_results(summary=None):
    paths = {
        'summary':'artifacts/generated-results/elliptic-curves/new_curve_conductor_screen_v1/summary.json',
        'submitted_summary':'artifacts/generated-results/elliptic-curves/submitted628_conductor_v2/summary.json',
        'small_rank22_certificate':'artifacts/generated-results/elliptic-curves/small_conductor_rank22_proof_v1.json'}
    if summary:
        paths['summary'] = summary
    elif CURRENT.exists():
        paths['summary'] = json.loads(CURRENT.read_text())['summary']
    manifest = {'schema':'elliptic-curves.current-conductor-screen.v1',**paths}
    for key, hashkey in [('summary','summary_sha256'),('submitted_summary','submitted_summary_sha256'),
                         ('small_rank22_certificate','small_rank22_sha256')]:
        manifest[hashkey] = hashlib.sha256((ROOT/paths[key]).read_bytes()).hexdigest()
    CURRENT.write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')


if __name__ == '__main__':
    import argparse
    from collections import Counter
    parser = argparse.ArgumentParser()
    parser.add_argument('--index',action='store_true')
    parser.add_argument('--summary',help='Repository-relative path of the certified summary to select')
    args = parser.parse_args()
    if args.index:
        index_current_results(args.summary)
    rows = load_conductor_inventory()
    print('CONDUCTOR DATABASE',len(rows),dict(Counter(r['conductor_information']['status'] for r in rows)))
