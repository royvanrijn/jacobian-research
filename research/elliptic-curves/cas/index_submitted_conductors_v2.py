#!/usr/bin/env python3
"""Current four-curve index, retaining and checking the original three proofs."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import index_submitted_conductors as original

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'artifacts/generated-results/elliptic-curves/submitted628_conductor_v2'
CURRENT = ROOT/'elliptic-curves/data/submitted_conductors_current.json'


def write(path,value):
    path.write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')


def run(check=False):
    prior = original.expected()
    if prior != original.read(original.SUMMARY):
        raise ArithmeticError('historic three-exact/one-unknown index changed')
    rows = []
    for old in prior['rows']:
        identifier = old['icarm_id']
        path = OUT/'curve628.json' if identifier == 628 else ROOT/old['certificate']
        d = original.read(path)
        script = ('certify_submitted628_nfs.sage' if identifier == 628
                  else 'certify_submitted627_630_conductors.sage')
        argv = ['sage','-python',str(ROOT/'elliptic-curves/cas'/script),'--check']
        if identifier != 628:
            argv += ['--id',str(identifier)]
        subprocess.run(argv,check=True)
        if d['icarm_id'] != identifier or d['inventory_id'] != old['inventory_id']:
            raise ArithmeticError('curve identity differs')
        rows.append({'icarm_id':identifier,'inventory_id':d['inventory_id'],'status':'EXACT',
            'conductor':d['conductor'],'bad_primes':d['bad_primes'],
            'certificate':str(path.relative_to(ROOT)),'certificate_sha256':original.hashed(path)})
    summary = {'schema':'elliptic-curves.submitted-conductor-index.v2',
        'exact_count':4,'unresolved_count':0,'rows':rows,
        'previous_summary':str(original.SUMMARY.relative_to(ROOT)),
        'previous_summary_sha256':original.hashed(original.SUMMARY),
        'checker_sha256':original.hashed(Path(__file__)),
        'claim_boundary':'Four complete prime lists and exact conductors; no rank change or universal record claim.'}
    summary_path = OUT/'summary.json'
    digest = hashlib.sha256((json.dumps(summary,sort_keys=True,indent=2)+'\n').encode()).hexdigest()
    manifest = {'schema':'elliptic-curves.current-submitted-conductors.v2',
        'summary':str(summary_path.relative_to(ROOT)),'summary_sha256':digest}
    if check:
        if original.read(summary_path) != summary or original.read(CURRENT) != manifest:
            raise ArithmeticError('current four-curve index differs')
    else:
        if summary_path.exists() and original.read(summary_path) != summary:
            raise FileExistsError('preserve v2 summary')
        write(summary_path,summary)
        write(CURRENT,manifest)
    print('SUBMITTED CONDUCTORS: 4 EXACT; 0 UNKNOWN')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args = parser.parse_args()
    run(args.check)
