#!/usr/bin/env python3
"""Read final checkpoint names only while the frozen census is running.

The original status glob can race workers' transient *.tmp.json files.
Keep this observer separate so the active arithmetic source hashes stay fixed.
"""
import argparse
from collections import Counter
import json
from pathlib import Path


def snapshot(out):
    rows=json.loads((out/'population.json').read_text())['curves']
    counts={}
    for mode in ['base','local','class']:
        c=Counter()
        for row in rows:
            p=out/mode/(row['curve_key']+'.json')
            if p.is_file():c[json.loads(p.read_text())['status']]+=1
        counts[mode]=dict(sorted(c.items()))
    return {'population_count':len(rows),'stage_status_counts':counts}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[2]/'artifacts/local/elliptic-curves/wide-arithmetic-profile-v1')
    args=ap.parse_args();print(json.dumps(snapshot(args.output),sort_keys=True,indent=2))


if __name__=='__main__':main()
