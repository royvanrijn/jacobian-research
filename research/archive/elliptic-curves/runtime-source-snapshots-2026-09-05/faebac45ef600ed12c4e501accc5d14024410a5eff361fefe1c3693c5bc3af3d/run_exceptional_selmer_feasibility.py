#!/usr/bin/env python3
"""Three equation-only descent probes; no point search or pairing attestation.

45 seconds and 1 GiB observed RSS per curve, one worker at a time. Each
worker gets only its equation and resource limit. Inputs, logs and partial
stages survive under --workdir. A completed ell2cover still needs independent
class/map binding before its quartics enter the mathematical panel.
"""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_elkies_2026_relative_2selmer_open import (
    DEFAULT_SAGE_PYTHON, file_sha256, supervise_source,
)

WORKER = r'''
import json
from pathlib import Path
from sage.all import EllipticCurve, QQ, pari
from sage.version import version
p = json.loads(Path(INPUT_PATH).read_text())
pari.allocatemem(p['stack_bytes'])
pari.setrand(1)
out = {'software': {'sage': version, 'pari': str(pari.version())},
       'full_selmer_dimension': None, 'covers': None, 'stages': []}
def stage(s):
    out['stages'].append(s)
    Path(OUTPUT_PATH).write_text(json.dumps(out, indent=2)+'\n')
    print(s, flush=True)
E = EllipticCurve(list(map(QQ, p['model'])))
stage('ellrankinit_start')
context = pari(E).ellrankinit()
stage('ellrankinit_complete')
stage('bnfcertify_start')
certificates = [int(pari.bnfcertify(b)) for b in context[2]]
assert certificates and all(c == 1 for c in certificates)
out['bnfcertify'] = certificates
stage('bnfcertify_complete')
stage('ell2cover_start')
cs = context.ell2cover()
out['covers'] = [{'quartic': str(c[0]), 'map': [str(v) for v in c[1]]} for c in cs]
out['full_selmer_dimension'] = len(cs)
stage('ell2cover_complete')
'''


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workdir', type=Path, required=True)
    args = ap.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    if (args.workdir / 'run.json').exists():
        ap.error('preserve the completed run; use a fresh workdir')
    # Only equations are extracted; no known point is sent to the worker.
    source = ROOT / 'artifacts/generated-results/elliptic-curves/exceptional_soluble_selmer_panel_v1.json'
    panel = json.loads(source.read_text())
    records = []
    for curve_id in (398, 400, 543):
        row = next(r for r in panel['curves'] if r['curve_id'] == curve_id)
        payload = {'model': row['model'], 'stack_bytes': 256_000_000}
        (args.workdir / f'input-{curve_id}.json').write_text(json.dumps(payload, indent=2)+'\n')
        result = args.workdir / f'worker-{curve_id}.json'
        measurement = supervise_source(str(DEFAULT_SAGE_PYTHON), WORKER, payload,
            result, args.workdir / f'worker-{curve_id}.log', timeout=45,
            rss_limit_bytes=1_073_741_824)
        records.append({'curve_id': curve_id, 'input': payload,
                        'measurement': measurement,
                        'worker': json.loads(result.read_text()) if result.exists() else None})
        out = {'schema': 'elliptic-curves.exceptional-selmer-feasibility.v1',
               'source_sha256': file_sha256(source), 'runner_sha256': file_sha256(Path(__file__)),
               'worker_source': WORKER, 'records': records,
               'point_searches': 0, 'maximum_total_worker_seconds': 135,
               'claim': 'Incomplete probes yield no complement, Selmer dimension or rank bound.'}
        (args.workdir / 'run.json').write_text(json.dumps(out, indent=2, sort_keys=True)+'\n')
        print(curve_id, measurement['outcome'], records[-1]['worker'], flush=True)


if __name__ == '__main__':
    main()
