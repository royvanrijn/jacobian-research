#!/usr/bin/env python3
"""Prepare cached arithmetic for the exceptional panel before optional descent.

45 seconds and 1 GiB observed RSS per curve, one worker at a time. Each
worker gets its equation, known discriminant support and resource limit.
Inputs, logs and partial stages survive under --workdir; completed arithmetic
also survives across workdirs in the shared content-addressed cache. Full
BNF/Selmer is requested only by --full-selmer. No point search is performed.
"""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_elkies_2026_relative_2selmer_open import DEFAULT_SAGE_PYTHON, file_sha256
from research_runtime.supervisor import supervise_source
from research_runtime.store import checkpoint

WORKER = r'''
import json
from pathlib import Path
from sage.all import pari
from sage.version import version
from research_runtime.sage_arithmetic import SageArithmetic
from research_runtime.store import FactStore, checkpoint
p = json.loads(Path(INPUT_PATH).read_text())
pari.allocatemem(p['stack_bytes'])
pari.setrand(1)
out = {'software': {'sage': version, 'pari': str(pari.version())},
       'full_selmer_dimension': None, 'covers': None, 'stages': []}
def stage(s):
    out['stages'].append(s)
    checkpoint(Path(OUTPUT_PATH), out)
    print(s, flush=True)
arithmetic = SageArithmetic(FactStore(Path(p['cache_dir'])))
stage('arithmetic_context_start')
context = arithmetic.prepare(p['model'], factor_primes=p['factor_hint_primes'], discover=True)
out['arithmetic_context'] = context.record()
out['arithmetic_context_key'] = context.key
stage('minimal_model_and_factorization_cached')
field = arithmetic.field(context.two_torsion, factor_primes=[2,*context.bad_primes], discover=True)
out['two_torsion_context_key'] = context.two_torsion.key
stage('maximal_order_and_polredabs_cached')
out['scheduling_features'] = arithmetic.scheduling_features(context.two_torsion)
out['full_selmer_requested'] = p['full_selmer']
if p['full_selmer']:
    stage('cached_field_complete_selmer_start')
    descent = arithmetic.full_selmer(context, requirement='complete-selmer', discover=True)
    out['descent'] = descent
    out['full_selmer_dimension'] = descent['full_selmer_dimension']
    stage('complete_selmer_basis_cached')
else:
    stage('arithmetic_prepared_no_bnf_requested')
'''


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workdir', type=Path, required=True)
    ap.add_argument('--cache-dir', type=Path, default=ROOT/'artifacts/local/arithmetic-cache')
    ap.add_argument('--full-selmer', action='store_true', help='request an unconditional complete Selmer basis after cached setup')
    ap.add_argument('--factor-support', type=Path, help='JSON mapping curve IDs to known discriminant prime supports')
    args = ap.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    if (args.workdir / 'run.json').exists():
        ap.error('preserve the completed run; use a fresh workdir')
    # Only equations are extracted; no known point is sent to the worker.
    source = ROOT / 'artifacts/generated-results/elliptic-curves/exceptional_soluble_selmer_panel_v1.json'
    panel = json.loads(source.read_text())
    from verify_icarm_curve398_rank30 import DISCRIMINANT_FACTORIZATION
    supports = {'398': [p for p, _ in DISCRIMINANT_FACTORIZATION]}
    if args.factor_support:
        supports.update(json.loads(args.factor_support.read_text()))
    records = []
    for curve_id in (398, 400, 543):
        row = next(r for r in panel['curves'] if r['curve_id'] == curve_id)
        payload = {'model': row['model'], 'stack_bytes': 256_000_000,
                   'factor_hint_primes': supports.get(str(curve_id), []),
                   'cache_dir': str(args.cache_dir.resolve()), 'full_selmer': args.full_selmer}
        (args.workdir / f'input-{curve_id}.json').write_text(json.dumps(payload, indent=2)+'\n')
        result = args.workdir / f'worker-{curve_id}.json'
        measurement = supervise_source(str(DEFAULT_SAGE_PYTHON), WORKER, payload,
            result, args.workdir / f'worker-{curve_id}.log', timeout=45,
            rss_limit_bytes=1_073_741_824)
        records.append({'curve_id': curve_id, 'input': payload,
                        'measurement': measurement,
                        'worker': json.loads(result.read_text()) if result.exists() else None})
        out = {'schema': 'elliptic-curves.exceptional-selmer-feasibility.v2',
               'source_sha256': file_sha256(source), 'runner_sha256': file_sha256(Path(__file__)),
               'worker_source': WORKER, 'records': records,
               'point_searches': 0, 'maximum_total_worker_seconds': 135,
               'claim': 'Incomplete probes yield no complement, Selmer dimension or rank bound.'}
        checkpoint(args.workdir / 'run.json', out)
        print(curve_id, measurement['outcome'], records[-1]['worker'], flush=True)


if __name__ == '__main__':
    main()
