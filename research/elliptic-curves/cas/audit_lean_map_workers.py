#!/usr/bin/env python3
"""Fixed retained-map equivalence and supervised cost audit; no point search."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import shutil

from research_runtime.supervisor import run, Limits
from research_runtime.store import checkpoint
from v3_warm_engine import certified_state
from pointed_quartic_search import PointedQuarticSearch
import pari_pointed_backend as backend

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(parent, output):
    output.mkdir(exist_ok=False)
    seed = json.loads((parent/'seed.json').read_text())
    state = certified_state(tuple(map(F, seed['curve'])),
        tuple(tuple(map(F, p)) for p in seed['points']), seed['proof'])
    worker = CAS/'lean_preconditioned_map_worker.py'
    inputs = {str((parent/'seed.json').relative_to(ROOT)): sha(parent/'seed.json')}
    rows = []
    # Fixed first/last four centres of the completed first 50-centre pass.
    for ci in (0, 1, 2, 3, 46, 47, 48, 49):
        for policy in ('preconditioned_full', 'factor_free'):
            old = parent/'epoch-00'/f'map-{ci:04d}-{policy}'
            wd = output/old.name
            wd.mkdir()
            for name in ('input.json', 'result.json', 'supervisor.json'):
                inputs[str((old/name).relative_to(ROOT))] = sha(old/name)
            shutil.copyfile(old/'input.json', wd/'input.json')
            supervisor = run([shutil.which('sage'), '-python', str(worker),
                '--input', str(wd/'input.json'), '--output', str(wd/'result.json')],
                limits=Limits(wall_seconds=5, rss_bytes=1024**3),
                log_path=wd/'worker.log', result_path=wd/'result.json')
            checkpoint(wd/'supervisor.json', supervisor)
            if supervisor['outcome'] != 'completed':
                raise ArithmeticError('lean worker did not complete')
            original = json.loads((old/'result.json').read_text())['mapping']
            actual = json.loads((wd/'result.json').read_text())['mapping']
            if actual != original:
                raise ArithmeticError('retained exact map changed')
            search = PointedQuarticSearch(state=state,
                centre={'coefficients': actual['centre']['representative']},
                coordinate_policy=actual['coordinate_policy'])
            backend.validate_map(search, actual)
            prior = json.loads((old/'supervisor.json').read_text())
            rows.append(dict(centre_index=ci, policy=policy,
                old_wall_seconds=prior['wall_seconds'],
                lean_wall_seconds=supervisor['wall_seconds'],
                exact_map_equal=True, identity_verified=True,
                input_sha256=sha(wd/'input.json'), result_sha256=sha(wd/'result.json'),
                supervisor_sha256=sha(wd/'supervisor.json')))
            print('EXACT_LEAN_MAP', ci, policy, supervisor['wall_seconds'], flush=True)
    sources = {str(p.relative_to(ROOT)):sha(p) for p in
        [Path(__file__), worker, CAS/'lean_factor_free_pari_mapping.sage',
         CAS/'lean_preconditioned_full_pari_mapping.sage']}
    if any(sha(ROOT/p) != h for p,h in inputs.items()):
        raise ArithmeticError('retained inputs changed')
    checkpoint(output/'report.json', dict(status='PASS_FIXED_EXACT_MAP_EQUIVALENCE',
        point_searches=0, rows=rows, inputs=inputs, sources=sources,
        old_map_wall_seconds=sum(r['old_wall_seconds'] for r in rows),
        lean_map_wall_seconds=sum(r['lean_wall_seconds'] for r in rows),
        claim_boundary='Fixed retained maps and exact identity checks only. Historical versus current worker timings; no whole-search speedup or universal equivalence claim.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.parent.resolve(), args.output.resolve())
