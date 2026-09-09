#!/usr/bin/env python3
"""Fixed prospective-coordinate exposure against completed parent boxes; no point search."""
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
from pointed_box_equivalence import box_key
import run_parent_seed_v3 as old_runner

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
    protocol = json.loads((parent/'protocol.json').read_text())
    terminal = json.loads((parent/'terminal.json').read_text())
    verified = json.loads((parent/'verified.json').read_text())
    old_runner.guard(protocol)
    if verified['status'] != 'PASS_INDEPENDENT_PRODUCTIVE_V3_REPLAY' or verified['terminal_sha256'] != sha(parent/'terminal.json'):
        raise ArithmeticError('verified parent required')
    selection = json.loads((parent/'epoch-00/landscape/selection.json').read_text())
    charts = [json.loads(x.read_text()) for x in sorted((parent/'epoch-00').glob('chart-*.json'))]
    completed = {}
    for chart in charts:
        if chart['search']['status'] == 'bounded_search_complete' and chart['search']['height_bound'] == 125000 and chart['search']['infinity_checked']:
            key = tuple(chart['mapping']['centre']['point'])
            completed.setdefault(key,set()).add(box_key(chart['mapping']['matrix']))
    inputs = {str((parent/'seed.json').relative_to(ROOT)): sha(parent/'seed.json')}
    rows = []
    for ci in range(16):
        for policy in ('preconditioned_full',):
            wd = output/f'map-{ci:04d}-{policy}'
            wd.mkdir()
            centre = selection['centres'][ci]
            checkpoint(wd/'input.json',dict(curve=seed['curve'],points=seed['points'],centre=centre,policy=policy))
            supervisor = run([shutil.which('sage'), '-python', str(worker),
                '--input', str(wd/'input.json'), '--output', str(wd/'result.json')],
                limits=Limits(wall_seconds=5, rss_bytes=1024**3),
                log_path=wd/'worker.log', result_path=wd/'result.json')
            checkpoint(wd/'supervisor.json', supervisor)
            if supervisor['outcome'] != 'completed':
                raise ArithmeticError('lean worker did not complete')
            actual = json.loads((wd/'result.json').read_text())['mapping']
            duplicate = box_key(actual['matrix']) in completed.get(tuple(centre['point']),set())
            search = PointedQuarticSearch(state=state,
                centre={'coefficients': actual['centre']['representative']},
                coordinate_policy=actual['coordinate_policy'])
            backend.validate_map(search, actual)
            rows.append(dict(centre_index=ci, policy=policy,
                lean_wall_seconds=supervisor['wall_seconds'], duplicate_completed_box=duplicate,
                identity_verified=True, input_sha256=sha(wd/'input.json'),
                result_sha256=sha(wd/'result.json'), supervisor_sha256=sha(wd/'supervisor.json')))
            print('EXACT_LEAN_MAP', ci, policy, supervisor['wall_seconds'], flush=True)
    sources = {str(p.relative_to(ROOT)):sha(p) for p in
        [Path(__file__), CAS/'pointed_box_equivalence.py', worker, CAS/'lean_factor_free_pari_mapping.sage',
         CAS/'lean_preconditioned_full_pari_mapping.sage']}
    for path in [parent/'protocol.json',parent/'terminal.json',parent/'verified.json',parent/'epoch-00/landscape/selection.json',*sorted((parent/'epoch-00').glob('chart-*.json'))]:
        inputs[str(path.relative_to(ROOT))] = sha(path)
    if any(sha(ROOT/p) != h for p,h in inputs.items()):
        raise ArithmeticError('retained inputs changed')
    checkpoint(output/'report.json', dict(status='PASS_FIXED_COORDINATE_EXPOSURE_AUDIT',
        point_searches=0, rows=rows, inputs=inputs, sources=sources,
        distinct_uncovered_boxes=sum(not r['duplicate_completed_box'] for r in rows),
        lean_map_wall_seconds=sum(r['lean_wall_seconds'] for r in rows),
        claim_boundary='Fixed first16 centres only, exact signed-permutation box comparison including infinity at height125000. No point search, existence claim, or claim about untested centres.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.parent.resolve(), args.output.resolve())
