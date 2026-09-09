"""Bounded map construction, immutable receipts, and arithmetic-only replay."""
import hashlib
import json
from pathlib import Path
import time
from research_runtime.store import checkpoint
from research_runtime.supervisor import run, Limits
from pointed_quartic_search import PointedQuarticSearch
import pari_pointed_backend as backend
from preconditioned_map_worker import MAPPERS

CAS = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def require(test, message):
    if not test:
        raise ArithmeticError(message)


def obtain(epoch, ci, policy, model, points, centre, state, protocol, *, replay=False):
    """Return (map or None, receipt hash, resource-limit reason or None).

    Replay never starts a worker, including for incomplete maps. Every returned
    map is checked by the independent pointed-quartic identity verifier.
    """
    base = epoch/f'map-{ci:04d}-{policy}'
    receipt_path = base.with_suffix('.json')
    payload = dict(curve=list(map(str, model)), points=[list(map(str, p)) for p in points],
                   centre=json.loads(json.dumps(centre)), policy=policy)
    worker = CAS/'preconditioned_map_worker.py'
    mapper = CAS/MAPPERS[policy]
    command = [protocol['map_python'], '-python', str(worker),
               '--input', str(base/'input.json'), '--output', str(base/'result.json')]
    if not receipt_path.exists():
        require(not replay, 'missing map receipt during replay')
        if base.exists():
            base.rename(base.with_name(base.name+f'-incomplete-{time.time_ns()}'))
        base.mkdir()
        checkpoint(base/'input.json', payload)
        supervisor = run(command,
                         limits=Limits(wall_seconds=protocol['seconds_per_map'],
                                       rss_bytes=protocol['map_rss_bytes']),
                         log_path=base/'worker.log', result_path=base/'result.json')
        checkpoint(base/'supervisor.json', supervisor)
        outcome = supervisor['outcome']
        require(outcome in ('completed', 'strict_wall_timeout', 'strict_rss_limit'),
                'map worker failed: '+outcome)
        checkpoint(receipt_path, dict(status=outcome, input_sha256=sha(base/'input.json'),
                    supervisor_sha256=sha(base/'supervisor.json'),
                    result_sha256=sha(base/'result.json') if outcome == 'completed' else None,
                    worker_sha256=sha(worker), mapper_sha256=sha(mapper)))
        print('BOUNDED_MAP', ci, policy, outcome, flush=True)
    receipt = read(receipt_path)
    require(receipt['input_sha256'] == sha(base/'input.json') and
            read(base/'input.json') == payload, 'map input differs')
    require(receipt['worker_sha256'] == sha(worker) and receipt['mapper_sha256'] == sha(mapper),
            'map source differs')
    require(receipt['supervisor_sha256'] == sha(base/'supervisor.json'), 'map supervisor differs')
    supervisor = read(base/'supervisor.json')
    require(supervisor['command'] == command and
            supervisor['limits']['wall_seconds'] == protocol['seconds_per_map'] and
            supervisor['limits']['rss_bytes'] == protocol['map_rss_bytes'] and
            supervisor['outcome'] == receipt['status'] and
            supervisor['log_sha256'] == sha(base/'worker.log'), 'map execution binding differs')
    if receipt['status'] in ('strict_wall_timeout', 'strict_rss_limit'):
        require(receipt['result_sha256'] is None, 'censored map claims a result')
        return None, sha(receipt_path), receipt['status']
    require(receipt['status'] == 'completed' and supervisor['returncode'] == 0 and
            receipt['result_sha256'] == sha(base/'result.json') == supervisor['worker_result_sha256'],
            'successful map result required')
    result = read(base/'result.json')
    require(result['status'] == 'EXACT_MAP_CONSTRUCTED' and
            result['input_sha256'] == receipt['input_sha256'] and
            result['worker_sha256'] == receipt['worker_sha256'] and
            result['mapper_sha256'] == receipt['mapper_sha256'], 'map result binding differs')
    mapping = result['mapping']
    require(mapping['centre'] == payload['centre'], 'map centre differs')
    search = PointedQuarticSearch(state=state, centre={'coefficients': centre['representative']},
                                 coordinate_policy=mapping['coordinate_policy'])
    backend.validate_map(search, mapping)
    return mapping, sha(receipt_path), None
