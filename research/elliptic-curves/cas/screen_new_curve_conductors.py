#!/usr/bin/env python3
"""Checkpointed factor discovery for a frozen post-discovery conductor screen."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT/'artifacts/local/elliptic-curves/conductor-record-screen-v1'
INPUT = ROOT/'artifacts/generated-results/elliptic-curves/inventory201_conductor_bounds_v1.json'
ECM = '/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/ecm'


def put(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')
    tmp.replace(path)


def prp(n):
    p = subprocess.run(['gp', '-fq'], input=f'print(ispseudoprime({n}));\n',
                       capture_output=True, text=True, timeout=10)
    if p.returncode or p.stdout.strip() not in ('0', '1'):
        raise ArithmeticError('probable-prime screening failed')
    return p.stdout.strip() == '1'


def prepare():
    data = json.loads(INPUT.read_text())
    rows = [r for r in data['rows'] if not r['catalogue_matches']]
    selected = set()
    for rank in range(22, 28):
        candidates = sorted((r for r in rows if r['rank_lower_bound'] >= rank),
                            key=lambda r: int(r['conductor_upper_bound']))
        selected.update(r['id'] for r in candidates[:3 if rank == 27 else 2])
    protocol = {
        'input': str(INPUT.relative_to(ROOT)),
        'input_sha256': hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'all_unpublished_ids': [r['id'] for r in rows],
        'priority_ids': sorted(selected), 'workers': 2,
        'pari_seconds_per_priority': 15, 'ecm_seconds_per_priority': 120,
        'ecm_schedule': [2000]*25 + [11000]*50 + [50000]*75 + [250000]*75 + [1000000]*50,
        'sigma': '2000000 + priority roster index*10000 + attempt',
        'ecm_binary': ECM, 'ecm_binary_sha256': hashlib.sha256(Path(ECM).read_bytes()).hexdigest(),
        'scope': 'After ICARM628, screen residual primality for all194 unpublished curves; '
                 'factor the two smallest conductor bounds at each rank22..27 and all three '
                 'unsubmitted rank27 curves. No point search or rank change. Discovery factors '
                 'are only inputs to a separate exact certificate checker.',
        'failure': 'Checkpoint partial factors; a PRP is not a proved prime, and incomplete '
                   'factorization is not an exact conductor or an exclusion by an upper bound.'}
    path = WORK/'factor_protocol.json'
    if path.exists():
        assert json.loads(path.read_text()) == protocol
    else:
        put(path, protocol)
    print('Prepared', len(rows), 'primality screens and', len(selected), 'factor priorities')


def run():
    protocol = json.loads((WORK/'factor_protocol.json').read_text())
    assert hashlib.sha256(INPUT.read_bytes()).hexdigest() == protocol['input_sha256']
    assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest() == protocol['script_sha256']
    rows = {r['id']: r for r in json.loads(INPUT.read_text())['rows']}
    for identifier in protocol['all_unpublished_ids']:
        path = WORK/'factors'/identifier/'initial.json'
        if not path.exists():
            n = int(rows[identifier]['remaining_cofactor'])
            put(path, {'id': identifier, 'factors': [str(n)] if n > 1 else [],
                       'residual_is_probable_prime': prp(n) if n > 1 else False})

    def worker(item):
        index, identifier = item
        directory = WORK/'factors'/identifier
        path = directory/'state.json'
        state = json.loads(path.read_text()) if path.exists() else {
            'id': identifier, 'factors': json.loads((directory/'initial.json').read_text())['factors'],
            'attempts': 0, 'ecm_seconds': 0, 'status': 'RUNNING'}
        if state['status'] != 'RUNNING':
            return
        factors = list(map(int, state['factors']))
        def checkpoint():
            assert math.prod(factors) == int(rows[identifier]['remaining_cofactor'])
            state['factors'] = list(map(str, sorted(factors)))
            put(path, state)
        if not (directory/'pari.json').exists():
            code = f'default(parisizemax,500000000);\nsetrand(628);\nprint(factorint({math.prod(factors)},9));\n'
            (directory/'pari.gp').write_text(code)
            start = time.monotonic()
            try:
                p = subprocess.run(['gp','-fq',str(directory/'pari.gp')],capture_output=True,
                                   text=True,timeout=protocol['pari_seconds_per_priority'])
                stdout, stderr, status = p.stdout, p.stderr, 'RETURNED'
            except subprocess.TimeoutExpired as e:
                stdout, stderr, status = e.stdout or b'', e.stderr or b'', 'TIMEOUT'
            stdout = stdout.decode() if isinstance(stdout,bytes) else stdout
            stderr = stderr.decode() if isinstance(stderr,bytes) else stderr
            (directory/'pari.stdout').write_text(stdout)
            (directory/'pari.stderr').write_text(stderr)
            parsed = [(int(a),int(b)) for a,b in re.findall(r'(\d+)\s*,\s*(\d+)',stdout)]
            if parsed and math.prod(a**b for a,b in parsed) == math.prod(factors):
                factors = [a for a,b in parsed for _ in range(b)]
            checkpoint()
            put(directory/'pari.json', {'status':status,'elapsed_seconds':time.monotonic()-start})
        for attempt in range(state['attempts'],len(protocol['ecm_schedule'])):
            composite = next((n for n in sorted(factors) if not prp(n)),None)
            if composite is None:
                state['status'] = 'PROBABLE_PRIME_FACTORIZATION'
                break
            if state['ecm_seconds'] >= protocol['ecm_seconds_per_priority']:
                break
            argv = [ECM,'-sigma',str(2000000+index*10000+attempt),'-maxmem','256',
                    str(protocol['ecm_schedule'][attempt])]
            start = time.monotonic()
            try:
                p = subprocess.run(argv,input=str(composite)+'\n',capture_output=True,text=True,
                    timeout=min(30,protocol['ecm_seconds_per_priority']-state['ecm_seconds']))
                stdout,stderr,status = p.stdout,p.stderr,'RETURNED'
            except subprocess.TimeoutExpired as e:
                stdout,stderr,status = e.stdout or b'',e.stderr or b'','TIMEOUT'
            stdout = stdout.decode() if isinstance(stdout,bytes) else stdout
            stderr = stderr.decode() if isinstance(stderr,bytes) else stderr
            elapsed = time.monotonic()-start
            (directory/f'ecm_{attempt:04d}.stdout').write_text(stdout)
            (directory/f'ecm_{attempt:04d}.stderr').write_text(stderr)
            put(directory/f'ecm_{attempt:04d}.json', {'argv':argv,'input':str(composite),
                'status':status,'elapsed_seconds':elapsed})
            for text in re.findall(r'Factor found in step \d+:\s*(\d+)',stdout):
                q = math.gcd(int(text),composite)
                if 1 < q < composite:
                    factors.remove(composite)
                    factors.extend([q,composite//q])
                    break
            state['attempts'] = attempt+1
            state['ecm_seconds'] += elapsed
            checkpoint()
        state['status'] = ('PROBABLE_PRIME_FACTORIZATION' if all(prp(n) for n in factors)
                           else 'BUDGET_EXHAUSTED')
        checkpoint()
        print(identifier,state['status'],[len(str(n)) for n in sorted(factors)],flush=True)
    with ThreadPoolExecutor(max_workers=protocol['workers']) as pool:
        list(pool.map(worker,enumerate(protocol['priority_ids'])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage',choices=['prepare','run'])
    args = parser.parse_args()
    prepare() if args.stage == 'prepare' else run()
