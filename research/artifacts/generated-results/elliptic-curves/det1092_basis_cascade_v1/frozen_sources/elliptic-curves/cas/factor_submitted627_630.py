#!/usr/bin/env python3
"""Checkpointed factor discovery for the four explicitly prioritized submissions."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / 'artifacts/local/elliptic-curves/submitted627-630-conductors-v1'


def run(row):
    identifier = row['icarm_id']
    prefix = WORK / f'pari_partial_{identifier}'
    if prefix.with_suffix('.json').exists():
        return json.loads(prefix.with_suffix('.json').read_text())
    script = ('default(parisizemax,536870912);setrand(627);\n'
              f'F=factorint({row["remaining_cofactor"]},9);\n'
              'for(i=1,matsize(F)[1],print(F[i,1]," ",F[i,2]," ",ispseudoprime(F[i,1])));\n')
    prefix.with_suffix('.gp').write_text(script)
    start = time.monotonic()
    try:
        p = subprocess.run(['gp', '-fq', '-s', '64000000'], input=script,
                           text=True, capture_output=True, timeout=60)
        stdout, stderr = p.stdout, p.stderr
        factors = [line.split() for line in stdout.splitlines() if len(line.split()) == 3]
        product = 1
        for q, e, _ in factors:
            product *= int(q)**int(e)
        if p.returncode or product != int(row['remaining_cofactor']):
            raise ArithmeticError('partial factorization reconstruction failed')
        result = {'status': 'FACTORS_RETURNED', 'icarm_id': identifier,
                  'factors': [{'factor': q, 'exponent': int(e), 'probable_prime': flag == '1'}
                              for q, e, flag in factors]}
    except subprocess.TimeoutExpired as e:
        stdout, stderr = e.stdout or b'', e.stderr or b''
        stdout = stdout.decode() if isinstance(stdout, bytes) else stdout
        stderr = stderr.decode() if isinstance(stderr, bytes) else stderr
        result = {'status': 'TIMEOUT', 'icarm_id': identifier, 'factors': []}
    result['wall_seconds'] = time.monotonic()-start
    prefix.with_suffix('.stdout').write_text(stdout)
    prefix.with_suffix('.stderr').write_text(stderr)
    prefix.with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n')
    print(identifier, result, flush=True)
    return result


if __name__ == '__main__':
    rows = json.loads((WORK/'protocol.json').read_text())['roster']
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(run, rows))
