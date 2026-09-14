#!/usr/bin/env python3
"""Test nonzero 2-torsion over every retained Q80 rational-bisection branch field.

This is a new finite arithmetic gate, not a new bisection or point search.
At a simple finite reduction root of q, a root of the monic parent cubic
over Q[t]/q would reduce to a root over F_p. A root-free fibre excludes it.
Only the 11952 atlas is used; no exceptional specialization enters selection.
"""
import argparse
from fractions import Fraction
import gzip
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import resource
import time

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / 'artifacts/generated-results/elkies-k3-q80-rational-branch-torsion-v1'
ATLAS = 'artifacts/generated-results/elkies-k3-r17-smooth-character-replay-inputs-v1.json.gz'
PARENTS = 'artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/input.json'
HELPER = 'elkies-k3/scripts/verify_r17_smooth_character_witness.py'
CHECKER = 'elkies-k3/scripts/verify_q80_rational_branch_torsion.py'
SPEC = importlib.util.spec_from_file_location('smooth_character', ROOT / HELPER)
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix == '.gz' else data)


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def primes():
    result, n = [], 101
    while len(result) < 64:
        if all(n % d for d in range(2, int(n**0.5)+1)):
            result.append(n)
        n += 1
    return result


def freeze(out):
    parent = next(r for r in read(ROOT/PARENTS)['parents'] if r['name'] == 'alternate-q80')
    paths = [ATLAS, PARENTS, HELPER, parent['source']]
    assert digest(ROOT/ATLAS) == V.PACKET_SHA256
    assert digest(ROOT/parent['source']) == parent['source_sha256']
    write_new(out/'input.json', {
        'schema': 'q80-rational-branch-torsion-input-v1',
        'bindings': {p: digest(ROOT/p) for p in paths},
        'producer_sha256': digest(Path(__file__)), 'checker_sha256': digest(ROOT/CHECKER),
        'parent': parent, 'chart': '11952', 'count': 39147,
        'primes': primes(), 'cpu_seconds_per_mode': 120, 'address_space_bytes': 4*1024**3,
        'selection': 'Every retained smooth rational-bisection translation class, in the existing atlas order; first64 primes at least101. No preview or exceptional-point selection.',
        'acceptance': 'A simple finite branch root at an integral parent prime with a smooth root-free monic cubic. Each class needs one witness; survivors are UNKNOWN.',
        'boundary': 'The branch-field gate and resulting common-quartic descent do not assert exact twist rank or exclude other branch fields, genera or parents.',
    })
    print('Frozen39147 branch fields and64 primes', flush=True)


def evaluate(coeff, t, p):
    result = 0
    for c in reversed(coeff):
        result = (result*t+c) % p
    return result


def run(out):
    start = time.process_time()
    packet = read(out/'input.json')
    assert digest(Path(__file__)) == packet['producer_sha256']
    assert digest(ROOT/CHECKER) == packet['checker_sha256']
    for path, expected in packet['bindings'].items():
        assert digest(ROOT/path) == expected, path
    atlas = read(ROOT/ATLAS)['atlases'][0]
    assert atlas['source'] == V.SPECS[0]
    V.verify_atlas(atlas, 39147)
    branches = [V.quadratic_atom(row[3], row[4]) for row in atlas['rows']]
    coeff = [[Fraction(c) for c in packet['parent'][key]] for key in ('A','B')]
    remaining, stages = set(range(39147)), []
    for p in packet['primes']:
        if not remaining:
            break
        if any(c.denominator % p == 0 for line in coeff for c in line):
            stages.append({'prime': p, 'status': 'DEFER_NONINTEGRAL_PARENT'})
            continue
        A, B = [[c.numerator*pow(c.denominator, -1, p) % p for c in line] for line in coeff]
        no_roots = []
        for t in range(p):
            a, b = evaluate(A,t,p), evaluate(B,t,p)
            if (4*a**3+27*b*b) % p and all((x**3+a*x+b) % p for x in range(p)):
                no_roots.append(t)
        bad = set(no_roots)
        sqrt = {x*x % p: x for x in range(p)}
        witnesses = []
        for i in sorted(remaining):
            c,b,a = [x % p for x in branches[i]]
            discriminant = (b*b-4*a*c) % p
            if not a or not discriminant or discriminant not in sqrt:
                continue
            s = sqrt[discriminant]
            roots = sorted({(-b+s)*pow(2*a,-1,p) % p, (-b-s)*pow(2*a,-1,p) % p})
            t = next((t for t in roots if t in bad), None)
            if t is not None:
                witnesses.append([i,t])
        remaining.difference_update(i for i,t in witnesses)
        stage = {'prime':p, 'input_sha256':digest(out/'input.json'),
                 'root_free_smooth_fibres':no_roots, 'witnesses':witnesses,
                 'remaining_indices':sorted(remaining)}
        name = 'prime-%d.json' % p
        write_new(out/name,stage)
        meta = {'prime':p, 'status':'TESTED', 'excluded':len(witnesses),
                'remaining':len(remaining), 'path':name, 'sha256':digest(out/name)}
        stages.append(meta)
        print(json.dumps({k:meta[k] for k in ('prime','excluded','remaining')}),flush=True)
    result = {'schema':'q80-rational-branch-torsion-result-v1',
              'status':'PASS' if not remaining else 'INCOMPLETE',
              'input_sha256':digest(out/'input.json'), 'stages':stages,
              'count':39147, 'excluded':39147-len(remaining), 'remaining_indices':sorted(remaining),
              'positive_target_complete':False, 'whole_family_genus_bound':'UNKNOWN',
              'cpu_seconds':round(time.process_time()-start,6)}
    write_new(out/'result.json',result)
    print(json.dumps({k:result[k] for k in ('status','excluded','cpu_seconds')}),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['freeze','run'])
    parser.add_argument('--output',type=Path,default=DEFAULT)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(120,125))
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    (freeze if args.mode == 'freeze' else run)(args.output)
