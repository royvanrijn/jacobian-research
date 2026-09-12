#!/usr/bin/env python3
"""Bounded exact full-frame shell traversal, preserving its incomplete prefix.

Reuses IntegerExactParity's rational LDL decomposition and integer branch
bounds. The adaptation retains signed vectors and a resumable frontier;
it never mistakes a parity-minimum bank for the requested full shells.
"""
import gzip
import hashlib
import json
from math import isqrt
from pathlib import Path
import resource
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'artifacts/local/elliptic-curves/marked-two-class-v1'
sys.path.insert(0, str(OUT/'source-snapshots/elliptic-curves/cas'))
from visibility_lattice_fast import IntegerExactParity


def write(name, data):
    path = OUT/name
    with path.with_suffix(path.suffix+'.tmp').open('w') as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write('\n')
    path.with_suffix(path.suffix+'.tmp').replace(path)


def main():
    protocol = json.loads((OUT/'protocol.json').read_text())
    assert hashlib.sha256((OUT/'geometry.json').read_bytes()).hexdigest() == protocol['geometry_sha256']
    own = 'elliptic-curves/rank-jump/enumerate_marked_two_class.py'
    assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest() == protocol['sources'][own]
    if (OUT/'enumeration.json').exists() or (OUT/'shell-prefix.jsonl.gz').exists():
        raise FileExistsError('retain the prior traversal; no automatic continuation')
    resource.setrlimit(resource.RLIMIT_CPU, (3600, 3605))
    geometry = json.loads((OUT/'geometry.json').read_text())
    G = geometry['frame_gram']
    oracle = IntegerExactParity(G)
    n = oracle.n
    Q = geometry['source_ns_gram']
    B = geometry['frame_basis_in_source_ns']
    F, O = geometry['fibre'], geometry['zero']
    walls = geometry['effective_wall_divisors']
    wall_duals = [[sum(Q[i][j]*s[j] for j in range(19)) for i in range(19)] for s in walls]
    words = [0]*n
    stack = []
    nodes = 0
    counts = {'10': 0, '14': 0}
    rejected = {'10': 0, '14': 0}
    best = {'10': [], '14': []}
    start = time.process_time()
    limit = 14*oracle.norm_scale
    max_nodes = protocol['enumeration']['total_nodes']
    stop = None

    def snapshot():
        return {'status': 'INCOMPLETE_PREFIX', 'nodes': nodes, 'counts': counts,
            'known_wall_rejections': rejected, 'cpu_seconds': time.process_time()-start,
            'active_coordinates': list(words), 'frontier': [list(s) for s in stack],
            'selected_prefixes': best, 'geometry_sha256': protocol['geometry_sha256']}

    def visit(i, used, stream):
        nonlocal nodes, stop
        if nodes >= max_nodes or time.process_time()-start >= 3590:
            stop = 'NODE_LIMIT' if nodes >= max_nodes else 'CPU_LIMIT'
            return False
        nodes += 1
        if nodes % 100000 == 0:
            write('enumeration-checkpoint.json', snapshot())
            stream.flush()
            print(json.dumps({'nodes': nodes, 'counts': counts, 'cpu': time.process_time()-start}), flush=True)
        if i < 0:
            norm, remainder = divmod(used, oracle.norm_scale)
            assert remainder == 0
            if norm not in (10, 14):
                return True
            key = str(norm)
            counts[key] += 1
            h = (norm-10)//4
            D = [2*O[j]+(h+4)*F[j]+sum(words[k]*B[k][j] for k in range(n)) for j in range(19)]
            wall = next((j for j, s in enumerate(wall_duals) if sum(a*b for a,b in zip(D,s)) < 0), None)
            row = {'norm': norm, 'reduced_frame_coordinates': list(words), 'native_ns_coordinates': D, 'negative_wall_index': wall}
            stream.write(json.dumps(row, separators=(',', ':'))+'\n')
            if wall is not None:
                rejected[key] += 1
                return True
            # Native coordinates use F_old, O_old+F_old, then MW17.
            order = [D[1], D[0]-D[1], sum(abs(x) for x in D), D]
            row['construction_cost_key'] = order
            best[key].append(row)
            best[key].sort(key=lambda r:r['construction_cost_key'])
            if len(best[key]) > 128:
                best[key].pop()
            return True
        b = oracle.shift_denominators[i]
        a = sum(c*words[j] for j,c in enumerate(oracle.shift_coefficients[i], i+1))
        weight = oracle.integer_weights[i]
        rad = isqrt((limit-used)//weight)
        lo, hi = -((rad+a)//b), (rad-a)//b
        stack.append([i, used, lo, hi])
        for x in range(lo, hi+1):
            words[i] = x
            stack[-1][2] = x
            if not visit(i-1, used+weight*(b*x+a)**2, stream):
                return False
        stack.pop()
        return True

    with gzip.open(OUT/'shell-prefix.jsonl.gz', 'xt') as stream:
        complete = visit(n-1, 0, stream)
    result = snapshot()
    result.update({'status': 'COMPLETE_SHELLS' if complete else 'INCOMPLETE_PREFIX',
        'stop_reason': stop, 'shells_complete': complete,
        'prefix_sha256': hashlib.sha256((OUT/'shell-prefix.jsonl.gz').read_bytes()).hexdigest(),
        'cpu_seconds': time.process_time()-start,
        'boundary': 'Every retained vector is exact. Known walls give exact individual rejections. Surviving vectors still require curve equations and irreducibility; no generic-translation deduplication was used.'})
    write('enumeration.json', result)
    write('frozen-candidates.json', {'geometry_sha256':protocol['geometry_sha256'],
        'enumeration_sha256':hashlib.sha256((OUT/'enumeration.json').read_bytes()).hexdigest(),
        'shells_complete':complete, 'rows':[row for key in ('10','14') for row in best[key]]})
    print(json.dumps({k:v for k,v in result.items() if k not in ('selected_prefixes','frontier','active_coordinates')}), flush=True)


if __name__ == '__main__':
    main()
