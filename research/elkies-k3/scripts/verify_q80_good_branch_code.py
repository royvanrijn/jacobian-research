#!/usr/bin/env python3
"""Independent short-word enumeration and two-dimensional subspace replay."""
import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import resource
import time

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / 'artifacts/generated-results/elkies-k3-q80-good-branch-code-v1/result.json'


def short_words(points):
    """Enumerate triples directly; a fourth entry is determined by its code."""
    by_code = {c: i for i, (_, _, c) in enumerate(points)}
    assert len(by_code) == len(points) and 0 not in by_code
    words = set()
    for i, j, k in combinations(range(len(points)), 3):
        if len({points[v][0] for v in [i, j, k]}) != 3:
            continue
        code = points[i][2] ^ points[j][2] ^ points[k][2]
        if code == 0:
            words.add((i, j, k))
        elif code in by_code:
            l = by_code[code]
            if l > k and points[l][0] not in {points[v][0] for v in [i, j, k]}:
                words.add((i, j, k, l))
    return sorted(words)


def double_supports(points, words):
    supports = [{points[v][0] for v in w} for w in words]
    return {tuple(sorted(a | b)) for a, b in combinations(supports, 2) if len(a | b) <= 4}


def planes(values):
    return {tuple(sorted((a, b, a ^ b))) for a, b in combinations(values - {0}, 2)}


def verify(path):
    started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
    resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
    result = json.loads(path.read_text())
    assert result['status'] == 'PASS'
    for name, expected in result['bindings'].items():
        assert sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
    base = ROOT / 'artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
    rat = json.loads((base / 'rational-fibres.json').read_text())['rows']
    quad = json.loads((base / 'quadratic-fibres.json').read_text())['rows']
    points = [(i, root, c) for i, row in enumerate(rat)
              for root, c in zip(row['rational_roots'], row['rational_codes'])]
    assert len(points) == 110
    words = short_words(points)
    assert len([w for w in words if len(w) == 3]) == 4
    assert len([w for w in words if len(w) == 4]) == 40
    assert not double_supports(points, words)
    retained = json.loads((ROOT / 'artifacts/generated-results/elkies-k3-q80-rational-branch-code-v1/audit.json').read_text())
    displayed = [[{'t': rat[points[v][0]]['t'], 'root': points[v][1], 'code': points[v][2]}
                  for v in word] for word in words]
    assert displayed == retained['words']
    rational_planes = set()
    for i, j in combinations(range(132), 2):
        left, right = {0, *rat[i]['rational_codes']}, {0, *rat[j]['rational_codes']}
        span = {a ^ b for a in left for b in right}
        assert len(span) == len(left) * len(right)
        rational_planes |= planes(span)
    quadratic_planes = set()
    for row in quad:
        if not row['smooth']:
            continue
        values = row['norm_codes']
        assert 0 not in values
        if len(values) == 3:
            assert values[0] ^ values[1] ^ values[2] == 0
            plane = tuple(sorted(values))
            assert len(set(values)) == 3 and plane not in quadratic_planes
            quadratic_planes.add(plane)
    assert len(quadratic_planes) == 1456
    assert not (quadratic_planes & rational_planes)
    k1 = ROOT / 'artifacts/generated-results/elkies-k3-q80-genus-one-k1-reciprocity-v1'
    k0 = ROOT / 'artifacts/generated-results/elkies-k3-q80-genus-one-k0-quartic-gate-v1'
    cubic = json.loads((k1 / 'cubic-census.json').read_text())
    assert cubic['orbits'] == 749320 and cubic['status'] == 'COMPLETE'
    rational_codes = {c for _, _, c in points}
    for row in cubic['zero_norm_rows']:
        if row['smooth']:
            assert row['norm_codes'].count(0) == 1
            assert not (set(row['norm_codes']) & rational_codes)
    assert json.loads((k0 / 'cubic-census.json').read_text())['signature_rows'] == []
    # Original hash-bound producer/replay checkpoints supply the degree-four input.
    compared = []
    for kind in ['producer', 'replay']:
        agg = json.loads((k0 / f'quartic-{kind}.json').read_text())
        pieces = []
        for name, expected in agg['checkpoints'].items():
            data = (k0 / name).read_bytes()
            assert sha256(data).hexdigest() == expected
            r = json.loads(data)['result']
            assert not r['candidate_rows'] and r['first_failed_character'][-1] == 0
            pieces.append(r)
        pieces.sort(key=lambda r: r['pair_begin'])
        assert [(r['pair_begin'], r['pair_end']) for r in pieces] == [(i, i+429) for i in range(0,4290,429)]
        assert sum(r['orbits'] for r in pieces) == 73620690
        compared.append([(r['orbits'],r['smooth_split'],r['first_failed_character'],r['split_character_checksum']) for r in pieces])
    assert compared[0] == compared[1]
    synthetic = [(s, j, c) for s, row in enumerate([[1,2,3],[4,8,12],[5,10,15]])
                 for j, c in enumerate(row)]
    assert double_supports(synthetic, short_words(synthetic)) == {(0,1,2)}
    assert tuple(sorted((1,4,5))) in planes({0,1,2,3,4,5,6,7})
    bounds = {'1+1':0,'2':0,'1+1+1+1':1,'1+1+2':1,'2+2':1,'1+3':1,'4':1}
    assert result['kernel_dimension_upper_bounds'] == bounds
    assert result['higher_degree_censuses_recomputed'] is False
    assert result['rank_bound_on_twist_claimed'] is False and result['positive_target_complete'] is False
    return {'status':'PASS','result_sha256':sha256(path.read_bytes()).hexdigest(),
            'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'independent_new_linear_algebra':True,'higher_degree_censuses_recomputed':False,
            'written_descent_formally_verified':False,'positive_target_complete':False,
            'rational_short_words':len(words),'rational_pair_planes':len(rational_planes),
            'quadratic_planes':len(quadratic_planes),'synthetic_rank_two_control':True,
            'elapsed_seconds':time.monotonic()-started}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=DEFAULT)
    parser.add_argument('--record', type=Path)
    args = parser.parse_args()
    answer = verify(args.input)
    if args.record:
        with args.record.open('x') as f:
            json.dump(answer, f, indent=2, sort_keys=True)
            f.write('\n')
    print(json.dumps(answer, sort_keys=True))
