#!/usr/bin/env python3
"""Bounded branch-character audit; necessary residue conditions, never points."""
import argparse
from collections import defaultdict, Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import resource
import time

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1/rational-fibres.json'


def audit():
    resource.setrlimit(resource.RLIMIT_CPU, (15, 15))
    resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
    start = time.monotonic()
    rows = json.loads(SOURCE.read_text())['rows']
    assert len(rows) == 132 and all(r['smooth'] for r in rows)
    assert [r['t'] for r in rows] == list(range(131)) + [None]
    points = []
    for i, r in enumerate(rows):
        codes = r['rational_codes']
        assert len(codes) in (0, 1, 3)
        assert len(codes) == len(r['rational_roots'])
        assert all(0 < c < 2**17 for c in codes)
        if len(codes) == 3:
            assert codes[0] ^ codes[1] ^ codes[2] == 0
        points.extend((i, root, code) for root, code in zip(r['rational_roots'], codes))
    assert len(points) == len({p[2] for p in points}) == 110
    by_code = {p[2]: j for j, p in enumerate(points)}
    pairs = defaultdict(list)
    triples = set()
    for i, j in combinations(range(len(points)), 2):
        if points[i][0] == points[j][0]:
            continue
        code = points[i][2] ^ points[j][2]
        pairs[code].append((i, j))
        k = by_code.get(code)
        if k is not None and len({points[v][0] for v in (i, j, k)}) == 3:
            triples.add(tuple(sorted((i, j, k))))
    quartics = set()
    for bucket in pairs.values():
        for left, right in combinations(bucket, 2):
            word = tuple(sorted(left + right))
            if len({points[v][0] for v in word}) == 4:
                quartics.add(word)
    words = sorted(triples | quartics)
    supports = [{points[v][0] for v in word} for word in words]
    rank_two = set()
    for i, j in combinations(range(len(words)), 2):
        union = supports[i] | supports[j]
        if len(union) <= 4:
            rank_two.add(tuple(sorted(union)))
    def present(word):
        return [{'t': rows[points[v][0]]['t'], 'root': points[v][1], 'code': points[v][2]} for v in word]
    return {
        'schema': 'q80-rational-branch-code-audit-v1',
        'scope': 'Four distinct rational reductions at p=131; necessary reciprocity patterns only',
        'source': str(SOURCE.relative_to(ROOT)),
        'source_sha256': sha256(SOURCE.read_bytes()).hexdigest(),
        'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'limits': {'cpu_seconds': 15, 'address_space_bytes': 1024**3},
        'rational_sites': 132, 'nonzero_root_codes': len(points),
        'weight_one_or_two_words': 0,
        'weight_three_words': len(triples), 'weight_four_words': len(quartics),
        'rank_two_supports_size_at_most_four': [list(s) for s in sorted(rank_two)],
        'words': [present(w) for w in words],
        'word_count_by_support_size': dict(Counter(map(len, supports))),
        'new_sections_constructed': 0, 'global_target_complete': False,
        'elapsed_seconds': time.monotonic() - start,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', type=Path, required=True)
    args = parser.parse_args()
    assert not args.record.exists(), 'Refusing to overwrite evidence'
    result = audit()
    args.record.parent.mkdir(parents=True, exist_ok=True)
    with args.record.open('x') as f:
        json.dump(result, f, indent=2, sort_keys=True)
        f.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'words'}, sort_keys=True))
