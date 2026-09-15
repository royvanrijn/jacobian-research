#!/usr/bin/env python3
"""New finite branch-code joins, reusing completed higher-degree censuses."""
import argparse
from collections import defaultdict
from hashlib import sha256
import importlib.util
from itertools import combinations
import json
from pathlib import Path
import resource
import time

ROOT = Path(__file__).resolve().parents[2]
DEG = 'artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1/'
K1 = 'artifacts/generated-results/elkies-k3-q80-genus-one-k1-reciprocity-v1/'
K0 = 'artifacts/generated-results/elkies-k3-q80-genus-one-k0-quartic-gate-v1/'
RATIONAL = 'artifacts/generated-results/elkies-k3-q80-rational-branch-code-v1/audit.json'
RATIONAL_SCRIPT = 'elkies-k3/scripts/audit_q80_rational_branch_code.py'


def read(path):
    return json.loads((ROOT / path).read_text())


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def audit():
    resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
    resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
    start = time.monotonic()
    inputs = [DEG + n for n in ['rational-fibres.json', 'quadratic-fibres.json', 'input.json',
                              'result.json', 'independent-replay.json']]
    inputs += [K1 + n for n in ['cubic-census.json', 'result.json', 'independent-replay.json']]
    inputs += [K0 + n for n in ['cubic-census.json', 'quartic-comparison.json', 'quartic-input.json',
                              'quartic-producer.json', 'quartic-replay.json', 'independent-replay.json',
                              'result.json']]
    inputs += [RATIONAL, RATIONAL_SCRIPT, str(Path(__file__).relative_to(ROOT))]
    # Authenticate the retained summaries against their original result records.
    for prefix, names in [(DEG, ['rational-fibres.json', 'quadratic-fibres.json']),
                          (K1, ['cubic-census.json']), (K0, ['cubic-census.json'])]:
        result = read(prefix + 'result.json')
        assert result['status'] == 'PASS'
        for name in names:
            assert digest(prefix + name) == result['records'][name]
        assert read(prefix + 'independent-replay.json')['status'] == 'PASS'
    rat = read(DEG + 'rational-fibres.json')['rows']
    quad_all = read(DEG + 'quadratic-fibres.json')['rows']
    assert len(rat) == 132 and len(quad_all) == 8515
    assert all(r['smooth'] for r in rat)
    quad = [q for q in quad_all if q['smooth']]
    assert len(quad) == 8514
    codes = {c for r in rat for c in r['rational_codes']}
    assert len(codes) == sum(len(r['rational_codes']) for r in rat) == 110
    assert 0 not in codes
    for row in rat + quad:
        values = row.get('rational_codes', row.get('norm_codes'))
        assert len(values) in (0, 1, 3) and len(values) == len(set(values))
        assert 0 not in values
        if len(values) == 3:
            assert values[0] ^ values[1] ^ values[2] == 0
    old = read(RATIONAL)
    assert old['source_sha256'] == digest(DEG + 'rational-fibres.json')
    assert old['script_sha256'] == digest(RATIONAL_SCRIPT)
    spec = importlib.util.spec_from_file_location('rational_branch_audit', ROOT / RATIONAL_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    fresh = mod.audit()
    for key in ['words', 'weight_three_words', 'weight_four_words',
                'rank_two_supports_size_at_most_four']:
        assert fresh[key] == old[key]
    assert not fresh['rank_two_supports_size_at_most_four']
    pair_spans = defaultdict(set)
    pairs = list(combinations(range(132), 2))
    for index, (i, j) in enumerate(pairs):
        left, right = [0] + rat[i]['rational_codes'], [0] + rat[j]['rational_codes']
        span = {a ^ b for a in left for b in right}
        assert len(span) == len(left) * len(right)  # no degree-two relation
        for code in span - {0}:
            pair_spans[code].add(index)
    signatures = set()
    mixed = []
    for q in quad:
        values = q['norm_codes']
        if len(values) != 3:
            continue
        signature = tuple(sorted(values))
        assert signature not in signatures
        signatures.add(signature)
        a, b = values[:2]
        for index in pair_spans[a] & pair_spans[b]:
            mixed.append({'rational_indices': list(pairs[index]), 'quadratic_t': q['t']})
    assert len(signatures) == 1456 and not mixed
    c1, c0 = read(K1 + 'cubic-census.json'), read(K0 + 'cubic-census.json')
    for c in [c1, c0]:
        assert c['status'] == 'COMPLETE' and c['orbits'] == 749320
        assert c['smooth_split'] == 125595 and c['nodal'] == 1
    zeros = [r for r in c1['zero_norm_rows'] if r['smooth']]
    assert len(zeros) == 7
    for row in zeros:
        values = row['norm_codes']
        assert values.count(0) == 1 and len(set(values)) == 2
        assert not (set(values) & codes)
    assert c0['signature_rows'] == []
    # Check every retained quartic checkpoint; do not repeat the 73-million-orbit run.
    aggregates = []
    for kind in ['producer', 'replay']:
        aggregate = read(K0 + 'quartic-' + kind + '.json')
        assert aggregate['orbits'] == 73620690 and aggregate['smooth_split'] == 12267623
        parts = []
        for name, expected_hash in sorted(aggregate['checkpoints'].items()):
            path = K0 + name
            assert digest(path) == expected_hash
            packet = read(path)
            assert packet['input_sha256'] == digest(K0 + 'quartic-input.json')
            part = packet['result']
            assert part['status'] == 'COMPLETE_RANGE'
            assert part['orbits'] == (part['pair_end'] - part['pair_begin']) * 131**2
            assert part['candidate_rows'] == [] and part['first_failed_character'][-1] == 0
            parts.append(part)
            inputs.append(path)
        parts.sort(key=lambda p: p['pair_begin'])
        assert [(p['pair_begin'], p['pair_end']) for p in parts] == [(i, i + 429) for i in range(0, 4290, 429)]
        assert sum(p['orbits'] for p in parts) == 73620690
        aggregates.append(parts)
    for a, b in zip(*aggregates):
        for key in ['orbits', 'smooth_split', 'first_failed_character', 'split_character_checksum', 'candidate_rows']:
            assert a[key] == b[key]
    return {
        'schema': 'q80-good-branch-code-v1', 'status': 'PASS',
        'scope': 'Reduced branch divisor of degree two or four at131, disjoint from the reduced parent discriminant',
        'bindings': {path: digest(path) for path in sorted(set(inputs))},
        'new_calculations': {'rational_short_words': [len(w) for w in fresh['words']],
                             'rational_pair_spans': len(pairs),
                             'smooth_full_split_quadratic_spans': len(signatures),
                             'mixed_rank_two_matches': mixed},
        'inherited_cubic_zero_rows': len(zeros),
        'inherited_cubic_signature_matches': 0, 'inherited_quartic_zero_pairs': 0,
        'kernel_dimension_upper_bounds': {'1+1': 0, '2': 0, '1+1+1+1': 1,
                                         '1+1+2': 1, '2+2': 1, '1+3': 1, '4': 1},
        'higher_degree_censuses_recomputed': False,
        'rank_bound_on_twist_claimed': False, 'positive_target_complete': False,
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
    print(json.dumps({k: v for k, v in result.items() if k not in ['bindings', 'new_calculations']}, sort_keys=True))
