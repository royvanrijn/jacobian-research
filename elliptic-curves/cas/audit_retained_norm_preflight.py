#!/usr/bin/env python3
"""Frozen retained-data evaluation of the reusable factor-free preflight."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from research_runtime import norm_ramification as gate

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
WORK = ROOT/'artifacts/local/elliptic-curves/retained-norm-preflight-v1'
INPUT = ART/'retained_norm_preflight_inputs_v1.json'
PROTOCOL = ART/'retained_norm_preflight_protocol_v1.json'
OUT = ART/'retained_norm_preflight_v1.json'
SOURCES = [ART/(x+'.json') for x in [
    'rank_jump_fresh_norm_projection_v1',
    'rank_jump_fresh_governing_panel_inputs_v1',
    'rank_jump_retained_norm_batch_capacity_inputs_v1']]


def read(path):
    return json.loads(path.read_text())


def hashed(path):
    return sha256(path.read_bytes()).hexdigest()


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def projection():
    panel, models, batch = map(read, SOURCES)
    models = {r['token']: r['model'] for r in models['cases']}
    cases = []
    for row in panel['rows']:
        assert row['status'] == 'PASS'
        a1, a2, a3, a4, a6 = models[row['token']]
        assert all(gate.rational(x) == 0 for x in [a1, a2, a3])
        cases.append({'id': row['token'], 'cubic_ascending': [a6, a4, '0', '1'],
                      'elements': [c['alpha_ascending'] for c in row['candidates']]})
    cases.append({'id': 'retained-reference-296', 'cubic_ascending': batch['cubic_ascending'],
                  'elements': [r['alpha_ascending'] for r in batch['relations']]})
    assert len(cases) == 12 and sum(len(c['elements']) for c in cases) == 428
    return {'schema': 'elliptic-curves.retained-norm-preflight-input.v1', 'cases': cases}


def freeze():
    projected = projection()
    write_new(INPUT, projected)
    paths = [*SOURCES, Path(__file__).resolve(), Path(gate.__file__).resolve(), INPUT]
    write_new(PROTOCOL, {
        'schema': 'elliptic-curves.retained-norm-preflight-protocol.v1',
        'bindings': {str(p.relative_to(ROOT)): hashed(p) for p in paths},
        'population': 'All132 original panel generators plus all296 fixed-box reference generators; no new candidates.',
        'rule': 'Peel only generators with an isolated nonsquare norm remainder. Repeat until stable. Retain shared-support cancellations.',
        'limits': {'wall_seconds': 60, 'cases': 12, 'generators': 428},
        'endpoint': 'Forced-zero coefficient coordinates and unresolved dimension cap, never curve exclusion or rational points.',
        'provenance': 'Retrospective retained dictionaries. Worker reads only cubics and element coefficients. Known outcomes informed motivation, not a blind discrimination test.',
        'following_campaign': None})


def calculate():
    protocol = read(PROTOCOL)
    for path, digest in protocol['bindings'].items():
        assert hashed(ROOT/path) == digest, path
    assert read(INPUT) == projection()
    rows = [{'id': c['id'], **gate.preflight(c['cubic_ascending'], c['elements'])}
            for c in read(INPUT)['cases']]
    return {'schema': 'elliptic-curves.retained-norm-preflight.v1', 'status': 'PASS',
            'protocol_sha256': hashed(PROTOCOL), 'rows': rows,
            'totals': {'generators': 428, 'forced_zero': sum(r['forced_zero_count'] for r in rows),
                       'unresolved': sum(len(r['unresolved_indices']) for r in rows)},
            'scope': 'Necessary dictionary-capacity restrictions only; UNKNOWN is not admissibility. Frozen point campaigns and scores are unchanged.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['freeze', 'build', 'check'])
    mode = parser.parse_args().mode
    if mode == 'freeze':
        freeze()
        print('FROZEN 12 cases / 428 generators')
    else:
        result = calculate()
        if mode == 'check':
            assert result == read(OUT), 'retained preflight differs'
        else:
            write_new(OUT, result)
        print(result['status'], result['totals'])
        for row in result['rows']:
            print(row['id'], row['forced_zero_count'], '/', row['input_generator_count'], 'forced zero')
