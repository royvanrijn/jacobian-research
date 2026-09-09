#!/usr/bin/env sage-python
"""Factor-free support audit of the immutable 567-reduction multiplier bank.

This does not rerun ideal reduction, enlarge the bank, or compute a class
group. Private unramified rational primes separate primitive ideal divisors
modulo rational principal ideals. Each process is capped at 25 seconds.
"""
import hashlib, json, signal
from pathlib import Path
from sage.all import QQ, ZZ, pari, matrix, gcd, prod

signal.alarm(25)
ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
OUT = ART / 'det1092_saved_unit_relations_v1'
BANK = ART / 'rank_jump_curve302_arithmetic_unit_class_v1.json'
ARITH = ART / 'rank_jump_curve302_strict_constructor_arithmetic_v1.json'

def read(path): return json.loads(path.read_text())
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def save(name, value):
    OUT.mkdir(exist_ok=True)
    path = OUT / name
    if path.exists(): assert read(path) == value
    else:
        with path.open('x') as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write('\n')

bank = read(BANK)
sources = [ROOT / entry['path'] for entry in bank['trials']]
assert len(sources) == 7
for path, entry in zip(sources, bank['trials']):
    assert sha(path) == entry['sha256']
save('protocol.json', {
    'classification': 'new bounded relation audit of frozen generic-only inputs',
    'inputs': {str(p.relative_to(ROOT)): sha(p) for p in [BANK, ARITH, *sources]},
    'script_sha256': sha(Path(__file__)),
    'limits': {'seconds': 25, 'saved_trials': 567, 'new_reductions': 0,
               'factorizations': 0, 'class_groups': 0, 'point_searches': 0},
    'test': 'Peel primitive reduced ideals having a private unramified rational prime. Duplicate ideals and residual starting-ideal relations are retained, not assumed trivial.'})

groups = {}
starting_norms = []
trials = []
for source_index, path in enumerate(sources):
    data = read(path)
    I = matrix(QQ, pari(data['starting_ideal_hnf']).sage())
    assert I.nrows() == I.ncols() == 3 and all(v in ZZ for v in I.list())
    starting_norms.append(ZZ(abs(I.det())))
    assert len(data['trials']) == 81
    for trial_index, item in enumerate(data['trials']):
        J = matrix(QQ, pari(item['reduced_ideal_hnf']).sage())
        assert J.nrows() == J.ncols() == 3 and all(v in ZZ for v in J.list())
        content = gcd(ZZ(v) for v in J.list())
        assert content > 0
        primitive = matrix(ZZ, J / content)
        key = tuple(primitive.list())
        if key not in groups:
            groups[key] = {'id': len(groups), 'matrix': primitive,
                           'norm': ZZ(abs(primitive.det())), 'occurrences': []}
        group = groups[key]
        group['occurrences'].append([source_index, trial_index])
        trials.append({'source': source_index, 'trial': trial_index,
                       'group': group['id'], 'content': str(content)})

ordered = sorted(groups.values(), key=lambda g: g['id'])
blocked_support = abs(ZZ(read(ARITH)['field_discriminant'])) * prod(starting_norms)
active = {g['id']: g['norm'] for g in ordered if g['norm'] != 1}
peeling = []
while active:
    product = prod(active.values())
    found = False
    for index in sorted(active):
        remaining = active[index]
        blocker = blocked_support * (product // remaining)
        removed = []
        while True:
            common = gcd(remaining, blocker)
            if common == 1: break
            removed.append(str(common))
            remaining //= common
        if remaining != 1:
            assert gcd(remaining, blocker) == 1
            peeling.append({'group': index, 'private_part': str(remaining),
                            'removed_gcds': removed})
            del active[index]
            found = True
            break
    if not found: break

report = {
    'status': 'PASS_FROZEN_BANK_SUPPORT_AUDIT',
    'classification': 'verified exact support computation; interpretation requires the primitive-ideal lemma',
    'saved_trials': len(trials), 'distinct_primitive_ideals': len(ordered),
    'primitive_unit_ideal_groups': [g['id'] for g in ordered if g['norm'] == 1],
    'cross_source_groups': [g['id'] for g in ordered
                            if len(set(pair[0] for pair in g['occurrences'])) > 1],
    'starting_norms': list(map(str, starting_norms)),
    'groups': [{'id': g['id'], 'norm': str(g['norm']),
                'matrix': [[str(v) for v in row] for row in g['matrix'].rows()],
                'occurrences': g['occurrences']} for g in ordered],
    'trials': trials, 'peeling': peeling, 'unpeeled_nonunit_groups': sorted(active),
    'boundary': 'No conclusion about units from unpeeled groups or cross-source repetitions without checking their exact ideal and multiplier relations.',
    'protocol_sha256': sha(OUT / 'protocol.json'),
}
save('support.json', report)
print(report['status'], 'groups', len(ordered), 'peeled', len(peeling),
      'unpeeled', len(active), 'cross_source', report['cross_source_groups'], flush=True)
