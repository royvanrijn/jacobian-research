"""Bounded reading of the complete retained degree-two character tables."""
import collections
import hashlib
import itertools
import json
from pathlib import Path
import resource

resource.setrlimit(resource.RLIMIT_CPU, (10, 15))
ROOT = Path(__file__).resolve().parents[4]
OLD = ROOT / 'artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
rational = json.loads((OLD / 'rational-fibres.json').read_text())['rows']
quadratic = json.loads((OLD / 'quadratic-fibres.json').read_text())['rows']
points = [(r['t'], e, c) for r in rational
          for e, c in zip(r['rational_roots'], r['rational_codes'])]
codes = {v[2]: v[:2] for v in points}
assert len(codes) == len(points) == 110 and 0 not in codes
groups = collections.defaultdict(list)
for a, b in itertools.combinations_with_replacement(points, 2):
    c = a[2] ^ b[2]
    if c in codes or c in (0, 13412):
        groups[c].append(['rational', a[:2], b[:2]])
for row in rational + quadratic:
    for e, c in zip(row['roots'], row['norm_codes']):
        if c in codes or c in (0, 13412):
            groups[c].append(['quadratic', row['t'], e, row['smooth']])
exceptional = []
for c, (t, e) in codes.items():
    for row in groups[c]:
        natural = (row[0] == 'rational' and row[1][0] == row[2][0] == t
                   or row[0] == 'quadratic' and row[1] == t)
        if not natural:
            exceptional.append([c, (t, e), row])
result = {
    'status': 'PREVIEW_ONLY',
    'bindings': {str((OLD/f).relative_to(ROOT)): hashlib.sha256((OLD/f).read_bytes()).hexdigest()
                 for f in ['input.json', 'result.json', 'rational-fibres.json', 'quadratic-fibres.json']},
    'rational_root_codes': len(codes),
    'quadratic_zero_rows': [r for r in quadratic if 0 in r['norm_codes']],
    'norm_13412_rows': groups[13412],
    'exceptional_single_code_matches': exceptional,
    'coefficient_solution': False,
}
out = Path(__file__).resolve().parent.parent / 'character-preview.json'
with out.open('x') as f:
    json.dump(result, f, indent=2, sort_keys=True)
    f.write('\n')
print(json.dumps({'rational_root_codes': len(codes), 'quadratic_zero_rows': len(result['quadratic_zero_rows']),
                  'norm_13412_rows': result['norm_13412_rows'], 'exceptional_matches': len(exceptional)}))
