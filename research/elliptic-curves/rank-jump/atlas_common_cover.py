#!/usr/bin/env python3
"""Group a frozen generic atlas by exact quadratic function field."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import reduce
import gzip
import json
from math import gcd
from pathlib import Path
import retrospective as r
from cover_experiment import sqrtq

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / 'ATLAS_COMMON_COVER_PROTOCOL.json'
GEOMETRY = r.OUT / 'rank_jump_solubility_first_geometry_v1.json.gz'
OUTPUT = r.OUT / 'rank_jump_atlas_common_cover_v1.json'
GEOMETRY_SHA = '3ada2afa96b4f63944b0fdeebc2a28910870ba0391eee29a8721588426cda450'


def analyse():
    raw = GEOMETRY.read_bytes(); assert r.digest(raw) == GEOMETRY_SHA
    covers = json.loads(gzip.decompress(raw)); assert len(covers) == 39119
    groups = defaultdict(list); degenerate = []
    normalized = []
    for cover in covers:
        f = cover['integer_quadratic']; label = cover['label']
        if len(f) != 3 or not f[2] or f[1] ** 2 == 4 * f[0] * f[2]:
            degenerate.append(label); continue
        scale = reduce(gcd, f)
        if f[2] < 0: scale = -scale
        primitive = tuple(x // scale for x in f)
        assert reduce(gcd, primitive) == 1 and primitive[-1] > 0
        row = {'label': label, 'scale': str(F(scale, cover['denominator_scale'] ** 2))}
        groups[primitive].append(row)
        normalized.append({'label': label, 'primitive': list(primitive), 'scale': row['scale']})
    same_branch = []; same_field = []; field_count = 0
    for polynomial, rows in sorted(groups.items()):
        classes = []
        for row in rows:
            matched = False
            for cls in classes:
                root = sqrtq(F(row['scale']) / F(cls[0]['scale']))
                if root is not None:
                    cls.append(row | {'square_ratio_root': str(root)}); matched = True; break
            if not matched: classes.append([row | {'square_ratio_root': '1'}])
        field_count += len(classes)
        if len(rows) > 1:
            same_branch.append({'primitive': list(polynomial), 'members': rows, 'field_classes': classes})
        for cls in classes:
            if len(cls) > 1: same_field.append({'primitive': list(polynomial), 'members': cls})
    assert not degenerate, 'Nonsquarefree rows require a separate squareclass audit'
    return {'schema': 'rank-jump.atlas-common-cover.v1', 'status': 'PASS', 'layer': 'incidence',
            'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes())
                         for p in (PROTOCOL, GEOMETRY, Path(__file__), HERE / 'retrospective.py', HERE / 'cover_experiment.py')},
            'equations': len(covers), 'distinct_branch_polynomials': len(groups),
            'distinct_quadratic_function_fields': field_count,
            'degenerate_labels': degenerate, 'same_branch_groups': same_branch, 'same_field_groups': same_field,
            'branch_group_size_counts': dict(sorted(Counter(str(len(v)) for v in groups.values()).items())),
            'normalized_roster_sha256': r.digest(json.dumps(normalized, sort_keys=True, separators=(',', ':')).encode()),
            'boundary': r.read(PROTOCOL)['boundary']}


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('mode', choices=['build', 'check']); mode = p.parse_args().mode
    result = analyse()
    if mode == 'build': r.write_new(OUTPUT, result)
    else: assert r.read(OUTPUT) == result
    print('PASS', result['equations'], 'equations;', result['distinct_branch_polynomials'], 'branch polynomials;',
          result['distinct_quadratic_function_fields'], 'quadratic fields;', len(result['same_field_groups']), 'repeated-field groups')
