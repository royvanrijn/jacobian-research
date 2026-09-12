#!/usr/bin/env python3
"""Freeze the existing point-free covers and reconstruct only generic points.

This preparation never opens a point-search or rational-lift certificate.
Run from any directory; all outputs must be new.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'artifacts/generated-results/elliptic-curves'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def square_root(x):
    x = F(x)
    a, b = isqrt(x.numerator), isqrt(x.denominator)
    assert a*a == x.numerator and b*b == x.denominator
    return F(a, b)


def write_new(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, sort_keys=True, indent=2)
        stream.write('\n')


def freeze(folder):
    inputs = [SOURCE/'rank_jump_reference_strict_class_construction_inputs_v1.json',
              SOURCE/'rank_jump_constructed_class_compaction_v1.json']
    ref, compact = [json.loads(p.read_text()) for p in inputs]
    f = list(map(F, ref['cubic_ascending']))
    assert f[2:] == [1, 1] and len(ref['generic_classes']) == 16
    assert compact['status'] == 'PASS' and [c['column'] for c in compact['cases']] == [6, 7]
    model = [F(1), F(0), F(0), f[1]/16, f[0]/64]
    generic = []
    for row in ref['generic_classes']:
        b = list(map(F, row['beta_ascending']))
        assert b[2] == 0 and b[1] < 0
        scale = square_root(-b[1])
        X = b[0]/(scale*scale)
        Y = square_root(sum(a*X**i for i, a in enumerate(f)))
        x, y = X/4, (Y-X)/8
        assert y*y+x*y == x**3+model[3]*x+model[4]
        assert F(row['norm']) == scale**6 * Y*Y
        generic.append([str(x), str(y)])
    # The integral short model is x_short=36*x+3, y_short=216*y+108*x.
    A, B = 1296*model[3]-27, 46656*model[4]-3888*model[3]+54
    short_points = []
    for x, y in generic:
        X, Y = 36*F(x)+3, 216*F(y)+108*F(x)
        assert Y*Y == X**3+A*X+B
        short_points.append([str(X), str(Y)])
    binding = {str(p.relative_to(ROOT)): digest(p) for p in inputs}
    binding[str(Path(__file__).resolve().relative_to(ROOT))] = digest(Path(__file__))
    covers = {'schema': 'constructed-class-blind-covers.v1',
        'cubic_ascending': ref['cubic_ascending'],
        'original_curve': list(map(str, model)),
        'coordinates': ['u', 'v', 'w', 's'],
        'cover_definition': 'beta*(u+v*theta+w*theta^2)^2=Q0+Q1*theta+Q2*theta^2; Q2=0, Q1+s^2=0',
        'map': 'X=Q0/s^2; Y=sqrt(Norm(beta))*Norm(u+v*theta+w*theta^2)/s^3; x=X/4, y=(Y-X)/8',
        'cases': [{k: c[k] for k in ['column', 'beta_ascending', 'norm', 'cover_quadrics', 'x_numerator', 'max_coefficient_bits']}
                  for c in compact['cases']],
        'boundary': 'Fixed class inputs only. No cover point, field square root, exceptional point, or point-combination input.'}
    for row in covers['cases']:
        row['positive_norm_square_root'] = str(square_root(row['norm']))
    subgroup = {'schema': 'constructed-class-generic-subgroup.v1',
        'original_curve': list(map(str, model)), 'original_points': generic,
        'curve': list(map(str, [0, 0, 0, A, B])), 'points': short_points,
        'generic_rank': 16,
        'transport': 'x_short=36*x+3; y_short=216*y+108*x',
        'provenance': 'Each x is reconstructed from its generic Kummer representative; y uses the positive cubic ordinate. Signs may differ from historical section ordering, but the generated subgroup is identical.'}
    folder.mkdir(parents=True, exist_ok=False)
    write_new(folder/'covers.json', covers)
    write_new(folder/'generic.json', subgroup)
    protocol = {'schema': 'constructed-class-blind-experiment.v1',
        'question': 'Recover points on two previously fixed extra classes without point-oracle inputs; compare ordinary V3 initialized only with the same generic subgroup.',
        'selection': 'User-fixed MW16-05 at t=3/17 and compaction columns 6 and 7; retrospective positive control, not prospective fibre selection.',
        'input_bindings': binding,
        'frozen_inputs': {n: digest(folder/n) for n in ['covers.json', 'generic.json']},
        'limits': {'computation_seconds_per_arm': 1200, 'rss_bytes_per_worker': 2*1024**3,
                   'maximum_parallel_workers_per_arm': 1, 'maximum_single_preparation_step_seconds': 60,
                   'v3_maximum_charts': 100, 'v3_chart_height': 125000,
                   'v3_seconds_per_chart': 10, 'v3_seconds_per_map': 5},
        'blindness': 'Separate fresh-context cover worker reads only covers.json, permitted algorithm sources/docs and its own generated files. Root handles generic-only V3 and independent rational verification. Shared filesystem separation is an instruction/dataflow boundary, not OS access isolation.',
        'cover_preparation': 'Use established genus-one minimisation/reduction where available; retain all transformations and coefficientwise identities. A bounded miss or unavailable backend cannot establish point absence.',
        'comparison': 'Report setup, reduction, search, verification and total costs separately; count recovered fixed covers separately from new V3 subgroup directions. No timing superiority claim from unequal or incomplete exposure.',
        'success_gate': 'Primitive rational cover coordinates, exact transport to the original elliptic equation, beta*xi^2=4*x-theta, and independent exact-arithmetic replay; no oracle-derived discovery input.',
        'forbidden': ['known exceptional points', 'oracle square roots in the cubic field', 'retrospective point combinations', 'old point-search output as execution input'],
        'prior_work': ['elliptic-curves/rank-jump/CONSTRUCTED_CLASS_BLOCK_AND_RATIONAL_LIFTS.md',
                       'elliptic-curves/rank-jump/TWO_CONSTRUCTED_STRICT_CLASSES_AND_302.md',
                       'elliptic-curves/notes/TWO_CLASS_RELATION_PILOT_2026-09-12.md'],
        'scope': 'One bounded paired experiment. No relation recollection, full class group, Selmer upper bound, parameter sweep, new curve, or rank upper bound.'}
    write_new(folder/'protocol.json', protocol)
    print(json.dumps({'output': str(folder), 'frozen_inputs': protocol['frozen_inputs']}))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    freeze(p.parse_args().output.resolve())
