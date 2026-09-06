#!/usr/bin/env python3
"""Complete the universal bad13 claim by checking the original residue cells."""
import argparse
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
import report_r17_13_scaling_geometry as previous
from research_runtime.store import checkpoint

ROOT = previous.ROOT
OUT = previous.ART / 'r17_13_scaling_geometry_v2.json'


def expected():
    prior = previous.expected()
    if cert.read(previous.OUT) != json.loads(json.dumps(prior)):
        raise ArithmeticError('retained v1 geometry proof differs')
    rows = []
    families = cert.read(previous.charts.INPUT)['families']
    if len(families) != 6 or len(prior['rows']) != 6:
        raise ArithmeticError('all six families required')
    for family, old in zip(families, prior['rows']):
        if family['family'] != old['family']:
            raise ArithmeticError('family order differs')
        a = list(map(int, family['A_coefficients_low_to_high']))
        b = list(map(int, family['B_coefficients_low_to_high']))
        if len(a) != 9 or len(b) != 13:
            raise ArithmeticError('homogeneous degrees8 and12 required')
        cells = []
        for residue in [*range(13), 'infinity']:
            n, d = (1, 0) if residue == 'infinity' else (residue, 1)
            A = sum(c * n**i * d**(8-i) for i, c in enumerate(a)) % 13
            B = sum(c * n**i * d**(12-i) for i, c in enumerate(b)) % 13
            delta = (-16 * (4*A**3 + 27*B**2)) % 13
            if delta != 0:
                raise ArithmeticError('original model has a good13 residue')
            cells.append({'residue': residue, 'A_mod13': A,
                          'B_mod13': B, 'discriminant_mod13': delta})
        if not old['scaled_models_still_bad_at13']:
            raise ArithmeticError('divided model bad reduction required')
        rows.append({**old, 'original_model_residues': cells,
                     'every_nonsingular_primitive_fibre_has_bad_reduction_at13': True})
    paths = [Path(__file__).resolve(), Path(previous.__file__).resolve(),
             previous.OUT, previous.charts.INPUT]
    return {'schema': 'elliptic-curves.r17-13-scaling-geometry.v2',
            'status': 'PASS',
            'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths},
            'retained_v1_certificate': str(previous.OUT.relative_to(ROOT)),
            'rows': rows, 'original_projective_discriminant_cells_checked': 84,
            'claim_boundary': prior['claim_boundary'] +
            ' Additionally, all84 original projective discriminant residues are checked directly. Outside the unique nonminimal cell, v1 proves the original model is13-minimal; inside it, v1 proves the once-divided model is13-minimal and bad. Thus every nonsingular fibre at a primitive rational parameter has bad reduction at13. This is a proof completion, not a new parameter or point search.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = expected()
    if args.check:
        if cert.read(OUT) != json.loads(json.dumps(result)):
            raise ArithmeticError('complete bad13 certificate differs')
    else:
        if OUT.exists():
            raise FileExistsError('preserve v2 complete bad13 certificate')
        checkpoint(OUT, result)
    print('ALL84 ORIGINAL RESIDUES AND SIX MINIMAL13 CLASSIFICATIONS PASS', flush=True)
