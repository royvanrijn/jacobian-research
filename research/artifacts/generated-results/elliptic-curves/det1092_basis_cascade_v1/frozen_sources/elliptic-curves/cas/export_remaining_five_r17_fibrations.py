#!/usr/bin/env python3
"""Standalone Sage equations and all generic sections for the five R17 search models."""
import argparse
import json
from pathlib import Path
import certify_compact_r17_candidates as cert

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
ATLAS = ART / 'compact_six_r17_atlas_v1.json'
OUT = ART / 'remaining_five_r17_fibrations.sage'
MANIFEST = ART / 'remaining_five_r17_fibrations_export_v1.json'


def expected():
    atlas = cert.read(ATLAS)
    for name, digest in atlas['checker_sources'].items():
        if cert.hashed(ROOT / name) != digest:
            raise ArithmeticError('original generic atlas checker changed')
    keep = ['103b2', '074d9', '07ca9', '08234', '08f72']
    families = []
    for f in atlas['families']:
        if f['family'] not in keep:
            continue
        if f['generic_rank_lower_bound'] != 17 or len(f['sections']) != 17:
            raise ArithmeticError('exact generic17 atlas required')
        families.append({k: f[k] for k in ['family', 'A_coefficients_low_to_high',
            'B_coefficients_low_to_high', 'sections', 'generic_height_gram',
            'generic_height_gram_determinant', 'base_matrix_a_b_c_d',
            'constant_scaling', 'parameter_convention']})
    if [f['family'] for f in families] != keep:
        raise ArithmeticError('complete remaining-five atlas order required')
    data = {'atlas_sha256': cert.hashed(ATLAS), 'families': families}
    script = '''# Existing compact R17 fibrations: exact equations and generic sections.
# These models are not newly discovered fibrations or new specialized rank results.
# Generic independence is inherited from the exact compact-atlas transport proof.
# This file independently checks all85 rational-function point identities and
# the supplied positive definite generic height forms; it does not recompute heights.
import json
from sage.all import QQ, FunctionField, EllipticCurve, matrix

DATA = json.loads(r\'''PAYLOAD\''')
K = FunctionField(QQ, 't')
t = K.gen()
curves = {}
generic_points = {}
height_grams = {}

def polynomial(coefficients, z):
    value = 0
    for c in reversed(coefficients):
        value = value*z + QQ(c)
    return value

def rational(record, z):
    return polynomial(record['numerator_coefficients_low_to_high'], z) / polynomial(record['denominator_coefficients_low_to_high'], z)

for f in DATA['families']:
    name = f['family']
    A = polynomial(f['A_coefficients_low_to_high'], t)
    B = polynomial(f['B_coefficients_low_to_high'], t)
    E = EllipticCurve(K, [A, B])
    points = [E(rational(P['X'], t), rational(P['Y'], t)) for P in f['sections']]
    G = matrix(QQ, f['generic_height_gram'])
    assert len(points) == 17 and G.nrows() == G.ncols() == 17
    assert G == G.transpose() and G.det() == QQ(f['generic_height_gram_determinant']) == 948
    assert all(G[:i,:i].det() > 0 for i in range(1,18))
    curves[name], generic_points[name], height_grams[name] = E, points, G

def fibre(family, parameter):
    """Return an integral Q fibre and its transported seventeen generic points.

    Singular fibres or poles of the generic sections raise an exception.
    Specialized independence must be checked separately before claiming rank17.
    """
    f = next(r for r in DATA['families'] if r['family'] == family)
    z = QQ(parameter)
    d = z.denominator()
    A = polynomial(f['A_coefficients_low_to_high'], z)*d**8
    B = polynomial(f['B_coefficients_low_to_high'], z)*d**12
    assert A.denominator() == B.denominator() == 1
    E = EllipticCurve(QQ, [A, B])
    points = [E(rational(P['X'], z)*d**4, rational(P['Y'], z)*d**6) for P in f['sections']]
    return E, points

assert len(curves) == 5 and sum(map(len,generic_points.values())) == 85
print('FIVE R17 FIBRATIONS:85 EXACT GENERIC POINT IDENTITIES PASS')
'''.replace('PAYLOAD', json.dumps(data, sort_keys=True, separators=(',', ':')))
    return script, {'schema': 'elliptic-curves.remaining-five-r17-export.v1',
        'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in [Path(__file__).resolve(), ATLAS]},
        'families': keep, 'generic_sections': 85,
        'claim_boundary': 'Standalone equations, generic points and transported generic height forms for five existing compact R17 fibrations. Extraction preserves the exact atlas values. Running the Sage file checks all rational-function point identities and positive definiteness of the supplied height forms; the original atlas proof supplies their geometric provenance. No new fibration, automatic rank17 claim on every specialization, exact rank, new point search or record.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--check', action='store_true'); args = p.parse_args()
    script, data = expected()
    if args.check:
        if OUT.read_text() != script:
            raise ArithmeticError('standalone generic export differs')
    else:
        if OUT.exists() or MANIFEST.exists():
            raise FileExistsError('preserve generic export')
        OUT.write_text(script)
    data['output'] = str(OUT.relative_to(ROOT)); data['output_sha256'] = cert.hashed(OUT)
    if args.check:
        if cert.read(MANIFEST) != data:
            raise ArithmeticError('generic export source bindings differ')
    else:
        cert.write(MANIFEST, data)
    print('EXACT FIVE R17 EXPORT EXTRACTION PASS', flush=True)
