#!/usr/bin/env sage -python
"""Exact generic height geometry of the nine retained quadratic covers.

No specialized numerical heights or public complementary points are used.
The old Gram and trace labels are bound to the specialization basis. Polynomial
identities check the section, its trace, and the hypotheses for h(T)=8.
"""
import argparse
import csv
import json
from hashlib import sha256
from pathlib import Path
import sys

from sage.all import QQ, ZZ, PolynomialRing, matrix, vector, pari

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'elliptic-curves/cas'))
from research_runtime.store import checkpoint

COVERS = ROOT / 'artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-covers-v1.json'
PRIORITY = ROOT / 'artifacts/generated-results/elkies-2026-bisection-equation-priority-full.tsv'
PUBLISHED = ROOT / 'elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json'


def read(path, inputs):
    inputs[str(path.relative_to(ROOT))] = sha256(path.read_bytes()).hexdigest()
    return json.loads(path.read_text())


def build(progress=None):
    inputs = {}
    covers = read(COVERS, inputs)
    R = PolynomialRing(QQ, 't'); t = R.gen(); K = R.fraction_field()
    def poly(xs): return R(list(map(QQ, xs)))
    def rat(row): return K(poly(row['numerator_coefficients_low_to_high']) / poly(row['denominator_coefficients_low_to_high']))
    roster = [(c, ch['direct_model']) for ch in covers['charts'] for f in ch['fibres'] for c in f['covers']]
    h = covers['historical_rank28_anchor']; roster.append((h, h['direct_model']))
    # Parse the published Gram constant without executing the historical search.
    import ast
    old_source = ROOT / 'elliptic-curves/cas/half_lattice_fake_descent_replay.sage'
    inputs[str(old_source.relative_to(ROOT))] = sha256(old_source.read_bytes()).hexdigest()
    tree = ast.parse(old_source.read_text())
    published_gram = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
                          and any(isinstance(x, ast.Name) and x.id == 'GENERIC_GRAM' for x in n.targets))
    inputs[str(PRIORITY.relative_to(ROOT))] = sha256(PRIORITY.read_bytes()).hexdigest()
    with PRIORITY.open() as stream:
        historical_word = next(list(map(int, r['published_basis_w'].split()))
                               for r in csv.DictReader(stream, delimiter='\t') if int(r['orbit_mask'], 0) == 88680)
    records = []
    for cover, model_path in roster:
        model = read(ROOT / model_path, inputs)
        equation = model.get('weierstrass_model', model)
        A, B = (poly(equation[k + '_coefficients_low_to_high']) for k in ('A', 'B'))
        delta = 4*A**3 + 27*B**2
        q = poly(cover['branch_quadratic_coefficients_low_to_high'])
        # Simple discriminant roots imply I1 fibres; infinity is smooth.
        if (A.degree(), B.degree(), delta.degree()) != (8, 12, 24):
            raise ArithmeticError('unexpected Weierstrass degrees')
        if delta.gcd(delta.derivative()).degree() != 0 or q.degree() != 2 or q.discriminant() == 0 or q.gcd(delta).degree() != 0:
            raise ArithmeticError('branch or reducible-fibre correction hypotheses failed')
        section = cover['eighteenth_section']
        x0, x1, y0, y1 = (poly(section[k + '_coefficients_low_to_high']) for k in ('x0', 'x1', 'y0', 'y1'))
        if any(p.degree() > bound for p, bound in zip((x0, x1, y0, y1), (4, 3, 6, 5))):
            raise ArithmeticError('section may intersect zero at infinity')
        if y0*y0 + q*y1*y1 != x0**3 + 3*x0*x1*x1*q + A*x0 + B or 2*y0*y1 != 3*x0*x0*x1 + q*x1**3 + A*x1:
            raise ArithmeticError('quadratic section equation failed')
        # Chord joining T and sigma(T): slope y1/x1; trace is rational.
        slope = K(y1/x1)
        X = slope*slope - 2*x0
        Y = -y0 + slope*(x0-X)
        if Y*Y != X**3 + A*X + B:
            raise ArithmeticError('trace equation failed')
        if 'sections' in model:
            G = matrix(QQ, model['sections']['height_gram'])
            word = vector(ZZ, cover['equation_basis_trace_word'])
            generic = [(rat(r['X']), rat(r['Y'])) for r in model['sections']['records']]
        else:
            G = matrix(QQ, published_gram); word = vector(ZZ, historical_word)
            source = read(PUBLISHED, inputs); generic = []
            for r in source['sections']:
                x = poly(r['x_coefficients_low_to_high'])
                if not generic: y = poly(r['y_coefficients_low_to_high'])
                else:
                    chord = r['chord']; ref = generic[chord['reference_basis_index']]
                    y = ref[1] + poly(chord['slope_coefficients_low_to_high'])*(x-ref[0])
                generic.append((K(x), K(y)))
        # Bind the trace to its exact basis word by function-field group law.
        from sage.all import EllipticCurve
        E = EllipticCurve(K, [A, B])
        # Add signed basis sections in small-height order. Computing 10*P
        # first creates huge rational functions which later cancel entirely.
        remaining = list(map(int, word)); accumulated = vector(ZZ, [0]*17)
        Q = E(0); max_intermediate_height = 0
        while any(remaining):
            choices = []
            for i, w in enumerate(remaining):
                if not w: continue
                step = 1 if w > 0 else -1
                proposed = vector(ZZ, list(accumulated)); proposed[i] += step
                choices.append((proposed*G*proposed, i, step, proposed))
            height, i, step, accumulated = min(choices, key=lambda r:(r[0],r[1]))
            Q += step*E(generic[i]); remaining[i] -= step
            max_intermediate_height = max(max_intermediate_height, int(height))
        if Q.is_zero() or Q[0] != X or Q[1] != Y or word*G*word != 10 or G.det() != 948:
            raise ArithmeticError('trace basis identity or inherited Gram failed')
        cross = G*word
        G18 = (2*G).augment(matrix(QQ, 17, 1, list(cross)))
        G18 = G18.stack(matrix(QQ, 1, 18, list(cross)+[8]))
        if not G18.is_positive_definite() or G18.det() != 2**17*948*3:
            raise ArithmeticError('height Gram or Schur complement failed')
        U = matrix(ZZ, pari(G18).qflllgram()).transpose()
        reduced = U*G18*U.transpose()
        if abs(U.det()) != 1: raise ArithmeticError('non-unimodular reduction')
        records.append({'cover_label': cover['label'], 'model_source': model_path,
            'old_gram': [[str(v) for v in r] for r in G.rows()],
            'trace_word': list(map(int, word)), 'trace_height_over_old_base': 10,
            'height_gram': [[str(v) for v in r] for r in G18.rows()],
            'section_height': 8, 'orthogonal_residual_height': 3, 'determinant': str(G18.det()),
            'intersection_checks': {'old_discriminant_degree': 24, 'old_discriminant_squarefree': True,
                'branch_disjoint_from_discriminant': True, 'branch_degree': 2,
                'section_polynomial_degrees': list(map(lambda p: int(p.degree()), (x0, x1, y0, y1))),
                'zero_section_intersection': 0, 'reducible_fibre_correction': 0,
                'chi_after_base_change': 4, 'section_and_trace_word_identities_exact': True},
            'reduced_basis_rows_in_original_basis': [list(map(int, r)) for r in U.rows()],
            'reduced_gram': [list(map(int, r)) for r in reduced.rows()]})
        print('MW18_GEOMETRY|' + cover['label'] + '|height=8|schur=3', flush=True)
        if progress:
            checkpoint(progress, {'status':'PARTIAL_CHECKPOINT','inputs':inputs,'covers':records})
    return {'schema': 'elliptic-curves.mw18-height-geometry.v1', 'status': 'PASS', 'inputs': inputs,
        'sources': {str(Path(__file__).resolve().relative_to(ROOT)): sha256(Path(__file__).read_bytes()).hexdigest()},
        'height_formula_reference': 'https://arxiv.org/abs/0907.0298, sections 11.8, 11.19 and 11.20',
        'basis': 'ordered specialized generic_R17, followed by cover_section', 'covers': records,
        'claim_boundary': 'Exact height Gram of the displayed generic rank-18 subgroup; no saturation or exact generic-rank claim.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    checkpoint(args.output, build(args.output.with_suffix('.partial.json')))
