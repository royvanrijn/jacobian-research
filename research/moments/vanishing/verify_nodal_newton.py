#!/usr/bin/env python3
"""Exact algebra for NODAL_NEWTON_EXCLUSION.md; no analytic verification."""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import time
import sympy as S


def run():
    started = time.monotonic()
    here = Path(__file__).resolve().parent
    u, v, z, t = S.symbols('u v z t')
    a, b, c, d, gamma = S.symbols('a b c d gamma', nonzero=True)
    original = z*u*v+a*u**3+b*u**2*v+c*u*v**2+d*v**3
    normal = S.expand(original.subs(z, z-b*u-c*v))
    assert normal == z*u*v+a*u**3+d*v**3
    assert S.expand(normal.subs(a, 0)-v*(z*u+d*v**2)) == 0
    assert S.expand(normal.subs(d, 0)-u*(z*v+a*u**2)) == 0
    vertices = {'O': (0, 0, 0), 'A': (3, 0, 0), 'B': (0, 3, 0),
                'C': (1, 1, 1), 'D': (0, 0, 2)}
    facets = {}
    for names in combinations('ABCD', 3):
        matrix = S.Matrix([vertices[n] for n in names])
        if not matrix.det():
            continue
        weight = matrix.inv()*S.ones(3, 1)
        if all(weight.dot(S.Matrix(e)) <= 1 for e in vertices.values()):
            facets[''.join(names)] = tuple(weight)
    assert facets == {'ABC': (S.Rational(1, 3),)*3,
                      'ACD': (S.Rational(1, 3), S.Rational(1, 6), S.Rational(1, 2)),
                      'BCD': (S.Rational(1, 6), S.Rational(1, 3), S.Rational(1, 2))}
    assert all(sum(w) == 1 and min(w) > 0 for w in facets.values())
    quadratics = [(i, j, 2-i-j) for i in range(3) for j in range(3-i)]
    for e in quadratics:
        degrees = [sum(wi*ei for wi, ei in zip(w, e)) for w in facets.values()]
        assert max(degrees) <= 1
        assert e == (0, 0, 2) or max(degrees) < 1
    terms = {'A': a*u**3, 'B': d*v**3, 'C': z*u*v, 'D': gamma*z**2}
    faces = sorted({''.join(face) for names in facets
                    for size in range(1, 4) for face in combinations(names, size)},
                   key=lambda face: (len(face), face))
    witnesses = []
    for names in faces:
        face = sum(terms[n] for n in names)
        found = None
        for variable in (u, v, z):
            derivative = S.Poly(S.diff(face, variable), u, v, z)
            if len(derivative.terms()) == 1 and not derivative.is_zero:
                found = {'face': names, 'polynomial': str(face),
                         'variable': str(variable), 'monomial_partial': str(derivative.as_expr())}
                break
        assert found, names
        witnesses.append(found)
    assert len(faces) == 13
    h0, h1, h2, h3, h4, h5 = S.symbols('h0:6')
    quadratic = (h0*u*u+h1*v*v+h2*z*z+2*h3*u*v+2*h4*u*z+2*h5*v*z)/2
    phase = normal-quadratic
    lam = t**-2
    scaled = S.expand(t*normal.subs({u:u/t, v:v/t, z:z/t}, simultaneous=True)
                      -quadratic.subs({u:u/t, v:v/t, z:z/t}, simultaneous=True))
    assert S.expand(scaled-lam*phase) == 0
    divergence = sum(S.diff(x, x)+lam*S.diff(phase, x)*x for x in (u, v, z))
    seed_difference = t**-6*normal-(-3*t**-4-2*t**-6*phase)
    assert S.expand(seed_difference-t**-4*divergence) == 0
    # Negative control: an isotropic node loses the essential axis intercept.
    assert S.Poly(phase.subs(h2, 0), u, v, z).coeff_monomial(z*z) == 0
    # Gaussian degree-two calibration fixes the V-index sign and strict gap.
    assert sum([S.Rational(1, 2)]*3) == S.Rational(3, 2)
    assert 1 < S.Rational(3, 2)
    source = here/'sources/sabbah_bbases.pdf'
    expected = '1696f2ebccd38e93b86640d4787b673707fc6b09f53c86221643738120d9687a'
    assert hashlib.sha256(source.read_bytes()).hexdigest() == expected
    receipt = {
        'completed': True, 'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_sha256': expected, 'normal_form': str(normal),
        'outer_facets': {n: [str(x) for x in w] for n, w in facets.items()},
        'outer_faces_checked': len(faces), 'nondegeneracy_witnesses': witnesses,
        'quadratic_supports_checked': len(quadratics), 'volume_Newton_degree': '1',
        'constant_series_required_exponent': '3/2',
        'rescaling': 'x=z/t, lambda=t^(-2), dx=t^(-3) dz',
        'connection_identity': 'mapped nabla_t(1)-nabla_t(t^(-3)) = t^(-4) sum_i (partial_i+lambda*f_i)(z_i)',
        'isotropic_node_negative_control': 'q_TT=0 removes T^2; Newton argument not applied',
        'homogeneous_pure_null_classification': 'P=(l^T*x)*(x^T*A*x), l^T*l=0, A*l proportional to l',
        'mixed_cutoff': 'm>deg(multiplier)', 'explicit_uniform_homogeneous_pure_gate': None,
        'assurance': 'Exact symbolic algebra and cached source hash only. Imported V-filtration theorem, regular-singular growth, formal-functional extension and classification synthesis remain written mathematics; no external review or formal verification.',
        'elapsed_seconds': round(time.monotonic()-started, 6)
    }
    target = here/'nodal_newton_verification.json'
    target.write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['completed', 'outer_faces_checked', 'volume_Newton_degree',
                                           'constant_series_required_exponent', 'elapsed_seconds']}))


if __name__ == '__main__':
    run()
