#!/usr/bin/env python3
"""Verify the full finite-cover degree calculation and length arithmetic.

Uses the retained exact msolve length9168 as a dependency, without rerunning
that CAS. No global A12 finiteness, actual d14 or full Image result is asserted.
"""
import hashlib
import json
from pathlib import Path

from anti_involution_cover import (
    S, a, c, lam, k, d, D0, FREE_CUBIC, ODD_TRIPLES,
    normalized_locus, chart_base, inverse_base, cubic_change, z, t, w,
)
from cyclic_frame_quotient import G, C, p, r, s0, s1, s2, U, normalized_relations
from verify_anti_involution_prefix import PRIME, STEM, leading_staircase


def zero(value):
    assert S.cancel(value) == 0, value


def main():
    root = Path(__file__).resolve().parent
    assert PRIME % 2 and PRIME % 3
    special = normalized_locus()
    for relation in normalized_relations():
        zero(relation.subs(special, simultaneous=True))
    gram = G.subs(special, simultaneous=True)
    companion = C.subs(special, simultaneous=True)
    rotation = S.diag(-1, 1, -1)
    assert rotation.det() == 1
    assert rotation.T*gram*rotation == gram
    assert rotation*companion*rotation == -companion
    for indices, coefficient in U.items():
        sign = S.prod(rotation[i, i] for i in indices)
        if indices in ODD_TRIPLES:
            assert sign == 1 and special[coefficient] == 0
        else:
            assert sign == -1

    # The determinant forces c to be a unit before solving for p.
    restricted_det = S.expand(G.det().subs({r: 0, s1: 0, s0: a, s2: c}))
    assert restricted_det == a*p*c*c-c**3
    zero(c*(c*c-a*p*c)-D0**2-( -restricted_det-D0**2))
    # c * ((c^2-a*p*c)/D0^2) = 1 modulo det(G)+D0^2.
    zero(c*(c*c-a*p*c)/D0**2-1+ (restricted_det+D0**2)/D0**2)

    source, inverse = chart_base(), inverse_base()
    for expression, expected in [(inverse['difference'], d),
                                 (inverse['bz'], source['bz']),
                                 (inverse['bw'], source['bw'])]:
        pulled = expression.subs({a: source[a], c: source[c]}, simultaneous=True)
        # The only remaining equation on the common cover fixes its orientation.
        zero(pulled.subs(d, 2*D0/(lam**6*k*k))-expected.subs(d, 2*D0/(lam**6*k*k)))
    zero((source['normalized_D']-D0).subs(d, 2*D0/(lam**6*k*k)))
    zero((inverse['eta']*2*k*inverse['difference'])-1)
    zero((inverse['k_square']-k*k).subs({a: source[a], c: source[c]}, simultaneous=True))
    zero((inverse['lambda_square']-lam*lam).subs(a, source[a]))
    zero((special[p].subs({a: source[a], c: source[c]}, simultaneous=True)-source[p])
         .subs(d, 2*D0/(lam**6*k*k)))

    # Both quadratic root extensions are monic with unit constant and derivative.
    # The amplitude root is likewise monic of degree6 over the original chart.
    X = S.Symbol('X')
    for degree, unit in [(6, D0/source['chart_D']), (2, a/2), (2, 4*c/a**2)]:
        f = X**degree-unit
        assert S.Poly(f, X).LC() == 1
        assert S.diff(f, X) == degree*X**(degree-1)
        assert degree % PRIME

    transform, original, y = cubic_change()
    determinant = S.factor(transform.det())
    assert determinant == S.Rational(9, 16)*d**6*k**6*lam**20
    assert transform.shape == (4, 4)
    reconstructed = (transform*S.Matrix(FREE_CUBIC))
    az, aw, ez, ew = reconstructed
    zero(original-(az*z**3+aw*w**3+ez*z*(t*t-z*w/2)+ew*w*(t*t-z*w/2)))
    zero(2*S.diff(original, z, w)+S.diff(original, t, 2))
    metric = S.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    v = S.Matrix([1, 0, 1])
    bz, bw = source['bz'], source['bw']
    operator = S.Matrix([[0, bw/2, 0], [bz/2, 0, bw/2], [0, bz/2, 0]])
    frame = S.Matrix.hstack(lam*v, lam**2*operator*v, lam**3*operator**2*v)
    physical = S.Matrix([z, t, w])
    assert (frame.T*metric*physical-S.Matrix(y)).applyfunc(S.expand).is_zero_matrix
    zero(frame.det()-source['normalized_D'])
    # Exact phase identity proves moment ideals pull back with factors lam^n
    # at every order, including nilpotent coefficient specializations.
    quad = (physical.T*metric*(lam*operator)*physical)[0]
    zero(quad-lam*t*(bz*z+bw*w))
    cubic_norm = lam*original
    phase_norm = y[0]+quad+cubic_norm
    zero(phase_norm-lam*(z+w+t*(bz*z+bw*w)+original))

    receipt_path = root/'anti_involution_prefix12_verification.json'
    receipt = json.loads(receipt_path.read_text())
    assert receipt['status'] == 'PASS' and receipt['prime'] == PRIME
    assert receipt['checker_sha256'] == hashlib.sha256(
        (root/'verify_anti_involution_prefix.py').read_bytes()).hexdigest()
    for name, digest in receipt['hashes'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest() == digest
    basis_size, bounds, count = leading_staircase(root/(STEM+'.out'))
    assert count == receipt['modular_length'] == 9168
    amplitude_degree, quotient_degree = 6, 4
    extended_length = amplitude_degree*count
    assert extended_length % quotient_degree == 0
    normalized_length = extended_length//quotient_degree
    assert normalized_length == 13752 and normalized_length+3*22 == 13818

    report = {
        'claim': 'SIC2B33CYCLICLENGTH20261116', 'status': 'PASS',
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'model_sha256': hashlib.sha256((root/'anti_involution_cover.py').read_bytes()).hexdigest(),
        'parent_receipt_sha256': hashlib.sha256(receipt_path.read_bytes()).hexdigest(),
        'prime': PRIME, 'input_modular_chart_length': count,
        'source_amplitude_cover_degree': amplitude_degree,
        'normalized_quotient_cover_degree': quotient_degree,
        'common_cover_modular_length': extended_length,
        'normalized_anti_halfturn_modular_length': normalized_length,
        'cubic_change_determinant': str(determinant),
        'cubic_change_matrix': [[str(value) for value in row] for row in transform.tolist()],
        'normalized_locus_substitution': {str(key): str(value) for key, value in special.items()},
        'inverse_base': {key: str(value) for key, value in inverse.items()},
        'DVR_flatness': 'Regular-sequence flatness at all special closed points '
            'kills Ann(ell) globally. Torsion-free basis lifting preserves all9168 '
            'special classes as Q-linearly independent. The rank6 common cover '
            'descends flatness through its rank4 cover of the normalized locus.',
        'normalized_characteristic_zero_dimension_lower_bound': normalized_length,
        'A12_dimension_lower_bound': normalized_length+3*22,
        'A12_finiteness_and_exact_length': 'UNCOMPUTED',
        'mu14_status_on_normalized_anti_halfturn_A12': 'unit by inherited SIC2B33AW14; '
            'therefore adding mu14 removes this full scheme, including nilpotents. '
            'No explicit inverse computed and no old moment14 CAS artifact replayed.',
        'next_global_gate': 'Full A14 with both mu13 and mu14, compared with51.',
        'd14': 17, 'd14_dependency': 'SIC2B33LOCALTHICK20261117',
        'sources': ['https://stacks.math.columbia.edu/tag/02JN',
                    'https://stacks.math.columbia.edu/tag/03PA'],
        'boundary': 'Exact symbolic cover and multiplicity-transfer verification, '
            'conditional on the retained exact msolve length9168 and local length22. '
            'No new CAS run, global A14 computation, local d14 recomputation, extra all-order pure '
            'component, independent Buchberger verification or full Image theorem.'}
    (root/'anti_involution_cover_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
