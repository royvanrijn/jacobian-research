#!/usr/bin/env python3
"""Exact complete-frame, moment-circuit and global residue controls.

Does not compute the global dimension or length of A12.
"""
from pathlib import Path
import gzip
import hashlib
import json
import math

from cyclic_frame_quotient import (
    S, p, r, s0, s1, s2, D, VARIABLES, U, G, C, DET,
    normalized_relations, trace_relations, moment_circuit,
    hilbert_residue, quadratic_moments, invariants_from_phase,
    krylov_stratum, rodrigues_point,
)
from conformal_euler_cubic import z, w, t, q, sphere_average, harmonic_basis


def main():
    root = Path(__file__).resolve().parent
    assert S.expand((C**3-p*C-r*S.eye(3))).is_zero_matrix
    assert C.T*G == G*C
    assert S.trace(C) == 0
    weighted_euler = sum(weight*variable*S.diff(DET, variable)
                         for variable, weight in [(p, 2), (r, 3), (s0, 2), (s1, 3), (s2, 4)])
    assert S.expand(weighted_euler-12*DET) == 0

    metric = S.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    B0 = S.Matrix([[0, 0, 0], [S.Rational(1, 2), 0, 0], [0, S.Rational(1, 2), 0]])
    original_v = S.Matrix([1, 0, 0])
    values, determinant = invariants_from_phase(original_v, B0, -z**3/2, (z, t, w))
    assert determinant == S.Rational(1, 8)
    assert values == rodrigues_point()
    assert all(S.simplify(f.subs(values)) == 0 for f in normalized_relations())
    reflection = S.diag(1, -1, 1)
    reflected, reflected_volume = invariants_from_phase(
        reflection*original_v, reflection*B0*reflection, -z**3/2, (z, t, w))
    assert reflected == values and reflected_volume == -S.Rational(1, 8)

    alpha = S.Symbol('alpha')
    cube_relation = S.Poly(alpha**3-1, alpha)
    three_points = rodrigues_point(alpha/4)
    for relation in normalized_relations():
        assert S.rem(S.Poly(relation.subs(three_points), alpha), cube_relation).is_zero
    assert S.gcd(alpha**3-1, 3*alpha*alpha) == 1
    assert S.expand((w+z*t-z**3/2).subs({z: -z, w: -w}, simultaneous=True)
                    +(w+z*t-z**3/2)) == 0

    progress = []
    def checkpoint(n, nodes, spatial):
        progress.append({'moment': n, 'nodes': nodes, 'spatial_terms': spatial})
        print('moment '+str(n)+': '+str(nodes)+' circuit nodes', flush=True)
    graph, outputs = moment_circuit(12, checkpoint)
    artifact = root/'cyclic_frame_moments_12.circuit.json.gz'
    if not artifact.exists():
        graph.save(artifact, outputs)
    with gzip.open(artifact, 'rt') as stream:
        saved = json.load(stream)
    assert saved['variables'] == list(map(str, VARIABLES))
    assert saved['moment_outputs'] == outputs
    assert len(saved['nodes']) == len(graph.nodes)
    assert all(tuple(raw) == generated for raw, generated in zip(saved['nodes'], graph.nodes))
    del saved
    assert all(v == 0 for v in graph.evaluate(values, outputs))

    full_cubic = z**3+w*(t*t-z*w/2)+t**3-3*z*w*t
    semisimple_B = S.Matrix([[0, 1, 4], [S.Rational(1, 2), 0, 1], [0, S.Rational(1, 2), 0]])
    repeated_B = S.Matrix([[0, 3, 8], [S.Rational(1, 2), 0, 3], [0, S.Rational(1, 2), 0]])
    tests = [
        ('nilpotent_B_nonisotropic_v', S.Matrix([1, 1, 1]), B0, full_cubic),
        ('nonnilpotent_B', original_v, semisimple_B, full_cubic),
        ('repeated_eigenvalue_B', original_v, repeated_B, full_cubic),
        ('mu2_zero_mu3_nonzero', S.Matrix([1, 1, -S.Rational(1, 2)]), B0, S.Integer(0)),
    ]
    generic_checks = []
    for label, vector, operator, cubic in tests:
        frame_values, volume = invariants_from_phase(vector, operator, cubic, (z, t, w))
        assert volume == S.Rational(1, 8)
        evaluated = graph.evaluate(frame_values, outputs[:4])
        physical = S.Matrix([z, t, w])
        phase = S.expand((vector.T*metric*physical)[0]
                         +(physical.T*metric*operator*physical)[0]+cubic)
        expected = [sphere_average(phase**n) for n in range(1, 5)]
        assert evaluated == expected
        if label == 'mu2_zero_mu3_nonzero':
            assert expected[1] == 0 and expected[2] == S.Rational(2, 5)
        generic_checks.append({'control': label, 'moments1to4': list(map(str, expected))})

    # The invariant harmonic norm is nondegenerate in every nonsingular
    # metric, by isometry. Check its seven-dimensional standard model.
    cubic_basis = tuple(harmonic_basis()[3].values())
    norm = S.Matrix([[sphere_average(f*g) for g in cubic_basis] for f in cubic_basis])
    assert norm.det() != 0 and norm.rank() == 7

    quadratic = quadratic_moments(6)
    assert S.expand(quadratic[1]-(4*p+5*s0)/15) == 0
    assert S.expand(quadratic[2]-2*(4*r+7*s1)/35) == 0
    basis = S.groebner(quadratic[1:], s0, s1, s2, p, r, order='lex')
    assert basis.reduce(DET)[1] == 0
    normalized = S.groebner([*quadratic[1:], DET+S.Rational(1, 64)], s0, s1, s2, p, r, order='lex')
    assert list(normalized) == [1]
    quadratic_point, volume = invariants_from_phase(original_v, semisimple_B, S.Integer(0), (z, t, w))
    quadratic_evaluated = graph.evaluate(quadratic_point, outputs[:6])
    assert quadratic_evaluated == [S.factor(f.subs(quadratic_point)) for f in quadratic]

    diagonal = S.diag(1, -2, 1)
    upper_nilpotent = S.Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    rank2_nondegenerate = S.Matrix([[S.Rational(1, 2), 0, S.Rational(1, 2)],
                                  [0, -1, 0], [S.Rational(1, 2), 0, S.Rational(1, 2)]])
    strata = [
        (S.zeros(3, 1), diagonal, (0, 0)),
        (S.Matrix([1, 0, 0]), diagonal, (1, 0)),
        (S.Matrix([0, 1, 0]), diagonal, (1, 1)),
        (S.Matrix([0, 1, 0]), upper_nilpotent, (2, 1)),
        (S.Matrix([1, 1, 1]), rank2_nondegenerate, (2, 2)),
        (original_v, B0, (3, 3)),
    ]
    for vector, operator, expected in strata:
        assert operator.T*metric == metric*operator and S.trace(operator) == 0
        assert krylov_stratum(vector, operator, metric) == expected

    contour, schur = hilbert_residue()
    assert contour == schur == S.Rational(317, 230400)
    residue = math.factorial(12)*contour
    assert residue == 659043 and 22*S.Rational(1, 2) == 11
    report = {
        'claims': ['SIC2B33CYCLIC20261116', 'SIC2B33HILBERT20261116'],
        'status': 'PASS',
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'model_sha256': hashlib.sha256((root/'cyclic_frame_quotient.py').read_bytes()).hexdigest(),
        'circuit_sha256': hashlib.sha256(artifact.read_bytes()).hexdigest(),
        'sympy_version': S.__version__,
        'normalized_quotient_variables': list(map(str, VARIABLES)),
        'normalized_quotient_relations': list(map(str, normalized_relations())),
        'circuit_prefix': 12, 'circuit_nodes': len(graph.nodes),
        'circuit_compressed_bytes': artifact.stat().st_size,
        'full_circuit_replay': True, 'generic_original_phase_controls': generic_checks,
        'orientation_sign_control': 'Same G,C,U with opposite D are retained as distinct SO(q) data.',
        'harmonic_norm_rank': 7,
        'quadratic_six_moment_gram_determinant_remainder': '0',
        'quadratic_normalized_prefix6_basis': ['1'],
        'quadratic_control_scope': 'Entire H3=0 regression only; not global A12 exhaustion.',
        'known_global_Krull_dimensions': {'A1': 11, 'A2': 10, 'A3': 9},
        'normalized_Rodrigues_quotient_points': 3,
        'local_first12_length_per_point': 22,
        'A12_finite_length_exhaustion_target': 66,
        'A12_global_Krull_dimension': 'UNCOMPUTED',
        'A12_global_vector_space_dimension': 'UNCOMPUTED',
        'hilbert_gamma0_residue': str(contour), 'hilbert_gamma0_Schur': str(schur),
        'hypothetical_regular_prefix12_residue': str(residue),
        'known_Rodrigues_component_contribution': '11',
        'D_zero_rank_Gram_rank_strata': [[0, 0], [1, 0], [1, 1], [2, 1], [2, 2]],
        'boundary': 'Complete quotient and exact moment input, with written low-prefix dimension '
                    'and Hilbert obstruction proofs. No global A12 dimension/length computation, '
                    'D-nonzero exhaustion, extra pure component, or full (3,3) Image theorem. '
                    'No independent specialist, external or formal verification.',
    }
    (root/'cyclic_frame_quotient_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in
                      {'normalized_quotient_relations', 'generic_original_phase_controls'}}, indent=2))


if __name__ == '__main__':
    main()
