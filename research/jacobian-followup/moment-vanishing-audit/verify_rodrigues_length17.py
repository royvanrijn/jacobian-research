#!/usr/bin/env python3
"""Audit the supplied exact local certificate, with optional full jet replay.

Preserves the uploaded bundle byte-for-byte. Rebuilds finite relation rows,
exports actual membership witnesses, checks multiplication and the all-order
shear bounds. The written Nakayama/completion argument is in the proof note.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile

import sympy as S

from conformal_euler_cubic import z, t, w, q, sphere_average

ARCHIVE_SHA256 = '920fc19dc1a43f89b7de809d0b3a3c1098e3472b973649e1d392561c8b291f87'
K_GENERATORS = ((5, 0), (4, 1), (3, 2), (2, 4), (1, 5), (0, 7))


def matrix(raw):
    return S.Matrix([[S.Rational(value) for value in row] for row in raw])


def encoded(values):
    return [str(value) for value in values]


def replay(root, bundle):
    import resource
    def limits():
        resource.setrlimit(resource.RLIMIT_CPU, (180, 180))
        resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    workspace = root.parents[2]
    destination = Path(tempfile.mkdtemp(prefix='rodrigues-length17-', dir=workspace/'tmp'))
    assert destination.resolve().is_relative_to(workspace.resolve())
    manifest = json.loads((bundle/'manifest.json').read_text())
    for name in ['manifest.json', *manifest['sha256']]:
        shutil.copyfile(bundle/name, destination/name)
    command = [sys.executable, 'verify_certificate.py']
    started = time.monotonic()
    with (root/'rodrigues_length17_replay.log').open('w') as log:
        result = subprocess.run(command, cwd=destination, stdout=log,
                                stderr=subprocess.STDOUT, timeout=180, preexec_fn=limits)
    record = {'archive_sha256': ARCHIVE_SHA256, 'returncode': result.returncode,
              'wall_seconds': time.monotonic()-started, 'command': command,
              'replay_directory': str(destination.relative_to(workspace)),
              'resource_limits': {'wall_seconds': 180, 'cpu_seconds': 180,
                                  'address_space_bytes': 2*1024**3},
              'uploaded_manifest_hashes_match': True,
              'generated_certificate_data_byte_equal':
                  (destination/'certificate_data.json').read_bytes()
                  == (bundle/'certificate_data.json').read_bytes()}
    (root/'rodrigues_length17_replay_run.json').write_text(json.dumps(record, indent=2)+'\n')
    assert result.returncode == 0 and record['generated_certificate_data_byte_equal']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    bundle = root/'rodrigues_local_length17_certificate'
    archive = root/'rodrigues_local_length17_certificate.zip'
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == ARCHIVE_SHA256
    with zipfile.ZipFile(archive) as zipped:
        original_manifest = zipped.read('rodrigues_local_length17_certificate/manifest.json')
    assert (bundle/'manifest.json').read_bytes() == original_manifest
    manifest = json.loads(original_manifest)
    for name, digest in manifest['sha256'].items():
        assert hashlib.sha256((bundle/name).read_bytes()).hexdigest() == digest
    if args.replay:
        replay(root, bundle)
    run = json.loads((root/'rodrigues_length17_replay_run.json').read_text())
    assert run['archive_sha256'] == ARCHIVE_SHA256
    assert run['returncode'] == 0 and run['generated_certificate_data_byte_equal']

    data = json.loads((bundle/'certificate_data.json').read_text())
    monomials = [tuple(exponents) for exponents in data['monomials']]
    assert monomials == [(i, degree-i) for degree in range(8) for i in range(degree+1)]
    residuals = {int(n): {tuple(exponents): S.Rational(value) for exponents, value in terms}
                 for n, terms in data['residuals'].items()}
    M14, M28 = matrix(data['first14_matrix']['entries']), matrix(data['first28_matrix']['entries'])
    for raw, observed in [(data['first14_matrix'], M14), (data['first28_matrix'], M28)]:
        reconstructed = S.Matrix([[residuals[n].get((i-u, j-v), S.S.Zero)
                                    for i, j in monomials] for n, u, v in raw['row_labels']])
        assert reconstructed == observed
    assert M14.shape == (32, 36) and M28.shape == (50, 36)
    assert M14.rank() == M28.rank() == 19
    truncated_lengths = {}
    for prefix, expected in [(12, 20), (13, 18), (14, 17)]:
        rows = [row for row, label in zip(M14.tolist(), data['first14_matrix']['row_labels'])
                if label[0] <= prefix]
        truncated_lengths[prefix] = 36-S.Matrix(rows).rank()
        assert truncated_lengths[prefix] == expected

    # Exact witnesses u*M14=e_(i,j), in the original32 row coordinates.
    independent_rows = M14.T.rref()[1]
    pivot_columns = M14.rref()[1]
    basis_rows = M14[list(independent_rows), :]
    inverse = basis_rows[:, list(pivot_columns)].inv()
    def membership(row):
        coefficients = row[:, list(pivot_columns)]*inverse
        assert coefficients*basis_rows == row
        full = [S.S.Zero]*32
        for index, value in zip(independent_rows, coefficients):
            full[index] = value
        assert S.Matrix(1, 32, full)*M14 == row
        return encoded(full)
    def monomial_row(exponent):
        row = S.zeros(1, 36)
        row[0, monomials.index(exponent)] = 1
        return row
    witnesses = {'row_labels': data['first14_matrix']['row_labels'],
                 'monomials': monomials, 'K_memberships': {}, 'tail_memberships': {}}
    for exponent in K_GENERATORS:
        witnesses['K_memberships'][str(exponent)] = membership(monomial_row(exponent))
    for exponent in monomials:
        if sum(exponent) == 7 or 3*exponent[0]+2*exponent[1] >= 13:
            membership(monomial_row(exponent))
    for n in range(15, 25):
        witnesses['tail_memberships'][str(n)] = membership(
            S.Matrix(1, 36, [residuals[n].get(exponent, 0) for exponent in monomials]))
    assert all(not residuals[n] for n in range(22, 29))
    assert all(any(i >= u and j >= v for u, v in K_GENERATORS)
               for i, j in monomials if i+j == 7)
    for i in range(14):
        for j in range(14):
            assert (3*i+2*j >= 13) == any(i >= u and j >= v for u, v in K_GENERATORS)

    # Quotient reduction and the actual exported multiplication operators.
    reduced, pivots = M14.rref()
    free = [index for index in range(36) if index not in pivots]
    basis = [tuple(exponents) for exponents in data['quotient_basis']]
    assert basis == [monomials[index] for index in free]
    reductions = {}
    for j, index in enumerate(free):
        reductions[monomials[index]] = S.eye(17)[:, j]
    for row, index in enumerate(pivots):
        reductions[monomials[index]] = S.Matrix([-reduced[row, j] for j in free])
    assert reductions[0, 6] != S.zeros(17, 1)
    Ma, Mb = matrix(data['multiplication_a']), matrix(data['multiplication_b'])
    for operator, step in [(Ma, (1, 0)), (Mb, (0, 1))]:
        expected = S.Matrix.hstack(*(reductions.get((i+step[0], j+step[1]), S.zeros(17, 1))
                                    for i, j in basis))
        assert operator == expected
    assert Ma*Mb == Mb*Ma
    unit = S.eye(17)[:, basis.index((0, 0))]
    for j, (i, k) in enumerate(basis):
        assert Ma**i*Mb**k*unit == S.eye(17)[:, j]
    dimensions = [S.Matrix.hstack(*(Ma**i*Mb**(degree-i) for i in range(degree+1))).rank()
                  for degree in range(8)]
    assert dimensions == [17, 16, 14, 11, 7, 4, 1, 0]
    assert Mb**6*unit != S.zeros(17, 1)
    hilbert = [left-right for left, right in zip(dimensions, dimensions[1:])]
    assert hilbert == [1, 2, 3, 4, 3, 3, 1]

    # Check the all-order support bound, including the actual quadratic IFT term.
    a, b = S.symbols('a b')
    F0 = w+z*t-z**3/2
    Qa = w*w-(t**3-3*z*w*t)/10+S.Rational(17, 70)*t
    Qb = w*t+S.Rational(2, 5)*z*(t*t-z*w/2)-z/5
    T1 = a*Qa+b*Qb
    T2 = (a*a*(-S.Rational(4, 637)*t*w+w**3/4-S.Rational(57, 6370)*z)
          +a*b*(-S.Rational(4, 385)*t*t+t*w*w/2+S.Rational(4, 385)*w*z)
          +b*b*(t*t*w/5-w*w*z/10+w/75))
    quadratic = S.Poly(T2, a, b)
    square = S.Poly(S.expand(T1*T1), a, b)
    for n in range(1, 11):
        for exponents in [(2, 0), (1, 1), (0, 2)]:
            lhs = n*sphere_average(F0**(n-1)*quadratic.coeff_monomial(exponents))
            if n >= 2:
                lhs += S.binomial(n, 2)*sphere_average(F0**(n-2)*square.coeff_monomial(exponents))
            assert lhs == 0
    inverse_shear = {z: z, t: t+z*z, w: w-z*t-z**3/2}
    pullback = lambda marker: S.expand(S.sympify(marker).subs(inverse_shear, simultaneous=True))
    assert pullback(q) == q and pullback(F0) == w
    def maximum_weight(marker):
        return max(i-k for (i, j, k), value in S.Poly(pullback(marker), z, t, w).terms() if value)
    tangent_weights = [maximum_weight(Qa), maximum_weight(Qb)]
    quadratic_weights = [maximum_weight(quadratic.coeff_monomial(exponents))
                          for exponents in [(2, 0), (1, 1), (0, 2)]]
    assert tangent_weights == [4, 3] and quadratic_weights == [9, 8, 7]
    harmonic_basis = [1, z, t, w, z*z, z*t, t*t-z*w, w*t, w*w,
                      z**3, z*z*t, z*(t*t-z*w/2), t**3-3*z*w*t,
                      w*(t*t-z*w/2), w*w*t, w**3]
    assert max(map(maximum_weight, harmonic_basis)) == 9
    assert all(5*i+4*j <= 2*(3*i+2*j) for i in range(14) for j in range(14))

    witness_path = root/'rodrigues_length17_row_witnesses.json'
    witness_path.write_text(json.dumps(witnesses, separators=(',', ':'))+'\n')
    # The settled actual index7 sharpens the earlier global membership gate.
    from cyclic_nilpotence_gate import LINEAR_TARGETS, BASE_VARIABLES, target_vectors
    from cyclic_frame_quotient import VARIABLES
    targets = target_vectors(7)
    assert len(targets) == len(LINEAR_TARGETS) == 14
    assert all(targets[i] == (LINEAR_TARGETS[i]**7, 0, 0) for i in range(13))
    generators = ['trace0', 'trace1', 'trace2']+['mu'+str(n) for n in range(2, 15)]
    assert len(generators) == 16
    circuit = root/'cyclic_frame_moments_14.circuit.json.gz'
    program = {
        'schema': 1, 'coefficient_field': 'QQ', 'prefix': 14,
        'source_circuit': circuit.name,
        'source_circuit_sha256': hashlib.sha256(circuit.read_bytes()).hexdigest(),
        'base_variables': list(map(str, BASE_VARIABLES)), 'basis': ['1', 's2', 's2^2'],
        'generators_after_monic_Gram_reduction': generators,
        'column_recipes': [{'generator': name, 's2_power': k}
                           for name in generators for k in range(3)],
        'target_exponent': 7, 'target_linear_forms': list(map(str, LINEAR_TARGETS)),
        'reduced_target_vectors': [[str(value) for value in vector] for vector in targets],
        'known_local_lengths': [17, 17, 17],
        'equivalent_global_length_target': 51,
        'memberships_verified': 0,
        'scope': 'Exact fourteen-seventh-power equivalent exhaustion target, with48 '
                 'module columns. No membership or global length51 is proved.'}
    assert len(program['column_recipes']) == 48
    (root/'cyclic_prefix14_nilpotence_module_program.json').write_text(json.dumps(program, indent=2)+'\n')
    report = {
        'claim': 'SIC2B33LOCALTHICK20261117', 'status': 'PASS',
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'archive_sha256': ARCHIVE_SHA256,
        'supplied_certificate_data_sha256': manifest['sha256']['certificate_data.json'],
        'row_witnesses_sha256': hashlib.sha256(witness_path.read_bytes()).hexdigest(),
        'uploaded_verifier_replayed': True, 'generated_data_byte_equal': True,
        'exact_transverse_lengths': {'12': 22, '13': 18, '14': 17, 'infinity': 17},
        'ordinary_truncation_degree': 8, 'ordinary_truncation_lengths': truncated_lengths,
        'relation_matrix_shapes': [[32, 36], [50, 36]], 'relation_ranks': [19, 19],
        'K_generators_with_exact_row_witnesses': K_GENERATORS,
        'finite_tail_indices_with_exact_row_witnesses': list(range(15, 25)),
        'inverse_shear_tangent_maximum_weights': tangent_weights,
        'inverse_shear_quadratic_maximum_weights': quadratic_weights,
        'higher_correction_marker_maximum_weight': 9,
        'weighted_filtration_vanishing': 'F^13 L14=0',
        'automatic_residual_tail_vanishing': 'n>24',
        'actual_multiplication_matrix_sizes': [[17, 17], [17, 17]],
        'maximal_ideal_power_dimensions': dimensions, 'hilbert_function': hilbert,
        'nilpotence_index': 7, 'minimal_consecutive_local_prefix': 14,
        'global_A14_exhaustion_target': 51,
        'A14_equivalent_nilpotence_target_exponent': 7,
        'A14_module_generator_columns': 48, 'global_memberships_verified': 0,
        'global_A14_finiteness_and_length': 'UNCOMPUTED',
        'msolve_dependency_for_this_local_result': False,
        'sources': ['https://stacks.math.columbia.edu/tag/07RC',
                    'https://stacks.math.columbia.edu/tag/00MC'],
        'boundary': 'Supplied exact rational Taylor verifier replayed byte-for-byte; '
            'separate finite-row, multiplication and shear-weight cross-check. '
            'Uses the written Nakayama and completion argument to identify the full '
            'ordinary local algebra. No new local branch search, global finiteness '
            'or length computation, independent specialist, formal verification '
            'or full (3,3) Image theorem.'}
    (root/'rodrigues_length17_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
