#!/usr/bin/env python3
"""Audit the supplied seven-target packet and export the graded module input.

The full supplied replay is separately retained. This checker validates its
hashes and rechecks local powers, cubic relation, exact rational Bezout
identities, augmented-ideal equality and every moment-circuit grading gate.
"""
from array import array
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile

from seven_target_gate import (
    S, p, r, s0, s1, s2, U, VARIABLES, K, EXPONENTS, PSI_COEFFICIENTS,
    PSI_TERMS, PSI, TARGETS, TARGET_DEGREES, TARGET_WEIGHTS, WEIGHTS,
    character, module_program,
)
from cyclic_frame_quotient import normalized_relations
from conformal_euler_cubic import z, t, w, sphere_average

ARCHIVE_SHA256 = '5ac529c2338c5404bace25091c2623e9a75710a5660455e1d58dc52b0361b020'


def replay(root, bundle, manifest):
    import resource
    def limits():
        resource.setrlimit(resource.RLIMIT_CPU, (180, 180))
        resource.setrlimit(resource.RLIMIT_AS, (3*1024**3, 3*1024**3))
    workspace = root.parents[2]
    destination = Path(tempfile.mkdtemp(prefix='rodrigues-seven-targets-', dir=workspace/'research/.scratch'))
    assert destination.resolve().is_relative_to(workspace.resolve())
    for name in ['manifest.json', *manifest['sha256']]:
        target = destination/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(bundle/name, target)
    started = time.monotonic()
    with (root/'rodrigues_seven_targets_replay.log').open('w') as log:
        result = subprocess.run([sys.executable, 'verify_gate.py'], cwd=destination,
                                stdout=log, stderr=subprocess.STDOUT, timeout=180,
                                preexec_fn=limits)
    outputs = ['global_target_nilindices.json', 'global_adapted_targets.json',
               'global_cubic_relation.json', 'odd_slice_support_certificate/odd_slice_support_data.json']
    record = {'archive_sha256': ARCHIVE_SHA256, 'returncode': result.returncode,
              'wall_seconds': time.monotonic()-started,
              'replay_directory': str(destination.relative_to(workspace)),
              'all_input_hashes_match': True,
              'output_bytes_equal': {name: (bundle/name).read_bytes() == (destination/name).read_bytes()
                                     for name in outputs},
              'resource_limits': {'wall_seconds': 180, 'cpu_seconds': 180,
                                  'address_space_bytes': 3*1024**3}}
    (root/'rodrigues_seven_targets_replay_run.json').write_text(json.dumps(record, indent=2)+'\n')
    assert result.returncode == 0 and all(record['output_bytes_equal'].values())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    bundle = root/'rodrigues_global_seven_target_certificate'
    archive = root/'rodrigues_global_seven_target_certificate.zip'
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == ARCHIVE_SHA256
    with zipfile.ZipFile(archive) as zipped:
        raw = zipped.read('rodrigues_global_seven_target_certificate/manifest.json')
    assert (bundle/'manifest.json').read_bytes() == raw
    manifest = json.loads(raw)
    for name, digest in manifest['sha256'].items():
        assert hashlib.sha256((bundle/name).read_bytes()).hexdigest() == digest
    if args.replay:
        replay(root, bundle, manifest)
    run = json.loads((root/'rodrigues_seven_targets_replay_run.json').read_text())
    assert run['archive_sha256'] == ARCHIVE_SHA256 and run['returncode'] == 0
    assert len(run['output_bytes_equal']) == 4 and all(run['output_bytes_equal'].values())
    old = root/'rodrigues_local_length17_certificate'
    assert all((bundle/'rodrigues_local_length17_certificate'/path.name).read_bytes()
               == path.read_bytes() for path in old.iterdir() if path.is_file())

    data = json.loads((old/'certificate_data.json').read_text())
    Ma = S.Matrix([[S.Rational(x) for x in row] for row in data['multiplication_a']])
    Mb = S.Matrix([[S.Rational(x) for x in row] for row in data['multiplication_b']])
    basis = [Ma**i*Mb**j for i, j in data['quotient_basis']]
    unit = S.eye(17)[:, 0]
    def operator(value):
        return sum((coefficient*basis[i] for i, coefficient in enumerate(value)), S.zeros(17))
    local = json.loads((bundle/'global_target_nilindices.json').read_text())
    vectors = {name: S.Matrix([S.Rational(x) for x in vector])
               for name, vector in local['normalized_target_vectors'].items()}
    normalized_s2 = S.Matrix([S.Rational(x) for x in local['normalized_s2']])
    # The adapter key U222 is the nilpotent T=U222+128*s2, not raw U222.
    values = {p: operator(vectors['p']), r: operator(vectors['r']),
              s0: operator(vectors['s0']), s1: operator(vectors['s1']),
              s2: operator(normalized_s2)}
    values.update({value: operator(vectors['U'+''.join(map(str, indices))])
                   for indices, value in U.items()})
    values[U[2, 2, 2]] -= 128*values[s2]
    def evaluate(expression):
        result = S.zeros(17)
        for powers, coefficient in S.Poly(expression, *VARIABLES).terms():
            term = coefficient*S.eye(17)
            for variable, exponent in zip(VARIABLES, powers):
                if exponent:
                    term = term*values[variable]**exponent
            result += term
        return result
    for value, exponent in zip(K, EXPONENTS):
        matrix = evaluate(value)
        assert matrix**exponent == S.zeros(17)
        assert matrix**(exponent-1) != S.zeros(17)
    five = S.Matrix.hstack(*(evaluate(term)*unit for term in PSI_TERMS))
    assert five.rank() == 4 and five*S.Matrix(PSI_COEFFICIENTS) == S.zeros(17, 1)
    assert evaluate(PSI) == S.zeros(17)
    for target in TARGETS:
        assert evaluate(target) == S.zeros(17)

    Z = (p, s0, s1, U[0, 0, 0], U[0, 0, 2], U[0, 2, 2], U[2, 2, 2]+128*s2)
    ideal_K = S.groebner([*normalized_relations(), *K], *VARIABLES, domain=S.QQ)
    ideal_Z = S.groebner([*normalized_relations(), *Z], *VARIABLES, domain=S.QQ)
    assert list(ideal_K) == list(ideal_Z)
    support = json.loads((bundle/'odd_slice_support_certificate/odd_slice_support_data.json').read_text())
    x = S.Symbol('x')
    for certificate in support['bezout'].values():
        first = S.Poly(S.sympify(certificate['first']), x, domain=S.QQ)
        second = S.Poly(S.sympify(certificate['second']), x, domain=S.QQ)
        left = S.Poly(S.sympify(certificate['left_over_Q']), x, domain=S.QQ)
        right = S.Poly(S.sympify(certificate['right_over_Q']), x, domain=S.QQ)
        assert left*first+right*second == S.Poly(1, x, domain=S.QQ)
    primitive_jacobian = S.Matrix(support['odd_jacobian'])
    assert primitive_jacobian.det() == 49411883664
    F0 = w+z*t-z**3/2
    raw_jacobian = S.Matrix([[n*sphere_average(F0**(n-1)*marker)
                              for marker in [w*w, t*w*w, t**3-3*z*t*w, z*z*t]]
                             for n in [3, 5, 7, 9]])
    row_scales = []
    for i in range(4):
        pivot = next(j for j in range(4) if raw_jacobian[i, j])
        scale = primitive_jacobian[i, pivot]/raw_jacobian[i, pivot]
        assert primitive_jacobian.row(i) == scale*raw_jacobian.row(i)
        row_scales.append(scale)

    assert [S.Poly(expr, *VARIABLES).total_degree() for expr in TARGETS] == list(TARGET_DEGREES)
    assert list(map(character, TARGETS)) == list(TARGET_WEIGHTS)
    assert [character(expr) for expr in normalized_relations()] == [0, 0, 11, 10]
    circuit = root/'cyclic_frame_moments_14.circuit.json.gz'
    circuit_sha256 = hashlib.sha256(circuit.read_bytes()).hexdigest()
    with gzip.open(circuit, 'rt') as stream:
        graph = json.load(stream)
    variables = {str(variable): weight for variable, weight in WEIGHTS.items()}
    grades = array('b')
    for index, node in enumerate(graph['nodes']):
        if node[0] == 'q':
            weight = 0 if node[1] else -1
        elif node[0] == 'v':
            weight = variables[node[1]]
        else:
            left, right = grades[node[1]], grades[node[2]]
            if node[0] == '+':
                assert left < 0 or right < 0 or left == right, ('grading mismatch', index)
                weight = right if left < 0 else left
            else:
                weight = -1 if left < 0 or right < 0 else (left+right) % 12
        grades.append(weight)
    assert all(grades[output] in (-1, n % 12)
               for n, output in enumerate(graph['moment_outputs'], 1))
    program = module_program(circuit_sha256)
    assert len(program['targets']) == 7 and len(program['source_column_recipes']) == 48
    assert program['memberships_verified'] == 0 and program['witness_multiplier_degree_bound'] is None
    program_path = root/'cyclic_prefix14_seven_target_module_program.json'
    program_path.write_text(json.dumps(program, indent=2)+'\n')
    report = {
        'claim': 'SIC2B33SEVENGATE20261117', 'status': 'PASS',
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'model_sha256': hashlib.sha256((root/'seven_target_gate.py').read_bytes()).hexdigest(),
        'archive_sha256': ARCHIVE_SHA256, 'four_JSON_outputs_replayed_byte_for_byte': True,
        'included_original_local_certificate_unchanged': True,
        'adapted_least_local_nilpotence_exponents': list(EXPONENTS),
        'cubic_five_vector_rank': 4,
        'Psi_primitive_integer_coefficients': list(PSI_COEFFICIENTS),
        'augmented_ambient_ideal_equality': '(K1,...,K7)=(Z) modulo full structural relations',
        'global_selected_locus_support_prefix': 6,
        'global_selected_locus_reduced_ideal_prefix': 9,
        'augmented_global_identity': 'I9+(Z)=I9+(K1,...,K7)=P',
        'primitive_integer_odd_moment_Jacobian_determinant': str(primitive_jacobian.det()),
        'raw_normalized_odd_moment_Jacobian_determinant': str(raw_jacobian.det()),
        'primitive_Jacobian_row_scales': list(map(str, row_scales)),
        'target_degrees': list(TARGET_DEGREES), 'target_characters_mod12': list(TARGET_WEIGHTS),
        'moment_circuit_nodes_graded_gate_for_gate': len(grades),
        'source_circuit_sha256': circuit_sha256,
        'module_program_sha256': hashlib.sha256(program_path.read_bytes()).hexdigest(),
        'source_column_recipes': 48, 'finite_ansatz_columns': 'UNCOMPUTED',
        'witness_degree_bound': 'UNKNOWN', 'global_memberships_verified': 0,
        'global_A14_finiteness_and_length': 'UNCOMPUTED',
        'msolve_dependency_for_this_reduction': False,
        'boundary': 'Proves a selected-locus augmented ideal identity and checks both '
            'equivalent seven-target gates, local powers and the exact mod12 exporter '
            'filter. No unrestricted I14 membership, witness-degree bound or global '
            'length51 is proved. The Jacobian integer value refers to primitive '
            'integer-scaled moment equations. No independent specialist or formal verification.'}
    (root/'seven_target_gate_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
