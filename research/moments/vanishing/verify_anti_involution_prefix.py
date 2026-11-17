#!/usr/bin/env python3
"""Exact input/staircase controls and optional bounded msolve replay.

The existence theorem uses a square complete intersection over a DVR.
It does not use an extracted modular point or an unverified Jacobian.
msolve's exact finite-field Groebner algorithm is part of the trusted base;
the leading-ideal receipt is not an independent Buchberger certificate.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile
import time

import numpy as np

from anti_involution_prefix import (
    PARAMETERS, az, aw, ez, ew, bz, bw, phase, exact_moment,
    primitive_moment, expression,
)
from conformal_euler_cubic import S, z, w, t, sphere_average
from cyclic_frame_quotient import (
    Circuit, invariants_from_phase, VARIABLES,
)

PRIME = 1073741827
STEM = 'anti_involution_prefix12_p1073741827_rabinowitsch'


def leading_staircase(path):
    text = path.read_text()
    assert '#Leading ideal data' in text
    assert '#field characteristic: '+str(PRIME) in text
    assert '#variable order:       inv, az, aw, ez, ew, bz, bw' in text
    assert '#monomial order:       graded reverse lexicographical' in text
    body = text[text.index('[')+1:text.rindex(']')]
    names = ['inv']+list(map(str, PARAMETERS))
    leading = []
    for monomial in body.split(','):
        powers = [0]*7
        for factor in monomial.strip().split('*'):
            variable, power = factor.split('^')
            powers[names.index(variable)] = int(power)
        leading.append(tuple(powers))
    assert (0,)*7 not in leading
    bounds = []
    for i in range(7):
        pure = [exponents[i] for exponents in leading
                if exponents[i] and sum(exponents) == exponents[i]]
        assert pure, 'missing pure leading power; zero-dimensionality not certified'
        bounds.append(min(pure))
    staircase = np.ones(bounds, dtype=bool)
    for exponents in leading:
        if all(a < b for a, b in zip(exponents, bounds)):
            staircase[tuple(slice(a, b) for a, b in zip(exponents, bounds))] = False
    count = int(staircase.sum())
    assert count == 9168 and staircase[(0,)*7]
    return len(leading), dict(zip(names, bounds)), count


def rational_evaluation(terms, values):
    return sum(coefficient*math.prod(value**exponent for value, exponent in zip(values, powers))
               for powers, coefficient in terms.items())


def replay(root, expected_source):
    """One bounded exact sparse run, retaining the primary log and leading ideal."""
    import resource
    def limits():
        resource.setrlimit(resource.RLIMIT_CPU, (220, 220))
        resource.setrlimit(resource.RLIMIT_AS, (4*1024**3, 4*1024**3))
    with tempfile.TemporaryDirectory(prefix='anti-prefix12-') as folder:
        destination = Path(folder)/'leading.out'
        command = ['msolve', '-f', str(root/(STEM+'.ms')), '-o', str(destination),
                   '-t', '2', '-v', '2', '-l', '2', '-m', '2048',
                   '-g', '1', '-P', '0', '-c', '0']
        started = time.monotonic()
        with (root/(STEM+'_replay.log')).open('w') as log:
            process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                     timeout=120, preexec_fn=limits)
        assert process.returncode == 0
        assert destination.read_bytes() == expected_source
        return {'same_leading_output': True, 'seconds': time.monotonic()-started,
                'linear_algebra': 'deterministic exact sparse (-l2)',
                'scope': 'Same-CAS replay, not independent Buchberger verification.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    assert S.isprime(PRIME)
    F = phase()
    involution = {z: -z, w: -w, t: t}
    assert S.expand(F.subs(involution, simultaneous=True)+F) == 0
    cubic = S.expand(F-z-w-t*(bz*z+bw*w))
    assert S.expand(2*S.diff(cubic, z, w)+S.diff(cubic, t, 2)) == 0
    metric = S.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    v = S.Matrix([1, 0, 1])
    B = S.Matrix([[0, bw/2, 0], [bz/2, 0, bw/2], [0, bz/2, 0]])
    assert B.T*metric == metric*B and S.trace(B) == 0
    determinant = S.Matrix.hstack(v, B*v, B*B*v).det().factor()
    assert S.expand(determinant-(bz-bw)*(bz+bw)**2/8) == 0
    assert (v.T*metric*v)[0] == 2

    saved = json.loads((root/'anti_involution_moments12.json').read_text())
    primitive, scales = {}, {}
    for n in range(2, 13, 2):
        primitive[n], scales[n] = primitive_moment(n)
        assert primitive[n] == {tuple(k): value for k, value
                                in saved['primitive_moments'][str(n)]}
        assert str(scales[n]) == saved['moment_scales'][str(n)]
        assert scales[n].numerator % PRIME and scales[n].denominator % PRIME
        assert all(value.denominator % PRIME for value in exact_moment(n).values())
    source = ('inv,'+','.join(map(str, PARAMETERS))+'\n'+str(PRIME)+'\n'
              +',\n'.join([expression(primitive[n], PRIME).replace(' ', '')
                            for n in range(2, 13, 2)]
                           +['inv*(bz-bw)*(bz+bw)-1'])+'\n')
    actual = (root/(STEM+'.ms')).read_text()
    assert re.sub(r'\s', '', actual) == re.sub(r'\s', '', source)

    # Independent Wick/covariance circuit controls for every even moment,
    # at two nonsingular cyclic frames of volume1/8.
    import gzip
    with gzip.open(root/'cyclic_frame_moments_12.circuit.json.gz', 'rt') as stream:
        saved_circuit = json.load(stream)
    graph = Circuit()
    graph.nodes = saved_circuit['nodes']
    controls = []
    for point in [(1, 2, 3, 4, 1, 0), (2, -1, -2, 1, 0, -1)]:
        substitution = dict(zip(PARAMETERS, point))
        values, volume = invariants_from_phase(v, B.subs(substitution),
                                              cubic.subs(substitution), (z, t, w))
        assert volume == S.Rational(1, 8)
        observed = graph.evaluate(values, saved_circuit['moment_outputs'])
        for n in range(1, 13):
            expected = (Fraction(0) if n % 2 else
                        scales[n]*rational_evaluation(primitive[n], point))
            assert observed[n-1] == expected
        controls.append({'parameters': list(point), 'moments1to12': 'all agree'})
    # A third direct integration control does not use the frame adapter.
    point = dict(zip(PARAMETERS, (1, 2, 3, 4, 5, 6)))
    for n in (2, 4, 6):
        assert sphere_average(F.subs(point)**n) == S.Rational(
            scales[n]*rational_evaluation(primitive[n], tuple(point.values())))

    leading_path = root/(STEM+'.out')
    basis_size, bounds, length = leading_staircase(leading_path)
    run = json.loads((root/(STEM+'_run.json')).read_text())
    assert run['status'] == 'completed' and run['returncode'] == 0
    log = (root/(STEM+'.log')).read_text()
    assert 'Dimension of quotient: 9168' in log
    assert '#equations                       7' in log
    assert '#variables                       7' in log
    replay_result = replay(root, leading_path.read_bytes()) if args.replay else None
    previous = root/'anti_involution_prefix12_verification.json'
    if replay_result is None and previous.exists():
        replay_result = json.loads(previous.read_text()).get('bounded_same_CAS_replay')
    files = ['anti_involution_prefix.py', 'anti_involution_moments12.json',
             STEM+'.ms', STEM+'.out', STEM+'.log', STEM+'_run.json']
    report = {
        'claim': 'SIC2B33CYCLICEXTRA20261116', 'status': 'PASS',
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'hashes': {name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in files},
        'prime': PRIME, 'variables': 7, 'equations': 7,
        'software': {'sympy': S.__version__, 'numpy': np.__version__,
                     'primary_msolve_package': 'Ubuntu msolve0.6.5-1build2'},
        'cyclic_inverse_equation': 'inv*(bz-bw)*(bz+bw)-1',
        'msolve_leading_basis_elements': basis_size, 'pure_leading_power_bounds': bounds,
        'modular_dimension': 0, 'modular_length': length,
        'moment_circuit_controls': controls,
        'odd_moments': 'identically zero by a determinant-one half-turn',
        'square_lifting_argument': 'At any special closed point the seven equations '
            'followed by the uniformizer form a system of parameters in an '
            'eight-dimensional regular local ring. Hence they are a regular sequence; '
            'the uniformizer is a nonzerodivisor and the characteristic-zero fibre '
            'is nonempty. No smoothness or unramified lift is asserted.',
        'characteristic_zero_chart_s0': 2,
        'normalized_A12_consequence': 's0^10 is NOT in I12; length66 is false.',
        'A12_vector_space_dimension_lower_bound': 69,
        'A12_global_finiteness': 'UNCOMPUTED',
        'A12_exact_length': 'UNCOMPUTED',
        'later_failure': 'mu14 is nonzero at every chart first12 solution by the '
                         'inherited anti-Weyl theorem SIC2B33AW14. Its historical '
                         'moment14 artifact is not reconstructed or replayed here.',
        'bounded_same_CAS_replay': replay_result,
        'trusted_computation': 'msolve exact finite-field Groebner/saturation algorithms; '
            'small explicit seven-equation presentation, retained primary leading ideal '
            'and bounded replay. No independent Buchberger or formal certificate.',
        'sources': ['https://stacks.math.columbia.edu/tag/02JN',
                    'https://stacks.math.columbia.edu/tag/00R4'],
        'boundary': 'Certifies existence of additional first12-moment candidates on '
            'D nonzero. Supplies neither explicit coordinates nor a nonsingular '
            'modular point, another all-order pure phase, an exact global length '
            'or the full (3,3) Image theorem.'}
    previous.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
