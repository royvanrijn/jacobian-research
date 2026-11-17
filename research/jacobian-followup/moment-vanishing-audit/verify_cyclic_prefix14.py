#!/usr/bin/env python3
"""Build and replay the full A14 input, retaining both mu13 and mu14.

This extends the exact first12 circuit with the same coefficient model.
It does not run a global Groebner calculation or compute local d14.
"""
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path

from cyclic_frame_quotient import (
    S, VARIABLES, normalized_relations, moment_circuit,
    invariants_from_phase, rodrigues_point,
)
from conformal_euler_cubic import z, t, w


def odd_factorial(value):
    return math.prod(range(value, 0, -2))


def direct_control(phase, maximum):
    """Direct spatial multiplication and integration, independent of frame Wick DAG."""
    terms = S.Poly(phase, z, w, t).terms()
    assert all(coefficient.is_Integer for powers, coefficient in terms)
    powers = {(0, 0, 0): 1}
    result = []
    for n in range(1, maximum+1):
        next_power = defaultdict(int)
        for left, coefficient in powers.items():
            for right, factor in terms:
                next_power[tuple(a+b for a, b in zip(left, right))] += coefficient*int(factor)
        powers = dict(next_power)
        moment = Fraction(0)
        for (i, j, height), coefficient in powers.items():
            if i == j and height % 2 == 0:
                moment += coefficient*Fraction(math.factorial(i)*odd_factorial(height-1),
                                               odd_factorial(2*i+height+1))
        result.append(moment)
    return result


def main():
    root = Path(__file__).resolve().parent
    old_path = root/'cyclic_frame_moments_12.circuit.json.gz'
    with gzip.open(old_path, 'rt') as stream:
        old = json.load(stream)
    progress = []
    def checkpoint(n, nodes, spatial):
        progress.append({'moment': n, 'nodes': nodes, 'spatial_terms': spatial})
        if n >= 12:
            print('moment '+str(n)+': '+str(nodes)+' nodes', flush=True)
    graph, outputs = moment_circuit(14, checkpoint, node_limit=1800000)
    assert outputs[:12] == old['moment_outputs']
    assert graph.nodes[:len(old['nodes'])] == [tuple(node) for node in old['nodes']]
    assert old['variables'] == list(map(str, VARIABLES))
    artifact = root/'cyclic_frame_moments_14.circuit.json.gz'
    if not artifact.exists():
        graph.save(artifact, outputs)
    with gzip.open(artifact, 'rt') as stream:
        saved = json.load(stream)
    assert saved['moment_outputs'] == outputs
    assert saved['variables'] == list(map(str, VARIABLES))
    assert len(saved['nodes']) == len(graph.nodes)
    assert all(tuple(node) == generated for node, generated in zip(saved['nodes'], graph.nodes))
    del saved, old

    assert graph.evaluate(rodrigues_point(), outputs) == [0]*14
    metric = S.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    v = S.Matrix([1, 1, 1])
    operator = S.Matrix([[0, 0, 0], [S.Rational(1, 2), 0, 0], [0, S.Rational(1, 2), 0]])
    cubic = w**3
    values, volume = invariants_from_phase(v, operator, cubic, (z, t, w))
    assert volume == S.Rational(1, 8)
    physical = S.Matrix([z, t, w])
    phase = S.expand((v.T*metric*physical)[0]
                     +(physical.T*metric*operator*physical)[0]+cubic)
    expected = direct_control(phase, 14)
    observed = graph.evaluate(values, outputs)
    assert observed == expected
    assert expected[12] != 0 and expected[13] != 0
    local = json.loads((root/'rodrigues_length17_verification.json').read_text())
    assert local['status'] == 'PASS' and local['exact_transverse_lengths']['14'] == 17
    report = {
        'status': 'PASS', 'scope': 'Full normalized D=1/8 A14 arithmetic input only.',
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'model_sha256': hashlib.sha256((root/'cyclic_frame_quotient.py').read_bytes()).hexdigest(),
        'old_prefix12_sha256': hashlib.sha256(old_path.read_bytes()).hexdigest(),
        'prefix14_sha256': hashlib.sha256(artifact.read_bytes()).hexdigest(),
        'prefix14_nodes': len(graph.nodes), 'prefix14_compressed_bytes': artifact.stat().st_size,
        'prefix12_gate_for_gate_agreement': True,
        'moment_outputs': {str(n): output for n, output in enumerate(outputs, 1)},
        'structural_equations': list(map(str, normalized_relations())),
        'moment_equations': ['mu'+str(n) for n in range(2, 15)],
        'total_global_equations': 17, 'global_variables': 15,
        'mu13_retained': True, 'mu14_retained': True,
        'Rodrigues_moments1to14': 'all0',
        'generic_direct_control_phase': str(phase),
        'generic_control_mu13': str(expected[12]), 'generic_control_mu14': str(expected[13]),
        'progress': progress,
        'global_A14_dimension_and_length': 'UNCOMPUTED',
        'actual_local_d14': 17,
        'local_length_dependency': 'SIC2B33LOCALTHICK20261117',
        'next_exhaustion_target': 'finite dim_C A14=51; three actual local contributions17',
        'boundary': 'Exact input extension and gate-for-gate replay with direct integration '
            'controls. No global ideal membership, finite length, local multiplication '
            'table, all-order component exhaustion or full Image theorem. The local '
            'value17 is read from its separately verified supplied certificate.'}
    (root/'cyclic_prefix14_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ['progress', 'structural_equations', 'moment_outputs']}, indent=2))


if __name__ == '__main__':
    main()
