#!/usr/bin/env python3
"""Replay the exact algebraic reduction; does not prove any target membership."""
import gzip
import hashlib
import json
from pathlib import Path

from cyclic_frame_quotient import (
    S, p, r, s0, s1, s2, U, VARIABLES, DET, trace_relations,
    rodrigues_point,
)
from cyclic_nilpotence_gate import (
    BASE_VARIABLES, LINEAR_TARGETS, C0, C1, C2, MONIC,
    MULTIPLICATION_S2, reduce_s2, vector, multiply, target_vectors,
    circuit_module_program,
)


def main():
    root = Path(__file__).resolve().parent
    assert len(LINEAR_TARGETS) == len(BASE_VARIABLES) == 14
    assert len(VARIABLES) == 15
    assert S.expand(MONIC+DET+S.Rational(1, 64)) == 0
    assert S.Poly(MONIC, s2).LC() == 1
    assert (MULTIPLICATION_S2**3-C2*MULTIPLICATION_S2**2
            -C1*MULTIPLICATION_S2-C0*S.eye(3)).applyfunc(S.expand).is_zero_matrix
    # Verify all basis products, including the degree-four reduction.
    for i in range(3):
        for j in range(3):
            left, right = [0]*3, [0]*3
            left[i] = right[j] = 1
            assert multiply(left, right, expand=True) == vector(s2**(i+j))
    # A nonlinear control exercises arbitrary coefficients in the base ring.
    a = p+r*s2+s0*s2*s2
    b = s1+U[(0, 1, 2)]*s2+(p*s0+r)*s2*s2
    assert multiply(vector(a), vector(b), expand=True) == vector(a*b)

    substitution = {p: 0, r: 0, s0: 0, s1: 0,
                    **{u: 0 for ijk, u in U.items() if ijk != (2, 2, 2)},
                    U[(2, 2, 2)]: -128*s2}
    assert S.expand((DET+S.Rational(1, 64)).subs(substitution)) == -s2**3+S.Rational(1, 64)
    assert all(S.expand(f.subs(substitution)) == 0 for f in trace_relations())
    assert S.gcd(s2**3-S.Rational(1, 64), 3*s2*s2) == 1
    assert all(ell.subs(rodrigues_point()) == 0 for ell in LINEAR_TARGETS)
    # The fourteen forms and s2 give an invertible affine-linear coordinate system.
    assert S.Matrix([*LINEAR_TARGETS, s2]).jacobian(VARIABLES).det() in (1, -1)
    targets = target_vectors()
    assert all(targets[i] == (LINEAR_TARGETS[i]**10, 0, 0) for i in range(13))
    raw = LINEAR_TARGETS[-1]**10
    reduced = sum(targets[-1][i]*s2**i for i in range(3))
    quotient, remainder = S.div(S.Poly(raw-reduced, s2), S.Poly(MONIC, s2))
    assert remainder.is_zero and quotient.degree() == 7
    source = root/'cyclic_frame_moments_12.circuit.json.gz'
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    with gzip.open(source, 'rt') as stream:
        circuit = json.load(stream)
    assert circuit['variables'] == list(map(str, VARIABLES))
    assert len(circuit['moment_outputs']) == 12
    program = circuit_module_program(source_sha)
    assert len(program['column_recipes']) == 42
    (root/'cyclic_nilpotence_module_program.json').write_text(json.dumps(program, indent=2)+'\n')
    report = {
        'claim': 'SIC2B33NILGATE20261116', 'status': 'PASS',
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'model_sha256': hashlib.sha256((root/'cyclic_nilpotence_gate.py').read_bytes()).hexdigest(),
        'moment_circuit_sha256': source_sha,
        'linear_targets': list(map(str, LINEAR_TARGETS)),
        'target_exponent': 10, 'known_local_lengths': [22, 22, 22],
        'known_local_maximal_ideal_power_zero': 10,
        'rank_three_module_base_variables': 14, 'module_columns': 42,
        'quotient_by_linear_targets': 'QQ[s2]/(s2^3-1/64)',
        'module_membership_checks_performed': 0,
        'A12_global_length': 'UNCOMPUTED',
        'boundary': 'Checks the explicit monic algebra, target reduction and support '
                    'description used in the written fourteen-certificate equivalence. '
                    'This checker does not verify a target membership or an outsider. '
                    'The separate half-turn square-system certificate disproves length66. '
                    'The equivalence uses the separately certified actual local index10 '
                    'and length22, not merely their weighted initial algebra.'}
    (root/'cyclic_nilpotence_gate_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
