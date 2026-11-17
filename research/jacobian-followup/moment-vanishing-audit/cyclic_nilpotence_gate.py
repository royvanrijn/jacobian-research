"""Exact rank-three module formulation of the fourteen A12 certificates.

No membership or global length is asserted by this module.  The saved moment
circuit may be interpreted in this algebra without expanding its polynomials.
"""
import sympy as S

from cyclic_frame_quotient import p, r, s0, s1, s2, U, VARIABLES, DET

BASE_VARIABLES = tuple(v for v in VARIABLES if v != s2)
LINEAR_TARGETS = (p, r, s0, s1,
                  *(v for ijk, v in U.items() if ijk != (2, 2, 2)),
                  U[(2, 2, 2)]+128*s2)
C0 = (-p*p*s0*s1*s1-2*p*r*s0*s0*s1-r*r*s0**3-r*s1**3
      +S.Rational(1, 64))
C1 = p*s1*s1+3*r*s0*s1
C2 = p*s0
MONIC = s2**3-C2*s2**2-C1*s2-C0
MULTIPLICATION_S2 = S.Matrix([[0, 0, C0], [1, 0, C1], [0, 1, C2]])


def reduce_s2(expression):
    """Canonical representative modulo the Gram equation, with no localization."""
    return S.rem(S.Poly(expression, s2), S.Poly(MONIC, s2)).as_expr().expand()


def vector(expression):
    reduced = S.Poly(reduce_s2(expression), s2)
    return tuple(reduced.nth(i) for i in range(3))


def multiply(left, right, expand=False):
    z0 = left[0]*right[0]
    z1 = left[0]*right[1]+left[1]*right[0]
    z2 = left[0]*right[2]+left[1]*right[1]+left[2]*right[0]
    z3 = left[1]*right[2]+left[2]*right[1]
    z4 = left[2]*right[2]
    result = (z0+C0*z3+C0*C2*z4,
              z1+C1*z3+(C0+C1*C2)*z4,
              z2+C2*z3+(C1+C2*C2)*z4)
    return tuple(S.expand(v) for v in result) if expand else result


def target_vectors(exponent=10):
    return tuple(vector(ell**exponent) for ell in LINEAR_TARGETS)


def circuit_module_program(circuit_sha256):
    """Compact recipes, not expanded moment columns or membership witnesses."""
    return {
        'schema': 1, 'coefficient_field': 'QQ',
        'base_variables': list(map(str, BASE_VARIABLES)),
        'basis': ['1', 's2', 's2^2'],
        'monic_relation': str(MONIC),
        'multiplication_s2': [[str(x) for x in row]
                              for row in MULTIPLICATION_S2.tolist()],
        'source_circuit': 'cyclic_frame_moments_12.circuit.json.gz',
        'source_circuit_sha256': circuit_sha256,
        'generators': ['trace0', 'trace1', 'trace2']
                      +['mu'+str(n) for n in range(2, 13)],
        'column_recipes': [{'generator': name, 's2_power': k}
                           for name in ['trace0', 'trace1', 'trace2']
                           +['mu'+str(n) for n in range(2, 13)]
                           for k in range(3)],
        'target_exponent': 10,
        'target_linear_forms': list(map(str, LINEAR_TARGETS)),
        'module_membership_checks_performed': 0,
        'external_negative_certificate': 'CYCLIC_NILPOTENCE_GATE_AND_PREFIX12_OUTSIDERS.md '
                                         'proves s0^10 is not in I12.',
        'scope': 'An exact rank-three module problem with42 columns over a '
                 'fourteen-variable polynomial ring. No successful membership '
                 'or global exhaustion certificate is asserted.'}
