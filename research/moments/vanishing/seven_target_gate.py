"""The seven equivalent A14 targets and their rational mod12 grading.

Only the targets and exact filtering rules are defined here. Witness degree
bounds, ideal memberships and global finiteness remain open.
"""
import sympy as S

from cyclic_frame_quotient import p, r, s0, s1, s2, U, VARIABLES
from cyclic_nilpotence_gate import BASE_VARIABLES, vector

K = (4*p+5*s0, s0, s1, U[0, 0, 0], U[0, 0, 2]-6*U[0, 1, 1],
     4*U[0, 2, 2]-13*U[1, 1, 2], U[2, 2, 2]+128*s2+S.Rational(896, 75)*p*p)
EXPONENTS = (3, 7, 5, 3, 2, 3, 2)
PSI_COEFFICIENTS = (3814472897542065, 9204109983515700, -9625759731192528,
                    -10061220068131439, -5712042093615)
PSI_TERMS = (s0**3, s1**2, s2*K[0], K[0]*s0**2, K[5])
PSI = S.expand(sum(coefficient*term for coefficient, term in zip(PSI_COEFFICIENTS, PSI_TERMS)))
POWER_TARGETS = tuple(S.expand(value**exponent) for value, exponent in zip(K, EXPONENTS))
TARGETS = (POWER_TARGETS[0], PSI, *POWER_TARGETS[2:])
TARGET_NAMES = ('K1^3', 'Psi_integer', 's1^5', 'u000^3', 'K5^2', 'K6^3', 'K7^2')
WEIGHTS = {p: 2, r: 3, s0: 2, s1: 3, s2: 4,
           **{value: (-sum(indices)-2) % 12 for indices, value in U.items()}}
TARGET_WEIGHTS = (6, 6, 3, 6, 4, 6, 8)
TARGET_DEGREES = (3, 3, 5, 3, 2, 3, 4)


def character(expression):
    terms = S.Poly(expression, *VARIABLES).terms()
    weights = {sum(exponent*WEIGHTS[variable] for variable, exponent in zip(VARIABLES, powers)) % 12
               for powers, coefficient in terms if coefficient}
    assert len(weights) <= 1, 'inhomogeneous expression: '+str(expression)
    return next(iter(weights)) if weights else None


def module_program(circuit_sha256):
    generators = [('trace'+str(k), (-k) % 12) for k in range(3)]
    generators += [('mu'+str(n), n % 12) for n in range(2, 15)]
    columns = [{'generator': name, 's2_power': k, 'character_mod12': (weight+4*k) % 12}
               for name, weight in generators for k in range(3)]
    targets = []
    for name, expression, weight, degree in zip(TARGET_NAMES, TARGETS, TARGET_WEIGHTS, TARGET_DEGREES):
        reduced = vector(expression)
        for i, value in enumerate(reduced):
            assert character(value) in (None, (weight-4*i) % 12)
        targets.append({'name': name, 'polynomial': str(expression), 'degree': degree,
                        'character_mod12': weight, 'reduced_vector': list(map(str, reduced)),
                        'component_characters_mod12': [(weight-4*i) % 12 for i in range(3)],
                        'multiplier_characters_by_source_column_mod12':
                            [(weight-column['character_mod12']) % 12 for column in columns]})
    return {'schema': 1, 'coefficient_field': 'QQ',
            'ambient_variables': list(map(str, VARIABLES)),
            'base_variables': list(map(str, BASE_VARIABLES)),
            'base_variable_characters_mod12': {str(v): WEIGHTS[v] for v in BASE_VARIABLES},
            'basis': ['1', 's2', 's2^2'], 'basis_characters_mod12': [0, 4, 8],
            'source_circuit': 'cyclic_frame_moments_14.circuit.json.gz',
            'source_circuit_sha256': circuit_sha256,
            'source_column_recipes': columns, 'targets': targets,
            'power_gate_local_exponents': list(EXPONENTS),
            'Psi_primitive_integer_coefficients': list(PSI_COEFFICIENTS),
            'global_length_equivalent_to_all_memberships': 51,
            'memberships_verified': 0, 'witness_multiplier_degree_bound': None,
            'finite_ansatz_column_count': None, 'grading_reduction_factor': None,
            'scope': 'Seven equivalent target memberships, all raw degrees<=5, with exact '
                     'mod12 multiplier filters. The48 source-generator recipes remain; '
                     'actual finite ansatz columns must be counted for the chosen degree '
                     'bound and target. No witness-degree bound, membership or global '
                     'length51 proof is asserted.'}
