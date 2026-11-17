#!/usr/bin/env python3
"""Export seven separate graded radical tests; do not solve them.

The mathematical cover and generic-specialization obligations are recorded in
FIELD_COVER_AND_GENERIC_COMPONENT_STRATEGY.md. Prefix14 is available, but a
certificate may use any later finite prefix and arbitrary radical exponents.
"""
import hashlib
import json
from pathlib import Path

from seven_target_gate import S, K, WEIGHTS, VARIABLES, character
from cyclic_nilpotence_gate import BASE_VARIABLES, MONIC, MULTIPLICATION_S2, vector


def main():
    root = Path(__file__).resolve().parent
    v = S.Symbol('v')
    generators = [('trace'+str(k), (-k) % 12) for k in range(3)]
    generators += [('mu'+str(n), n % 12) for n in range(2, 15)]
    sources = [{'generator': name, 's2_power': k, 'character_mod12': (weight+4*k) % 12}
               for name, weight in generators+[('unit_relation', 0)] for k in range(3)]
    assert len(sources) == 51
    problems = []
    for i, function in enumerate(K, 1):
        weight = character(function)
        inverse_weight = (-weight) % 12
        unit_relation = S.expand(1-v*function)
        terms = S.Poly(unit_relation, *VARIABLES, v).terms()
        all_weights = [WEIGHTS[variable] for variable in VARIABLES]+[inverse_weight]
        assert {sum(exponent*grade for exponent, grade in zip(powers, all_weights)) % 12
                for powers, coefficient in terms if coefficient} == {0}
        reduced = vector(unit_relation)
        assert sum(reduced[k]*S.Symbol('s2')**k for k in range(3)) == unit_relation
        problems.append({'name': 'invert_K'+str(i), 'function': str(function),
                         'function_degree': S.Poly(function, *VARIABLES).total_degree(),
                         'function_character_mod12': weight,
                         'auxiliary_variable': str(v), 'auxiliary_character_mod12': inverse_weight,
                         'unit_relation': str(unit_relation),
                         'unit_relation_remainder_vector': list(map(str, reduced)),
                         'unit_target_vector': ['1', '0', '0'],
                         'source_column_multiplier_characters_mod12':
                             [(-column['character_mod12']) % 12 for column in sources],
                         'verified': False, 'certificate_prefix': None,
                         'witness_degree_bound': None})
    assert [problem['function_degree'] for problem in problems] == [1, 1, 1, 1, 1, 1, 2]
    old_program = json.loads((root/'cyclic_prefix14_seven_target_module_program.json').read_text())
    circuit = root/old_program['source_circuit']
    assert hashlib.sha256(circuit.read_bytes()).hexdigest() == old_program['source_circuit_sha256']
    program = {'schema': 1, 'coefficient_field': 'QQ',
               'purpose': 'Rodrigues-only reduced coverage option, not a necessary condition '
                          'for Image when other pure components are safe.',
               'available_prefix': 14, 'later_prefixes_allowed': True,
               'moment_extension_performed': False,
               'source_circuit': circuit.name,
               'source_circuit_sha256': old_program['source_circuit_sha256'],
               'base_variables_per_problem': list(map(str, BASE_VARIABLES))+[str(v)],
               'base_variable_characters_mod12': {str(var): WEIGHTS[var] for var in BASE_VARIABLES},
               'Gram_basis': ['1', 's2', 's2^2'], 'Gram_basis_characters_mod12': [0, 4, 8],
               'monic_Gram_relation': str(MONIC),
               'multiplication_s2': [[str(value) for value in row]
                                     for row in MULTIPLICATION_S2.tolist()],
               'source_column_recipes_at_prefix14': sources,
               'source_recipe_count_at_prefix_N': '3*(N+3), after removing monic Gram',
               'problems': problems, 'tests_are_separate': True,
               'unit_certificates_verified': 0, 'radical_exponents_prescribed': False,
               'finite_witness_column_counts': None,
               'scope': 'Seven separate 1 in I_N+(1-v*K_i) tests; no unit or radical '
                        'membership is proved. At fixed N>=14 Rodrigues-only coverage '
                        'is equivalent to length51 because the local lengths are known. '
                        'The51 source recipes at prefix14 are bookkeeping, not an '
                        'exhaustion-length computation. Generic safety of additional '
                        'pure components and D=0 coverage are separate routes/obligations.'}
    (root/'cyclic_field_radical_gate_program.json').write_text(json.dumps(program, indent=2)+'\n')
    print(json.dumps({'input_export': 'PASS', 'separate_problems': len(problems),
                      'function_degrees': [problem['function_degree'] for problem in problems],
                      'later_prefixes_allowed': True, 'verified_unit_certificates': 0}, indent=2))


if __name__ == '__main__':
    main()
