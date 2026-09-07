#!/usr/bin/env python3
"""Finite-presentation regressions; never certifies an arithmetic presentation."""
import argparse
import itertools
from pathlib import Path
import retrospective as r

OUT = r.OUT / 'rank_jump_class_creation_relation_gate_v1.json'
CAP = r.OUT / 'rank_jump_strict_cover_creation_capacity_v1.json'
IDEALS = r.OUT / 'rank_jump_generic_sunit_carrier_v1.json'
PROTOCOL = Path(__file__).with_name('CLASS_CREATION_RELATION_GATE_PROTOCOL.json')


def invariants(R, A):
    from sage.all import GF, ZZ
    assert R.base_ring() == ZZ and R.rank() == R.nrows()
    assert A.nrows() == R.nrows()
    D, L, V = R.smith_form()
    assert L * R * V == D and abs(L.det()) == abs(V.det()) == 1
    orders = [abs(int(D[i, i])) for i in range(R.nrows())]
    B = L * A
    assert all(int(2 * B[i, j]) % orders[i] == 0
               for i in range(R.nrows()) for j in range(A.ncols()))
    c = sum(d % 2 == 0 for d in orders)
    rho4 = sum(d % 4 == 0 for d in orders)
    detected = R.augment(A).change_ring(GF(2)).rank() - R.change_ring(GF(2)).rank()
    return {'quadratic_character_dimension': c, 'four_rank': rho4,
            'elementary_factor_dimension': c - rho4,
            'generic_ideal_image_mod_2_dimension': int(detected)}


def calculate():
    from sage.all import ZZ, matrix, diagonal_matrix
    checks = 0
    # Independent enumeration of Hom(C,Z/2), Hom(C,Z/4), and characters
    # restricted to C[2]; explicitly abstract groups, not new number fields.
    for ds in itertools.product([1, 2, 3, 4, 6, 8, 12], repeat=2):
        R = diagonal_matrix(ZZ, ds)
        A = diagonal_matrix(ZZ, [d // 2 if d % 2 == 0 else 0 for d in ds])
        got = invariants(R, A)
        counts = [sum(all(d * x % q == 0 for d, x in zip(ds, v))
                      for v in itertools.product(range(q), repeat=2)) for q in (2, 4)]
        assert counts[0] == 2 ** got['quadratic_character_dimension']
        assert counts[1] == 2 ** (got['quadratic_character_dimension'] + got['four_rank'])
        restrictions = {tuple((x * int(A[i, i])) % 2 for i, x in enumerate(v))
                        for v in itertools.product(range(2), repeat=2)
                        if all(d * x % 2 == 0 for d, x in zip(ds, v))}
        assert len(restrictions) == 2 ** got['elementary_factor_dimension']
        assert got['generic_ideal_image_mod_2_dimension'] == got['elementary_factor_dimension']
        U = matrix(ZZ, [[1, 3], [0, 1]])
        V = matrix(ZZ, [[1, 0], [2, 1]])
        assert got == invariants(U * R * V, U * A)
        checks += 1
    models = []
    for label, ds, av in [('ten_C4', [4]*10, [2]*10),
                          ('ten_C2', [2]*10, [1]*10),
                          ('missing_relation', [2], [1]),
                          ('relation_completed', [1], [0])]:
        models.append({'label': label, 'provenance': 'ABSTRACT_REGRESSION_ONLY',
                       **invariants(diagonal_matrix(ZZ, ds), diagonal_matrix(ZZ, av))})
    assert models[0]['quadratic_character_dimension'] == models[1]['quadratic_character_dimension'] == 10
    assert models[0]['generic_ideal_image_mod_2_dimension'] == 0
    assert models[1]['generic_ideal_image_mod_2_dimension'] == 10
    panel = []
    for row in r.read(CAP)['rows']:
        panel.append({'token': row['token'], 'family': row['family'],
            'oracle_rank_derived_strict_dimension_lower_bound': row.get('total_strict_rational_dimension_lower_bound', 'UNKNOWN'),
            'generic_strict_dimension': row.get('generic_strict_dimension', 'UNKNOWN'),
            'equation_only_class_two_rank': 'UNKNOWN',
            'equation_only_class_four_rank': 'UNKNOWN',
            'certified_two_primary_relation_presentation': 'ABSENT_IN_THIS_INPUT',
            'paired_absolute_class_creation': 'UNKNOWN'})
    retro = [{'id': row['id'], 'case_index': row['case_index'],
              'elementary_factor_dimension_lower_bound': row['generic_detected_rank'],
              'provenance': 'RETROSPECTIVE_ORACLE_CHARACTERS',
              'absolute_two_rank': 'UNKNOWN', 'absolute_four_rank': 'UNKNOWN'}
             for row in r.read(IDEALS)['rows']]
    return {'schema': 'rank-jump.class-creation-relation-gate.v1',
        'status': 'PASS_ABSTRACT_REGRESSIONS_AND_SOURCE_AUDIT',
        'exhaustive_two_factor_cases': checks, 'abstract_models': models,
        'fresh_panel': panel, 'earlier_ideal_certificates': retro,
        'new_arithmetic_classes_constructed': 0,
        'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes())
                     for p in [Path(__file__), PROTOCOL, CAP, IDEALS]}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['capture', 'check'])
    args = parser.parse_args()
    result = calculate()
    if args.mode == 'capture': r.write_new(OUT, result)
    else: assert r.read(OUT) == result
    print(result['status'], result['exhaustive_two_factor_cases'])
