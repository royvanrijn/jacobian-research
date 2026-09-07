#!/usr/bin/env python3
"""Scalar gluing and simultaneous elliptic norm obstructions, retrospectively."""
import argparse
from itertools import product
from math import isqrt
from pathlib import Path
import retrospective as r
import production_twist_blocks as blocks
import scalar_cup

PROTOCOL = Path(__file__).with_name('QUADRATIC_NORM_BLOCK_PROTOCOL.json')
MODELS = r.OUT/'rank_jump_production_minus_twist_inputs_v1.json'
INPUT = r.OUT/'rank_jump_quadratic_norm_block_inputs_v1.json'
OUTPUT = r.OUT/'rank_jump_quadratic_norm_blocks_v1.json'


def mm(A, B):
    return tuple(sum(A[2*i+k]*B[2*k+j] for k in range(2)) % 2
                 for i in range(2) for j in range(2))


def centralizer(generators):
    gl = [M for M in product(range(2), repeat=4) if (M[0]*M[3]-M[1]*M[2]) % 2]
    return [list(M) for M in gl if all(mm(M, g) == mm(g, M) for g in generators)]


def export():
    models = {x['case_index']: x for x in r.read(MODELS)['cases']}
    cups = {x['case_index']: x for x in r.read(scalar_cup.OUTPUT)['production_cases']}
    rows = []
    for b in r.read(blocks.OUTPUT)['rows']:
        i = b['case_index']
        rows.append({
            'case_index': i, 'id': b['id'],
            'cubic_ascending': models[i]['integral_cubic_ascending'],
            'generic_rank': b['generic_rank'], 'known_independent_rank': b['known_independent_rank'],
            'observed_quotient_rank': b['observed_quotient_rank'],
            'generic_strict_dimension': b['generic_strict_dimension'],
            'CT_matrix': cups[i]['scalar_cup_matrix'],
            'hyperbolic_pairs': b['hyperbolic_pairs'], 'radical_basis': b['radical_basis']
        })
    r.write_new(INPUT, {
        'schema': 'rank-jump.quadratic-norm-block-inputs.v1',
        'source_hashes': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes())
                          for p in (MODELS, scalar_cup.OUTPUT, blocks.OUTPUT)},
        'rows': rows, 'small_control_cubic_ascending': [-1, -14, -11, 1],
        'boundary': 'Retrospective CT matrices and witness masks, not prospective features.'
    })


def cubic_certificate(coeff, bound):
    d, c, b, a = map(int, coeff)
    assert a == 1
    disc = b*b*c*c-4*a*c**3-4*b**3*d-27*a*a*d*d+18*a*b*c*d
    assert disc
    square = disc > 0 and isqrt(disc)**2 == disc
    witness = next((p for p in r.primes(bound)
                    if disc % p and all((x**3+b*x*x+c*x+d) % p for x in range(p))), None)
    assert witness is not None, 'UNKNOWN: bounded irreducibility witness absent'
    return {'discriminant': str(disc), 'discriminant_square': square,
            'irreducible_reduction_prime': witness,
            'group': 'C3' if square else 'S3'}


def compute():
    inp = r.read(INPUT)
    for name, digest in inp['source_hashes'].items():
        assert r.digest((r.ROOT/name).read_bytes()) == digest
    cycle, swap = (0, 1, 1, 1), (0, 1, 1, 0)
    s3, c3 = centralizer([cycle, swap]), centralizer([cycle])
    assert s3 == [[1, 0, 0, 1]] and len(c3) == 3
    bound = r.read(PROTOCOL)['limits']['maximum_modular_prime']
    rows = []
    for old in inp['rows']:
        cert = cubic_certificate(old['cubic_ascending'], bound)
        assert cert['group'] == 'S3'
        M = list(map(r.pack, old['CT_matrix']))
        n = len(M)
        pairs = [[v['strict_basis_mask'] for v in pair] for pair in old['hyperbolic_pairs']]
        flat = [v for pair in pairs for v in pair]
        radical = [v['strict_basis_mask'] for v in old['radical_basis']]
        assert r.rank(flat+radical) == n
        for i, x in enumerate(flat+radical):
            for j, y in enumerate(flat+radical):
                expected = int(i < len(flat) and j < len(flat) and i//2 == j//2 and i != j)
                assert blocks.pairing(M, x, y) == expected
        obstruction = len(flat)
        assert r.rank(M) == obstruction
        rows.append({
            'id': old['id'], 'case_index': old['case_index'], 'cubic_certificate': cert,
            'generic_rank': old['generic_rank'], 'known_independent_rank': old['known_independent_rank'],
            'observed_quotient_rank': old['observed_quotient_rank'],
            'retained_strict_dimension': n, 'generic_strict_dimension': old['generic_strict_dimension'],
            'nondegenerate_block_masks': flat, 'CT_rank': obstruction,
            'local_global_elliptic_norm_defect_dimension_at_least': obstruction,
            'local_global_elliptic_norm_defect_order_at_least': 2**obstruction,
            'retained_norm_kernel_dimension_at_most': n-obstruction,
            'independent_twist_Sha_classes_killed_by_Qi_at_least': obstruction,
            'field': 'Q(i)', 'nonzero_selected_Sha_classes_period_and_index': 2,
            'rational_2_torsion_gluings_to_any_quadratic_twist': 1,
            'smooth_genus_two_Jacobian_from_this_2_gluing': False,
            'full_norm_defect_dimension': 'UNKNOWN'
        })
    small = cubic_certificate(inp['small_control_cubic_ascending'], bound)
    assert small['group'] == 'C3'
    return {
        'schema': 'rank-jump.quadratic-norm-blocks.v1',
        'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes())
                     for p in (Path(__file__), PROTOCOL, INPUT, Path(blocks.__file__))},
        'status': 'PASS', 'GL2_centralizers': {'S3': s3, 'C3': c3},
        'small_control': small, 'rows': rows,
        'geometric_object': 'Res_(Q(i)/Q)(E over Q(i)), with its geometrically product principal polarization',
        'norm_defect_identity': 'E(Q)_locally_norm / N E(Q(i)) = ker(Sha(E^(-1)/Q) -> Sha(E^(-1)/Q(i)))',
        'boundary': 'Inherited CT evidence proves lower bounds for one simultaneous solubility obstruction. The norm index is finite and does not add free rank. Original ranks, twist ranks, and full norm defects remain UNKNOWN.'
    }


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=['export', 'build', 'check'])
    args = p.parse_args()
    if args.mode == 'export':
        export()
    else:
        data = compute()
        if args.mode == 'check':
            assert r.read(OUTPUT) == data
            print('PASS quadratic norm blocks and scalar-gluing gate')
        else:
            r.write_new(OUTPUT, data)
            print(data['status'], [(x['case_index'], x['CT_rank']) for x in data['rows']])
