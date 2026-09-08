#!/usr/bin/env python3
"""Refine the class-rank target using proved2-Selmer parity over Q."""
import argparse
import json
from math import prod
from pathlib import Path
import subprocess

import certify_small_conductor_class_rank_target as base

ROOT, ART = base.ROOT, base.ART
OUTPUT = ART / 'small_conductor_selmer_rank_target_v1.json'
GP = Path('/usr/bin/gp')


def expected():
    original = base.expected()
    if json.loads(json.dumps(original)) != base.original.cert.read(base.OUTPUT):
        raise ArithmeticError('base Brumer--Kramer proof differs')
    proof = base.original.cert.read(base.PROOF)
    primes = [0] + [int(p) for p, _ in proof['discriminant_factorization']]
    program = ('E=ellinit([' + ','.join(proof['integral_model']) + ']);print(['
               + ','.join(f'ellrootno(E,{p})' for p in primes) + ']);quit\n')
    process = subprocess.run([str(GP), '-q', '-f'], input=program, text=True,
                             capture_output=True, timeout=10, check=True)
    if process.stderr:
        raise ArithmeticError(process.stderr)
    signs = json.loads(process.stdout)
    if len(signs) != len(primes) or any(s not in (-1, 1) for s in signs) or prod(signs) != 1:
        raise ArithmeticError('fixed global sign+1 required')
    offset = original['brumer_kramer']['offset']
    rows = [{'class_two_rank_upper_bound': g, 'raw_selmer_upper_bound': g+offset,
             'parity_refined_selmer_and_rank_upper_bound': 2*((g+offset)//2)}
            for g in range(15, 21)]
    return {
        'schema': 'elliptic-curves.small-conductor-selmer-parity-target.v1',
        'status': 'PASS', 'sources': {str(p.relative_to(ROOT)): base.original.cert.hashed(p)
                                    for p in [Path(__file__).resolve(), base.OUTPUT]},
        'gp_binary': str(GP), 'gp_binary_sha256': base.original.cert.hashed(GP),
        'root_number_program': program, 'root_number_stdout': process.stdout,
        'local_root_numbers': dict(zip(map(str, primes), signs)), 'global_root_number': 1,
        'theorem': 'For elliptic curves overQ, the parity of dim_F2 Sel_2(E/Q)-dim_F2 E(Q)[2] equals the root number. This follows from the proved2-parity theorem for the2-infinity Selmer corank and the alternating Cassels--Tate pairing on the finite quotient. It does not assert algebraic rank parity.',
        'references': [
            'https://annals.math.princeton.edu/2010/172-1/p11',
            'https://arxiv.org/html/1606.07178#A1'],
        'rational_two_torsion_dimension': 0, 'selmer_dimension_parity': 0,
        'brumer_kramer_offset': offset,
        'unconditional_class_two_rank_lower_bound': 15,
        'upper_bound_implications': rows,
        'sufficient_class_two_rank_upper_bound_for_exact_rank22': 16,
        'necessary_class_two_rank_lower_bound_for_rank_at_least23': 17,
        'current_class_two_rank_upper_bound': None, 'current_rank_upper_bound': None,
        'argument': 'An upper bound g<=16 gives Sel_2 dimension<=23. Its proved even parity improves this to22, matching the certified22 points. Conversely rank>=23 requires even Sel_2 dimension>=24, forcing g>=17.',
        'claim_boundary': 'An exact theorem-directed target, not a completed class-group bound or descent. No GRH, BSD or finiteness of Sha is assumed for this implication. If a future class-group bound depends on GRH, the resulting rank statement inherits that assumption.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    data = expected()
    if args.check:
        if json.loads(json.dumps(data)) != base.original.cert.read(OUTPUT):
            raise ArithmeticError('Selmer-parity target differs')
    else:
        if OUTPUT.exists():
            raise FileExistsError('preserve parity target')
        base.original.cert.write(OUTPUT, data)
    print('PROVED IMPLICATION: class2-rank <=16 suffices for exact rank22; rank>=23 requires class2-rank>=17')
