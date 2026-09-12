#!/usr/bin/env sage-python
"""BNF-free cubic principal relations with exact witnesses and mod-two algebra.

Small KSW-style smooth-norm prototype, not the production NFS/Julia-reduction
implementation. No global factor-base generation theorem is asserted.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
from math import prod
from pathlib import Path
import time

from sage.all import QQ, ZZ, PolynomialRing, pari, prime_range, RealIntervalField
import two_class_relation_core as core
import wide_arithmetic_profile_core as storage


def encode_element(nf, alpha):
    return [str(v) for v in nf.nfalgtobasis(alpha)]


def index_form(nf, f):
    """Delone-Faddeev multiplication-table form of the maximal cubic order."""
    zk = nf.nf_get_zk()
    assert zk[0] == 1 and len(zk) == 3
    product_coordinates = nf.nfalgtobasis(zk[1]*zk[2])
    omega = pari.Mod(zk[1]-product_coordinates[2], f)
    theta = pari.Mod(zk[2]-product_coordinates[1], f)
    w2, t2 = nf.nfalgtobasis(omega**2), nf.nfalgtobasis(theta**2)
    a, b, c, d = -ZZ(w2[2]), ZZ(w2[1]), -ZZ(t2[2]), ZZ(t2[1])
    assert omega*theta == -a*d
    assert omega**2 == -a*c+b*omega-a*theta
    assert theta**2 == -b*d+d*omega-c*theta
    initial = [int(d), int(c), int(b), int(a)]
    R = PolynomialRing(QQ, 'z')
    assert R(initial).discriminant() == nf.disc()
    alpha = -omega/a
    assert sum(initial[i]*alpha**i for i in range(4)) == 0
    reduced, matrix, steps, capped = core.reduce_form(initial)
    p, q, r, s = matrix
    beta = (s*alpha-q)/(p-r*alpha)
    assert sum(reduced[i]*beta**i for i in range(4)) == 0
    assert R(reduced).discriminant() == nf.disc()
    assert all(QQ(v).denominator() == 1 for v in nf.nfalgtobasis(reduced[3]*beta))
    return beta, {'initial_cubic_ascending': initial, 'reduced_cubic_ascending': reduced,
                  'GL2_matrix': matrix, 'height_descent_steps': steps, 'height_descent_capped': capped,
                  'reduction_method': 'finite exact height descent, not certified Julia reduction',
                  'generator_in_integral_basis': encode_element(nf, beta),
                  'integral_basis_polynomials': [str(v) for v in nf.nf_get_zk()],
                  'form_discriminant_equals_field_discriminant': True}


def collect(input_path, output):
    start = time.monotonic()
    inp = json.loads(input_path.read_text())
    protocol = inp['protocol']
    R = PolynomialRing(QQ, 'x')
    f = pari(R(list(map(QQ, inp['integral_monic_cubic_ascending']))))
    hints = list(map(ZZ, inp['certified_bad_rational_primes']))
    assert all(p.is_prime(proof=True) for p in hints)
    pari.addprimes(hints)
    nf = pari.nfinit([f, sorted(set([ZZ(2)]+hints))])
    assert pari.nfcertify(nf) == []
    assert str(nf.disc()) == inp['field_discriminant']
    assert list(map(int, nf.nf_get_sign())) == inp['field_signature']
    beta, form = index_form(nf, f)
    c = form['reduced_cubic_ascending']
    denominator_factors = [(int(p), int(e)) for p, e in ZZ(c[3]).factor(proof=True)]
    assert prod(p**e for p, e in denominator_factors) == c[3]
    rational_support = sorted(set(map(int, prime_range(protocol['factor_base_bound']+1))) |
                              {p for p, e in denominator_factors})
    primorial = prod(rational_support)
    ideals, factor_base, above = [], [], {}
    for p in rational_support:
        indices = []
        for P in sorted(nf.idealprimedec(p), key=lambda P: str(nf.idealhnf(P))):
            i = len(ideals)
            ideals.append(P)
            indices.append(i)
            factor_base.append({'column': i, 'rational_prime': p, 'ramification_index': int(P[2]),
                'residue_degree': int(P[3]), 'norm': str(ZZ(p)**int(P[3])),
                'hnf': [[str(v) for v in row] for row in nf.idealhnf(P).sage().rows()],
                'mandatory_denominator_prime': p in {q for q, e in denominator_factors}})
        assert sum(factor_base[i]['ramification_index']*factor_base[i]['residue_degree'] for i in indices) == 3
        above[p] = indices
    relations, rows = [], []
    for p in rational_support:
        exponents = [[i, factor_base[i]['ramification_index']] for i in above[p]]
        row = sum((e % 2)*(1 << i) for i, e in exponents)
        rows.append(row)
        relations.append({'kind': 'canonical_rational_prime', 'rational_prime': p,
                          'generator_in_integral_basis': [''+str(p), '0', '0'],
                          'signed_ideal_exponents': exponents, 'parity_hex': hex(row)})
    canonical_count = len(rows)
    # Deliberately include all residue degrees above the rational support.
    # Canonical (p) rows effect the cubic degree-2/3 eliminations mod 2.
    bach = 12*RealIntervalField(128)(abs(ZZ(nf.disc()))).log()**2
    setup = {'schema': 'elliptic-curves.two-class-relation-setup.v1', 'status': 'PASS_CERTIFIED_ORDER_AND_FACTOR_BASE',
        'input_sha256': storage.sha256_file(input_path), 'form': form,
        'denominator_factorization': denominator_factors, 'factor_base': factor_base,
        'rational_support': rational_support, 'canonical_count': canonical_count,
        'baseline': core.matrix_stats(rows, len(ideals), canonical_count),
        'bach_GRH_bound_interval': str(bach), 'factor_base_generation': 'UNKNOWN_NOT_CERTIFIED',
        'global_class_2rank_upper_bound': None, 'global_class_2rank_estimate': None,
        'setup_wall_seconds': round(time.monotonic()-start, 6)}
    storage.write_immutable(output/'setup.json', setup)
    print(f"MOD2_SETUP|key={inp['curve_key']}|columns={len(ideals)}|canonical_rank={canonical_count}|deficiency={len(ideals)-canonical_count}", flush=True)
    checkpoints, smooth, norm_bits, residual_bits = [], 0, Counter(), Counter()
    trial_hash = sha256()
    for attempted, (a, b) in enumerate(core.primitive_pairs(protocol['candidate_budget']), 1):
        numerator = core.form_value(c, a, b)
        residual = core.strip_support(numerator, primorial)
        norm_bits[abs(numerator).bit_length()] += 1
        residual_bits[residual.bit_length()] += 1
        trial_hash.update(f'{a},{b},{numerator},{residual}\n'.encode())
        if residual == 1:
            gamma = a+b*beta
            assert QQ(nf.nfeltnorm(gamma)) == QQ(numerator)/c[3]
            exponents = []
            norm = QQ(1)
            # Since a*beta is integral and its norm has known rational
            # support, no unrecorded ideal can divide gamma. Signed
            # valuations are checked individually and by ideal identity.
            ideal_product = nf.idealhnf(1)
            for p in rational_support:
                if numerator % p and c[3] % p:
                    continue
                for i in above[p]:
                    e = int(nf.idealval(gamma, ideals[i]))
                    if e:
                        exponents.append([i, e])
                        norm *= QQ(p)**(factor_base[i]['residue_degree']*e)
                        ideal_product = nf.idealmul(ideal_product, nf.idealpow(ideals[i], e))
            assert norm == abs(QQ(nf.nfeltnorm(gamma)))
            assert nf.idealhnf(ideal_product) == nf.idealhnf(gamma)
            row = sum((e % 2)*(1 << i) for i, e in exponents)
            rows.append(row)
            smooth += 1
            relations.append({'kind': 'smooth_binary_cubic', 'attempt': attempted, 'a': a, 'b': b,
                'norm': str(QQ(numerator)/c[3]), 'generator_in_integral_basis': encode_element(nf, gamma),
                'signed_ideal_exponents': exponents, 'parity_hex': hex(row)})
        if attempted in protocol['candidate_prefixes'] or attempted == protocol['candidate_budget']:
            point = {'attempted': attempted, 'smooth_relations': smooth,
                     **core.matrix_stats(rows, len(ideals), canonical_count),
                     'trial_sha256': trial_hash.hexdigest()}
            checkpoints.append(point)
            storage.write_immutable(output/f'prefix-{attempted:06}.json', point)
            storage.write_immutable(output/f'relations-{attempted:06}.json', {'rows': relations})
            print(f"MOD2_PREFIX|key={inp['curve_key']}|attempts={attempted}|smooth={smooth}|deficiency={point['deficiency']}", flush=True)
    result = {'schema': 'elliptic-curves.two-class-relation-result.v1',
        'status': 'PASS_BOUNDED_RELATION_PROTOTYPE_NOT_CLASS_RANK', 'curve_key': inp['curve_key'],
        'input_sha256': storage.sha256_file(input_path), 'setup_sha256': storage.sha256_file(output/'setup.json'),
        'candidate_prefixes': checkpoints, 'norm_numerator_bit_histogram': dict(sorted(norm_bits.items())),
        'residual_bit_histogram': dict(sorted(residual_bits.items())), 'trial_sha256': trial_hash.hexdigest(),
        'final_matrix': core.matrix_stats(rows, len(ideals), canonical_count),
        'canonical_relation_count': canonical_count, 'noncanonical_relation_count': smooth,
        'factor_base_generation': 'UNKNOWN_NOT_CERTIFIED', 'global_class_2rank_upper_bound': None,
        'global_class_2rank_estimate': None, 'full_class_group_computed': False,
        'measurement_gate': 'NO_NONCANONICAL_RELATIONS' if smooth == 0 else 'PARTIAL_RELATIONS_ONLY',
        'boundary': 'Deficiency bounds only the image of the chosen factor base in Cl(K)/2Cl(K). It is neither an estimate nor an upper/lower bound on the full class 2-rank without a generation argument. Untested principal relations can lower the deficiency.',
        'wall_seconds': round(time.monotonic()-start, 6)}
    storage.write_immutable(output/'result.json', result)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    collect(args.input, args.output)
