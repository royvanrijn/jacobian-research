#!/usr/bin/env sage-python
"""One fixed-budget extension of the frozen two-field relation pilot.

The worker receives only the prior equation input/setup, never MW coordinates.
Every relation has an exact principal-ideal witness and independent GF(2) rank.
"""
import argparse
from collections import Counter
from hashlib import sha256
from math import prod
from pathlib import Path
import time
from sage.all import QQ, ZZ, GF, PolynomialRing, matrix, pari

import two_class_skew_core as skew
import two_class_relation_core as core
import wide_arithmetic_profile_core as storage


def run(input_path, setup_path, output):
    start = time.monotonic()
    import json
    inp, setup = [json.loads(p.read_text()) for p in [input_path, setup_path]]
    c = setup['form']['reduced_cubic_ascending']
    box = skew.rectangle(c, 18)
    protocol = {'input_sha256': storage.sha256_file(input_path), 'setup_sha256': storage.sha256_file(setup_path),
                'sources': {str(p): storage.sha256_file(p) for p in [Path(__file__), Path(skew.__file__), Path(core.__file__)]},
                'rectangle': box, 'selection': 'minimum exact coefficient triangle bound among A*B=2^18 dyadic boxes',
                'rational_support_unchanged': True, 'max_wall_seconds': 180,
                'candidate_rule': 'all primitive pairs in [-A,A] x [1,B], ordered by short axis then long axis',
                'filter': 'exact modular-root progression division by every support prime, including all powers; no score cutoff',
                'new_point_searches': False, 'full_class_group': False,
                'serious_density_gate': 'at least100 independent noncanonical rows and at least5 percent of initial deficiency removed',
                'global_class_2rank_estimate': None, 'factor_base_generation': 'UNKNOWN_NOT_CERTIFIED'}
    storage.write_immutable(output/'plan.json', protocol)
    R = PolynomialRing(QQ, 'x')
    f = pari(R(list(map(QQ, inp['integral_monic_cubic_ascending']))))
    hints = list(map(ZZ, inp['certified_bad_rational_primes']))
    assert all(p.is_prime(proof=True) for p in hints)
    pari.addprimes(hints)
    nf = pari.nfinit([f, sorted(set([ZZ(2)]+hints))])
    assert pari.nfcertify(nf) == [] and str(nf.disc()) == inp['field_discriminant']
    assert [str(v) for v in nf.nf_get_zk()] == setup['form']['integral_basis_polynomials']
    beta = nf.nfbasistoalg(pari(list(map(QQ, setup['form']['generator_in_integral_basis']))).Col())
    assert sum(c[i]*beta**i for i in range(4)) == 0
    assert str(R(c).discriminant()) == inp['field_discriminant']
    assert all(QQ(v).denominator() == 1 for v in nf.nfalgtobasis(c[3]*beta))
    fb, support = setup['factor_base'], setup['rational_support']
    ideals, above, masks = [], {}, []
    for p in support:
        above[p] = []
        for P in sorted(nf.idealprimedec(p), key=lambda P: str(nf.idealhnf(P))):
            i = len(ideals)
            assert fb[i]['column'] == i and fb[i]['rational_prime'] == p
            assert fb[i]['hnf'] == [[str(v) for v in r] for r in nf.idealhnf(P).sage().rows()]
            assert (fb[i]['ramification_index'], fb[i]['residue_degree']) == (int(P[2]), int(P[3]))
            above[p].append(i)
            ideals.append(P)
        masks.append(sum((fb[i]['ramification_index'] % 2)*(1 << i) for i in above[p]))
    assert len(ideals) == len(fb)
    canonical = len(masks)
    fa, fbeta = R([-c[0], c[1], -c[2], c[3]]), R([c[3], -c[2], c[1], -c[0]])
    roots = []
    for p in support:
        ra, rb = [sorted(map(int, pari.polrootsmod(poly, p))) for poly in [fa, fbeta]]
        assert all(int(fa(r)) % p == 0 for r in ra) and all(int(fbeta(r)) % p == 0 for r in rb)
        roots.append((p, ra, rb))
    storage.write_immutable(output/'roots.json', roots)
    print('SKEW_SETUP', inp['curve_key'], box, flush=True)
    relations, bits, residual_bits = [], Counter(), Counter()
    digest = sha256()
    attempted = 0
    for line, values in enumerate(skew.sieve_lines(c, box, roots), 1):
        for a, b, numerator, residual in values:
            attempted += 1
            digest.update(f'{a},{b},{numerator},{residual}\n'.encode())
            bits[abs(numerator).bit_length()] += 1
            residual_bits[residual.bit_length()] += 1
            if residual != 1:
                continue
            gamma = a+b*beta
            assert QQ(nf.nfeltnorm(gamma)) == QQ(numerator)/c[3]
            exponents, product = [], nf.idealhnf(1)
            for p in support:
                if numerator % p and c[3] % p:
                    continue
                for i in above[p]:
                    e = int(nf.idealval(gamma, ideals[i]))
                    if e:
                        exponents.append([i, e])
                        product = nf.idealmul(product, nf.idealpow(ideals[i], e))
            assert nf.idealhnf(product) == nf.idealhnf(gamma)
            mask = sum((e % 2)*(1 << i) for i, e in exponents)
            masks.append(mask)
            relations.append({'a': a, 'b': b, 'attempt': attempted, 'norm': str(QQ(numerator)/c[3]),
                              'generator_in_integral_basis': [str(v) for v in nf.nfalgtobasis(gamma)],
                              'signed_ideal_exponents': exponents, 'parity_hex': hex(mask)})
        storage.write_immutable(output/f'line-{line:05}.json', {'attempted': attempted, 'trial_sha256': digest.hexdigest(),
                                'relations': relations, 'matrix': core.matrix_stats(masks, len(fb), canonical)})
    # Independent complete smoothness replay: a primorial-GCD algorithm, not
    # modular root completeness or its progression implementation.
    rehash = sha256()
    smooth = []
    primorial = prod(support)
    A, B = box['A'], box['B']
    pairs = ((a, b) for b in range(1, B+1) for a in range(-A, A+1)) if A >= B else (
             (a, b) for a in range(-A, A+1) for b in range(1, B+1))
    from math import gcd
    replayed = 0
    for a, b in pairs:
        if gcd(a, b) != 1:
            continue
        n = sum(c[i]*a**i*(-b)**(3-i) for i in range(4))
        residual = abs(n)
        while True:
            common = gcd(residual, primorial)
            if common == 1:
                break
            residual //= common
        rehash.update(f'{a},{b},{n},{residual}\n'.encode())
        replayed += 1
        if residual == 1:
            smooth.append((a, b))
    assert replayed == attempted and rehash.hexdigest() == digest.hexdigest()
    assert smooth == [(r['a'], r['b']) for r in relations]
    M = matrix(GF(2), len(masks), len(fb), {(j, i): 1 for j, v in enumerate(masks)
               for i in range(len(fb)) if (v >> i) & 1}, sparse=False)
    stats = core.matrix_stats(masks, len(fb), canonical)
    assert M.rank() == stats['rank_mod2']
    gate = stats['rank_gain_beyond_canonical'] >= 100 and 20*stats['rank_gain_beyond_canonical'] >= len(fb)-canonical
    result = {'status': 'PASS_EXACT_BOUNDED_SKEW_SIEVE_NOT_CLASS_RANK', 'attempted': attempted,
              'rectangle': box, 'noncanonical_relations': len(relations), 'relations': relations,
              'final_matrix': stats, 'trial_sha256': digest.hexdigest(),
              'independent_primorial_candidate_replay': 'PASS', 'independent_sage_GF2_rank': int(M.rank()),
              'norm_numerator_bit_histogram': dict(sorted(bits.items())), 'residual_bit_histogram': dict(sorted(residual_bits.items())),
              'density_gate': 'PASS_ENGINEERING_GATE_NOT_AUTHORIZATION' if gate else 'NOT_SERIOUS_RELATION_DENSITY',
              'factor_base_generation': 'UNKNOWN_NOT_CERTIFIED', 'global_class_2rank_estimate': None,
              'global_class_2rank_upper_bound': None, 'wall_seconds': round(time.monotonic()-start, 6),
              'bindings': {p.name: storage.sha256_file(p) for p in sorted(output.iterdir()) if p.is_file()}}
    storage.write_immutable(output/'result.json', result)
    print('SKEW_DONE', inp['curve_key'], 'attempted', attempted, 'relations', len(relations),
          'deficiency', stats['deficiency'], result['density_gate'], flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--setup', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    run(args.input, args.setup, args.output)
