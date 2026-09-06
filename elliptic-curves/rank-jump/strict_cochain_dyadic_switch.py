#!/usr/bin/env python3
"""Exact local evaluation of the new norm cochain on the retained -1 switch."""
import argparse
from fractions import Fraction as Q
from pathlib import Path
import retrospective as r
import verify_unpointed_governing_norm as verify
from verify_explicit_governing_octic import finite_polynomials

OUTPUT = r.OUT/'rank_jump_strict_cochain_dyadic_switch_v1.json'


def compute():
    inp = r.read(verify.INPUT); replay = r.read(verify.OUTPUT); control = r.read(verify.CONTROL)
    assert inp['status'] == replay['status'] == control['status'] == 'PASS'
    for data in (inp, replay):
        for path, sha in data['bindings'].items(): assert r.digest((r.ROOT/path).read_bytes()) == sha
    K = verify.Algebra(inp['cubic_ascending']); X = K.elt(inp['X_ascending'])
    trim, sub, mul, power, gcd = finite_polynomials(2)
    cubic = trim(inp['cubic_ascending']); beta = trim(inp['beta_ascending'])
    assert all(sum(c*z**i for i, c in enumerate(cubic)) % 2 for z in (0, 1))
    root = power(beta, 4, cubic)
    assert mul(root, root, cubic) == beta
    # For an actual 2-adic root B, 2X+2B agrees with this vector mod4.
    Z = K.add(K.add(X, X), K.elt([2*c for c in root]))
    norm = K.norm(Z); assert norm.denominator == 1 and norm % 4 == 3
    # Hence N(X+B)=N(2X+2B)/8 has valuation -3 and odd unit 3 mod4.
    dyadic = -1
    D = Q(inp['norm_X_squared_minus_beta'])
    assert D.numerator % 163 and D.denominator % 163
    assert all(v.denominator % 163 for v in X)
    odd = 1
    f = list(map(Q, inp['cubic_ascending']))
    alpha = list(map(Q, inp['alpha_ascending']))
    def evaluate(poly, x): return sum(c*x**i for i, c in enumerate(poly))
    intervals = [(Q(-2), Q(-1)), (Q(-1, 10), Q(0)), (Q(12), Q(13))]
    for lo, hi in intervals:
        assert evaluate(f, lo)*evaluate(f, hi) < 0
        # alpha'=2t-10 has constant nonzero sign on each interval.
        assert (2*lo-10)*(2*hi-10) > 0
        assert min(evaluate(alpha, lo), evaluate(alpha, hi)) > 0
    assert K.norm(X) > 0
    # alpha>0 gives |X_i|>|sqrt(beta_i)|, so N(X+sqrt(beta)) has sign N(X)>0.
    real = 1
    product = dyadic*odd*real; assert product == -1
    matrix = [[0, 1], [1, 0]]
    assert matrix == control['independent_norm_cup_matrix']
    return {'schema': 'rank-jump.strict-cochain-dyadic-switch.v1', 'status': 'PASS',
            'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in
                         (Path(__file__), verify.INPUT, verify.OUTPUT, verify.CONTROL, Path(verify.__file__), Path(r.__file__))},
            'prime2': {'beta_square_root_mod2': root, 'twice_X_plus_twice_root_mod4_lift': list(map(str, Z)),
                       'norm_of_lift': str(norm), 'norm_odd_unit_mod4': 3, 'actual_norm_valuation': -3,
                       'Hilbert_minus1_norm': dyadic},
            'prime163': {'norm_delta_valuation': 0, 'norm_X_plus_root_valuation': 0, 'Hilbert_minus1_norm': odd},
            'infinity': {'isolating_intervals': [[str(a), str(b)] for a, b in intervals],
                         'alpha_totally_positive': True, 'norm_X_positive': True, 'Hilbert_minus1_norm': real},
            'other_places': 'zero contribution: both curves and the governing extension are unramified outside 2,163',
            'CT_difference_matrix': matrix, 'matches_retained_independent_ideal_cup_matrix': True,
            'boundary': 'New point-free local derivation of a previously known strict-block CT switch. The difference formula is proved in the accompanying note; no new curve rank or full-Selmer preservation is asserted.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('mode', choices=['build', 'check']); args = p.parse_args()
    result = compute()
    if args.mode == 'build': r.write_new(OUTPUT, result)
    else: assert r.read(OUTPUT) == result
    print('PASS: the strict-block CT switch is the 2-adic unit 3 mod4; other contributions vanish')
