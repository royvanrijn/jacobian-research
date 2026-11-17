"""Exact actual-direction rank-one test; no imported producer code."""
from pathlib import Path
import hashlib
import json
from math import factorial
import sympy as S

ROOT = Path(__file__).resolve().parent


def rank_one_data(phi, x, t):
    poly = S.Poly(S.expand(phi), x)
    d = poly.degree()
    assert d >= 2
    top = poly.nth(d)
    second = poly.nth(d - 1)
    slope = S.cancel(S.diff(top, t) / (d * top))
    shift = S.cancel((S.diff(second, t) - (d - 1) * slope * second) / (d * top))
    residual = S.cancel(S.diff(phi, t) - (slope * x + shift) * S.diff(phi, x) - slope)
    is_one = S.cancel(S.diff(residual, x)) == 0
    if is_one:
        assert S.cancel(S.diff(phi, t) - residual - (slope + (slope * x + shift) * S.diff(phi, x))) == 0
    return is_one, slope, shift, residual


def run():
    x, t = S.symbols('x t', nonzero=True)
    alpha, beta, gamma = [S.Function(name)(t) for name in ['alpha', 'beta', 'gamma']]
    y = alpha * x + beta
    is_one, A, B, lam = rank_one_data(y ** 4 + y + gamma, x, t)
    assert is_one and S.cancel(lam - (S.diff(gamma, t) - S.diff(alpha, t) / alpha)) == 0
    assert S.cancel(A - S.diff(alpha, t) / alpha) == 0
    assert S.cancel(B - S.diff(beta, t) / alpha) == 0
    assert rank_one_data(-x ** 2 / 2 + t * x ** 4, x, t)[0] is False
    assert rank_one_data(t * (x ** 5 + x), x, t)[0] is False
    assert rank_one_data((x - t) ** 5 + (x - t), x, t) == (True, 0, -1, 0)
    assert rank_one_data(-x ** 2 / 2 + t * x - t ** 2 / 2, x, t) == (True, 0, -1, 0)
    c = S.symbols('c')
    is_one, _, _, lam = rank_one_data(t * (x ** 6 + c), x, t)
    assert is_one and S.cancel(lam - c + 1 / (6 * t)) == 0

    a, c, k, p, q = S.symbols('a c k p q')
    phase_repair = t * (p * x + q) + k * t ** 2 * (a + c * x)
    linear = S.expand(phase_repair).coeff(x, 1)
    defect = S.expand(phase_repair.subs(x, 0) + linear ** 2 / 2)
    assert defect.coeff(t, 4) == k ** 2 * c ** 2 / 2
    assert defect.subs(k, 0).coeff(t, 2) == p ** 2 / 2
    assert defect.subs({k: 0, p: 0}).coeff(t, 1) == q

    def wick(z_power, w_power, real_power):
        if z_power != w_power or real_power % 2:
            return 0
        half = real_power // 2
        return factorial(z_power) * factorial(real_power) // (2 ** half * factorial(half))

    assert wick(0, 2, 2) == 0  # Rejected Q=WT, P=WT, m=1.
    assert wick(2, 2, 2) == 2  # Correct Q=Z^2, P=WT, m=2.

    # Independent rank-two module check: quotient (1,0) is horizontal.
    connection = S.Matrix([[0, 0], [1, 1 / t]])
    quotient = S.Matrix([[1, 0]])
    assert quotient * connection == S.zeros(1, 2)
    assert S.Matrix.hstack(S.Matrix([1, 0]), connection[:, 0]).det() == 1
    return {'completed': True, 'generic_affine_shape_and_scalar': True,
            'translation_direction_rank': 1, 'scaling_x5_plus_x_rank_one': False,
            'Gaussian_quartic_rank_one': False, 'Gaussian_translation_rank_one': True,
            'constant_b_repair_constraints': True, 'rank_two_trivial_quotient': True,
            'rejected_sharpness_moment': 0, 'accepted_sharpness_moment': 2,
            'producer_package_executed': False,
            'universal_written_proofs_certified_by_checker': False,
            'independent_replay': False, 'external_review': False, 'formal_verification': False}


if __name__ == '__main__':
    result = run()
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (ROOT / 'actual_direction_rank_one_verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
