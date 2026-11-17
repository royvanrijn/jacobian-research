"""Exact algebraic boundary checks; not a contour/asymptotic certificate."""
from pathlib import Path
import hashlib
import json
import sympy as S

ROOT = Path(__file__).resolve().parent


def run():
    x, y, t = S.symbols('x y t', nonzero=True)
    a, b, c, d, p = S.symbols('a b c d p', nonzero=True)
    mobius = (a * y + b) / (c * y + d)
    preimage = (d * p - b) / (a - c * p)
    assert S.cancel(S.diff(mobius, y) / (mobius - p) - (1 / (y - preimage) - 1 / (y + d / c))) == 0
    # Pullback through a ramification point preserves a simple log pole.
    for e in (2, 3, 5):
        assert S.cancel(S.diff(y ** e, y) / y ** e - e / y) == 0

    # Moving-pole derivative must be included in the actual connection.
    omega = 1 / (x - t) + 2 / (x + t ** 2)
    phi = t / (x - t) - (x + t ** 2) ** 2 / t
    rho = 1 / t
    actual = S.diff(omega, t) + S.diff(phi, t) * omega + S.diff(rho, t) / rho * omega
    full = S.diff(omega * S.exp(phi) * rho, t) / (S.exp(phi) * rho)
    assert S.cancel(S.expand_power_exp(full - actual)) == 0
    assert S.cancel(S.diff(omega, t) - (1 / (x - t) ** 2 - 4 * t / (x + t ** 2) ** 2)) == 0
    assert S.diff(omega, t) != 0

    # A differential residue checks the rational-extension obstruction.
    assert S.residue(1 / t, t, 0) == 1
    assert S.diff(S.log(t), t) == 1 / t
    # Boundary failures: growing residues and absent normalization.
    assert S.residue(t * S.exp(x) / x, x, 0) == t
    assert S.residue(S.exp(x) / x, x, 0) == 1
    # The higher-pole Gaussian counterexample after x=t*y, t>0.
    gaussian_phase = -x ** 2 / (2 * t ** 2)
    assert gaussian_phase.subs(x, t * y) == -y ** 2 / 2
    assert S.diff(t * y, y) / t == 1
    # dx has a double pole in y=1/x, whereas dx/x has a simple one.
    assert S.diff(1 / y, y) == -1 / y ** 2
    assert S.cancel(S.diff(1 / y, y) / (1 / y)) == -1 / y
    return {'completed': True, 'mobius_log_pullback': True,
            'ramification_indices_checked': [2, 3, 5],
            'moving_pole_full_connection': True,
            'rational_extension_obstruction_residue': 1,
            'boundary_examples': ['higher-pole Gaussian scaling', 'growing residue', 'absent decay normalization'],
            'universal_contour_proof_certified_by_checker': False,
            'period_comparison_certified_by_checker': False,
            'independent_replay': False, 'external_review': False, 'formal_verification': False}


if __name__ == '__main__':
    result = run()
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (ROOT / 'logarithmic_curve_criterion_verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
