"""Exact all-index algebraic identities; not an analytic proof verifier."""
from pathlib import Path
import hashlib
import json
import sympy as S

ROOT = Path(__file__).resolve().parent


def run():
    m, j = S.symbols('m j', integer=True, nonnegative=True)
    assert S.combsimp(S.binomial(m, j) * S.factorial(j) / S.factorial(m) - 1 / S.factorial(m - j)) == 0
    k = S.symbols('k', integer=True, positive=True)
    c, u = S.symbols('c u', nonzero=True)
    # General multiplicity-k leading root model, t=1/(c*u^k).
    parameter = 1 / (c * u ** k)
    residue = -1 / (parameter * k * c * u ** (k - 1))
    assert S.simplify(residue + u / k) == 0

    x, t = S.symbols('x t', nonzero=True)
    d = S.symbols('d', integer=True, positive=True)
    H, r = S.Function('H')(x), S.Function('r')(x)
    rho = t ** (-1 / d)
    v_original = 1 / (1 - t * H)
    v = v_original / rho
    assert S.simplify(rho * v - v_original) == 0
    assert S.simplify(S.diff(rho, t) / rho + 1 / (d * t)) == 0
    # The full normalized tangent includes the derivative of the moving form.
    tangent = S.diff(v, t) + r * v + S.diff(rho, t) / rho * v
    assert S.simplify(rho * tangent - (S.diff(v_original, t) + r * v_original)) == 0
    assert S.diff(v_original, t) == H / (1 - t * H) ** 2

    z, w, T = S.symbols('z w T')
    A = S.Function('A')(z, T)
    h = z * A
    assert S.diff(h, z).subs(z, 0) == A.subs(z, 0)
    # Reflection preserves the neutral radial monomial.
    assert (z * w).xreplace({z: w, w: z}) == z * w
    return {'completed': True, 'all_index_resolvent_coefficient_identity': True,
            'arbitrary_multiplicity_root_residue': '-u/k',
            'full_moving_form_gauge_identity': True,
            'trapped_neutral_coefficient': True,
            'positive_negative_reflection': True,
            'universal_analytic_input_certified_by_checker': False,
            'independent_replay': False, 'external_review': False, 'formal_verification': False}


if __name__ == '__main__':
    result = run()
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (ROOT / 'neutral_face_verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
