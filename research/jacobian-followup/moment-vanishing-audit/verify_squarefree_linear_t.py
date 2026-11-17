"""Universal formal-coordinate identities, not an analytic certificate."""
from pathlib import Path
import hashlib
import json
import sympy as S

ROOT = Path(__file__).resolve().parent


def run():
    z, t, T = S.symbols('z t T', nonzero=True)
    a, b = S.Function('a')(z), S.Function('b')(z)
    F = (z / t - a) / b
    D = 1 - t * (S.diff(a, z) + S.diff(b, z) * F)
    assert S.cancel(S.diff(F, z) - D / (t * b)) == 0
    g = S.Function('g')(t, T)
    bg = S.Function('b')(g)
    w = S.diff(g, T) / (t * bg)
    velocity = S.diff(g, t) / S.diff(g, T)
    assert S.simplify(S.diff(w, t) + w / t - S.diff(w * velocity, T)) == 0
    Q, B = S.Function('Q')(t, T), S.Function('B')(t, T)
    actual = S.diff(w, t) * Q + w * S.diff(Q, t) + w * Q * S.diff(B, t)
    tangent = w * (S.diff(Q, t) - S.diff(Q, T) * velocity +
                   (S.diff(B, t) - (S.diff(B, T) - T) * velocity - 1 / t) * Q)
    ibp = S.diff(w * velocity, T) * Q + w * velocity * (S.diff(Q, T) + (S.diff(B, T) - T) * Q)
    assert S.simplify(actual - tangent - ibp) == 0
    for roots in [[-1, 1], [0, 1, 2], [-2, -1, 1, 3]]:
        poly = S.prod(z - r for r in roots)
        residues = sum(1 / (S.diff(poly, z).subs(z, r) * (z - r)) for r in roots)
        assert S.cancel(1 / poly - residues) == 0
    lam, mu = S.symbols('lambda mu', nonzero=True)
    X = T + lam / mu - 1 / (mu * t)
    trapped_D = 1 - t * (lam + mu * T)
    assert S.cancel(trapped_D + t * mu * X) == 0
    assert S.cancel(t * (-mu) * X / trapped_D) == 1
    assert S.limit(t * X, t, 0) == -1 / mu
    return {'completed': True, 'arbitrary_a_b_jacobian': True,
            'arbitrary_formal_g_conservation': True, 'full_tangent_identity': True,
            'partial_fraction_root_counts': [2, 3, 4], 'trapped_component_retained': True,
            'universal_analytic_input_certified_by_checker': False,
            'independent_replay': False, 'external_review': False, 'formal_verification': False}


if __name__ == '__main__':
    result = run()
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (ROOT / 'squarefree_linear_t_verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
