"""Generic symbolic identities and small exact normal-form examples."""
from pathlib import Path
import hashlib
import json
import sympy as S

ROOT = Path(__file__).resolve().parent


def bezout_vector(polys, z):
    g = S.Integer(0)
    vector = [S.Integer(0)] * len(polys)
    for i, b in enumerate(polys):
        if b == 0:
            continue
        if g == 0:
            leading = S.Poly(b, z).LC()
            g = b / leading
            vector[i] = 1 / leading
        else:
            u, v, new_g = S.gcdex(g, b, z)
            vector = [S.expand(u * a) for a in vector]
            vector[i] += v
            g = S.expand(new_g)
    assert S.expand(sum(a * b for a, b in zip(vector, polys)) - g) == 0
    return g, vector


def run():
    z, w, T1, T2 = S.symbols('z w T1 T2')
    Ts = [T1, T2]
    us = [S.Function('u1')(z), S.Function('u2')(z)]
    vs = [S.Function('v1')(z), S.Function('v2')(z)]

    def shear(vec):
        return {z: z, w: w - sum(u * T for u, T in zip(vec, Ts)) - z * sum(u ** 2 for u in vec) / 2,
                **{T: T + z * u for T, u in zip(Ts, vec)}}

    phi = shear(us)
    quad = 2 * z * w + sum(T ** 2 for T in Ts)
    assert S.expand(quad.subs(phi, simultaneous=True) - quad) == 0
    assert S.det(S.Matrix([phi[v] for v in [z, T1, T2, w]]).jacobian([z, T1, T2, w])) == 1
    for var in [z, w, T1, T2]:
        composed = phi[var].subs(shear(vs), simultaneous=True)
        assert S.expand(composed - shear([u + v for u, v in zip(us, vs)])[var]) == 0
        assert S.expand(phi[var].subs(shear([-u for u in us]), simultaneous=True) - var) == 0

    f = S.Function('f')(z, w, T1, T2)
    deriv = lambda f: sum(u * (z * S.diff(f, T) - T * S.diff(f, w)) for u, T in zip(us, Ts))
    lap = lambda f: 2 * S.diff(f, z, w) + sum(S.diff(f, T, 2) for T in Ts)
    expected = 2 * sum(S.diff(u, z) * (z * S.diff(f, T, w) - T * S.diff(f, w, 2)) for u, T in zip(us, Ts))
    assert S.expand(lap(deriv(f)) - deriv(lap(f)) - expected) == 0
    assert S.expand(lap(phi[w] ** 2) + 4 * sum(S.diff(u, z) * (T + z * u) for u, T in zip(us, Ts))) == 0
    for u, T in zip(us, Ts):
        assert S.diff(u * S.diff(f, T), w) - S.diff(u * S.diff(f, w), T) == 0

    b, u, rem = [S.Function(name)(z) for name in ['b', 'u', 'R']]
    G = S.Function('G')
    a = z * b * u + rem
    direct = (w + u * T1 - z * u ** 2 / 2) * (a + b * (T1 - z * u)) + G(z, T1 - z * u)
    normal = w * (rem + b * T1) + (u * T1 - z * u ** 2 / 2) * (rem + b * T1) + G(z, T1 - z * u)
    assert S.expand(direct - normal) == 0

    examples = []
    for bs, a in [([(z - 1) ** 2, z ** 2], z ** 8 + 3 * z ** 2 + 2 * z + 5),
                  ([z ** 2 * (z - 1), z * (z - 1) ** 2], z ** 7 + 2 * z ** 3 + 3 * z + 4)]:
        g, vec = bezout_vector(bs, z)
        q, rem = S.div(a, z * g, z)
        chosen = [S.expand(q * v) for v in vec]
        assert S.expand(a - z * sum(b * u for b, u in zip(bs, chosen)) - rem) == 0
        assert S.degree(rem, z) < S.degree(g, z) + 1
        examples.append({'gcd': str(g), 'remainder': str(rem)})
    image = (w - z * T1 - z ** 3 / 2) ** 2
    assert S.expand(lap(image) + 4 * (T1 + z ** 2)) == 0
    return {'completed': True, 'generic_spectators_checked': 2,
            'quadratic_and_jacobian': True, 'group_and_inverse': True,
            'gaussian_ibp_identity': True, 'full_laplacian_commutator': True,
            'finite_laplacian_obstruction': True,
            'generic_normal_form_identity': True, 'bezout_examples': examples,
            'laplacian_counterexample': '-4*(T1+z**2)',
            'all_dimension_written_proof_certified_by_checker': False,
            'independent_replay': False, 'external_review': False, 'formal_verification': False}


if __name__ == '__main__':
    result = run()
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (ROOT / 'gaussian_shear_normal_form_verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
