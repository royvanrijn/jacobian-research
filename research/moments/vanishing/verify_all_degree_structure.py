"""Small exact transcription checks; universal arguments are in the proof note."""
from pathlib import Path
from math import comb, factorial
import hashlib
import json
import sympy as S

ROOT = Path(__file__).resolve().parent


def real_moment(n):
    return 0 if n % 2 else factorial(n) // (2 ** (n // 2) * factorial(n // 2))


def wick(z, w, t=0):
    return factorial(z) * real_moment(t) if z == w else 0


def run():
    checks = 0
    first = []
    for d in range(2, 8):
        for m in range(2 * (d + 1) + 1):
            direct = sum(comb(m, b) * wick(m - b, d * b) for b in range(m + 1))
            closed = factorial(m) // factorial(m // (d + 1)) if m % (d + 1) == 0 else 0
            assert direct == closed
            checks += 1
        nu = 2 * d - 1 if d % 2 == 0 else 4 * d - 2
        nonzero = []
        for m in range(2 * nu + 1):
            direct = sum(comb(m, b) * wick(d * (m - b), (d - 1) * b, b) for b in range(m + 1))
            k = m // (2 * d - 1)
            closed = comb(m, d * k) * factorial(d * (d - 1) * k) * real_moment(d * k) if m % (2 * d - 1) == 0 and d * k % 2 == 0 else 0
            assert direct == closed
            if m and direct:
                nonzero.append(m)
            checks += 1
        assert nonzero[0] == nu
        first.append({'degree': d, 'binary_first': d + 1, 'homogeneous_first': nu})

    s, alpha, beta, gamma = S.symbols('s alpha beta gamma')
    for k in (1, 2, 3):
        exponent = alpha * (s + beta) + gamma ** 2 / 2 + gamma * s * (s + beta) ** k - beta * s * (s + beta) ** (2 * k) / 2
        zero = {alpha: 0, beta: 0, gamma: 0}
        assert exponent.subs(zero) == 0
        assert S.diff(exponent, gamma).subs(zero) == s ** (k + 1)
        assert S.diff(exponent, beta).subs(zero) == -s ** (2 * k + 1) / 2

    # Coefficient polarization makes distinct occupations independent.
    # Balanced safe shear: pure moments vanish despite a nonzero universal M_4.
    c0, c1, c2 = S.symbols('c0 c1 c2')
    # Support W, Z*T, Z^3 for k=1: balanced occupations first occur at m=4.
    support = [(0, 1, 0), (1, 0, 1), (3, 0, 0)]
    universal = []
    for m in range(1, 5):
        value = 0
        for b in range(m + 1):
            for c in range(m - b + 1):
                a = m - b - c
                z, w, t = [sum(n * v[i] for n, v in zip((a, b, c), support)) for i in range(3)]
                value += factorial(m) // (factorial(a) * factorial(b) * factorial(c)) * wick(z, w, t) * c0 ** a * c1 ** b * c2 ** c
        universal.append(S.expand(value))
    assert universal[:3] == [0, 0, 0]
    assert universal[3] == 12 * c0 ** 2 * c1 ** 2 + 24 * c0 ** 3 * c2
    assert universal[3].subs({c0: 1, c1: 1, c2: S.Rational(-1, 2)}) == 0
    return {'completed': True, 'exact_moment_comparisons': checks,
            'first_nonzero_examples': first, 'shear_source_checks': 9,
            'polarized_support_first': 4, 'polarized_M4': str(universal[3]),
            'universal_proofs_certified_by_checker': False,
            'independent_replay': False, 'formal_verification': False, 'external_review': False}


if __name__ == '__main__':
    result = run()
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    target = ROOT / 'all_degree_structure_verification.json'
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
