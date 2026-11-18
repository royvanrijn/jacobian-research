#!/usr/bin/env python3
"""Check the two fixed Laurent minors and their sole exceptional specialization."""
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json

from verify_uniform_fox_ideal_models import free_action, reconstruct, red

HERE = Path(__file__).resolve().parent


def laurent_row(word):
    run = 0
    row = [Counter() for _ in range(4)]
    for letter in word:
        slot = abs(letter) - 1
        if letter < 0:
            run -= 1
            row[slot][run] -= 1
        else:
            row[slot][run] += 1
            run += 1
    assert run == 0
    return row


def product(a, b):
    result = Counter()
    for i, u in a.items():
        for j, v in b.items():
            result[i + j] += u * v
    return result


def normalized_minor(a, b, i, j):
    c = product(a[i], b[j])
    c.subtract(product(a[j], b[i]))
    support = [i for i, x in c.items() if x]
    lo, hi = min(support), max(support)
    return lo, [c[i] for i in range(lo, hi + 1)]


def determinant(matrix):
    """Fraction-free elimination with row pivots."""
    a = [r[:] for r in matrix]
    n = len(a)
    previous = 1
    sign = 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[k][k] * a[i][j] - a[i][k] * a[k][j]
                a[i][j], remainder = divmod(numerator, previous)
                assert remainder == 0
            a[i][k] = 0
        previous = a[k][k]
    return sign * a[-1][-1]


def resultant(f, g):
    m, n = len(f) - 1, len(g) - 1
    a = [[0] * i + f[::-1] + [0] * (n - i - 1) for i in range(n)]
    b = [[0] * i + g[::-1] + [0] * (m - i - 1) for i in range(m)]
    return determinant(a + b)


def trim(f):
    while f and not f[-1]:
        f.pop()
    return f


def gcd_mod(a, b, p):
    a, b = trim([x % p for x in a]), trim([x % p for x in b])
    while b:
        while len(a) >= len(b):
            shift = len(a) - len(b)
            factor = a[-1] * pow(b[-1], -1, p)
            for i, x in enumerate(b):
                a[i + shift] = (a[i + shift] - factor * x) % p
            trim(a)
        a, b = b, a
    u = pow(a[-1], -1, p)
    return [u * c % p for c in a]


def specialize(row, t, p):
    return [sum(v * pow(t, e, p) for e, v in entry.items()) % p for entry in row]


def positive_crossing(a, b, t, p):
    u = pow(t, -1, p)
    return b, [(u * x + (1 - u) * y) % p for x, y in zip(a, b)]


def verify():
    literal, _, _, _ = reconstruct()
    relators = literal["all_four_closure_relators"]
    h0, h3 = map(laurent_row, (relators[0], relators[3]))
    s, f = normalized_minor(h0, h3, 0, 3)
    u, g = normalized_minor(h0, h3, 1, 3)
    assert (s, u) == (-5, -5)
    assert f == [1, -6, 15, -21, 18, -8, 0, 4, -3, 1]
    assert g == [-2, 10, -25, 40, -45, 36, -22, 9, -2]
    r = resultant(f, g)
    assert r == -83 and gcd_mod(f, g, 83) == [24, 1]
    # Corrupt one coefficient to test the actual resultant computation.
    damaged = f[:]
    damaged[0] += 1
    assert resultant(damaged, g) != -83
    p, t = 83, 59
    w = specialize(h0, t, p)
    y0 = specialize(laurent_row(relators[1]), t, p)
    bottom = free_action(literal["literal_braid_word"] + [2] * 6)
    y1 = specialize(laurent_row(red(bottom[1] + [-2])), t, p)
    assert w == [10, 12, 60, 1] and y0 == [51, 66, 7, 42] and y1 == [1, 26, 55, 1]
    base_bottom_rows = [specialize(laurent_row(relators[i]), t, p) for i in (1, 2)]
    for row, i in zip(base_bottom_rows, (1, 2)):
        row[i] = (row[i] + 1) % p
    a, b = base_bottom_rows
    for _ in range(6):
        a, b = positive_crossing(a, b, t, p)
    a[1] = (a[1] - 1) % p
    assert a == y1
    def minors(v):
        return [(w[0] * v[j] - w[j] * v[0]) % p for j in (1, 2)]
    assert minors(y0) == [48, 81] and minors(y1) == [82, 75]
    combo = lambda v: (77 * v[0] - 34 * v[1]) % p
    assert combo(minors(y0)) == combo(minors(y1)) == 29
    assert pow(t, -6, p) == 81 != 1
    return {
        "status": "PASS exact scalar determinantal unit over all parameters",
        "source_sha256": {n: sha256((HERE / n).read_bytes()).hexdigest() for n in
                           ("uniform_fox_ideal_inputs.json", "verify_uniform_fox_ideal_models.py")},
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "fixed_rows_zero_based": [0, 3],
        "normalized_pairs": [[0, 3], [1, 3]], "polynomials_ascending": [f, g],
        "resultant": r, "exceptional_characteristic": p, "exceptional_t": t,
        "mod_83_gcd_ascending": [24, 1], "R0": w, "R1_base": y0, "R1_next": y1,
        "mod_83_minors_base_and_next": [minors(y0), minors(y1)],
        "minor_combination": {"weights": [77, -34], "constant": 29, "tail_nontrivial_eigenvalue": 81},
        "all_n_reason": "The normalized two-strand derivative has distinct eigenvalues 1 and -t^{-1} at t=59. Six-step powers are affine linear in 81^n, so the verified two values force the minor combination to be the constant 29.",
        "rejected_control": "Change constant term of the first normalized minor.",
        "assurance": {"whole_proof_independent_replay": False, "formal_verification": False,
                      "external_specialist_review": False}}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
