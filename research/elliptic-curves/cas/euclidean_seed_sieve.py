"""Exact, Sage-free Euclidean bisections. Incidence is NOT rank admission.

Coefficient lists are low-to-high over Q. No floating point, factorization,
record parameters, or exceptional point coordinates enter this module.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt


def poly(values):
    p = tuple(map(F, values)) or (F(0),)
    while len(p) > 1 and not p[-1]:
        p = p[:-1]
    return p


def add(a, b):
    return poly((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                for i in range(max(len(a), len(b))))


def scale(a, c):
    return poly(x * F(c) for x in a)


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return poly(out)


def divrem(a, b):
    a, b = poly(a), poly(b)
    if b == (0,):
        raise ZeroDivisionError('zero polynomial divisor')
    q = [F(0)] * max(1, len(a)-len(b)+1)
    while a != (0,) and len(a) >= len(b):
        k, c = len(a)-len(b), a[-1]/b[-1]
        q[k] += c
        a = sub(a, (F(0),)*k + scale(b, c))
    return poly(q), a


def exact_div(a, b):
    q, r = divrem(a, b)
    if r != (0,):
        raise ArithmeticError('non-exact polynomial division')
    return q


def inverse_mod(a, modulus):
    r0, r1, s0, s1 = poly(modulus), poly(a), (F(0),), (F(1),)
    while r1 != (0,):
        q, r = divrem(r0, r1)
        r0, r1, s0, s1 = r1, r, s1, sub(s0, mul(q, s1))
    if len(r0) != 1 or not r0[0]:
        raise ArithmeticError('trace numerator is not a unit modulo h^2')
    return divrem(scale(s0, 1/r0[0]), modulus)[1]


def evaluate(a, t):
    t, out = F(t), F(0)
    for c in reversed(a):
        out = out*t + F(c)
    return out


def square_certificate(value):
    """An exact square/nonsquare receipt, including negative and zero values."""
    value = F(value)
    result = {'value': str(value), 'numerator': str(value.numerator),
              'denominator': str(value.denominator)}
    if value < 0:
        return {**result, 'status': 'NONSQUARE', 'negative': True}
    a, b = isqrt(value.numerator), isqrt(value.denominator)
    square = a*a == value.numerator and b*b == value.denominator
    return {**result, 'status': 'SQUARE' if square else 'NONSQUARE',
            'numerator_floor_sqrt': str(a), 'denominator_floor_sqrt': str(b),
            'root': str(F(a, b)) if square else None}


def conic_from_trace(A, B, h, nx, ny):
    """Checked application of the norm-ten formula, not a new family theorem.

    Recheck the trace equation, every division, degree bound, smoothness, and
    the lift identically in Q[t,W]/(W^2-q). Exceptional charts are not guessed.
    """
    A, B, h, nx, ny = map(poly, (A, B, h, nx, ny))
    if not (len(A) <= 9 and len(B) <= 13 and len(h) == 4 and h[-1] == 1
            and len(nx) <= 11 and len(ny) <= 16):
        raise ArithmeticError('outside certified degree-three trace chart')
    h2 = mul(h, h)
    if mul(ny, ny) != add(add(mul(mul(nx, nx), nx), mul(mul(A, nx), mul(h2, h2))),
                         mul(B, mul(mul(h2, h2), h2))):
        raise ArithmeticError('trace is not on the parent')
    m = divrem(scale(mul(ny, inverse_mod(nx, h2)), -1), h2)[1]
    g = exact_div(add(mul(m, nx), ny), h2)
    b = exact_div(sub(mul(m, m), nx), h2)
    k = exact_div(sub(mul(m, b), scale(g, 2)), h2)
    q = exact_div(sub(sub(scale(mul(m, k), 4), scale(mul(b, b), 3)), scale(A, 4)), h2)
    if not all(len(p) <= bound+1 for p, bound in ((m, 5), (g, 9), (b, 4), (k, 3))):
        raise ArithmeticError('Euclidean coefficient degree exceeded')
    if len(q) not in (2, 3) or (len(q) == 3 and q[1]**2 == 4*q[0]*q[2]):
        raise ArithmeticError('not a smooth nonsplit degree-one/two cover')
    x0, x1, y0, y1 = scale(b, F(1, 2)), scale(h, F(1, 2)), scale(mul(h, k), F(-1, 2)), scale(m, F(-1, 2))
    const = add(mul(mul(x0, x0), x0), add(scale(mul(mul(x0, mul(x1, x1)), q), 3), add(mul(A, x0), B)))
    linear = add(scale(mul(mul(x0, x0), x1), 3), add(mul(mul(mul(x1, x1), x1), q), mul(A, x1)))
    if add(mul(y0, y0), mul(mul(y1, y1), q)) != const or scale(mul(y0, y1), 2) != linear:
        raise ArithmeticError('polynomial map replay failed')
    values = {'h': h, 'nx': nx, 'ny': ny, 'm': m, 'g': g, 'b': b, 'k': k, 'q': q}
    return {**{name: list(map(str, p)) for name, p in values.items()},
            'maps': [list(map(str, p)) for p in (x0, x1, y0, y1)]}


def incidence(conic, parameter):
    """One branch per split: its companion differs by an inherited trace.

    Ramification gives 2Q=P_w and cannot add a rational MW direction. A
    nonzero split still needs exact independence admission on its fibre.
    """
    t = F(parameter)
    cert = square_certificate(evaluate(conic['q'], t))
    base = {'parameter': str(t), 'square_certificate': cert}
    if cert['status'] != 'SQUARE':
        return {**base, 'status': 'EXACT_NONSPLIT'}
    W = F(cert['root'])
    if not W:
        return {**base, 'status': 'RAMIFIED_NO_NEW_DIRECTION'}
    x0, x1, y0, y1 = (evaluate(p, t) for p in conic['maps'])
    return {**base, 'status': 'SPLIT_REQUIRES_ADMISSION',
            'point': [str(x0+x1*W), str(y0+y1*W)]}


def address(index):
    """Signed Calkin-Wilf addresses; unique, O(log index), no target data."""
    if isinstance(index, bool) or not isinstance(index, int) or index < 0:
        raise ValueError('nonnegative integer address required')
    a = b = 1
    for bit in bin(index//2+1)[3:]:
        if bit == '0':
            b += a
        else:
            a += b
    return str(F(a if index % 2 == 0 else -a, b))


def ordered_orbits(rows, offset, count):
    if offset < 0 or count < 1:
        raise ValueError('invalid orbit window')
    selected = [r for r in rows if r['category'] == 'rational' and F(r['minimum_norm']) == 10]
    selected.sort(key=lambda r: (sha256(('euclidean-seed-v1/'+r['orbit_mask']).encode()).hexdigest(), int(r['orbit_mask'])))
    if len({r['orbit_mask'] for r in selected}) != len(selected):
        raise ValueError('duplicate orbit IDs')
    if offset+count > len(selected):
        raise ValueError('orbit window exceeds the frozen rational norm-ten atlas')
    return [{'mask': int(r['orbit_mask']), 'word': list(map(int, r['parent_MW17_w'].split()))}
            for r in selected[offset:offset+count]]
