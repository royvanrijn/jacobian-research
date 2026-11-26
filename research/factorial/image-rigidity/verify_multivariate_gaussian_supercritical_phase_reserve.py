"""Certify (7) in the supercritical Gaussian phase-reserve proof."""

from collections import Counter

S = 1 << 45

class I:
    """Closed interval with endpoints in units of 1/S."""

    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi

    def __add__(self, other):
        other = box(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -box(other)

    def __rsub__(self, other):
        return box(other) + -self

    def __mul__(self, other):
        other = box(other)
        if self.lo >= 0 and other.lo >= 0:
            return I(self.lo * other.lo // S, ceil(self.hi * other.hi, S))
        endpoints = [x * y for x in (self.lo, self.hi)
                     for y in (other.lo, other.hi)]
        return I(min(endpoints) // S, ceil(max(endpoints), S))

    __rmul__ = __mul__

def box(value):
    return value if isinstance(value, I) else I(value * S, value * S)

def ceil(value, divisor):
    return -((-value) // divisor)

def classify(d, r, v):
    e = 1 - r
    b = e * e * v
    a = r - b
    if a.hi < 0 or 2 * b.lo > S:
        return "domain"
    # Tighten the enclosure using constraints that define the domain.
    a = I(max(a.lo, 0), a.hi)
    b = I(b.lo, min(b.hi, S // 2))
    da = d * a
    if (4 * v * (1 + da) * r * r - 1 - 2 * r).lo > 0:
        return "domain"

    db = d * b
    f = 1 + da - db
    c = 1 + 2 * e * v
    dc = d * c
    one_d = 1 + d
    ub = 1 - db
    j = f * c - 2 * d * one_d * a * e * v
    a_term = a * (1 + da)
    y_test = dc * (j + a_term * dc) - v * ub * one_d * one_d
    if y_test.lo > 0:
        return "height"

    h = (19 - 10 * a) * (3 - 2 * db)
    g = 24 * f * (3 - 2 * db + b) - h * one_d * a
    hc = h * c
    t_test = j * g * hc + one_d * (a_term * hc * hc * ub - v * g * g)
    if t_test.lo > 0:
        return "phase"

initial = (I(0, S), I(0, S), I(0, 2 * S))
stack = [(initial, 0)]
counts = Counter()
while stack:
    coordinates, depth = stack.pop()
    verdict = classify(*coordinates)
    if verdict:
        counts[verdict] += 1
        counts["depth"] = max(counts["depth"], depth)
        continue
    if depth >= 65:
        raise ArithmeticError(tuple((x.lo, x.hi) for x in coordinates))
    axis = max(range(3), key=lambda i: coordinates[i].hi - coordinates[i].lo)
    x = coordinates[axis]
    mid = (x.lo + x.hi) // 2
    for part in I(x.lo, mid), I(mid, x.hi):
        stack.append((coordinates[:axis] + (part,) + coordinates[axis + 1:], depth + 1))
expected = {"domain": 1691, "height": 6084, "phase": 5162, "depth": 24}
if counts != expected:
    raise ArithmeticError((dict(counts), expected))
print(dict(counts))
