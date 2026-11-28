"""Exact rational Young and polynomial checks for the six-cell (2,2,1,1) gap."""

from functools import lru_cache
from itertools import combinations
from sympy import Matrix, Poly, Rational as Q, symbols

# The tableau order is ABCD, ACBD, ACDB, CABD, CADB, CDAB.
a, b, c, d, a_box, b_box, c_box, d_box = symbols("a b c d A B C D")
boxes = {a_box: (0, 1), b_box: (1, 1), c_box: (2, 0), d_box: (3, 0)}
order = (
    (a_box, b_box, c_box, d_box),
    (a_box, c_box, b_box, d_box),
    (a_box, c_box, d_box, b_box),
    (c_box, a_box, b_box, d_box),
    (c_box, a_box, d_box, b_box),
    (c_box, d_box, a_box, b_box),
)
tab = [((0, 0), (1, 0)) + tuple(boxes[x] for x in row) for row in order]
identity = Matrix.eye(6)

def tableau_inversions(tableau):
    reading = sorted(range(6), key=lambda i: tableau[i])
    return sum(i > j for i, j in combinations(reading, 2))

def adjacent(which):
    matrix = Matrix.zeros(6)
    for col, tableau in enumerate(tab):
        i, j = tableau[which]
        k, l = tableau[which + 1]
        axial = (l - k) - (j - i)
        matrix[col, col] = Q(1, axial)
        switched = list(tableau)
        switched[which], switched[which + 1] = switched[which + 1], switched[which]
        switched = tuple(switched)
        if switched in tab:
            row = tab.index(switched)
            matrix[row, col] = (
                Q(1) if tableau_inversions(switched) > tableau_inversions(tableau)
                else 1 - Q(1, axial**2)
            )
    return matrix

gen = [adjacent(i) for i in range(2, 5)]
for g in gen:
    assert g**2 == identity
assert gen[0] * gen[2] == gen[2] * gen[0]
for g, h in zip(gen, gen[1:]):
    assert g * h * g == h * g * h

@lru_cache(None)
def representation(permutation):
    if permutation == (0, 1, 2, 3):
        return identity
    for j in range(3):
        if permutation[j] > permutation[j + 1]:
            swapped = list(permutation)
            swapped[j], swapped[j + 1] = swapped[j + 1], swapped[j]
            return representation(tuple(swapped)) * gen[j]
    raise ValueError("Invalid permutation")

t_first = 2 * adjacent(1)
assert t_first == Matrix.diag(1, 1, 1, -2, -2, -2)
transports = []
for i in range(4):
    swap = list(range(4))
    swap[i], swap[0] = swap[0], swap[i]
    action = representation(tuple(swap))
    transports.append(action * t_first * action)

means = symbols("m1:5")
last = identity[:, 5]
field = Matrix.zeros(6, 1)
for i, j in combinations(range(4), 2):
    k, l = [index for index in range(4) if index not in (i, j)]
    field += (representation((i, j, k, l)) * last
              * (means[i] - means[j]) * (means[k] - means[l]))

rescale = Matrix.diag(Q(5, 9), Q(10, 27), Q(10, 9), Q(10, 9), Q(5, 12), Q(5, 3))
t = [b + c + d, c + d, d, 0]
t0 = a + b + c + d
gap_values = dict(zip(means, t))
small = 2*b*c + 3*b*d + 2*c**2 + 2*c*d
large = 4*b*c + 3*b*d + 4*c**2 + 4*c*d
displayed = Matrix([b*d, small, c*(b+c+d), c*(b+c+d), large, b*d])
field_on_gaps = rescale.inv() * field.subs(gap_values, simultaneous=True)
assert (field_on_gaps - displayed).expand() == Matrix.zeros(6, 1)

expected = [
    Matrix([0, 0, 0, 3*c*t0, 3*(4*a*c+3*a*d+large), 3*d*(a+b)]),
    Matrix([0, 8*(b+c)*t0,
            4*a*b+4*a*c+3*a*d+4*b**2+8*b*c+4*b*d+4*c**2+4*c*d,
            (b+4*c)*t0,
            4*a*b+16*a*c+12*a*d+4*b**2+20*b*c+13*b*d+16*c**2+16*c*d,
            3*b*d]),
    Matrix([4*b*t0, 4*(3*b+2*c)*t0, c*(a+4*b+4*c+4*d), 4*c*t0,
            3*a*b+4*a*c+3*b**2+19*b*c+12*b*d+16*c**2+16*c*d,
            b*(a+b+c+4*d)]),
    4*displayed,
]
prefix = Matrix.zeros(6, 1)
counts = []
for i, transport in enumerate(transports):
    replaced = t.copy()
    replaced[i] = t0
    prefix += ((identity - transport)
               * field.subs(dict(zip(means, replaced)), simultaneous=True))
    coefficients = rescale.inv() * prefix
    assert (coefficients - expected[i]).expand() == Matrix.zeros(6, 1)
    assert all(all(coefficient >= 0 for coefficient in Poly(term, a, b, c, d).coeffs())
               for term in coefficients)
    counts.append(sum(len(Poly(term, a, b, c, d).terms())
                      for term in coefficients if term.expand() != 0))

sample = dict(zip(means, [-1, -2, -3, -4]))
adjacent_sample = sample | {means[0]: 0}
witness = (Q(108, 5)
           * (field.subs(sample) - transports[0] * field.subs(adjacent_sample)))
assert list(witness) == [-12, -40, -24, 264, 531, 180]
print("Positive field polynomials and all four exact prefix vectors verified.")
print("Nonzero monomials in the four prefix vectors:", counts)
print("Direct cutoff obstruction:", list(witness))

