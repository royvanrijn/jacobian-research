#!/usr/bin/env sage
"""Smith-saturate the rejected intermediate CM-43 q=60 RR module over GF(p).

For ``D=Q79+4O-43F`` use the smooth-fiber basis
``{1,x,y,x^2,h*z_Q}`` and clear ``D0=h^2*x-Nx``.  In numerator rows
``[1,x,x^2,x^3,y,xy]`` the five columns form the polynomial matrix

    [-Nx,   0,   0,   0, Ny]
    [ h^2,-Nx,   0,   0,  0]
    [   0, h^2,  0,-Nx,  0]
    [   0,   0,  0, h^2, 0]
    [   0,   0,-Nx,  0, h^3]
    [   0,   0, h^2, 0,  0].

The Smith left transformation saturates the rank-five submodule inside the
six-row ambient polynomial module.  A shifted weak-Popov reduction with row
weights ``[0,4,8,12,6,10]`` then gives the finite global degree bounds.
Imposing coefficients T^0,...,T^42=0 is the exact -43F gate.

The resulting nullity six is a certified boundary computation, not the
generic q=60 pencil.  The full CM-43 chamber audit shows that this q=60 class
has fixed CM-only sections ``4*(P1-P2)+(P3-P2)`` and collapses to the old
fiber.  Thus the four excess dimensions record boundary degeneration; they
must not be removed by invented II* gates.  Recovering the generic pencil
requires deformation away from CM-43 or a non-collapsed boundary anchor.
"""

import argparse
from pathlib import Path

from sage.all import *


parser = argparse.ArgumentParser()
parser.add_argument("--p", type=int, default=101)
args = parser.parse_args()

load(str(Path(__file__).resolve().parent / "verify_cm43_humbert8_anchor.sage"))

p = ZZ(args.p)
assert p.is_prime()
field = GF(p)
R = PolynomialRing(field, "T")
Tmod = R.gen()

Q79 = 4*point1-5*point2+point3
Nxq = RT(Q79[0].numerator())
Nyq = RT(Q79[1].numerator())
hq = RT(Q79[0].denominator().sqrt())
assert Q79[0].denominator() == hq**2
assert Q79[1].denominator() == hq**3


def reduce_polynomial(poly):
    assert all(QQ(value).denominator() % p for value in poly)
    return R([
        field(QQ(value).numerator())/field(QQ(value).denominator())
        for value in poly
    ])


Nx, Ny, h = map(reduce_polynomial, (Nxq, Nyq, hq))
h2 = h**2
h3 = h**3
module_matrix = matrix(R, (
    (-Nx, 0, 0, 0, Ny),
    (h2, -Nx, 0, 0, 0),
    (0, h2, 0, -Nx, 0),
    (0, 0, 0, h2, 0),
    (0, 0, -Nx, 0, h3),
    (0, 0, h2, 0, 0),
))

smith, left, right = module_matrix.smith_form()
assert smith == left*module_matrix*right
assert left.is_invertible() and right.is_invertible()
invariant_degrees = tuple(
    ZZ(smith[index, index].degree()) for index in range(5)
)
assert invariant_degrees == (0, 0, 0, 0, 290)

# If S=L*M*R, then the first five columns of L^{-1} generate the saturation
# of the image of M in R^6.
saturated = left.inverse()[:, :5]
assert saturated.nrows() == 6 and saturated.ncols() == 5

row_shifts = (0, 4, 8, 12, 6, 10)


def column_pivot(column):
    candidates = []
    for row, entry in enumerate(column):
        if entry:
            candidates.append((ZZ(entry.degree())+row_shifts[row], row))
    assert candidates
    # Highest row breaks a shifted-degree tie, the standard weak-Popov rule.
    shifted_degree, row = max(candidates)
    return shifted_degree, row


def shifted_weak_popov(basis):
    basis = matrix(R, basis)
    while True:
        pivots = [column_pivot(basis.column(index)) for index in range(5)]
        collision = None
        for left_index in range(5):
            for right_index in range(left_index+1, 5):
                if pivots[left_index][1] == pivots[right_index][1]:
                    collision = (left_index, right_index)
                    break
            if collision is not None:
                break
        if collision is None:
            return basis
        first, second = collision
        if pivots[first][0] > pivots[second][0]:
            first, second = second, first
        pivot_row = pivots[first][1]
        reducer = basis[pivot_row, first]
        target = basis[pivot_row, second]
        degree_delta = ZZ(target.degree()-reducer.degree())
        assert degree_delta >= 0
        factor = target.leading_coefficient()/reducer.leading_coefficient()
        basis.set_column(
            second,
            basis.column(second)-factor*Tmod**degree_delta*basis.column(first),
        )


popov = shifted_weak_popov(saturated)
pivots = tuple(column_pivot(popov.column(index)) for index in range(5))
assert len(set(row for _, row in pivots)) == 5
generator_weights = tuple(degree-120 for degree, _ in pivots)
assert sorted(generator_weights) == [-45, -44, -43, -41, -39]

# All polynomial multiples whose infinity weight is nonpositive.
global_columns = []
column_labels = []
for generator, weight in enumerate(generator_weights):
    assert weight < 0
    for degree in range(-weight+1):
        global_columns.append(Tmod**degree*popov.column(generator))
        column_labels.append((generator, degree))
assert len(global_columns) == 217

# D=A-43F at T=0: kill all six numerator coefficients through T^42.
entries = {}
for column_index, column in enumerate(global_columns):
    for row, polynomial in enumerate(column):
        for degree in range(43):
            value = polynomial[degree]
            if value:
                entries[(43*row+degree, column_index)] = value
fiber_gate = matrix(field, 6*43, len(global_columns), entries, sparse=True)
rank = fiber_gate.rank()
kernel = fiber_gate.right_kernel_matrix()
assert fiber_gate.dimensions() == (258, 217)
assert rank == 211 and kernel.nrows() == 6

print(
    f"CM43Q60SMITH|p={p}|invariant_degrees={invariant_degrees}"
    f"|generator_weights={generator_weights}|pivots={pivots}",
    flush=True,
)
print(
    f"CM43Q60SMITH|gate_dimensions={fiber_gate.dimensions()}|rank={rank}"
    f"|nullity={kernel.nrows()}|interpretation=CM_boundary_degeneration",
    flush=True,
)
print("CM43Q60SMITH|status=BOUNDARY_SIX_DIMENSIONAL_KERNEL", flush=True)
