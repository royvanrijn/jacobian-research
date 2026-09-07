#!/usr/bin/env sage
"""Recover the Humbert-8 summand intrinsically in the rootless rank-17 NS.

This is a verifier, not a search.  It checks one short primitive class in
``U + (-rank17_gram)`` and identifies its orthogonal complement with the
indefinite genus of the fixed part of the marked Kumar H2 frame.  It also
proves that a q=8 elliptic-fiber presentation cannot be orthogonal to this
Humbert-8 class.
"""

from pathlib import Path

from sage.all import (
    QQ,
    ZZ,
    Genus,
    block_diagonal_matrix,
    gcd,
    identity_matrix,
    matrix,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
U = matrix(ZZ, [[0, 1], [1, 0]])


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


rank17 = load_matrix(ROOT / "elkies-k3/data/lattice/rank17_gram.txt")
kumar = load_matrix(
    ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_2.txt"
)
ns = block_diagonal_matrix(U, -rank17)

# The frame part belongs to {x : rank17*x is divisible by 4}.  Its norm 36
# makes (4,4,x) have square -4 in U+(-rank17).
x = vector(ZZ, (-1, 0, -3, 0, 2, -2, 1, -2, 1, 1, 0, 1, 0, 0, -2, -2, 2))
h = vector(ZZ, (4, 4) + tuple(x))
assert x * rank17 * x == 36
assert all(value % 4 == 0 for value in rank17 * x)
assert h * ns * h == -4
assert gcd([abs(ZZ(value)) for value in ns * h]) == 4
assert gcd([abs(value) for value in h]) == 1

# Reflection in h: v |-> v + (v.h)h/2.  Divisibility four makes it integral.
reflection = identity_matrix(QQ, ns.nrows())
reflection += (ns * h.column() / 2) * matrix(QQ, [h])
assert all(value.denominator() == 1 for value in reflection.list())
reflection = reflection.change_ring(ZZ)
assert reflection * reflection == identity_matrix(ZZ, ns.nrows())
assert reflection * ns * reflection.transpose() == ns
assert reflection.det() == -1

# In Smith coordinates the cyclic discriminant group has order 948.  The
# reflection acts by 475: -1 on the 4-part, +1 on the 3- and 79-parts.
smith, left, _ = ns.smith_form()
assert tuple(smith.diagonal()) == (1,) * 18 + (948,)
smith_action = left * reflection * left.inverse()
assert all(smith_action[-1, index] % 948 == 0 for index in range(18))
multiplier = ZZ(smith_action[-1, -1]) % 948
assert multiplier == 475
assert (multiplier % 4, multiplier % 3, multiplier % 79) == (3, 1, 1)

# Since div(h)=|h^2|=4, h splits off integrally.  The complement is the
# determinant-237 fixed lattice expected from H237.
orthogonal_basis = matrix(ZZ, [list(h * ns)]).right_kernel_matrix()
orthogonal_gram = orthogonal_basis * ns * orthogonal_basis.transpose()
split_basis = matrix(ZZ, [list(h)] + orthogonal_basis.rows())
assert abs(split_basis.det()) == 1
assert split_basis * ns * split_basis.transpose() == block_diagonal_matrix(
    matrix(ZZ, [[-4]]), orthogonal_gram
)
assert orthogonal_gram.det() == -237
orthogonal_smith = tuple(orthogonal_gram.smith_form()[0].diagonal())
assert orthogonal_smith == (1,) * 17 + (237,)

# Delete the orthogonal height-four row from the Kumar H2 frame.  Including
# its hyperbolic plane gives the other rank-18 determinant-237 lattice.
fixed_indices = list(range(15)) + [16]
kumar_fixed_frame = kumar.matrix_from_rows_and_columns(
    fixed_indices, fixed_indices
)
kumar_fixed_ns = block_diagonal_matrix(U, -kumar_fixed_frame)
assert kumar_fixed_ns.det() == -237
assert Genus(orthogonal_gram) == Genus(kumar_fixed_ns)

# There is no q=8 fiber in h^perp.  If f=(a,b,v), then ab=8,
# v*R*v=16 and ell.v=a+b for ell=x*R/4.  Cauchy--Schwarz gives
# v*R*v >= (ell.v)^2/(ell*R^-1*ell).  The two factor types have sums 9
# and 6.  Sum 9 is already impossible; sum 6 forces equality and hence
# v=2*x/3, which is not integral.
ell = vector(QQ, x * rank17 / 4)
dual_norm = ell * rank17.inverse() * ell
assert dual_norm == QQ(9) / 4
assert QQ(9) ** 2 / dual_norm > 16
assert QQ(6) ** 2 / dual_norm == 16
extremal = QQ(6) / dual_norm * ell * rank17.inverse()
assert extremal == QQ(2) / 3 * x
assert any(value.denominator() != 1 for value in extremal)

print(
    "R17H8SPLIT|h={}|norm=-4|divisibility=4|reflection_multiplier={}".format(
        tuple(h), multiplier
    ),
    flush=True,
)
print(
    "R17H8SPLIT|complement_rank=18|det=-237|smith=1^17,237|"
    "genus=KUMAR_H2_FIXED",
    flush=True,
)
print(
    "R17H8SPLIT|q8_orthogonal_fiber=IMPOSSIBLE|"
    "factor_sums=6,9|status=PASS",
    flush=True,
)
