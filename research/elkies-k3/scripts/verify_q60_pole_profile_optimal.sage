#!/usr/bin/env sage
"""Prove the 0,2,4 section-pole profile of the q=60 frame is optimal."""

from sage.all import Matrix, QQ, QuadraticForm, ZZ, matrix, vector


height = Matrix(QQ, [[4, 0, 0], [0, QQ(20) / 3, 1], [0, 1, 12]])
scaled = (3 * height).change_ring(ZZ)


def correction(v):
    """E6 local correction; only the middle generator has nonzero class."""
    return QQ(0) if ZZ(v[1]) % 3 == 0 else QQ(4) / 3


def zero_intersection(v):
    h = vector(QQ, v) * height * vector(QQ, v)
    return (h - 4 + correction(v)) / 2


# If P.O <= 2, then h(P) <= 8.  In the integral form 3H this is norm <= 24,
# or quadratic-form value <= 12.  Fincke--Pohst enumeration is exhaustive.
form = QuadraticForm(ZZ, scaled)
shells = form.short_vector_list_up_to_length(13, up_to_sign_flag=True)
short = []
for shell in shells:
    for row in shell:
        v = vector(ZZ, row)
        if v and zero_intersection(v) <= 2:
            short.append(v)

short_module = matrix(ZZ, [list(v) for v in short]).row_module()
assert short_module.rank() == 2
assert all(v[2] == 0 for v in short)

# The reduced basis attains P.O = 0,2,4, so the lower bound is sharp.
basis = matrix.identity(ZZ, 3)
profile = tuple(ZZ(zero_intersection(row)) for row in basis.rows())
assert abs(basis.det()) == 1
assert profile == (0, 2, 4)

print(
    "Q60POLES|vectors_with_PO_le_2={}|span_rank=2|"
    "every_integral_basis_max_PO_ge_4".format(2 * len(short))
)
print("Q60POLES|attained_profile=0,2,4|status=PASS")
