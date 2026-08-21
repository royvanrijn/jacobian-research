#!/usr/bin/env sage
"""Certify that the short CM-43 q=8 fiber is not a generic NS class.

The generic determinant-948 Kumar frame is enlarged at CM discriminant -43
by a rank-one vector of norm 40764 and the marked primitive-closure glue
``k=211``.  In the marked closure the q=8 witness has nonzero coordinate
along that added rank-one factor, whereas the q=60 witness does not.

The same obstruction is visible in the explicit Mordell--Weil quotient.  If
``P3=(0,0,1)`` and ``Q79=(4,-5,1)``, then the horizontal part
``R=P1-2*P2`` of q=8 has a nonzero projection onto the primitive CM-only
direction ``W=(116,92,29)``.  Consequently q=8 is a CM-boundary
factorization/normal-form device, not a fibration class on the generic
rank-19 family.  The q=60 fiber, whose horizontal part is Q79, remains the
first generic neighbor.
"""

from pathlib import Path

from sage.all import *


BASE = Path(__file__).resolve().parents[1]


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ],
    )


def discriminant_generator(gram):
    smith, left, _ = gram.smith_form()
    diagonal = list(smith.diagonal())
    assert diagonal[:-1] == [1] * (gram.nrows() - 1)
    order = abs(ZZ(diagonal[-1]))
    last = vector(ZZ, [0] * (gram.nrows() - 1) + [1])
    integral_covector = left.inverse() * last
    dual_vector = gram.inverse() * integral_covector
    assert order * dual_vector in ZZ**gram.nrows()
    return order, dual_vector


def primitive_closure(gram, norm, disc_coefficient, rank1_coefficient):
    _, generator = discriminant_generator(gram)
    ambient = block_diagonal_matrix(gram, matrix(ZZ, [[norm]]))
    glue = vector(
        QQ,
        list(disc_coefficient * generator)
        + [QQ(rank1_coefficient) / norm],
    )
    assert (glue * ambient * glue) % 2 == 0
    denominator = ZZ(glue.denominator())
    generators = denominator * identity_matrix(ZZ, gram.nrows() + 1)
    generators = generators.stack(
        matrix(ZZ, 1, gram.nrows() + 1, [
            ZZ(denominator * entry) for entry in glue
        ])
    )
    basis = generators.row_module().basis_matrix().change_ring(QQ) / denominator
    extension = basis * ambient * basis.transpose()
    assert all(entry.denominator() == 1 for entry in extension.list())
    return extension.change_ring(ZZ), basis


generic_frame = load_gram(
    BASE / "data/fibrations/kumar_e7e8_mw2_frame_2.txt"
)
marked_frame = load_gram(
    BASE / "data/fibrations/kumar_cm43_marked_e7e8_mw3_frame.txt"
)
assert generic_frame.det() == 948 and marked_frame.det() == 43

# The rows of closure_basis express the marked closure basis in the orthogonal
# ambient basis [generic frame(17), CM rank-one vector].
closure, closure_basis = primitive_closure(generic_frame, 40764, 211, 43)
assert closure == marked_frame
assert ZZ(1 / abs(closure_basis.det())) == 948

q8_frame = vector(ZZ, (
    156, -78, 0, 0, -78, 0, -78, 0, 0,
    0, 0, 0, 0, 0, 0, -1, -155, -32,
))
q8_ambient = q8_frame * closure_basis
expected_q8_ambient = vector(QQ, (
    QQ(52)/79, QQ(78)/79, QQ(104)/79, QQ(156)/79,
    QQ(130)/79, QQ(104)/79, QQ(78)/79,
    0, 0, 0, 0, 0, 0, 0, 0,
    -QQ(1)/4, QQ(27)/79, QQ(1)/316,
))
assert q8_ambient == expected_q8_ambient
assert q8_ambient[-1] != 0

q60_frame = vector(ZZ, (
    0, 0, -1, -1, -1, -1, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0,
))
q60_ambient = q60_frame * closure_basis
assert q60_ambient == q60_frame
assert q60_ambient[-1] == 0

# Independent check in the explicit CM-43 Mordell--Weil quotient.
height = matrix(QQ, (
    (QQ(5)/2, -QQ(1)/2, -1),
    (-QQ(1)/2, QQ(5)/2, 0),
    (-1, 0, 4),
))
p3 = vector(QQ, (0, 0, 1))
q79 = vector(QQ, (4, -5, 1))
r = vector(QQ, (1, -2, 0))
assert matrix(QQ, (p3, q79, r)).det() == -3

generic_rows = matrix(QQ, (p3, q79))
generic_gram = generic_rows * height * generic_rows.transpose()
projection_coefficients = (
    vector(QQ, (r*height*p3, r*height*q79)) * generic_gram.inverse()
)
assert projection_coefficients == vector(QQ, (-QQ(1)/4, QQ(27)/79))
residual = r - projection_coefficients * generic_rows
primitive_cm_direction = vector(QQ, (116, 92, 29))
assert primitive_cm_direction*height*p3 == 0
assert primitive_cm_direction*height*q79 == 0
assert primitive_cm_direction*height*primitive_cm_direction == 40764
assert residual == -primitive_cm_direction/316
assert residual*height*residual == QQ(129)/316 == 9*QQ(43)/948

print(
    "CM43Q8MEMBERSHIP|q8_cm_coordinate={}|generic=0|"
    "q60_cm_coordinate={}|generic=1".format(q8_ambient[-1], q60_ambient[-1]),
    flush=True,
)
print(
    "CM43Q8MEMBERSHIP|R=-1/4*P3+27/79*Q79-1/316*W|"
    "W=(116,92,29)|height_W=40764|residual_height=129/316",
    flush=True,
)
print("CM43Q8MEMBERSHIP|status=PASS", flush=True)
