#!/usr/bin/env sage
"""Derive the resolved D7 rows for the CM24 third q12 pencil.

The unique effective lattice lift gives

    D = Q + 2*O + 4*F + R,

with one selected component from the quadratic I2 pair, the full affine E6
vector, and the displayed D7 correction at infinity.  The ambient coefficient
space is

    a(W) + b*X + c(W)*z_Q,
    deg(a)<=4, deg(c)<=2,

of dimension nine.  This verifier constructs the CM compositum, checks the
selected I2 node and the two singular-point values of z_Q, and emits the seven
exact component conditions:

* c(r)=0 and a(r)+b*Qx(r)=0 at one root r of the quadratic I2 factor;
* a(-27/2)+162*s*c(-27/2)=0 at the IV* cusp;
* four rows at infinity obtained from the complete D7 ideal.

For the canonical rational double point

    Y^2 + U^2*Z + Z^6 = 0,

the correction cycle has coefficients ``(3,2,5,4,3,3,6)`` in the resolution
graph order used below.  These are exactly the divisorial valuations of Y.
The corresponding complete ideal is

    (Y, U^2, Z*U, Z^3),

and its quotient has basis ``1,Z,Z^2,U``.  In the actual infinity chart the
3-jet is ``y^2-u^2(u+9t)`` and there is no degree-four term, so modulo this
length-four quotient one may take ``Z=-(u+9t)``, ``U=u``.  Higher terms in the
analytic D7 coordinate change lie in the complete ideal.  Reducing the nine
ambient sections in this quotient gives the required four resolved rows.

The conjugate choice of the quadratic I2 component is equivalent.  The
resulting 7-by-9 matrix has rank seven and a two-dimensional kernel.
"""

from sage.all import (
    Matrix, PolynomialRing, QQ, QuadraticField, vector,
)


K = QuadraticField(-6, "s")
s = K.gen()
extension_ring = PolynomialRing(K, "J")
J = extension_ring.gen()
L = K.extension(J**2+3, "j")
j = L.gen()
s = L(s)

base = PolynomialRing(L, "W")
W = base.gen()
A = (
    -27*W**6+59049*W**4+L(13286025)/8*W**3
    + L(129140163)/8*W**2+L(1162261467)/32*W
    - L(10460353203)/64
)
B = (
    54*W**9-177147*W**7-L(97253703)/8*W**6
    - L(7360989291)/16*W**5-L(331244518095)/32*W**4
    - L(4487491524087)/32*W**3-L(144886352214753)/128*W**2
    - L(1303977169932777)/256*W-L(5147278302366225)/512
)
Qx = (
    -L(8)/27*W**4+22*W**3-L(243)/2*W**2+729*W
    - L(492075)/8
)
Qy = s*(
    L(16)/243*W**6-L(22)/3*W**5+L(333)/2*W**4
    - L(2025)/4*W**3+L(190269)/4*W**2
    - L(177147)/16*W+L(199290375)/32
)
assert Qy**2 == Qx**3+A*Qx+B

quadratic_i2 = W**2+L(27)/2*W+L(729)/4
selected_root = L(27)/4*(-1+j)
conjugate_root = L(27)/4*(-1-j)
assert quadratic_i2(selected_root) == quadratic_i2(conjugate_root) == 0


def verify_i2_node(root):
    local_A = A(root)
    local_B = B(root)
    node_x = -3*local_B/(2*local_A)
    assert Qx(root) == node_x
    assert Qy(root) == 0
    assert node_x**3+local_A*node_x+local_B == 0
    assert 3*node_x**2+local_A == 0


verify_i2_node(selected_root)
verify_i2_node(conjugate_root)

e6_point = -L(27)/2
assert A(e6_point) == B(e6_point) == 0
e6_z_value = -Qy(e6_point)/Qx(e6_point)
assert e6_z_value == 162*s

# At infinity use t=1/W, xbar=t^4*X, ybar=t^6*Y and zbar=t^2*z.
# The cusp is xbar=ybar=0.  The leading Q coordinates give zbar=2*s/9.
qx_infinity = Qx[4]
qy_infinity = Qy[6]
infinity_z_value = -qy_infinity/qx_infinity
assert infinity_z_value == L(2)/9*s


def evaluation_row(point, include_b_qx=False, c_scale=0):
    row = [point**degree for degree in range(5)]
    row.append(Qx(point) if include_b_qx else L(0))
    row.extend(c_scale*point**degree for degree in range(3))
    return vector(L, row)


def finite_rows(root):
    c_at_root = vector(
        L, [0, 0, 0, 0, 0, 0, 1, root, root**2]
    )
    m_at_root = evaluation_row(root, include_b_qx=True)
    e6_cusp = evaluation_row(e6_point, c_scale=e6_z_value)
    return [c_at_root, m_at_root, e6_cusp]


# Certify the tangent D7 coordinates directly from the local equation.  Here
# t=1/W, xbar=t^4*X, ybar=t^6*Y, and u=xbar-3*t.
local_ring = PolynomialRing(L, names=("tloc", "uloc", "yloc"))
tloc, uloc, yloc = local_ring.gens()
local_A = sum(L(A[index])*tloc**(8-index) for index in range(7))
local_B = sum(L(B[index])*tloc**(12-index) for index in range(10))
local_equation = (
    yloc**2-(uloc+3*tloc)**3-local_A*(uloc+3*tloc)-local_B
)
homogeneous_parts = local_equation.homogeneous_components()
assert homogeneous_parts[2] == yloc**2
assert homogeneous_parts[3] == -uloc**3-9*tloc*uloc**2
assert homogeneous_parts.get(4, local_ring(0)) == 0

# Canonical D7 graph order: the long arm is 2-1-4-3-7 and the two short
# arms are 5-7 and 6-7.  The coordinate valuations come from the explicit
# minimal resolution of Y^2+U^2*Z+Z^6.
d7_intersection = Matrix(
    QQ,
    [
        [2, -1, 0, -1, 0, 0, 0],
        [-1, 2, 0, 0, 0, 0, 0],
        [0, 0, 2, -1, 0, 0, -1],
        [-1, 0, -1, 2, 0, 0, 0],
        [0, 0, 0, 0, 2, 0, -1],
        [0, 0, 0, 0, 0, 2, -1],
        [0, 0, -1, 0, -1, -1, 2],
    ],
)
d7_correction = vector(QQ, (3, 2, 5, 4, 3, 3, 6))
assert d7_intersection*d7_correction == vector(QQ, (0, 1, 0, 0, 0, 0, 1))
d7_valuations = {
    "Z": vector(QQ, (2, 2, 2, 2, 1, 1, 2)),
    "U": vector(QQ, (2, 1, 4, 3, 3, 3, 5)),
    "Y": vector(QQ, (3, 2, 5, 4, 3, 3, 6)),
}
assert d7_valuations["Y"] == d7_correction
for valuation in (
    d7_valuations["Y"],
    2*d7_valuations["U"],
    d7_valuations["Z"]+d7_valuations["U"],
    3*d7_valuations["Z"],
):
    assert min(valuation-d7_correction) >= 0

# Work only to length four.  The allowed residues are 1,Z,Z^2,U; all other
# monomials lie in (Y,U^2,ZU,Z^3).  Since the coordinate change has no
# quadratic correction, t=-(Z+U)/9 and xbar=(2U-Z)/3 in this quotient.
residue_ring = PolynomialRing(L, names=("Z", "U"))
Z, U = residue_ring.gens()
local_t = -(Z+U)/9
local_x = (2*U-Z)/3
allowed_exponents = ((0, 0), (1, 0), (2, 0), (0, 1))


def residue_vector(polynomial):
    coefficients = residue_ring(polynomial).dict()
    return vector(
        L,
        [
            sum(
                (coefficient for exponent, coefficient in coefficients.items()
                 if tuple(exponent) == allowed),
                L(0),
            )
            for allowed in allowed_exponents
        ],
    )


# zbar=t^2*z_Q=(ybar+Qybar)/(xbar-Qxbar).  Its denominator is a unit at the
# D7 point.  Expand its ybar=0 restriction through degree two.
qx4, qx3, qx2 = (L(Qx[index]) for index in (4, 3, 2))
qy6, qy5, qy4 = (L(Qy[index]) for index in (6, 5, 4))
denominator_constant = -qx4
denominator_linear = local_x-qx3*local_t
denominator_quadratic = -qx2*local_t**2
inverse_denominator = (
    1/denominator_constant
    - denominator_linear/denominator_constant**2
    + denominator_linear**2/denominator_constant**3
    - denominator_quadratic/denominator_constant**2
)
zbar_residue = residue_ring(
    (qy6+qy5*local_t+qy4*local_t**2)*inverse_denominator
)
assert residue_vector(zbar_residue)[0] == infinity_z_value

# Coefficient order: a0,...,a4,b,c0,c1,c2.  After multiplication by t^4
# these nine local sections are t^4,...,1,xbar,t^2*zbar,t*zbar,zbar.
local_ambient_basis = (
    [local_t**(4-index) for index in range(5)]
    + [local_x]
    + [local_t**(2-index)*zbar_residue for index in range(3)]
)
infinity_rows = list(
    Matrix(L, [residue_vector(value) for value in local_ambient_basis]).transpose().rows()
)
assert Matrix(L, infinity_rows).rank() == 4

rows = Matrix(L, finite_rows(selected_root)+infinity_rows)
conjugate_rows = Matrix(L, finite_rows(conjugate_root)+infinity_rows)
assert rows.rank() == conjugate_rows.rank() == 7
assert rows.ncols()-rows.rank() == 2
kernel = rows.right_kernel().basis_matrix()
assert kernel.nrows() == 2 and kernel*rows.transpose() == 0

print(
    f"Q80THIRDLOCALGATES|field=QQ(sqrt(-6),sqrt(-3))|"
    f"selected_I2_root={selected_root}|conjugate={conjugate_root}|"
    "ambient=1,W,W2,W3,W4,X,z,Wz,W2z",
    flush=True,
)
print(
    "Q80THIRDLOCALGATES|D7_normal_form=Y2+U2Z+Z6|"
    "D7_complete_ideal=(Y,U2,ZU,Z3)|D7_quotient_basis=1,Z,Z2,U|"
    f"D7_residue_rows={tuple(tuple(row) for row in infinity_rows)}",
    flush=True,
)
print(
    f"Q80THIRDLOCALGATES|rows={tuple(tuple(row) for row in rows.rows())}|"
    f"rank={rows.rank()}|nullity={rows.right_kernel().dimension()}|"
    f"kernel={tuple(tuple(row) for row in kernel.rows())}|"
    "status=PASS_RESOLVED_D7_KERNEL",
    flush=True,
)
