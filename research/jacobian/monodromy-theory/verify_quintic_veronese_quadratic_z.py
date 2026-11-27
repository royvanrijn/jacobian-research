"""Check the polynomial identities in QuinticVeroneseQuadraticZExclusion."""

from sympy import Function, Matrix, expand, symbols


def zero(expression):
    assert expand(expression) == 0, expression


def zero_matrix(matrix):
    for entry in matrix:
        zero(entry)


def check():
    x, y, z, a, b, l1, l2, c = symbols("x y z a b l1 l2 c")
    N = Matrix([x**2, x*y, y**2])
    nu = Matrix([y**2, -2*x*y, x**2])
    ell = l1*x + l2*y

    g_coeffs = symbols("g0:7")
    g = ell*z**2 + sum(g_coeffs[i]*term for i, term in enumerate(
        (z*x**2, z*x*y, z*y**2, x**3, x**2*y, x*y**2, y**3)
    ))
    H5 = g*N
    jac = H5.jacobian((x, y, z))
    zero_matrix(
        jac.adjugate()
        - g*(Matrix([0, 0, 5*g])-g.diff(z)*Matrix([x, y, z]))*nu.T
    )

    A = ell*N
    B = a*N.diff(x) + b*N.diff(y)
    mixed = (
        3*A.diff(x).dot(B.diff(y).cross(B))
        + 3*B.diff(x).dot(A.diff(y).cross(B))
        + 2*B.diff(x).dot(B.diff(y).cross(A))
    )
    zero(mixed + 14*ell*(b*x-a*y)**2)

    h = ell + c
    A2 = h*N
    L = Matrix([Function(f"L{i}")(x, y) for i in range(3)])
    S = nu.dot(L)
    fourth_power = (
        A2.diff(x).dot(A2.diff(y).cross(L))
        + 2*L.diff(x).dot(A2.diff(y).cross(A2))
        + 2*A2.diff(x).dot(L.diff(y).cross(A2))
    )
    zero(fourth_power - h*((7*ell+6*c)*S-2*h*(x*S.diff(x)+y*S.diff(y))))


if __name__ == "__main__":
    check()
    print("Verified the leading adjugate and both mixed determinant identities.")
