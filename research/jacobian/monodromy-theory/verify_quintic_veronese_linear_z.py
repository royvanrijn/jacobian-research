"""Replay the determinant identities and resonance in the linear-in-z reduction."""

from sympy import Function, Matrix, expand, symbols


def zero(expression):
    assert expand(expression) == 0, expression


def check():
    x, y, a, b = symbols("x y a b")
    q0, q1, q2 = symbols("q0 q1 q2")
    N = Matrix([x**2, x*y, y**2])
    nu = Matrix([y**2, -2*x*y, x**2])
    q = q0*x**2 + q1*x*y + q2*y**2

    A = q*N
    B = a*N.diff(x) + b*N.diff(y)
    mixed = (
        3*A.diff(x).dot(B.diff(y).cross(B))
        + 3*B.diff(x).dot(A.diff(y).cross(B))
        + B.diff(x).dot(B.diff(y).cross(A))
    )
    zero(mixed + 22*q*(b*x-a*y)**2)

    r = symbols("r")
    M = Matrix([Function(f"M{i}")(x, y) for i in range(3)])
    mixed = (
        (6-r)*N.diff(x).dot(N.diff(y).cross(M))
        + 6*M.diff(x).dot(N.diff(y).cross(N))
        + 6*N.diff(x).dot(M.diff(y).cross(N))
    )
    zero(mixed - ((12-2*r)*nu.dot(M) - 6*nu.dot(x*M.diff(x)+y*M.diff(y))))
    assert [(r, j) for r in range(1, 6) for j in range(3) if 12-2*r-6*j == 0] == [(3, 1)]

    b1, b2, b3, d1, d2, d3, s = symbols("b1 b2 b3 d1 d2 d3 s")
    H = N + x*Matrix([b1,b2,b3]) + y*Matrix([d1,d2,d3])
    combination = (H.diff(x) + s*H.diff(y)).subs(
        {x: -(b1+s*d1)/2, y: -(b3/s+d3)/2}
    )
    zero(combination[0])
    zero(combination[2])
    zero(2*s*combination[1] - (-d1*s**3+(2*d2-b1)*s**2+(2*b2-d3)*s-b3))


if __name__ == "__main__":
    check()
    print("Verified the weighted reduction determinant identities and resonance.")
