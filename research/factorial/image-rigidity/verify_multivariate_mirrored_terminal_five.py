from math import comb

from sympy import Poly, symbols


def verify_corner(v_num, v_den, y_num, y_den, u, w, y):
    d = 1 + 4*u - 5*w
    common_den = v_den*y_den
    c = (1-5*w)*common_den+10*y_num*v_den-6*v_num*y_den
    numerator = Poly(
        c**2+4*u*c*d*common_den+(6*v_num-2*v_den)*v_den*y_den**2*d**2,
        w,
    )

    # Projectivize the interval for w, with y in [0, +infinity].
    low_num, low_den = u**2*(u+5), 5*u+1
    high_num, high_den = 15*u, 16-u
    w_num = low_num*high_den+y*high_num*low_den
    w_den = (1+y)*low_den*high_den
    m = numerator.degree()
    poly = Poly(sum(
        coefficient*w_num**j*w_den**(m-j)
        for (j,), coefficient in numerator.terms()
    ), u, y)
    reduced = poly.exquo(Poly((1-u)**4, u, y))

    # Homogenize u=x/(1+x) via the binomial theorem.
    degree = reduced.degree(u)
    y_degree = reduced.degree(y)
    coefficients = {
        (i, j): 0 for i in range(degree+1) for j in range(y_degree+1)
    }
    for (i, j), coefficient in reduced.terms():
        for k in range(degree-i+1):
            coefficients[i+k, j] += coefficient*comb(degree-i, k)
    if not coefficients or max(coefficients.values()) >= 0:
        raise ArithmeticError("A transformed coefficient is nonnegative")
    return len(coefficients)


def main():
    u, w, y = symbols("u w y")
    v_bounds = ((u**2*(u+5), 5*u+1), (u*(1+2*u), 2+u))
    y_bounds = ((w**2*(w+3), 3*w+1), (w*(1+3*w), 3+w))
    for i, (v_num, v_den) in enumerate(v_bounds):
        for j, (y_num, y_den) in enumerate(y_bounds):
            print(f"corner {i},{j}: {verify_corner(v_num, v_den, y_num, y_den, u, w, y)} strictly negative coefficients")


if __name__ == "__main__":
    main()
