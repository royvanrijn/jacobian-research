"""Check the finite coefficient claims in the doubled-square j exclusion.

From the repository root run this file with SymPy and Singular on PATH.
The proof uses each calculation only over characteristic zero.
"""

from fractions import Fraction
from itertools import product
import subprocess

from sympy import Matrix, Poly, diff, expand, symbols


x, y, z, Y, P = symbols("x y z Y P")
a, b, c, d, e, f, g, h, i, ell = symbols("a b c d e f g h i ell")

INITIAL = (
    {(4, 0, 2), (3, 2, 1), (2, 4, 0)},
    {(5, 0, 1), (4, 2, 0)},
    {(6, 0, 0), (4, 0, 1), (3, 2, 0)},
)
AXIAL = ((4, 0, 2), (5, 0, 1), (4, 0, 1))


def weight(monomial, weights):
    return sum(exponent * w for exponent, w in zip(monomial, weights))


def cut(supports, weights, highest):
    return tuple(
        {m for m in allowed if weight(m, weights) < weight(next(iter(top)), weights)}
        | top
        for allowed, top in zip(supports, highest)
    )


def competing_monomials(supports, base, axis, lower, upper, reference):
    """List ties in the indicated interval of one varying weight."""
    answer = []
    for allowed, ref in zip(supports, reference):
        ties = {}
        for m in allowed - {ref}:
            offset = weight(m, base) - weight(ref, base)
            slope = m[axis] - ref[axis]
            if slope:
                tie = Fraction(-offset, slope)
                if lower < tie <= upper:
                    ties.setdefault(tie, set()).add(m)
            else:
                assert offset <= 0, (m, ref, base)
        answer.append(ties)
    return tuple(answer)


def check_support():
    monomials = {m for m in product(range(7), repeat=3) if sum(m) <= 6}
    supports = cut((monomials,) * 3, (2, 3, 4), INITIAL)
    assert all(
        weight(m, (1, 1, 2)) <= weight(ref, (1, 1, 2))
        for allowed, ref in zip(supports, AXIAL) for m in allowed
    )
    assert competing_monomials(supports, (1, 1, 0), 2, 2, 4, AXIAL) == (
        {3: {(1, 0, 3), (0, 1, 3)}, 4: {(0, 0, 3)}},
        {Fraction(5, 2): {(0, 0, 3)}, 3: {(2, 0, 2), (1, 1, 2)},
         4: {(1, 0, 2), (0, 1, 2)}},
        {3: {(1, 0, 2), (0, 1, 2)}, 4: {(0, 0, 2)}},
    )
    at_three = (
        INITIAL[0],
        INITIAL[1] | {(2, 0, 2), (1, 2, 1), (0, 4, 0)},
        INITIAL[2] - {(6, 0, 0)},
    )
    assert competing_monomials(
        cut(supports, (1, 2, 3), at_three), (0, 2, 3), 0, 0,
        1, ((2, 4, 0), (0, 4, 0), (3, 2, 0)),
    ) == (
        {1: {(3, 2, 1), (4, 0, 2)}, Fraction(1, 2): {(0, 0, 3), (0, 3, 1)}},
        {1: {(1, 2, 1), (2, 0, 2), (4, 2, 0), (5, 0, 1)}},
        {1: {(4, 0, 1)}, Fraction(2, 3): {(0, 0, 2), (0, 3, 0)},
         Fraction(1, 2): {(1, 1, 1)}, Fraction(1, 3): {(0, 1, 1)}},
    )
    removed = (
        {(1, 0, 3), (0, 1, 3)},
        {(0, 0, 3), (2, 0, 2), (1, 1, 2), (0, 1, 2)},
        {(1, 0, 2), (0, 1, 2)},
    )
    axial_four = tuple(allowed - excluded for allowed, excluded in zip(supports, removed))
    assert all(
        weight(m, (1, 1, 4)) <= weight(ref, (1, 1, 4))
        for allowed, ref in zip(axial_four, AXIAL) for m in allowed
    )
    assert competing_monomials(
        axial_four, (1, 0, 4), 1, 1, Fraction(5, 2), AXIAL,
    ) == (
        {2: {(0, 2, 2), (2, 1, 2)}, Fraction(7, 3): {(1, 3, 1)},
         Fraction(12, 5): {(0, 5, 0)}, Fraction(5, 2): {(2, 4, 0), (3, 2, 1)}},
        {Fraction(5, 3): {(0, 3, 1)}, 2: {(1, 2, 1), (3, 1, 1)},
         Fraction(9, 4): {(0, 4, 0)}, Fraction(7, 3): {(2, 3, 0)},
         Fraction(5, 2): {(0, 2, 1), (4, 2, 0)}},
        {2: {(0, 2, 1), (2, 1, 1)}, Fraction(7, 3): {(1, 3, 0)},
         Fraction(5, 2): {(3, 2, 0)}},
    )
    final = cut(axial_four, (2, 5, 8), (INITIAL[0], INITIAL[1], at_three[2]))
    assert tuple({m for m in allowed if m[2] > ref[2]}
                 for allowed, ref in zip(final, AXIAL)) == (set(), {(0, 0, 2)}, set())


def jacobian(outputs, variables):
    return Poly(Matrix(outputs).jacobian(variables).det(), *variables)


def check_determinants():
    w_three = jacobian(
        (x**4*z**2 + (a*x+b*y)*z**3,
         x**5*z + (c*x**2+d*x*y)*z**2,
         x**4*z + (e*x+f*y)*z**2), (x, y, z),
    )
    selected = [(8, 0, 4), (4, 1, 5), (0, 2, 6)]
    assert [w_three.coeff_monomial(mon) for mon in selected] == [
        -b-4*d+6*f, -5*b*d+5*b*f+2*d*f, b*d*f,
    ]

    first = x**2*y**4+a*z**3+b*y**3*z
    third = x**3*y**2+c*x*y*z
    half = Poly(diff(first, x)*diff(third, z)-diff(first, z)*diff(third, x), z)
    assert [expand(half.nth(n) - expected) for n, expected in zip(
        (2, 1, 0), (-9*a*x**2*y**2, -b*c*y**4, x**2*y**5*(-3*b+2*c)),
    )] == [0, 0, 0]

    seven_thirds = jacobian(
        (x**4*z**2+a*z**3+d*x*Y*z,
         x**5*z+b*x*z**2+e*x**2*Y,
         x**4*z+c*z**2+f*x*Y), (x, Y, z),
    )
    selected = [(9, 0, 2), (6, 1, 1), (3, 2, 0)]
    assert [seven_thirds.coeff_monomial(mon) for mon in selected] == [
        -d-4*e+6*f, -5*d*e+5*d*f+2*e*f, d*e*f,
    ]

    five_halves = jacobian(
        (x**2*P**2+a*z**3, x**4*P+b*x*z**2+g*P*z, x**3*P+c*z**2), (x, P, z),
    )
    selected = [(6, 1, 2), (5, 2, 1), (4, 3, 0), (2, 1, 3)]
    assert [
        expand(five_halves.coeff_monomial(mon) - expected) for mon, expected in
        zip(selected, (3*a, 4*(2*b-3*c), 4*g, -9*a*g-4*b*c))
    ] == [0]*4


def check_weight_two_radical():
    determinant = jacobian(
        (z**2*(x**4+d*y**2+f*x**2*y+a*z),
         x*z*(x**4+g*y**2+h*x**2*y+b*z),
         z*(x**4+i*y**2+ell*x**2*y+c*z)), (x, y, z),
    )
    first = (d, f, g, h, i, ell)
    second = (a, b, c, g, i, ell-h, f-2*h, d-h**2)
    singular = [
        'LIB "primdec.lib";',
        'ring R=0,(a,b,c,d,f,g,h,i,ell),dp;',
        'ideal I=' + ','.join(str(coef).replace('**', '^') for coef in determinant.coeffs()) + ';',
        'ideal G=std(radical(I));',
    ]
    for left, right in product(first, second):
        expression = str(expand(left*right)).replace('**', '^')
        singular.append(f'if(reduce({expression},G)!=0){{print("Nonzero reduction");}}')
    singular.extend(('print("All products vanish modulo the radical");', 'quit;'))
    run = subprocess.run(
        ['Singular', '-q'], input='\n'.join(singular)+'\n', check=True,
        capture_output=True, text=True,
    )
    assert run.stdout.strip() == 'All products vanish modulo the radical', run.stdout
    assert not run.stderr, run.stderr


if __name__ == "__main__":
    check_support()
    check_determinants()
    check_weight_two_radical()
    print("Doubled-square third mixed correction: support, determinants and radical verified.")
