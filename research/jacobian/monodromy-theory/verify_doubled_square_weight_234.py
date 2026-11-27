"""Replay the doubled-square-member weight-(2,3,4) coefficient reduction."""

from fractions import Fraction
from subprocess import run

from sympy import Matrix, Poly, expand, symbols


def tie_terms():
    ordinary = [
        {(4, 0, 2), (3, 2, 1), (2, 4, 0)},
        {(5, 0, 1), (4, 2, 0)},
        {(6, 0, 0)},
    ]
    ties = {}
    for component, target_weight in enumerate((8, 7, 6)):
        for j in range(1, 7):
            for i in range(7 - j):
                for l in range(7 - i - j):
                    monomial = (i, j, l)
                    if monomial not in ordinary[component]:
                        if i + j + l == 6 or 3 * (i + j) + 5 * l >= (22, 20, 18)[component]:
                            continue
                    weight = Fraction(target_weight - i - 2 * l, j)
                    if 1 < weight <= Fraction(3, 2):
                        ties.setdefault(weight, set()).add((component, i, j, l))
    return ties


assert tie_terms() == {
    Fraction(6, 5): {(2, 0, 5, 0)},
    Fraction(5, 4): {(1, 0, 4, 1), (2, 1, 4, 0)},
    Fraction(4, 3): {(0, 0, 3, 2), (1, 1, 3, 1), (2, 0, 3, 1), (2, 2, 3, 0)},
    Fraction(7, 5): {(1, 0, 5, 0)},
    Fraction(3, 2): {
        (0, 0, 4, 1), (0, 1, 2, 2), (0, 2, 4, 0), (0, 3, 2, 1),
        (1, 0, 2, 2), (1, 1, 4, 0), (1, 2, 2, 1), (1, 4, 2, 0),
        (2, 0, 4, 0), (2, 1, 2, 1), (2, 3, 2, 0),
    },
}

def first_axial_ties():
    ties = set()
    for component, (i0,l0) in enumerate(((4,2),(5,1),(6,0))):
        for l in range(l0+1,6):
            for i in range(6-l):
                for j in range(6-i-l):
                    w = Fraction(i0-i-j,l-l0)
                    if Fraction(3,2) < w < 2:
                        ties.add((w,component,i,j,l))
    return ties

assert first_axial_ties() == {
    (Fraction(5,3),1,0,0,4), (Fraction(5,3),2,1,0,3), (Fraction(5,3),2,0,1,3),
}

x, P, z = symbols("x P z")
p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14 = symbols("p1:15")
D = symbols("D")
A1 = x**2*P**2 + p1*P**2*z + p2*x*P*z**2 + p3*x**2*z**3 + p4*z**4
B1 = x**4*P + p5*x*P**2 + p6*x**2*P*z + p7*P*z**2 + p8*x**3*z**2 + p9*x*z**3
C1 = x**6 + p10*x**3*P + p11*x**4*z + p12*P**2 + p13*x*P*z + p14*x**2*z**2 + D*z**3
det = Poly(Matrix([A1, B1, C1]).jacobian([x, P, z]).det(), x, P, z)

steps = [
    ((11,0,2),p3,0,-18),
    ((10,1,1),p2,2*p8,-12),
    ((9,0,3),p4,p8**2,-24),
    ((9,2,0),p11,2*p6-p1,-6),
    ((8,1,2),p9,p8*(p6-p1),36),
    ((6,3,0),p13,(p1*p10-14*p1*p5+4*p10*p6+4*p5*p6)/6,-6),
    ((7,2,1),p14,(p1-p6)**2+(p10-4*p5)*p8+2*p7,-12),
    ((5,2,2),D,4*p1*p5*p8-4*p5*p6*p8+2*p6*p7
        -p1*((p1-p6)**2+(p10-4*p5)*p8+2*p7)
        +(p1*p10-14*p1*p5+4*p10*p6+4*p5*p6)*p8/6,-18),
]
substitutions = {}
for term, eliminated, value, factor in steps:
    coefficient = det.coeff_monomial(term).subs(substitutions)
    assert expand(coefficient - factor * (eliminated - value)) == 0
    substitutions[eliminated] = value

equations = sorted({
    expand(c.subs(substitutions, simultaneous=True)) for c in det.coeffs()
} - {0}, key=str)
generators = ",\n".join(
    str(Poly(c).clear_denoms()[1].as_expr()).replace("**","^") for c in equations
)

singular = """
ring R = 0,(p1,p5,p6,p7,p8,p10,p12),dp;
ideal I =
""" + generators + """;
ideal claims = (p1*p10)^6,(p1*p8)^2,(p1*p5)^4,
(p1^2-p1*p6-p5*p8+p7)^3,(p12*(p1*p6+p8*p10-2*p7))^3,
(p6*(p10-2*p5))^6,(p6*p10^2+8*p1*p12-4*p6*p12)^6,
(p8*(p5*p10-2*p12))^3,(p8*(p5^2-p12))^3;
matrix witnesses = lift(I,claims);
if (matrix(I)*witnesses == matrix(claims)){
    print("Verified all nine membership identities over QQ.");
} else {
    ERROR("The coefficient ideal identities failed.");
}
quit;
"""
out = run(["Singular", "-q"], input=singular, text=True, capture_output=True, check=True)
assert out.stdout.strip() == "Verified all nine membership identities over QQ." and not out.stderr.strip(), out.stdout + out.stderr
print("Verified the two weight tie lists, eight eliminations, and nine rational ideal identities.")
