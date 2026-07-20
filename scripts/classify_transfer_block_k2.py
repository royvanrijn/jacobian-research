"""Exact local-algebra classification of the two-transfer block Z_2."""

import sympy as sp

Z = sp.symbols("Z")
a, b, c, d = sp.symbols("a b c d")
u5, u4, u3, u2, u1, u0 = sp.symbols("u5 u4 u3 u2 u1 u0")

V = Z**4+a*Z**3+b*Z**2+c*Z+d
U = Z**6+u5*Z**5+u4*Z**4+u3*Z**3+u2*Z**2+u1*Z+u0
difference = sp.Poly(sp.expand(U**2-V**3), Z)

# The six highest coefficient equations are triangular and eliminate U.
solution = {}
for degree, variable in zip(range(11, 5, -1), (u5, u4, u3, u2, u1, u0)):
    equation = sp.expand(difference.coeff_monomial(Z**degree).subs(solution))
    roots = sp.solve(equation, variable, dict=False)
    assert len(roots) == 1
    solution[variable] = sp.factor(roots[0])

remaining = [
    sp.factor(difference.coeff_monomial(Z**degree).subs(solution))
    for degree in range(5, -1, -1)
]

# Adapt V to the reduced parameter S=Z^2+pZ+q:
# V=S^2+X*Z+Y.
p, q, X, Y = sp.symbols("p q X Y")
adapted = {
    a: 2*p,
    b: p**2+2*q,
    c: 2*p*q+X,
    d: q**2+Y,
}
remaining = [sp.factor(equation.subs(adapted)) for equation in remaining]

candidate = [X**3, 2*X*Y-p*X**2, Y**2-q*X**2]
computed_groebner = sp.groebner(remaining, Y, X, q, p, order="lex", domain=sp.QQ)
candidate_groebner = sp.groebner(candidate, Y, X, q, p, order="lex", domain=sp.QQ)

# Equality of ideals, checked in both directions.
assert all(candidate_groebner.reduce(equation)[1] == 0
           for equation in remaining)
assert all(computed_groebner.reduce(equation)[1] == 0
           for equation in candidate)

# Monic leading terms Y^2, X*Y, X^3 give a free k[p,q]-basis.
basis = (sp.Integer(1), X, Y, X**2)
assert all(candidate_groebner.reduce(element)[1] == element for element in basis)

# The coincident-root fiber p=q=0 has Hilbert function (1,2,1), nilpotency
# index three, and two-dimensional socle span(X^2,Y); it is not the
# Gorenstein tensor product k[e1,e2]/(e1^2,e2^2).
special_relations = [relation.subs({p: 0, q: 0}) for relation in candidate]
special = sp.groebner(special_relations, Y, X, order="lex", domain=sp.QQ)
assert all(special.reduce(element)[1] == element for element in basis)
assert special.reduce(X**2)[1] == X**2
assert special.reduce(X**3)[1] == 0
assert special.reduce(X*Y)[1] == 0
assert special.reduce(Y**2)[1] == 0

print("PASS: the six high coefficients eliminate U triangularly")
print("PASS: Z_2 has ideal (X^3, 2XY-pX^2, Y^2-qX^2)")
print("PASS: Z_2 is finite flat of rank 4 over k[p,q]")
print("PASS: the coincident-root fiber is k[X,Y]/(X^3,XY,Y^2)")
