"""Replay the exact coefficient and collision identities in the conductor-ten proof."""

from sympy import Poly, Rational, diff, expand, factor, symbols


z, a, b, c, k = symbols("z a b c k")
da, db, dc, dk = symbols("da db dc dk")
S, P, A0, B0 = symbols("S P A0 B0")
U = z**4 + a*z**2 + b*z + c
D = z**2 + a/2
V = 3*U*D/2 - D**3/2 + 3*b**2/8 + k*z
delta = a**2 - 4*c
A = b*delta - 4*k*a/3
B = delta**2 - 8*a*b**2 + 32*k*b/3
G = b**3/4 - 3*A*a/16 - 3*k*a**2/8 + k*delta/4


def derivative(f):
    return sum(diff(f, x)*dx for x, dx in zip((a, b, c, k), (da, db, dc, dk)))


jac = Poly(expand(diff(U, z)*derivative(V) - derivative(U)*diff(V, z)), z)
assert Poly(expand(V**2-U**3), z).degree() == 7
assert Poly(expand(V**2-U**3), z).LC() == 2*k
assert jac.degree() == 4 and jac.nth(4) == 4*dk and jac.nth(3) == 0
assert expand((jac.nth(2) - 3*derivative(A)/4).subs(dk, 0)) == 0
assert expand((jac.nth(1) + 3*derivative(B)/32).subs(dk, 0)) == 0
assert expand((jac.nth(0)-derivative(G)-3*a*derivative(A)/16).subs(dk, 0)) == 0

# The divided difference of z^n is the complete symmetric polynomial
# of degree n-1 in its two arguments.
differences = [0, 1]
for _ in range(2, 7):
    differences.append(expand(S*differences[-1] - P*differences[-2]))


def divided_difference(f):
    return sum(Poly(f, z).nth(i)*differences[i] for i in range(7))


assert expand(divided_difference(U) - (S**3-2*S*P+a*S+b)) == 0
quintic = 2*S**5 + 4*a*S**3 - 2*b*S**2 + 3*delta*S - 8*k
product = (S**3+a*S+b)/(2*S)
assert factor(divided_difference(V).subs(P, product) + quintic/8) == 0
gap = -S**2 - 2*a - 2*b/S
gap_reduced = -3*b/S + 3*delta/(2*S**2) - 4*k/S**3
assert factor(gap - gap_reduced + quintic/(2*S**3)) == 0

# Use b=w^(-6) so all five candidate scales are Laurent monomials.
w, r, d = symbols("w r d")
quintic_d = 2*S**5+4*a*S**3-2*b*S**2+3*d*S-8*k
small = {a:4*k/(3*b), d:A0/b+16*k**2/(9*b**2)}
large = {a:9*b**4/(2*k**2), d:6*b**3/k}


def leading_weight(expr):
    terms = expand(expr).as_ordered_terms()
    order = min(term.as_powers_dict().get(w, 0) for term in terms)
    return order, expand(sum(term/w**order for term in terms
                             if term.as_powers_dict().get(w, 0) == order))


for branch, power, order, expected in (
    (small, -2, -10, 2*r**2*(r**3-1)),
    (small, 3, 0, -2*r**2-8*k),
    (large, -12, -60, 2*r**3*(r**2+9/k**2)),
    (large, 3, -15, 18*r*(r**2+k)/k**2),
    (large, 18, 0, 18*r/k-8*k),
):
    actual_order, initial = leading_weight(
        quintic_d.subs(branch).subs({b:w**-6,S:r*w**power}))
    assert actual_order == order and expand(initial-expected) == 0

assert [Rational(10-2*i, 2*i+3) for i in range(1,5)] == [
    Rational(8,5), Rational(6,7), Rational(4,9), Rational(2,11)]
print("Verified degree-ten parameter, collision, and root-scale identities.")
