# Generated Sage input for the common-arithmetic-fibers explicit example.
Q.<S> = PolynomialRing(QQ)
R.<x,y,z> = PolynomialRing(QQ, order='degrevlex')
factor_coefficients_ascending = [['-19', '0', '0', '1'], ['1', '1', '1']]
factors = [
    sum(QQ(c) * S^i for i, c in enumerate(coefficients))
    for coefficients in factor_coefficients_ascending
]
P5 = prod(factors)
G5 = P5 - P5(0)
t = 1 + x*y
q = t^2*z - 19*y^2*(3*t + 1)
F = vector(R, [
    q*t,
    -5*q^5*t^2*x^3 - 4*q^4*t^2*x^2 + 38*q*t - 3*q*x + 19*y,
    3*q^5*x^5 + 2*q^4*x^4 + x^3*z + 19*x*(5 - 3*t),
])
normalized_target = vector(QQ, ['1', '0', '-2'])
integral_scaling = diagonal_matrix(QQ, [1, 19, 19])
integral_target = integral_scaling * normalized_target
inverse_polynomial = G5 - G5[1] * (
    normalized_target[1]*S^2 + normalized_target[2]
) / 2
assert P5 == S^5 + S^4 + S^3 - 19*S^2 - 19*S - 19
assert inverse_polynomial == P5
assert integral_target == vector(QQ, ['1', '0', '-38'])
assert det(jacobian(F, (x, y, z))) == QQ("-722")
print("PASS: generated Sage paper example")
