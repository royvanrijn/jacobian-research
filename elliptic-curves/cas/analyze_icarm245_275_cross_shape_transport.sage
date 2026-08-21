#!/usr/bin/env sage
"""Classify the smallest cross-shape transports between ICARM 245 and 275.

Both curves are specializations of Fermigier's two-parameter six-root
construction.  Exact finite-reduction bases decompose each rank-20 subgroup
as twelve generic directions plus eight exceptional directions.  This script
joins every exceptional abscissa on curve 245 to every exceptional abscissa
on curve 275 by the affine line in the full ``(u,v,t,x)`` space.

For each of the 64 lines it evaluates Fermigier's quartic exactly, removes
the complete rational square factor, and records the degree (and hence genus)
of the residual square condition.  This is a finite classification of this
specific affine ansatz; it is not an exclusion of nonlinear transports.
"""

from collections import Counter

from sage.all import PolynomialRing, QQ


parameter_ring = PolynomialRing(QQ, "k")
k = parameter_ring.gen()
function_field = parameter_ring.fraction_field()
x_ring = PolynomialRing(function_field, "X")
X = x_ring.gen()


def fermigier_roots(u, v):
    """Return Fermigier's six labelled roots over an arbitrary QQ-algebra."""

    alpha1 = (
        -v + v**2 - v**3 + v**4 + 2*u*v + v**2*u - 2*v**3*u
        + v**4*u + u**2 + u**2*v - 2*v**2*u**2 - 2*v**3*u**2
        + u**3*v**2 - u**4
    )
    alpha2 = (
        v**3*u**2 - 2*v**2*u**2 + v**2*u - 2*u**3*v**2
        - 2*u**3*v + u**2*v + u**4*v + 2*u*v - v**4 - u
        + v**2 - u**3 + u**2 + u**4
    )
    alpha3 = (
        -v**4*u + 2*v**3*u + v**3*u**2 + v**2*u**2 + v**2*u
        - u**3*v**2 - 2*u**3*v - 2*u**2*v + u**4*v - v + v**3
        - 2*u**3 + u**2 + u**4
    )
    alpha4 = (
        v - v**2 + v**3 - v**4 + u - 2*u*v + v**2*u + 2*v**3*u
        - 2*u**2 - 2*u**2*v + v**2*u**2 + v**3*u**2 + u**3
        - u**4*v
    )
    alpha5 = (
        v**4*u - 2*v**3*u - v**3*u**2 + v**2*u**2 - 2*v**2*u
        + u**3*v**2 + 2*u**3*v + u**2*v + v**4 - u - 2*v**3
        + v**2 + u**3 - u**4*v
    )
    alpha6 = (
        v - 2*v**2 + v**3 + u - 2*u*v - 2*v**2*u - v**4*u
        - u**2 + u**2*v + v**2*u**2 + u**3 + 2*u**3*v
        + u**3*v**2 - u**4
    )
    return alpha1, alpha2, alpha3, alpha4, alpha5, alpha6


def mestre_quartic(roots, parameter):
    """Construct ``(g^2-q(X-t)q(X+t))/t^2`` exactly."""

    product = x_ring(1)
    for root in roots:
        product *= X-(root+parameter)
        product *= X-(root-parameter)
    approximant = [function_field(0)]*7
    approximant[6] = function_field(1)
    for index in range(5, -1, -1):
        polynomial = sum(
            approximant[degree]*X**degree for degree in range(7)
        )
        square = polynomial**2
        approximant[index] = (
            product[6+index]-square[6+index]
        )/2
    polynomial = sum(
        approximant[degree]*X**degree for degree in range(7)
    )
    remainder = polynomial**2-product
    assert remainder.degree() <= 4
    return remainder/parameter**2


def rational(value):
    return QQ(value)


# Canonical quartic abscissas selected by exact mod-3 rank-20 certificates.
curve245_canonical = tuple(map(rational, (
    "-1069/530", "1107/2", "-1217/10", "-2901/10",
    "3085/46", "-5773/2", "9517/170", "9933/10",
)))
curve275_canonical = tuple(map(rational, (
    "-94898/11", "79501/73", "92446/307", "245933/341",
    "136343/349", "-250232/583", "454737/605", "379655/803",
)))

# Native Fermigier coordinates.  The canonical normalizations are
# x245=16*(x+375/16) and x275=(16/3)*(x+1671/16).
curve245_native = tuple(value/16-QQ(375)/16 for value in curve245_canonical)
curve275_native = tuple(QQ(3)/16*value-QQ(1671)/16 for value in curve275_canonical)

u = function_field((3-9*k)/2)
v = function_field((4-5*k)/2)
t0 = QQ(5801)/160
t1 = QQ(10239)/176
t = function_field(t0+k*(t1-t0))
quartic = mestre_quartic(
    tuple(function_field(root) for root in fermigier_roots(u, v)), t
)


def odd_squareclass(polynomial):
    polynomial = parameter_ring(polynomial)
    answer = parameter_ring(1)
    for factor, exponent in polynomial.squarefree_decomposition():
        if exponent % 2:
            answer *= factor
    return answer.monic()


records = []
for left_index, left in enumerate(curve245_native, 1):
    for right_index, right in enumerate(curve275_native, 1):
        x_value = function_field(left+k*(right-left))
        value = function_field(quartic(x_value))
        assert value(0).is_square() and value(1).is_square()
        numerator_kernel = odd_squareclass(value.numerator())
        denominator_kernel = odd_squareclass(value.denominator())
        kernel = (numerator_kernel*denominator_kernel).monic()
        degree = kernel.degree()
        genus = max(0, (degree-1)//2)
        records.append((degree, genus, left_index, right_index, kernel))

degree_counts = Counter(record[0] for record in records)
minimum_degree = min(degree_counts)
minimum_records = [record for record in records if record[0] == minimum_degree]
minimum_factor_patterns = Counter(
    tuple((factor.degree(), exponent) for factor, exponent in kernel.factor())
    for _, _, _, _, kernel in minimum_records
)

print(
    "ICARM245275XSHAPE|pairs=64|path=affine(u,v,t,x)|"
    f"kernel_degree_counts={tuple(sorted(degree_counts.items()))}|"
    f"minimum_degree={minimum_degree}|minimum_genus={(minimum_degree-1)//2}",
    flush=True,
)
print(
    "ICARM245275XSHAPE|"
    f"minimum_factor_patterns={tuple(sorted(minimum_factor_patterns.items()))}",
    flush=True,
)
assert minimum_degree > 4
print(
    "ICARM245275XSHAPE|status=PASS_NO_GENUS_AT_MOST_ONE_AFFINE_TRANSPORT",
    flush=True,
)
