# Differential order restrictions for polynomial Weyl pairs

**Theorem PolynomialGcdTwoNondividingOrdersExclusion.** Over a
characteristic-zero field, polynomial differential operators of orders
\(2a,2b\), with \(a,b\ge2\) coprime, cannot have a nonzero scalar
commutator.

**Corollary PolynomialWeylReducedOrderBounds.** Suppose
\(P,Q\in A_1(k)=k[x][D;\partial_x]\) have nonzero scalar commutator
and generate a proper subalgebra. Each differential order is composite
and at least six. By successive replacements of one member by itself
minus a polynomial in the other, the pair can be brought to a
presentation with orders at least six, neither dividing the other,
and with gcd at least three. If one member has order six, it can be
kept fixed, and the reduced partner has order at least nine and
congruent to three modulo six.

## Proof

For the theorem the
[higher-root argument](POLYNOMIAL_GCD_TWO_HIGHER_ROOT_EXCLUSION.md)
forces the common factor of the leading coefficients to have degree
at most one. The
[infinity argument](POLYNOMIAL_GCD_TWO_LOW_DEGREE_LEADING_FACTOR.md)
excludes both constant and linear common factors.

To prove the corollary, take \(P\) of order \(n>0\), \(Q\) of order
\(m>0\), and let \(p,q\in k[x]\) be their respective leading
coefficients. If \(m=jn\) for a positive integer \(j\), the
order-\((n+m-1)\) commutator coefficient gives
\[
                    n p q'-m q p'=0.
\]
Thus \(q/p^j\in k^\times\). Replacing \(Q\) by \(Q-(q/p^j)P^j\)
strictly decreases its order and preserves both the commutator
and the generated subalgebra. This replacement is induced by a
triangular automorphism of the presenting Weyl algebra after scaling
its generators to the relation with commutator one.

By [QLO1](LOW_ORDER_DC1_AUTOMORPHISMS.md), a pair generating a proper
subalgebra has both orders at least four. If either order is prime,
the [original-field prime-order normal form](ALL_COPRIME_ALGEBRAIC_STRING_OBSTRUCTION.md)
subtracts a polynomial in that member from its mate to produce a
mate of order one. It then generates by the same low-order result.
This rules out prime orders in every polynomial presentation.

If \(P\) has order four, successively subtract powers of this fixed
member from \(Q\) while its order is divisible by four. The residual
partner still has order at least four. An odd order is excluded by
the [coprime scalar theorem](ALL_COPRIME_ALGEBRAIC_STRING_OBSTRUCTION.md),
and an even order congruent to two modulo four is excluded by the
theorem above. This rules out order four in every polynomial
presentation. The minimum composite order still possible is six.

For a general pair, whenever either order divides the other make the
same subtraction, interchanging their roles if needed. The sum of the
orders strictly decreases at each step and stays at least twelve, so
the process stops with neither order dividing the other. For those
orders, gcd one is excluded by the coprime theorem and gcd two by
the displayed theorem. If one member has order six, subtract only
from its partner, since order less than six is already impossible.
The final partner order is not divisible by six and has a nontrivial
common divisor with six other than two. It is therefore congruent to
three modulo six and at least nine.

All changes here take place in the one-variable polynomial Weyl
algebra. In particular, the order bounds concern globally polynomial
presentations. No normalization into a polynomial presentation is
asserted for a pair over a larger differential coefficient field.
