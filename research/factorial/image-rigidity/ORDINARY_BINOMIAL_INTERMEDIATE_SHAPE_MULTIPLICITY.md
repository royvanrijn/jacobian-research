# Intermediate-shape root multiplicity for ordinary binomial moments

For integers \(m\geq0\), \(d\geq1\) and \(\gamma>0\), use the
monic reciprocal moment polynomial
\[
  Q_{m,\gamma,d}(x)=\sum_{k=0}^m \binom{m}{k}(\gamma)_{dk}x^{m-k},
  \qquad (\gamma)_j=\prod_{i=0}^{j-1}(\gamma+i).
\]
For ordinary factorial moments and support \(X^a(A+B X^d)\), the two
actual reciprocal polynomials at terminal \(N\) are
\(Q_{N-1,\gamma,d}\) and \(Q_{N,\gamma+a,d}\), with
\(\gamma=a(N-1)+1\). Fix this \(\gamma\), set \(P=Q_{N,\gamma,d}\),
and write \(E=x\,d/dx\).

## The derivative and shape-shift identities

For \(N\geq2\) one has
\[
  Q_{N-1,\gamma,d}=\frac{P'}{N},\qquad
  (\gamma)_a Q_{N,\gamma+a,d}
   =\prod_{j=0}^{a-1}(\gamma+j+dN-dE)\,P .                 \tag{1}
\]
The empty product is the identity. For the first equation, compare
coefficients using \((N-k)\binom Nk=N\binom{N-1}k\).
For the second equation each factor multiplies the coefficient of
\(x^{N-k}\) by \(\gamma+j+dk\). Their product gives
\((\gamma)_{dk}(\gamma+dk)_a=(\gamma)_a(\gamma+a)_{dk}\).
The operators commute since each is a polynomial in \(E\).

Thus the actual consecutive coprimality assertion is exactly the
coprimality of \(P'\) and the shifted operator expression on the
right side of (1). The rescaling factors are nonzero rational
numbers. Over any characteristic-zero field these rational
polynomial assertions are equivalent to coprimality over
\(\mathbb C\), by Bézout and scalar extension.

## Uniform multiplicity bound

The [all-order ordinary binomial theorem](ORDINARY_BINOMIAL_RESIDUAL_REGION.md)
implies
\[
  \operatorname{mult}_x Q_{N,a(N-1)+1,d}\leq\max(1,a)
\]
at every root, for all \(a\geq0\), \(d\geq1\) and \(N\geq2\).
Indeed, if a root has multiplicity at least \(\max(2,a+1)\), then
\(P'\) and every derivative of \(P\) of order at most \(a\) vanish
there. The shifted operator in (1) has differential order at most
\(a\), so this would give a common root of the two actual moment
polynomials. Since \(P\) has degree \(N\), one can also take the
minimum of this bound and \(N\).

At \(a=0\), equation (1) gives the exact equivalence between
ordinary constant-plus-power consecutive coprimality and
squarefreeness of \(Q_{N,1,d}\). At \(a=1\), at any root of
\(P'\), equation (1) reduces to
\[
  \gamma Q_{N,\gamma+1,d}=(\gamma+dN)P .
\]
Thus at \(\gamma=N\) the ordinary lower-one consecutive assertion
is equivalent to squarefreeness of \(Q_{N,N,d}\). Both implications
hold separately for every \(N,d\). In terms of the reciprocal
\(T_{N,N,d}\), the same equivalence follows because \(Q_{N,N,d}\)
has nonzero constant and leading coefficients. The older
[lower-one theorem](ORDINARY_LOWER_ONE_STRONG_FACTORIAL.md) therefore
gives this squarefreeness at every gap, including gap six.

At \(a=2\), for example, write \(c=\gamma+dN\). At every critical
point of \(P\) the shifted operator in (1) reduces to
\[
  c(c+1)P+d^2x^2P'' .
\]
The all-order ordinary theorem excludes its vanishing at any of
those critical points. The uniform multiplicity statement concerns
the fixed-intermediate-shape polynomial \(Q_{N,a(N-1)+1,d}\).
Squarefreeness of the actual upper polynomial \(Q_{N,aN+1,d}\)
for positive \(a\) remains a separate assertion.

## Ordinary lower-one Hurwitz determinants

The [uniform discriminant identity](gap_three_hurwitz_structure.tex)
at `thm:discriminant` and the sign argument at
`cor:top-nonnegative` give
\(\mathcal D_{N,d}(N)>0\) if and only if \(T_{N,N,d}\) is squarefree.
The adjoint-polynomial recursion in that source has integer
coefficients at \(\gamma=N\) starting from one. Its cyclic matrix
minors, the coefficients of \(W^{[d]}_{N,N}\), and the resulting
Hurwitz determinant are integers. Consequently, for all
\(d\geq1\) and \(N\geq2\),
\[
    \mathcal D_{N,d}(N)\geq1 .
\]
At gap three, the same argument and the two positive tail relations
give \(\Delta_{2N-2},\Delta_{2N-1},\Delta_{2N}\geq1\) at the
ordinary lower-one shape. Conversely, the determinant bound for
all \(d,N\) implies the ordinary lower-one consecutive theorem
by the discriminant identity and (1). The positivity for all
positive Gamma shapes and the sign pattern of other Hurwitz minors
have their own quantifiers.

