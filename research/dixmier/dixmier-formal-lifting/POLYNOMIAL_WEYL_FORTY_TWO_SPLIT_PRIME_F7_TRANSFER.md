# Split-prime transfer of the degree-forty-two Weyl residual to planar F7

**Theorem PolynomialWeylFortyTwoSplitPrimeF7Transfer.**
An actual proper polynomial Weyl pair over \(\mathbb C\) with a
member of Bernstein degree \(42\) would give a characteristic-zero
commutative Keller pair satisfying every hypothesis of Standing
Data 1.1 in the \(F_7\) case of
[Niedbala Giraudin, Hexagon 2611.00896v1](https://hexagonmath.org/2611.00896v1).

**Corollary PolynomialWeylFortyTwoExclusionUnderPlanarF7.**
The normalized planar exclusion in Theorem 1.2 of that
preprint excludes the degree-forty-two polynomial Weyl residual.
The implication uses the pinned \(F_7\) theorem, without an
exhaustiveness assumption on any planar degree census.

## Algebraic coefficients and split primes

Take \(P,Q\) in the centered Fourier presentation of the
[degree-forty-two Weyl residual](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md),
the [finite-pole theorem](POLYNOMIAL_WEYL_FORTY_TWO_FINITE_POLE.md),
and the
[regular finite-pole theorem](POLYNOMIAL_WEYL_FORTY_TWO_REGULAR_FINITE_POLE.md).
Use \([Y,X]=1\) and \([Q,P]=1\).
The support rectangles have corners \((12,30),(42,105)\), with
those corners the only monomials of maximal \(X\)-degree or
maximal \(Y\)-degree. The lowest Euler grades are \(-18,-63\),
where the Euler grade of \(X^iY^j\) is \(i-j\).
For \(C=XY\) their components are \(X^{-18}f(C),X^{-63}g(C)\)
with \(\deg f=30,\deg g=105\).
At weight \((3,-1)\), their leading forms are
\[
     \alpha \bigl(X h(XY^3)\bigr)^2,\qquad
     \beta \bigl(X h(XY^3)\bigr)^7,\qquad
     h(z)=(z-u)^4(z-v),\qquad \alpha\beta uv(u-v)\ne0.       \tag{1}
\]
Here \(\alpha,\beta\) also give the single-corner coefficients.

We can take these data over a number field. Keep the two support
rectangles, the pure maximal-coordinate faces, the two lower
Euler bounds, the equations (1), and the equation \([Q,P]=1\).
These are finitely many polynomial equations over \(\mathbb Q\)
in the coefficients and \(\alpha,\beta,u,v\), with nonzero
conditions imposed by inverse variables. Their complex solution
gives a solution over \(\overline{\mathbb Q}\) by the
Nullstellensatz. The argument requires only these data from the
algebraic solution. In particular it imposes the support bounds
before specializing coefficients.

Choose a finite Galois extension of \(\mathbb Q\) containing
these algebraic coefficients and all roots of \(f,g\).
Enlarge the finite set of excluded primes to accommodate
denominators, the required nonzero coefficients, and primes
\(p\le147\). There are infinitely many remaining rational
primes splitting completely in this Galois extension.
One elementary proof uses a primitive integral element.
Infinitely many primes divide the nonzero values of its
minimal polynomial. Indeed, if there were a finite list
of such primes, choose a nonzero value \(A\) and restrict
the input to its original argument plus multiples of
\(A\) times the product of the listed primes.
The resulting polynomial values divided by \(A\) are
integers congruent to one modulo each listed prime,
and a value of absolute size greater than one gives a
contradiction. Away from the finite exceptional primes,
a root of the minimal polynomial modulo \(p\) gives a
degree-one prime. Normality then gives complete splitting.
Thus we may reduce the Weyl data with both \(f\) and \(g\)
split into linear factors over \(\mathbb F_p\).

## What the center remembers

For each such prime, write \(x=X^p,\ y=Y^p\) for the two
central variables of the reduced Weyl algebra. Denote
\(P^p,Q^p\), expressed as polynomials in \(x,y\), by
\(U_p,V_p\), respectively. They are central by Theorem
1.7(iii)(a) of
[Lauritzen and Thomsen, arXiv:1508.02499v1](https://arxiv.org/pdf/1508.02499v1).
Their ordinary Jacobian determinant is a nonzero scalar
by Proposition 1.12 and the standard center Poisson
bracket of Proposition 1.11 of the same paper.
These results apply to the reduction of the actual Weyl
endomorphism over its finitely generated coefficient
ring of characteristic zero.

At \((1,0),(0,1),(3,-1)\) the contraction loss is strictly
positive. Consequently the leading symbols of \(P^p,Q^p\)
are the \(p\)-th powers of the reduced leading symbols of
\(P,Q\) at these weights. Rewriting in \(x,y\) gives the
same two support rectangles, the same unique maximal
\(x\)-degree and maximal \(y\)-degree corners, and
\[
  \operatorname{in}_{3,-1}(U_p)=\overline\alpha^{\,p}
        \bigl(x\,\overline h^{(p)}(xy^3)\bigr)^2,\qquad
  \operatorname{in}_{3,-1}(V_p)=\overline\beta^{\,p}
        \bigl(x\,\overline h^{(p)}(xy^3)\bigr)^7.            \tag{2}
\]
Here the superscript on \(\overline h^{(p)}\) applies
Frobenius to coefficients, so the multiplicities and
distinctness of the two nonzero roots are unchanged.

The Euler grading has zero contraction loss, but it is
an exact grading. The lowest components of the powers
are the exact products of the lowest components.
For \(a=18\) or \(63\) and the corresponding split
coefficient polynomial \(F=f\) or \(g\), one has in the
localized Weyl algebra
\[
 (X^{-a}F(C))^p=X^{-ap}\prod_{j=0}^{p-1}F(C-aj).
\]
Because \(p\nmid a\), the orbit of each reduced root under
the \(a\)-shift is the prime field. Thus each root
contributes \(C^p-C\) to the product. The falling-factorial
identity \(X^pY^p=(C)_p=C^p-C\) now gives
\[
     \operatorname{in}_{-1,1}(U_p)
          =\overline\alpha^{\,p} x^{12}y^{30},\qquad
     \operatorname{in}_{-1,1}(V_p)
          =\overline\beta^{\,p} x^{42}y^{105}.              \tag{3}
\]
In particular there is a single maximizing point at the
weight \((-1,1)\) for each center polynomial. Splitting
the indicial roots in the prime field is the step that
renders these two initial forms monomial.

## Returning the finite support system to characteristic zero

Now consider independent commutative polynomials \(U,V\).
Impose the two rectangles and pure maximal-coordinate
faces above, the two initial forms in the shape (2) with
distinct nonzero roots and nonzero scales, the two pure
initial forms in (3), and a nonzero constant Jacobian.
These conditions are a finite polynomial system over
\(\mathbb Z\) in the coefficients and finitely many
parameters, with the required nonzero products inverted
by introducing inverse variables. In particular the pure
forms in (3) impose the bounds \(j-i\le18,63\) and set
the other coefficients at equality to zero.
The center pair supplies a solution for infinitely many
prime characteristics. If the system had no solution
over \(\overline{\mathbb Q}\), its coefficient equations
would generate the unit ideal over \(\overline{\mathbb Q}\),
and hence over \(\mathbb Q\) by faithful flatness.
Clearing denominators in a certificate would rule out
solutions in every characteristic outside a finite list.
Thus there is a characteristic-zero pair satisfying the
whole system.

## Exact matching to Standing Data 1.1

We check the normalization using
[Guccione, Guccione and Valqui, arXiv:1401.1784v3](https://arxiv.org/pdf/1401.1784v3).
The characteristic-zero pair \(U,V\) is polynomial and
has nonzero scalar Jacobian. Its Bernstein degrees
\(42,147\) and maximal \(x\)-degrees \(12,42\) have the
ratios required by Definition 4.3 for a \((2,7)\)-pair.
The start and end of its \((1,0)\) leading part for \(U\)
are both \((12,30)\), of negative diagonal difference,
so the pair is standard by that same definition.

The \((3,-1)\) edge of \(U\) has normalized endpoints
\((6,15),(1,0)\). The first has larger \(y\)-exponent
and negative diagonal difference. Thus it gives a
regular corner as in Definition 5.5. The leading powers
of the common form commute, and the lower endpoint has
positive diagonal difference and the edge has more
than one root. This is type II.b in the classification
following Remark 5.8, so the starting triple in Definition
6.2 is \(((6,15),(1,0),(3,-1))\), by Proposition 6.1
and Remark 6.3.

The regularity condition for the *left* successor is
also literal. Both rectangles have a single common
normalized corner at the positive coordinate weights.
At every weight \((-s,1)\) with \(0\le s\le1\), this is
still the unique maximizer by (3) and the pure maximal
\(y\)-degree face: its weight is the combination
\((1-s)j+s(j-i)\).
Hence there is no edge from \((1,0)\) through
\((-1,1)\) in the angular order used for the successor.
Since the Jacobian at the origin is nonzero, \(U\)
has a nonzero linear coefficient. The associated
monomial wins against \((12,30)\) when \(s\) is
sufficiently large, so there is an edge at finite
\(s>1\). This proves the strict inequality
\[
                    (-1,1)<\operatorname{Succ}_{U}(1,0)
                              <(-1,0).
\]
Finally the root pattern in (2) is the required \(4+1\)
pattern of \(F_7\). Thus all Standing Data 1.1 hold for
\((U,V)\). Assuming the \(F_7\) exclusion, the pinned Theorem 1.2
excludes the pair.

The finite transfer above checks the hypotheses for this
particular shape and appeals to Theorem 1.2 only at the
last step. The proof in that source uses its two Newton-side
reductions and the exact unit-ideal certificate of
Proposition 5.2. The
[repository replay receipt](../../jacobian/planar-degree-frontier/SHAPE_EXCLUSION_REPLAY.json)
records the internal replay of that certificate.
