# Consecutive binomial moment schemes and Gaussian transfer

The [all-order ordinary binomial theorem](ORDINARY_BINOMIAL_RESIDUAL_REGION.md#ordinary-binomial-strong-factorial)
applies to every integer support \(a\geq0,d\geq1\) and terminal \(M\geq2\).
The paper gives the same statement as its [main theorem](../../../results/strong-factorial-binomials/main.tex#L51).
The [fixed-pair coefficient and radial argument](SUPPORT_UNIFORM_COEFFICIENT_NILPOTENCE_AND_RADIAL_GAUSSIAN.md)
started from a support-uniform eventual theorem. The stronger input now gives its
conclusions at every terminal. All rings below are commutative unital
\(\mathbb Q\)-algebras.

## The universal two-moment quotient

For \(F(Z)=A Z^a+B Z^{a+d}\), let \(\mu_m(A,B)=\Lambda(F^m)\), where
\(\Lambda(Z^k)=k!\). For fixed \(a,d,M\), write
\[
 R_M=\mathbb Q[A,B]/(\mu_{M-1},\mu_M),\qquad J=(A,B)R_M.
\]
The ordinary theorem says the two forms have no common projective zero over
\(\mathbb C\). The fixed-pair argument in §1 of the cited note then gives a
homogeneous regular sequence of degrees \(M-1,M\). Its Hilbert series, length
and ordinary nilpotence index are respectively
\[
 (1+\cdots+t^{M-2})(1+\cdots+t^{M-1}),\qquad M(M-1),\qquad 2M-2.
\]
The top graded piece \(J^{2M-3}\) is a line over \(\mathbb Q\). This proves
that in every coefficient ring the two equations force \((A,B)^{2M-2}=0\).
The coefficients of \(F^r\) are unit rational multiples of the generators
\(A^{r-j}B^j\) of \((A,B)^r\) at distinct exponents. It follows that
\((A,B)^r=0\) if and only if \(F^r=0\). In a reduced coefficient ring the two
equations force \(A=B=0\).

There is a reverse implication for the universal assertion. Any nonzero
common coefficient vector over \(\mathbb C\) would give a specialization
in which \((A,B)\) is the unit ideal and the two moments vanish. Applying
the universal ideal bound excludes such vectors. The argument applies
at each pair \(M-1,M\), so quantifying it over all supports and all terminals
recovers the full ordinary theorem. The same reverse test works directly
with polynomial nilpotence or with reduced coefficient rings.

## Two literal Gaussian laws

Let \(X,Y,Z\) be independent standard normal variables and \(G\) their
normalized Wick functional. For \(T=(X^2+Y^2)/2\) one has
\(G(T^k)=k!\). The exact factorial-to-Gaussian substitution in §§3–4 of
the cited fixed-pair note is injective over every coefficient ring and
preserves the coefficient ideal. It gives both directions of the
universal conclusion for \(H=T^a(A+B T^d)\) in dimension two.

A second realization uses independent radial variables
\(S=X^2/2\), \(T'=(Y^2+Z^2)/2\) of shapes \(1/2,1\), and
\(U=4ST'=X^2(Y^2+Z^2)\). For each integer \(k\geq0\) one has
\[
 G(U^k)=4^k(1/2)_k(1)_k=(2k)!.
\]
The first equality follows by independence and the normalized Gaussian
moments. The second follows by splitting \((2k)!\) into its even and odd
factors. Thus if \(H=U^a(A+B U^d)\), its two Gaussian moment polynomials
are exactly those for the ordinary support \((2a,2d)\) in the
lower-exponent and gap coordinates. Both substitutions have injective
powers over every coefficient ring. Distinct powers have distinct
homogeneous degrees and each power of \(U\) has coefficient one on
\(X^{2k}Y^{2k}\). In particular the coefficient ideals and literal
nilpotence exponents are preserved.

In each Gaussian family over a characteristic-zero field, a nonzero phase
has a nonzero unmarked moment in every consecutive positive pair.
Over a \(\mathbb Q\)-algebra the two equations force \(H^{2M-2}=0\), hence
every original polynomial marker \(Q\) satisfies \(G(QH^r)=0\) for
\(r\geq2M-2\). This marked cutoff is sharp for each support and pair.
In the corresponding universal quotient, take \(r=2M-3\) and write
the nonzero polynomial \(H^r=e h\), where \(e\) spans the top graded
line over \(\mathbb Q\) and \(h\) is a nonzero polynomial with rational
coefficients in the original variables. Use \(h\) itself as the original
marker. Then \(G(hH^r)=eG(h^2)\ne0\), since the real standard Gaussian
has a strictly positive density. Both prescribed unmarked moments still
vanish in that quotient.

## Independent common factors and ordinary restrictions

If two multivariate ordinary factorial monomials differ in exactly one
coordinate, write their binomial as \(V F(Z)\) with \(V\) the common
monomial in the other coordinates. At each order \(m\) its moment is
\(\Lambda(V^m)\mu_m(A,B)\) by the product formula for the factorial
functional. The first factor is a product of factorials and a unit in
every coefficient ring considered here. The one-active-coordinate binomial
family therefore has the same consecutive nonvanishing theorem and exact
universal coefficient ideal bound. Taking \(V=1\) recovers the univariate
statement. The argument applies to the independently distributed variables
in this construction, including any independent Gaussian radial factors
whose pure moments are nonzero rationals.

The ordinary theorem also directly covers its former restrictions on
\(a,d,M\), including lower-exponent slices, finite boxes, total-degree
bounds, square-root-gap and integer-cusp supports. For the normalized
polynomials at a common nonzero affine or multiplicative change of the
root variable, an alleged common zero maps to a common zero of the same
ordinary consecutive pair. The old sufficient analytic, local-arithmetic
and irreducibility questions can ask for stronger assertions than this
exclusion. For instance individual squarefreeness for general positive
lower exponent and consecutive coprimality at arbitrary positive Gamma
shape remain separate.

