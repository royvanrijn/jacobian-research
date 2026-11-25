# Root and coefficient consequences of eventual two-monomial separation

Fix distinct exponent vectors \(\alpha,\beta\) with a private variable
on each side as in the
[eventual separation theorem](MULTIVARIATE_TWO_MONOMIAL_EVENTUAL_BORDER.md).
Write \(F=AX^\alpha+BX^\beta\), \(M_j(A,B)=\Lambda(F^j)\), and let
\(m_0\) be a cutoff for the theorem. Thus \(M_m(1,t)\) and
\(M_{m+1}(1,t)\) are coprime for each \(m\ge m_0\).

## Eventual squarefreeness from a varying shared row

Suppose at least one shared coordinate has \(\alpha_i\ne\beta_i\).
Use the positive normalization \(P_m(z)=L_m(z)+z^m U_m(1/z)\)
and the fixed annulus from the proof of the separation theorem.
The functions \(L_m\) and \(U_m\) are bounded away from zero and
their first derivatives are uniformly bounded on the fixed disk
containing \(z\) and \(1/z\) for every root \(z\). At a root,
\(z^m U_m(1/z)=-L_m(z)\), and differentiation gives
\[
 \frac{zP_m'(z)}{L_m(z)}
   =\frac{zL_m'(z)}{L_m(z)}
      +\frac{z^{-1}U_m'(1/z)}{U_m(1/z)}-m .
\]
The two fractions are bounded in absolute value independently of
the root and of \(m\). The derivative is therefore nonzero for
all roots for sufficiently large \(m\). Positive rescaling preserves
multiplicities. Since \(M_m(1,t)\) has rational coefficients, its
Bézout identity with its derivative transfers squarefreeness to
every characteristic-zero field.

## Coefficient rings and the finite range of remaining orders

For \(m\ge m_0\), the homogeneous forms \(M_m,M_{m+1}\) have no
common projective zero. At \([0:1]\) their pure factorial coefficients
are nonzero, and on \(A\ne0\) this is the separation theorem. They form
a regular sequence in \(\mathbb Q[A,B]\) with degrees \(m,m+1\).
The two-form Hilbert series is
\[
 (1+t+\cdots+t^{m-1})(1+t+\cdots+t^m).
\]
The quotient has dimension \(m(m+1)\) over \(\mathbb Q\), has a
nonzero component in degree \(2m-1\), and vanishes in all higher
degrees. These are the same homogeneous two-form consequences as in
[the ordinary binomial argument](ORDINARY_BINOMIAL_MOMENT_SCHEMES_AND_GAUSSIAN_TRANSFER.md).

After base change to any commutative \(\mathbb Q\)-algebra, vanishing
at orders \(m,m+1\) therefore forces \((A,B)^{2m}=0\) and
\(F^{2m}=0\). Both exponents are sharp in the universal quotient.
For the phase itself, the monomials
\(X^{(2m-1-k)\alpha+k\beta}\) for \(0\le k\le 2m-1\)
are distinct, so a nonzero coefficient of degree \(2m-1\) survives.
Conversely a nilpotent coefficient ideal makes all sufficiently
high moments vanish. Thus eventual vanishing of the moments is
equivalent to nilpotence of \((A,B)\). Over a characteristic-zero
field the only phase in this fixed family with all powers in
\(\ker\Lambda\) is zero.

For a fixed support and any certified cutoff \(m_0\), checking the
remaining coprimality statements reduces to rational resultants at
\(1\le m<m_0\). The separation theorem supplies a support-dependent
existence cutoff.

## The universal coefficient scheme at terminal four

The [terminal-four theorem](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#terminal-4-for-arbitrary-supports)
gives the same homogeneous complete intersection at \(m=3\)
for every pair of distinct monomials, including supports without
private variables. For every fixed tuple of positive real Gamma
shapes, work over \(\mathbb R\) with the two moment forms at orders
three and four. Their pure coefficients are nonzero and their
dehomogenizations are coprime over \(\mathbb C\). The Hilbert-series
argument above gives dimension twelve, highest nonzero degree five,
and sharp coefficient-ideal and polynomial nilpotence at power six
after evaluation in a commutative \(\mathbb R\)-algebra. At ordinary
shapes the construction and the same conclusions hold over
\(\mathbb Q\) and arbitrary commutative \(\mathbb Q\)-algebras.

## Individual roots at fixed order and small curvature

Fix \(n\ge1\). In the
[small-curvature limit](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#gaussian-limit-of-the-centered-roots),
the normalized order-\(n\) polynomial tends coefficientwise, uniformly
over the supports, to an imaginary rescaling of \(\mathrm{He}_n\).
Every root of \(\mathrm{He}_n\) is simple. Indeed a repeated root
would also annul \(\mathrm{He}_{n-1}\) by the derivative formula
\(\mathrm{He}'_n=n\mathrm{He}_{n-1}\). The Hermite recurrence
\(\mathrm{He}_{j+1}=x\mathrm{He}_j-j\mathrm{He}_{j-1}\) would descend
to \(\mathrm{He}_0=1\), a contradiction. Both identities follow
by differentiating the generating function in its two variables.

The limiting leading coefficient and discriminant are nonzero.
By continuity, the normalized order-\(n\) polynomial is squarefree
for sufficiently small positive curvature, with a cutoff depending
only on \(n\). The real affine normalization preserves multiplicities,
and the original ordinary moment polynomial has rational coefficients.
Consequently the squarefreeness also holds over every
characteristic-zero field.
