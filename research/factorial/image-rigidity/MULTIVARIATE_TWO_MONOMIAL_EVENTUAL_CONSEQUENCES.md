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
[the ordinary binomial argument](ORDINARY_BINOMIAL_RADIAL_TRANSFERS.md).

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
