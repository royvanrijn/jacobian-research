# A cross-parity link for zero-constant binary quadratics

Let \(f=\alpha X^2/2+XY+\gamma Y^2/2+dX+eY\) and
\(P_n=\Lambda_2(f^n)/n!\), where
\(\Lambda_2(X^iY^j)=i!j!\). Work in a finite extension of
\(\mathbb Q_2\), with \(v(2)=1\) and \(v(0)=+\infty\). Suppose
\[
             \tfrac12<s=\min(v(\alpha),v(\gamma),v(d),v(e))
                         <\tfrac34 .
\]
Choose \(c\) to be a coefficient attaining \(s\), and use hats for
the reductions of the four coefficients divided by \(c\). Define
\[
D=(1+\widehat\alpha T+\widehat d^{\,2} T^2)
  (1+\widehat\gamma T+\widehat e^{\,2} T^2)
  =1+l_1T+l_2T^2+l_3T^3+l_4T^4 .
\]
A nonzero cross coefficient of minimal valuation can be normalized
to one before applying this result. The hypotheses refer to
relative valuations after that normalization.

In this sector, every four consecutive positive moment indices
contain a nonzero moment in any of the following cases.

- \(\deg D\le2\).
- \(\deg D=3\) and \(l_1\ne0\).
- \(\deg D=4\) and \(l_1l_3\ne0\).

For \(\deg D=3\) with \(l_1=0\), the four-window still holds at
even starts, and the five-window holds at every start. For
\(\deg D=4\) with exactly one of \(l_1,l_3\) nonzero, the
five-window holds at odd starts. The argument works over any
residue field of characteristic two and permits fractional
valuations.

## The even coefficient layer

Use the cross expansion and digit identity from
[weighted matching controls](BINARY_ZERO_CONSTANT_QUADRATIC_WEIGHTED_MATCHING_LOCAL.md#exact-cross-weights-below-the-unit-valuation).
The \(j\)-cross part of \(\sum_n P_n(T/c)^n\) is
\(j!c^{-j}T^j J_{j,X}J_{j,Y}\). Here the vertex coefficients
inside each \(J_j\) are scaled by \(c\), and
\[
\overline{J_{j,X}}=
  (1+\widehat d T)^{1-(j\bmod2)}
   (1+\widehat\alpha T+\widehat d^{\,2}T^2)^{-1-\lfloor j/2\rfloor}.
\]
The same identity holds for \(Y\). The multiplier has valuation
\(\phi(j)=v_2(j!)-js\). Its binary bit of place \(b\) has cost
\((1-s)2^b-1\). The bits at places zero and one have negative cost
and every higher bit has positive cost. Thus \(j=3\) uniquely
minimizes the unrestricted cost. At an odd index, a zero moment
forces the corresponding coefficient of \(T^3/D(T)^2\) to
vanish, as in the leading matching control.

At an even index the smallest possible contribution comes
uniquely from \(j=1\). Indeed, every even \(j\) has cost
at least \(\phi(2)=1-2s>\phi(1)\). Every odd \(j\)
whose second bit is zero has cost at least \(\phi(1)\),
strictly except at \(j=1\). Finally, for odd \(j\) whose
second bit is one, the displayed digit identity is even in
\(T\) at both vertices. The coefficient of odd degree in
\(J_{j,X}J_{j,Y}\) has integer-polynomial coefficients
divisible by \(2\): the digit formula holds in the polynomial
ring in the four scaled parameters over \(\mathbb F_2\).
The valuation gain is therefore at least \(1\), including in
a ramified extension. Such a cross part at an even index
has valuation at least \(\phi(3)+1=2-3s>\phi(1)\).
It follows that at any positive even index, a zero moment
forces the corresponding coefficient of \(T/D(T)\) to
vanish.

Put
\[
  \sum_{k\ge0} a_k U^k
       =\frac1{1+l_1^2U+l_2^2U^2+l_3^2U^3+l_4^2U^4},
  \qquad a_k=0\quad(k<0).
\]
Since \(1/D(T)^2=\sum a_k T^{2k}\) and
\(1/D(T)=D(T)\sum a_kT^{2k}\), the two moment tests are
\[
 P_{2r+1}=0\ \Longrightarrow\ a_{r-1}=0,\qquad
 P_{2r+2}=0\ \Longrightarrow\ l_1a_r+l_3a_{r-1}=0.          \tag{1}
\]

## Forced boundary values

Write \(M=\sum_n P_nT^n\), and let \(A,B\) be the normalized
moment series of the restrictions of \(f\) to the \(X,Y\)
axes, respectively. Put
\[
 L(T)=\big(1-\alpha+(d-\alpha e)T\big) A(T)
      +\big(1-\gamma+(e-\gamma d)T\big) B(T).
\]
With \(\Delta=\alpha\gamma-1\) and
\(C=2de-\gamma d^2-\alpha e^2\), the integration-by-parts
identity for the factorial functional gives
\[
 2\Delta T^2 M'
   =(\alpha+\gamma-2-2\Delta T+CT^2)M+L.                  \tag{2}
\]
To see this directly, take \(U_X=\Lambda_2(X\exp(Tf))\) and
\(U_Y=\Lambda_2(Y\exp(Tf))\). The identity
\(\Lambda_2(\partial_X g)=\Lambda_2(g)-\Lambda_1(g(0,Y))\)
and its \(Y\)-version give
\[
 T(\alpha U_X+U_Y+dM)=M-B,\qquad
 T(U_X+\gamma U_Y+eM)=M-A .
\]
Applying the same identity to \(X\exp(Tf)\) and
\(Y\exp(Tf)\) and using
\(2f=X f_X+Y f_Y+dX+eY\) gives
\(2TM'=(1+dT)U_X+(1+eT)U_Y-2M\).
Substitution yields (2). In particular, a block of zero
moments at \(m,\ldots,m+k-1\) forces \(L_n=[T^n]L=0\) at
\(m+2,\ldots,m+k-1\).

The vertex identity at \(j=0\) shows that \(L(T/c)\) is
integral and has residue
\[
 \frac{(1+\widehat dT)^2}{1+\widehat\alpha T+\widehat d^{\,2}T^2}
 +\frac{(1+\widehat eT)^2}{1+\widehat\gamma T+\widehat e^{\,2}T^2}
 = T D'(T)/D(T).
\]
Its coefficient at \(2r+3\) is
\[
 l_1 a_{r+1}+(l_3+l_1l_2)a_r
       +(l_3l_2+l_1l_4)a_{r-1}+l_3l_4 a_{r-2}.           \tag{3}
\]
A vanishing exact \(L_{2r+3}\) forces this expression to
vanish.

## Combining the two parities

Let \(h=\deg D\). The reciprocal series defining the \(a_k\)
has degree-\(h\) recurrence with nonzero trailing coefficient
\(l_h^2\) and \(a_0=1\). A block of \(h\) successive
nonnegative indices with zero coefficients propagates
backward to \(a_0=0\). A block used below that starts at a
negative index contains \(a_0\).

The leading matching bound already gives the four-window
when \(h\le2\). Consider an even start \(2r+2\) for a
four-zero block. By (1), \(a_r=a_{r+1}=0\) and
\(l_3a_{r-1}=0\). Equation (3) at \(2r+5\), a forced
boundary index, gives \(l_1a_{r+2}+l_3l_4a_{r-1}=0\).
When \(h=4\) and \(l_1l_3\ne0\), this yields the four
consecutive zeros \(a_{r-1},\ldots,a_{r+2}\). When \(h=3\),
one has \(l_3\ne0\), so \(a_{r-1},a_r,a_{r+1}\) already
contradict the degree-three recurrence.

For an odd start \(2r+1\), use \(r\ge1\), since \(P_1\)
is a unit. Four zero moments give \(a_{r-1}=a_r=0\)
and \(l_1a_{r+1}=0\). At \(2r+3\), (3) then gives
\(l_3l_4a_{r-2}=0\). If \(h=4\) and \(l_1l_3\ne0\),
this yields the four consecutive zeros
\(a_{r-2},\ldots,a_{r+1}\). If \(h=3\) and \(l_1\ne0\),
the three consecutive zeros starting at \(a_{r-1}\)
suffice.

Five zero moments at the same odd start already give
\(a_{r-1}=a_r=a_{r+1}=0\) by (1), which suffices
for \(h=3\). For \(h=4\), the boundary test at
\(2r+3\) forces \(l_3l_4a_{r-2}=0\), while that at
\(2r+5\) forces \(l_1a_{r+2}=0\). Hence either
\(l_1\ne0\) or \(l_3\ne0\) gives the fourth zero
coefficient.

The remaining quartic cases in this valuation band
are even starts when exactly one of \(l_1,l_3\) is zero,
and all starts when both are zero. The second-layer and
boundary tests alone can permit zero blocks there.
For instance, take \(D=1+T+T^4\) over \(\mathbb F_2\)
and split it into two vertex quadratics over a finite
extension of the residue field. Here \(a_{12}=a_{13}
=a_{14}=0\), by the recurrence
\(a_k=a_{k-1}+a_{k-4}\) with \(a_0=1\).
Thus the tests at moment indices \(26,\ldots,30\)
and boundary indices \(28,\ldots,30\) all vanish.
Those are residue conditions at two specified precisions.
The test does not establish an exact zero moment over
the coefficient field.
