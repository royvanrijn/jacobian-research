# Integer anchors in the deep pure-cross residue

Let \(f=xy+\alpha x^2/2+\gamma y^2/2+dx+ey\) and set
\(R_n=\Lambda_2(f^n)/(n!)^2\). Work over a finite extension of
\(\mathbb Q_2\), with \(v(2)=1\) and \(v(0)=+\infty\).
Suppose all four noncross natural coefficients have valuation
strictly above one. Define
\[
                  C=\gamma d^2-2de+\alpha e^2 .
\]

**Theorem.** If \(C=0\), every \(R_n\) is a unit. If \(C\ne0\),
fix an integer \(a\ge-1\) and \(B\ge1\). For all sufficiently
large integers \(L\), any positive integer \(n\ne a\) with
\(n\equiv a\bmod 2^L\) and \(n\le 2^{BL}\) satisfies
\((R_n,R_{n+1})\ne(0,0)\).
The threshold can depend on \(a,B\) and the coefficients.

## The normalized recurrence

The [interior and boundary identity](ZERO_CONSTANT_BINARY_QUADRATIC_MATCHING_AND_BOUNDARY_LIFTS.md#the-exact-interior-and-boundary-recurrence)
gives
\[
 D R_n=2\Delta R_{n-1}+\frac{C}{n(n-1)}R_{n-2}
              +\frac{T_n}{n!},\qquad
 D=\alpha+\gamma-2,\quad \Delta=\alpha\gamma-1
 \quad(n\ge2).
\]
Here \(v(D)=1\), \(v(\Delta)=0\).
If at least one noncross coefficient is present, let \(s>1\)
be the least of their valuations. The axis recurrence gives
\(v(A_n),v(B_n)\ge ns\) for the normalized axis moments.
The weights defining \(T_n\) therefore give \(v(T_n)\ge ns\).
Legendre's formula yields
\[
                 v(T_n/n!)\ge(s-1)n+s_2(n),
\]
where \(s_2\) is the binary digit sum. In the case \(C=0\),
the recurrence reduces modulo the maximal ideal to
\(R_n=(2\Delta/D)R_{n-1}\), with unit multiplier and \(R_1\)
a unit. The phase \(xy\) alone gives \(R_n=1\) directly.

## The two-window at integer anchors

Take the distinguished column \((1,1)\) in the
[endpoint normalization](FACTORIAL_ENDPOINT_NORMALIZATION_AND_INTEGER_ANCHORS.md).
Its slope preference holds: the quadratic coefficients have
valuation at least \(s-1>0\), and the linears lose one in
the slope and still have positive excess. At every negative
integer \(a\) write \(V_a\) for the convergent anchor series.
At nonnegative \(a\), write \(G(a)\) for the endpoint series
after multiplication by \(n-a\).
Indeed \(C\ne0\) implies two rows can decrease in one count
direction. This holds if \(de\ne0\) by using both linears.
If \(\gamma d^2\ne0\), using one \(y^2\) and two \(x\)
terms decreases both rows, and the other case is symmetric.

We first show that \(V_{-1}\) is a unit. In the exact count
formula at \(n=-1\), any positive row shift gives a zero factor
at the start of its rising product. For a surviving count
with \(h\) replacements let \(L_x,L_y\) be the decreases.
Their sum is the total number of linears in the count.
Up to sign the rational multiplier of the coefficient ratios
is the ordinary integer count multinomial divided by
\(L_x!L_y!\). The loss from these factorials is at most
\(L_x+L_y\). Each linear thus contributes valuation at least
\(s-1\) after this loss, and each quadratic ratio has the
same lower bound. Every count with \(h>0\) contributes
valuation at least \((s-1)h\). The empty count is one.

For fixed \(a\), approach it by integers \(n\equiv a\bmod 2^L\)
as in the anchor theorem. To obtain the relations it suffices
to take \(n=a+2^L\) and use the anchor theorem with any fixed
exponent bound above one. Multiplying the normalized recurrence
by \(n-a\) and passing to the corresponding limits gives
\[
\begin{split}
     D G(0)&=-C V_{-2},\\
     D G(1)&=2\Delta G(0)+C V_{-1},\\
     D G(a)&=2\Delta G(a-1)+\frac{C}{a(a-1)}G(a-2)
                           \qquad(a\ge2).
\end{split}
\]
The boundary term goes to zero by the valuation bound above.
At \(a=0\) the factor \(n\) cancels and leaves \(1/(n-1)\)
in the term with \(R_{n-2}\), whose limit is \(V_{-2}\).
At \(a=1\) the corresponding cancellation leaves \(1/n\)
and the limit \(V_{-1}\).
The second identity and \(C V_{-1}\ne0\) exclude \(G(0)=G(1)=0\).
Since every second-predecessor coefficient in the third identity
is nonzero, two adjacent zero values would force \(G(0)=G(1)=0\)
by backward propagation. For \(a=-1\) use \(V_{-1}\) alone.
The uniform error estimate at each fixed anchor now transfers the
nonzero value to at least one of the two exact moments in the
claimed range, using a slightly larger exponent bound for \(n+1\)
if needed.
