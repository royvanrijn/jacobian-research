# Singleton-slot windows and exact interpolation for factorial moments

Write
\[
 f=\sum_{\nu}\frac{u_\nu}{\nu!}X^\nu,\qquad
 P_n=\Lambda(f^n)/n!,\qquad E(T)=\sum_{n\ge0}P_nT^n .
\]
The sum has finite support, no constant term, and integral
natural coefficients in a finite extension of \(\mathbb Q_p\).
A singleton is an exponent vector \(p^b e_i\). Throughout
the residue argument, assume that each coordinate has at most
one unit singleton. Let \(I\) be the set of coordinates with
one, denote it by \(p^{b_i}e_i\), write \(a_i=u_{p^{b_i}e_i}\),
and put \(q=|I|\).

The factorial insertion at a singleton depends modulo the
maximal ideal only on the part of the inserted exponent at
or above its digit position. This gives an actual-support
window at every order for an explicit dyadic class of
arbitrary degree and dimension. When all the unit
singletons are linear, the same decomposition gives local
analytic interpolation with no initial transient if each
other monomial meets a unit-linear coordinate.

## 1. The residue series

Use a bar for reduction to the residue field. Let \(H\) be
the set of unit nonsingletons. For a count vector \(m\) on
\(H\), put
\[
 h=\sum_{\nu\in H}m_\nu,\qquad L_i=\sum_{\nu\in H}m_\nu\nu_i,\qquad
 C(m)=
 \frac{\prod_iL_i!}{\prod_{\nu\in H}m_\nu!\prod_{\nu\in H}(\nu!)^{m_\nu}},
 \qquad
 \ell_i(m)=\left\lfloor\frac{L_i}{p^{b_i}}\right\rfloor\quad(i\in I).
\]
Then the residue series is
\[
 \overline E(T)=
 \sum_{\substack{0\le m_\nu<p\\ p\nmid C(m)}}
 \frac{\overline{C(m)}\prod_{\nu\in H}(\bar u_\nu T)^{m_\nu}}
      {\prod_{i\in I}(1-\bar a_iT)^{1+\ell_i(m)}} .          \tag{1}
\]
In the retained terms the coordinate additions and the
single-column digit multiplications are carry-free, so
\[
       \ell_i(m)=
        \sum_{\nu\in H}m_\nu\left\lfloor\frac{\nu_i}{p^{b_i}}\right\rfloor .
                                                               \tag{2}
\]

To prove (1), first fix the nonsingleton counts. By the
[factorial coefficient valuation formula](FACTORIAL_MOMENT_SERIES_PRIME_POWER_RATIONALITY.md#degree-bounds-in-characteristic-p),
\(C(m)\) is integral and is divisible by \(p\) if any
nonsingleton count is at least \(p\). The same formula gives
the carry-free condition whenever \(p\nmid C(m)\).
Terms using a nonunit coefficient vanish in the residue
field. For each coordinate \(i\in I\), write \(Q=p^{b_i}\).
The remaining count \(N\) of its singleton brings the
factor
\[
   (\bar a_i T)^N\,
   \overline{\frac{(QN)!}{N!\,(Q!)^N}}\,
   \overline{\binom{QN+L_i}{L_i}} .                           \tag{3}
\]
The normalized pure singleton ratio is congruent to one
modulo \(p\). Its valuation is zero by Legendre's formula,
and the digit factorials of \(QN\) equal those of \(N\),
so the unit-factorial formula gives the stated residue.
Writing \(L_i=Q\ell_i+r\) with \(0\le r<Q\), Lucas's
formula makes the last factor of (3) equal to
\(\overline{\binom{N+\ell_i}{\ell_i}}\).
Its generating series is \((1-\bar a_iT)^{-1-\ell_i}\).
Multiplication over the coordinates and summation over
the retained counts proves (1).

## 2. All-order windows from the slot loads

Suppose \(q\ge1\), and over the retained count vectors
define
\[
  M_i=\max_m\ell_i(m),\qquad
  R=\sum_{i\in I}M_i,\qquad
  B=\max_m\!\left(h-\sum_{i\in I}\ell_i(m)\right).
\]
The empty vector is retained. Formula (1) has common
denominator \(\prod_{i\in I}(1-\bar a_iT)^{1+M_i}\),
of degree \(q+R\), and numerator of degree at most \(R+B\).
In particular, if \(B<q\), this is a nonzero proper
rational series because its constant coefficient is one.
After any cancellation the denominator still has positive
degree at most \(q+R\). Its leading coefficient is nonzero,
so the coefficient recurrence propagates any zero block
of its degree back to the constant coefficient. Therefore
every window of length \(q+R\) contains a unit normalized
moment. When \(q+R\) is at most the actual support size,
this proves the actual-support Strong Factorial window.

At \(p=2\) the retained count vectors are the subsets of
\(H\) with disjoint binary digit slots, and
\(\overline{C(m)}=1\). For example, suppose every unit
nonsingleton has
\[
 \left\lfloor\nu_i/2^{b_i}\right\rfloor=1
 \quad\hbox{at exactly one coordinate }i\in I,
\]
with all the other quotients at coordinates in \(I\)
equal to zero. Then \(B=0\). Two selected
nonsingletons that hit the same coordinate would both
occupy digit \(b_i\) there. Thus \(M_i\le1\), with equality
exactly when some unit nonsingleton hits that coordinate.
The resulting window length is \(q+\#\{i:M_i=1\}\),
at most \(q+|H|\). Coefficients may be arbitrary units
and the exponents below the singleton positions and in
coordinates outside \(I\) are unrestricted. In more
general profiles the maxima \(M_i\) use the actual
carry-free subsets, so high slot loads that compete for
a digit do not have to be added across incompatible
subsets.

## 3. Exact interpolation with unit linears

Assume now that each unit singleton is \(e_i\) at a
coordinate \(i\in I\). Keep every other nonzero column
in an index set \(H'\), including nonunit singletons and
nonunit nonsingletons. For any count vector on \(H'\),
use \(h,L_i,C(m)\) as above with \(H'\) in place of \(H\).
Grouping the moment series by these counts gives the exact
coefficientwise identity
\[
 E(T)=\sum_{m\ge0}
           C(m)\prod_{\nu\in H'}(u_\nu T)^{m_\nu}
             \prod_{i\in I}(1-a_iT)^{-1-L_i}.                \tag{4}
\]
Indeed the singleton ratio in (3) is now exactly one,
and the insertion is exactly \(\binom{N+L_i}{L_i}\).

At precision \(p^k\) only counts with \(h=O(k)\)
are retained. Each nonsingleton contributes at least
\(v_p(m_\nu!)\ge m_\nu/p-1\) to \(v(C(m))\)
before coefficient valuations. Each nonunit singleton
contributes at least \(m_\nu v(u_\nu)\) from its
coefficient. More explicitly, take
\(\lambda\) to be the minimum of \(1/p\) and all the
\(v(u_\nu)>0\) for nonunit singletons, and let \(g\)
be the number of nonsingletons. The valuation
of the scalar prefactor in (4) is at least
\(\lambda h-g\).
Every term with \(h\ge(k+g)/\lambda\) can be discarded
modulo \(p^k\). Since \(L_i\le h\max_{\nu\in H'}\nu_i\),
the common denominator degree at this precision is linear
in \(k\). Clearing that denominator gives a numerator
of degree at most its degree plus the maximum retained
\(h\), also linear in \(k\). The reciprocal roots of the
denominator are fixed units. The
[finite-difference interpolation argument](QUADRATIC_DYADIC_ANALYTIC_UNIT_ENVELOPE.md#finite-differences-and-local-analyticity)
therefore applies with \(p\) in place of two, giving
locally analytic unit envelopes and a threshold linear
in \(k\). For the analytic convergence, on a disk
\(z=b+p^j y\) the Gauss loss for \(\binom zr\)
is at most \(r/(p^j(p-1))\), by summing the
factorial valuations at levels above \(j\).
Choosing \(j\) large enough makes this loss smaller
than the linear growth of the Mahler coefficient
valuations.

Suppose in addition that \(q\ge1\) and every
\(\nu\in H'\) has \(\sum_{i\in I}\nu_i\ge1\).
Then \(h\le\sum_{i\in I}L_i\), so every rational
summand in (4) is proper. The common numerator stays
proper at each precision. Consequently the recurrence
already holds on the entire nonnegative sequence and
its compatible unit envelopes interpolate \(P_n\)
**exactly at every nonnegative integer**.
The finitely many residue classes are those for a
fixed step \(L\) making every \(a_i^L-1\) have
valuation at least one. The same recurrences extend
the sequence to negative integers.

This exact interpolation has a forced negative block.
In the covered case just stated, if \(q\ge2\), the
interpolated values at \(-1,\ldots,-(q-1)\) vanish.
Indeed a proper rational summand in (4) is
\(O(T^{-q})\) at infinity, since its denominator
degree minus its numerator degree is at least \(q\).
For each precision let \(N(T)/D(T)\) be the proper
rational approximation. The backward values of its
unit recurrence are the negatives of the coefficients
of its expansion at infinity. To see this, combine the
series at zero with minus the series at infinity as a
bilateral coefficient sequence. Multiplication by \(D\)
gives \(N-N=0\), and the unit leading coefficient of
\(D\) gives the unique backward recurrence. The
coefficients of \(T^{-1}\) through \(T^{-q+1}\)
are zero at each precision, proving the assertion.

The exact interpolation and the residue window apply to
the stated singleton configurations. Profiles without
unit singletons, or with several different singleton
positions at one coordinate, can require other arguments.

## 4. Higher-precision windows from suppressed multiple insertions

Choose a unit linear of the phase \(X_i\). Scale its
coefficient to one and assume all remaining natural coefficients
\(u_\nu\) are integral. The following argument separates
the arrays with zero or one insertion from those with several.
It applies even if those other columns include unit singletons.

For \(j\ge2\) put
\[
     C_\nu(j)=\frac{\prod_i(j\nu_i)!}{j!(\nu!)^j},\qquad
     K_{\mu,\nu}=\sum_i v_p\binom{\mu_i+\nu_i}{\mu_i}.
\]
For an integer \(k\ge1\), impose the following inequalities:
\[
\begin{array}{ll}
 v_p(C_\nu(j))+j\,v(u_\nu)\ge k
   &\text{for every }\nu\text{ and }2\le j<p,\\
 w_\nu-1+p\,v(u_\nu)\ge k
   &\text{for every }\nu,\\
 K_{\mu,\nu}+v(u_\mu)+v(u_\nu)\ge k
   &\text{for every distinct pair }\mu,\nu.
\end{array}                                                    \tag{5}
\]
Here \(w_\nu\) is the sum of all base-\(p\) exponent digits.
At \(p=2\), the first line is empty, and \(K_{\mu,\nu}\)
counts the binary carries with their multiplicities.
For unit coefficients the second line then asks for weight
at least \(k+1\).

Under (5), let \(A_1,\ldots,A_t\) be the distinct **positive**
exponents \(\nu_i\) among the remaining columns. Set
\[
 B_0=1+\sum_{\nu_i=0}u_\nu,\qquad
 B_j=\sum_{\nu_i=A_j}u_\nu\quad(1\le j\le t).
\]
If
\[
                 v(B_0)+\sum_{j=1}^t v_p(A_j)<k,               \tag{6}
\]
then every window of length \(t+1\) at a positive index contains
a nonzero exact moment. This is at most the actual support size,
and can be shorter when several columns share the same
exponent along \(X_i\).

To prove the suppression, group the exact expansion by the
counts of all columns other than \(X_i\), keeping only \(X_i\)
as the distinguished linear. The identity (4) still holds
with this single distinguished column. Factor the scalar
\(C(m)\) into the product of the single-column factors
\(C_\nu(m_\nu)\) and the coordinate multinomial factors
\(\prod_i\binom{\sum_\nu m_\nu\nu_i}{(m_\nu\nu_i)_\nu}\).
Every factor is integral. If some \(m_\nu=j\ge2\), the
corresponding scalar weighted by \(u_\nu^j\) has valuation at
least \(k\). For \(j<p\) this follows from the first line
of (5). For \(j\ge p\), use
\(v_p(C_\nu(j))\ge (w_\nu-1)v_p(j!)\ge w_\nu-1\) and the
second line. If all counts are zero or one and at least two
are one, merge any selected pair first in each coordinate
multinomial. This extracts the pair of binomial factors
whose total valuation is \(K_{\mu,\nu}\).
The last line of (5) suppresses that array.
The remaining generating series therefore gives, at \(n\ge1\),
\[
         P_n\equiv B_0+
                \sum_{j=1}^t B_j\binom{n+A_j-1}{A_j}
                                                    \pmod{p^k}. \tag{7}
\]

For the window argument let \(\Delta\) be the forward
difference in \(n\), let \(Z^{\underline r}\) denote the
falling factorial, and write
\[
       H(Z)=\prod_{j=1}^t(Z-A_j)=\sum_{r=0}^t\gamma_r
                                               Z^{\underline r}.
\]
The \(\gamma_r\) are integers because the falling factorial
basis is monic and triangular. With \(n^{\overline r}\) the
rising factorial, repeated Pascal gives the exact identity
\[
 n^{\overline r}\Delta^r\binom{n+A-1}{A}
                =A^{\underline r}\binom{n+A-1}{A}.
\]
It holds also for \(r>A\), when both sides vanish.
Consequently the integral operator
\(\sum_{r=0}^t\gamma_r n^{\overline r}\Delta^r\) applied to (7)
gives \(B_0H(0)\) modulo \(p^k\). The operator uses only
\(P_n,\ldots,P_{n+t}\), and (6) makes its value nonzero.

More precisely put \(m_n=\min_{0\le r\le t}
v_p(\gamma_r n^{\overline r})\), taking \(v_p(0)=+\infty\).
The error after applying the operator has valuation at least
\(k+m_n\). Thus a proposed exact zero window must satisfy
\(k+m_n\le v(B_0H(0))\).
For \(t=1\) and \(B_0=1\) this implies
\(v_p(A_1)-v_p(n)\ge k\), since
\(m_n=\min(v_p(n),v_p(A_1))\).
In particular a dyadic binomial with one linear monomial,
positive shared exponent \(A\) on the other monomial, and
integral relative natural coefficient has all consecutive
nonvanishing if the other monomial has binary weight
\(w_\nu>v_2(A)+1\). The pair condition is vacuous for this
binomial. For many columns sharing \(A\), the same
two-window holds under (5) with \(k>v_2(A)\).
For example, if \(A\equiv2\bmod4\), binary weight at least
three is sufficient for the binomial.
