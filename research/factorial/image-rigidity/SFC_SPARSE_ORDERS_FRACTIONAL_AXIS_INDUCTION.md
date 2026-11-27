# Fractional axis induction for prescribed factorial orders

Write \(L(x^i y^j)=i!j!\). An axis profile of size \(s\) fixes
offsets \(u_i\in\mathbb Z_{\ge0}\) and assigns each column to
one of two sides, with exponents \((u_i,D-u_i)\) on the first
side and \((D-u_i,u_i)\) on the second. Offsets on each side
are distinct. Use only ordinary integer degrees \(D\) at
which all exponents are nonnegative and the columns are
distinct.

**Fixed-axis theorem.** For any \(n\ge1\) and any finite list
of size-\(s\) axis profiles with \(s\ge2\), one can choose the
same \(s\) orders for all the profiles, beginning with
\(n,n+1\) and followed by increasing powers of two. The common
factorial zeros of those orders can occur only at a computable
finite set of degrees for each profile. The result holds over
every characteristic-zero field. The same theorem holds
in any number of variables for profiles with at most two
columns at each total fixed offset, assigned to different
growing coordinates whenever there are two.

**Consecutive trinomial case.** For an axis trinomial with
two columns on one side and one on the other, suppose the
unique largest fixed offset belongs to the two-column side.
Then for every \(n\ge1\), the actual consecutive orders
\(n,n+1,n+2\) have only a computable finite set of possibly
exceptional degrees. At every start, condition (2) also
covers a sufficiently separated unique maximum over a
same-side survivor pair. Successive separation also gives
consecutive support-sized orders for multivariable axes. At even starts this degree sieve
also holds when the mark leaves a same-side pair with
gap larger than the binary length bound below. The
same bound, increased by one, also covers a tied maximum.
The two even marks narrow the remaining wide-gap baselines
to (3). Condition (4) also treats deep upper-endpoint
baselines at arbitrary gaps.

**Independent heights.** Fix a finite list of bivariate
homogeneous offset supports, each with \(s\ge2\) columns.
For each start \(n\ge1\), at most \(2s-2\) computable orders,
including \(n,n+1\), restrict all common factorial zeros on
all independent nonnegative integer translates of those
supports to a computable finite set of height pairs. The
same result holds for inhomogeneous offset supports if each
exponent row has distinct entries and the Gaussian weights
specified in Section 4 are distinct. The exact exceptional
height pairs for the selected orders are computable. At
start two and \(s\ge3\), homogeneous supports admit at most
\(2s-3\) orders including \(2,3,4\). For \(s\ge3\) the
\(2s-2\) orders can include \(n+2\) at even starts and
\(n+3\) at odd starts.

Edo and van den Essen formulate the Strong Factorial
Conjecture with consecutive orders in
[Section 2.1 of their paper](https://arxiv.org/abs/1304.3956).
Liu and Sun's
[Theorem 2.6](https://doi.org/10.1017/S0004972719000546)
uses all positive orders to prove homogeneous FC(2).
The fixed-axis theorem interpolates between these scopes by
prescribing the first adjacent pair and allowing the later
orders to be chosen. The height theorem gives four orders for
trinomials and does not decide three consecutive moments.

## 1. Fractional elimination on fixed axes

For a count vector \(\alpha\) with \(|\alpha|=m\) supported
on one side, set \(U=\sum_i\alpha_i u_i\). Define the
decoupled form \(P_m(D)\) by retaining only these count
vectors, with coefficient
\[
                 \binom m\alpha\frac{U!}{(mD)_{\underline U}},
\]
where \((T)_{\underline U}=T(T-1)\cdots(T-U+1)\).
At ordinary degrees these are the same-side terms of
\(L(f^m)/(mD)!\). This is the rational continuation in the
[initial axis sieve](SFC_HOMOGENEOUS_TRINOMIAL_AXIS_DEGREE_SIEVE.md).

Set \(M=\max_i u_i\). When \(M\ge1\), a terminal order
\(q\) gives a pole at \(D_*=M-1/q\). The denominator at
order \(q\) has a zero only for \(U=qM\), since the
vanishing factor \(qD-(qM-1)\) is last in the largest
denominator. Within a side this occurs only on the pure
term of its maximal column, and only if that column has
offset \(M\). The residues coincide on the two sides.
Thus \((D-D_*)P_q(D)\) specializes to a nonzero multiple
of the sum of the pure \(q\)-th powers of the maximal
columns. All retained coefficients at each \(m<q\) are
regular and strictly positive there, since each denominator
factor is at least \(1-m/q\).

For two columns on the same side, use \(q=n+1\). The larger
offset is unique. Its coefficient variable must vanish at
the mark, and the \(n\)-th equation forces the other to
vanish. The decoupled resultant in \(D\) is therefore a
nonzero rational function. For two columns on opposite
sides, at every real \(D>M\) the two forms are diagonal
with strictly positive coefficients. No ratio can have
both its \(n\)-th and \((n+1)\)-st powers negative real.
Their quotient would make the ratio positive real.
This proves the two-column assertion even when \(M=0\).

Inductively choose orders for the finite family of all
subprofiles, keeping the same previous orders whenever
a column is deleted. Suppose the \(r-1\) selected forms
have a nonzero rational resultant in \(D\) for every
\((r-1)\)-column subprofile. Consider an \(r\)-column
profile and a prospective power of two \(q\) higher than
all the selected orders. There are finitely many values
\(D_*=M-1/q\) at which the lower resultant on any specified
subprofile vanishes or has a pole. Indeed the values are
distinct and the relevant rational functions are nonzero.

If the maximum is unique, the terminal mark deletes its
column and the lower equations then exclude a common
direction outside those finitely many \(q\). In the
consecutive trinomial case the two remaining columns are
on opposite sides. Their lower equations are diagonal
and strictly positive at \(D_*\) itself, so \(q=n+2\)
works directly.

If the maximum is tied, use variables \(Y,Z\) for the two
maximal columns. They are on opposite sides. The terminal
equation gives \(Y^q+Z^q=0\). At \(Y=0\) it forces \(Z=0\),
and the lower equations already exclude a common direction
on \(Y=0\) for all but finitely many \(q\), by deleting the
\(Y\) column. For \(Y\ne0\), set \(Z=\xi Y\). The resultant
\(R(\xi)\) of the \(r-1\) lower equations in the \(r-1\)
remaining variables is a rational-coefficient polynomial
in \(\xi\). It is nonzero. To see this homogenize the
substitution as \(Y=\eta V,\ Z=\xi V\), with \(V\) a
coefficient variable. At \((\eta,\xi)=(0,1)\) the
resultant is that of the subprofile with \(Y\) deleted
and is nonzero. Dehomogenizing at \(\eta=1\) cannot
annihilate the whole polynomial: by the resultant
change-of-variables law, scaling \((\eta,\xi)\) scales
the resultant by the same power as scaling the single
variable \(V\). Equivalently its degree in \((\eta,\xi)\)
is \(\prod_{j<r}m_j\), where \(m_j\) are the lower
orders. For the exclusion it suffices to use the crude
degree bound
\[
                   \deg R\le(r-1)\prod_{j<r}m_j
\]
from coefficient degrees and resultant multihomogeneity.
Every solution of \(\xi^q=-1\) has minimal polynomial
\(T^q+1\) over \(\mathbb Q\) when \(q\) is a power of two:
it is a primitive \(2q\)-th root of unity. Taking \(q\)
larger than the displayed bound excludes all such roots.

This proves simultaneous nonvanishing at the mark for
every sufficiently large power of two at the \(r\)-th
step. All requirements can be tested with exact rational
arithmetic. In particular one can specialize the marked
forms and compute their projective common-zero test
while searching the powers of two. A nonvanishing mark
proves that the unspecialized decoupled resultant is
nonzero, since multiplication by \(D-D_*\) is invertible
over \(\mathbb Q(D)\).

### Every start with a separated unique maximum

Suppose the unique maximum is \(M\), the two survivors
are on one side with offsets \(a<b\), and set
\(g=M-b\), \(\delta=b-a\). For every start \(n\ge1\)
put \(q=n+2\). If
\[
        \#\{1\le h\le g:\gcd(h,q)=1\}>\delta n(n+1),
        \tag{2}
\]
the decoupled resultant at the consecutive orders
\(n,n+1,n+2\) is nonzero. In particular it suffices
that \(g\ge q(\lfloor\delta n(n+1)/\varphi(q)\rfloor+1)\).
The bound does not depend on the common offset \(a\).

At \(D_h=M-h/q\) the order-\(q\) denominator has a
vanishing factor at \(qD_h=qM-h\). The pure maximal
column has a simple pole. Every other same-side term
has \(U\le qM-g\le qM-h\), so its denominator stops
before that factor. Multiplication by \(D-D_h\)
therefore deletes the maximal column at each of
these marks. The lower orders are regular since
\(h/q\) has denominator \(q\).

For \(m=n,n+1\), divide the surviving binary form by
its constant pure coefficient and clear the remaining
denominator \((m(D-a))_{\underline{m\delta}}\).
In its coefficient of \(X^{m-k}Y^k\) the degree in
\(D\) is at most \(\delta(m-k)\). The resultant of
the two binary forms is nonzero as a polynomial in
\(D\). Indeed, for the pair itself the order-\((n+1)\)
mark at \(D=b-1/(n+1)\) deletes the higher column by
its pole, and the lower pure coefficient is nonzero
there. The normalization and denominator clearing
are by nonzero rational functions in \(D\), which
preserve generic nonvanishing. The degree of the
cleared resultant is at most \(\delta n(n+1)\):
scaling \(X\) scales the binary resultant with
weight \(n(n+1)\), so every coefficient monomial in
it has that total weight when the coefficient of
\(X^{m-k}Y^k\) has weight \(m-k\).
On the marks \(D_h\) we have \(D_h>b\), so all these
normalizations are regular and invertible. The
distinct marks in (2) cannot all be zeros of the
resultant. A nonzero mark gives the claimed
consecutive-order degree sieve.

The same argument gives support-sized consecutive
orders in any dimension when the total fixed
offsets are successively separated. Use the
multivariable axis profiles of Section 2 and write
the distinct total fixed offsets as
\(u_1<\cdots<u_s\). For \(r\ge3\), put
\(W=u_{r-1}-u_1\), \(g=u_r-u_{r-1}\), and
\(q=n+r-1\). Suppose at each such \(r\) that
\[
 \#\{1\le h\le g:\gcd(h,q)=1\}
       >(r-1)W\prod_{m=n}^{q-1}m .
\]
Induct on the number of columns, using the
two-column result of Section 1. For the first
\(r-1\) columns, multiply every order-\(m\) form by
\((mD)_{\underline{mu_1}}\) and then clear the
remaining denominator with
\((m(D-u_1))_{\underline{mW}}\). Each coefficient
now has degree at most \(mW\), independently of
the individual fixed exponents. Resultant
multihomogeneity bounds the degree of this
nonzero lower resultant by the right side of
the inequality. At \(D=u_r-h/q\), we have
\(D>u_{r-1}\), and the terminal residue again
deletes only the maximal column. One of the
marks avoids all roots of the lower resultant.
Section 2 gives the effective ordinary degree
cutoff for the \(s\) consecutive orders
\(n,\ldots,n+s-1\).

### Even starts and wide offsets

For an even start \(n\), mark the order \(q\) divisible
by four in \(\{n,n+1,n+2\}\). Write \(m_e,m_o\) for the
other even and odd orders, so \(v_2(m_e)=1\) and
\(v_2(m_o)=0\). Put \(e=v_2(q)\ge2\). At
\(D_*=M-1/q\) those two orders are regular even if they
exceed \(q\), since \(q\) divides neither. The residue
calculation at order \(q\) is unchanged.

For a same-side pair of offsets \(a<b\), set
\(\delta=b-a\), \(B(t)=\lceil\log_2(t+1)\rceil\), and
\(K_m=2B(m)+B(mb)\). After dividing by the pure moment
of the lower offset, the coefficient of the \(k\)-th
power of the coefficient ratio in order \(m\) is
\[
 c_{m,k}=\binom mk\frac{(ma+\delta k)!}{(ma)!}
                \frac{1}{(m(D_*-a))_{\underline{\delta k}}}.
\]
Each denominator factor has valuation \(v_2(m)-e\).
By the binary digit sum formula for factorials, the
coefficient valuation equals
\[
     (1+e-v_2(m))\,\delta k+\epsilon_{m,k},
     \qquad \max_k\epsilon_{m,k}-\min_k\epsilon_{m,k}\le K_m .
\]
Indeed the error is the binomial valuation plus
\(s_2(ma)-s_2(ma+\delta k)\), with oscillations at most
\(2B(m)\) and \(B(mb)\), respectively. For coefficients
of valuation \(\sigma k+\epsilon_k\) with error
oscillation \(K\), every root has valuation in
\([-\sigma-K,-\sigma+K]\). Outside this interval
the constant or leading term is uniquely least.
This is the slope comparison of the
[row-translation sieve](SFC_ROW_TRANSLATION_DIVISOR_SIEVE.md#even-starts-with-arbitrary-unmarked-rows)
applied to the falling factorial. Here the two
interval centers differ by \(\delta\). Thus when a
unique maximum leaves a same-side pair, the bound
\(\delta>K_{m_e}+K_{m_o}\) proves exclusion at the mark.

For a tied maximum put \(b=M\), take \(a\) to be the
smaller offset on the two-column side and use \(Y,Z\)
for the tied columns. The terminal equation gives
\(Y^q+Z^q=0\). If \(Y=0\), both maximal variables are
zero and a lower equation excludes the third alone.
Otherwise put \(Z=\xi Y\). The lower equations then
have the same binary forms as above except that their
leading coefficients are multiplied by \(1+\xi^m\).
At a place over two, decompose the root of unity
\(\xi^m\) into its odd-order and power-of-two parts.
The power-of-two part has order \(2^{e+1-v_2(m)}\)
because \(\xi^q=-1\). If the odd-order part is
nontrivial, its residue is nontrivial and \(1+\xi^m\)
is a unit. If it is trivial, \(1+\xi^m\) has valuation
\(1/2^{e-v_2(m)}\). This follows by applying the
Eisenstein polynomial \((1+T)^{2^{e-v_2(m)}}+1\)
to \(-\xi^m-1\). In both cases the multiplier is
nonzero and adds at most \(1/2\) to the error
oscillation. The two intervals are therefore
disjoint when \(\delta>K_{m_e}+K_{m_o}+1\).
The coordinate at the smaller offset cannot vanish
on a common direction since the two leading
coefficients are nonzero.

The binary lengths depending on \(b\) can be removed
away from explicit residue classes of \(a\). Fix
\(R\ge1\) and suppose the even-order interval
\(m_ea+1,\ldots,m_ea+m_e\delta\) contains no multiple of
\(2^{R+1}\). The contribution to the valuation of
the factorial ratio in \(c_{m,k}\) from levels
\(r\le R\) is
\[
 \sum_{r=1}^{R}
      \left(\left\lfloor\frac{ma+\delta k}{2^r}\right\rfloor
           -\left\lfloor\frac{ma}{2^r}\right\rfloor\right)
           =(1-2^{-R})\delta k+E_{m,k},
 \qquad |E_{m,k}|<R.
\]
The even order has no higher-level contribution.
For the odd order higher levels add a nondecreasing
function of \(k\). Its center is already the lower
one, so only its upper endpoint is needed, and the
constant-term comparison preserves that upper bound.
The centers differ by \(\delta\), with the required
error bounds \(2B(m)+2R\), before the possible tied
correction. Thus both cases are excluded if
\[
               \delta>4R+2B(m_e)+2B(m_o)+1 .
\]
The omitted residue classes modulo \(2^{R+1}\) have
proportion at most \(m_e\delta/2^{R+1}\).
Indeed, for each \(1\le j\le m_e\delta\), the
congruence \(m_ea+j=0\bmod 2^{R+1}\) has two
solutions when \(j\) is even and none otherwise.
Taking the largest permissible \(R\) gives an
exponentially decreasing bound on this proportion
as \(\delta\) grows at fixed even start. The estimate
does not depend on the absolute size of \(a\).

### The two even marks and high factorial valuations

Use the notation of the preceding subsection. At an
even start \(n\ge4\), put \(Q=n+2\), \(u=n/2\) and
\(v=(n+2)/2\). The same-side pair after deletion, or the
pair on the two-column side at a tie, has offsets
\(a,a+\delta\). The maximum is \(M\) as before. Fix
\(R\ge1\) and assume
\[
  2^R\ge uv\delta,\qquad
  K=2B(Q)+2R+1,\qquad
  \delta>2K+1.
\]
If the decoupled resultant vanishes identically at
the three prescribed orders, then for some integer
\(1\le h\le\delta\),
\[
 R\le v_2(a+h)\le R+Q^2(\delta+2K+3).
 \tag{3}
\]
Thus the obstruction for wide gaps is supported on
explicit residue classes near \(-1,\ldots,-\delta\).
The very deep classes are excluded at the original
consecutive orders as well.

First mark either even order \(q\in\{n,n+2\}\),
and write \(m_e\) and \(m_o\) for the other even and
odd orders. Both are regular at \(D=M-1/q\).
If \(4\mid q\), the preceding subsection applies.
If \(v_2(q)=1\), the falling factorial in the
\(m_e\)-th equation has integral factors at two,
so its reciprocal adds a nonincreasing function of
\(k\) to the coefficient valuations. In the odd
equation each reciprocal denominator factor has
valuation one. Truncating the numerator factorial
valuations at level \(R\) therefore gives centers
\((1-2^{-R})\delta\) and
\((2-2^{-R})\delta\) for the coefficient slopes in
the even and odd equations. We need a lower bound
on the even root valuations and an upper bound on
the odd ones. The nonincreasing even correction
preserves the former by the leading-term comparison.
Higher numerator levels in the odd equation add a
nondecreasing correction, which preserves the latter
by the constant-term comparison.

At a tie the last coefficient is again modified by
\(1+\xi^m\), with \(\xi^q=-1\). For the odd order its
valuation is at most \(1/2\). For the even order
with \(v_2(q)=1\), the two-primary part of \(\xi^{m_e}\)
is trivial. Its multiplier is a unit if that root
of odd order is nontrivial and has valuation one
otherwise. Thus the errors in the required endpoint
comparisons are at most \(K\) for either choice of
mark. With gap greater than \(2K\), a mark can fail
only if the corresponding even numerator interval
contains a multiple of \(2^{R+1}\).

If both marks fail, there are integers
\(1\le j_0\le u\delta\), \(1\le j_1\le v\delta\)
such that \(ua+j_0\) and \(va+j_1\) are divisible
by \(2^R\). The integer \(vj_0-uj_1\) is divisible
by \(2^R\) and has absolute value strictly smaller
than \(uv\delta\). Hence \(j_0=uh,\ j_1=vh\) since
\(\gcd(u,v)=1\). The odd member of \(u,v\) now gives
\(a+h=0\bmod 2^R\).

To exclude high values \(H=v_2(a+h)\), use the
mark \(q\) divisible by four and put \(e=v_2(q)\).
Here \(2^{R+1}>Q\delta\), so for each remaining
order \(m\) only the factor \(ma+mh\) can contribute
above level \(R\). With \(j_m=\lceil mh/\delta\rceil\)
and \(L_m=m-j_m+1\), the coefficient valuations are
\[
 \sigma_m k+\epsilon_{m,k}
       +(H+v_2(m)-R)\mathbf 1_{k\ge j_m},
 \qquad
 \sigma_m=(1-2^{-R}+e-v_2(m))\delta .
\]
The oscillation of \(\epsilon_{m,k}\) is at most \(K\),
including a last-coefficient correction at a tie.
For \(H\) larger than the upper bound in (3),
the segment from \(j_m-1\) to \(m\) is the last
Newton edge. Indeed its rise beyond \(\sigma_m L_m\)
is \(H+v_2(m)-R+\epsilon_{m,m}-\epsilon_{m,j_m-1}\).
The segment lies strictly below its interior points
because the extra jump divided by \(L_m\) exceeds
\(2K\). The same inequality places its backward
extension below all earlier points. All earlier
edge slopes belong to \([\sigma_m-K,\sigma_m+K]\).
These intervals for the two orders are disjoint.

The slope of the last edge in each order exceeds
both earlier slope intervals. If \(L_{m_e}=L_{m_o}\),
their last slopes differ by at least
\(\delta-(1+2K)/L_{m_e}>0\). If the lengths are
different, the contributions of \(H-R\) to their
slopes differ by at least \((H-R)/Q^2\). All remaining
contributions to the slope difference have total
absolute value at most \(\delta+2K+1\).
The two Newton polygons have disjoint slopes in
either case, proving (3).

### A narrow-gap consequence at the upper endpoint

Even when \(\delta\) is small, the class \(a\) very
close to \(-\delta\) at two can be excluded for every
even start \(n\ge2\). In fact write \(s=s_2(\delta-1)\)
and \(H=v_2(a+\delta)\). It is enough that
\[
                   2^H> (n+2)\delta,\qquad
                   H>\delta+1+2s .
 \tag{4}
\]
Use the mark \(q\) divisible by four. In the factorial
ratio for order \(m\), the factors before \(ma+m\delta\)
have exactly the valuations of the descending integers
from \(m\delta-1\) to \(1\), since
\(ma+j=m(a+\delta)-(m\delta-j)\). The last factor
has valuation \(H+v_2(m)\).

For an order \(m\), set aside its last coefficient.
With \(l=m-1-k\ge1\), the rise from degree \(k\) to
degree \(m-1\) is
\[
 (e-v_2(m))\delta l+v_2(m)-v_2\binom mk
       +\delta l+s-s_2(\delta(l+1)-1).
\]
For the even order, all slopes of this truncated
polygon are at most \(e\delta+s\), by comparison
with its leading term. If \(\delta\) is a power of
two, the bound improves to \(e\delta\), since
\(s_2(\delta(l+1)-1)=s+s_2(l)\ge s+1\).
In the truncated odd order the first \(\delta k\)
factors contain at least \(\delta k/2\) even
integers. For even \(\delta\) this follows by pairing.
For odd \(\delta\), \(a\) is odd and the first factor
is even. Comparing with the constant term makes
every odd truncated slope at least \(e\delta+\delta/2\).
These bounds are disjoint: \(s_2(\delta-1)<\delta/2\)
for all positive \(\delta\) except \(2,4\), which
are powers of two. To see the digit inequality, for
odd \(\delta\) use \(s_2(2k)=s_2(k)\le k\).
For even \(\delta=2k+2\) use
\(s_2(2k+1)=1+s_2(k)<k+1\) when \(k\ge2\).

All truncated slopes in either order are at most
\(e\delta+\delta+s\) by the same degree-\(m-1\)
comparison. The rise to the last coefficient is
\[
                  (e-v_2(m))\delta+H+v_2((\delta-1)!)+\tau_m ,
\]
where \(\tau_m=0\) in the unique-maximum case and
\(0\le\tau_m\le1/2\) at a tie. The \(v_2(m)\) from
the last numerator factor cancels the binomial
valuation at degree \(m-1\). Condition (4) makes this
rise larger than both truncated slope bounds, so
the final edges have width one. Their slopes
differ by at least \(\delta-1/2>0\).
No slope is common. This includes the gap-one case
with \(v_2(a+1)>\max\{2,\log_2(n+2)\}\).

### Exact certificates with the offsets left free

For start \(3\) and \(1\le\delta\le5\), and for starts
\(4,5,6,7\) with \(\delta=1\), the consecutive degree sieve
holds for every integer \(a\) and \(M\) in both the
same-side survivor and tied-maximum cases.

Write \(C=M-a\). For each \(q\in\{n,n+1,n+2\}\),
mark \(D=M-1/q\). The other two orders are regular
because \(n\ge3\) and \(q\) divides neither. The
coefficient of \(z^k\) in the binary form after
clearing the denominator is
\[
  \binom mk
  \prod_{j=1}^{\delta k}(ma+j)
  \prod_{j=\delta k}^{m\delta-1}(m(C-1/q)-j).
  \tag{5}
\]
At a tie set \(C=\delta\) and multiply the leading
coefficient by \(1+t^m\), where \(t^q=-1\).
A common direction at the mark forces the binary
resultant of (5) to vanish. At a tie it forces the
resultant of that eliminant with \(t^q+1\) to vanish.
These implications include the binary chart at infinity,
since we specialize the Sylvester determinants with
their original binary sizes. The case where both
maximal variables vanish is excluded by a lower
pure moment.

For gap one, before denominator clearing the binary
polynomial is a Gauss hypergeometric polynomial with
Jacobi parameter \(-m(C-1/q)-1<-1\). The hypotheses in
[Arvesú, Driver and Littlejohn, Theorem 1](https://arxiv.org/pdf/2009.10196)
require both a positive-weight parameter range and
coprimality of the pair.

The checker removes factors from the three eliminants
by exact gcd division with
\[
\begin{split}
 & (q a+qC-1)(q(C-a-\delta)-1)
       \prod_{k=1}^{\delta}(a+k)(qC-1-q(k-1))
          &&\text{for the unique maximum},\\
 & (q a+q\delta-1)\prod_{k=1}^{q\delta}(a+k)
          &&\text{for the tied maximum}.
\end{split}
\]
All displayed factors are nonzero at the admissible
integer parameters. The quotients are cleared to
integer polynomials. If all three marks fail at one
integer profile, then all three stripped eliminants
vanish at its simultaneous residues. This is impossible
if their reductions generate the unit ideal over a
prime field, or if a complete search modulo a specified
integer has no surviving residues.

The [certificate checker](verify_sfc_consecutive_axis_marks.py)
uses the unit ideal over \(\mathbb F_{65537}\) except
for the unique maximum at starts \(4\) and \(6\),
gap one. The complete dyadic searches have 128 and
524288 survivors at levels 18 and 20, respectively,
and zero at the next level in each case. Each step
lifts both coordinates of all surviving residues and
reevaluates the eliminants. Internal replay used
Python 3.12.3 and Singular 4.3.2.

~~~sh
python3 research/factorial/image-rigidity/verify_sfc_consecutive_axis_marks.py 3 1 2 3 4 5
python3 research/factorial/image-rigidity/verify_sfc_consecutive_axis_marks.py 4 1 --mode tied
python3 research/factorial/image-rigidity/verify_sfc_consecutive_axis_marks.py 4 1 --mode unique --dyadic-depth 19
python3 research/factorial/image-rigidity/verify_sfc_consecutive_axis_marks.py 5 1
python3 research/factorial/image-rigidity/verify_sfc_consecutive_axis_marks.py 6 1 --mode tied
python3 research/factorial/image-rigidity/verify_sfc_consecutive_axis_marks.py 6 1 --mode unique --dyadic-depth 22
python3 research/factorial/image-rigidity/verify_sfc_consecutive_axis_marks.py 7 1
~~~

## 2. Returning to ordinary degrees

Let \(Q\) be the largest selected order, which may
differ from the marked order. At \(D\ge2QM\), any
mixed-side count vector gives
two factorial arguments each at least \(D/2\). With
\(l=\lfloor D/2\rfloor\), its coefficient in the moment
divided by \((mD)!\) is at most \(Q!2^{-l}\). Indeed it
is the multinomial coefficient divided by a binomial
coefficient whose smaller index is at least \(l\);
\(\binom Nk\ge2^k\) for \(0\le k\le N/2\), by the
product expression for the binomial coefficient.
The actual and decoupled coefficients are bounded by
\(Q!\) at the ordinary large degrees. Consequently
the actual resultant differs from the decoupled one
by at most \(C2^{-l}\) for an effectively computable
constant \(C\): expand the fixed resultant polynomial
and telescope each monomial of it. If \(M=0\) for
opposite-side pairs, also require \(D\ge1\).

The decoupled resultant is a nonzero rational function
of \(D\). Beyond a computable threshold its absolute
value is bounded below by \(cD^{-b}\) for a computable
\(c>0\) and \(b\ge0\), by bounding the lower coefficients
of its numerator and using an upper bound for the
denominator. Comparison with the exponential bound gives
an effective ordinary cutoff. For example, one can
bound the ratio of successive values of \(l^b2^{-l}\)
above by a fixed number less than one for all sufficiently
large \(l\), and then search for a strict inequality.
Test the exact factorial forms at the finitely many
remaining admissible degrees. On one side the rational
resultant is already the ordinary normalized resultant
and its numerator suffices.

The projective common-zero tests use rational
coefficients. If there is a common direction over any
characteristic-zero field, a dehomogenized coefficient
chart has a zero over \(\overline{\mathbb Q}\) by the
Nullstellensatz, and therefore over \(\mathbb C\). A
failure in the computed list refers to a geometric
common zero. An individual field can have fewer failures.

### More than two variables

In \(d\) variables, assign each column a growing
coordinate and fix the other \(d-1\) exponents. Its
total fixed offset is the sum of those exponents.
The growing coordinate has exponent \(D\) minus
the total fixed offset. Suppose an offset occurs at
most twice, and if twice its columns grow in
different coordinates. The decoupled forms retain
only count vectors whose columns have the same
growing coordinate. Their numerators are now the
products of off-coordinate factorials and the
denominator remains \((mD)_{\underline U}\), with
\(U\) the count-weighted total fixed offset.

Section 1 applies up to the tied maximum. There
the terminal mark now gives \(\xi^q=-c\), with
\(c>0\) rational because the two residues are
positive multiples of their factorial numerators.
Their ratio need not be one. If \(q\) is a power
of two, the ratio \(\xi/\bar\xi\) is a primitive
\(q\)-th root of unity: the argument of \(\xi\)
is an odd multiple of \(\pi/q\). Complex
conjugation preserves the rational degree of
\(\xi\). If that degree is \(h\), the compositum
of \(\mathbb Q(\xi)\) and \(\mathbb Q(\bar\xi)\)
has degree at most \(h^2\), so \(h^2\ge q/2\).
Choosing sufficiently large \(q\) still exceeds the
fixed resultant degree bound. The subprofile at
\([\eta:\xi]=[0:1]\) and the two-column base use the
same proofs as before. For the ordinary comparison,
a mixed count vector still has at least two factorial
arguments at least \(D/2\). Its multinomial
denominator is at least the binomial coefficient
formed by choosing either of these arguments out
of \(mD\). The same exponential estimate proves
the degree cutoff in any dimension.

## 3. Translations in one coordinate on distinct offsets

Fix arbitrary nonnegative exponents \((u_i,v_i)\) with
distinct \(v_i\), and add the integer \(T\ge0\) to every
\(v_i\). The \(u_i\) are unrestricted. The ordinary
normalized forms are polynomials in \(T\) with
coefficients
\[
                 \binom m\alpha
                 \Bigl(\sum_i\alpha_i u_i\Bigr)!
                 \prod_{j=1}^{\sum_i\alpha_i v_i}(mT+j).
\]
Subtract the least fixed \(v_i\) from every \(v_i\) and
add it to \(T\). For terminal order \(q\), put \(T=-1/q\)
in this rebased family. The terminal form retains only
the pure term at the unique least \(v_i=0\), because
every other term contains \(qT+1\). The restricted
lower forms are polynomial in the rebased \(T\), and
their resultants are nonzero by induction. Rebasing a
subprofile by its own minimum only shifts its parameter
and multiplies each form by a nonzero order-dependent
polynomial, so generic nonvanishing is unaffected.
For two columns \(q=n+1\) works directly: the surviving
pure coefficient in order \(n\) is nonzero since \(n<q\).
At each subsequent step sufficiently large powers of two
avoid the finitely many bad marks on the subprofiles.
Thus on every fixed line of this kind there are finitely
many ordinary failures for the selected \(s\) orders,
and the same selection works for any prescribed finite
list of lines. The exact forms are polynomial, so the
resultant gives an effective degree bound or exact list.

## 4. Both independent heights

The
[Gaussian normalization for independent heights](SFC_BIVARIATE_TRANSLATION_GAUSSIAN_JETS.md#1-exact-graded-normalization)
reduces escape with both heights large to a common
direction of
\[
 G_m(a)=\mathbb E\left[\left(\sum_{g,k}a_{g,k}m^{d_g}
                   Z_m^k\right)^m\right],
 \qquad Z_m\sim \mathcal N(0,1/m).
\]
Here \(d_g\) is a total-degree class, and \(0\le k<r_g\)
for its multiplicity \(r_g\). Choose the adjacent pair
\(n,n+1\) and \(s-2\) additional even orders with
pairwise distinct two-adic valuations also distinct from
that of the even member of \(n,n+1\). For \(s\ge3\), take
the first additional even order to be \(n+2\) if \(n\) is
even and \(n+3\) if \(n\) is odd. The two even orders are
two apart and have distinct valuations. Any further even
orders can be higher powers of two with new valuations.
This gives one odd order and \(s-1\) even orders with
distinct positive two-adic valuations. For a homogeneous support
the \(s\) real numbers \(\lambda_{g,k}=d_g-k/2\) are
distinct. In the inhomogeneous assertion their
distinctness is the Gaussian weight condition.

Suppose there is a common direction. Take algebraic
coefficients using the Nullstellensatz and a place over
two. For an even order of valuation \(e=v_2(m)>0\), the
nonzero terms in the moment expansion have valuations
\[
      v_2\binom m\alpha+
      \sum_{g,k}\alpha_{g,k}\bigl(v_2(a_{g,k})+e\lambda_{g,k}\bigr).
\]
The double factorial for an even Gaussian power is odd,
and every pure term in an even order occurs. Thus if
one of the affine weights is uniquely least, its pure
term is uniquely least in the expansion and gives
nonvanishing. A zero forces a tie on the lower envelope.
The slope drops strictly at each tie. With at most \(s\)
distinct slopes, \(s-1\) distinct tie points exhaust
all transitions and require all \(s\) coefficients to
be nonzero. The largest slope is uniquely least before
the first tie, in particular at \(e=0\). The largest
slope occurs at \(k=0\) in the highest total-degree
class, so its pure term also occurs at the odd order.
It is then uniquely least, a contradiction. This extends
the [three-column valuation argument](SFC_COLLINEAR_TERNARY_EVEN_START_AND_DYADIC_REDUCTION.md#even-orders-and-distinct-valuations).

By the exact Gaussian normalization and its compactness
bound, these orders leave a computable finite list of
horizontal and vertical height lines. On a homogeneous
line the varying degree is an axis profile with all
columns on one side. On an inhomogeneous line in the
stated class, its varying row has distinct offsets and
Section 3 applies. Choose \(s-2\) additional powers of
two simultaneously for the finitely many lines, using
\(n,n+1\) as the first two orders. Sections 1–3 give
finitely many computable candidates on those lines.
Testing the full selected set of at most \(2s-2\)
moments gives the exact exceptional height list.

At start two, use \(2,3,4\) and \(s-3\) higher powers
of two for the Gaussian orders. These still have the
required distinct even valuations. On the remaining
homogeneous lines, the
[orders \(2,3,4\) axis theorem](SFC_HOMOGENEOUS_AXIS_WINDOWS_2_3_Q.md)
gives nonzero rational resultants on all three-column
subprofiles. Begin the induction in Section 1 with that
base to add at most \(s-3\) orders. This gives the bound
\(2s-3\).

The bound uses two choices of powers of two. The Gaussian
orders first identify the lines, and a single later
selection works on all those lines. To get a consecutive
theorem, one would need an exclusion at the prescribed
terminal order. For trinomials it is enough on a fixed
axis to settle both the same-side adjacent pair at
\(D=M-1/(n+2)\) after a unique maximum is deleted and
the tied case when the cyclotomic degree is small.
Condition (2) handles separated maxima at every start.
For even starts, (3) bounds the depth of the remaining
wide-pair obstruction, and (4) covers the deep
upper-endpoint cases at every gap.
Moving the pole to \(M-h/q\) can retain mixed terms:
a count vector with deficit \(qM-U<h\) also has the
vanishing denominator factor when \(1\le h<q\).
For independent heights there is also a Gaussian
obligation when the prescribed consecutive orders
supply too few distinct even valuations.
