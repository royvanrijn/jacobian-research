# Endpoint normalization and integer anchors for factorial moments

Let
\[
       f=a_0X^\sigma\left(1+\sum_{\nu\ne\sigma}c_\nu X^{\nu-\sigma}\right)
\]
be a finite polynomial over a finite extension of \(\mathbb Q_p\),
with \(a_0\ne0\), nonzero displayed coefficients and
nonnegative ordinary exponent columns. The expression in
parentheses can be a Laurent polynomial. Normalize the
ordinary moments at the distinguished column:
\[
       R_n=\frac{\Lambda(f^n)}{a_0^n\prod_i(n\sigma_i)!}.
\]
Assume a strict effective slope preference,
\[
 \tau=\min_{\nu\ne\sigma}
    \left(v(c_\nu)+\frac{|\nu|-|\sigma|}{p-1}\right)>0,
         \qquad v(p)=1.                                     \tag{1}
\]
This normalization applies when all columns are
nonsingletons. The result below gives integer-anchor limits
for arbitrary gap orientations. Under a condition on the
decreasing exponent rows it gives locally analytic
interpolation on the full order domain.

For a vector \(m=(m_\nu)\) of nonnegative counts denote
\[
 h=\sum_\nu m_\nu,\qquad D_i=\sum_\nu m_\nu(\nu_i-\sigma_i).
\]
Write \(z^{\underline h}\) for the falling factorial.
Associate to the vector \(m\) the rational function
\[
 A_m(z)=
 \frac{z^{\underline h}\prod_\nu c_\nu^{m_\nu}}{\prod_\nu m_\nu!}
 \frac{\prod_{D_i>0}\prod_{j=1}^{D_i}(z\sigma_i+j)}
      {\prod_{D_i<0}\prod_{j=0}^{-D_i-1}(z\sigma_i-j)}.        \tag{2}
\]
At each nonnegative integral \(n\) the moments have the exact
finite expansion \(R_n=\sum_{h\le n}A_m(n)\). The factors
in its denominator are positive integers because
\(-D_i\le h\sigma_i\le n\sigma_i\).

## Integer anchors for arbitrary gap orientations

Let \(d_-\) be the largest number of negative coordinates of
\(D\) among nonnegative integral count vectors. Fix an integer
anchor \(a\). Put \(\lambda=0\) when \(a<0\), and
\(\lambda=\max(0,d_--1)\) when \(a\ge0\).
Then each
\[
                         B_m(z)=(z-a)^\lambda A_m(z)
\]
is regular at \(a\), and the series \(V_a=\sum_m B_m(a)\)
converges. For any fixed real \(b\ge1\), there is a constant
\(C\) such that, for all sufficiently large integers \(s\),
\[
 (n-a)^\lambda R_n\equiv V_a\pmod{p^{\,s-C}}
 \quad\left(\begin{array}{l}
 n\ge1,\ n\ne a,\ n\equiv a\bmod p^s,\\
 n\le p^{bs}
 \end{array}\right).                                       \tag{3}
\]
The congruence is in the sense of the normalized valuation.
Thus \(V_a\ne0\) excludes exact moment vanishing in these
order classes for large \(s\). A finite collection of
successive integer anchors gives the same conclusion for
windows when at least one anchor value in the collection
is nonzero. The upper bound in (3) can be enlarged by a
fixed factor to account for the window offsets.

Here is the removal of poles. If \(a<0\), an index with
\(D_i<0\) has \(\sigma_i>0\), so \(a\sigma_i-j\ne0\).
If \(a\ge0\), each decreasing coordinate supplies at most
one zero denominator factor, of the form
\(\sigma_i(z-a)\). Such a zero requires \(-D_i>a\sigma_i\),
hence \(h>a\). The numerator then supplies the factor
\(z-a\) from \(z^{\underline h}\). After that cancellation,
\(\lambda\) covers all remaining powers of \(z-a\).

For estimates, a nonzero interval of \(H\) consecutive
integers of absolute values at most \(M\) has product
valuation at least \(v_p(H!)\) and at most
\(H/(p-1)+\lfloor\log_p M\rfloor\), by counting multiples of
each \(p^j\). The factorial lower bound
\[
             v_p(H!)\ge H/(p-1)-1-\log_p(1+H)
\]
follows from Legendre's formula. For actual positive \(n\)
and \(h\le n\), the multinomial factor in (2) is integral.
The positive rows therefore give their lengths divided by
\(p-1\), up to a logarithmic loss in \(1+h\). Each negative
row loses its length divided by \(p-1\), plus at most a
logarithmic term in \(1+n\). In particular there is a constant
\(C_0\) depending on the support and prime such that
\[
 v(A_m(n))\ge
 \tau h-C_0\bigl(1+\log_p(1+h)\bigr)-d\log_p(1+n),           \tag{4}
\]
where \(d\) is the number of variables.

For the anchor estimate, choose an integer
\(J_h=\lceil\log_p(C_1(1+h))\rceil+1\) with \(C_1\)
large enough that \(p^{J_h}\) exceeds every factor length
and the absolute value of every nonzero affine factor at
\(a\). This is possible because the support and anchor
are fixed. After the pole removal, consider
\(v(z-a)\ge s\ge J_h\). Divide each affine factor by
\(p^{\min(v(g(a)),J_h)}\), with the minimum \(J_h\)
if the factor vanishes at \(a\). The resulting factors
are integral on that domain. The nonzero denominator
factors are units, and each factor changes by an amount
of valuation at least \(s-J_h\) between \(z\) and \(a\).
This property is preserved by products and unit inverses.
The residual nonnegative power of \(z-a\) is integral and
has the same comparison bound.

The stripped valuations in a numerator interval of length
\(H\) sum to at least \(v_p(H!)\), because at each level up
to \(J_h\) the original interval has at least
\(\lfloor H/p^j\rfloor\) multiples. If we removed one zero
from the falling factorial, this lower bound loses at most
\(J_h\). Stripped valuations in each denominator interval
sum to at most its length divided by \(p-1\) plus \(2J_h\),
including the possible replacement of a vanishing affine
factor by the constant \(\sigma_i\). Using the scalar
\(\prod m_\nu!^{-1}\) and (1), the net stripped valuation
is consequently at least \(\tau h-C_2J_h\) for a fixed
\(C_2\). We have proved
\[
 \begin{aligned}
       v(B_m(a))&\ge \tau h-C_2J_h,\\
 v(B_m(z)-B_m(a))&\ge
       s+\tau h-(C_2+1)J_h\qquad(v(z-a)\ge s\ge J_h).
 \end{aligned}                                               \tag{5}
\]
The first bound follows at \(a\) directly from the stripped
factors. It proves convergence of the anchor series.

Choose \(H\) sufficiently large that the bound in (4) is
at least \(s-C\) whenever \(h\ge Hs\) and \(n\le p^{bs}\).
For example, absorb the logarithm of \(1+h\) into
\((\tau/2)h\) plus a constant, and take
\(H>2(db+1)/\tau\). Multiplication by \((n-a)^\lambda\)
does not reduce the valuation. The first bound in (5)
also disposes of the anchor tail past \(Hs\).
For \(h<Hs\), eventually \(J_h\le s\), and
\(\tau h-(C_2+1)J_h\) is bounded below independently of \(h\).
The second bound in (5) compares all head terms
at precision \(s-C\). The head is present in the actual
moment because \(n\ne a\) and \(n\equiv a\bmod p^s\)
give \(n\ge p^s-|a|\) for large \(s\).
Taking the two tails and the head together proves (3).

## A finite negative-anchor witness

Suppose \(\sigma_i>0\) in every variable and the distinguished
column has strictly least total degree. Let
\(\delta=\min_{\nu\ne\sigma}(|\nu|-|\sigma|)>0\).
At \(a=-r\), \(r\ge1\), any count with \(D_i\ge r\sigma_i\)
has a zero positive-row factor in (2). Every potentially
nonzero term therefore satisfies
\[
       h\delta\le\sum_i D_i\le r|\sigma|-d.                  \tag{6}
\]
Thus \(V_{-r}\) is a finite rational polynomial in the
coefficient ratios, with constant coefficient one.
In particular if \(\delta>|\sigma|-d\), then \(V_{-1}=1\)
for every choice of coefficients. At any prime satisfying
(1), (3) gives the nonzero-moment witness at sufficiently
deep orders \(n\equiv-1\bmod p^s\) in the stated range.
For example, the fully shared columns \((2,2,2)\) and
\((1,1,12)\) have two decreasing rows and satisfy this
strict gap bound.
Good primes over a number field, where the coefficient
ratios are units, satisfy (1). The exponent gaps in
individual rows may have either sign.

## Full analytic interpolation with at most one decreasing row per count

Suppose every count vector has at most one negative \(D_i\),
and every row that can decrease has
\(\sigma_i=q=p^e\) for some \(e\ge0\). The row can depend
on the count. There is then a locally analytic function
\(F\) on \(\mathbb Z_p\) and constants \(c>0,C'\)
such that
\[
                     v(R_n-F(n))\ge c n-C'\qquad(n\ge1).    \tag{7}
\]
If no row can decrease then \(R_n=F(n)\) for every
nonnegative integer \(n\), and \(F(0)=1\).
Whenever \(a<0\), the analytic function satisfies
\(F(a)=V_a\). In particular the finite witness
\(V_{-1}=1\) under (6) certifies all sufficiently large
orders in some fixed class \(n\equiv-1\bmod p^s\)
without the upper bound in (3).
The same finite-precision window criterion applies to
\(F(z),\ldots,F(z+w-1)\) on \(\mathbb Z_p\): absence of a
common zero certifies all sufficiently late \(w\)-windows.

For the proof, if the negative row has length \(L=-D_i>0\)
put \(r=\lceil L/q\rceil\le h\). The multiples of \(q\)
in its falling product give the factor \(q^r z^{\underline r}\).
Cancel it from (2), writing the remaining quotient as
\[
 \frac{z^{\underline h}}{(qz)^{\underline L}}
    =\frac{(z-r)^{\underline{h-r}}}{Q}
       \prod_{\substack{1\le j<L\\q\nmid j}}
                               (1-qz/j)^{-1},
 \qquad Q=q^r\!\!\prod_{\substack{1\le j<L\\q\nmid j}}(-j).
                                                               \tag{8}
\]
The scalar satisfies
\[
 v_p(Q)=\sum_{\ell=1}^e \left\lceil L/p^\ell\right\rceil
        \le \frac{L(1-1/q)}{p-1}+e .
\]
The remaining rational factors in (8) are analytic units
on the closed unit disk, since \(v_p(q/j)\ge1\).
When there is no negative row take \(L=r=0,Q=1\).
In all cases the falling factorial in (8) is divisible by
\((h-r)!\) on \(\mathbb Z_p\). The positive row factors
are divisible by \(D_i!\) there. Counting valuations gives
\[
  v(A_m(z))\ge
  \sum_\nu m_\nu v(c_\nu)+\frac{\sum_{D_i>0}D_i-L}{p-1}
     -C_3\bigl(1+\log_p(1+h)\bigr)\qquad(z\in\mathbb Z_p),    \tag{9}
\]
because \(r\le L/q+1\), and the linear valuation costs of
\(r\) and \(Q\) together are at most \(L/(p-1)\) plus a
fixed constant. Here and below \(A_m\) denotes its
pole-cancelled expression.

To justify analytic convergence, on a disk \(z=a+p^s y\)
with \(a\) an integer and \(|y|\le1\), each product of
\(H\) consecutive affine factors in the numerator has
Gauss valuation at least
\(\sum_{\ell\le s}\lfloor H/p^\ell\rfloor\).
Indeed each affine coefficient of \(y\) is divisible by
\(p^s\), and the constant factors form a consecutive
integer interval. Compared with \(v_p(H!)\) the loss is
at most \(Hp^{-s}/(p-1)\).
The total numerator lengths are bounded by a support
constant times \(h\). Choose fixed \(s\) so that their total
linear loss is below \((\tau/2)h\). The unit inverses in
(8) introduce no Gauss loss. Thus
\(\sum_m A_m(z)\) converges on every such disk as a
restricted power series. Define \(F\) by these sums.
By (9) its terms with \(h>n\) have valuation at least
\(cn-C'\) uniformly on \(\mathbb Z_p\), for some \(c>0\),
proving (7). If no row decreases, \(z^{\underline h}\)
vanishes at \(z=n\) for all \(h>n\), which makes the
interpolation exact. For negative integers there was no
pole in (2), proving \(F(a)=V_a\).

Finally, the functions on \(\mathbb Z_p\) are continuous
and the domain is compact. The absence of a common zero
gives a positive lower bound on the maximum absolute
value in the window vector. The bound (7) then transfers
it to sufficiently late windows. The finite tail bounds
and the finite head continuity estimates allow this
criterion to be checked at a fixed precision.

The coordinatewise nondecreasing case also has a direct
unit test. Suppose all relative coefficients are integral.
For \(n\) divisible by \(p\), any count with a component
\(m_\nu\ge p\) has a positive rising interval of length
at least \(p\). That interval product is divisible by
\(p\), and the multinomial factor is integral. If every \(m_\nu<p\)
and \(h>0\), the falling factorial in (2) supplies \(n\)
and the count factorials are units. Thus \(R_n\equiv1\bmod p\)
at all positive multiples of \(p\). At \(p=2\) this gives
an all-order two-window even with arbitrary support size.
The underlying division by the endpoint factorial is
the one used at prime orders in
[van den Essen, Wright and Zhao, Proposition 4.7](https://arxiv.org/pdf/1008.3962).
The conclusions here include every positive multiple of
the fixed prime, exact interpolation in the ordered case,
and the integer-anchor comparison for mixed gaps.
