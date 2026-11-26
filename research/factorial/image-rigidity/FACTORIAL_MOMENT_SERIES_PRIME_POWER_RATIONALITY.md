# Prime-power rationality of factorial moment series

Let \(\nu\) range over a finite set of distinct nonzero exponent
vectors in \(\mathbb Z_{\ge0}^{d}\), and put
\[
 f=c+\sum_\nu \frac{u_\nu}{\nu!}X^\nu,\qquad
 \nu!=\prod_i\nu_i!,\qquad
 P_n=\Lambda_d(f^n)/n!,\qquad E(T)=\sum_{n\ge0}P_nT^n.
\]
Work in a finite extension of \(\mathbb Q_p\), with valuation \(v(p)=1\)
and ring of integers \(\mathcal O\). Assume \(u_\nu\in\mathcal O\).

**Theorem.** If \(c=0\) or \(v(c)>1/(p-1)\), then \(E(T)\) has integral
coefficients and its reduction modulo \(p^a\mathcal O\) is rational for
every \(a\ge1\), with a denominator of constant coefficient one.
For \(c=0\) the assertion holds universally over
\((\mathbb Z/p^a)[(u_\nu)_\nu]\).
A finite matrix given below produces the entire reduced sequence
in the moment index \(n\).

For fixed local coefficients, the matrix gives an exact two-sided
decision as to whether this modulus certifies a nonzero moment in
every positive-index window of a given length. The decision concerns
this one modulus. A surviving zero window modulo \(p^a\) need not
consist of exact zero moments.

The main distinction is between the moment order \(n\) and the
exponent vector of each term of \(f^n\). A partition recurrence in
the exponent vector has periodic coefficients modulo \(p^a\).
Summing the recurrence by residue classes then gives a constant
matrix in the moment order.

## The integral partition recurrence

First let \(c=0\). For an exponent vector \(N\ge0\) define
\[
       W_N(T)=N!\,[X^N]\exp(Tf),\qquad W_0=1.
\]
If a component of \(N\) is negative set \(W_N=0\).
For \(N\ne0\) choose the least coordinate \(i\) with \(N_i>0\).
Differentiation in \(X_i\), followed by exact coefficient comparison,
gives
\[
 W_N=T\sum_{\nu_i>0} u_\nu
           \binom{N_i-1}{\nu_i-1}\prod_{j\ne i}\binom{N_j}{\nu_j}
                  W_{N-\nu}.                                      \tag{1}
\]
Indeed the derivative of \(X^\nu/\nu!\) in \(X_i\) is
\(X^{\nu-e_i}/(\nu-e_i)!\) if \(\nu_i>0\). The remaining factors
are the coefficient factorial of \(N-e_i\) divided by
\((\nu-e_i)!(N-\nu)!\).
This is also the count obtained by specifying the block containing
one distinguished label of color \(i\). Each variable is a color
and a monomial \(\nu\) is a block with \(\nu_j\) labels of color \(j\).
Induction on \(\sum N_i\) proves that \(W_N\) has coefficients
in \(\mathbb Z[(u_\nu)_\nu]\). We have \(E=\sum_N W_N\);
the sum is formal, since an order \(n\) uses only exponent vectors
with \(\sum N_i\le n\max_\nu\sum\nu_i\).

In one color, (1) is the usual restricted set partition recurrence
with each block weighted by \(T\) and its coefficient. For comparison,
[Rocks, Theorem 3.10](https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf)
uses the factorial map modulo \(p\) for sums of \(p\)-power monomials
on separate variables. The recurrence above keeps an arbitrary
collection of multicolored blocks and works at every fixed precision.

## The finite matrix at precision \(p^a\)

Choose \(M\ge1\) at least the largest exponent entry, and take
\(K=p^{a+\lfloor\log_p M\rfloor}\). For \(0\le b\le M\),
\(\binom{t}{b}\bmod p^a\) is periodic in \(t\ge0\) with period \(K\).
Indeed for \(1\le j\le b<K\),
\[
 v_p\binom Kj\ge v_p(K)-v_p(j)\ge a
\]
by \(j\binom Kj=K\binom{K-1}{j-1}\), and the claim follows from
Vandermonde's identity.

A state records \(N_j\bmod K\) and whether \(N_j=0\), for each
coordinate \(j\). There are at most \((K+1)^d\) states. The distinction
between zero and a positive multiple of \(K\) determines the
coordinate used in (1). For each state \(b\) write \(S_b(T)\) for
the sum of the \(W_N\) in that state. Adding \(\nu\) to a vector in
state \(b\) determines a unique target state \(b'\): add its residues
and take the union of the positive-coordinate sets. Let \(i\) be
the first positive coordinate of \(b'\). The transition from \(b\)
to \(b'\) for \(\nu\) is allowed only if \(\nu_i>0\), and then has
weight
\[
      u_\nu\binom{N_i-1}{\nu_i-1}\prod_{j\ne i}\binom{N_j}{\nu_j}
                          \pmod{p^a},                           \tag{2}
\]
where the upper entries can be replaced by their representatives
in \([0,K-1]\). Here \(N\) is the target. Set \(C_{b',b}\) to the
sum of these weights.

Summing (1) by states gives
\[
 S=e_0+TCS,\qquad
 E=\boldsymbol 1^{\mathsf T}(I-TC)^{-1}e_0
                                      \pmod{p^a}.              \tag{3}
\]
Every term on the right of (1) comes from a nonnegative predecessor,
and each allowed transition enumerates exactly those predecessors.
The unique vector of total exponent zero supplies \(e_0\).
Thus (3) has no additional boundary terms. The inverse is formal
and the adjugate gives a denominator \(\det(I-TC)\) with constant
coefficient one, even over the ring with nilpotents
\(\mathcal O/(p^a)\).

After specializing the coefficients, the vector sequence
\(C^ne_0\) takes values in a finite set, since \(\mathcal O/(p^a)\)
is finite. Iterate to its first repetition to obtain the preperiod
and the cycle. For a proposed length \(s\), test the scalar sequence
\(\boldsymbol1^{\mathsf T}C^ne_0\) for zero blocks through the
preperiod and one period of starting indices, using periodic
continuation at the end. This is an exact finite procedure and
includes all positive starts. The [linear precision theorem](FACTORIAL_LINEAR_PRECISION_ANALYTIC_INTERPOLATION.md)
replaces the growing residue matrix by a fixed differential
matrix and obtains joint local analytic interpolation in order
and coefficients.

For \(c\ne0\), multiply (3) by \(\exp(cT)\). If \(v(c)>1/(p-1)\),
the bound
\[
 v(c^k/k!)\ge k\bigl(v(c)-1/(p-1)\bigr)
\]
makes that factor a polynomial at every fixed precision.
The larger class with equality need not be rational at all
precisions. For example, at \(p=2,\ c=2,\ f=c\), the coefficient
\(2^n/n!\) has valuation \(s_2(n)\), where \(s_p(n)\) is digit sum.
Modulo \(4\), its positive-order support is exactly the powers
of two. This infinite support is not eventually periodic,
so the series is not rational with a denominator of constant
coefficient one over \(\mathbb Z/4\).

## Degree bounds in characteristic \(p\)

At residue precision, the finite matrix has a more explicit
formula that avoids enumerating its states. Call \(\nu\) a singleton
if it is \(p^b e_i\), where \(e_i\) is a coordinate vector.
Write \(d_{\nu,i,t}\) for digit \(t\) of \(\nu_i\), and
\(w_\nu=\sum_{i,t}d_{\nu,i,t}\). Thus the singletons have
\(w_\nu=1\). Write a bar for reduction to the residue field.
For each coordinate and position \(t\ge0\), put
\[
 A_{i,t}(T)=1-\sum_{\nu=p^b e_i,\ b\le t}
                         (\bar u_\nu T)^{p^{t-b}} .
\]
Let \(h_i\) be the largest position with a nonzero exponent digit
in coordinate \(i\), with \(h_i=-1\) if there is none. Then, for
zero constant term, \(E\) reduces to the following rational function:
\[
 \overline E(T)=
 \frac{\displaystyle
   \sum_{m}
     \prod_{\nu:\ w_\nu>1}
       \frac{(\bar u_\nu T)^{m_\nu}}{m_\nu!\prod_{i,t}(d_{\nu,i,t}!)^{m_\nu}}
     \prod_{i,\ 0\le t\le h_i} l_{i,t}!\,A_{i,t}(T)^{p-1-l_{i,t}}
      }{\displaystyle \prod_i A_{i,h_i+1}(T)} .
                                                               \tag{4}
\]
Here the finite sum is over the nonsingleton counts
\(0\le m_\nu<p\) with
\(l_{i,t}=\sum_{\nu:\,w_\nu>1}d_{\nu,i,t}m_\nu\le p-1\)
at every slot. Denominators inside the finite sum are field
units. A sufficiently divisible constant, with
\(v(c^k/k!)>0\) for each \(k\ge1\), leaves (4) unchanged.

To prove (4) consider the integral native coefficient of
\(\prod u_\nu^{m_\nu}\) in \(P_{\sum m_\nu}\):
\[
    C(m)=\frac{\prod_i(\sum_\nu m_\nu\nu_i)!}
                   {\prod_\nu m_\nu!\prod_{\nu,i}(\nu_i!)^{m_\nu}} .
\]
Factoring first by monomial, then by coordinate, expresses this
as the product of single-monomial coefficients and ordinary
multinomial coefficients. Legendre's formula gives the valuation
of the single-monomial coefficient as
\[
 \frac{(w_\nu-1)m_\nu+s_p(m_\nu)-\sum_i s_p(m_\nu\nu_i)}{p-1}
 \ \ge\ (w_\nu-1)v_p(m_\nu!).                                  \tag{5}
\]
The inequality uses \(s_p(m_\nu\nu_i)\le
s_p(m_\nu)s_p(\nu_i)\), obtained by multiplying digits and
then performing carries. Thus a nonsingleton count must be
less than \(p\) for its term to survive. For such counts,
the single-monomial valuation vanishes exactly when the
digit multiplication is carry-free. For a singleton,
multiplication by \(p^b\) just shifts the count digits.
The coordinate multinomial coefficients require that the
combined slots also be carry-free.

For a factorial \(M!\), its unit part modulo \(p\) is
\((-1)^{v_p(M!)}\prod_t M_t!\), where \(M_t\) are the base-\(p\)
digits of \(M\). This follows by separating the factors divisible
by \(p\) and applying Wilson's theorem to complete blocks.
In a surviving \(C(m)\) the signs cancel because the total
valuation is zero. The resulting coefficient at a slot loaded
by \(l\) nonsingleton digits, after summing the allowed
singleton digits at that slot, is
\[
     l!\sum_{j=0}^{p-1-l}\binom{l+j}{j}(1-A_{i,t})^j
                           =l!\,A_{i,t}^{p-1-l}.
\]
This is the multinomial theorem for the singleton digits,
followed by the binomial theorem in characteristic \(p\).
For \(t>h_i\) there is no nonsingleton load, and successive
\(A_{i,t}\) are Frobenius powers of \(A_{i,h_i+1}\).
The tail product of their \((p-1)\)-st powers telescopes to
\(1/A_{i,h_i+1}\). This proves (4).
For \(p=2\) each nonsingleton count is zero or one, and the
capacity test is disjointness of the selected binary digit
slots. For zero-constant quadratics this specializes to the
[ordinary matching identity](QUADRATIC_STRONG_FACTORIAL_WINDOWS_FROM_MATCHINGS.md#2-the-universal-matching-generating-function).

## Which poles can occur at higher precision

The denominator from (3) can be much larger than necessary. In fact
all precisions admit a rational expression whose denominator reduces
to a product of powers of fixed polynomials. At coordinate \(i\),
retain just its unit singleton monomials \(p^b e_i\). If there are
any, let \(h\) be their largest position and put
\[
                U_i(T)=1-\sum_{p^b e_i\ {\rm retained}}
                                   (u_{p^b e_i}T)^{p^{h-b}} .
\]
A coordinate with no unit singleton contributes \(U_i=1\).
Then \(E\) modulo \(p^a\) has an expression with denominator
\(\prod_i U_i^{D_i(a)}\) for computable exponents \(D_i(a)\).
Thus higher-precision terms can change the multiplicities but
require no new residue pole locations.

To see this, remove the retained singleton columns from the support
temporarily. For every other column its count is bounded modulo
\(p^a\): either (5) tends to infinity with the count, or the column
has a nonunit coefficient whose powers tend to zero. At fixed
counts of these columns put \(L_i=\sum_\nu m_\nu\nu_i\).
Their contribution to \(E\) is their integral native coefficient
times \(\prod_\nu(u_\nu T)^{m_\nu}\), multiplied by the product
of the row series
\[
                 \sum_{N\ge0}\binom{N+L_i}{L_i}\,W_N^{(i)}(T).
\]
Here \(W_N^{(i)}\) is the one-color partition polynomial of only
the retained singletons, and the binomial factor is the ratio
\((N+L_i)!/(N!L_i!)\). Since only finitely many \(L_i\) occur,
one can choose \(K=p^r>p^h\) large enough for all these binomials
and the one-color recurrence coefficients to be periodic modulo
\(p^a\). In one color, aggregate by \(N\bmod K\) alone. The recurrence
adds a block of positive size, so after aggregation the \(N=0\)
term is just the initial vector. The shifted row series is a
linear combination of these \(K\) components.

Let \(D(T)\) be the matrix denominator for this one-color recurrence.
Modulo \(p\), the transition adding \(p^b\) has weight \(u_{p^b e_i}\)
if its starting residue is divisible by \(p^b\), and zero
otherwise. This follows by expanding \((1+Z)^{N-1}\) by base-\(p\)
digits in the binomial \(\binom{N-1}{p^b-1}\).
Every edge increases the representative in \([0,K-1]\) unless
it lands at residue zero. Indeed a permitted step starts at a
multiple of its length and cannot pass \(K\). Thus every directed
cycle contains zero. In the determinant expansion, at most one
cycle appears in any contributing permutation. Consequently
\(\overline D=1\) minus the total weight of paths from \(0\) to
\(K\) with aligned steps \(p^b\), with each step weighted by
\(\bar u_{p^b e_i}T\).

For an aligned interval of length \(p^\ell\), its path polynomial
equals the weight of a single step of that length, if allowed,
plus the \(p\)-th power of the path polynomial for length
\(p^{\ell-1}\). In the second case every smaller step lies in a
single one of the \(p\) equal aligned subintervals.
The path weight at \(K\) is therefore
\(\sum_b(\bar u_{p^b e_i} T)^{K/p^b}\), by Frobenius.
It follows that \(D\equiv U_i^{K/p^h}\pmod p\).
For the inverse denominator use the finite geometric expansion
of \(1/D\) around \(1/U_i^{K/p^h}\): their difference in
denominator polynomials is divisible by \(p\), so truncating
after the \(a\)-th term is exact modulo \(p^a\).
This gives the choice \(D_i(a)=aK/p^h\).
A row with no retained singleton has \(W_0^{(i)}=1\) and no
other terms. Taking the finite sum over the bounded counts
proves the claimed common denominator.

## A support window for shared variables of arbitrary degree

Suppose there is at most one singleton in each coordinate, and
every nonsingleton exponent in a coordinate with singleton
\(p^b e_i\) is below \(p^b\). Singleton digits then use only
slots above all nonsingleton loads in that coordinate. Formula
(4), or the carry calculation before it, gives
\[
                  \overline E(T)
                     =\frac{V(T)}{\prod_{\rm singletons}(1-\bar u_\nu T)} ,
\]
where \(V(0)=1\) is a polynomial of degree at most
\((p-1)h\) for \(h\) nonsingletons. In particular, at \(p=2\)
let \(q\) be the number of unit singletons and let \(h\) count
only unit nonsingletons. If the displayed rational function
is nonpolynomial, the [residue denominator argument](QUADRATIC_MATCHING_RESIDUE_DENOMINATOR_WINDOWS.md#proof-of-the-rational-series-criterion)
gives a window of length at most \(\max(q,h)\) at every
positive start. The nonpolynomial condition is automatic
when \(q>h\). This includes arbitrary independent coefficients
and arbitrarily large exponents under the stated local
normalization, and the two sorts of monomial may use the
same variable. The general rationality statement by itself
does not bound the reduced denominator by support size.

## An exact window from a shared digit position

Here is a valuation criterion that also allows the singleton's
coefficient to be slightly smaller in magnitude than the other
natural coefficients. Assume \(c=0\) and that the support has
exactly one singleton \(\sigma=p^b e_i\), with coefficient \(u\ne0\).
Suppose every other column \(\nu\) has a nonzero digit at position
\(b\) in coordinate \(i\). Write \(s_\nu=v(u)-v(u_\nu)\).
If \(s_\nu<1/p\) for every other column, then for every positive
multiple \(n\) of \(p\),
\[
                             v(P_n)=n v(u).                   \tag{6}
\]
Consequently every window of length \(p\) has a nonzero moment.
This is the actual-support bound when the support has at least
\(p\) terms. In particular, at \(p=2\) it gives a two-window
for every support size greater than one under these hypotheses.

To prove (6), use a count vector \(m\) of total size \(n\).
By Legendre's formula the exact valuation from the factorial
coefficient, beyond the contribution
\(\sum_\nu(w_\nu-1)v_p(m_\nu!)\), is
\[
 \frac1{p-1}\sum_j\left(\sum_\nu s_p(\nu_j)s_p(m_\nu)
                -s_p\Big(\sum_\nu m_\nu\nu_j\Big)\right).      \tag{7}
\]
This counts carry operations when the digit products are written
as separate summands. For each nonsingleton write
\(m_\nu=pq_\nu+r_\nu\), \(0\le r_\nu<p\), and put
\(R=\sum_{\nu\ne\sigma}r_\nu\). The initial digit occupancy at
slot \((i,b)\) is at least \(R+(m_\sigma\bmod p)\).
Since \(p\mid n\), this is \(p\lceil R/p\rceil\).
This slot alone must produce at least \(\lceil R/p\rceil\)
carry operations. Incoming carries can only increase that count.
The valuation of the term with count vector \(m\), minus
\(n v(u)\), is therefore bounded below by
\[
             \sum_{\nu\ne\sigma}q_\nu(w_\nu-1-p s_\nu)
                 +\lceil R/p\rceil-\sum_{\nu\ne\sigma} r_\nu s_\nu .
\]
Here \(v_p(m_\nu!)\ge q_\nu\). Since every \(w_\nu\ge2\) and
\(s_\nu<1/p\), the bound is strictly positive as soon as any
nonsingleton is used. The all-singleton term has precisely
valuation \(n v(u)\), as (5) is an equality of value zero
for a singleton. It is the unique least-valuation term.

With just one nonsingleton \(\nu\), the same proof gives the
stronger threshold
\[
                  s_\nu<\min\bigl((w_\nu-1)/p,\ 1/(p-1)\bigr).
\]
Indeed \(\lceil r_\nu/p\rceil=1\) for \(r_\nu>0\).
At \(p=2\) this is an actual two-window for a pure dyadic monomial
\(X_i^{2^b}\) and any nonsingleton whose \(i\)-th exponent has
digit \(b\) set, throughout the stated relative-valuation range.
The other exponent rows are unrestricted.
