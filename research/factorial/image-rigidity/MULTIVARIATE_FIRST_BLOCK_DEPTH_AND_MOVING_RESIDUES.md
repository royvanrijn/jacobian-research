# First-block prime depths with arbitrary moving rows

At any two adjacent orders an ordinary first factorial block
gives a prime-depth obstruction in the physical coefficient
coordinate. It applies at primes that divide several rows.
The same inequality leads to all-degree exclusion for an
unbounded set of other row types at arbitrarily large heights
when their residues modulo some first-block integers are
controlled.

Write \(n,N=n+1\) for the orders and select \(r\ge1\) copies
of a row with smaller exponent \(L\) and positive absolute gap
\(d\). Reverse both columns if needed so this row is increasing.
For another row type \(\rho\) let \(\ell_\rho,d_\rho,r_\rho\)
be its smaller exponent, absolute gap and multiplicity, and let
\(L_\rho\) be the least common multiple of the two full ramps
\(j\ell_\rho+t\), \(j=n,N,\ 1\le t\le jd_\rho\).
Equal rows can be dropped. There is no restriction on the
number of other types, their heights, or their signs.

## The exact depth inequality

Choose one of the orders \(m\), and let \(m'\) be the other.
At a prime \(p\) put
\[
 b_p=\max_{1\le t\le m'd}v_p(m'L+t),\qquad
 R_p=\sum_\rho r_\rho v_p(L_\rho),\qquad
 W_p=\sum_{t=1}^d\bigl(v_p(mL+t)-b_p\bigr)_+ .
\]
A physical common root requires, at every prime,
\[
             rW_p\ \le\ (m+1)(r b_p+R_p)+m v_p(m'!).            \tag{1}
\]
The first block can use several factors at once. For the
individual factor \(F=mL+h,\ 1\le h\le d\), (1) in particular
implies
\[
 r v_p(F)\ \le\ (m+2)r b_p+(m+1)R_p+m v_p(m'!).                 \tag{2}
\]

Strip the factorial at the smaller endpoint in each order,
which changes the coefficients only by an overall scalar.
At index \(k\) the selected row uses the first \(kd\) factors.
At each \(e\le b_p\) the valuation layer satisfies
\[
 \left\lfloor\frac{jL+kd}{p^e}\right\rfloor
      -\left\lfloor\frac{jL}{p^e}\right\rfloor
      -\frac{kd}{p^e}
   =\left\{\frac{jL}{p^e}\right\}
      -\left\{\frac{jL+kd}{p^e}\right\}.
\]
Its oscillation is at most one, and the subtracted linear
slope is the same at the two orders. Remove this common
slope for the \(b_p\) layers of the selected copies. By the
[common slope removal for the other rows](MULTIVARIATE_RETAINED_COMMON_SLOPE_BUDGET.md#arbitrary-retained-rows-at-one-prime),
one can also remove their common linear slopes, with total
oscillation at most \(R_p\) at each order. These changes move
both Newton polygons by the same slope.

At order \(m'\) there are no higher layers in the selected
row. With the binomial coefficients included, the reduced
ordinates have oscillation at most
\(E+v_p(m'!)\), where \(E=rb_p+R_p\).
Their lower convex minorant takes values inside that same
range. Every unit slope is consequently at most
\(E+v_p(m'!)\), since it is the difference between minorant
values at consecutive integer positions.

At order \(m\), the higher layers contribute zero at
index zero and at least \(rW_p\) at every positive index.
The binomial coefficient has valuation zero at index zero
and nonnegative valuation elsewhere. Subtract any candidate
common slope \(\lambda\le E+v_p(m'!)\) times \(k\) from the
reduced ordinates at order \(m\). The difference of a
positive-index ordinate from index zero is at least
\(rW_p-E-m(E+v_p(m'!))\).
If (1) fails, index zero is the unique minimizer, so there
is no face of that slope. The Newton polygon product rule
then excludes every common factor. This proves (1) and (2).

## A joint test with no chosen row type

At fixed \(p,m,m'\), replace \(b_p\) by the separate unmarked
maximum \(b_{\rho,p}\) for each row type, including the
previously selected type. Set \(E=\sum_\rho r_\rho b_{\rho,p}\).
For each row form the marked ramp array
\[
 v_{\rho,t}=(v_p(m\ell_\rho+t)-b_{\rho,p})_+,\qquad
                     1\le t\le m d_\rho .
\]
The net higher-layer change from index zero to \(k\) is
\[
 U(k)=
 \sum_{\rho\ {\rm increasing}}r_\rho\sum_{t=1}^{kd_\rho}v_{\rho,t}
 -
 \sum_{\rho\ {\rm decreasing}}
          r_\rho\sum_{t=(m-k)d_\rho+1}^{m d_\rho}v_{\rho,t}.
\]
If
\[
              U(k)>E+k(E+v_p(m'!))\qquad(1\le k\le m),        \tag{3}
\]
there is no common factor. One can also reverse both columns
and use (3) on the reversed profile. The proof removes the
first \(b_{\rho,p}\) layers for every row by the same floor
identity, with signs reversed in the affine slope for
decreasing rows. The error has oscillation at most \(E\),
and all remaining layers vanish at the unmarked order.
At the marked order, their exact change from index zero
is \(U(k)\). The unique-minimizer argument above now
gives (3). Positive contributions may come from many
different first blocks at the same prime.

For example, take two increasing unit-gap rows at smaller
exponent \(L=(2^{20}-1)/3\) and two decreasing unit-gap
rows at smaller exponent one. At orders \(2,3\) use
\(m=3,p=2\). The selected first block is \(2^{20}\),
while \(b_2=2\) and \(R_2=4\). Inequality (1) would
require \(36\le35\). The sole prime of the selected
first block is shared by the other rows.

## An integer certificate without prime factorization

For one factor \(F=mL+h\) in the first block, choose a nonzero
integer \(B_u\) for each factor \(u\) in the selected row's
order-\(m'\) ramp and in all the other rows' two ramps. Require
\(\gcd(F,u)\mid B_u\). These certificates can be found from
small congruence combinations or integer polynomial Bezout
identities. Write \(D_0\) for their least common multiple
in the selected order-\(m'\) ramp, \(D_\rho\) for the least
common multiple in row type \(\rho\), and \(D\) for the least
common multiple of \(D_0\) and all \(D_\rho\). Absolute values
are understood in each least common multiple.

Let \(F_{\rm deep}=\prod_{v_p(F)>v_p(D)}p^{v_p(F)}\).
A common root requires the integer divisibility
\[
       F_{\rm deep}^{\,r}\mid
       D_0^{(m+2)r}(m'!)^{m}
                 \prod_\rho D_\rho^{(m+1)r_\rho}.              \tag{4}
\]
Indeed at each prime in this product defining \(F_{\rm deep}\),
every certified ramp factor has \(v_p(u)<v_p(F)\) and hence
\(v_p(u)\le v_p(B_u)\). Insert these upper bounds in (2).
Moreover \(F/\gcd(F,D)\mid F_{\rm deep}\), so replacing
\(F_{\rm deep}\) in (4) by \(F/\gcd(F,D)\) gives a test without
prime factorization. At one prime outside \(D\) with \(p>m'\),
(2) already excludes a common factor independently of all
multiplicities, as soon as \(p\mid F\).

Only the selected row's order-\(m'\) ramp occurs among its
certificates. Arbitrarily high prime powers at its other
order-\(m\) factors cause no loss. The bound on another row
uses one maximal depth at each prime, independent of the
number of times a small prime divides its ramp.

## All-degree exclusion for moving residues

Fix positive integers \(d,d_0\). Let \(J>3\) be an integer
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Keep a selected row with smaller exponent \(L\), gap \(d\)
and arbitrary positive count \(r\). Write
\[
             Q(L)=\prod_{2\le m<J}(mL+1).
\]
For every other row type, allow any smaller exponent of
the form
\[
       \ell_\rho=c_\rho L+b_\rho+q_\rho Q(L),\qquad
       |c_\rho|,|b_\rho|\le W,\quad 1\le d_\rho\le d_0,         \tag{5}
\]
with arbitrary integer \(q_\rho\) subject to \(\ell_\rho\ge0\).
The signs, multiplicities and set of types can all depend on
\(L\). The alternative condition giving exactly the same
proof is a potentially different representation
\(\ell_\rho\equiv c_{\rho,m} L+b_{\rho,m}\pmod{mL+1}\)
for each \(2\le m<J\), with the same bounds. Use the
parameters for the current mark \(m\) in that case.

For each adjacent pair \(2\le n<N<J\) and each
\(m\in\{n,N\}\), assume the noncoincidence inequalities
\[
        m(jb_\rho+t)-jc_\rho\ne0
        \quad(j=n,N,\ 1\le t\le jd_\rho).                     \tag{6}
\]
For instance all these inequalities hold if \(c_\rho=1\)
and \(1\le b_\rho\le W\). A prescribed fixed core is also
allowed by using \(c_\rho=q_\rho=0,\ b_\rho=\ell_\rho\) and
increasing the bounds \(d_0,W\) to cover it.

There is a constant \(C\), depending only on \(J,d,d_0\),
such that every sufficiently large \(L\) satisfies
SFC(2) in every degree in either of the following regimes.

1. \(W=o((\log\log L)^2/\log\log\log L)\), with arbitrary
   independent multiplicities.
2. \(W\le c\log L\), where \(c>0\) depends only on
   \(J,d,d_0\), and \(\sum_\rho r_\rho\le A r\) for any
   fixed \(A\). The height threshold is computable from
   \(A,J,d,d_0\).

In both regimes classification of the other row types
takes place at the actual integer \(L\). In particular
it permits a growing number of residues and arbitrary
additional height scales through \(q_\rho\).

To prove this, take \(F=mL+1\). A factor \(u=j\ell_\rho+t\)
in another row has
\[
                   m u\equiv m(jb_\rho+t)-jc_\rho\pmod F.
\]
The nonzero integer on the right is a certificate of
absolute value at most \(C(W+1)\). A selected-row factor
\(u=m'L+t\) has certificate \(mt-m'\), which is nonzero because
\(m,m'\ge2\) are consecutive and coprime. The certificate
bound holds uniformly over both orders, all types and all
\(q_\rho\).

For regime 1, apply
[Pasten, Corollary 1.5](https://arxiv.org/pdf/2312.03566)
to the coprime positive integers \(1,n(NL+1)\), whose
sum is \(N(nL+1)\). Their product with the sum equals
\(nN(nL+1)(NL+1)\). Its greatest prime factor is bounded
below by a positive constant times
\((\log\log L)^2/\log\log\log L\), uniformly over the
finitely many pairs. Hence above a height threshold a
prime \(p>J+C(W+1)\) divides one of the two first factors
\(nL+1, NL+1\). Choose its corresponding order as \(m\).
It divides none of the certificates, and (2) excludes.
With the elementary depth test alone, any value of \(L\)
with such a prime has the same conclusion for that pair.

For regime 2, take \(m=N\), \(F=NL+1\) at every pair.
Since \(D\) divides the least common multiple of all
positive integers up to \(C(W+1)\),
\(\log D\le 4C(W+1)\log 2\). To see the elementary bound,
at most one power of each prime lies in \((s,2s]\), and
each such power contributes one factor of its prime to
\(\binom{2s}{s}\). This bounds the increment of the
logarithm of the least common multiple by \(2s\log2\).
Sum over dyadic intervals up to the least power of two
above \(C(W+1)\).
Choose \(c\) small enough that \(\log D\le \tfrac12\log L\)
for large \(L\). The log of the left side of (4), divided
by \(r\), is then at least \(\tfrac12\log L\), using
\(F/\gcd(F,D)\). On the right \(\log D_0\) and
\(\log(m'!)\) are bounded in terms of \(J,d\), and
\[
      \log D_\rho\le (n+N)d_0\log(C(W+1)).
\]
Its total log divided by \(r\) is \(O_{J,d,d_0,A}(\log(W+2))\).
For all sufficiently large \(L\) these bounds contradict
(4), with a computable threshold.

The first moment has no common root with the second by
even positivity on the real axis. Terminal orders at or
above \(J\) are covered by the absolute cutoff. Thus the
same height threshold covers every degree and gives
coprimality over any characteristic-zero field. The residue
conditions are sufficient for all-degree closure. For arbitrary
profiles the paired-order conditions (1), (3) and (4) apply
directly to their integer factorial blocks.
