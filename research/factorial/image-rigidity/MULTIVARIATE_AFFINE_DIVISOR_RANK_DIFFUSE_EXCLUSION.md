# Exact diffuse exclusion from affine factorial divisor rank

Use the [signed-balanced affine template](MULTIVARIATE_BALANCED_AFFINE_JET_DIFFUSE_EXCLUSION.md)
with fixed ordinary core \(C\), adjacent orders \(n,N=n+1\), \(n\ge2\),
and appended types
\[
              (\beta_i M+a_i,\ \beta_i M+b_i)
\]
of multiplicities \(r c_i+s_i\). Here the slopes \(\beta_i>0\) and
offsets are rational, the gaps \(d_i=b_i-a_i\) are nonzero integers,
\(c_i\ge0\) and \(s_i\) are integers, and
\(\sum_i c_i d_i=0\), with at least one positive \(c_i\).
Admissible \(M>0,r\) give ordinary endpoints and multiplicities.
Write \(u=1/M\), \(A=\sum_i c_i d_i^2/\beta_i\),
\(B=\sum_i c_i d_i/\beta_i\), \(h=\sum_i s_i d_i\),
\(\tau=Aru\) and \(\Lambda=N/n\). The positive diffuse reference
is along the signed-bias slope \(B/A\).

There is a finite divisor matrix for this template. Its rows are
indexed by coefficient positions \((n,1),\ldots,(n,n-1)\) and
\((N,1),\ldots,(N,N)\). Its columns record the exponent of each
distinct linear factorial factor in the rational functions
\[
 \mathcal F_{j,k}(M)=
 \prod_i
 \frac{
      \bigl(\exp R^{(i)}_{j,k}(1/M)\bigr)^{n c_i}}
      {\bigl(\exp R^{(i)}_{n,n}(1/M)\bigr)^{k c_i}},
 \quad
 R^{(i)}_{j,k}(u)=\sum_{\ell=1}^{k d_i}
                 \log\!\left(1+\frac{a_i+\ell/j}{\beta_i}u\right).
                                                               \tag{1}
\]
Negative upper summation limits use the signed sum, so the
exponentials in (1) are literal finite rational products. The
column divisors are taken over \(\mathbb Q(M)\), disregarding
nonzero constants and counting a possible divisor at \(M=0\)
in the same way as the others.

**Theorem.** If the matrix in (1) has full row rank, then every
positive physical Gaussian reference contact has a computable
ordinary exclusion neighborhood for this template. That is, for
each positive exceptional strength \(\tau_*\) there exist
\(\varepsilon,\delta>0\) such that the original consecutive
moments are coprime whenever \(0<u<\varepsilon\) and
\(|\tau-\tau_*|<\delta\). The conclusion applies to multiple
common roots and nonisolated loci of the two first limiting
directions. Fixed offsets \(s_i\) in the row counts are allowed.

The full-rank test can be replaced by the proper projection test
below. At all paired orders below the absolute cutoff, if such
tests pass and the zero-strength pairs
\(Q_n^C(z),Q_N^C(\Lambda^h z)\) are coprime, the template is
regular in every degree for every admissible multiplicity
parameter at all sufficiently large \(M\).

## Algebraic projection of a common factor

Let \(L>0\) clear the denominators of
\[
 f_{j,k}^0=\frac{k(k-j)}{2j}
                +\frac BA\,\frac{k}{2}\left(\frac1j-\frac1n\right).
\]
Set \(b_{j,k}=[z^k]Q_j^C\,\Lambda^{hk\boldsymbol 1_{\{j=N\}}}\).
Let \(\boldsymbol Y\) range over the coefficient multiplier
torus with the positions of (1) free and
\(Y_{n,0}=Y_{N,0}=Y_{n,n}=1\). For \(X\ne0\), consider
\[
                   \sum_{k=0}^j b_{j,k}Y_{j,k}
                                          X^{L f_{j,k}^0}\,z^k .
                                                               \tag{2}
\]
Let \(Z\) be the Zariski closure of the projection onto
\(\boldsymbol Y\) from the locus where (2) has a common factor
of degree at least two. Then \(Z\) is a proper algebraic subset.
Indeed, fixing \(X\) and the lower coefficients leaves finitely
many possible monic quadratic divisors of the lower polynomial,
including multiplicity. For each such divisor \(H\), requiring
the upper polynomial to vanish modulo \(H\) imposes two
independent linear conditions in its free coefficients.
The classes of \(z,z^2\) are independent in \(\mathbb C[z]/H\)
because \(H\) has nonzero constant term. The incidence therefore
has codimension at least two for fixed \(X\). Allowing the one
coordinate \(X\) makes its image of codimension at least one.

Write \(T_{\rm div}\) for the subtorus that is the image of the
monomial map with columns given by (1). The relative version of
the test is
\[
                              T_{\rm div}\not\subset Z.        \tag{3}
\]
This is finite algebraic elimination of the Sylvester minors in
the Laurent variables, followed by the known monomial substitution.
Full row rank makes \(T_{\rm div}\) the entire multiplier torus,
so it implies (3). In particular \(Z\) depends on the core and
the template through \(f^0,h\), with no choice of exceptional
strength needed.

For later use, the resultant in \(X\) of (2) is nonzero for
every fixed \(\boldsymbol Y\) in the complex torus, even at
nonpositive or nonreal coefficient tuples. Put
\(q=X^{-L/(2nN)}\) on the positive real axis. Apart from the
upper argument scale \(q^{B/A}\), the Gaussian exponents are
\((nN/j)k(j-k)\). If \(B\ne0\), their polynomials before that
scale tend to endpoint binomials as \(q\downarrow0\) and all
root moduli are bounded above and away from zero. The upper
scale then separates the moduli. If \(B=0\), let \(q\to+\infty\).
The [two-term asymptotics at large clock](MULTIVARIATE_GAUSSIAN_TWO_CLOCK_RESULTANT_AT_INFINITY.md#exact-leading-terms-and-the-two-clock-degrees)
place the degree-\(j\) roots at the distinct modulus scales
\[
                  q^{(nN/j)(2k-j-1)},\qquad 1\le k\le j.
\]
That proof uses a simple nonzero binomial root after rescaling
each two-term face and works with the present nonzero complex
coefficients. The exponent sets in the two orders are disjoint.
One order is even and gives odd exponents, because its partner
is odd. The other gives even exponents. Thus again there is
coprimality for a real interval of \(X>0\).

## Prime specialization of the exact template

For every distinct linear divisor of (1) choose a separate prime
that divides none of the finitely many constant multipliers,
its leading coefficient or any cross determinant with the other
linear divisors. At each chosen prime impose a residue on \(M\)
at which its own affine factor has valuation one and the other
factors are units. This is possible by choosing first its
unique zero modulo the prime and then one of the lifts
modulo its square where it is nonzero. The Chinese remainder
theorem gives arbitrarily large positive integers \(M\)
satisfying all the conditions. The rational functions (1) are
positive for such \(M\) large enough. Their prime valuation
vectors at the chosen primes are the respective divisor
columns. The specialization here is of rational functions;
it is permitted even if one of those sampling values \(M\)
gives nonintegral formal row endpoints.

Put \(f_{j,k}(u)=\log\mathcal F_{j,k}(1/u)/(nAu)\), continued
analytically at zero. It is exactly the repeated-block multiplier
after a common lower endpoint normalization, as in the
[balanced affine normalization](MULTIVARIATE_BALANCED_AFFINE_JET_DIFFUSE_EXCLUSION.md#exact-physical-normalization-and-its-first-jet).
For each of the sampled \(u=1/M\), the complex exponential
curve \(\bigl(\exp(y f_{j,k}(u))\bigr)_{j,k}\) has Zariski closure
containing \(T_{\rm div}\). In fact monomials with the same
exponential frequency give a multiplicative relation among
the positive rationals (1), and its exponent difference
pairs to zero with every selected prime valuation column.
Hence both monomials agree on \(T_{\rm div}\). If a Laurent
polynomial vanishes on the exponential curve, the coefficients
sum to zero within each group of equal frequencies by the
linear independence of distinct exponentials. The polynomial
then vanishes on \(T_{\rm div}\).

## Exact elimination at positive strengths

At a positive physical reference contact \(\tau_*\), the clock
\(q_*=\exp[-\tau_* /(2nN)]\) is algebraic by the nonzero
Laurent reference resultant, so \(\tau_*\) is transcendental.
Use the analytic arrays
\[
 b_{j,k}\exp\left(\tau_* f_{j,k}^0+G_{j,k}(u)
                  +v f_{j,k}(u)
                  +t\frac{f_{j,k}(u)-f_{j,k}^0}{u}\right),      \tag{4}
\]
where \(G_{j,k}(u)\) is the known log multiplier of the
fixed additive offsets \(s_i\) and \(G_{j,k}(0)=0\).
The actual tail uses \(v=\tau-\tau_*,\ t=\tau_*u\).
The quotient in (4) is analytic. All Taylor coefficients
of (4) are algebraic by the rational row log series and
the algebraic clock. On \((u,t)=(0,0)\), the physical
resultant is a nonzero function of \(v\), as proved above.

Apply the [contour eliminant](MULTIVARIATE_ISOLATED_TWO_CLOCK_DIFFUSE_EXCLUSION.md#exact-elimination-of-an-accumulating-contact)
to the resultant in \(v\) and the Sylvester minors detecting
common factor degree at least two. A physical common root
of the ordinary polynomials is nonreal by strict even
positivity and supplies its conjugate as well. The
contour eliminant has algebraic Taylor coefficients in
\((u,t)\). If actual common roots accumulated as \(u\to0\)
and \(\tau\to\tau_*\), this eliminant would vanish identically
on the germ \(t=\tau_* u\), and hence identically in
\((u,t)\). The last step follows coefficient by coefficient
from the transcendence of \(\tau_*\).

An identically zero eliminant would give, for every small
\((u,t)\), a corank-two contact at some \(v\) in the fixed
bounded contour disk. Fix one of the sampled positive
\(u=1/M\) sufficiently close to zero. At this \(u\), let
\(y=v+t/u\) and \(X=\exp[(\tau_*-t/u)/L]\). The arrays
(4) are exactly (2) with
\[
                      Y_{j,k}=\exp G_{j,k}(u)\,
                                              \exp(y f_{j,k}(u)).
\]
Take a Laurent polynomial \(P\) vanishing on \(Z\) whose
restriction to \(T_{\rm div}\) is nonzero, as permitted by
(3). This nonzero restriction persists after translation
by \(\exp G(u)\) for all sufficiently small \(u\), by
evaluating at a fixed witness point of \(T_{\rm div}\).
The specialization result implies
\(P(\exp G(u)\exp(y f(u)))\) is a nonzero entire function of
\(y\). It has finitely many zeros in any compact disk.
As \(t\) runs over a small real interval, all the candidate
\(v\)'s, hence the \(y\)'s, stay in a compact disk at the
fixed \(u\). One \(y\) must therefore work at infinitely
many different \(t\), hence at infinitely many different
nonzero \(X\). This contradicts the nonzero resultant
in \(X\) for each fixed coefficient multiplier tuple.

The nonzero germ furnished by this contradiction gives
computable neighborhoods by the same finite Taylor search
and analytic remainder bounds as in the
[isolated-template proof](MULTIVARIATE_BALANCED_AFFINE_JET_DIFFUSE_EXCLUSION.md#algebraic-isolation-and-exact-exclusion).
The all-degree assertion uses bounded \(\tau\) at any
physical failure by the absolute curvature theorem,
finitely many orders below the absolute angular cutoff,
and the assumed zero-strength tests. These are the same
compact limit and coefficient estimates as for the
balanced affine template.

## An infinite family beyond the single-mark test

Fix any \(r_0\ge2\), let \(K=r_0+1\), and take equal main
multiplicities of the four types
\[
 (M,M+r_0),\quad(M+1,M),\quad
 (M+K+r_0,M+K),\quad(M+K,M+K+1).
\]
This is the [signed-balanced obstruction to the single-mark test](MULTIVARIATE_INDEPENDENT_HEIGHT_WINDOW_CLOSURE.md#a-balanced-unit-gap-obstruction-to-the-single-mark-count-test)
at orders \(2r_0,2r_0+1\). The divisor test here has full rank
at every consecutive pair \(n,N=n+1\), \(n\ge2\). Its limiting
signed bias slope and its tail signed gap are both zero
when the four multiplicities are equal. Hence above a computable
height the nonempty template by itself, or the template appended
to any fixed ordinary core already satisfying SFC(2), is regular
in every degree.
For an empty core, the zero-strength limit is covered by the
[uniform weak-curvature theorem](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md).

To verify the rank, suppose a linear combination of
the logarithmic divisor rows vanishes. Write its
coefficients on the raw logarithmic arrays at each order
as \(\alpha_{j,k}\), absorbing the common normalization into
the coefficient at \((n,n)\), and divide out the common
factor \(n\) from (1). Put
\(\alpha_{j,0}=-\sum_{k\ge1}\alpha_{j,k}\), which has no
effect on the raw logarithms since the index-zero multiplier
is one. Define \(A_{j,l}=\sum_{k\ge l}\alpha_{j,k}\)
for \(1\le l\le j\).

The two offset windows are disjoint. Set
\(l_0=\lfloor j/r_0\rfloor+1\). In the first window, for each
\(l\ge l_0\) choose
\((l-1)r_0<h'\le lr_0\) with \(h'>j\) and \(j\nmid h'\).
At \(l=l_0\) the choice \(h'=j+1\) works, and later blocks
contain at least two consecutive integers above \(j\).
The factor \(jM+h'\) has no contribution from the unit-gap
row, and the other order cannot share its divisor since
\(h'/j\) is nonintegral and the orders are coprime.
Its exponent is therefore \(A_{j,l}=0\).

The second window is the reflection of the first. Reflection
at the constant-normalized raw coefficients replaces the
array \(\alpha_{j,k}\) by \(\alpha_{j,j-k}\): the subtraction
of the index-\(j\) log factor produces the specified
\(\alpha_{j,0}\) term. Its tail sum at \(l\) is
\(-A_{j,j-l+1}\). Applying the same divisors to this window
forces \(A_{j,k}=0\) for \(k\le j-\lfloor j/r_0\rfloor\).
Because \(r_0\ge2\), the two ranges cover \(1,\ldots,j\).
This holds at both orders, so all \(\alpha_{j,k}\) are zero.
