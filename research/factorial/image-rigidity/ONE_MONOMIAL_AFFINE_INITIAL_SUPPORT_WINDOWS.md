# Initial support-size windows for one monomial and an affine phase

Let \(F=aX^\beta+c+\sum_{i=1}^t b_i X_i\) over a
characteristic-zero field, with \(a b_1\cdots b_t\ne0\),
\(D=|\beta|\ge2\), and any additional monomial variables allowed.
Coordinates with \(\beta_i=0\) among the linear terms are allowed.
The support size is \(t+1\) if \(c=0\) and \(t+2\) otherwise.

**Theorem.** The first actual-support-size factorial moments cannot
all vanish in either of these cases.

1. \(c=0\), and all \(\beta_i\) at the nonzero linear terms equal
   the same positive integer \(d\). The other exponents of \(X^\beta\)
   are arbitrary.
2. \(D\ge D_*(t)\) for a computable threshold depending only on \(t\).
   The affine coefficients and all monomial exponents are arbitrary.

The first case has no total-degree restriction. By independent
linear adjunction, the second also excludes a shifted support-size
zero window starting at \(n\) whenever
\(D\ge\max_{t\le j\le t+n-1}D_*(j)\).
The support-size window for unrestricted starting orders and
monomial degrees remains open in this argument.

## Exact exponential-moment composition

Scale \(F\) by \(a^{-1}\) and absorb that factor into its affine
coefficients. Abbreviate
\[
 A_r=\frac{\prod_i(r\beta_i)!}{r!},\quad
 A(z)=\sum_{r\ge0}A_rz^r,\quad
 P(z)=\prod_{i=1}^t(1-b_i z),\qquad
 Q(z)=\prod_{i=1}^t(1-b_i z)^{\beta_i}.
\]
Then the normalized moment series is
\[
         B(z)=\sum_{n\ge0}\frac{\Lambda(F^n)}{n!}z^n
               = \frac{e^{cz}}{P(z)} A(z/Q(z)).                \tag{1}
\]
Indeed \(\Lambda(x^r e^{b z x})=r!/(1-bz)^{r+1}\) follows
coefficientwise from the binomial series. Expanding the exponential
of the monomial term and multiplying these individual identities
gives (1). When all the linear terms are independent of \(X^\beta\),
this specializes to the independent linear adjunction identity.
[Rocks, Theorem 3.12 and Corollary 3.13](https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf)

The multinomial allocation among \(e\) separate groups of size \(p\)
gives, for \(n_1+\cdots+n_r=p\), \(n_i>0\) and \(e\ge0\),
\[
             \binom{ep}{en_1,\ldots,en_r}
                      \ \ge\ \binom{p}{n_1,\ldots,n_r}^{e}.
\]
Applying this to all variable exponents shows
\[
          \frac{\prod_i A_{n_i}}{A_p}
                  \le \binom{p}{n_1,\ldots,n_r}^{\,1-D}.        \tag{2}
\]

## Equal attached exponents

The case \(t=0\) is immediate, so put \(p=t+1\ge2\).
Now \(Q=P^d\), and if the first \(p\) moments vanish,
\(P(z)\equiv A(z/P(z)^d)\bmod z^{p+1}\).
The unique formal solution \(P_*(z)\) with constant term one
has \(P_*(z)=A(u)\) for \(u=zA(u)^{-d}\).
Uniqueness follows coefficient by coefficient, since each
unknown coefficient in \(P\) affects the right-hand side only
at later orders. Lagrange inversion gives, with \(k=dp-1>0\),
\[
  [z^p]P_*=\frac1p[u^{p-1}]A'(u)A(u)^{-dp}
                       =-\frac1k[u^p]A(u)^{-k}.               \tag{3}
\]
The coefficient identity can also be obtained by changing the
formal residue variable \(z=uA(u)^d\) and integrating by parts.
The degree of \(P\) is \(p-1\), so a nonzero coefficient in
(3) proves the assertion.

The single-factor contribution to \([u^p]A^{-k}\) is \(-k A_p\).
For \(j=1,\ldots,p-1\) the absolute contribution from \(j+1\)
factors, divided by \(k A_p\), is bounded using (2) by
\[
 T_j=\frac{\binom{k+j}{j+1}}{k}\binom{p-1}{j}
                        \big((p)_{\downarrow j}\big)^{\,1-D}.
 \tag{4}
\]
In fact the multinomial for \(j+1\) positive parts is least when
all but one part equal one, by repeatedly moving one unit from
a smaller part greater than one to a largest part. That value
is \((p)_{\downarrow j}\), and the composition count is
\(\binom{p-1}{j}\).

For \(p\ge4\) the bound (4) is largest at \(d=1,D=p-1\).
First replace \(D\) by \(d(p-1)\le D\).
An increment of \(d\) multiplies the numerator factor
\((dp)_j/(j+1)!\) by at most \(2^j\) and multiplies the
denominator by \(((p)_{\downarrow j})^{p-1}\ge2^{j(p-1)}\).
For \(p=4\) the three bounds in this largest case are
\(3/8,5/72,5/576\), with sum less than one.
For \(p\ge5,d=1,D=p-1\), \(T_1=(p-1)/(2p^{p-3})\le1/10\).
Successive ratios are
\[
 \frac{T_{j+1}}{T_j}
       =\frac{p+j}{(j+1)(j+2)}
                    \frac{p-j-1}{(p-j)^{p-2}} .
\]
For \(1\le j\le p/2\) this is at most
\(2p/(6(p/2)^{p-3})\le4/(3p)\).
For \(j>p/2\), \(j\le p-2\), it is at most
\(8/(p 2^{p-3})\le2/5\).
So the sum of (4) is less than one.

For \(p=3,d\ge2\), the same monotonicity gives the bounds at
\(d=2,D=4\), namely \(2/9,7/216\).
For \(p=3,d=1,D\ge3\) the bounds are at most \(1/3,1/18\).
At \(p=3,d=1,D=2\) the values \(A_1,A_2,A_3=1,2,6\)
give \([u^3]A^{-2}=-12+12-4=-4\).
For \(p=2\) the only correction relative to \(k A_2\) is exactly
\(2d/\prod_i\binom{2\beta_i}{\beta_i}\).
It is less than one since the factor for the attached variable
already exceeds \(2d\) for \(d\ge2\); for \(d=1\) there is
another positive exponent because \(D\ge2\).
This proves the nonzero coefficient in (3) in every case.

## Arbitrary affine coefficients at large total degree

Use complex absolute values after embedding the finitely generated
coefficient field in \(\mathbb C\). If \(t=0\), the monomial
has a nonzero first moment. With a nonzero constant, a zero
first moment forces that constant to be minus the mean of the
ordinary real monomial after normalization. The second moment
is then its strictly positive variance.

For \(t\ge1\), the following power-sum estimate allows a constant
term. There is an explicit \(\epsilon_t>0\) such that complex
numbers \(u_1,\ldots,u_t\) with \(\max |u_i|=1\) have a power
sum of absolute value greater than \(\epsilon_t\) in each of the
ranges \(1,\ldots,t\) and \(2,\ldots,t+1\). One can use
\[
              \epsilon_t=\frac{1}{4t\,3^t(t+1)!(8t)^{t+1}}.
\]
To prove the second assertion, suppose \(|p_j|\le\epsilon_t\)
for \(2\le j\le t+1\), where \(p_j=\sum_i u_i^j\).
The degree-\(t\) polynomial \(\prod_i(1-u_i z)\) agrees through
degree \(t+1\) with
\(\exp(-p_1 z-\sum_{j=2}^{t+1}p_j z^j/j)\).
Its coefficient of \(z^{t+1}\) is zero. On the unit circle
the contribution from changing the exponential \(e^{-p_1 z}\)
is bounded by \(3^t(e^{t\epsilon_t}-1)\), since \(|p_1|\le t\).
Cauchy's coefficient estimate gives
\[
                |p_1|^{t+1}\le (t+1)!3^t\,2t\epsilon_t
                                  <(8t)^{-(t+1)} .
\]
Now \(\sum_{j=1}^t |p_j|/j\le1/(4t)\).
On the unit circle the exponential truncated in its exponent
at degree \(t\) differs from 1 by at most \(1/(2t)\).
The polynomial is its Taylor truncation at degree \(t\);
each of its nonconstant coefficients has modulus at most
\(1/(2t)\). It cannot have a unit-circle zero, contrary
to \(\max|u_i|=1\). The first assertion follows by the same
last argument with all \(p_1,\ldots,p_t\) bounded by \(\epsilon_t\).

Write \(L=t\) for \(c=0\), and \(L=t+1\) otherwise.
Let \(T=\max |b_i|>0\), and \(\kappa_r=[z^r]\log A(z)\).
If the first \(L+1\) native moments vanish, the logarithm of
(1) gives
\[
 -\frac{\sum_i b_i^n}{n}
     =\sum_{r=1}^n \kappa_r[z^{n-r}]Q(z)^{-r}
 \quad (2\le n\le L+1),
 \tag{5}
\]
also at \(n=1\) when \(c=0\).
If \(D\ge L+3\), (2) shows that
\[
                     \tfrac12 A_r\le|\kappa_r|\le 2A_r
                  \qquad(1\le r\le L+1).
 \tag{6}
\]
Indeed every product with at least two factors in the logarithm
has ratio at most \(2^{1-D}\) to \(A_r\), and the total number
of ordered positive compositions of \(r\) is \(2^{r-1}\).

For each nonnegative integer \(e\), the sequence
\((en)!/(n!)^e\) is positive log-convex with first value at
\(n=0\) equal to one. Its next-term ratio is
\(\prod_{j=1}^e (e-(e-j)/(n+1))\), which increases in \(n\).
Consequently \((A_r/(r!)^{D-1})^{1/r}\) is nondecreasing.
For \(j\ge2\),
\[
 \frac{\log(j!)}{j}-\frac{\log((j-1)!)}{j-1}\ \ge\ \frac{1}{4j}.
\]
To see this, the difference is \(1/j\) times
\(\log j-\log((j-1)!)/(j-1)\), and at least half the numbers
in \(1,\ldots,j-1\) are at most \(j/2\).
In particular
\[
 \frac{A_r^{1/r}}{A_L^{1/L}}\le e^{-(D-1)/(4L)}
               \ (r<L),\qquad
 \frac{A_{L+1}}{A_L^{(L+1)/L}}\ge e^{(D-1)/4}.                \tag{7}
\]

Fix \(K=16 L/\epsilon_t\). We claim
\(T<K A_L^{1/L}\) when \(D\) is large enough in terms of \(t\).
If the bound fails, apply the power-sum estimate to \(b_i/T\).
For one of the available orders \(n\le L\) in (5),
the modulus of its left side divided by \(T^n\) is
greater than \(\epsilon_t/L\).
On the right,
\[
 |[z^h]Q^{-r}|\le \binom{D r+h-1}{h}T^h
\]
by the binomial series and \(\sum_{i=1}^t\beta_i\le D\).
The term \(r=L\), if present, is at most \(2/K^L\).
All other terms together are at most
\(2L((D+1)(L+1))^L e^{-(D-1)/(4L)}\),
by (6)–(7). For large \(D\) their sum is less than
\(\epsilon_t/(2L)\), a contradiction.

Finally use (5) at \(n=L+1\), move its term \(\kappa_{L+1}\)
to the left, and use the bound on \(T\).
After dividing by \(A_L^{(L+1)/L}\), (6)–(7) bound the
left side below by \(\tfrac12 e^{(D-1)/4}\).
The right side is at most
\((t+2L)K^{L+1}((D+1)(L+1))^{L+1}\).
The exponential exceeds this polynomial for all sufficiently
large \(D\), giving the contradiction.
Both tests for sufficiently large \(D\) are effective:
the relevant polynomial times decreasing exponential
is decreasing beyond a computable threshold, and rational
enclosures find a first strict verified bound on that range.
Taking the larger choice for \(L=t,t+1\) gives \(D_*(t)\).

The shifted consequence uses the
[exact initial-window lift](SFC_SHIFT_TO_INITIAL_LINEAR_ADJUNCTION.md).
It adds at most \(n-1\) independent linear terms, preserves
the monomial and its degree, and preserves whether the constant
term is present. The maximum of the corresponding initial
thresholds therefore suffices for the original shifted window.

## Nonvanishing under compatible affine tilts

For arbitrary total degree \(D\ge2\), (1) also rules out eventual
vanishing for one nonlinear monomial plus any affine phase.
In fact \(A_r\ge(r!)^{D-1}\) by the same multinomial comparison
with \(r\) parts equal to one. Thus \(A\) has zero radius of
convergence. The substitution \(z/Q(z)\) has derivative one at
the origin and has an analytic inverse there. If \(B\) converged
near zero then so would \(A\), by (1). In particular \(B\) is
not polynomial. The
[uniform propagation theorem](FACTORIAL_INHOMOGENEOUS_UNIFORM_PROPAGATION.md)
gives a coefficient-uniform finite nonvanishing window for
this whole class at each fixed \(D,t\), since there are at most
\(D+t\) active variables.

The same substitution works for a homogeneous nonlinear phase
\(H\) of degree \(D\) in place of the monomial, provided every
monomial of \(H\) has the same total degree \(d_g\) within each
group of variables whose nonzero linear coefficients coincide.
Each group's contribution to \(Q\) in (1) is then
\((1-b_g z)^{d_g}\), while \(A_r=\Lambda(H^r)/r!\).
This follows by applying the variablewise weighted exponential
identity from (1) to each monomial of \(H^r\).
If \(H\) has infinitely many nonzero moments, \(A\) has
zero radius of convergence: otherwise the additional
factorial denominators in the simplex normalization
\(\Lambda(H^r)/(Dr+m-1)!\) would make its generating series
entire. That normalized period is regular singular at
infinity by the
[homogeneous propagation theorem](FACTORIAL_HOMOGENEOUS_UNIFORM_PROPAGATION.md).
An entire solution with a regular singularity at infinity has
polynomial growth and is polynomial, by the local Frobenius
expansion and Cauchy estimates. This contradicts the assumed
infinitely many nonzero moments. The same analytic inverse
argument now gives non-eventual vanishing and a finite
window for \(H\) plus the compatible affine tilt.
