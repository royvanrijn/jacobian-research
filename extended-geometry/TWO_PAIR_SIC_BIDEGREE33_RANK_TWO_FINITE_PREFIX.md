# The isolated rank-two cubic finite-prefix component

## 1. The exact point

Work on the normalized anti-Weyl chart

\[
s_4=-s_2,\quad s_5=s_1,\quad s_6=-s_0,\qquad
t_2=0,\quad t_3=t_1,\quad t_4=-t_0.
\tag{1.1}
\]

All odd moments vanish identically.  On this chart coefficient rank at
most two is cut out by

\[
\begin{aligned}
R={}&-9s_0s_2+9s_1^2-9s_2^2-9s_3^2-8s_3-t_0^2+4t_1^2+1,\\
I={}&-9s_0s_3+s_0+18s_1s_2+9s_2s_3+9s_2-4t_0t_1.
\end{aligned}                                             \tag{1.2}
\]

Let \(c\) be the rational vector represented by the terminating decimals

\[
(-0.8702368803753021,-1.8595344299498098,
-0.3028762025151954,1.5192025917402650,
-2.3287623389711780,-1.5350489449499409)
\tag{1.3}
\]

in the coordinate order \((s_0,s_1,s_2,s_3,t_0,t_1)\), and put

\[
X=\prod_{j=1}^6[c_j-10^{-10},c_j+10^{-10}].              \tag{1.4}
\]

The exact representation of the point is the system

\[
R=I=\mu_2=\mu_4=\mu_6=\mu_8=0                           \tag{1.5}
\]

together with the rational isolating box \(X\).  No floating-point root
acceptance is part of the proof.

## 2. Certificate

The checker constructs the primitive anti-Weyl moments directly from the
contraction formula and performs a Krawczyk test using rational interval
arithmetic.  If \(J\) is the Jacobian of (1.5), it takes the exact rational
matrix \(A=J(c)^{-1}\) and verifies

\[
c-Af(c)+(1-AJ(X))(X-c)\subset\operatorname{int}(X).       \tag{2.1}
\]

Thus (1.5) has a unique real zero \(\xi\) in \(X\).

In the row-major coefficient matrix, all sixteen \(3\)-by-\(3\) minors
reduce to zero modulo \((R,I)\).  The minor on rows and columns
\(\{1,2\}\) is bounded away from zero on \(X\).  Hence

\[
\operatorname{rank}C(\xi)=2.                              \tag{2.2}
\]

On the smooth twelve-dimensional determinantal chart using the same
minor, rational interval evaluation exhibits a nonzero \(8\)-by-\(8\)
minor of the differential of \((\mu_1,\ldots,\mu_8)\).  Its rank is at
least eight.  At a moment-zero point, the three infinitesimal diagonal
\(\mathrm{SL}_2\) directions and radial scaling lie in the kernel.  They
are independent here: the normalized non-null
\(\operatorname{Sym}^2\) component has one-dimensional stabilizer, and
the nonzero-weight sextic coordinates break that stabilizer; radial
scaling cannot be an infinitesimal \(\mathrm{SL}_2\) direction because
the quadratic discriminant is nonzero.  Therefore the rank is exactly
eight and these four directions exhaust the tangent kernel.  The induced
projective GIT point is reduced and isolated.

The same interval evaluation gives, after removal of the positive integer
contents of the moments,

\[
\boxed{\mu_{10}(\xi)>0,\qquad \mu_{12}(\xi)<0,\qquad
\mu_{14}(\xi)>0.}                                        \tag{2.3}
\]

The non-null normalized \(\operatorname{Sym}^2\) discriminant is a
nonzero invariant, so \(\xi\) is semistable.  We have proved:

> **Theorem 2.1.** The rank-at-most-two nine-moment fiber contains an
> exact-rank-two semistable reduced isolated projective point.  Its first
> missing moment is already nonzero.  Consequently this component does
> not lift to either \(Z(\mu_1,\ldots,\mu_8,\mu_{12})\) on the rank-two
> locus or \(Z(\mu_1,\ldots,\mu_{12},\mu_{14})\) in the ambient space.

This realizes the semistable point previously forced only by the
rank-stratified Hilbert-series calculation.  It is a finite-prefix
component, not an all-order moment-zero point.

## 3. Smaller corrected-component charts

For the remaining anti-Weyl classification set

\[
x=s_1,\quad y=s_3,\quad A=s_0^2,\quad B=s_2^2,\quad
C=t_0^2,\quad D=t_1^2,\quad p=s_0s_2.
\]

The equations \(\mu_2=R=0\) give

\[
\begin{aligned}
A={}&-15B+\frac{14}{3}C+\frac{56}{3}D-6x^2-10y^2-\frac{70}{3},\\
p={}&-B-\frac19C+\frac49D+x^2-y^2-\frac89y+\frac19.
\end{aligned}                                             \tag{3.1}
\]

The remaining rank equations are

\[
p^2=AB                                                     \tag{3.2}
\]

and

\[
16CD=(1-9y)^2A+81(2x+y+1)^2B
 +18(1-9y)(2x+y+1)p.                                     \tag{3.3}
\]

After (3.1), the primitive corrected moments have the following sparse
profiles in \((x,y,B,C,D)\):

| moment | terms | degree |
| --- | ---: | ---: |
| \(\mu_4\) | 39 | 4 |
| \(\mu_6\) | 119 | 6 |
| \(\mu_8\) | 294 | 8 |
| \(\mu_{12}\) | 1218 | 12 |

Exact characteristic-zero Gröbner calculation is now immediate on this
quotient:

\[
\bigl(p^2-AB,\ (3.3),\ \mu_4,\mu_6,\mu_8,\mu_{12}\bigr)
=\mathbb Q[x,y,B,C,D].                                   \tag{3.4}
\]

The `msolve` replay reduces 955 pairs and returns the one-element basis
\([1]\).  Because the quotient equations are necessary for every original
anti-Weyl point, possible squaring artifacts do not weaken the exclusion.

> **Proposition 3.1.** The corrected exact-rank-two system
> \(Z(\mu_1,\ldots,\mu_8,\mu_{12})\) has no point on the normalized
> anti-Weyl chart.

Off the exceptional Hurwitz boundary, a rank-two cubic pencil can instead
be normalized to

\[
A_1=X^2(X+T),\qquad A_2=T^2(\lambda X+T),                 \tag{3.5}
\]

with an arbitrary \(2\)-by-\(4\) channel matrix.  Projective normalization
and elimination of \(\mu_1\) leave seven variables.  The coincident-critical-
point and coincident-critical-value divisors must be treated separately;
(3.5) is not a boundary classification.  This generic seven-variable chart
and its exceptional divisors are the remaining rank-two component-
extraction frontier.  Section 4 closes one complete parameter fibre and
one doubly exceptional elimination branch inside this chart.

## 4. Two exact exclusions on the generic Hurwitz chart

Dehomogenize (3.5), normalize the first channel coefficient, and write

\[
\begin{aligned}
F(x,y)&=(1+x)B(y)+x^2(\lambda+x)D(y),\\
B(y)&=1+a_1y+a_2y^2+a_3y^3,\\
D(y)&=b_0+b_1y+b_2y^2+d_3y^3.
\end{aligned}                                             \tag{4.1}
\]

The first moment is linear with nonzero constant pivot and gives

\[
d_3=-1-\frac{a_1}{3}-\frac{\lambda b_2}{3}.              \tag{4.2}
\]

The calculation below is localized at the nonzero quadratic
discriminant

\[
\begin{aligned}
\Delta={}&(9+2a_1+\lambda b_2)^2\\
 &+(3+2\lambda b_1+3b_2)
 \bigl((3-\lambda)a_1+2a_2-3\lambda-\lambda^2b_2\bigr)
\end{aligned}                                             \tag{4.3}
\]

and at the channel minor

\[
M_{01}=b_1-a_1b_0.                                       \tag{4.4}
\]

These two localizations are part of the statement: neither their zero
divisors nor the other projective channel charts are silently included.

There are two characteristic-zero unit calculations.

> **Proposition 4.1 (the complete \(\lambda=0\) fibre).**  In
> \(\mathbb Q[a_1,a_2,a_3,b_0,b_1,b_2,\Delta^{-1},M_{01}^{-1}]\), the
> moments \(\mu_2,\ldots,\mu_8\) at \(\lambda=0\) generate the unit
> ideal.

Thus the entire declared \(\lambda=0\) fibre of the generic Hurwitz
chart contains no all-order pure-moment point.

The coefficient of \(a_3\) in the primitive second moment is

\[
P_1=3(3b_0\lambda+3b_0+4b_1).                            \tag{4.5}
\]

On \(P_1=0\), substitute

\[
b_1=-\frac34b_0(\lambda+1).                              \tag{4.6}
\]

The remaining second moment is linear in \(a_2\).  Its primitive pivot
is

\[
P_2=-\frac34
 \bigl(9b_0\lambda^2+2b_0\lambda+9b_0-16b_2-16\bigr).
\tag{4.7}
\]

On the secondary boundary \(P_2=0\), substitute

\[
b_2=\frac{9b_0\lambda^2+2b_0\lambda+9b_0-16}{16}.        \tag{4.8}
\]

Exact rational elimination after (4.6)--(4.8), retaining the equation
\(\mu_2=0\) and the open condition \(\Delta M_{01}\ne0\), gives

\[
\bigl(\mu_2,\mu_3,\ldots,\mu_8\bigr)=(1).                \tag{4.9}
\]

> **Proposition 4.2 (doubly exceptional pivot branch).**  The locus
> \(P_1=P_2=0\) in the declared generic Hurwitz chart has no common zero
> of \(\mu_2,\ldots,\mu_8\), and hence no all-order pure-moment point.

Both calculations construct the moments of \(3F\) over the
\(521\)-bit Mersenne prime \(2^{521}-1\).  Every integer coefficient of
the order-\(m\) contraction is bounded in absolute value by
\((3m)!\,52^m\), and the prime is larger than twice this bound through
order eight.  Symmetric representatives therefore recover the exact
integer polynomials, rather than a heuristic modular image.  Scaling
\(F\) by \(3\) does not change their zero sets.  `msolve` over
\(\mathbb Q\) returns the one-element basis \([-1]\) in both cases.

The square system \(\mu_2,\ldots,\mu_8\) in the original seven
parameters has mixed volume \(74\,144\).  This explains the size of the
unreduced elimination, but is only a numerical-algebraic complexity
measurement.  Generic finite-field elimination, native saturation, and
the reduced branches below exceeded the recorded memory limits; those
failures have no mathematical force.

On the tertiary open inside \(P_1=0,\ P_2\ne0\), eliminating \(a_2\)
and then \(a_3\) reveals exact powers \(P_2^2,P_2^3,P_2^3,P_2^4\) in
the cleared moments \(\mu_4,\ldots,\mu_7\), respectively.  Cancelling
them leaves \(4\,493,11\,332,17\,566,40\,418\) terms over
\(\mathbb F_{1\,073\,741\,827}\).  Even this four-parameter square
system produced no basis under either `msolve`'s exact-sparse or
probabilistic-sparse kernel with a 5 GB address-space limit; the native
saturation kernel also terminated without a basis.  These figures record
the current computational ceiling only.  In particular, they are neither
a finite-field exclusion nor evidence that the branch is nonempty.

The surviving pieces of this one channel chart are now explicit:

1. \(P_1\ne0\), where \(\mu_2\) eliminates \(a_3\);
2. \(P_1=0,\ P_2\ne0\), where \(\mu_2\) eliminates \(a_2\); and
3. the zero divisors of \(\Delta M_{01}\).

On the first piece, the reduced \(\mu_3\) is quadratic in \(a_2\); on
the second, the reduced \(\mu_3\) is linear in \(a_3\).  These are exact
triangular reductions, but their remaining unspecialized ideals have
not been proved unit.  The projective \(B(0)=0\) chart, the other
channel-minor opens, and the exceptional cubic-pencil divisors also
remain separate.  Propositions 4.1--4.2 are therefore genuine
characteristic-zero exclusions, not a global rank-two classification.

## 5. Consequence for the Picard--Fuchs work

The fixed integral fiber used in the existing order-eight/order-fourteen
calculation has raw moments

\[
2502,\quad120004752,\quad24426240171840,\quad
13040211749112437760.                                    \tag{5.1}
\]

It is therefore a generic benchmark, not a point of any relevant
finite-prefix zero fiber.  Its exact reconstructed operators, Ore division,
forward-coefficient audit, and fifty-value bridge remain valuable regression
tests.  They do not settle the corrected component locus.

The characteristic-zero gap is specifically the relative boundary identity.
Closed-path Griffiths--Dwork reduction produces Picard--Fuchs operators for
closed periods, while relative-period theory obtains an inhomogeneous
equation from the exact boundary form.  Reduction-based creative telescoping
can intentionally avoid certificates, which is why an exact operator alone
does not identify the interval endpoint trace.  See
[Lairez](https://arxiv.org/abs/1404.5069),
[Bostan--Lairez--Salvy](https://arxiv.org/abs/1301.4313), and
[Li--Lian--Yau](https://arxiv.org/abs/0910.4215).  The recent twisted-period
extension likewise formulates the problem through a Griffiths--Dwork
\(D\)-module
[on relative twisted forms](https://arxiv.org/abs/2604.09129).

For the compact chart, both interval boundaries \(y=0\) and \(y=\infty\)
must be retained.  The next invariant proof object is the
\(14_{\rm interior}+2_{0}+2_{1}\) relative Gauss--Manin connection.  It
should be constructed only after an exact corrected-system component has
survived the missing-moment audit.  At an algebraic survivor the natural
coefficient field is its residue field \(K\); descent to \(\mathbb Q\)
requires a separate argument and restriction of scalars may increase
scalar order.
