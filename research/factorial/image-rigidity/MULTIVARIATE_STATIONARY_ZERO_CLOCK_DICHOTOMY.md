# Stationary Gaussian contacts with zero clock response

Fix an ordinary finite two-monomial core and adjacent orders \(n,N=n+1\),
\(n\ge2\). Take a simple physical reference contact at positive Gaussian
strength \(W_*\), signed bias \(\lambda_*\), and integer tail gap \(h\).
Use the [signed jet notation](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md).
The phase difference has expansion
\(g(W_*+\delta)=a\delta^{2r}+O(\delta^{2r+1})\), where \(a\ne0\)
and \(r\ge1\). Write \(A(W)=\Im\Delta\kappa(W)\) and
\(B(W)=\Im\Delta\gamma(W)\) for the signed cubic and clock responses.
Assume \(B(W_*)=0\).

The earlier [curvature triangle](MULTIVARIATE_STATIONARY_BIAS_CURVATURE_TRIANGLE.md)
uses the phase divided by the inverse-height moment \(Z\). When the
clock response vanishes, its limiting affine expression is zero at
the triangle's origin. A further uniform exact expansion, at the
scale of the quartic moment, yields a different compact optimization
problem. Its ordinary reciprocal gaps remain discrete.

## The refined phase

For the unequal rows of an ordinary shared tail let \(m_i\) be the
midpoint, \(d_i\) the nonzero signed integer gap, and put
\[
 x_i=d_i/m_i,\quad w_i=d_i^2/m_i,\quad
 \eta=\max_i|x_i|,\quad W=\sum_iw_i,\quad \lambda=\sum_ix_i .
\]
Use \(K=\sum_iw_i x_i,\ Z=\sum_i x_i^2\), and \(B_4=\sum_iw_i x_i^2\).
Two additional sums are
\[
                    V=\sum_i x_i^3,\qquad U=\sum_i x_i^2/m_i .
\]
The normalized variables of a nonempty tail are
\[
 (P,Q,\mu,\nu)
       =\left(\frac{K}{\sqrt{B_4}},\frac{Z}{\sqrt{B_4}},
              \frac{V}{B_4},\frac{U}{B_4}\right).                 \tag{1}
\]
The sign of \(P\) is relevant.

At either order \(j\), the zero-endpoint tail log multiplier on its
coefficient grid is
\[
  W F_j+K H_j+Z G_j+B_4 I_j
             -\frac{V}{j}H_j+\frac{U}{6j^2}F_j+O_j(\eta B_4).   \tag{2}
\]
Here \(F_j=-j\phi,\ H_j=j\Psi,\ G_j=-F_j/(2j)\), and \(I_j=j\Xi\)
where the argument of the zero-endpoint shapes is \(k/j\).
To verify (2), with \(t=s-1/2\), the \(i\)-th contribution to the
second derivative of the log symbol per degree is
\[
 j d_i^2\psi_1(1+j(m_i+t d_i)).
\]
Writing \(y=m_i+t d_i\), the first three terms are
\(d_i^2(y^{-1}-(2j)^{-1}y^{-2}+(6j^2)^{-1}y^{-3})\).
The trigamma Laplace expansion has error bounded by
\(C_jd_i^2/y^4\) on \(y>0\), as follows by integrating the bounded
Taylor remainder in \(u/(e^u-1)=1-u/2+u^2/12+O(u^3)\)
as in the [signed center estimate](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md#a-matched-expansion-with-a-second-order-ordinary-remainder).
Expand the three displayed powers of \(y\) through degrees \(2,1,0\)
in \(t x_i\), respectively. For \(\eta\le1/2\) their sum over the
tail is
\[
            W-tK+t^2B_4-Z/(2j)+tV/j+U/(6j^2)+O_j(\eta B_4).
\]
The bound follows term by term from
\(\sum_i|x_i|^3\le B_4,\ U\le B_4\), and \(1/m_i\le |x_i|\le\eta\).
Integrating with zero endpoints uses
\(\phi''=-1,\ \Psi''=-t,\ \Xi''=t^2\) and gives (2).

Let \(C=\Im\Delta\sigma\) be the quartic response and
\(D_{11}=\Im\Delta\beta^H\) the second cubic response, evaluated at
\(W_*\). The additional coefficients at \(W_*\) are
\[
 \begin{aligned}
 E&=\Im\Delta(-\kappa_j/j),&
 F&=\Im\Delta(\alpha_j/(6j^2)),\\
 D_{12}&=\Im\Delta(-\kappa'_j/(2j)),&
 D_{22}&=\Im\Delta(\alpha'_j/(8j^2)).
 \end{aligned}                                                    \tag{3}
\]
The last two formulas follow from the exact clock shift
\(W\mapsto W-Z/(2j)\) in the presence of the factor
\(\exp(KH_j)\). In particular \(D_{12}\) multiplies \(KZ\),
and \(D_{22}\) multiplies \(Z^2\). Here \(F=0\), since
\(g'(W_*)=0\) and
\(B(W_*)=\Im(\alpha_n/(2n)-\alpha_N/(2N))=0\)
force both individual angular velocities to vanish.

The angular expansion keeps \(g(W)\), \(A(W)K\) and \(B(W)Z\),
and its other terms through scale \(B_4\) are
\[
 CB_4+EV+FU+D_{11}K^2+D_{12}KZ+D_{22}Z^2+o(B_4)                \tag{4}
\]
as \((W,\eta)\to(W_*,0)\), uniformly in the tail. For the proof,
first keep the coefficients in (3) and (4) as functions of \(W\)
on a fixed complex disk of simple log-root branches. Cauchy–Schwarz
gives \(K^2\le WB_4\) and \(Z^2\le WU\le WB_4\).
The coefficient vector's linear perturbation starting with \(KH_j+ZG_j\)
is \(O_W(\sqrt{B_4})\). Every omitted nonlinear term is
\(O_W(B_4^{3/2})=O_W(\eta B_4)\), since \(B_4\le W\eta^2\).
The error \(O_W(\eta B_4)\) with \(W\)-dependent coefficients
holds with one \(W\) derivative by Cauchy's estimate on a smaller
disk. Replacing these coefficients by their values at \(W_*\)
costs \(O_W(|W-W_*|B_4)\) in the function. Its derivative
costs \(O_W(B_4)\). This also supplies the derivative error
after rescaling around \(W_*\).

## The exact compact moment set

For \(q\in\{-1,0,1\}\) use three variables \(A_q,B_q,C_q\)
with
\[
   A_q,C_q\ge0,\qquad q B_q\ge0,\qquad B_q^2\le A_q C_q.
\]
For \(W_0>0,\ |\lambda_0|\le W_0\), let
\(\mathcal C(W_0,\lambda_0)\) be the set of \((P,Q,\mu)\) such that
\[
 \sum_q A_q=W_0,\quad \sum_q q A_q=\lambda_0,\quad
 \sum_q C_q=1,\qquad
 P=\sum_q B_q,\quad Q=\sum_q q B_q,\quad \mu=\sum_q q C_q.          \tag{5}
\]
This is a compact semialgebraic feasibility problem in nine
variables.

The set (5) is exactly the set of limiting \((P,Q,\mu)\) along
ordinary tails with \((W,\lambda,\eta)\to(W_0,\lambda_0,0)\)
and any fixed integer signed gap. For necessity put
\(q_i=1/d_i\) and \(z_i=x_i/\sqrt{B_4}\).
The \(i\)-th row contributes a rank-one triple
\((A,B,C)=(w_i,w_i z_i,w_i z_i^2)\) at \(q_i\).
The six totals in (5) sum to the required values because
\(q_i w_i=x_i\). Every total is affine in \(q_i\).
Interpolate each triple between \(q=0\) and
\(q=\operatorname{sgn}(q_i)\) using weights \(1-|q_i|,|q_i|\).
Both new triples have the required signs. At a fixed \(q\),
sums of the associated positive semidefinite two-by-two
matrices remain positive semidefinite. This proves the finite
constraints for every ordinary profile and for all limits.

For each of the three \(q\), its feasible triple is a sum of at
most two rank-one triples. If \(B_q\ne0\), use
\((B_q^2/C_q,B_q,C_q)\) and \((A_q-B_q^2/C_q,0,0)\).
If \(B_q=0\), use \((A_q,0,0)\) and \((0,0,C_q)\).
Omit any zero triple. The sign condition persists.
Assign one group of ordinary rows to each resulting triple.
Write \(A,B,C\) for its coordinates. When \(A>0\), put
\(z=B/A\), so \(C=Az^2\).
To obtain the limits use a common positive parameter
\(\varepsilon\to0\) and make the combined \(B_4/\varepsilon^2\)
tend to one.

- If \(A>0\) and \(q\ne0\), take the fixed signed gap
  \(d=1/q\). For \(z\ne0\) let \(x\) be asymptotic to
  \(\varepsilon z\), using a midpoint nearest \(d/(\varepsilon z)\)
  in \(\mathbb Z+d/2\). If \(z=0\), use target
  \(x=\operatorname{sgn}(d)\varepsilon^2\) in the same rounding.
  Repeat this row \(\lfloor A/(dx)\rfloor\) times, so the
  group curvature tends to \(A\).
- If \(A>0\) and \(q=0\), use one row. When \(z\ne0\)
  choose signed gap \(d\) nearest \(A/(\varepsilon z)\) and
  midpoint nearest \(d/(\varepsilon z)\) in the same lattice.
  When \(z=0\) use the positive target \(x=\varepsilon^2\)
  in place of \(\varepsilon z\) in both choices.
  Then \(1/d\to0\), \(dx\to A\), and \(x/\varepsilon\to z\).
- If \(A=0\), then \(C>0\) and \(B=0\). For \(q\ne0\)
  take one row with \(d=1/q\). For \(q=0\) take one row
  with \(d=\lceil\varepsilon^{-1/2}\rceil\). Use target
  \(x=\operatorname{sgn}(d)(C\varepsilon^2/|d|)^{1/3}\)
  and the midpoint rounding above. Now \(dx\to0\),
  \(dx^2/\varepsilon\to0\), and
  \(dx^3/\varepsilon^2\to C\).

Each midpoint tends to infinity with \(|d|/m\to0\), so both
endpoints are positive integers. The displayed group limits give
every required coordinate, first with \(\varepsilon\) in place of
\(\sqrt{B_4}\). Their third coordinates sum to one, allowing
that substitution. To fix the signed gap, let \(D\) be the
remaining integer correction. Append \(|D|\) unit rows of
the correct sign at a smaller exponent exceeding
\((1+|D|)\varepsilon^{-4}\). Their total curvature is
\(O(\varepsilon^4)\) and their relative gaps are \(O(\varepsilon^4)\),
so none of the limits change.

## Exclusion and transverse unfolding

Assume both \(A(W_*+\delta)\) and \(B(W_*+\delta)\) vanish
through order \(r-1\), and let \(b,c\) be their coefficients
of order \(r\), respectively. For \(r=1\) the assumptions
are just vanishing of the two responses at the contact.
Define the quadratic
\[
 \begin{aligned}
 f(Y,P,Q,\mu)=\operatorname{sgn}(a)\bigl[
     aY^2+(bP+cQ)Y+C+E\mu
       +D_{11}P^2+D_{12}PQ+D_{22}Q^2\bigr].
 \end{aligned}                                                     \tag{6}
\]
Minimize it on \((P,Q,\mu)\in\mathcal C(W_*,\lambda_*)\)
and \(Y\in\mathbb R\) if \(r\) is odd, \(Y\ge0\) if \(r\) is even.
The minimum exists because the moment set is compact and the
coefficient of \(Y^2\) is strictly positive.

**Theorem.** If the minimum in (6) is strictly positive,
there are a neighborhood of the reference and a relative-gap
tolerance on which no ordinary shared tail can give an exact
angular contact of the indicated branches. The conclusion is
uniform in the tail's row count, gaps and heights. It includes
saturated bias. If the minimum is strictly negative and
\(|\lambda_*|<W_*\), the reference is a limit of simple transverse
interior Gaussian contacts of finite ordinary cores. The adjacent
orders and integer total gap are preserved, and every strict
curvature cap persists. The two strict tests concern this branch
pair at the specified positive contact.

For exclusion, the angular equation, the response vanishing
orders, and \(K,Z=O_W(\sqrt{B_4})\) first give
\[
              |\delta|^{2r}
                      =O_W(\sqrt{B_4}\,|\delta|^r+B_4).
\]
Thus \(Y=\delta^r/\sqrt{B_4}\) is bounded at any candidate
contact. Divide by \(B_4\) and use (4) to get \(0=f+o(1)\).
The compact moment set contains every limiting normalized
vector. A sequence approaching the reference through arbitrarily
diffuse angular contacts would contradict the positive minimum.

For unfolding, choose a point with \(f<0\).
In each rank-one triple replace \((A,B,C)\) by
\((tA,\sqrt t B,C)\), with \(0<t<1\) close enough to one to
preserve the strict negative value.
The new moment vector has \(W_0=tW_*\), \(\lambda_0=t\lambda_*\),
\(P,Q\) multiplied by \(\sqrt t\), and \(\mu\) unchanged.
Realize it by exact ordinary anchors. Keep these anchors exact
and use remaining Gaussian strength \(T=W-W_{\rm A}\).
At \(W=W_*+B_{4,{\rm A}}^{1/(2r)}X\), the signed angular
difference divided by \(B_{4,{\rm A}}\) converges with one
real derivative to the polynomial (6) at \(Y=X^r\) for the
scaled moments. This follows from (2)–(4) and the Taylor
expansions of \(g,A,B\). Its quadratic in \(Y\) is negative
at an admissible point, so it has a simple nonzero zero in
the admissible range. Composing with \(Y=X^r\) keeps this
zero simple. The actual angular gap has a nearby simple zero.

Take the remaining integer gap to be the original \(h\) minus
the signed gap of the anchors. Solve the radial equality for
the remaining bias. The [exact center correction](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md#a-matched-expansion-with-a-second-order-ordinary-remainder)
for the anchors is \(-\lambda_{\rm A}/(2nN)+O(Z_{\rm A})\)
after stripping their integer gap. The normalized root
moduli therefore force \(\lambda_{\rm rem}\to(1-t)\lambda_*\).
The remaining strength tends to \((1-t)W_*\), which is
strictly greater than the limiting absolute remaining bias.
The root scales and the combined coefficients converge, including
when the anchor signed gap grows. The larger-exponent curvature
of the anchors differs from their midpoint curvature by
\(O(\eta W_{\rm A})\), which gives the cap assertion.

The strict conditions admit finite certificates by [real quantifier elimination](https://arxiv.org/pdf/1609.02879).
Use (5) and rational outer intervals shrinking to \(W_*,\lambda_*\)
and to the coefficients of (6). Bound \(Y\) outside a fixed
compact interval by its positive quadratic term. Compactness and
the strict margin ensure that a finite real algebraic positivity
test succeeds. The root and coefficient bounds in (2)–(4) and
this certificate give a computable local tolerance.
For strict negative minimum at interior bias, first scale the
moment vector as in the proof, leaving strict slack in the
remaining bias bound. Approximate its finitely many rank-one
triples and \(X\) by algebraic data, preserving the signs and
normalizing the third coordinate. The two strict inequalities
persist. Enumerating such data with enclosures for the reference
constants therefore certifies strict negativity whenever it occurs.

For a finite ordinary core the minimum is nonzero
when \(\lambda_*/W_*\) is rational and
\(\operatorname{sgn}(a)C\ne |E|\).
Indeed put \(\ell=\lambda_*/W_*\) and use the semialgebraic set
\(\mathcal C(1,\ell)\). Its \(P,Q\) scale by \(\sqrt W\)
in \(\mathcal C(W,\ell W)\), while \(\mu\) stays unchanged.
Write \(Y=\sqrt W Y_0\). The expression minimized in (6) is
then \(g_0+W R_0\), where
\(g_0=\operatorname{sgn}(a)(C+E\mu)\) and \(R_0\) is a fixed
quadratic in \(Y_0,P,Q\) at mass one. If
\(\operatorname{sgn}(a)C<|E|\), take \(Y_0=P=Q=0\) and
all the third-coordinate weight at the suitable \(q=\pm1\).
This gives a negative minimum at every \(W>0\).
Otherwise \(g_0>0\) uniformly. If the minimum vanishes at \(W_c>0\),
the minimizing point at \(W_c\) is strictly negative for
\(W>W_c\). For \(0<W<W_c\), every value is a convex combination
of its value at \(W_c\) and a strictly positive value at zero,
so the minimum is positive. Thus \(W_c\) is an isolated zero.
It is real algebraic by quantifier elimination, since \(\ell\)
is rational and the response coefficients are algebraic at the
exceptional clock. The actual \(W_*=-2nN\log q_*\) is
transcendental by Hermite–Lindemann. Therefore it is not \(W_c\).
The sign searches apply at the actual contact. At interior bias
this gives the exclusion or unfolding conclusion. At saturated
bias only a positive minimum gives the stated exclusion.

## A lower-order clock response

If the order \(s\) of \(B(W_*+\delta)\) satisfies \(1\le s<r\),
write \(c_s\) for its nonzero leading coefficient. At an
interior contact there is a further transverse unfolding test
when \(s\) is odd or \(a c_s<0\). Use reflected pairs of unit-gap
rows with a common midpoint tending to infinity and total
curvature tending to any fixed
\(0<L<W_*-|\lambda_*|\). Their signed gap and bias vanish.
They have \(K=V=0\), \(B_4,U=O(Z^2)\).
At \(W-W_*=Z^{1/(2r-s)}X\), divide (4) by
\(Z^{2r/(2r-s)}\). The limit with one derivative is
\(aX^{2r}+c_sX^s\), since \(2r/(2r-s)<2\).
The sign condition gives a nonzero simple real root of this
polynomial. The physical scale is recovered as before, with
limiting remaining strength \(W_*-L>|\lambda_*|\).

If instead \(s\) is even and \(a c_s>0\), the clock term in
the phase has the same sign as \(g\) near the contact.
Under the same order-\(r\) vanishing hypothesis on \(A\),
the strict positive certificate (6) with \(c=0\) still excludes
ordinary angular contact. Indeed the sign of \(B(W)Z\) can
only increase the signed phase. The angle equation first
bounds \(\delta^r/\sqrt{B_4}\) using the other terms, and then
(6) with \(c=0\) supplies the positive margin after division
by \(B_4\). A negative minimum of this latter test has no
corresponding general unfolding assertion.

The zero minima and the low-order clock cases not covered by
these strict tests remain open. The negative branch produces
Gaussian contacts with exact ordinary anchors. Exclusion of exact
ordinary completions near the resulting transverse contacts is a
separate arithmetic question.
