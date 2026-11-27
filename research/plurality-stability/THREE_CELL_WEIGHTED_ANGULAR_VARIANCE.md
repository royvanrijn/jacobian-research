# Weighted variance penalties for the sharp circular arc comparison

The sharp arc slack absorbs a weighted radial variance whenever the
weight is a positive semidefinite kernel with two pointwise upper
bounds. In a one-mode comparison, the same argument allows a larger
weight at the cost of a loss involving only the median label of the
trial direction. That loss is small when the profile remains close to
the trial line. The actual first angular moment deficit is available in
both comparisons.

## Weighted cancellation of the dominant-label terms

Use the [sharp arc coefficients](THREE_CELL_SHARP_POSITIVE_ARC_ENVELOPE.md)
and [variance absorption notation](THREE_CELL_SHARP_VARIANCE_ABSORPTION.md).
Let $p=p(t)$ and $q=p(s)$ be two radial mass profiles in the
three-label simplex, $D=|p-q|^2$, $F(D)=D(1-D/2)$ and
$u_i(t)=(p_i(t)-1/2)_+$. Define
\[
 Q(t,s)=d^{-1}e^{-\kappa(t+s)},\qquad a=\gamma\sqrt{ts},\qquad
 n(a)=\frac{M(a)-I_0(a)}4,\qquad b(a)=2(P_\sharp(a)-N(a)),
\]
where $N(a)=2\int_0^1 v e^{a\cos(\pi v)}dv$ as in (SV3).
Write \(\mathscr L_a(p,q)\) for the difference between the sharp
scalar upper expression and the aligned-arc pairing. The positive
and negative Taylor bounds in (SV3) and (SR2) give
\[
 \mathscr L_a(p,q)\ge b(a)\sum_i u_i(t)u_i(s)+n(a)F(D).        \tag{WV1}
\]

Let $W(t,s)$ be a symmetric nonnegative integrable kernel such that
its quadratic form on bounded functions in $L^2(\nu)$ is positive
semidefinite. Suppose $L>2$ and
\[
    W\le Qn(a),\qquad LW\le Q b(a)\quad\text{pointwise}.        \tag{WV2}
\]
Then for each $j\in\{1,2,3\}$,
\[
\begin{split}
 \iint Q\mathscr L_a\,d\nu d\nu
  \ \ge\ &\frac14\iint W D\,d\nu d\nu
 -\left(\frac{8}{L-2}-L\right)_+\langle u_j,Wu_j\rangle \\
 &+\iint (Q n(a)-W)F(D)\,d\nu d\nu .                         \tag{WV3}
\end{split}
\]
In particular no dominant-label error is needed for $L\ge4$.
The kernel can be a sum of localized positive product kernels. When
both inequalities hold for a function of $a$ times $Q$, a series
with nonnegative coefficients in $a$ gives the required positive
form by expanding the finite partial sums into rank-one forms and
passing to the integrable limit.

To prove the cancellation, use the pointwise estimate in (SV1),
\[
 D(D-3/2)_+\le4\sum_{i\ne j}u_i(t)u_j(s).
\]
Indeed if $D>3/2$, both profiles have dominant coordinates with
different labels. With $u=\sum_i u_i(t)$ and $v=\sum_i u_i(s)$,
their norms give $D\le1+2u^2+2v^2\le3/2+2uv$, using
$u^2+v^2-uv\le1/4$ on $[0,1/2]^2$. Also $D\le2$.
Since $F(D)=D/4-D(D-3/2)/2$, subtract the off-diagonal bound in
the part weighted by $W$. Let $x_i=\langle u_i,Wu_i\rangle^{1/2}$.
Positive semidefiniteness implies
\(\langle u_i,Wu_k\rangle\le x_i x_k\). For the two indices
$i,k\ne j$,
\[
 \begin{split}
 L\sum_i x_i^2-2\sum_{\ell\ne m}x_\ell x_m
 &\ge (L-2)(x_i^2+x_k^2)-4x_j(x_i+x_k)+Lx_j^2\\
 &\ge \big(L-8/(L-2)\big)x_j^2 .
 \end{split}
\]
The last step completes a square for each index distinct from $j$.
This proves (WV3). The first angular deficit is subtracted from the
paired angular integral before applying the aligned-arc Taylor
comparison, as in (SR4)–(SR5). The use of its aligned-arc slack in
(WV3) therefore leaves that deficit available.

## Explicit linear and cubic weights

The following functions satisfy (WV2) with $W=Qw(a)$.
\[
 \begin{array}{c|c}
 L & w(a)\\ \hline
 4 & \dfrac{24}{25\pi^2}a+\dfrac{a^3}{10000}\\[3pt]
 3 & \dfrac{a}{\pi^2}+\dfrac{7a^3}{1000}.
 \end{array}                                                    \tag{WV4}
\]
The nonnegative odd series for $n(a)$ in (SR7) starts with
\(a/\pi^2+7a^3/(54\pi^2)\), which bounds both functions above.
For the other constraint, the following three lower bounds hold for
every $a\ge0$.
\[
 P_\sharp(a)-N(a)\ \ge\
 \max\left\{
 \frac{48a}{25\pi^2}+\frac{a^3}{5000},
 \frac32\left(\frac a{\pi^2}+\frac{7a^3}{1000}\right),
 \frac75\left(\frac a{\pi^2}+\frac{a^3}{125}\right)
 \right\}.                                                       \tag{WV5}
\]
For $0\le a\le8$ expand the exponential through degree 16.
The remaining odd-even pairs are nonnegative by the signed moment
inequalities in (SV4), since each even coefficient is nonpositive
and $a/(j+1)\le1$ for the remaining odd degrees $j\ge17$.
The Fourier identity for the signed density in (SV4) yields exact
rational coefficients for the resulting polynomial after division
by $a$ and multiplication by $\pi^2$. Use $\pi^2<10$ in the
cubic term of each proposed bound. The
[rational checker](check_three_cell_weighted_angular_variance.py)
verifies positive Bernstein coefficients on all 128 subintervals
of width $1/16$, including the extension at $a=0$. In the order
displayed in (WV5), the floors of one million times the smallest
coefficients are 2634, 14054, and 66136.

For $a\ge8$, subtract the threshold $e^{-a/\sqrt2}$ from the
signed integral for $P_\sharp-N$, which has mass zero and a
single sign change at $v=3/4$. The remaining integrand is
nonnegative, and restricting to $v\in[0,1/4]$ gives
\((e^{7a/10}-e^{-7a/10})/32\). At $a=8$, this is larger than the
third bound in (WV5) by the exponential series through degree 12
and $\pi^2>9$, also checked rationally. After division by $a^3$
the exponential lower bound increases, because its logarithmic
derivative before division is at least $7/10>3/a$. The third
bound divided by $a^3$ decreases. The third bound exceeds the
other two for $a\ge8$ by direct coefficient comparison. This
proves (WV5) and (WV4).

## Correction on the radial operator and the median-label loss

For the first row of (WV4), the Dirichlet term of (WV3) subtracts
a variance in each of the two measures \(\mu_n\), $n=1,3$, of
(SR7). More generally if $w(a)=\sum_n c_n a^n$ has nonnegative
coefficients, write \(\omega_n=c_n(2\rho)^n d\,\Gamma(1+n/2)^2\).
Then
\[
 \frac14\iint W D\,d\nu d\nu=\frac12\sum_n\omega_n\operatorname{Var}_{\mu_n}(p),
 \qquad \langle u_i,Wu_i\rangle=\sum_n\omega_n(\mathbb E_{\mu_n}u_i)^2.  \tag{WV6}
\]
This follows by taking the radial factors of $a^n Q$ and normalizing
each to $\mu_n$. The first weight in (WV4) keeps the variance term
of (SV1) and adds a positive cubic-feature variance. Equivalently
the corrected radial kernel is $K+W/2$, and its row weight is
$J+\tfrac14\int W(t,s)d\nu(s)$. This form also applies to any
other $W$ satisfying (WV2) with $L\ge4$.

For the second row, the coefficient on the negative term in (WV3)
is five. Suppose a feasibility estimate uses $p-p_*=\alpha\phi+g$
and let $j$ be the median coordinate of $\alpha$. At any radial
point choose $i\ne j$ with $\alpha_i\phi\ge\alpha_j\phi$, using
the two extreme coordinates. If $p_j>1/2$, then
$p_j-p_i\ge2u_j$ since $p_i\le1-p_j$. Thus
\[
             u_j\le |g|/\sqrt2 .                                \tag{WV7}
\]
One can improve this bound by the size of the trial. A centered
three-coordinate vector has $|\alpha_j|\le|\alpha|/\sqrt6$.
two coordinates of the median's sign have magnitude at least its
magnitude, and the remaining coordinate has at least twice that
magnitude. For $0<\delta<1/6$ and
$E=\{|\alpha\phi|>\sqrt6(1/6-\delta)\}$, the predicted median
mass is at most $1/2-\delta$ on the complement of $E$.
Writing $r_*=\sqrt{2/3}$ and using $|g_j|\le r_*|g|$, we get
$u_j\le(r_*|g|-\delta)_+\le r_*^2 |g|^2/(4\delta)$ there.
If $G^2=\int J_*|g|^2d\nu$ with the positive trial row weight
$J_*$ and $d\mu_n\le B_n J_*d\nu$, then
\[
 \mathbb E_{\mu_n}u_j\le \frac{r_*^2 B_n G^2}{4\delta}
                       +\sqrt{\frac{B_n\mu_n(E)}2}\,G .        \tag{WV8}
\]
The last term uses (WV7) on $E$. The nonnegative last line of
(WV3) further retains at least the purity-times-variance reserve
of (SR7) in $\mu_3$, with coefficient
\((7/(54\pi^2)-7/1000)(2\rho)^3d\,\Gamma(5/2)^2\).
This median-label comparison and the unrestricted comparison with
$L=4$ are two lower bounds on the same slack. Either can be used
at each feasible amplitude. An all-correlation norm or amplitude
bound closing the remaining gap is still needed.
