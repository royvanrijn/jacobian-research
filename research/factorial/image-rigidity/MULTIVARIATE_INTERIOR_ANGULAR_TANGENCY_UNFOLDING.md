# Ordinary rows unfold stationary interior Gaussian contacts

Use the [compactified interior reference](MULTIVARIATE_COMPACTIFIED_INTERIOR_CONTACT_FINITE_REDUCTION.md)
at adjacent orders \(n,N=n+1\), with finite or countable retained
ordinary core \(R\) of summable curvature and total signed degree
difference \(\Delta_0\). Suppose there is a physical contact at
strength \(T_*>0\) and bias \(|\lambda_*|<T_*\), simple in both
polynomials. The unscaled endpoint-normalized root branches
\(z_j(T)\) are nonreal near \(T_*\). Choose their phase difference
\(\theta(T)=\arg z_N(T)-\arg z_n(T)\) with value zero at the
contact.

Assume \(\theta(T_*+x)=c x^p+O(x^{p+1})\), with \(c\ne0\) and
\(p\) even. The order is finite by the angular quotient
resultant at the endpoint-binomial limit. Odd \(p\) is already
handled by ordinary finite truncation. Put
\[
 \phi(s)=\frac{s(1-s)}{2},\qquad
 \Psi(s)=\frac{s(1-s)(2s-1)}{12},\qquad
 \Xi(s)=\frac{(s-1/2)^4-1/16}{12}.
\]
As in the compactified reference, let \(\mathcal L_T(H)\)
denote the first derivative of the phase difference on
multiplying each order's coefficient at \(k\) by
\(\exp(u j H(k/j))\) at \(u=0\). Write
\[
   \mathcal L_{T_*+x}(\Psi)=b x^s+O(x^{s+1})
\]
when it has finite order \(s\), with \(b\ne0\), and take \(s=\infty\)
for the zero germ. Let
\[
 \begin{split}
  C_2&=\tfrac12\left(\frac{(\arg z_n)'(T_*)}{n}
                        -\frac{(\arg z_N)'(T_*)}{N}\right),\\
  C_4&=\mathcal L_{T_*}(\Xi),\\
  H_2&=\tfrac12\left.\partial_u^2
                      \theta_{\Psi}(T_*,u)\right|_{u=0}.
 \end{split}
\]
Here \(\theta_{\Psi}(T,u)\) is the phase difference of the
polynomials with the same \(\exp(u j\Psi(k/j))\) multiplier.

**Theorem.** Each of the following strict tests gives a sequence
of finite ordinary cores with simple transverse physical
interior contacts. The consecutive orders and the integer
\(\Delta_0\) are the original ones.

1. \(s<p/2\). The finite contact strengths tend to \(T_*\).
   This includes stationary cubic response at \(T_*\) if its
   order of vanishing is low enough.
2. For some real \(0<w<T_*-|\lambda_*|\), put
   \(C(w)=C_2/w^2+C_4/w+H_2\). If \(s=p/2\), require
   \(b^2>4c C(w)\). If \(s>p/2\), require \(c C(w)<0\).
   The new contact strengths tend to \(T_*-w\).
3. For some real \(0<2w<T_*-|\lambda_*|\), require
   \(c(C_2+w C_4)<0\). There is no hypothesis on \(s\).
   The new strengths tend to \(T_*-2w\).

In every case the normalized reference pairs, relative physical
scales and signed biases converge to their original values.
The added ordinary rows in cases 2 and 3 have summed limiting
curvature \(w\) and \(2w\), respectively. Their variance replaces
an equal amount of Gaussian strength. If the stated sign
condition holds for arbitrarily small \(w\), one can also
choose the strengths to tend to \(T_*\). The total curvature
cap is preserved whenever the original cap has strict slack.

## First-order response along a higher-order contact

If a row has midpoint \(m\), signed gap \(d\), midpoint
curvature \(w_d=d^2/m\), \(A_d=d^2/m^2\), and
\(\eta_d=|d|/m\), its normalized log contribution after replacing
\(w_d\) of Gaussian strength is \(jR_j(k/j)\), where
\[
 R_j(s)=A_d d\,\Psi(s)
         +\frac{A_d}{2j}\phi(s)+A_d w_d\Xi(s)
         +O_j\bigl(A_d(w_d\eta_d+\eta_d+m^{-1})\bigr).          \tag{1}
\]
This is the [exact compensated chord expansion](MULTIVARIATE_COMPACTIFIED_INTERIOR_CONTACT_FINITE_REDUCTION.md#ordinary-row-perturbations-at-fixed-effective-strength).
It follows by expanding the trigamma function in the second
derivative of the row chord and integrating with zero endpoint
values. The error is uniform at all the sample points of the
two fixed orders. The Gaussian comparison parameter \(T\)
in this proof is the effective strength before subtracting
the compensating variance. The added log perturbations are
independent of \(T\).

For case 1 put \(\nu=p/(p-s)\in[1,2)\). Choose a positive
integer \(L\) such that \(2L+2>(L+2)\nu\), and let
\[
 m=t^{2L+1},\quad d=\pm2t^L,\quad w_d=4/t,\quad
                         u=8t^{-(L+2)}
\]
at large integers \(t\). Both endpoints \(m\pm d/2\) are
ordinary positive integers, and (1) becomes
\(R_j=\pm u\Psi+o(u^\nu)\), because the largest displayed
error after the cubic term is \(O(t^{-(2L+2)})\).
The squared norm of the whole log perturbation is also
\(o(u^\nu)\).

Choose a finite truncation with all coefficient and center
errors \(o(u^\nu)\). Such truncations exist by summability,
and the implicit function theorem for simple root branches
makes the first-order phase response uniform in \(T\).
It also gives uniform bounds for any fixed number of \(T\)
derivatives of the Taylor remainders in the log coefficients.
At \(T=T_*+u^{1/(p-s)} y\), the new phase difference divided
by \(u^\nu\) converges on bounded \(y\) intervals to
\[
                              c y^p\pm b y^s .                \tag{2}
\]
The convergence holds with its first \(y\) derivative.
Choose the row orientation so the nonzero \(c,b\) terms in
(2) have opposite signs. Then (2) has a positive simple root,
because its nonzero factor is a binomial in \(y^{p-s}\).
The perturbed phase difference has a nearby simple real zero.
Its physical Gaussian strength is \(T-w_d\to T_*\).
This proof uses ordinary rows even when the original
\(\mathcal L_{T_*}(\Psi)\) vanishes.

## Second-order and reflected adjustments

In case 2 take \(d=\pm 2t\) and choose an integer midpoint
\(m\) nearest to \(d^2/w\). Then \(w_d\to w\),
\(\eta_d\to0\), and \(u=A_d|d|\to0\). The identities
\[
                 \frac{A_d}{u^2}=\frac1{w_d^2},\qquad
                 \frac{A_d w_d}{u^2}=\frac1{w_d}
\]
show which terms of (1) survive at order \(u^2\). The
remainder in (1) is \(o(u^2)\). On finite truncations of
accuracy \(o(u^2)\), Taylor expansion of the analytic
phase response therefore gives the perturbed difference
\[
 \theta(T)\ \pm u\,\mathcal L_T(\Psi)
        +u^2\bigl(C(w)+o(1)\bigr)
\]
as \(T\to T_*\). The constant \(H_2\) is the quadratic
response to the leading log perturbation. The clock term
has response \(C_2/w^2\), since a log direction
\(\phi/(2j)\) in order \(j\) equals a clock change of
\(-1/(2j)\) in that order. As before, expansion of the
analytic root logarithm also gives derivative convergence
after the rescaling \(T=T_*+u^{2/p}y\). The limiting phase
polynomial after division by \(u^2\) is
\[
 \begin{cases}
   c y^p\ \pm b y^{p/2}+ C(w),&s=p/2,\\
   c y^p+C(w),&s>p/2.
 \end{cases}
\]
In the first line the discriminant in \(y^{p/2}\) is
strictly positive. At least one real zero in that variable
is nonzero. By choosing the sign of \(b\), one can make
a nonzero simple zero in \(y^{p/2}\) positive, which gives
a positive simple zero in \(y\). The sign condition in the
second line gives a positive simple zero directly. Simple
zeros persist by the uniform first derivative and opposite
endpoint signs on a small fixed interval around the zero.

For case 3 add a reflected pair of the same row as in
case 2, compensating both variances. The leading cubic
terms cancel exactly. The summed perturbation has scale
\(2A_d\) and, after division by that scale, tends to
\(\phi/(2j)+w\Xi\) at order \(j\). Its squared norm is
\(o(A_d)\). Use the rescaling
\(T=T_*+(2A_d)^{1/p}y\). The normalized phase difference
converges with one derivative to \(c y^p+C_2+w C_4\),
which has a positive simple zero under the stated sign.
The summed compensating strengths are \(2w_d\to2w\).

## Physical scale and exact sign tests

In each construction choose \(\Delta_0\) as the unchanged
total integer degree difference, with the Gaussian part's
integer gap adjusted relative to the finite retained core.
The center correction for a new row tends to zero, including
when \(w_d\) tends to \(w>0\). In fact the
[exact row center estimate](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#uniform-estimates-for-rows-of-large-size)
writes it as \(-\log(b/a)/(2nN)\) with an error bounded by
\(|b-a|/(4n^3\min(a,b)^2)\), and both terms tend to zero
since \(|d|/m\to0\). The correction from the retained
truncation converges by summability. Thus \(\Gamma_{\rm fin}\),
the retained center at fixed \(\Delta_0\), tends to \(\Gamma\).

At every constructed angular zero set
\[
  \lambda_{\rm fin}=2nN\left(\Gamma_{\rm fin}
              -\log |z_{N,\rm fin}/z_{n,\rm fin}|\right).
\]
This makes the physical root moduli agree and the angular
zero makes the roots equal. The normalized coefficients
and roots converge by the compensated perturbation estimates.
Hence \(\lambda_{\rm fin}\to\lambda_*\). The respective strict
conditions on \(w\) preserve \(|\lambda_{\rm fin}|<T_{\rm fin}\).
Angular transversality in effective \(T\) is the same as in
the new Gaussian strength as the added rows are fixed in
each comparison. The endpoint and midpoint curvatures of
the added rows differ by \(o(1)\). This proves preservation
of any strict curvature cap.

For a finite ordinary core at an algebraic exceptional clock,
the phase derivatives used here are algebraic. They are
obtained by differentiating the simple root equation, then
taking imaginary parts of logarithmic root derivatives.
The jets for \(p,s\), when finite, can be tested with exact
algebraic coefficients. For the second-order signs one can
take a rational \(w\) in the indicated open interval.

In all three tests it suffices to inspect the cubic response
through derivative order \(p/2\), including when its germ is
zero. For a simple even contact of an infinite retained core
with unbounded gaps, the earlier [signed unfolding tests](MULTIVARIATE_COMPACTIFIED_INTERIOR_CONTACT_FINITE_REDUCTION.md#signed-unfolding)
and the tests here give explicit residual conditions. A
contact admitting none of these reductions to a finite
sign-changing interior contact has
\[
                C_2=C_4=0,\qquad s\ge p/2 .
\]
If \(s=p/2\), it must satisfy \(b^2\le4cH_2\). If
\(s>p/2\), it must satisfy \(cH_2\ge0\). Indeed the earlier
tests give \(C_2=0\) from the common angular velocity and
\(C_4=0\) from the unbounded gaps. The remaining inequalities
then follow by the first and second tests above.
