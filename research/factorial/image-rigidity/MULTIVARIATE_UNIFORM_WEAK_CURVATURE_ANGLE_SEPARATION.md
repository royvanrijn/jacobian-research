# Exact angular separation at a fixed weak-curvature cap

For ordinary distinct exponent columns \(a,b\), set
\[
    W=\sum_{a_i\ne b_i}\frac{(b_i-a_i)^2}{(a_i+b_i)/2}.
\]
There is a computable positive rational \(\eta_*\) such that for every
terminal \(N=n+1\ge2\), the two factorial moment polynomials have
disjoint root arguments if \(W\le\eta_*\).
In particular, they are coprime over every characteristic-zero field.
The cap is uniform in the exponent rows and dimension.
For fully shared profiles and \(n\ge2\), it is also enough that
\(nX\le\eta_*\), where \(X\) is the order-\(n\) Gauss reciprocal sum.
The proof specifies a finite choice of the cap and compares roots across orders.

The analytic estimates are those of the
[effective vanishing-curvature comparison](MULTIVARIATE_EFFECTIVE_VANISHING_CURVATURE.md).
We use equality of root arguments to compare the exact solutions of
the modulus equations in its regular charts. In the fold chart,
the modulus localization already controls the possible difference
of real fields tightly enough to identify the Airy label.

## Field normalization and the common argument

Take \(W\) small enough that all unequal rows are shared, and write
\(\varepsilon=\sqrt{W/2}\). For order \(j\) the normalization in
the effective comparison has
\[
        t=-\exp\!\left(-\sum_i d_i\log(jm_i)+\varepsilon\xi\right),
           \qquad m_i=(a_i+b_i)/2,\quad d_i=b_i-a_i .
\]
The displayed sum is real. Thus equal physical arguments at orders
\(n,N\) allow representatives with \(\Im\xi_n=\Im\xi_N=y\) and
\(|\varepsilon y|\le\pi\) for any positive rescaling of either
polynomial. Conjugating both roots reduces to
\(y\ge0\). Negative real roots use \(y=0\), and positive real
arguments cannot be roots because the polynomial coefficients are
positive. The real part bound \(|\Re \xi_j|\le8\varepsilon\) from
the effective comparison holds at every root separately and is
independent of any relation between the two real parts.

The large-field and one-saddle regions there exclude single-order
roots. Its regular two-saddle region and fixed small fold region
exhaust the remaining bounded \(y\)-interval. Both fields enter
the same region since \(y\) is common and \(|\Re\xi_j|\) is
uniformly small. Shrink the physical curvature cap as needed
after fixing the two regions.

## Regular modulus comparison

Use the two retained saddle packets below the fold, with
\(D=P_R-P_L\) their stationary phase difference.
Here \(\Re D'\ge c_D>0\) on each horizontal joining segment,
and \(|e^D-1|\ge g_D>0\) on the regular field region.
These bounds are those of the
[regular comparison](MULTIVARIATE_EFFECTIVE_VANISHING_CURVATURE.md#uniform-regular-comparison-off-a-fixed-fold-neighborhood).
Near a retained saddle the phase may use any continuous local
logarithm of the nonzero base, because \(e^{j\log g}=g^j\)
for every integer \(j\). In particular at \(y=0\) the argument
\(\pi\) of the base at the left saddle is part of \(P_L\).
Its order parity enters the exponential. The leading
amplitude convention stays the same at both consecutive orders.
The possible vanishing base between the saddles only enters the
strict-loss estimate for the modulus off the local packets.

At an exact root in the chart, the Gaussian comparison gives
\[
        \Gamma(\xi_j)e^{jD(\xi_j)}=-(1+e_j),\qquad
                    |e_j|\le K/\sqrt n ,
 \tag{1}
\]
where the analytic nonvanishing leading amplitude ratio
\(\Gamma\) is independent of \(j\) for a fixed profile.
Indeed the effective density representation has
\(p_j(s)=\sqrt j\,e^{jA(s)}(a_0(s)+O(1/n))\), where the shared
profile determines the same \(A,a_0\) and moment base
\(g_\varepsilon\) at each order.
Saddle locations, Hessians and leading local contour factors
are consequently independent of \(j\), and the two local Gaussian
powers cancel in the ratio. The relative errors in (1) follow from
the fixed quadratic losses and the strict tail and connector losses
for each packet. An \(O(n^{-1/2})\) error suffices here.

On the compact regular chart fix effective upper bounds \(M_\Gamma\)
for \(|\log|\Gamma||\), \(L_\Gamma\) for \(|\Gamma'/\Gamma|\),
and \(L_D\) for \(|D'|\). Choose \(n\) large enough that
\(K/\sqrt n\le1/4\). Take absolute values and logarithms in (1):
\[
 \left|j\Re D(\xi_j)+\log|\Gamma(\xi_j)|\right|
                                      \le 2 K/\sqrt n.
 \tag{2}
\]
In particular \(|\Re D(\xi_j)|\le(M_\Gamma+1)/n\).
Let \(\delta=\xi_N-\xi_n\), which is real.
The real derivative lower bound for \(D\) first gives
\[
       |\delta|\le A/n,\qquad A=2(M_\Gamma+1)/c_D.
\]
Subtract (2) in its signed form at the two orders, using this
first bound and the derivative bound for \(\log|\Gamma|\).
One obtains
\[
    N|\delta|\le B/\sqrt n,\qquad
    B=\bigl(4K+M_\Gamma+1+L_\Gamma A\bigr)/c_D .
 \tag{3}
\]
In detail, the subtracted leading real phase is
\(N(\Re D(\xi_N)-\Re D(\xi_n))+\Re D(\xi_n)\).
The extra real phase is \(O(1/n)\), and the difference of log
amplitudes is at most \(L_\Gamma|\delta|\). Thus (3) bounds
the real separation without requiring an integer real drift.

Now divide the two complex equations (1). The ratio of amplitudes
can be written using the integral of \(\Gamma'/\Gamma\) along
the horizontal segment. Its logarithm together with
\(N(D(\xi_N)-D(\xi_n))\) has modulus at most
\(B(L_\Gamma+L_D)/\sqrt n\). For \(n\) large enough that this
is at most \(1/2\), the division gives the explicit bound
\[
        |e^{D(\xi_n)}-1|
             \le \bigl(6K+2B(L_\Gamma+L_D)\bigr)/\sqrt n.
\]
Taking the right side smaller than \(g_D\) excludes equal
arguments throughout the regular region.

## Airy labels at the same imaginary field

Use the cubic chart and analytic finite integrals of the
[effective fold comparison](MULTIVARIATE_EFFECTIVE_VANISHING_CURVATURE.md#effective-airy-residual-and-linked-roots).
Write \(\Lambda(\xi)\) for its physical Airy parameter and
\(z_j=j^{2/3}\Lambda(\xi_j)\) for its scaled coordinate.
There are effective constants \(C_0>0,c_0>0,L_0>0\) on a fixed
convex fold disk such that at any true zero
\[
 z_j=a_{k_j}+O(C_0n^{-1/3}),\qquad
           -\Im\Lambda'(\xi)\ge c_0,\qquad |\Lambda'(\xi)|\le L_0 .
\]
Here \(a_k\) is a negative real Airy zero. Equal imaginary
physical fields and \(|\Im\Lambda(\xi_j)|\le C_0/n\) give
\[
             |\delta|\le2C_0/(c_0n).
\]
Consequently, if \(|\Lambda|\le1\) in the admitted disk,
\[
 |z_N-z_n|
 \le N^{2/3}L_0|\delta|
       +(N^{2/3}-n^{2/3})|\Lambda(\xi_n)|
 \le \bigl(4C_0 L_0/c_0+1\bigr)n^{-1/3}.
\]
Thus the two labels are at distance at most \(C_1 n^{-1/3}\)
with \(C_1=2C_0+4C_0L_0/c_0+1\).
The Airy equation at a simple zero \(a\) excludes another zero
within \(c_1(1+|a|)^{-1/2}\) for a positive effective \(c_1\),
by the normalized integral equation in that comparison.
Choose the admitted bound \(\mu\) on \(|\Lambda|\) after the
localization constants, and then take \(n\) large enough that
\[
       \mu+n^{-2/3}+C_0/n < \bigl(c_1/(2C_1)\bigr)^2 .
\]
Since \(1+|a_{k_n}|\le1+n^{2/3}\mu+C_0 n^{-1/3}\), this makes
another label impossible. Both fields therefore have the
same Airy label.

The common-label branch in that comparison has
\[
 \xi'(j)=(\Lambda^{-1})'
             \bigl(-\tfrac23 a_k j^{-5/3}+O(n^{-2})\bigr),
                                      \qquad n\le j\le N .
\]
Its inverse derivative has imaginary part uniformly positive
on the fold chart. The Airy zero \(a_k<0\) has a uniform
positive lower bound on \(|a_k|\), so for sufficiently large
\(n\) the imaginary motion over one order is strictly
positive and bounded below by a positive multiple of
\(|a_k|N^{-5/3}\). The actual integer roots differ from the
interpolating branch by exponentially small errors, which
cannot remove that imaginary displacement. The negative
imaginary fold follows by conjugation. The two real parts may differ in either direction within the common chart.

## Finite hierarchy and the all-order cap

Fix the large-field exclusion, a fold coordinate chart and its
finite contours and Airy constants as in the effective comparison.
The constants \(C_0,c_0,L_0,c_1\) are taken on those fixed domains.
Shrink the admitted fold field radius, if necessary, for the label
inequality and common-label branch above. Next choose a covered
interval of imaginary fields and take the effective constants
\(c_D,g_D,M_\Gamma,L_\Gamma,L_D,K\) on the remaining regular
regions. This uses compact domains bounded away from the fold.
The amplitude estimates are uniform on those domains because the
density and phase bounds depend on the rows only through \(W\).

Decrease \(\varepsilon\) to ensure the real field bound enters
these charts while preserving all strict contours. This gives a
computable positive rational curvature cap \(W_{\rm a}\).
Choose an integer \(N_{\rm a}\ge3\) large enough to satisfy the
regular inequalities, the Airy label inequality, and the derivative
and residual bounds at the fold. These are strict inequalities
with the constants already fixed. The resulting angular separation
holds for every \(N\ge N_{\rm a}\) at \(0<W\le W_{\rm a}\).

An absolute all-order cap can now be taken as
\[
                 \eta_*=\min(W_{\rm a},10^{-18}N_{\rm a}^{-3}).
\]
For \(3\le N<N_{\rm a}\), use \(S\le W\) and the
[explicit matched Gauss remainder](MULTIVARIATE_MATCHED_GAUSS_UNIFORM_REMAINDER.md)
at the smaller curvature cap. For \(N=2\) the first moment has one
negative real root. The second moment is strictly positive on
real arguments since it is the expectation of a square and the
two monomials are distinct.
Finally, for \(n\ge2\) and fully shared profiles, a row of smaller
exponent \(a\) and gap \(d\) gives
\[
   n\sum_{h=1}^d\frac{d}{na+h}
              \ge\frac{d^2}{a+d/n}\ge\frac{d^2}{a+d/2}.
\]
Summing yields \(nX\ge W\), which proves the additional cap
formulation.
