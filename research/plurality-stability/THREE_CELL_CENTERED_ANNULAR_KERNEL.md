# A centered radial kernel for the annular comparison

Balance permits an exact subtraction from the raw radial mass
kernel while retaining the centered angular terms. In the annular
proof of (JA1), this improves the inner mixed bound from
\(O(\sqrt{r_0})\) to \(O(r_0^{5/2})\) and removes the inverse inner
cutoff from the exceptional-set mean estimate. The remaining
quantitative local amplitude threshold still depends on the angular
comparison. These estimates apply to ordinary balanced planar profiles
at every compact interval of positive correlations.
[The annular comparison (JA1)](THREE_CELL_JOINT_ANGULAR_ENVELOPE.md)

## Exact centering and a weighted kernel bound

Use \(d=1-\rho^2\), \(d\mu(r)=r e^{-r^2/2}dr\), and
\(b(r)=p(r)-\mathbf1/3\). Write
\[
 \alpha=\frac{\rho^2}{2d},\quad \lambda=\frac{\rho}{d},\qquad
 P(r,s)=d^{-1}e^{-\alpha(r^2+s^2)}I_0(\lambda rs).
\]
In (DG2) the raw mass term is
\(\iint P(r,s)b(r)\cdot b(s)d\mu(r)d\mu(s)\).
Balance \(\int b\,d\mu=0\) permits replacement of \(P\) there by
\[
 \begin{split}
 P^c(r,s)&=P(r,s)-P(r,0)-P(0,s)+P(0,0)\\
 &=d^{-1}\big[(1-e^{-\alpha r^2})(1-e^{-\alpha s^2})
       +e^{-\alpha(r^2+s^2)}(I_0(\lambda rs)-1)\big].
 \end{split}                                                \tag{CA1}
\]
This is a replacement in the raw mass term. The joint measure
for the angular-centered term remains \(P\,d\mu\,d\mu\).
[Polar stability identities (DG1)–(DG2)](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md#exact-separation-of-radial-masses-and-angular-changes)

The Bessel series gives \(I_0(z)-1\le z^2 I_0(z)/4\). Thus
\[
 0\le P^c(r,s)\le r^2s^2\left(\frac{\alpha^2}{d}
                                    +\frac{\lambda^2}4 P(r,s)\right).
                                                               \tag{CA2}
\]
Let \(J=[\rho_-,\rho_+]\Subset(0,1)\),
\(d_*=1-\rho_+^2\), \(\alpha_*=\rho_+^2/(2d_*)\),
\(\lambda_*=\rho_+/d_*\), and
\(c_*=(1-\rho_+)/[2(1+\rho_+)]\).
Put \(w(r)=\min(r,1)\). For \(0<r_0\le1\), (CA2) implies
\[
 \left\|1_{\{\min(r,s)<r_0\}}\frac{P^c(r,s)}{\sqrt{w(r)w(s)}}\right\|_{\rm HS}
 \le \sqrt{\frac25}\,r_0^{5/2}
 \left[
   \frac{\alpha_*^2}{d_*}\Psi(1/2)
       +\frac{\lambda_*^2}{4d_*}\Psi(c_*)\right],
 \qquad
 \Psi(c)=\left(c^{-3}+\frac{3\sqrt\pi}{8c^{5/2}}\right)^{1/2}.
                                                               \tag{CA3}
\]
The norm uses \(d\mu(r)d\mu(s)\). In fact \(I_0(\lambda rs)\le
e^{\lambda rs}\) and \(2rs\le r^2+s^2\) give
\[
 P(r,s)^2d\mu(r)d\mu(s)
      \le d_*^{-2}rs\,e^{-c_*(r^2+s^2)}\,dr\,ds .
\]
For \(r<r_0\) we have \(r^4/w(r)=r^3\). The squared kernel
bound from either term of (CA2) then integrates in \(r\) with
at most the factor \(r_0^5/5\). The integral in \(s\) is bounded
using \(s^5/w(s)\le s^5+s^4\), whose Gaussian integrals give
\(\Psi(c)^2\). Add the two Hilbert–Schmidt bounds and account
for both orders of the inner endpoint to obtain (CA3).

For middle radii \(M=[r_0,R]\), the simpler pointwise bound is
\[
 \frac{P^c(r,s)}{\sqrt{rs}}\le D_{J,R}(rs)^{3/2},\qquad
 D_{J,R}=\frac{\alpha_*^2}{d_*}
      +\frac{\lambda_*^2}{4d_*}
                      \exp\!\left(\frac{\rho_+R^2}{1+\rho_+}\right).
                                                               \tag{CA4}
\]
Indeed \(P\le d_*^{-1}\exp(\rho_+R^2/(1+\rho_+))\) on
\([0,R]^2\) by its Gaussian exponent and \(2rs\le r^2+s^2\).
For pairs in \(M^2\) touching a bad set \(B\subset M\), the
Hilbert–Schmidt norm in (CA4) is bounded by
\[
 \sqrt2 D_{J,R}
   \left(\int_B r^3d\mu\int_M s^3d\mu\right)^{1/2}.
                                                               \tag{CA5}
\]
This gives an \(O_{J,R}(\sqrt{\mu(B)})\) middle mixed cost
with no singular constant as \(r_0\downarrow0\). It is
stronger when the bad set lies near the origin.

## The middle quadratic and the mass error

Set \(G=M\setminus B\) and \(X_G=\int_G r|b|^2d\mu\).
After switching the raw mass term to \(P^c\), the quadratic
on \(G^2\) still uses the expansion (JA4). Its difference
from the ordinary quadratic with kernel \(P\) is exactly
\[
 -2M_G\cdot L_G,\qquad
 M_G=\int_G b\,d\mu,\qquad
 L_G=\int_G\big(P(r,0)-P(0,0)/2\big)b(r)\,d\mu(r).
                                                               \tag{CA6}
\]
In fact the two removed \(P(r,0)\) terms each contribute
\(-M_G\) times a single integral, and the added \(P(0,0)\)
term is \(P(0,0)|M_G|^2\). Since
\(|P(r,0)-P(0,0)/2|\le1/(2d_*)\), Cauchy–Schwarz
gives \(|L_G|\le(\pi/2)^{1/4}\sqrt{X_G}/(2d_*)\).
For any \(\eta>0\) the correction in (CA6) costs at most
\(\eta X_G+(\pi/2)^{1/2}|M_G|^2/(4\eta d_*^2)\).
This can be included with the usual mean correction in (JA6).

Write \(E=[0,\infty)\setminus G\), \(Y_E=\int_E w|b|^2d\mu\),
and \(\delta=\mu(B)\). For the annular bad set, rearrangement
of the decreasing weight \(1/r\) gives the exact bound
\[
 \int_B \frac{d\mu(r)}r
 \le \int_0^{r_\delta}e^{-r^2/2}dr\le r_\delta,\qquad
 r_\delta=\sqrt{-2\log(1-\delta)}.                            \tag{CA7}
\]
If \(\delta\le 1-e^{-1/2}\), then \(r_\delta\le1\). The same
bound holds for \(1/w\) by rearranging \(1/w\) onto
\([0,r_\delta]\). Therefore balance and Cauchy–Schwarz give
\[
 |M_G|^2\le
       [r_0+r_\delta+\mu((R,\infty))]Y_E .                    \tag{CA8}
\]
In particular a bad-set contribution can be controlled by
\(r_\delta\) rather than \(\mu(B)/r_0\), without using the
position of the bad set in the annulus.

## The terms that still set the cutoff

Apply the centered angular bound of Heilman's Lemma 6.1 on pairs
touching \(E\). For the normalized kernel of (DG1), the lemma
applies to arbitrary circular partitions of either width
profile. Its penalty is at least
\(c_0\min(\lambda rs,1)(|b(r)|^2+|b(s)|^2)\), where one may take
\(c_0=(.109/2)(1-e^{-\pi/6})\). For the partner radius conditioned
on \(r\), the event that the independent Gaussian component
perpendicular to the mean has standard coordinate of absolute
value at least one has probability \(2\overline\Phi(1)\)
and implies \(s\ge\sqrt d\). On this event
\(\min(\lambda rs,1)\ge\rho_- w(r)\).
Both orders of the pair then give a negative term at least
\(2c_0\,2\overline\Phi(1)\rho_-Y_E\). This choice of weight makes
the bad-set diagonal uniform in \(R\). The mixed term on those
pairs is the centered raw kernel, controlled by (CA3), (CA5),
and the Gaussian tail bound following (CA3). The latter remains
integrable after multiplication by the polynomial factors in
(CA2), so it gives a uniformly vanishing outer mixed cost as
\(R\to\infty\).
[Heilman, Lemma 6.1 and normalization (7)](https://arxiv.org/html/2306.03312v3)

The inner raw mixed cost is now \(O_J(r_0^{5/2})\) in the
normalized energy, but \(O_J(r_0)\) remains in the good-set
mean correction (CA8). The missing good-set quadratic diagonal
is still estimated as in (JA6), with costs from small partner
radii, tails and \(\mu(B)\). On the good pairs, JA4 still
requires an angular parameter cutoff and a mass threshold that
ensures the joint envelope expansion. Its proof of (JA2) gives
that mass threshold by compactness, and its cubic remainder
constant is unspecified. Thus these estimates quantify several
radial errors without giving an explicit admissible amplitude
for (JA1). A quantitative angular envelope estimate and its
remainder are still needed for that application.

The bounds are analytic with internal checking. They do not
change the proven correlation interval.

## An exponential reserve for two central profiles

There is an independent way to strengthen the angular reserve
when two profiles are central, without using the shape of
maximizing circular pairs. Use \(t=r^2/2\) in this section,
\(d\nu(t)=e^{-t}dt\), and the parameters
\(\kappa_\rho=\rho^2/d\), \(\gamma=2\rho/d\),
\(a=\gamma\sqrt{ts}\). Write \(z(t)=p(t)-\mathbf1/3\).
For \(1/3<h\le5/12\) let \(G_h\) be the central set of
(GC1), so every coordinate of its profiles lies between
\(2/3-h\) and \(h\). Write
\[
 d_h=2(h-1/3),\qquad m_h=-\cos\big(\pi(2/3-d_h)\big),\qquad
 u=\frac{3}{2\pi^2}.
\]
For \(t,s\in G_h\) each displacement of \(p(t)+p(s)\) from
\(2\mathbf1/3\) has magnitude at most \(d_h\le1/6\).
The positive Taylor ratio of its coordinate is an average
of \(e^{a\cos(\pi v)}\) over values
\(v\in[2/3-d_h,2/3+d_h]\). It is therefore at most
\(e^{-m_h a}\). In the normalization of (SR3), the
positive Taylor slack on these pairs is at least
\[
       \frac{P_\sharp(a)-e^{-m_h a}}4\,|z(t)+z(s)|^2 .
                                                               \tag{CA9}
\]
[Sharp positive coefficient (SP1)–(SP4)](THREE_CELL_SHARP_POSITIVE_ARC_ENVELOPE.md#sharp-angular-inequality)
[Central profiles (GC1)–(GC3)](THREE_CELL_GRADED_CENTRAL_EXTENSION.md#one-central-partner-and-two-central-partners)

Under the probability density of \(P_\sharp(a)\) in (SP4), let
\(X=\cos(\pi v)\). Then \(\mathbb E X=-u\). The tangent
identity gives
\[
 P_\sharp(a)-e^{-ua}
 =\mathbb E\big[e^{aX}-(1+(X+u)a)e^{-ua}\big].
\]
The integrand is nonnegative by convexity. On \(X\ge0\), the
subtracted term is at most \(1+Xa\), since
\((1+Xa)e^{ua}\ge(1+Xa)(1+ua)\ge1+(X+u)a\).
Hence
\[
 P_\sharp(a)\ge e^{-ua}
             +\sum_{n\ge2}\frac{a^n}{n!}M_n^+,\qquad
 M_n^+=\int_0^{1/2}3v\cos^n(\pi v)\,dv .
                                                               \tag{CA10}
\]
All the series terms are nonnegative. For \(m_h>u\), the
exponential difference satisfies
\[
 e^{-ua}-e^{-m_h a}
       \ge(m_h-u)a\exp\!\left(-\frac{m_h+u}{2}a\right)
\]
by Jensen on the integral of \(a e^{-va}\) for
\(v\in[u,m_h]\). Thus (CA9) gives a linear reserve with
damping parameter \((m_h+u)/2\), as well as the undamped
terms \(M_n^+a^n/(4n!)\) for every \(n\ge2\). For example
\(h=7/20\) gives damping less than \(.28\). The undamped
degree-two coefficient is
\[
      \frac{M_2^+}{8}=\frac{3}{128}-\frac{3}{32\pi^2}>0 .
\]
These are parts of the positive Taylor slack itself. The
same-dominant-label reserve vanishes when both profiles lie
in \(G_h\). If the one-central GC3 coefficient is used
throughout the union of central events, its both-central
portion must be deducted before adding (CA9). Alternatively
use the one-central estimate only on the pairs with exactly
one central profile. The first angular moment deficit is
separate from these scalar Taylor terms.

For any of the nonnegative terms proportional to \(a^ne^{-c a}\)
with \(c\ge0\), the planar radial factor yields the rank-one
minorant
\[
 \frac1d e^{-\kappa_\rho(t+s)} a^n e^{-c a}
 \ge\frac{\gamma^n}{d}
     t^{n/2}e^{-(\kappa_\rho+c\gamma/2)t}
     s^{n/2}e^{-(\kappa_\rho+c\gamma/2)s}.
\]
The single-radius normalization is
\(Z_{n,c}=\Gamma(1+n/2)[d/(1+c\rho)]^{1+n/2}\).
Let \(\mu_{n,c}\) be its normalized measure relative to
\(\nu\). Restricted to \(G_h^2\), the integral of
\(|z(t)+z(s)|^2\) against \(\mu_{n,c}\otimes\mu_{n,c}\)
is at least
\(2\mu_{n,c}(G_h)\int_{G_h}|z|^2d\mu_{n,c}\).
In addition, weighted orthogonality in a trial decomposition
\(z=\alpha\phi+g\), \(\|\phi\|_J=1\), improves the selected
raw norm estimate. Let \(T\) be the multiplication form for
integration against \(\mu_{n,c}\) on any selected subset of
\(G_h\), with \(0\preceq T\preceq B I\) in the weighted
trial space. If \(U=\langle\phi,T\phi\rangle\) and
\(r=\|g\|_J/|\alpha|\), then
\[
 \|T^{1/2}z\|_J^2\ge |\alpha|^2
                     (\sqrt U-\sqrt{B-U}\,r)_+^2 .
                                                               \tag{CA11}
\]
For \(U>0\) project onto
\((\alpha/|\alpha|)T^{1/2}\phi/\sqrt U\), using
\(\|P_{\phi^\perp}T\phi\|_J^2\le U(B-U)\);
the case \(U=0\) is immediate.
Thus (CA9)–(CA11) supply a finite-amplitude reserve in terms
of tilted central mass and norm. They supply no finite
correlation endpoint by themselves.
