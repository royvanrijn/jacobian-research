# A signed force expansion for two matched Gegenbauer slots

The angular force at a root of a two-slot circular convolution
admits a uniform signed expansion from the exact differential
equation. It compares a root with the positive characteristic
exponent obtained by freezing the equation at that root.
For the ordinary doubled shared unit pair
\((1,2),(2,1),(1,2),(2,1)\), the first correction vanishes at
the limiting first angular edge. The resulting finite-order
reserve proves strict consecutive interlacing for all sufficiently
large degrees by an argument within the exact ordinary array.

## Force comparison from the right boundary condition

Use the notation of the
[two-slot Sonin bound](MULTIVARIATE_TWO_SLOT_SONIN_POLAR_BARRIER.md#two-upper-bounds-on-the-force).
In degree \(n\), put \(K=n/2\), \(A_j=K+\lambda_j\), and
\[
 B=K+\lambda_1+\lambda_2,\quad
 C=K^2-\lambda_1\lambda_2,\quad
 \Omega^2=K A_1 A_2/B,\quad D=\Omega^2-C.
\]
Suppose \((\lambda_1/n,\lambda_2/n)\) ranges in a compact subset
of positive pairs with product greater than \(1/4\).
Write \(f\) for the lower centered polynomial. Let \(\theta\)
be one of its roots in \(0<\theta\le\theta_0<\pi\) and normalize
\(f'(\theta)=1\). In displacement \(s\) from that root put
\(g=f'\), \(w=f''+\Omega^2 f\), and
\[
 a(s)=B\tan((\theta+s)/2),\qquad a=a(0),\qquad
 \rho=-f''(\theta)>0 .
\]
The exact equation is \(w'=a(s)w+Dg\). In particular \(w(\pi)=0\)
in the angular variable. The Sonin comparison gives
\(|g(\theta+s)|\le1\) and hence \(|f(\theta+s)|\le s\).
Its functional is nonnegative to the right because it decreases
and its terminal value is \(g(\pi)^2+\Omega^2f(\pi)^2\).
Therefore \(|w(\theta+s)|\le\sqrt D\sqrt{1+\Omega^2s^2}\).
The integrated equation also gives
\(a(s)|w(\theta+s)|\le D\) for \(0<s<\pi-\theta\).

Let \(\lambda\) be the positive root of
\[
         \lambda^3-a\lambda^2+C\lambda-a\Omega^2=0.
\]
It is unique by Descartes' sign rule when \(a>0\), with the
continuous value \(\sqrt{-C}\) at \(a=0\). Set
\[
 b=\lambda-a,\qquad d=a\Omega^2/\lambda,\qquad
 H=\lambda^2+b\lambda+d .
\]
Then \(b>0\), \(C=d-b\lambda\), and the function
\[
 {\cal L}=w+b g-b\Omega^2 f/\lambda=f''+b f'+d f
\]
obeys the exact first-order equation
\[
       {\cal L}'-\lambda{\cal L}=(a(s)-a)w.                       \tag{1}
\]
The claimed expansion, uniform over the indicated shapes and
angles, is
\[
       \rho=b+a'(0)\frac{b(3\Omega^2-\lambda^2)}{H^2}
                  +O(n^{-1}).                                   \tag{2}
\]
Here the \(O(n^{-1})\) constant depends on the shape compact
set and on \(\theta_0\).

To verify the error, write \(S=\pi-\theta\), \(h(s)=e^{-\lambda s}\).
All functions inside the following integrals are evaluated at
\(\theta+s\). Integration of (1) backwards gives
\[
 \int_0^S h(s)s^j{\cal L}(\theta+s)\,ds
   =\frac{S^{j+1}}{j+1}h(S){\cal L}(\pi)
      -\frac1{j+1}\int_0^S h(s)s^{j+1}(a(s)-a)w(\theta+s)\,ds .
                                                               \tag{3}
\]
The terminal term is exponentially small since \(S\) is bounded
below and \(\lambda/n\) is bounded below. On
\(0\le s\le(\pi-\theta_0)/2\), we have \(a(s)-a=O(ns)\)
and \(|w|=O(n(1+ns))\). The complementary integral has
\(|(a(s)-a)w|\le D=O(n^2)\) and an exponentially small
weight. Thus the integrals of \(h{\cal L}\) and \(hs{\cal L}\)
are \(O(n^{-1})\) and \(O(n^{-2})\), respectively.

Let \(J_j=\int_0^S h(s)s^j f(\theta+s)\,ds\) and
\(L_j=\int_0^S h(s)s^j{\cal L}(\theta+s)\,ds\).
Integrating by parts in \({\cal L}=f''+bf'+df\), and using
\(f(\theta)=0\), \(f'(\theta)=1\), gives up to exponentially
small terminal terms
\[
 H J_0=1+L_0,\qquad H J_1=(2\lambda+b)J_0+L_1 .
\]
A second integration by parts gives
\[
 \begin{aligned}
 \int_0^S h(s)s w(\theta+s)\,ds
  &=(\lambda^2+\Omega^2)J_1-2\lambda J_0
       +O(e^{-c n})\\
  &=\frac{b(3\Omega^2-\lambda^2)}{H^2}(1+L_0)
       +\frac{\lambda^2+\Omega^2}{H}L_1+O(e^{-c n})\\
  &=\frac{b(3\Omega^2-\lambda^2)}{H^2}+O(n^{-2})
 \end{aligned}                                                    \tag{4}
\]
for a constant \(c>0\). Here \(H\) is bounded above and below
by positive multiples of \(n^2\), and
\(\Omega^2-d=b\Omega^2/\lambda\) yields the numerator in
the second line.

At the root, \({\cal L}(\theta)=b-\rho\), so (1) gives
\(\rho-b=\int_0^S h(s)(a(s)-a)w(\theta+s)\,ds-h(S){\cal L}(\pi)\).
Replacing \(a(s)-a\) by \(a'(0)s\) costs \(O(n^{-1})\).
On the initial half interval the Taylor remainder is \(O(ns^2)\),
and the bound on \(w\) gives this error after integration.
The complementary integral is exponentially small, using
\(|(a(s)-a)w|\le D\) and the Sonin bound for the subtracted term.
Now (4) proves (2).

## The ordinary doubled shared unit pair

Here \(\lambda_1=\lambda_2=n+1\). The
[exact ordinary upper sign](MULTIVARIATE_TWO_SLOT_SONIN_POLAR_BARRIER.md#an-angular-range-for-two-equal-shared-unit-pairs)
reduces to \(\rho<P(t)/Q(t)\), \(t=\tan(\theta/2)\), with
\(P,Q\) as displayed there. That result covers \(t\ge1/5\)
in every degree. Write \(\delta=1/n\), \(r=P/(nQ)\), and
use bars for normalization by \(n\) for \(a,b,\lambda,B\), and
by \(n^2\) for \(C,d,\Omega^2\). Near the small-angle edge
these are smooth functions of \(t,\delta\), determined by
\[
 \bar a=(5/2+2\delta)t,\quad
 \bar C=1/4-(1+\delta)^2,\quad
 \bar\Omega^2=\frac{(3/2+\delta)^2}{5+4\delta},\quad
 F(\bar b)=(\bar b+\bar a)^3
    -\bar a(\bar b+\bar a)^2+\bar C(\bar b+\bar a)
       -\bar a\bar\Omega^2=0.
\]
Define \(s_f=(4\bar d-\bar b^2-\delta^2)/3\).
By the exact root force identity, at a lower root
\[
 0\le\frac1{n^2}\sum_{j\ne i}\csc^2((\theta_i-\theta_j)/2)
    =(\rho/n)^2+\frac43\bar a\,(\rho/n)
                                +\frac{4\bar C+\delta^2}{3}.     \tag{5}
\]
Replacing \(\rho/n\) by \(\bar b\) in the last expression gives
\(s_f\).

At \(\delta=0\), \(s_f(t,0)\) is negative on \(0\le t<t_*\)
and positive on \(t>t_*\), where \(t_*=\sqrt{15}/25\).
Indeed a zero requires \(\bar b^2=4\bar d\). With
\(u=\bar a/\bar b\), the identities for \(C,d\) give
\(\bar\Omega^2/(-\bar C)=(1+u)/(u(3+4u))=3/5\).
The unique positive solution is \(u=1/2\), yielding the stated
\(t_*\), and the signs follow from the sign at zero and the
positive derivative below.

At \(\delta=0\) and \(t>t_*\), \(r>\bar b\). In fact exact
substitution into the cubic gives
\[
 F(r)=
 \frac{972t(1+t^2)^2(125t^2-3)(125t^4+302t^2+225)}
            {(125t^4+234t^2+45)^3}>0 ,
\]
and the positive root in \(\lambda=\bar b+\bar a\) is unique.
For a sequence of lower roots with \(0\le t\le1/5\) and
\(n\to\infty\), (2) and (5) therefore exclude a limit below
\(t_*\). A limit strictly above \(t_*\) gives the desired
strict sign directly.

At \((t_*,0)\), \(\bar\lambda^2=3\bar\Omega^2\). The leading
correction in (2) then implies
\[
                  \rho-b=O(|t-t_*|+\delta).
\]
The exact derivative values from the displayed cubic and the
rational \(P/Q\) are
\[
 \begin{array}{c|cc}
     & \partial_t & \partial_\delta\\ \hline
       s_f & \sqrt{15}/3 & -32/75\\
       r-\bar b & 50/11 & -16\sqrt{15}/165
 \end{array}
\]
with \(s_f=r-\bar b=0\) there. In particular the local curve
\(s_f(t_e(\delta),\delta)=0\) satisfies
\(t_e=t_*+O(\delta)\). At a root with \(t\to t_*\),
(5) and the improved error give
\(s_f\ge-O(\delta(|t-t_*|+\delta))\). The positive derivative
of \(s_f\) in \(t\) and \(|t-t_*|\le|t-t_e|+O(\delta)\)
then imply \(t\ge t_e-O(\delta^2)\).

At the curve \(t_e\) the derivative table gives
\[
 n(r-\bar b)(t_e(\delta),\delta)
   =-\frac{16\sqrt{15}}{165}
        +\frac{50}{11}\frac{32/75}{\sqrt{15}/3}
        +O(\delta)
   =\frac{16\sqrt{15}}{55}+O(\delta).
\]
The \(t\)-derivative of \(r-\bar b\) remains positive nearby,
so at a lower root approaching this edge
\[
                  \frac P Q-\rho\ \ge\ \frac{16\sqrt{15}}{55}-o(1).
\]
The sign holds for all sufficiently large \(n\) by compactness
of \(0\le t\le1/5\). The existing larger-angle bound, reflection,
and the distinct central-root parities finish eventual strict
interlacing for this ordinary profile. Formula (2) gives an
analytic route at the edge. An explicit cutoff and the remaining
small angles at bounded degrees would be needed for an
every-degree conclusion by this argument.
