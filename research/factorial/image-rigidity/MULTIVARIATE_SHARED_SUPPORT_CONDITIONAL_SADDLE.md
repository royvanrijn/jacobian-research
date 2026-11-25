# Conditional saddle integrals for monomials with shared support

Let \(a,b\in\mathbb Z_{>0}^r\) be distinct. After removing
coordinates with equal exponents, set \(d_i=b_i-a_i\ne0\),
\(A=\sum_i a_i\), and \(\Delta=\sum_i d_i\). Work with the
ordinary factorial functional and
\[
 M_n(t)=\Lambda((X^a+tX^b)^n),\qquad
 I_n(\zeta)=n^{-nA-r}M_n(n^{-\Delta}\zeta).
\]
This reduction applies with mixed signs of \(d_i\) and with
arbitrarily many varying coordinates.

The integral for \(I_n\) reduces to one complex variable on
the strip \(|\Im s|<\pi/2\). Its stationary equation has
exactly two critical points left after counting the exterior
real critical points of the rational parameter map. Every
nonzero real scale other than an interior cusp has a neighborhood
where consecutive moments have no simultaneous zero for all
sufficiently large orders. This includes real coexistence with two
dominant saddles and equal total degrees. At an interior cusp
the argument instead covers every fixed bounded quartic-scaled
field window, again including equal total degrees.

## An exact conditional integral

Put \(x_i(q)=a_i+d_iq\) and let
\[
 {\cal V}=\{q\in\mathbb C:\Re x_i(q)>0\ \text{for all }i\},\qquad
 s(q)=\sum_i d_i\log x_i(q).
\]
Use the right-half-plane logarithm. The vertical strip
\({\cal V}\) contains \((0,1)\) and can be unbounded on one
side. On this convex domain
\[
 \Re s'(q)=\Re\sum_i\frac{d_i^2}{x_i(q)}>0.
\]
Integrating the derivative on the segment between any two
points proves injectivity. On the real interval
\({\cal V}\cap\mathbb R\), the function \(s(q)\) increases
from \(-\infty\) to \(+\infty\). At finite endpoints use
the logarithm of a vanishing \(x_i\). At an infinite endpoint
all the \(d_i\) have the same sign and
\(\Re s(q)=\Delta\log|q|+O(1)\).

The image contains the whole horizontal strip
\(\Sigma=\{s:|\Im s|<\pi/2\}\). To prove this, every term
\(d_i\Im\log x_i(q)\) has the sign of \(\Im q\).
At a finite vertical side approached with a nonzero imaginary
part, a coordinate touching that side contributes a limiting
absolute value \(\pi|d_i|/2\ge\pi/2\).
At a real endpoint the real part of \(s(q)\) escapes, since
the vanishing logarithms have coefficients of the same sign.
At \(|q|\to\infty\) one has
\(\Re s(q)=\Delta\log|q|+O(1)\). If \(\Delta=0\),
the signs are mixed, so \(\Re q\) is bounded in
\({\cal V}\), and
\[
 \lim_{|q|\to\infty}|\Im s(q)|
    =\frac{\pi}{2}\sum_i |d_i|\ge\pi .
\]
These exclusions show that the image is relatively closed
in \(\Sigma\): a sequence with images converging in
\(\Sigma\) has a subsequence converging inside \({\cal V}\).
The image is open and contains the real line, proving the
claim. Write \(q=q(s)\) for the holomorphic inverse on
\(\Sigma\), and use \(x_i=x_i(q(s))\) below.

Scale the original exponential integral by \(X_i=ny_i\).
For real \(s\) use logarithmic coordinates \(y_i=x_i e^{u_i}\)
with \(u\) on the Euclidean hyperplane
\(W=\{u:\sum_i d_i u_i=0\}\).
The equation \(d\cdot\log y=s\) gives a Jacobian factor
\(1/\|d\|\) in logarithmic coordinates: the \(s\)
column has scalar product one with \(d\).
On \(W\), \(\sum_i(a_i-x_i)u_i=0\).
Consequently the exact formula is
\[
 \begin{aligned}
 I_n(\zeta)&=\int_{\mathbb R}(1+\zeta e^s)^n F_n(s)\,ds,\\
 F_n(s)&=e^{n a(s)}\frac{\prod_i x_i}{\|d\|}
 \int_W \exp\!\left(\sum_i u_i-n\sum_i x_i(e^{u_i}-1-u_i)\right)du,\\
 a(s)&=\sum_i a_i\log x_i-\sum_i x_i .
 \end{aligned}                                                     \tag{1}
\]
The right side defining \(F_n\) extends holomorphically to
\(\Sigma\). On a compact subset of \(\Sigma\) the real parts
\(\Re x_i\) are bounded below by a positive constant. For
fixed \(n\) the real part of the integrand exponent is
\(\sum_i[(1+n\Re x_i)u_i-n\Re x_i(e^{u_i}-1)]\),
which is bounded above by \(C-c\sum_i|u_i|\) there
for constants \(C,c>0\). This justifies holomorphy by
dominated integration.

There are analytic functions \(B_0,B_1\) on \(\Sigma\),
with \(B_0\) nonvanishing and positive on the real line,
such that uniformly on compact subsets,
\[
 F_n(s)=e^{n a(s)}n^{-(r-1)/2}
               \left(B_0(s)+\frac{B_1(s)}n+O(n^{-2})\right).
 \tag{2}
\]
Choose an orthonormal real coordinate matrix \(E\) for \(W\)
and put \(H(s)=E^{\mathsf T}\operatorname{diag}(x_i)E\).
The Hermitian real part of \(H\) is positive definite,
since for each complex vector \(\xi\ne0\),
\(\Re(\bar\xi^{\mathsf T}H\xi)=\sum_i\Re x_i |(E\xi)_i|^2>0\).
Its determinant is therefore nonzero. Use the square
root of the determinant that is positive for real \(s\).
The Gaussian integral with this quadratic form is
\((2\pi)^{(r-1)/2}(\det H)^{-1/2}\), by analytic
continuation from real \(s\). In particular
\[
 B_0(s)=\frac{\prod_i x_i}{\|d\|}
            \frac{(2\pi)^{(r-1)/2}}{\sqrt{\det H(s)}} .
\]

For (2), \(e^u-1-u\ge u^2/[2(1+|u|)]\) for real \(u\).
For \(u<0\) this follows by integrating
\(1-e^{-v}\ge v/(1+v)\); for \(u\ge0\) it follows from
the series. Thus outside a small fixed ball the integral
in (1) has an exponentially small relative tail, uniformly
on a compact parameter set. Inside the ball first put
\(u=z/\sqrt n\). Expand the phase through degree five
and the amplitude through degree three. The resulting
terms at orders \(n^{-1/2}\) and \(n^{-3/2}\) are odd in
\(z\) and have zero Gaussian integral. The remaining
error has integral \(O(n^{-2})\) after extracting the
Jacobian \(n^{-(r-1)/2}\). To see that estimate, write
the phase remainder after its quadratic part as
\(P_3(u)+P_4(u)+P_5(u)+O(\|u\|^6)\) and use the
exponential Taylor remainder through degree three.
After rescaling, all discarded terms are bounded by
\(n^{-2}\) times a fixed polynomial in \(\|z\|\) times
a Gaussian with slightly smaller quadratic coefficient.
The ball can be chosen uniformly small by the compact
lower bound on the real part of \(H\).
The \(n^{-1}\) coefficient is a finite sum of Gaussian
moments of even polynomials and is analytic. This also
covers \(r=1\) with the zero-dimensional Gaussian.

## Stationary equation and its critical points

The scalar phase is
\(\Phi_\zeta(s)=\log(1+\zeta e^s)+a(s)\)
where the first logarithm is defined. Differentiating gives
\(a'(s)=-q(s)\), since
\[
 \frac{d}{dq}\left(\sum_i a_i\log x_i-\sum_i x_i\right)
        =-q\sum_i d_i^2/x_i .
\]
Thus at a stationary point
\[
 \frac{\zeta e^s}{1+\zeta e^s}=q,\qquad
 \zeta(q)=\frac{q}{1-q}\prod_i x_i(q)^{-d_i}.
 \tag{3}
\]
For real \(q\in(0,1)\), let
\[
 \rho(q)=q(1-q)\sum_i\frac{d_i^2}{x_i(q)}
       =\sum_i\left(b_i-d_iq-\frac{a_i b_i}{x_i(q)}\right).
\]
It vanishes at both endpoints, is positive in the
interior, and is strictly concave, because
\[
 \rho''(q)=-2\sum_i\frac{a_i b_i d_i^2}{x_i(q)^3}<0.
\]
Moreover \((\log\zeta(q))'=(1-\rho(q))/(q(1-q))\).
At a stationary point,
\(\partial_s^2\Phi_\zeta=(\rho(q)-1)/s'(q)\).

There is also an exact global root count for the rational
critical equation. Write \(p_i=-a_i/d_i\), group equal
poles, and let \(m_0\) be the number of distinct poles.
Define the rational function
\[
 {\cal L}(q)=\sum_i \frac{d_i}{q-p_i}-\frac1q+\frac1{q-1}.
 \tag{4}
\]
This is minus the logarithmic derivative of (3).
The following count concerns the zeros of (4), away from
the poles. The rational map (3) can also ramify at its own
zeros or poles.

The residue at a pole left of zero is positive and at a
pole right of one it is negative. If \({\cal L}(q)=0\)
at a real \(q<0\), there must be a pole in \((q,0)\),
since otherwise every term of the grouped pole sum and
the difference \(-1/q+1/(q-1)\) are positive.
At such a zero,
\[
 (q{\cal L})'(q)=\sum_i\frac{a_i}{(q-p_i)^2}
                     -\frac1{(q-1)^2}>0
\]
because \(a_i\ge1\) at a pole in \((q,0)\).
Every zero on \(q<0\) is therefore simple with
\({\cal L}'(q)<0\). Similarly a zero at \(q>1\) requires
a pole in \((1,q)\), and
\[
 ((q-1){\cal L})'(q)=\sum_i\frac{b_i}{(q-p_i)^2}
                         -\frac1{q^2}>0 .
\]
Thus every such zero is simple with positive derivative.
There is exactly one zero in each gap between successive
poles on the same side, by the pole signs and the
direction of crossing. There is no zero between the
closest poles and the interval \([0,1]\).
If \(\Delta<0\), the far left tail has one additional zero;
if \(\Delta>0\), the far right tail has one additional
zero, whenever that side has poles. The sign tests
here use \({\cal L}(q)\sim\Delta/q\) at infinity.
For \(\Delta=0\) the leading term is \((1-A)/q^2<0\),
since both signs occur and \(A\ge2\), and neither tail
adds a zero. If all signs agree, the applicable far tail
has no extra zero by the same tests.

After clearing denominators in (4), the numerator has
degree \(m_0+1\) if \(\Delta\ne0\) and \(m_0\) otherwise.
The exterior zeros just counted number \(m_0-1\) and
\(m_0-2\), respectively. Thus the two remaining zeros
of (4), counted with multiplicity, form a nonreal conjugate
pair, a double zero inside \((0,1)\), or two simple zeros
inside \((0,1)\), according as \(\max\rho\) is below,
equal to, or above one.

## Comparing the two terminal orders

At a nondegenerate saddle, differentiating the stationary
value with respect to \(\log\zeta\) gives \(q\).
For consecutive orders at the same \(t\),
\(\zeta'=(1+1/n)^\Delta\zeta\).
If \(s_j(\zeta)\) is an analytic saddle branch and
\(\Phi_j(\zeta)=\Phi_\zeta(s_j(\zeta))\), then
\[
 (n+1)\Phi_j(\zeta')-n\Phi_j(\zeta)
          =\Phi_j(\zeta)+\Delta q_j(\zeta)+O(n^{-1}).
 \tag{5}
\]
The exponential of this leading increment is
\[
 e^{\Phi_j+\Delta q_j}=e^{-A}\Psi_a(q_j),\qquad
 \Psi_a(q)=\frac{\prod_i x_i(q)^{a_i}}{1-q},\qquad
 \Psi_b(q)=\frac{\prod_i x_i(q)^{b_i}}q,\qquad
 \zeta(q)=\frac{\Psi_a(q)}{\Psi_b(q)}.
 \tag{6}
\]
Thus once two nondegenerate saddle contributions are known
to exhaust the integral with nonzero analytic amplitudes
and uniform relative \(o(1)\) errors, simultaneous
cancellation requires limiting equality of \(\Psi_a\) at
the two saddles, and hence also of \(\Psi_b\).
Equal dominance means equal real parts of \(\Phi_j\).
Combining this with (6) forces equal \(\Re q_j\) when
\(\Delta\ne0\). The functions in this condition come from
the stationary equation and the change in order, independent
of which contour establishes the two-contribution expansion.

The order link also gives a monotonicity along individual
regular level components of the stationary parameter.
On the upper part of \({\cal V}\) choose
\[
 X(q)=s(q)+\log(1-q)-\log q,\qquad
 F(q)=\log\Psi_a(q),
\]
using the right-half-plane logarithms for \(x_i\) and
the usual branches for \(q,1-q\) there. Thus \(e^X=1/\zeta(q)\)
and
\[
 F'=\Delta-qX',\qquad
 (F-\Delta q)'=-qX'.
\]
On a connected regular component of \(\Re X=\text{constant}\),
\(\Im X\) is a strictly monotone coordinate. Its differential
is nonzero on the tangent line, by the Cauchy–Riemann equations,
and orients the component consistently. A closed component is
excluded by this same nonvanishing differential of a
single-valued function. Along such a component in the upper
half-plane,
\[
 \frac{d\Re(F-\Delta q)}{d\Im X}=\Im q>0 .
\]
Thus two points on the same regular component cannot have
equal saddle modulus \(|e^{\Phi_j}|\). For \(\Delta=0\)
this also gives strict monotonicity of \(|\Psi_a|\).
The comparison leaves pairs on distinct components
and passages through critical points.

## A local exclusion around the positive real axis

For real \(y_i>0\) and \(\zeta_0>0\), write
\[
 G_{\zeta_0}(y)=\log(y^a+\zeta_0y^b)-\sum_i y_i .
\]
Its upper level sets are compact in the positive orthant.
If \(\sum y_i\to\infty\), the logarithm grows at most
logarithmically in that sum. If \(\sum y_i\) stays
bounded and a coordinate tends to zero, both monomial
terms tend to zero. On each level set
\(s=d\cdot\log y\), the unique transverse maximum is
\(y_i=x_i(q(s))\), by strict concavity in the logarithmic
coordinates. Hence the real maxima are precisely the
global maxima of \(\Phi_{\zeta_0}\) in (1).

If \(\max\rho<1\), each positive \(\zeta_0\) has exactly
one nondegenerate real maximum, by (3). If \(\max\rho=1\),
the same holds except at the parameter \(\zeta_0\) of the
unique double critical point. If \(\max\rho>1\), there
are two turning points. Between their two parameter
values the three stationary points have signs of
\(\partial_s^2\Phi\) in the order negative, positive,
negative. The two maxima exchange dominance once:
the derivative of their phase difference with respect to
\(\log\zeta\) is \(q_+-q_->0\). At the two ends of the
interval the incipient or disappearing maximum has the
same value as the adjoining minimum. The other maximum
is strictly higher, because the real phase is strictly
monotone between the stationary points. There is
therefore one positive coexistence scale \(\zeta_*\)
with two distinct nondegenerate global maxima. Away
from it the single global maximum is nondegenerate.

Fix a positive \(\zeta_0\) where every real global
maximum is nondegenerate. The real contour yields
a uniform local saddle expansion on a small complex
disk about \(\zeta_0\). To justify localization, exclude
fixed small intervals in \(s\) around the real global
maxima. Compact upper level sets give a strict gap
\(\eta>0\) on the remaining region of the original
orthant. Its integral has the bound
\(e^{(n-1)(\max G_{\zeta_0}-\eta)}\int e^{G_{\zeta_0}(y)}dy\).
The last integral is finite. For complex \(\zeta\) nearby,
the estimate persists with a slightly smaller gap by
\[
 |y^a+\zeta y^b|
 \le (1+|\zeta-\zeta_0|/\zeta_0)(y^a+\zeta_0y^b).
\]
On each retained interval, deform inside \(\Sigma\) to
a horizontal segment through the corresponding analytic
saddle. The two short connecting segments retain
the original gap by continuity, and a sufficiently small
parameter disk keeps the shifted saddle nondegenerate.
Equation (2) and the one-dimensional Gaussian expansion
then give the corresponding local contribution as
\[
 n^{-r/2}e^{n\Phi_j(\zeta)}{\cal A}_j(\zeta)
      \left(1+\frac{{\cal B}_j(\zeta)}n+O(n^{-2})\right).
 \tag{7}
\]
The same Taylor argument as in (2) gives the error
uniformly. The amplitudes are analytic, and
\({\cal A}_j(\zeta_0)>0\) since the contour is traversed
from left to right through a real maximum. In the
two-maximum case the phase gap on the omitted region
remains bounded away from zero relative to both local
contributions after shrinking the disk. Its contribution
can be absorbed into one of their remainder terms.

A single maximum in (7) excludes zeros for all large
orders in a smaller disk. For two maxima set
\(\delta(\zeta)=\Phi_+(\zeta)-\Phi_-(\zeta)\);
then \(\delta(\zeta_*)=0\) and
\(\zeta_*\delta'(\zeta_*)=q_+(\zeta_*)-q_-(\zeta_*)>0\).
At any zero in the disk, (7) gives
\[
 e^{n\delta(\zeta)}
 =-\frac{{\cal A}_-(\zeta)}{{\cal A}_+(\zeta)}
           \left(1+\frac{{\cal B}_-(\zeta)-{\cal B}_+(\zeta)}n
                       +O(n^{-2})\right).
 \tag{8}
\]
If \(\Delta\ne0\), divide this equation by its
counterpart at order \(n+1\), parameter \(\zeta'\).
A simultaneous zero requires
\(\exp(\delta+\Delta\zeta\delta')=1+O(n^{-1})\)
by (5). At \(\zeta_*\) the left side is a positive real
number distinct from one. Shrinking the disk excludes
this possibility for large \(n\).
If \(\Delta=0\), the two parameters agree, and (8)
at the consecutive orders gives
\(e^{\delta(\zeta)}=1+O(n^{-2})\).
Near its simple zero \(\delta(\zeta_*)=0\), this
forces \(\delta(\zeta)=O(n^{-2})\) and \(\zeta\to\zeta_*\).
Then the left side of (8) tends to one while the right
side tends to a negative number. This completes the
local consecutive exclusion for both cases.

## Bounded scaled windows at the cusp

Suppose \(\max\rho=1\), attained at \(q_c\in(0,1)\).
Put \(s_c=s(q_c)\), \(\zeta_c=\zeta(q_c)\),
\(\beta_c=\log\zeta_c\) and \(K=q_c(1-q_c)\).
For every \(H<\infty\) there is \(n_0(a,b,H)\) such that
the consecutive orders have no common zero \(t\) in
\[
 \left|\log(n^\Delta t)-\beta_c\right|\le H n^{-3/4},
 \tag{9}
\]
using the branch near \(\zeta_c\). This includes
equal-degree integer cusps, for example \(a=(1,3)\)
and \(b=(3,1)\), where \(\rho\) attains one at \(1/2\).

To obtain the jet, write \(\beta(q)=\log\zeta(q)\) and
\(Q(s,\beta)=e^{\beta+s}/(1+e^{\beta+s})\).
The two log odds at the same \(s\) differ by
\(\beta-\beta(q(s))\).
At the cusp \(\beta'(q_c)=\beta''(q_c)=0\), and
\(\beta'''(q_c)=-\rho''(q_c)/K>0\).
Since \(q'(s_c)=K\), the scalar phase
\(\phi(s,\beta)=a(s)+\log(1+e^{\beta+s})\) has
\[
 \begin{aligned}
 \phi_s(s,\beta_c)&=Q(s,\beta_c)-q(s)
       =\frac{\rho''(q_c)K^3}{6}(s-s_c)^3
            +O((s-s_c)^4),\\
 \phi_{s\beta}(s_c,\beta_c)&=K .
 \end{aligned}
\]
Thus its quartic coefficient is
\(g=\rho''(q_c)K^3/24<0\).
Choose \(\ell>0\) with \(g\ell^4=-1/12\)
and put \(b_*=(K\ell)^{-1}>0\).

Let \(\epsilon=n^{-1/4}\), allowing positive real \(n\)
in the integral (1), and set
\(\beta=\beta_c+b_*\epsilon^3 h\).
After removing the nonzero factor
\(\exp(\epsilon^{-4}\phi(s_c,\beta))\epsilon^{2r-1}\)
and a fixed positive constant, the resulting function
\(J_\epsilon(h)\) has the limit
\[
 Z(h)=\int_{\mathbb R}e^{-y^4/12+hy}\,dy .
 \tag{10}
\]
The convergence is locally uniform in complex \(h\),
and \(J_\epsilon\) extends \(C^1\) to \(\epsilon=0\)
on the real half-line, locally uniformly with its
holomorphic \(h\) derivatives. Here the real-\(n\)
extension of the integral uses the analytic logarithm
near positive \(\zeta_c\) and agrees with \(I_n\) at
integer \(n\).

For these claims use (1) before integrating \(u\).
The complete real phase has a unique maximum at
\((s_c,0)\). Choose a smooth cutoff supported in a
small real neighborhood of this point and equal to
one in a smaller neighborhood. Put
\(s=s_c+\ell\epsilon y\) and \(u=\epsilon^2 v\) there.
The Jacobian has power \(\epsilon^{2r-1}\).
The scalar exponent after removing the value at
\(s_c\) expands as
\[
 \epsilon^{-4}[\phi(s,\beta)-\phi(s_c,\beta)]
          =-y^4/12+hy+o(1).
\]
For fixed \(y\) this is differentiable at \(\epsilon=0\):
the base difference in \(s\) has the factor
\((s-s_c)^4\), and the parameter-dependent difference
has the factor \((s-s_c)(\beta-\beta_c)\).
Both factors come with analytic quotients. Uniformly on
the cutoff support and bounded \(h\), the real part is at
most \(-c y^4+C|y|\); its \(\epsilon\) derivative is
bounded by a fixed polynomial in \(|y|\). The latter
follows by differentiating the two factored analytic
quotients on their compact domain.

In the transverse exponent,
\[
 \epsilon^{-4}(e^{\epsilon^2v_i}-1-\epsilon^2v_i)
   =v_i^2/2+\epsilon^2 v_i^3 E_i(\epsilon^2v_i)
\]
with analytic \(E_i\). Its real part after multiplication
by the positive real \(x_i(s)\) is bounded below by a
fixed positive multiple of \(v_i^2\) on the support.
The \(\epsilon\) derivatives of the exponent and amplitude
are bounded by polynomials in \(|y|+\|v\|\).
This gives integrable domination by a polynomial times
\(e^{-c y^4+C|y|-c\|v\|^2}\) for the scaled local integral
and its first derivative. At \(\epsilon=0\) the transverse
integral is the positive Gaussian at \(s_c\), which
supplies the fixed normalization constant in (10).
The derivative of the transformed cutoff is covered
by the same domination.

On the complement there is a strict phase gap at
\(\beta_c\) by compact upper level sets. Changing \(\beta\)
by \(O(\epsilon^3)\) changes the phase uniformly by
\(O(\epsilon^3)\) on the whole real orthant: the ratio
\(y^b/(y^a+e^\beta y^b)\) is uniformly bounded there
for \(\beta\) near \(\beta_c\). Differentiating the
normalized remainder contributes only powers of
\(\epsilon^{-1}\) and an integrable factor
\(1+|G_{\zeta_c}(y)|\) against \(e^{G_{\zeta_c}(y)}\).
For the latter integrability one can replace
\(e^{G_{\zeta_c}}\) by \(e^{G_{\zeta_c}/2}\), which is
integrable by the explicit monomial formula.
The remainder and its derivative thus vanish
exponentially. Cauchy's formula on a larger compact
\(h\)-disk gives the same convergence for the
\(h\) derivatives.

The zeros of \(Z\) are imaginary and simple. The
imaginary-zero assertion, also with a real quadratic
term \(\tau u^2/2\) in the exponent, is the empty-product
case of Newman's Theorem 2, taking the quartic density
with quadratic coefficient \(\tau/2+1\) and the Gaussian
multiplier \(e^{-u^2}\). [Newman, Fourier Transforms with Only Real Zeros](https://sites.math.northwestern.edu/~auffing/papers/Newman.pdf)
For simplicity write \(f(x,\tau)=Z(ix,\tau)\) with
that quadratic term. Integration by parts gives
\(f_{xxx}+3\tau f_x-3x f=0\) and
\(f_\tau=-f_{xx}/2\). A triple zero would force zero
initial data in the third-order equation. At a double
zero with \(c=f_{xx}/2\ne0\),
\(\eta^{-2}f(x_0+\eta z,\tau_0-\eta^2)\to c(z^2+1)\),
locally uniformly. Rouché's theorem would produce
nonreal \(x\)-zeros at real \(\tau_0-\eta^2\), a
contradiction. Also \(Z(0)>0\).

In any bounded \(h\)-window, enlarge the boundary to a
compact disk avoiding model zeros. Rouché's theorem and
the simplicity just proved put all zeros into disjoint
small disks about the finitely many relevant model zeros
\(h_j\), one per disk. The real \(C^1\) implicit-function
theorem gives branches \(h_j(\epsilon)=h_j+O(\epsilon)\)
with bounded derivatives, using the \(C^1\) extension
across zero by its first Taylor jet.
The corresponding log-parameter branches are
\[
 \beta_j(n)=\beta_c+b_* n^{-3/4}h_j(n^{-1/4}),\qquad
 \beta_j'(n)=-\tfrac34 b_*h_j n^{-7/4}+O(n^{-2}).
 \tag{11}
\]
If the same \(t\) is a zero at two consecutive orders
within (9), then the change in the \(h\)-coordinate
is \(O(n^{-1/4})\), since the exact change in \(\beta\)
at fixed \(t\) is \(\Delta\log(1+1/n)\).
Thus the two zeros must lie on the same branch.
For \(\Delta\ne0\) the exact change is of order \(n^{-1}\),
whereas (11) gives order \(n^{-7/4}\).
For \(\Delta=0\) the exact change is zero, whereas (11)
gives a nonzero consecutive change since \(h_j\ne0\).
This proves (9).

## Local exclusion around every negative real scale

Fix a negative real scale \(\zeta_0=-u\), \(u>0\). For
all sufficiently large \(n\), consecutive orders have no common
zero \(t\) with \(n^\Delta t\) in a fixed complex neighborhood of
\(\zeta_0\). The neighborhood and terminal threshold depend on
the fixed exponent profile and on \(\zeta_0\).

To prove this, split the real conditional contour at
\(s_0=-\log u\). Use the real phases
\[
 \phi_L(s)=a(s)+\log(1-u e^s)\quad(s<s_0),\qquad
 \phi_R(s)=a(s)+\log(u e^s-1)\quad(s>s_0).
\]
The critical equation on both intervals is \(Q(s)=q(s)\),
where \(Q(s)=-u e^s/(1-u e^s)\) and
\(Q'(s)=Q(s)(1-Q(s))<0\). Since \(q'(s)>0\), the
second derivative of each phase is strictly negative
throughout its interval. Each phase attains a maximum in the
interior. Indeed, its exponential is
\(|y^a-u y^b|e^{-\sum_i y_i}\) at \(y_i=x_i(q(s))\).
At \(s=s_0\) this exponential tends to zero. As \(s\) tends to
either infinite endpoint, \(y\) leaves every compact set in the
positive orthant since \(d\cdot\log y=s\), and the exponential
again tends to zero. Here all exponents are positive and the
exponential dominates the polynomial at infinity. Thus there
is exactly one nondegenerate maximum on each interval.

Near either of these points we can use an analytic branch of
the corresponding logarithm in complex \(s,\zeta\), using
\(1+\zeta e^s\) on the left and \(-1-\zeta e^s\) on the right.
They define analytic stationary values
\(\Phi_L,\Phi_R\), real at \(\zeta_0\). The same local contour
deformation and Gaussian expansion as in (7) apply. On the right
the integral of the original integer power has the extra factor
\((-1)^n\). The required localization holds even across the cut
in the real phases: set
\[
 K_\zeta(y)=|y^a+\zeta y^b|\,e^{-\sum_i y_i}.
\]
On the complement of the retained saddle intervals it has a strict
gap below its global maximum at \(\zeta_0\). The gap persists near
\(\zeta_0\) by the uniform estimate
\[
 |K_\zeta(y)-K_{\zeta_0}(y)|
  \le|\zeta-\zeta_0|\,y^b e^{-\sum_i y_i}.
\]
The right-hand factor is bounded. Also \(\int K_\zeta(y)\,dy\)
is uniformly bounded for \(\zeta\) in a compact disk, so
bounding \(K_\zeta^{n-1}\) by the gap controls the omitted
integral. In particular one dominant real maximum gives a
zero-free disk for sufficiently large \(n\).

When the two real maxima have equal height, both saddle
expansions (7) hold with uniform relative errors and the sign
\((-1)^n\) on the right. Write
\(D(\zeta)=\Phi_R(\zeta)-\Phi_L(\zeta)\). A zero in a sufficiently
small disk must satisfy
\[
 (-1)^n e^{nD(\zeta)}
 =-\frac{{\cal A}_L(\zeta)}{{\cal A}_R(\zeta)}
                  (1+O(n^{-1})).
\]
The two amplitudes remain nonzero on that disk. For a common
zero divide the relation at \(n+1,\zeta'=(1+1/n)^\Delta\zeta\)
by the relation at \(n,\zeta\). Formula (5) gives
\[
 -\exp\bigl(D(\zeta)+\Delta\zeta D'(\zeta)\bigr)=1+O(n^{-1}).
\]
At \(\zeta_0\) both \(D\) and \(D'\) are real, since the
stationary values are real along the nearby negative axis.
The left side is a strictly negative real number. A fixed
smaller disk excludes this equality at all sufficiently large
orders. In particular no term of relative order \(1/n\) needs
to be matched between the two expansions here.

## A fixed annulus for possible common roots

At \(\zeta_0=0\) the scalar phase is \(a(s)\), which has a
unique nondegenerate real maximum where \(q(s)=0\), since
\(a'(s)=-q(s)\) and \(a''(s)=-q'(s)<0\). Equivalently the
original real phase is maximized at \(y_i=a_i\). The uniform
bound with \(K_\zeta\) in the preceding section gives the
strict localization gap for all sufficiently small
\(|\zeta|\). The one-saddle version of (7) therefore gives a
fixed zero-free disk about zero for all sufficiently large
orders. This argument applies to the reversed profile as
well. For any \(t\ne0\),
\[
 M_n^{a,b}(t)=t^n M_n^{b,a}(1/t),\qquad
 n^{-\Delta}/t=1/(n^\Delta t).
\]
Consequently the individual scaled polynomials have all roots
in a common fixed annulus \(r_*\le |n^\Delta t|\le R_*\)
for all sufficiently large orders. The radii depend only on
the fixed profile.

A finite subcover of the negative real part of the annulus
then gives an angular neighborhood of the negative ray that
excludes common roots at all sufficiently large consecutive
orders. If \(\max\rho\ne1\), the same is true of the positive
ray by the positive real exclusion. If \(\max\rho=1\),
every possible sequence of common roots approaching the
positive ray must approach the unique cusp scale. The bounded
quartic window excludes a fixed bounded rescaling there.

## Nearby-prime Newton polygons

For the ordinary factorial coefficients there is also a
sparse exact arithmetic reduction. Write \(n=p+u\) with
\(p\) prime and \(u\ge1\), and suppose
\[
 E u<p,\qquad En<p^2,\qquad E=\max_i\max(a_i,b_i).
\]
For both orders \(m=p+w\), \(w=u,u-1\), the coefficient of
\(t^k\) in the unnormalized \(M_m(t)\) has valuation
\[
 v_p([t^k]M_m)=\mathbf1_{\{w<k<p\}}+
                \sum_i\left\lfloor\frac{a_i m+d_i k}{p}\right\rfloor .
 \tag{12}
\]
This is Legendre's formula with all arguments below \(p^2\).
The binomial valuation follows since
\(v_p(m!)=1\), \(v_p(k!)=\mathbf1_{\{k\ge p\}}\),
and \(v_p((m-k)!)=\mathbf1_{\{k\le w\}}\).

Each floor in (12) is monotone in \(k\) with net change
\(d_i\). The lower hull uses only the two endpoints of
each constant plateau, since the other points of a plateau
lie on their connecting segment. There are at most
\(1+2+\sum_i|d_i|\) plateaus, so the hull can be built
from a list whose size is independent of \(n\).
A nonhorizontal edge with rise \(h\) and width \(b\) has
at most \(1+\gcd(|h|,b)\le1+|h|\) eligible residual
coefficient indices. For each such index compute the
unit of each factorial by
\[
 p^{-v_p(l!)}l!\equiv
 (-1)^{\lfloor l/p\rfloor}\lfloor l/p\rfloor!\,(l\bmod p)!
 \pmod p,\qquad l<p^2,
\]
which follows from Wilson's theorem. The binomial unit
is the ratio of the unit of \(m!\) to the units of
\(k!\) and \((m-k)!\).
Matching lower slopes across the two orders, the degree
of the gcd of the corresponding residual polynomials
bounds the width of a common factor in units of the
reduced slope denominator. To justify this, at slope
\(j/e\) substitute \(t=\pi^{-j}Y\) with \(\pi^e=p\),
remove the minimum valuation and reduce. Multiplicativity
passes to the edge polynomials in \(Y^e\) after their
lowest monomials are removed. A positive common degree
must be even, because the even-order moment is strictly
positive on the real line.

The initial plateau \(0\le k\le w\) has valuation \(A\),
and the terminal plateau \(p\le k\le p+w\) has valuation
\(\sum b_i\). Their units, after fixed nonzero factors,
are copies of the coefficients of \(M_w(t)\) and
\(t^p M_w(t)\), respectively. This follows from Lucas's
congruence for the binomial coefficient and the displayed
factorial-unit formula. If all \(d_i>0\), the initial
plateau is the full zero-slope residual, since each
interior ordinate exceeds \(A\) and the terminal plateau
has larger valuation. This is the small-order residual
transfer used for the constant boundary.

With mixed signs that transfer can fail. For
\(a=(1,2,2), b=(4,1,1)\), \(m=p+1\) and any prime \(p>10\),
the initial plateau has valuation five. At \(k=3\) the
factorial arguments are \(p+10,2p-1,2p-1\), and (12)
gives valuation four. Thus the minimum is below the
initial plateau despite \(\sum b_i>\sum a_i\). A horizontal
edge at the minimum can require the interior coefficients
instead. The sparse hull and the nonhorizontal residual
tests remain valid for these mixed profiles.

The neighborhoods, scaled windows, and thresholds above
depend on the fixed profile and the chosen real scale. The global
all-order problem still requires control of the complex
saddle contributions away from these neighborhoods,
including transitions where the dominant saddles
coalesce or where other contributions have the same
exponential order. A contour-independent criterion
such as (6) restricts consecutive cancellation once the
contributors are identified.
