# An obstruction to uniform Gaussian half logarithms near critical coupling

For \(a>0\) and \(z\ne0\) use
\[
 \begin{split}
 G_a(t,z)&=e^{-a/4-t^2/(4a)}(e^{-t/2}+ze^{t/2}),\\
 I_j(z;a)&=\sqrt{\frac{j}{4\pi a}}\int_{-\infty}^0 G_a(t,z)^j\,dt .
 \end{split}
\]
Then \(I_j\) is a polynomial with positive coefficients and
\(B_j(z;2a)=I_j(z;a)+z^j I_j(1/z;a)\), where \(B_j\) is the
[Gaussian reference](MULTIVARIATE_DIFFUSE_LOW_TEMPERATURE_INTERLACING.md).
We distinguish this half polynomial, obtained by ending the integral
at the real midpoint, from the analytic moving-endpoint halves in
the supercritical saddle proof.

**Theorem.**

1. At each fixed \(0<a<2\), the canonical \(I_j\) has zeros
   approaching \(z=1\) with \(\log z\) on scale \(j^{-1/2}\).
   This also holds on scale \(j^{-3/4}\) when
   \(a=a_j=2+\nu_j/\sqrt j\), \(\nu_j\to\nu\in\mathbb R\).
2. Fix \(0<a<2\). For any nonnegative-coefficient polynomial
   split \(B_j(z;2a)=L_j(z)+z^j L_j(1/z)\), suppose \(L_j\)
   is zero-free on a fixed annulus about \(|z|=1\) along an
   infinite subsequence. The fraction of its zeros inside the
   inner boundary, relative to \(j\), must tend to \(1/2\).
   In particular no such sequence of splits has single-valued
   analytic half logarithms on a degree-uniform annulus.
3. There is \(\varepsilon>0\) such that for
   every fixed \(2<a<2+\varepsilon\), the canonical halves \(I_j\)
   have zeros in every fixed annulus about the circle for all
   sufficiently large \(j\).

These obstructions concern the analytic halves used to transfer
coefficient multipliers. The full Gaussian polynomials may still
interlace in these regimes.

## Scaling limits of the canonical half

For \(t\) real and \(v\) complex near zero,
\[
 \frac{G_a(t,e^v)}{G_a(0,e^v)}
   = e^{-t^2/(4a)}
        \frac{\cosh((t+v)/2)}{\cosh(v/2)} .
\]
The addition formula gives a bound on the modulus of the
cosh quotient by \(\cosh(t/2)e^{C|v||t|}\), with a fixed
constant near \(v=0\). Indeed it is
\(\cosh(t/2)|1+\tanh(v/2)\tanh(t/2)|\), and
\(|\tanh(t/2)|\le |t|/2\) for real \(t\).

At fixed \(0<a<2\), set \(c=(2-a)/(8a)>0\),
\(t=-y/\sqrt j\), and \(v=w/\sqrt j\). The normalized half has
locally uniform limit in complex \(w\)
\[
 \frac{I_j(e^{w/\sqrt j};a)}{G_a(0,e^{w/\sqrt j})^j}
       \longrightarrow
 \frac1{\sqrt{4\pi a}}\int_0^\infty e^{-cy^2-wy/4}\,dy.       \tag{1}
\]
The quadratic Taylor term of \(\log\cosh(\cdot/2)\) gives the
exponent. For \(|w|\le M\), the inequality
\(\log\cosh(t/2)\le t^2/8\) and the perturbation bound above
majorize the transformed integrand by \(e^{-cy^2+CM y}\).
This proves local uniform convergence. Completing the square
identifies the limit with a nonvanishing factor times
\(\operatorname{erfc}(w/(8\sqrt c))\).
[DLMF, complementary error function zeros](https://dlmf.nist.gov/7.13#ii)

For \(a_j=2+\nu_j/\sqrt j\) with \(\nu_j\to\nu\),
set \(t=-y/j^{1/4}\), \(v=w/j^{3/4}\), and use
\(\log\cosh(t/2)=t^2/8-t^4/192+O(t^6)\).
Then locally uniformly in complex \(w\),
\[
 j^{-1/4}\frac{I_j(e^{w/j^{3/4}};a_j)}
                   {G_{a_j}(0,e^{w/j^{3/4}})^j}
 \longrightarrow
 \frac1{\sqrt{8\pi}}\int_0^\infty
            e^{-y^4/192+\nu y^2/16-wy/4}\,dy .                \tag{2}
\]
The quadratic correction is
\(\sqrt j(1/8-1/(4a_j))\to\nu/16\).
To dominate the integrals, strictness of \(\tanh(t/2)<t/2\)
for \(t>0\), the quartic expansion at zero, and the large-\(|t|\)
growth of \(\cosh\) give a constant \(c_0>0\) with
\(\log\cosh(t/2)-t^2/8\le-c_0\min(t^4,t^2)\).
The scaled real exponent is bounded above by
\[
 -c_0\min(y^4,\sqrt j\,y^2)+C_\nu y^2+CM y .
\]
For \(y\le j^{1/4}\) the quartic controls the quadratic and
linear terms. For \(y\ge j^{1/4}\) and sufficiently large \(j\),
the expression is at most \(-y^2+CM y\).
These two bounds give a fixed integrable majorant on every
bounded complex \(w\)-set.

Both limit functions have zeros. To see this directly, omit their
positive constants and write the integral as
\(L(w)=\int_0^\infty \exp(-p(y)-wy/4)\,dy\), with \(p(y)=cy^2\)
or \(p(y)=y^4/192-\nu y^2/16\).
This entire function is bounded in modulus by
\(C\exp(C|w|^\rho)\) with \(\rho=2\) or \(4/3\), respectively,
by absorbing the linear term into half the leading power of \(p\).
If \(L\) were zero-free, it would have an entire logarithm \(P\).
For each integer \(m\ge1\), the coefficient of \(w^m\) in \(P\)
is twice the \(m\)-th Fourier coefficient of
\(\Re P(Re^{i\theta})\), divided by \(R^m\).
Subtract a constant upper bound \(C+C R^\rho\)
from this real part. The mean of the real part is \(\Re P(0)\),
so that coefficient has modulus \(O(R^{\rho-m})\).
Thus \(P\) would be a polynomial.
On the positive real line, substituting \(y=4t/w\) and using
dominated convergence gives \(L(w)\sim4/w\).
This implies \(\Re P(w)=-\log w+O(1)\), impossible for a
polynomial on that line. Finally, locally uniform convergence
and Rouché on a small circle about a limit zero give a sequence
of zeros with scaled logarithms tending to that zero, in both (1)
and (2). The denominators in (1)–(2) are nonzero near \(z=1\).


## Obstruction for arbitrary subcritical polynomial halves

Fix \(0<a<2\). The full polynomial \(B_j(z;2a)\) is
zero-free on a fixed complex disk about \(z=1\) for all
sufficiently large \(j\). Indeed at \(z=1\) the full integral
has a unique strict real maximum at \(t=0\), by
\(\log\cosh(t/2)\le t^2/8\), and its second logarithmic
derivative there is \(-2c<0\).
The saddle continues analytically to a neighborhood of \(z=1\).
Shift the real contour horizontally through this saddle.
The strict modulus loss away from a fixed local interval persists
by continuity on bounded pieces, and the Gaussian term gives
uniform decay on the tails. In the local interval Taylor expansion
has strictly negative real quadratic term and a uniformly bounded
cubic remainder. Gaussian integration after removing the stationary
exponential gives a nonzero analytic amplitude and an
\(O(j^{-1/2})\) error, uniformly on a smaller disk about \(1\).
This proves the zero-free assertion.

Now suppose \(L_j\) is a nonnegative-coefficient polynomial half
with no zeros in a fixed annulus \(r<|z|<R\) about \(|z|=1\),
along an infinite subsequence. Let \(m_j\) count the zeros
inside the inner boundary. Then \(z^{-m_j}L_j(z)\) admits an
analytic logarithm there by the argument principle. Put
\(F_j=j^{-1}\log(z^{-m_j}L_j)\), using the branch real at \(z=1\).
The coefficient bound \(0\le [z^k]L_j\le[z^k]B_j\) gives
\(\Re F_j\le\log(1+R)+|\log r|\), and
\(L_j(1)=B_j(1)/2\ge1\) bounds \(\Re F_j(1)\) below.
Harnack's inequality, applied to a constant minus \(\Re F_j\),
bounds the real parts locally uniformly from below.
Local Cauchy estimates for the harmonic real parts then bound
the derivatives of \(F_j\) on compact subsets. Together with
the normalization of the imaginary part at \(1\), this makes
\(F_j\) a normal family.

Pass to a convergent further subsequence \(m_j/j\to m\) and
\(F_j\to F\). For \(\theta\) in the fixed positive angular gap
of \(B_j\), the exact symmetric half decomposition gives
\[
             |(j-2m_j)\theta-2j\,\Im F_j(e^{i\theta})|<\pi .
\]
Indeed the quantity is continuous, vanishes at zero, and
an odd multiple of \(\pi\) at a point would give a root there.
Divide by \(j\). The resulting equality
\((1-2m)\theta=2\Im F(e^{i\theta})\) on the gap extends to
all real \(\theta\) by real analyticity. Periodicity then
forces \(m=1/2\). Every convergent subsequence of \(m_j/j\)
has this limit. A single-valued analytic logarithm of \(L_j\)
on the fixed annulus would require \(m_j=0\), a contradiction.

## Canonical halves at fixed \(a>2\) near critical

Use the [supercritical left saddle construction](MULTIVARIATE_DIFFUSE_LOW_TEMPERATURE_INTERLACING.md#uniform-eventual-interlacing-for-sge4delta).
At \(z=e^{i\theta}\) the saddle has \(t_L=-r+iy\), where
\(u=\theta+y\) and
\[
   a=\frac{r(\cosh r+\cos u)}{\sinh r},\qquad
   y=\frac r{\sinh r}\sin u .
\]
The real logarithmic saddle-over-endpoint difference
\(\Xi(\theta)=\log|G_a(t_L,e^{i\theta})/G_a(0,e^{i\theta})|\)
is real analytic for \(0<\theta<\pi\).
Stationarity gives
\[
              \Xi'(\theta)=\tfrac12\tan(\theta/2)-y/(2a).
\]
At \(a=2\), the same construction gives the simple negative-real-part
saddle for all positive \(\theta\), with \(r\to0\) as \(u\to0\).
Then \(y/u\to1\) and \(\theta/u\to0\).
Thus \(\Xi(0)=0\) by continuity, and \(\Xi'<0\) on a short
positive interval at \(a=2\).
Fix one \(u\) in that interval. The equations above imply
that \(\Xi<0\) at its associated angle for all \(a>2\)
sufficiently close to \(2\).
On the other hand \(\Xi(0)>0\) for every \(a>2\): the left
saddle is the strict real half-line maximum, and its
modulus exceeds that at the midpoint.
There are therefore two nonempty open angular arcs of opposite sign.

Deform the canonical half integral to the horizontal line
through \(t_L\), ending at \(iy\), and then the connector
from \(iy\) to zero. On the horizontal half-line the modulus
has its unique strict maximum at \(t_L\), with negative
logarithmic curvature, as in the cited construction.
Thus its logarithmic rate is \(\log|G_a(t_L,z)|\) by
the local Gaussian estimate, which has a nonzero amplitude.
On the connector the phase of \(G_a(iv,e^{i\theta})\) is
\(\theta/2\) for \(0\le v\le y\) and \(0<\theta<\pi\).
The derivative of its log modulus there is
\[
 \frac{v}{2a}-\tfrac12\tan((\theta+v)/2)<0.
\]
For \(a\ge2\) this follows from
\(\tan((\theta+v)/2)\ge(\theta+v)/2\).
The connector integral has no cancellation and its rate is
\(\log |G_a(0,z)|\), by the strict endpoint maximum.
The two rates are unequal on the two chosen arcs.
Consequently \(j^{-1}\log|I_j(e^{i\theta};a)|\)
converges to the saddle rate on the arc with \(\Xi>0\)
and to the endpoint rate on the arc with \(\Xi<0\).

If there were a fixed annulus and an infinite zero-free
subsequence, the normal-family argument above would give a
real-analytic limit of \(j^{-1}\log|I_j(e^{i\theta};a)|\)
around the circle. It would match the saddle rate on one
open arc and the endpoint rate on another. Both rates are
real analytic on \(0<\theta<\pi\). Equality to the first
rate on a nonempty open arc forces equality to it everywhere
on that interval, contradicting the second arc.
This proves the fixed supercritical obstruction.

