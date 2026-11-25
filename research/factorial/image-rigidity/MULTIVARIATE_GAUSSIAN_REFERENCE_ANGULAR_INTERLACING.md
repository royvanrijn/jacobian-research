# Consecutive angular interlacing for the quadratic Gaussian reference

For \(\tau>0\) and \(j\ge1\) put
\[
 R_j(z;\tau)=\sum_{k=0}^j\binom jk e^{-\tau k(j-k)}z^k .
\]
These are the endpoint-normalized Gaussian reference polynomials
for factorial binomials. At degree \(j\ge2\), an outer loss
\(-\log u_j\) corresponds to \(\tau=(-\log u_j)/(j-1)\).

**Theorem.** At equal \(\tau>0\), the roots of \(R_n\) and
\(R_{n+1}\) strictly interlace in argument, for every \(n\ge1\).
Write \(x_i(\tau)\) and \(y_i(\tau)\) for the positive angular
displacements of their roots from the negative real axis,
in decreasing order. For every \(0<\tau'\le\tau\),
\[
                       x_i(\tau)>y_{i+1}(\tau')
\]
whenever both heights exist. Thus the inner gaps stay strictly separated
also at the two physical Gaussian reference temperatures of an ordinary
factorial profile. The outer gaps \(y_i(\tau')-x_i(\tau)\) at distinct
physical temperatures and the transfer from the reference to the exact
factorial coefficients are separate obligations.

## Root motion on the upper semicircle

These polynomials are the unitary Hermite polynomials under
\(z\mapsto-z\) and \(\sigma^2=2\tau\). The circle location also
follows from the Lee–Yang theorem. [Kabluchko, Lemma 2.1](https://arxiv.org/pdf/2203.05533v3)

Define
\[
 f_j(x,\tau)=(2i)^{-j}(-1)^j e^{\tau j^2/4}e^{-ijx/2}
                                      R_j(-e^{ix};\tau).
\]
Its Fourier multiplier gives \(\partial_\tau f_j=-\partial_x^2 f_j\),
with \(f_j(x,0)=\sin^j(x/2)\). At a simple root \(x\) its
speed is \(f_j''(x,\tau)/f_j'(x,\tau)\).
The circle factorization of \(R_j\) writes \(f_j\) as a
nonzero constant times \(\prod_\alpha\sin((x-\theta_\alpha)/2)\),
where \(\theta_\alpha\) are the signed heights. Consequently the
speed at a positive root \(x=x_i\) is
\[
 \dot x=\cot x+\epsilon\cot(x/2)+
                     \sum_{k\ne i}F(x,x_k),\qquad
 F(x,a)=\frac{2\sin x}{\cos a-\cos x}.                              \tag{1}
\]
Here \(\epsilon=1\) in odd degree and zero in even degree.
The self-reflected root contributes \(\cot x\).
The paired force \(F(x,a)\) is strictly increasing as a function
of \(a\) on each of \([0,x)\) and \((x,\pi)\).
It is positive on the former and negative on the latter, and
\(F(x,0)=2\cot(x/2)\).

For each fixed \(j\), the heights divided by \(\sqrt{2\tau}\)
converge to the roots of the probabilists' Hermite polynomial
\(\mathrm{He}_j\) as \(\tau\downarrow0\). In fact the coefficient of
\(x^p\tau^m\) in the Taylor series of \(f_j\) is
\((-1)^m f_j^{(p+2m)}(0,0)/(p!m!)\), with superscript denoting
the \(x\) derivative. The initial sine power vanishes to order \(j\)
and begins with \(2^{-j}x^j\). Upon substituting
\(x=\sqrt{2\tau}\,v\), the terms \(p+2m=j\) give
\(2^{-j}(2\tau)^{j/2}\mathrm{He}_j(v)\), and higher terms are
\(o(\tau^{j/2})\) locally uniformly in \(v\).
The scaled function is real analytic in \(v\) and \(s=\sqrt{2\tau}\)
after division by \(2^{-j}s^j\), including at \(s=0\).
The real analytic implicit-function theorem gives one simple real
height branch per Hermite root. The Hermite roots strictly
interlace by their Jacobi recurrence, which gives strict
interlacing at sufficiently small positive \(\tau\).

As long as all heights stay real and distinct, the velocities of the
signed heights are \(\dot\theta_\alpha=\sum_{\beta\ne\alpha}
\cot((\theta_\alpha-\theta_\beta)/2)\), by the same factorization.
This is twice the gradient of
\[
        \mathcal E=\sum_{\alpha<\beta}
                    \log|\sin((\theta_\alpha-\theta_\beta)/2)|.
\]
Thus \(\mathcal E\) is nondecreasing. Its summands are nonpositive,
and a collision would send it to \(-\infty\).
A simple real zero of the real analytic trigonometric polynomial
stays real locally, so all heights remain real and simple at every
positive \(\tau\). Positivity of \(R_j(1;\tau)\) excludes heights
\(\pm\pi\). Real reciprocity pairs the positive heights with their
negatives, with the additional simple height zero in odd degree.

## First contact across adjacent degrees

Consider a first contact at a positive height \(x\).
The same-set heights remain distinct. At a contact \(y_i=x_i=x\)
each other same-index pair satisfies \(y_k\ge x_k\) on the same
side of \(x\). By (1) their force differences are nonnegative.
For even \(n\) the extra center in degree \(n+1\) contributes
\(\cot(x/2)>0\). For odd \(n\) the extra smallest positive height
in degree \(n+1\) supplies \(F(x,y_{\lfloor n/2\rfloor+1})\),
which exceeds \(2\cot(x/2)\) and outweighs the lost central force
\(\cot(x/2)\). Thus \(\dot y_i-\dot x_i>0\).

At a contact \(y_{i+1}=x_i=x\), pair \(y_{k+1}\) with \(x_k\)
where both exist. Now \(y_{k+1}\le x_k\) on the same side of \(x\),
so the paired force differences are nonpositive. The extra
largest height \(y_1>x\) contributes negatively.
For odd \(n\) the lost central force subtracts \(\cot(x/2)\).
For even \(n\), the smallest positive height
\(x_{\lfloor n/2\rfloor}<x\) has no pair in this comparison,
so its force of more than \(2\cot(x/2)\) is subtracted and the
new center adds only \(\cot(x/2)\).
Hence \(\dot y_{i+1}-\dot x_i<0\). Each sign contradicts being
a first loss of the indicated strict inequality. Symmetry and the
unique structural zero at the odd degree give full interlacing.

Each positive height moves strictly outward under its own heat flow.
Indeed its analytic Hermite branch has the form \(s(v_0+O(s))\)
with \(v_0>0\) and \(s=\sqrt{2\tau}\), so its speed starts positive.
Differentiate (1) in \(\tau\). The speed vector satisfies a linear
system whose off-diagonal entries are \(\partial_a F(x,a)>0\).
On any compact interval of positive times all the entries are
continuous. The integrating-factor expression in each component
keeps that component strictly positive up to any proposed first
zero, excluding such a zero. Hence \(y_{i+1}(\tau')\le
y_{i+1}(\tau)<x_i(\tau)\) for \(\tau'\le\tau\).
This is the complete-graph Lee–Yang monotonicity result under
\(\tau=2t\), since \(k\) minus spins give interaction sum
\(\binom j2-2k(j-k)\). [Hou, Jiang and Newman, Theorem 2](https://arxiv.org/pdf/2208.00917v1)

For the ordinary factorial reference pair with lower and upper
outer losses \(U,W\), the times are \(\tau=U/(n-1)\) and
\(\tau'=W/n\). The
[outer-ratio inequality](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md)
gives \(\tau'/\tau>n/(n+1)\). Each matched reciprocal satisfies
\(x'\le x\), and concavity gives
\(\log(1+nx')\le n\log(1+(n-1)x)/(n-1)\), so \(\tau'\le\tau\).


## Exact finite controls

With \(q\in\{1/4,3/4,9/10,49/50,999/1000\}\) and \(2\le n\le18\),
use \(e^{-\tau}=q\) at both degrees for the equal-time controls.
For the lower physical time ratio use \(e^{-\tau}=q^{n+1}\)
and \(e^{-\tau'}=q^n\). After removing the odd degree's factor
\(1+z\), express each reciprocal polynomial in \(z+z^{-1}\).
Exact rational root isolation gave 405 strict outer comparisons
and 360 strict inner comparisons across 85 pairs at each of these
settings. The physical lower-endpoint controls test the stronger
outer inequality at those inputs.

## The first-contact sign at unequal heat parameters

The sign argument at equal heat parameters does not extend to arbitrary
height configurations under the physical lower clock ratio \(n/(n+1)\).
For \(n=4\), take a boundary configuration of the interlacing inequalities
with
\[
 x_1=y_1=\frac{\pi}{2},\qquad
 x_2=\arccos(1/10),\qquad y_2=\arccos(99/1000).
\]
The other inequalities are strict. At the contact, the degree-four speed
in (1) is \(20\), and the degree-five speed at its own time is
\(1+2000/99\). When the upper clock runs at \(4/5\) of the lower clock,
the speed of \(y_1-x_1\) is
\[
           \frac45(1+2000/99)-20=-\frac{1504}{495}.
\]
This boundary point gives an outward crossing for the cone of arbitrary
height configurations. The Hermite-initialized trajectories impose
additional relations between the heights.
