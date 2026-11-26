# Gaussian root comparison at physical heat times

Let
\[
R_j(z;\tau)=\sum_{k=0}^j\binom jk e^{-\tau k(j-k)}z^k
\]
be the Gaussian reference for an order-\(j\) factorial moment. Its roots
are simple and on the unit circle for \(\tau>0\). Write
\(x_{j,i}(\tau)>0\) for the positive angular heights from \(-1\), in
decreasing order, and put
\(\lambda_{j,i}(\tau)=\sin^2(x_{j,i}(\tau)/2)\).
An odd degree has the additional structural height zero.
The [equal-clock theorem](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md)
gives strict consecutive interlacing and gives the inner gaps whenever the
upper clock is slower. This note gives a root-motion bound, an explicit
late-clock outer-gap theorem, and an exact first-contact test for the
remaining physical times.

## The squared-height equation

For \(j=2m+\epsilon\), \(\epsilon\in\{0,1\}\), the monic polynomial in
\(s=\sin(x/2)\) for the centered Gaussian roots is
\[
 P_j(s;\tau)=s^\epsilon
              \prod_{i=1}^m(s^2-\lambda_{j,i}(\tau)).
\]
The centered Fourier heat equation gives
\[
 \partial_\tau P_j
   =\tfrac14\big((s^2-1)\partial_s^2+s\partial_s-j^2\big)P_j,
 \qquad P_j(s;0)=s^j .
\]
In particular, if \(b_{j,r}\) is the elementary symmetric function of
order \(r\) in the squared heights, then
\[
 \dot b_{j,r}=-r(j-r)b_{j,r}
       +\frac{(j-2r+2)(j-2r+1)}4 b_{j,r-1}.                  \tag{1}
\]
Here \(b_{j,0}=1\) and \(b_{j,r}(0)=0\) for \(r>0\). Thus each squared-height
power sum has an exact description by a triangular linear equation
and Newton's identities. In particular,
\[
 \sum_i\lambda_{j,i}(\tau)=\frac j4(1-e^{-(j-1)\tau}).         \tag{2}
\]

At a positive root, differentiation and factorization give
\[
 \dot\lambda_i=F_{j,i}(\boldsymbol\lambda)
  =a-b\lambda_i+2\lambda_i(1-\lambda_i)
                         \sum_{k\ne i}\frac1{\lambda_i-\lambda_k},
 \qquad a=\frac{1+2\epsilon}{2},\quad b=1+\epsilon .          \tag{3}
\]
For example, at a root \(s_i\) the ratio
\(P_j''(s_i)/P_j'(s_i)\) is
\((1+2\epsilon)/s_i+4s_i\sum_{k\ne i}(s_i^2-s_k^2)^{-1}\).
Equation (3) is cooperative on the chamber of ordered distinct roots
in \((0,1)\), since every off-diagonal derivative is positive.

The circle location and the angular heat dynamics follow from
[Kabluchko, Section 2](https://arxiv.org/abs/2203.05533) and
[Hou, Jiang and Newman, Theorem 2](https://arxiv.org/abs/2208.00917).
The comparison below follows directly from (3).

## Motion and comparison with the equilibrium roots

Every squared height satisfies
\[
 0<\dot\lambda_i(\tau)\le
       \frac{\lambda_i(\tau)}{2(e^{\tau/2}-1)}.              \tag{4}
\]
Consequently \(\lambda_i(\tau)/(1-e^{-\tau/2})\) is nonincreasing.
For the proof, the vector \(F\) is strictly positive along the physical
trajectory by the cooperative heat-flow argument. Let \(J\) denote its
Jacobian. Differentiating (3) under a common scaling of the coordinates
and using \(b-a=1/2\) gives
\[
 (J\boldsymbol\lambda)_i-F_i
       =-\frac{F_i+\lambda_i/2}{1-\lambda_i}.
\]
At a maximum \(M=\max_i F_i/\lambda_i\), cooperativity gives
\((JF)_i\le M(J\boldsymbol\lambda)_i\) in every maximizing coordinate.
Its upper Dini derivative is therefore at most \(-M(M+1/2)\).
The Hermite limit at time zero has \(\tau M\to1\).
Integrating the inequality for \(1/M\) proves (4).

There is also an exponential barrier depending on the root index. Put
\[
 \mu_i=\cos^2\frac{(2i-1)\pi}{2j},\qquad
 e=e^{-(j-1)\tau},\qquad
 B_i(u)=\frac{\mu_i(1-u)(1-e)}{1-u+e(u-\mu_i)}.
\]
For \(j\ge2\) and every finite positive time,
\[
 B_i(\max_k\mu_k)\ \le\ \lambda_{j,i}(\tau)\
                   \le\ B_i(\min_k\mu_k).                   \tag{5}
\]
Indeed, \(\boldsymbol\mu\) is the equilibrium in (3), as follows by
letting \(\tau\) tend to infinity in \(R_j\) and factoring \(1+z^j\).
For \(\phi_s(u)=su/(1+(s-1)u)\) put \(d_i=1+(s-1)\lambda_i\). Then
\[
 F_{j,i}(\phi_s(\boldsymbol\lambda))
   =F_{j,i}(\boldsymbol\lambda)
           +(j-1)(1-s)\frac{\lambda_i(1-\lambda_i)}{d_i}.    \tag{6}
\]
To check this, use
\(\phi_s(\lambda_i)-\phi_s(\lambda_k)
=s(\lambda_i-\lambda_k)/(d_i d_k)\) in (3), then use
\(d_k=d_i+(s-1)(\lambda_k-\lambda_i)\) in its sum.
Along \(\phi_s(\boldsymbol\mu)\), one has
\(\partial_s\phi_s(\mu_i)=\mu_i(1-\mu_i)/[1-(1-s)\mu_i]^2\).
Take \(s(0)=0\) and
\(\dot s=(j-1)(1-s)[1-(1-s)u]\).
The trajectory is a componentwise subsolution for \(u=\max_k\mu_k\)
and a supersolution for \(u=\min_k\mu_k\). To compare with the singular
initial data, delay the subsolution by any positive time and advance
the supersolution by any positive time. Ordered distinct roots then
permit cooperative comparison starting just after the delay for the
subsolution or just after zero for the supersolution. Send these shifts
to zero. Solving the scalar equation gives (5).

**Late-clock theorem.** Let \(N=n+1\), \(n\ge2\), and
\(0<\tau'\le \tau\). If
\[
 e^{-(N-1)\tau'} <
    \sin^2\!\frac{\pi}{2N}\,(1-e^{-2/n}),                      \tag{7}
\]
then \(R_n(\,\cdot\,;\tau)\) and \(R_N(\,\cdot\,;\tau')\)
strictly interlace in angle. A simpler sufficient condition is
\((N-1)\tau'>\log(nN^2)\).

To prove the outer inequalities, let
\(\alpha=(2i-1)\pi/(2N)\) and \(\beta=(1+1/n)\alpha\) for a common
positive index. Writing \(\mu=\cos^2\alpha\), we have
\[
 \log\frac{\mu}{\cos^2\beta}
       =2\int_\alpha^\beta\tan v\,dv
       \ge \frac2n(1-\mu).
\]
Here \(\alpha\tan\alpha\ge\sin^2\alpha\). Concavity of
\(1-e^{-2v/n}\) on \([0,1]\) gives
\((\mu-\cos^2\beta)/\mu\ge(1-\mu)(1-e^{-2/n})\).
By (5), the upper trajectory's deficit from \(\mu\), divided by \(\mu\),
is at most
\(e^{-(N-1)\tau'}(1-\mu)/[1-\max_k\mu_k]\).
Condition (7) makes the upper squared height larger than the lower
equilibrium \(\cos^2\beta\), which bounds the lower squared height
at every time. The inner gaps follow from the equal-clock comparison
and monotonicity. The simpler sufficient condition uses
\(\sin(\pi/(2N))\ge 1/N\) and \(1-e^{-2/n}\ge1/n\).

## A root-moment test for the outer gaps at the physical clock ratio

For the full interval \(\tau'/\tau\in[n/N,1]\), it suffices to prove the
outer gaps at the minimum ratio \(n/N\), since the upper heights increase
with their own clock. The [scaled Hermite comparison](MULTIVARIATE_WEAK_CURVATURE_ALL_TERMINALS.md)
gives strict outer gaps at this ratio for sufficiently small positive
time at each degree.

Suppose there is a first outer contact. At its common squared height
\(\lambda=x_i=y_i\), write the lower and upper ordered squared-height
arrays as \((x_k)\), \((y_k)\). When \(n\) is odd append a zero to the
lower array to pair with the extra positive upper height. The
intervals \([x_k,y_k]\) are nonnegative-oriented, and the strict inner
inequalities put each nondegenerate interval on one side of
\(\lambda\). Define the positive measure
\[
 d\nu(s)=\sum_{k\ne i} {\bf1}_{[x_k,y_k]}(s)\,ds .
\]
The paired forces satisfy the exact identity
\[
 F_{N,i}(\boldsymbol y)-F_{n,i}(\boldsymbol x)
 = (1-\lambda)
       \left[1+2\lambda\int\frac{d\nu(s)}{(\lambda-s)^2}\right].       \tag{8}
\]
For even \(n\), the change in the drift of (3) is \(1-\lambda\).
For odd \(n\) it is \(-(1-\lambda)\), and the force of the appended
zero contributes \(2(1-\lambda)\) to the upper drift after the added
smallest-root force is split into its zero value and its increment.
In both parities the remaining paired forces are the integrals in
(8), by integrating the derivative of \(1/(\lambda-s)\).

At \(\tau'=n\tau/N\), the sign of the contact velocity difference
is the sign of
\[
 n(1-\lambda)
     \left[1+2\lambda\int\frac{d\nu(s)}{(\lambda-s)^2}\right]
                                      -F_{n,i}(\boldsymbol x).       \tag{9}
\]
A first contact requires (9) to be nonpositive. The matched
moments of the contact measure are
\[
 \int s^r\,d\nu(s)=
      \frac{\sum_k y_k^{r+1}-\sum_k x_k^{r+1}}{r+1},
       \qquad r\ge0.                                                \tag{10}
\]
These differences are generated by (1) at the two prescribed times,
without locating the individual \([x_k,y_k]\). Any real polynomial
\(q\) gives
\[
 \int\frac{d\nu(s)}{(\lambda-s)^2}\ \ge\
    \frac{\big(\int q\,d\nu\big)^2}
         {\int(\lambda-s)^2q(s)^2\,d\nu(s)}                          \tag{11}
\]
when the denominator is positive. This is Cauchy–Schwarz. The
best lower bounds over polynomials of increasing degree converge to
the left side, since the compact support avoids \(\lambda\) and
polynomials uniformly approximate \((\lambda-s)^{-2}\) there.

Equations (8)–(11) give an all-order sufficient test for the missing
outer sign using differences of the linked squared-height power sums
and the lower root speed. A degree-uniform positive bound for (9)
throughout the intermediate clock region remains open. The equations
here apply to the Gaussian reference and leave the angular comparison
from that reference to ordinary factorial roots as a separate step.

## A birth–death comparison for the Gaussian squared heights

For each \(j=2m+\epsilon\), the polynomial
\[
 G_j(z;\tau)=\prod_{i=1}^m(1-\lambda_{j,i}+\lambda_{j,i}z)
           =\sum_{\ell=0}^m p_{j,\ell}(\tau)z^\ell
\]
is the generating function of a birth–death chain begun at state zero.
In the raw clock, its rates at state \(\ell\) are
\[
 a_{j,\ell}=(m-\ell)(m+\epsilon-\tfrac12-\ell),\qquad
 d_{j,\ell}=\ell(\ell-\tfrac12).
\]
In fact \(G_j=\sum_r b_{j,r}(z-1)^r\). Substitute (1) and collect
powers of \(z\) to get
\[
 \dot G_j(z)=(z-1)\sum_\ell p_{j,\ell}
          \big[a_{j,\ell}z^\ell-d_{j,\ell}z^{\ell-1}\big].
\]
The term with \(\ell=0\) in the death sum vanishes. This is the
forward equation on \(\{0,\ldots,m\}\), with \(G_j(z;0)=1\).

At the common normalized time \(t=j\tau\), the rates are divided by
\(j\). In consecutive degrees \(n,N\), at each common state the upper
birth rate is at least the lower rate and the upper death rate is at
most the lower rate. For even \(n=2m\) the normalized birth difference
is \((m-\ell)(n+\tfrac12-(m-\ell))/(nN)\). For odd \(n=2m+1\)
it is \((m+\tfrac12-\ell)(n-(m-\ell))/(nN)\). The normalized death
comparison uses the same raw death rate at a common state.
Couple reactions at equal states using the smaller rate in each
direction. The extra allowed reactions there are an upper birth
and a lower death. At distinct ordered states, nearest-neighbor
jumps cannot reverse the order. The count at degree \(N\) therefore
stochastically dominates that at degree \(n\) at the physical lower
clock ratio \(n/N\). In particular
\(b_{N,r}(t/N)\ge b_{n,r}(t/n)\) for every \(r\), by taking the
expectation of \(\binom{\ell}{r}\).

This comparison is weaker than the outer angle inequalities.
The [unequal-clock boundary example](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md#the-first-contact-sign-at-unequal-heat-parameters)
already has coordinatewise unequal squared heights apart from the
contact and gives strict count dominance, while the outer contact
velocity has the wrong sign. Applying the birth–death comparison to
(9) would require further information about the actual trajectories.

## Connection with the matched ordinary tanh moments

The sharp lower Gaussian clock ratio occurs as an actual diffuse
ordinary limit at every \(c>0\). For a positive integer \(\ell\)
take \(\lfloor c\ell\rfloor\) copies of the mirrored unit pair
\((\ell,\ell+1),(\ell+1,\ell)\). At fixed \(j\) every Gauss factor
has
\[
 T_j=\frac{j-1}{j(2\ell+1)+1},\qquad
 U_j=2\lfloor c\ell\rfloor\operatorname{artanh}(T_j).
\]
Thus the outer time \(U_j/(j-1)\) tends to \(c/j\). The signed
moments cancel and every unsigned moment of order \(p>1\) is
\(O(\ell^{1-p})\). The endpoint-normalized log multiplier for
one mirrored pair at \(k\) is
\(-k(j-k)/(j\ell)+O_j(\ell^{-2})\), by expanding the logarithms
of the factors \(j\ell+1+h\) and subtracting the endpoint chord.
Consequently the exact ordinary polynomials converge to
\(R_j(\cdot;c/j)\). The matched tanh error bounds can therefore
tend to zero even at the clocks in (9).

There is a direct conversion to squared-height power sums for any
mirrored ordinary profile. If \(S_{j,k}\) is the \(k\)-th Newton
sum of its normalized degree-\(j\) polynomial, then for \(r\ge1\)
\[
 \sum_i\lambda_{j,i}^r
 =\frac{j\binom{2r}{r}
         +2\sum_{k=1}^r\binom{2r}{r+k}S_{j,k}}{2\cdot4^r}.         \tag{12}
\]
Indeed the signed-height roots satisfy
\(\lambda=(2+z+z^{-1})/4\), conjugate pairs give two equal
contributions and the structural root \(-1\) contributes zero.
Expanding \((2+z+z^{-1})^r=z^{-r}(1+z)^{2r}\) and using
reciprocity proves (12). Only the first \(r\) endpoint
coefficients enter its Newton sums. Thus one can obtain the
moments in (10) from endpoint coefficients assembled by the
[matched tanh kernels](MULTIVARIATE_MATCHED_TANH_KERNEL_TRANSPORT.md),
without enclosing independent errors in the middle of the
coefficient array. A force comparison for an exact factorial
path would still be needed to use those moments in a factorial
analogue of (9).
