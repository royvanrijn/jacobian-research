# Absolute imbalance twenty five for private unit pairs

Consider an ordinary two-monomial profile with exactly one private
unit row on each side and arbitrary shared rows.
Write \(\Delta=\sum_i(b_i-a_i)\) for the signed degree
difference. If \(|\Delta|\ge25\), the consecutive physical
factorial moments are coprime in every degree. The bound
uses no eventual terminal cutoff.

More precisely, when the shared curvature is below \(15/8\),
a physical common root forces \(|\Delta|\le24\). The proof
bounds the log radius of each normalized root by \(11/m\)
for \(m\ge7\), and \(10/m\) for \(4\le m\le7\).
It refines the
[uniform endpoint packet argument](MULTIVARIATE_PRIVATE_UNIT_ABSOLUTE_SIGNED_IMBALANCE.md)
by using the imaginary resolvent below the concentration
threshold and a sharper real-kernel bound above it.

## Shared curvature and the endpoint laws

Use the notation
\[
       V=\sum_{\rm shared}d_i^2/\ell_i,\qquad
       S=\sum_{\rm shared}d_i^2/(\ell_i+d_i),
\]
where \(d_i>0\) is the absolute gap and \(\ell_i>0\)
the smaller exponent of each unequal shared row.
The [curvature fifteen-eighths theorem](MULTIVARIATE_BOTH_PRIVATE_UNIT_CURVATURE_15_OVER_8.md)
covers \(S\ge15/8\), so take \(0<S<15/8\). The case \(S=0\)
has \(\Delta=0\). The ordinary packing gives
\[
                               V<61/12<6.
                                                               \tag{1}
\]
Indeed a row with \(d_i/\ell_i>1\) and \(d_i\ge4\) costs
more than two in \(S\). The only remaining possibilities
with this ratio and \(S<15/8\) are \((\ell_i,d_i)=(1,2),(2,3)\).
At most one can occur because their costs are at least \(4/3\).
Every other row contributes at most twice its \(S\) to \(V\).
The three resulting bounds are \(2(15/8)\),
\(4+2(15/8-4/3)=61/12\), and
\(9/2+2(15/8-9/5)\).

After the private unit factorials cancel, normalize both
endpoint coefficients to one by a positive argument scaling
and a positive scalar. Let \(\mu\) be its positive endpoint
moment law in degree \(m\), and \(\mu^\vee\) its order-\(m\)
tilted reciprocal law. For either law let \(B=\log T\)
and \(\mu_B=\mathbb E B\). The log Mellin function \(f\)
has zero endpoints and
\[
 -V/2\le\mu_B\le-3V/16,\qquad
 \operatorname{Var} B\le V/m,\qquad 0\le f''(x)\le V/m
                                                 \quad(0\le x\le m).
                                                               \tag{2}
\]
The upper trigamma integral bound gives the second derivative
estimate. Integration between the endpoints gives the slope
and variance upper bounds. For the negative mean, the triangular
Jensen estimate in the [endpoint means proof](MULTIVARIATE_PRIVATE_UNIT_ABSOLUTE_SIGNED_IMBALANCE.md#endpoint-means-and-concentration)
gives \(-f'(0)\ge6V/31>3V/16\) for \(m\ge4\).
Here \(d_i/\ell_i\le2\) by (1) and its packing argument.
On \([-1,1]\), the same trigamma upper bound gives
\(f''\le2V/m\), since \(m\ell_i-d_i\ge m\ell_i/2\).
In particular
\[
 \mathbb E e^{s(B-\mu_B)}\le e^{s^2 V/(2m)}\quad(0\le s\le m),
 \qquad
 \mathbb E\cosh(B+t)\le e^{V/m}\cosh(|t|+V/2).
                                                               \tag{3}
\]
These estimates apply to either law separately.

## Radial walls by imaginary-part bounds

At \(z=e^{t+i\theta}\), \(\pi/m\le\theta\le\pi\), put
\(X=B+t\) and \(E=\mathbb E(1-e^{X+i\theta})^{-1}\).
For \(0<\theta\le\pi/2\), the pointwise upper bound
\(1/\sin\theta\) and Jensen's inequality in the imaginary
part bound the ratio of any two such packet moduli by
\[
                   \frac{2(\mathbb E\cosh X-\cos\theta)}
                        {\sin^2\theta},
                                                               \tag{4}
\]
where the expectation is for the packet in the denominator.
If \(m^2(\mathbb E\cosh X-1)\le K\), then (4) is at most
\(K/2+2\), since \(\sin\theta\ge2\theta/\pi\ge2/m\).
For \(\theta\ge\pi/2\) the kernel's modulus is at most one.
Its real part is at least \((1+e^X)^{-2}\).
Jensen and \(\mathbb E T\le1\) bound the packet modulus ratio
above by \((1+e^{|t|})^2\).

Take \(mV\le512\) and \(|t|=11/m\). Uniformly for \(m\ge4\),
\[
                 m^2(\mathbb E\cosh X-1)
                       \le K=86^2(86/80)\,12.                \tag{5}
\]
For \(10\le m\le86\), use (3) and bound the left side by
\(m^2e^{6/m}\cosh(3+11/m)\). This is increasing in \(m\),
since its logarithmic derivative is at least
\((2m-17)/m^2\).
At \(86\), use \(e^{6/86}\le86/80\) and
\(\cosh(3+11/86)<12\).
For \(4\le m\le9\), the upper bound
\(81\cdot5\cdot[(11/4)^6+1]/2\) is less than \(K\).
For \(m\ge86\), use \(mV\le512\) in (3) to obtain
\(m^2[e^{512/m^2}\cosh(267/m)-1]\).
This expression decreases with \(m\) by the positive power
series in \(1/m^2\). Its value at \(86\) is bounded
by the same \(K\).
The comparisons \(\cosh(3+11/86)<12\) and \(e^{6/4}<5\)
follow from \(27/10<e<11/4\) and
\(e^{11/86}\le86/75\).
Thus (4) is at most \(K/2+2<(27/10)^{11}<e^{11}\).
The bound at \(\theta\ge\pi/2\) is smaller as well.

For \(4\le m\le7\) instead use \(|t|=10/m\). Equation (3)
gives \(m^2(\mathbb E\cosh X-1)\le49\cdot5\cdot140\),
using \(\cosh(11/2)<140\).
Indeed \(e^{11/2}<(11/4)^5(7/4)\) and \(e^{-11/2}<1\).
Then (4) is less than \(e^{10}\) and the large-angle bound
is also less than \(e^{10}\). This needs no restriction
on \(mV\).

## The concentration regime

Take \(mV\ge512\), \(m\ge4\), and \(|t|=11/m\le3V/128\).
Equation (2) gives \(\mathbb E X\le-21V/128\).
Chebyshev's inequality gives probability at least \(1/2\)
to the event \( |B-\mu_B|\le V/16\).
On this event
\[
                         -75V/128\le X\le-V/16 .
\]
The Chernoff estimates from (3) are
\[
 \mathbb P(X>0)\le e^{-55mV/4096},\qquad
 \mathbb P(X>-V/32)\le e^{-9mV/1024},
\]
obtained with \(s=5m/32\) and \(s=m/8\), respectively.
At \(mV=512\) these bounds are less than \(1/900\) and \(1/85\).
The functions \(u e^{-55u/4096}\) and \(u e^{-9u/1024}\)
decrease for \(u\ge512\).

For \(\theta\le V\) and \(\theta\le\pi/2\), on \(X\le-V/32\)
the kernel is bounded above by
\(\min(64/V,2/\theta)\le66/(V+\theta)\).
On the complement its modulus is at most \(m/2\), using
\(\theta\ge\pi/m\). The second tail estimate therefore adds
at most \(256/(85V)\) to the upper bound, giving
\(73/(V+\theta)\).

To bound the real part below, let \(\sigma=1-e^X\) on the
Chebyshev event. There \(V/32\le \sigma\le75V/128\), and
the real part of the kernel is at least
\(\sigma/(\sigma^2+\theta^2)\). All negative contributions
to the real part are on \(X>0\).
For \(\theta\le V/32\), the main event yields
\(32/(75V)\) and the negative tail costs at most
\(256/(900V)\). For \(V/32\le\theta\le V\), one has
\(\theta\sigma/(\sigma^2+\theta^2)\ge1/33\).
To see this, take the ratio \(\theta/\sigma\) in
\([4/75,32]\) and minimize \(x/(1+x^2)\) at the endpoints.
The main contribution is thus at least \(1/(66\theta)\),
and the negative tail costs at most \(1/(450\theta)\)
by the upper \(2/\theta\) on the kernel.
Both cases give \(|E|\ge1/[512(V+\theta)]\).

When \(V\le\theta\le\pi/2\), the imaginary part everywhere
has positive sign. On the main event \(e^X>1/64\) and the
kernel denominator is at most \(V^2+\theta^2\le2\theta^2\).
Together with \(\sin\theta\ge\theta/2\) and event
probability at least \(1/2\), this gives
\(\operatorname{Im}E\ge1/(512\theta)\).
The upper bound \(2/\theta\) suffices here.
For \(\theta\ge\pi/2\) the upper bound is one and the
main event gives real part at least \(1/8\).
Thus in every angle under consideration the packet ratio
is at most \(73\cdot512<e^{11}\).
All bounds used on averages apply to mixtures of the two
endpoint laws as well.

## Physical separation

The exact geometric sum at \(z=e^{t+i\theta}\) equals the
left endpoint packet plus \(z^m\) times the right packet
evaluated at \(z^{-1}\).
The bounds above hold for both signs of \(t\), so the two
terms cannot cancel when \(|z^m|=e^{11}\) or \(e^{-11}\),
and at orders \(4\) through \(7\) they cannot cancel on
the tighter walls \(e^{\pm10}\).
Interpolation of \(\mu\) to \((\mu+\mu^\vee)/2\) preserves
the packet estimates and the endpoint normalization.
The reciprocal mixture has unit-circle roots by the
[geometric angle theorem](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md).
Positive coefficients exclude the angles below \(\pi/m\)
throughout the interpolation, since the terms there lie in
a sector of width less than \(\pi\).
Root continuity therefore gives the claimed normalized radial
bounds. Conjugation handles negative angles.

At a common physical root of orders \(n,N=n+1\), the
[endpoint drift](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md#an-explicit-linked-radial-condition)
bounds the difference from \(\Delta\log(N/n)\) by \(V/(nN)\).
For \(n\ge7\) the radial bounds and \(\log(N/n)>1/N\)
give \(|\Delta|<11(15/7)+6/7<25\).
For \(4\le n\le6\), the tighter walls give
\(|\Delta|<10(2+1/4)+6/4=24\).
At \(n=2,3\), Eneström–Kakeya using \(|f'|\le V/2\)
bounds each normalized log radius by \(V/2\).
This gives \(|\Delta|<(61/12)(4+1/3)<25\).
At \(n=1\), the second moment is strictly positive
on the real line and the first root is real.
This proves the assertion. The complex coprimality descends
to every characteristic-zero field by rational Bézout identities.
