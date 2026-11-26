# An explicit cutoff for physical coprimality with two private units

Consider an ordinary two-monomial profile with one private unit
row on each side and arbitrary shared rows.

**Theorem.** At every lower order \(n\ge2048\), the consecutive
physical factorial moments are coprime over every characteristic-zero
field. The bound is an explicit number independent of the dimension,
the shared exponents and the signed degree difference. It does
not use the unevaluated absolute ordinary cutoff.

The [high curvature theorem](MULTIVARIATE_BOTH_PRIVATE_UNIT_CURVATURE_15_OVER_8.md)
covers shared absolute curvature \(S\ge15/8\).
The [weak curvature theorem](MULTIVARIATE_PRIVATE_UNIT_WEAK_THREE_TENTHS.md)
covers the shared smaller-endpoint curvature \(V\le3/10\).
The [higher-label theorem](MULTIVARIATE_PRIVATE_UNIT_HIGHER_LABEL_LINKED_RESOLVENT.md)
excludes labels \(k\ge2\) at the orders considered here.
It suffices to treat the first common positive angle box
with \(S<15/8\) and \(V>3/10\).

## Endpoint law at the first angular label

Use \(N=n+1\). In each order \(m=n,N\), cancel the private
factorials and normalize both endpoints to one. Let \(G_m\)
be the argument scale of this normalization, and let
\(T_m>0\) give the normalized coefficients as its moments through
degree \(m\). The independent Gamma representation gives
\(T_m=\prod_i X_i^{d_i}/G_m\) from the shared rows with
\(d_i=b_i-a_i\), where \(X_i\) has shape \(m a_i+1\).
There is also the reversed endpoint law. Write \(B_m=\log T_m\)
for the log of either endpoint variable, in its own orientation.

The [ordinary packing and absolute radial walls](MULTIVARIATE_PRIVATE_UNIT_IMBALANCE_TWENTY_FIVE.md)
give \(V<61/12<6\) and log radii \(|x_m|\le11/m\) at the
normalized roots. The associated triangular Jensen estimate
for the negative log endpoint mean gives
\[
 -\mathbb E B_m\ge \frac{V}{2(1+2\min(V,2)/3+1/m)}
                        \ge \frac9{40}\sqrt V.                \tag{1}
\]
For the first inequality use \(|d_i|/\ell_i\le\min(V,2)\)
in the denominator \(\ell_i+2|d_i|/3+1/m\) of the individual
row bounds. For the second use \(1/m<1/2000\). The ratio
\(\sqrt V/[2(1+2V/3+1/2000)]\) has a single maximum on
\([3/10,2]\), so checking both endpoints suffices there.
Above \(2\) the corresponding ratio increases.
The same endpoint law has
\[
  \operatorname{Var}(B_m)\le V/m,\qquad
  \mathbb E e^{t(B_m-\mathbb E B_m)}\le e^{t^2 V/(2m)}
                                      \quad(0\le t\le m).
\]

At a putative common root in the first box,
\(\pi/n<\theta<2\pi/(n+2)\). In particular
\(\theta<7/n\), \(\theta<\pi/2\) and \(\sin\theta\ge2/n\).
In all estimates below the real shift \(x\) can be any point
in \([-11/n,11/n]\). Put \(s=B_m+x\). Since \(n\ge2048\)
and \(V>3/10\), both \(|x|<\sqrt V/100\) and \(|x|<V/50\).
By (1), \(\mathbb E s\le-43\sqrt V/200\).
The Chernoff bound with \(t=(23/200)m/\sqrt V\) gives
\[
  \mathbb P(s>-\sqrt V/10)
              \le e^{-529m/80000}\le e^{-m/152}.
\]
The function \(m e^{-m/152}\) decreases for \(m\ge2048\),
and \((27/10)^{13}>6\cdot2048\cdot30\). Thus for either
endpoint law
\[
             nV\,\mathbb P(s>-\sqrt V/10)<1/30.             \tag{2}
\]

Set \(F(s)=(1-e^{s+i\theta})^{-1}\) and \(E_m(x)=\mathbb E F(s)\).
The bound \(|F|\le1/\sin\theta\le n/2\) holds everywhere.
With probability at least \(3/4\), Chebyshev gives
\(|B_m-\mathbb E B_m|\le V/12\), since \(mV>576\).
The weaker endpoint mean bounds
\(-V/2\le\mathbb E B_m\le-3V/16\) give
\(-5V/8\le s\le-V/12\) on this event.
Put \(a=1-e^s\). Then \(V/18\le a\le5V/8\), by
\(1-e^{-u}\ge u/(1+u)\) for \(u\ge0\) and \(V<6\).
Since \(\theta\le7V/(nV)<a/4\), the real part of the kernel
is at least \(a/(a^2+\theta^2)\ge(16/17)/a\).
All negative contributions to its real part occur at \(s>0\)
and cost at most \(1/(60V)\) by (2). Consequently
\[
                             |E_m(x)|\ge 1/V.                 \tag{3}
\]

For \(w<0\), the derivative obeys
\[
        |F'(w)|\le \frac{e^w}{(1-e^w)^2}
                      =\frac1{4\sinh^2(w/2)}\le w^{-2}.
\]
Its unrestricted bound is \(1/[4\sin^2(\theta/2)]\le n^2/4\).
Using the first estimate on \(s\le-\sqrt V/10\) and the second
on its complement gives, together with (2) and (3),
\[
  |\partial_x\log |E_m(x)||\le 100
                              +\frac{Vn^2}{4}\mathbb P(s>-\sqrt V/10)
                           < \frac{3n}{50}.                  \tag{4}
\]
The kernels are bounded at the fixed positive angle, so
differentiation under the expectation is valid. Conjugating
the angle leaves all modulus estimates unchanged.

## The adjacent packet difference

Use the
[exact signed Beta noise](MULTIVARIATE_PRIVATE_UNIT_HIGHER_LABEL_LINKED_RESOLVENT.md#exact-signed-log-noise)
to write \(B_n\) in law as \(B_N+Z\), with \(Z\)
independent of \(B_N\). In either orientation its log moment
generating function \(\Lambda\) satisfies
\[
 0\le \Lambda'(0)\le V/(4nN),\qquad
 0\le\Lambda''(t)\le V/[nN(1-2|t|/n)^2]\quad (|t|<n/2).
\]
In particular \(\mathbb E Z^2\le2V/n^2\).
At \(t=\pm n/4\), Taylor's integral bound gives
\(\Lambda(t)\le V/(16N)+V/8<1\), so \(\mathbb E e^{tZ}<3\).
Thus for \(p_Z=\mathbb P(|Z|>\sqrt V/20)\),
\[
                  p_Z\le6e^{-n\sqrt V/80},\qquad
                  nV p_Z<1/60.
\]
Indeed \(n\sqrt V>1120\), the function \(u e^{-u/80}\)
decreases there, and \((27/10)^{14}>6(5/2)1120\cdot60\).
This also gives \(\sqrt{p_Z}<1/100\).

Use \(s=B_N+x\). On \(s\le-\sqrt V/10\) and
\(|Z|\le\sqrt V/20\), the entire interval from \(s\) to \(s+Z\)
lies below \(-\sqrt V/20\). In this region
\[
   |F''(w)|\le \frac2{(1-e^w)^3}
                    \le\frac{23000}{V^{3/2}}.
\]
Here \(1-e^{-\sqrt V/20}\ge (8/9)\sqrt V/20\), since
\(\sqrt V/20<1/8\).
Expand \(F(s+Z)\) to first order on this event.
Independence bounds the absolute first-order term by
\[
 \frac{100}{V}\bigl(\mathbb E Z+\sqrt{\mathbb E Z^2\,p_Z}\bigr).
\]
The remainder costs at most
\(11500\,\mathbb E Z^2/V^{3/2}\).
On the complementary event use \(2\sup |F|\le n\)
and the two tail bounds. Dividing by (3) proves
\[
 \frac{|E_n(x)-E_N(x)|}{|E_N(x)|}
     \le \frac1{30}+\frac1{60}+\frac5n+\frac{57650}{n^2}
                                                    <\frac7{100}.     \tag{5}
\]
For the numerical bound use \(\sqrt V<5/2\),
\(V<6\), \(\sqrt{\mathbb E Z^2}\le2\sqrt V/n\),
and \(n\ge2048\).
This compares the packets at the same shift and in either
orientation.

## Integer drift

Let \(K_m(x)\) be the log modulus of the forward packet
at \(e^{x+i\theta}\) minus that of the reverse packet
at \(e^{-x-i\theta}\).
The finite geometric sum gives \(m x_m=K_m(x_m)\) at the
normalized roots. Equation (4) bounds \(|K_m'|\) by \(3n/25\)
between the two radii. Equation (5) bounds \(|K_n-K_N|\)
at a fixed radius by \(2[-\log(1-7/100)]\).
Subtracting the two root equations gives
\[
   \frac{22}{25}N|x_N-x_n|
      \le\frac{11}{n}+2[-\log(1-7/100)]<\frac15,
\]
using \(-\log(1-u)\le u+u^2/[2(1-u)]\).
At a common physical root \(x_N-x_n=\log(G_N/G_n)\).
The [exact center drift bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md#an-explicit-linked-radial-condition)
would then give \(N|x_N-x_n|>1-V/n\) if the
integer \(\Delta=\sum_i d_i\) were nonzero. The last
two bounds conflict for \(n\ge2048\).
When \(\Delta=0\), the
[balanced Euclidean theorem](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md)
excludes the physical common root. This completes the
explicit eventual proof in the remaining first-label regime.
