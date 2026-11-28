# Exact angular grid witnesses for the Gaussian agreement triangle relaxation

The radial-cut counterexample has an exact counterpart on
a polar grid with \(N>12\) equal angular bins, where \(N\)
is a multiple of three and the radial cells refine the
median-radius cut. The witness improves on the true
Gaussian simplex value at sufficiently high correlation
if \(N\ge93\). For \(15\le N\le90\) this direction instead
loses near correlation one. In particular, it loses on
grids of sizes \(24\) and \(36\), and wins at sizes \(93\),
\(96\) and \(120\). The witness uses the integrated Gaussian
objective and has all diagonal entries equal to one.
[Continuum radial-cut counterexample](GAUSSIAN_AGREEMENT_TRIANGLE_RELAXATION.md#a-balanced-mixed-radius-counterexample)
[The integrated grid relaxation](GAUSSIAN_AGREEMENT_TRIANGLE_DISCRETIZATION.md#integrated-finite-values-converge-to-the-continuum-value)

## Exact row balance

Set \(a=2\pi/3\), \(m=N/3\), \(\Delta=2\pi/N\), and
\(c_N=1-1/(2m)\). For \(x,y\) let \(t_N(x,y)\) denote
the circular geodesic distance between their angular
bin centers. Let \(B(x,y)\) be one when their radial
cells lie in opposite median halves and zero otherwise.
For sufficiently small \(\epsilon>0\), put
\[
 F_{N,\epsilon}(x,y)=
 \left(1-(1-c_N\epsilon)\frac{t_N(x,y)}a
                      -\epsilon B(x,y)\right)_+.             \tag{AG1}
\]
Here \(1-F_{N,\epsilon}\) is the truncation at one of the
sum of two nonnegative pseudometrics, namely the scaled
bin-center angular distance and the radial cut metric.
Truncation preserves the triangle inequality. The box,
symmetry and diagonal-one conditions also hold.

For a fixed \(x\), the angular bins at distance
\(j\Delta\) have multiplicity two for \(1\le j\le m\),
and multiplicity one at zero. The two radial halves
each have probability \(1/2\), independently of the
angular bin of \(Y\). If
\(0<\epsilon<\min\{1/2,[c_N(m+1)]^{-1}\}\), both radial
types in (AG1) have positive argument for \(j<m\).
At \(j=m\) only the same-half term is positive, and
for \(j>m\) both are zero. Thus the row mean is affine
in \(\epsilon\). At zero it is
\[
        \frac{1+2\sum_{j=1}^{m-1}(1-j/m)}{N}
        =\frac mN=\frac13 .
\]
Its derivative, multiplied by \(N\), is
\[
 -\frac12+2\sum_{j=1}^{m-1}\left(c_N\frac jm-\frac12\right)
         +c_N
             =m c_N-m+\frac12=0 .
\]
Moreover \(F_{N,0}\) is the average of the \(N\) grid
rotations of an ordinary tripod obtained by grouping
\(m\) consecutive bins per label. Two bins at distance
\(j\Delta\) have average agreement \((1-j/m)_+\).
Every one of these partitions is a regular Gaussian
simplex, so the integrated value at \(\epsilon=0\)
is the true \(S_\rho(T)\).

## Objective slope and the grid threshold

For correlated standard planar Gaussians \(X,Y\) at
correlation \(\rho\), put \(\alpha=\arccos\rho\),
\(\sigma^2=1-\rho^2\), and write \(t\) for their
unquantized geodesic angle. The slope of the integrated
objective of (AG1), which is affine throughout the
feasible interval above, is
\[
 L_{\rho,N}=
 \mathbb E\!\left[\left(c_N\frac{t_N}{a}-B\right)
                         1_{\{t_N<a\}}\right]
                  +c_N\mathbb P(t_N=a,\ B=0).              \tag{AG2}
\]
Its difference from \(c_N\mathbb E t_N/a-\mathbb E B\)
has magnitude at most \((3/2)\mathbb P(t_N\ge a)\).
Since \(|t_N-t|\le\Delta\) and \(a-\Delta>\pi/2\),
this probability is \(O(1-\rho)\) as \(\rho\uparrow1\).
Indeed \(t\ge\pi/2\) requires
\(\rho |X|+\sigma Z_1\le0\), where \(Z_1\) is standard
Gaussian independent of \(|X|\). Integration in the
Rayleigh radius gives probability at most
\(\sigma^2/(4\rho^2)\).

Common rotations of \(X,Y\) leave their law invariant.
For a fixed unquantized gap \(t\le\pi-\Delta\), averaging
over the common rotation makes the directed bin-center
gap take its two neighboring grid values with mean
\(t\). The values are both at most \(\pi\), so the
geodesic gap has the same mean. For \(t>\pi-\Delta\)
the error is bounded by \(\Delta\), and the same
Gaussian probability bound applies. Therefore
\(\mathbb E t_N=\mathbb E t+O(1-\rho)=\alpha+O(1-\rho)\).
The identity \(\mathbb E t=\alpha\) follows by
averaging sign disagreement over linear directions.
[The angle identity and large-angle bound](GAUSSIAN_AGREEMENT_TRIANGLE_RELAXATION.md#a-strict-gain-on-an-explicit-correlation-band)

For the exact radial crossing limit, let \(R=\sqrt{2\log2}\).
By symmetry \(\mathbb P(B=1)=2\mathbb P(|X|\le R<|Y|)\).
Condition on \(|X|=r\) and substitute \(r=R-\sigma u\)
in its Rayleigh integral. The conditional probability
is
\[
 \mathbb P\bigl(|\rho(R-\sigma u)e_1+\sigma Z|>R\bigr)
\]
with \(Z\) standard planar Gaussian. For each fixed
\(u\ge0\) it tends to \(\mathbb P(Z_1>u)\) as
\(\sigma\to0\), by expanding the squared radius.
On \(0\le u\le R/\sigma\) it is at most
\(\mathbb P(|Z|>u)=e^{-u^2/2}\), because
\(\rho(R-\sigma u)\le R-\sigma u\). The Rayleigh
density is bounded, so dominated convergence yields
\[
 \lim_{\rho\uparrow1}\frac{\mathbb P(B=1)}{\sigma}
     =2R e^{-R^2/2}\int_0^\infty\mathbb P(Z_1>u)\,du
     =\sqrt{\frac{\log2}{\pi}} .
\]
As \(\sigma/\alpha\to1\), (AG2) consequently gives
\[
 \lim_{\rho\uparrow1}\frac{L_{\rho,N}}{\arccos\rho}
    =\frac{c_N}{a}-\sqrt{\frac{\log2}{\pi}}.                \tag{AG3}
\]
The right side increases with \(N\). It is positive
precisely when
\[
            N>\frac{3}{2(1-\frac23\sqrt{\pi\log2})}.
\]
The threshold lies strictly between \(90\) and \(93\).
For instance \(3.1415<\pi<3.1416\) and
\(0.6931<\log2<0.6932\) put \(\pi\log2\) strictly
between \(9(59/60)^2/4\) and \(9(61/62)^2/4\).
These decimal bounds follow from the alternating
series in Machin's formula for \(\pi\) and the positive
series \(2\sum_{j\ge0}3^{-(2j+1)}/(2j+1)\) for \(\log2\).
This proves the stated signs for each grid size and
sufficiently high correlation. A negative slope
only rules out this witness direction on that grid.
