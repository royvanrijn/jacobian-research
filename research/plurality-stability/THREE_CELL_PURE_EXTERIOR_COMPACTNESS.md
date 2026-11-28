# Compactness of a pure exterior for Gaussian maximizers

Fix \(J=[a,b]\subset(0,1)\) and a compact set of positive
three-cell mass vectors. For each \(R>0\) there exist
\(\delta>0\) and \(R_0>R\) with the following property.
If \(A\) is an actual planar Gaussian noise-stability
maximizer in these parameter sets and for some label \(i\)
\[
                 \gamma_2(A_i\setminus B_R)\le\delta,
\]
then \(A_i\) is contained in \(B_{R_0}\) up to null sets.
In fact the adjusted-score winner set of \(i\), including its
ties with maximal scores, is contained in \(B_{R_0}\).
The constants depend on \(R\) and the parameter sets alone.
Consequently any unbounded winning label has a uniform
positive Gaussian mass outside each fixed disk.

## A forbidden annulus and its contraction

Let \(U_i=T_t1_{A_i}\) at the maximizing correlation \(t\).
By (PP1), on the closed adjusted-score winner set of \(i\)
one has \(U_i\ge\theta_*>0\), where the same \(\theta_*\)
works for every label and maximizer in the parameter sets.
By (PP2), there is \(M<\infty\) such that from any winning
point \(x\) of label \(i\) one can select a strict winning
point \(x'\) of the same label with \(|x'-tx|\le M\).
Indeed the indicated source ball has positive \(i\)-area
and the pair tie sets are null.
[Uniform score threshold and contracted source recurrence, PP1--PP2](THREE_CELL_PRESSURE_CONTRACTION_PROPAGATION.md)

Choose \(R_0>R\) so large that
\(\mathbb P(|Z|\ge aR_0-R)<\theta_*/4\) for a standard planar
Gaussian \(Z\). For \(|x|\ge R_0\) the contribution to \(U_i(x)\)
from \(A_i\cap B_R\) is then below \(\theta_*/4\): a sample
\(tx+\sqrt{1-t^2}Z\) in \(B_R\) must satisfy
\(|Z|\ge t|x|-R\ge aR_0-R\).

Fix an outer radius \(L\) with
\[
               L>R_0,\qquad bL+M<L,\qquad aL-M>R_0.
\]
For all \(t\in J\), \(|x|\le L\), and \(y\in\mathbb R^2\),
the relative Mehler kernel satisfies
\[
 M_t(x,y)=\frac{1}{1-t^2}
  \exp\!\left(\frac{|x|^2}{2}
      -\frac{t^2}{2(1-t^2)}|y-x/t|^2\right)
   \le C_L:=\frac{e^{L^2/2}}{1-b^2}.
\]
If \(\delta=\theta_*/(4C_L)\), the contribution from
\(A_i\setminus B_R\) to \(U_i(x)\) is at most \(\theta_*/4\)
on this disk. In particular \(U_i<\theta_*\) in the closed
annulus \(R_0\le |x|\le L\), which is free of winning
points of label \(i\).

If \(x_0\) were an \(i\)-winning point with \(|x_0|>L\),
iterate the contracted source recurrence to obtain \(i\)-winning
points with \(|x_{n+1}-t x_n|\le M\). The bound
\(|x_n|\le b^n|x_0|+M/(1-b)\) makes the sequence enter the
closed disk of radius \(L\), since \(M/(1-b)<L\). At the
first entry, the previous point has radius greater than \(L\),
so the new radius exceeds \(aL-M>R_0\). This contradicts
the forbidden annulus. The winner rule proves the claimed
containment.

## Uniform confinement near reference partitions

Suppose a reference partition \(A^*\) equals label \(k\)
outside a disk of radius \(R\). The preceding constants work
for every actual maximizer in the parameter sets whose
summed Gaussian label distance from \(A^*\) is less than
\(\delta\). For \(i\ne k\), its tail mass outside \(B_R\)
is bounded by the corresponding indicator difference from
\(A^*\), and it has no winning points outside \(B_{R_0}\).
No optimality assumption on \(A^*\) is needed for this
comparison.

In particular, if actual maximizers converge strongly in
Gaussian label distance to a partition with a single exterior
label, that same label fills the exterior of a common disk
for all sufficiently late maximizers. The parameter and mass
sequences need only lie in the stated compact positive sets.
This also gives a compact-family version. Cover any compact
family of reference partitions with single exterior labels
in Gaussian label distance by these neighborhoods, and take
a finite subcover. The largest of its finitely many output
radii is an exterior radius for every actual maximizer in
their union. Within each neighborhood the exterior label
is the one of its reference partition, which may differ
across the family.

## A first-moment lower bound on the exterior radius at a right tie

Let \(A\) be an ordinary balanced planar partition consisting of
one exterior label outside \(B_R\), with \(R<\infty\). Define
\[
 M_R=\int_0^R s^2e^{-s^2/2}\,ds,\qquad
 M=M_\infty=\sqrt{\pi/2},\qquad
 \mathcal D(R)=\frac{9}{8\pi}(M_R/M)^2.
\]
Then \(W_1(A)<\mathcal D(R)\), without a maximizing or stationary
hypothesis. In particular if \(A\) ties the centered simplex
\(T\) at \(0<r<1\) and simplex-optimal correlations accumulate
at \(r\) from the right, then
\[
       w(r):=2\bigl(S_r(T)-1/3\bigr)/r-S'_r(T)
                     \ \le\ W_1(A)\ <\ \mathcal D(R).           \tag{PE1}
\]
Whenever \(w(r)>0\), this forces a strictly larger exterior
radius than the unique solution of \(\mathcal D(R)=w(r)\).
The condition is a necessary radius bound rather than an
upper bound on the support of the bounded labels.

To prove the first-moment bound, put
\(m_i=\int x1_{A_i}\,d\gamma_2\) and \(W=\sum_i|m_i|^2=W_1(A)\).
Contributions to all three \(m_i\) from outside \(B_R\)
vanish: the two non-exterior labels are absent there and the
exterior integral is over a centered radial set. Hence
\[
 W=\int_{B_R}m_{A(x)}\cdot x\,d\gamma_2(x)
     \le \int_{B_R}\max_i m_i\cdot x\,d\gamma_2(x)
     =\frac{M_R}{2\pi}\sum_{i<j}|m_i-m_j|
     \le\frac{3M_R}{2\pi}\sqrt W .
\]
For the equality, the angular average of the maximum is
half the average of the maximum minus the minimum, by
inversion symmetry. For three real values the latter
difference is half the sum of the three pairwise absolute
differences. The circular average of
\(|(m_i-m_j)\cdot\omega|\) is \(2|m_i-m_j|/\pi\).
For the last inequality use Cauchy--Schwarz and
\(\sum_{i<j}|m_i-m_j|^2=3W\), since \(\sum_i m_i=0\).
Thus \(W\le\mathcal D(R)\).
If equality held, \(W>0\) and the three moment vectors
would form an equilateral triangle centered at zero.
Equality in the pointwise step would assign each point
of \(B_R\) to its largest linear score, up to null sets.
Each of the two bounded labels would then have exactly
\(\gamma_2(B_R)/3<1/3\) mass, a contradiction.

For (PE1), the right-contact hypothesis gives
\(S'_r(A)\le S'_r(T)\) by the derivative of the deficit for
the fixed partition \(A\). At a tie the positive-correlation
Hermite expansion gives
\[
 S'_r(A)=\frac{2}{r}\bigl(S_r(T)-1/3\bigr)-W_1(A)
                   +\sum_{j\ge3}(j-2)r^{j-1}W_j(A).
\]
The last sum is nonnegative, proving the inequality.
Also \(w(r)<W_1(T)=9/(8\pi)\): otherwise all simplex
Hermite weights of degree at least three would vanish
and its nonconstant indicators would be quadratic
polynomials almost everywhere, which bounded indicators
cannot be. The function \(\mathcal D\) is strictly increasing
from zero to \(W_1(T)\), giving the stated inverse.
[Hermite slope comparison at a right tie, NS4](THREE_CELL_NO_COMMON_CIRCLE_MOMENT_SLOPE.md#comparing-derivatives-at-a-tie)

If only one label is known to be contained in \(B_R\),
a related necessary bound is
\[
 W_1(A)<\left(\frac{3M_R}{2\pi}
                      +\frac{\sqrt2(M-M_R)}{\pi}\right)^2.
\]
The same angular comparison applies inside \(B_R\); on
each outer circle use only the two possible labels.
Their angular maximum has average
\(|m_i-m_j|/\pi\le \sqrt{2W}/\pi\).
Equality in that last bound requires \(m_i=-m_j\)
and the third moment zero, incompatible with the
nondegenerate equilateral equality inside \(B_R\).
