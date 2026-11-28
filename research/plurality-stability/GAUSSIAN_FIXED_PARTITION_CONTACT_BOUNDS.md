# Finite contact bounds for fixed Gaussian competitors

For every fixed cell count, a first-moment deficit controls the
number and order of possible contacts of a fixed balanced
partition with the simplex. The resulting estimates are
independent of ambient dimension and hold for arbitrary
measurable cells. Under the first-moment theorem used for the
small-correlation interval, the estimate is uniform over
partitions a fixed positive Gaussian distance from the simplex
orbit. The perimeter limit further excludes either endpoint
for each fixed nonsimplex shape.

## Holomorphic stability difference

Let \(A\) be balanced with \(k\) cells in dimension \(d\ge k-1\)
and let \(T\) denote the reference simplex. Write
\[
 W_m(A)=\sum_{i=1}^k\|\Pi_m1_{A_i}\|_2^2,\qquad
 H_A(z)=\sum_{m\ge1} \bigl(W_m(T)-W_m(A)\bigr)z^{m-1}.
\]
Here \(\Pi_m\) is orthogonal projection onto Hermite degree
\(m\). The constant noise-stability terms agree by balance.
For \(0<\rho<1\),
\[
              \rho H_A(\rho)=S_\rho(T)-S_\rho(A).
\]
Parseval gives \(\sum_{m\ge1}W_m(A)=1-1/k\). Thus \(H_A\)
is holomorphic on the unit disk, continuous on its closure,
and bounded in modulus by \(M=2(1-1/k)\) there. Its value
at zero is the first-moment deficit
\[
                   \Delta_A=W_1(T)-W_1(A).
\]
All of the analytic estimates below use just
\(\Delta_A\ge\lambda>0\). If the Gaussian label distance
\(d_k(A)\) to the simplex orbit is at least \(\eta>0\),
the dimension-independent small-correlation estimate
(L3), divided by \(\rho\) and taken to zero, gives
\(\Delta_A\ge C_k\eta^2\). This application uses the same
global first-moment input as (L3).
[Small-correlation stability and its first-moment input, L3](SMALL_POSITIVE_CORRELATION_INTERVAL.md#the-remaining-global-implication)

## A uniform bound on contacts and near contacts

Fix \(0<b<r<1\). Let \(Z_A\) be the zeros of \(H_A\) in
\(\{|z|\le r\}\), with multiplicities. Then
\[
  |Z_A|\le N_r=\left\lfloor \frac{\log(M/\lambda)}{\log(1/r)}\right\rfloor.
                                                               \tag{FCB1}
\]
Indeed for each finite multiset of zeros \(a\) divide
\(H_A\) by the product of the factors
\[
                         \frac{z-a}{1-\overline a z}.
\]
The quotient has removable singularities at the chosen
zeros and modulus at most \(M\) by the maximum principle,
since the factors have modulus one on the unit circle.
At zero the value has modulus at least
\(\lambda r^{-N}\) if the multiset has \(N\) elements
of modulus at most \(r\). This proves (FCB1).

Divide by all the factors for \(Z_A\), and call the
quotient \(G\). Then \(\log(M/|G|)\) is nonnegative
harmonic in \(|z|<r\), with value at zero at most
\(\log(M/\lambda)\). Harnack's inequality gives on
\(|z|\le b\)
\[
                |G(z)|\ge C,\qquad
                C=M(\lambda/M)^{(r+b)/(r-b)}.
\]
For real \(0\le\rho\le b\) the factored terms each have
modulus at least \(|\rho-a|/(1+br)\). Hence
\[
 |H_A(\rho)|\ \ge\
 C\left(\frac{\min\{1,\operatorname{dist}(\rho,Z_A)\}}
                {1+br}\right)^{N_r},                       \tag{FCB2}
\]
with the factor taken to be one when \(Z_A\) is empty.
When \(N_r>0\) and the threshold \(\epsilon\) is sufficiently
small, the set \(\{|H_A|\le\epsilon\}\cap[0,b]\) lies within
distance \((1+br)(\epsilon/C)^{1/N_r}\) of those complex
zeros. The length of the real sublevel set is at most
\(2N_r(1+br)(\epsilon/C)^{1/N_r}\), with the bound
interpreted as an empty sublevel set for \(N_r=0\) and
\(\epsilon<C\). The same estimate for the unnormalized
stability difference holds on \([a,b]\) for \(a>0\) after
replacing \(\epsilon\) by \(\epsilon/a\).

There is also a finite Taylor test with an explicit
margin. At \(0\le\rho\le b\), let \(a_j(\rho)\) be the Taylor
coefficients at zero of
\[
                         H_A\!\left(\frac{\rho+w}{1+\rho w}\right).
\]
These are linear combinations of the derivatives of
\(H_A\) at \(\rho\) up to order \(j\). Cauchy's inequality
gives \(|a_j(\rho)|\le M\). Evaluating at \(w=-\rho\)
and choosing \(n\) with \(M b^n/(1-b)\le\lambda/2\)
gives
\[
                 \max_{0\le j<n}|a_j(\rho)|
                              \ge \lambda(1-b)/2 .         \tag{FCB3}
\]
Thus a fixed distant partition cannot match a sufficiently
long finite correlation jet of the simplex. This result also
applies at a correlation where both partitions are global
maximizers. It does not assert that stationarity or global
maximality forces those jets to match.

## Each fixed shape has finitely many positive contacts

For a measurable partition \(A\), write
\(P_\gamma(A)=\tfrac12\sum_i P_\gamma(A_i)\) for the
Gaussian cluster perimeter. Put \(t=-\log\rho\).
The Ornstein–Uhlenbeck characterization of Gaussian
perimeter gives
\[
 \lim_{\rho\uparrow1}\frac{1-S_\rho(A)}{\sqrt{-\log\rho}}
       =\frac{2}{\sqrt\pi}P_\gamma(A),                     \tag{FCB4}
\]
including value \(+\infty\). In fact the numerator is
\(\sum_i\int_{A_i^c}T_\rho1_{A_i}\,d\gamma\).
Each term divided by \(\sqrt t\) has limit
\(P_\gamma(A_i)/\sqrt\pi\) if this perimeter is finite,
and infinite lower limit otherwise by the finite
lower-limit criterion.
[Ambrosio, Miranda and Pallara, Section 3, equation (3.1) and its converse](https://cvgmt.sns.it/media/doc/paper/2601/finepropertiesfin.pdf)

The Gaussian multi-bubble theorem gives
\(P_\gamma(A)>P_\gamma(T)\) for a balanced nonsimplex
partition, by equality rigidity. A translated standard
simplicial cluster with all volumes equal is centered.
One way to see this is to write the cell probabilities
as winners of independent Gaussian scores with additive
shifts. Unequal shifts give unequal winning probabilities
by interchange of the two coordinates.
[Milman and Neeman, Theorems 1.1 and 1.2](https://annals.math.princeton.edu/2022/195-1/p02)
[Strict comparison of translated simplex masses](GAUSSIAN_POSTERIOR_PARTITION_COMPATIBILITY.md#the-conditional-comparison)

Consequently \(S_\rho(A)<S_\rho(T)\) for all \(\rho\)
sufficiently close to one for each **fixed** nonsimplex \(A\).
If \(\Delta_A>0\), strict comparison also holds sufficiently
close to zero, and (FCB1) bounds the number of contacts in
between. With the existing first-moment theorem this covers
every fixed balanced nonsimplex. In particular the set where
one such partition beats the simplex is a finite union of
open intervals with closures in \((0,1)\). The exclusion
near one depends on the partition. These bounds do not
control an optimizing partition that changes with the
correlation.
