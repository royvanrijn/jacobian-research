# Quantitative bounds for three-cell maximizers near one half

On the correlation interval
\[
 J=[249/500,251/500]
\]
every balanced one-dimensional three-cell partition has stability at
least \(19/10000\) below the planar regular simplex. This gives an
explicit lower bound on the smaller singular value of the planar
cell-moment matrix of every global maximizer. There is also an explicit
sublevel estimate for its adjusted pair scores on each bounded disk.
These are necessary conditions for a competing maximizer. Global
simplex optimality on \(J\) remains open.

## An explicit one-dimensional gap

Write \(\phi,\Phi\) for the standard Gaussian density and distribution.
Put \(q=\Phi^{-1}(2/3)\), \(b=\Phi^{-1}(5/6)\), and \(I=\phi(q)\).
For a balanced one-dimensional partition let \(W_m\) denote its total
degree-\(m\) Gaussian Hermite weight. Balance and Parseval give
\(W_0=1/3\) and \(\sum_m W_m=1\). Consequently
\[
 S_{1/2}\le \frac5{12}+\frac38 W_1+\frac18 W_2 .              \tag{NH1}
\]

For any increasing statistic \(g\) on an atomless line, the largest
possible norm of its three cell integrals, at fixed equal masses, is
obtained by assigning three consecutive quantile intervals of \(g\).
Indeed, for a given vector of integrals \(m=(m_i)\), maximize its scalar
product with the new vector of integrals. Sort the coefficients \(m_i\).
Writing them as a constant and two nonnegative increments reduces the
maximization to integrating \(g\) over upper tails of masses \(1/3\) and
\(2/3\). These are simultaneously attained by the quantile assignment.
Cauchy–Schwarz then bounds \(|m|\) by the norm of that assignment's
integral vector. Apply this to \(g=x\) and to \(g=(x^2-1)/\sqrt2\),
sorting by \(|x|\) in the latter case. It gives
\[
 W_1\le M:=2I^2<265/1000,\qquad
 W_2\le u^2+(u-v)^2+v^2<182/1000,
 \quad u=\sqrt2\,qI,\quad v=\sqrt2\,b\phi(b).                  \tag{NH2}
\]
The second formula uses
\(\int_{|x|>t}(x^2-1)d\gamma_1=2t\phi(t)\).

Put \(\delta=M-W_1\). If \(\delta\ge1/50\), (NH1) and (NH2)
give \(S_{1/2}\le12751/24000\). For the remaining case, sort the three
first moments as \(s, t-s,-t\), where \(s,t\ge0\),
\(s\le 2t\), and \(t\le 2s\). Write \(d_1=I-s\) and \(d_3=I-t\).
These are nonnegative by the halfline bound. Their feasible polygon
has vertices \((0,0),(I/2,0),(0,I/2),(I,I)\). The convex quadratic
\(d_1^2+d_3^2-d_1d_3-I(d_1+d_3)/2\) is nonpositive at these
vertices, hence everywhere in the polygon. Expanding the squared
first moments gives
\[
 \delta\ge I(d_1+d_3),\qquad d_1+d_3\le\delta/(.36),         \tag{NH3}
\]
using \(I>36/100\).

For each extreme cell compare its second moment against its
corresponding outer halfline. Write \(k=qI/\sqrt2\) for the degree-two
moment of each halfline and \(\Delta_j\) for the cell's excess over
this moment. For the right cell \(A_1\), with \(H=(q,\infty)\),
\[
 d_1=\int |x-q|\,|1_{A_1}-1_H|d\gamma_1,\qquad
 \Delta_1=-\frac1{\sqrt2}\int (x+q)|x-q|
                         |1_{A_1}-1_H|d\gamma_1.              \tag{NH4}
\]
The first identity follows by subtracting the equal masses from the
first-moment difference. For the second, also subtract the equal
masses from the second-moment difference, and factor \(x^2-q^2\).
Reflection gives the same estimates for the left cell.

In (NH4) the positive part of the integrand, including the minus
sign, arises from \(x<-q\). Bound \(-x-q\) there by
\((3-q)+(-x-3)_+\), and bound \(x+q\) on the negative part
by \((3+q)+(x-3)_+\). With \(\bar\Phi(3)=1-\Phi(3)\), the
respective Gaussian tail integrals before division by \(\sqrt2\) are
\[
 \begin{split}
 \int_3^\infty (x-3)(x+q)\phi(x)dx&=(1-3q)\bar\Phi(3)+q\phi(3),\\
 \int_3^\infty (x-3)(x-q)\phi(x)dx&=(1+3q)\bar\Phi(3)-q\phi(3).
 \end{split}
\]
The first is at most \(q\phi(3)\), since \(q>1/3\). The second
is at most \(\phi(3)/3\), by Mills' bound
\(\bar\Phi(3)\le\phi(3)/3\). The scalar bounds below then give
\[
 \Delta_j\le (182/100)d_j+14/10000,\qquad
 |\Delta_j|\le (243/100)d_j+14/10000.                         \tag{NH5}
\]
The middle cell's degree-two moment is \(-2k-\Delta_1-\Delta_3\).
Put \(D=\delta/(.36)\). By (NH3)–(NH5) and \(k<111/1000\),
\[
 \begin{split}
 W_2&=6k^2+6k(\Delta_1+\Delta_3)
                  +2(\Delta_1^2+\Delta_3^2+\Delta_1\Delta_3)\\
 &\le 6(.111)^2+6(.111)(1.82D+.0028)+2(2.43D+.0028)^2.
 \end{split}
\]
Insert this bound and \(W_1<.265-\delta\) into (NH1). The resulting
upper bound is a convex quadratic in \(\delta\), so its maximum on
\([0,1/50]\) is at an endpoint. The values at the two endpoints are
less than \(.526\) and \(.5312\), respectively. In particular
\(12751/24000\) bounds every balanced one-dimensional three-cell
partition.

Here are elementary enclosures for the constants. All finite decimals
used as bounds are rational. For \(0\le x\le1\) define
\[
 H_n(x)=\sum_{j=0}^n\frac{(-1)^j x^{2j+1}}{(2j+1)2^j j!}.
\]
The alternating exponential series gives
\(H_5(x)\le \int_0^x e^{-y^2/2}dy\le H_6(x)\).
Using \(.39894<1/\sqrt{2\pi}<.39895\) gives the following
rational enclosures by substitution.

| \(x\) | Lower bound for \(\Phi(x)\) | Upper bound for \(\Phi(x)\) |
|---:|---:|---:|
| .430 | .66640 | .66641 |
| .431 | .66676 | .66678 |
| .967 | .83322 | .83324 |
| .968 | .83347 | .83349 |

Thus \(.430<q<.431\) and \(.967<b<.968\). Alternating
exponential sums through degrees five and six also give
\(.3635<I<.364\) and \(.2497<\phi(b)<.2500\). Using
\(1.4142<\sqrt2<1.4143\), these imply \(.221<u<.222\),
\(.341<v<.343\), and \(k<.111\). The bound for \(W_2\) in
(NH2) follows from \((.222)^2+(.343)^2+(.122)^2<.182\).
For (NH5), \(\phi(3)<.0045\) follows from \(e>2.71\),
\(e^{1/2}>1+1/2+1/8+1/48\), and \(1/\sqrt{2\pi}<.399\).
The bounds for \(1/\sqrt{2\pi}\) follow by squaring and using
\(3.14159<\pi<3.14160\), which in turn follows, for example,
from the alternating arctangent series in Machin's identity
\(\pi=16\arctan(1/5)-4\arctan(1/239)\) through order eleven.

The [exact simplex value](THREE_CELL_DEFICIT_AND_PERIMETER.md#exact-correlation-derivative-for-the-three-cell-simplex)
gives
\[
 S_{1/2}(T)=\frac23-\frac34
              \left(\frac12-\frac{\arcsin(1/4)}\pi\right)^2
             >2673/5000.
\]
For the last inequality use \(\pi<22/7\) and the first three
positive terms \(x+x^3/6+3x^5/40\) of \(\arcsin x\) at \(x=1/4\).
On \([.49,.51]\), the Hermite series gives
\(0\le S'_\rho\le(2/3)\sup_{m\ge1}m(.51)^{m-1}\le17/25\)
for every balanced partition. Termwise differentiation is justified
by the summable Hermite weights and the uniform derivative bound.
The maximum of \(m(.51)^{m-1}\) is at \(m=2\), since the ratios
for \(m\ge2\) are at most \(1.5(.51)\). Both derivatives thus lie
in the same interval of length \(17/25\), so their difference has
absolute value at most \(17/25\). In particular, for every
\(\rho\in J\) and every balanced one-dimensional partition \(D\),
\[
 S_\rho(T)-S_\rho(D)\ge 19/10000.                             \tag{NH6}
\]

## Directional perimeter and moment rank

Let \(A\) be a balanced planar global maximizer at \(0<\rho<1\).
Use the normal, score-gradient, and volume-flux conventions of the
[global interface bounds](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#a-global-matrix-bound-on-volume-flux).
In particular, \(P_v=\sum_e\int|v\cdot N_e|d\nu_e\) for a unit
vector \(v\), and let \(C v\) be the vector of projected cell
moments, with the common minus sign if using \(C=BL\) of (IB4).
Let \(G=L^*L\) in the \(L^2(a\,d\nu)\) norm, so
\(v^{\mathsf T}Gv=\sum_e\int a_e(v\cdot N_e)^2d\nu_e\).
The full-rank result for actual maximizers supplies invertibility of
\(C\) onto the centered label space and of \(G\). The input (IB4)
therefore applies.

For any centered label vector \(\lambda\), apply Cauchy–Schwarz
to the three interfaces and use (G3):
\[
 \lambda^{\mathsf T}BB^*\lambda
 \ge\frac{(\sum_{e=(i,j)}|\lambda_i-\lambda_j|\nu_e(\Sigma_e))^2}
          {\sum_e\int a_e\,d\nu_e}
 \ge\frac{3I^2}{2\rho S'_\rho(A)}|\lambda|^2.                \tag{NH7}
\]
For the second inequality, the sum in the numerator before squaring
is the integral of the Gaussian perimeters of the superlevel sets
of the function equal to \(\lambda_i\) on cell \(i\). This follows
by integrating the length \(|\lambda_i-\lambda_j|\) of the level
interval seen by each interface. For almost every level strictly between
the smallest and largest label value the volume is \(1/3\) or \(2/3\),
so Gaussian isoperimetry bounds the perimeter below by \(I\).
For three centered numbers the squared range is at least \(3/2\)
times their squared norm. Indeed, fixing a range \(r\) and letting
the middle value vary between the endpoints gives a convex quadratic
for the sum of squared centered values. Its maximum is \(2r^2/3\)
when the middle value is an endpoint.

By (IB4) and (NH7), inversion on the two-dimensional spaces gives
\(C^{\mathsf T}C\succeq[3(1-\rho)I^2/(2\rho S'_\rho(A))]G\).
A second Cauchy–Schwarz inequality and \(\rho\sum_e\int a_e^{-1}d\nu_e<1\)
from (IB1) yield
\[
 P_v\le \sqrt{\frac{2S'_\rho(A)}{3(1-\rho)I^2}}\,|Cv|,
 \qquad S'_\rho(A)\le \frac{4\phi(0)^2}{1-\rho^2}.           \tag{NH8}
\]
The last bound follows from (IB2) and (IB5). In particular the
factor in the directional inequality is explicit on compact
positive-correlation intervals.

For completeness, a small directional perimeter gives an actual
one-dimensional competitor in the following quantitative sense. Let
\(p_i(x)\) be the Gaussian masses of the cells on the fibers parallel
to \(v\), and assign each fiber to a label of maximum mass, using a
fixed tie rule. If \(D_0\) is this cylindrical partition, its total
symmetric-difference mass with \(A\) is
\(2\int(1-\max_i p_i(x))d\gamma_1(x)\). This is at most
\(\sum_i\int2p_i(x)(1-p_i(x))d\gamma_1(x)\), since
\(\sum_i p_i(x)^2\le\max_i p_i(x)\). On each fiber the Gaussian
isoperimetric profile \(\mathcal I(p)=\phi(\Phi^{-1}(p))\)
satisfies \(\mathcal I(p)\ge2\phi(0)\min(p,1-p)\) by concavity.
Thus \(2p(1-p)\le\mathcal I(p)/\phi(0)\), which is bounded
by the fiber perimeter divided by \(\phi(0)\). The one-dimensional
BV slicing formula, with Gaussian weights in both coordinates,
integrates the sum of the three fiber perimeters to \(2P_v\).
Consequently the symmetric-difference mass is at most \(2P_v/\phi(0)\).
By reallocating surplus mass on the atomless base line, make a
balanced cylindrical partition \(D\). This changes \(D_0\) by a
total symmetric-difference mass at most the sum of the three
original cell mass errors, which is itself at most the preceding
symmetric-difference mass with \(A\). We obtain
\[
 \sum_i\gamma_2(A_i\triangle D_i)\le 4P_v/\phi(0).
\]
Changing a partition from \(A\) to \(D\) can change agreement
only if one of the two Gaussian points changes label. Hence the
absolute stability difference is at most the total symmetric-difference
mass. For \(\rho\in J\), maximality and (NH6) give
\(P_v\ge 19\phi(0)/(4\cdot10000)\). Combining with (NH8) proves
\[
 \min_{|v|=1}|Cv|\ge\frac{19I}{8\cdot10000}
             \sqrt{\frac{3(1-\rho)(1-\rho^2)}2}
               >\frac5{10^5},\qquad \rho\in J.              \tag{NH9}
\]

## Uniform sublevel estimates for the maximizing scores

Write \(F=T_\rho f-c\) for any adjusted pair contrast of a balanced
planar global maximizer, where \(f=1_{A_i}-1_{A_j}\) and \(|c|\le1\).
The bound on the price difference follows because each label wins on
a set of positive mass and the scores lie in \([0,1]\). A product
comparison gives an explicit nonzero lower bound on the smoothed
contrast. Keep \(A_k\) fixed and use the sign of an independent Gaussian
coordinate to split \(A_i\cup A_j\) into its two equal-mass labels.
Dimension independence gives the same optimum in the product space.
With \(g=1_{A_i\cup A_j}\), the contributions of the two labels are
\((\langle g,T_\rho g\rangle+\langle f,T_\rho f\rangle)/2\) for
the original partition and
\((1/2+\arcsin(\rho)/\pi)\langle g,T_\rho g\rangle\) for the
independent split. Therefore
\[
 \langle f,T_\rho f\rangle\ge
             \frac{8}{9\pi}\arcsin \rho.                    \tag{NH10}
\]
Here \(\langle g,T_\rho g\rangle\ge(2/3)^2\) by its constant
projection. This is the independent split comparison of the
[product analysis](THREE_CELL_PRODUCT_RECOMBINATION.md#transition-matrix-restrictions-from-explicit-product-rules).
For \(\rho\in J\) the right side exceeds \(7/50\), using
\(\arcsin \rho\ge\rho\) and \(\pi<22/7\). As \(f\) has
mean zero and \(\|f\|_2=\sqrt{2/3}\), this implies
\(\|F\|_2>1/6\). Also \(|F|\le2\) on the real plane. Its squared
norm outside \(B_4\) is at most \(4e^{-8}<1/100\), so some
\(x_0\in B_4\) has \(|F(x_0)|\ge1/10\).

Fix \(R\ge1\), and set
\[
 L=R+5,\qquad M_R=1+e^{3L^2},\qquad
 N_R=\max\left(1,\left\lceil\frac{\log(10M_R)}{\log(4/3)}
                          \right\rceil\right).
\]
For any real unit vector \(e\), the function
\(F_e(z)=F(x_0+ze)\) extends to an entire function of the complex
variable \(z\). In the Gaussian conditional density write the
quadratic form without conjugation. Its modulus gives
\(|F_e(z)|\le1+\exp(\rho^2|\operatorname{Im}z|^2/(2(1-\rho^2)))\),
so \(|F_e|\le M_R\) on \(|z|\le4L\) for \(\rho\in J\).
Jensen's formula gives at most \(N_R\) zeros, with multiplicity,
in \(|z|<3L\). Remove these zeros by the finite Blaschke product
on the disk of radius \(3L\). The remaining analytic function
has modulus at most \(M_R\) by the maximum-modulus principle,
has modulus at least \(1/10\) at zero, and has no zeros in the
open disk. Harnack's inequality for the nonnegative harmonic
logarithm of \(M_R\) divided by its modulus gives modulus at
least \(1/(100M_R)\) on \(|z|\le L\).

On \(|z|\le L\), each Blaschke factor for a zero \(a\) has
modulus at least \(|z-a|/(4L)\). Thus, for
\(0<u<1/(100M_R)\), a real point with \(|z|\le L\) and
\(|F_e(z)|\le u\) is within \(4L(100M_Ru)^{1/N_R}\) of at
least one zero. The total length of such real points is at
most \(8LN_R(100M_Ru)^{1/N_R}\). Integrating in polar
coordinates about \(x_0\), with radial Jacobian at most \(L\),
and bounding the Gaussian density by \(1/(2\pi)\), proves
\[
 \gamma_2\{x\in B_R:|F(x)|\le u\}
 \le 8L^2N_R(100M_Ru)^{1/N_R}.                               \tag{NH11}
\]
The bound is uniform over all three pairs and all planar global
maximizers on \(J\). It requires neither gradient transversality
there nor a bound on the number of interfaces.

One application is to finite score approximations. If all three
adjusted scores are uniformly approximated on \(B_R\) to error
\(\varepsilon<1/(200M_R)\), then changing winners inside \(B_R\)
is confined to the union of the three sublevel sets in (NH11)
with \(u=2\varepsilon\). The additional disagreement mass outside
the disk is at most \(e^{-R^2/2}\). Twice the total disagreement
mass is the sum of the cell symmetric-difference masses. By
conditional expectation it bounds the integrated sum of absolute
circular mass errors.
There are explicit finite approximations of the scores. If
\(\eta\in(251/500,1)\), truncating the Hermite expansion of
\(T_\rho1_{A_i}\) at degree \(n\) gives uniform error on the disk
at most
\[
 \frac1{\sqrt3}\left(\frac{251}{500\eta}\right)^{n+1}
                  \sup_{x\in B_R}\sqrt{M_{\eta^2}(x,x)} .
\]
This follows by Cauchy–Schwarz from the Mehler expansion and
\(\|1_{A_i}\|_2^2=1/3\). The finite coefficients and centered
prices range over bounded sets. These estimates give an explicit
error for approximating maximizer profiles by polynomial winners
on a bounded disk. Ruling out all resulting profiles outside the
known local comparison neighborhood still requires a comparison.

The evidence consists of the analytic derivations and an internal
audit of the rational bounds. The derivations use the existing
dimension reduction, stationarity, and interface bounds for actual
global maximizers.
