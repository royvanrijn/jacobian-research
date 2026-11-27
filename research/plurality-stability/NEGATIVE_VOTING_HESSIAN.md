# The simplex interface gap across the negative voting interval

The first-moment interface form at the regular \(k\)-simplex has an
optimal gap of order \(\ell\log(k)/k\) modulo rotations. Here
\(\ell\) is the distance between two cell first moments. The full
negative paired Gaussian Hessian has a gap of order at least
\(\log(k)/k\) in the norm weighted by the exact local score derivative,
uniformly for \(0<r\le1/(k-1)\) and every ambient dimension \(d\ge k-1\)
when \(k\) is sufficiently large. Each partition in a pair has the same
mass vector as its counterpart, with the common vector free.

The [weighted interface formulation](NEGATIVE_SIMPLEX_INTERFACE_LIFT.md#weighted-hessian-and-compactness)
is the specialization to the regular simplex of the paired variational
problem in Heilman's *Hyperstable Sets with Voting and Algorithmic
Hardness Applications*, Problem 7.1 and Lemma 7.6. Those formulas
describe the form whose sign is estimated here.
[Heilman, arXiv:2209.11216](https://arxiv.org/html/2209.11216#S7).

## The exact two-row variance

Use the active label space \(E=\mathbf1^\perp\subset\mathbb R^k\).
Let \(\nu_{12}\) be Gaussian surface measure on the facet \(F_{12}\)
where the first two coordinates are equal and are maximal, and put
\(v_k=\nu_{12}(F_{12})\) and \(\ell=kv_k\). For \(k\ge4\) write
\(\gamma_k\) for the variance of a unit linear functional in
\(H^0=\{x\in E:x_1=x_2=0\}\) under the normalized facet law. In the
[balanced first-moment block decomposition (S82)](SIMPLEX_INTERFACE_SYMMETRY.md#an-explicit-first-moment-gap-from-convex-covariance),
the two-row interaction has its only nonzero multiplicity eigenvalue
\((k-1)v_k\gamma_k\). The standard interaction eigenvalue is at most
\(\ell/2\). The rotation interaction has eigenvalue \(\ell\) on the
geometric rotations and is zero on their orthogonal complement in its
isotypic space. All other blocks have zero first-moment interaction.

With \(\phi,\Phi\) the standard Gaussian density and CDF, define
\[
 I_j(k)=\int_{\mathbb R}\phi(y)^j\Phi(y)^{k-j}\,dy .
\]
The exact facet and variance formulas are
\[
 v_k=\sqrt2 I_2(k),\qquad
 1-\gamma_k=\frac{k}{3}\frac{I_4(k)}{I_2(k)}.
 \tag{VH1}
\]

To see this, add the free common-coordinate direction to the active
space and sample independent standard Gaussians \(Z_1,\dots,Z_k\) on
the resulting tie plane. Write \(Y=Z_1=Z_2\), which has variance
\(1/2\) under the normalized Gaussian law on that plane. The facet
condition is precisely that the other independent coordinates are at
most \(Y\). Restricting surface measure to the plane contributes the
factor \(\phi(0)\) from the orthogonal standard Gaussian coordinate.
The product of \(\phi(0)\) with the density of \(Y\) is
\(\sqrt2\,\phi(y)^2\). This gives \(v_k\) in (VH1) and the probability
density \(\phi(y)^2\Phi(y)^{k-2}/I_2(k)\) for \(Y\) on the facet.

The added common coordinate does not affect any linear functional in
\(H^0\). Given \(Y=y\) in the facet model, \(Z_3,\ldots,Z_k\) are
independent standard Gaussians truncated above at \(y\).
A unit contrast of them has mean zero since its coefficients sum
to zero. Its variance is the common truncated variance
\(1-y\lambda(y)-\lambda(y)^2\), where
\(\lambda(y)=\phi(y)/\Phi(y)\). This follows from integration
by parts in the first two moments of the truncated Gaussian.
Averaging yields
\[
 1-\gamma_k
   =\frac{\int y\phi(y)^3\Phi(y)^{k-3}\,dy + I_4(k)}
         {I_2(k)}.
\]
The integral of \((\phi^3\Phi^{k-3})'\) vanishes by the Gaussian
decay at both ends, so the numerator is
\((1+(k-3)/3)I_4(k)\). This proves (VH1) for every \(k\ge4\).

## The asymptotic gap

Set \(z_k=\Phi^{-1}(1-1/k)\). For each fixed positive integer \(j\),
\[
 I_j(k)\sim \frac{\Gamma(j)z_k^{j-1}}{k^j}
 \quad (k\to\infty).
 \tag{VH2}
\]
Indeed, put \(Q=1-\Phi\) and use \(u=kQ(y)\) as the integration
variable. With \(y_k(u)=Q^{-1}(u/k)\), the rescaled integral is
\[
 \frac{k^j I_j(k)}{z_k^{j-1}}
   =\int_0^k
      \left(\frac{k\phi(y_k(u))}{z_k}\right)^{j-1}
      (1-u/k)^{k-j}\,du .
\]
For each fixed \(u>0\), the Mills bounds give
\(y_k(u)/z_k\to1\) and
\(\phi(y_k(u))/(y_k(u)Q(y_k(u)))\to1\).
The integrand therefore tends to \(u^{j-1}e^{-u}\).

For a dominating function, suppose first \(u\le k/2\), so
\(y_k(u)\ge0\). The bounds
\(\phi(y)\le C(1+y)Q(y)\) and \(Q(y)\le e^{-y^2/2}\) for
\(y\ge0\), together with \(z_k\asymp\sqrt{\log k}\), give
\[
 \frac{k\phi(y_k(u))}{z_k}
   \le C u\bigl(1+\sqrt{\max\{0,\log(1/u)\}}\bigr)
\]
for all sufficiently large \(k\). For \(k/2<u<k\) the same
bound follows from \(\phi(y)\le\phi(0)\) and \(k<2u\).
The other factor in the integrand is at most \(e^{-u/2}\)
when \(k\ge2j\). This majorant is integrable at zero and infinity,
and the integrand can be set to zero on \(u\ge k\).
Dominated convergence proves (VH2).

Let \(D\) be the first-moment map on signed facet speeds as in (S82).
Equip the balanced speeds with their Gaussian surface \(L^2\) norm
and let \(\mathcal R_k\) be the geometric rotation space.
Combining (VH1)–(VH2) and \(z_k^2\sim 2\log k\) gives
\[
 1-\gamma_k=(2+o(1))z_k^2/k.
\]
For sufficiently large \(k\), the two-row interaction eigenvalue
\((k-1)v_k\gamma_k\) exceeds the standard bound \(\ell/2\).
The exact balanced first-moment gap is therefore
\[
 \inf_{\substack{f\perp\mathcal R_k,\ f\ {\rm balanced}\\ \|f\|=1}}
  \bigl(\ell\|f\|^2-\|Df\|_F^2\bigr)
  =\ell\bigl(1-(1-1/k)\gamma_k\bigr)
  =(4+o(1))\,\ell\frac{\log k}{k}.
 \tag{VH3}
\]
The two-row field from (S82) attains the infimum. This improves the
first-moment lower bound \(\ell/k\) by a factor asymptotic to
\(4\log k\).

## The full negative paired gap

Use \(D_r,K_r,C,\mathcal R_k\) from
[(NL1)–(NL2)](NEGATIVE_SIMPLEX_INTERFACE_LIFT.md#weighted-hessian-and-compactness),
with \(D_r\) the form defined by the positive local weight \(a_r\)
and \(K_r\) the interaction at correlation \(-r\).
There are absolute \(c>0\) and \(k_0\) such that for all \(k\ge k_0\)
and \(0<r\le1/(k-1)\), the active symmetric and antisymmetric forms
satisfy

\[
\begin{aligned}
 (D_r+K_r)(f,f)&\ge c\frac{\log k}{k}
      \operatorname{dist}_{D_r}(f,\mathcal R_k)^2,\\
 (D_r-K_r)(f,f)&\ge \tfrac12D_r(f,f)\qquad(Cf=0).
\end{aligned}
\tag{VH4}
\]

In all asymptotic estimates below, constants and errors are uniform
for \(0<r\le1/(k-1)\).

### The reciprocal local weight

Let \(\mu\) be normalized Gaussian measure on \(F_{12}\), and let
\(U=p\cdot x\), where \(p\) is the unit vector in the direction
\(p_0=(e_1+e_2)/2-\mathbf1/k\). Put \(m_k=\mathbb E_\mu U\).
The [standard row identity (S49)](SIMPLEX_INTERFACE_SYMMETRY.md#every-cell-count-the-standard-sector-is-a-scalar-positive-kernel)
holds at negative correlation, with the row equal to \(a_r/r\) by
(NL6). Absorb the Gaussian square in the input variable by the
dilation \(y=\sigma Y\), \(\sigma=\sqrt{1-r^2}\). The result is

\[
 \frac{a_r(x)}{r\ell}
  =\frac1{\sigma}e^{-\beta|x|^2}\mathbb E
                   \exp\bigl(-(r/\sigma)x\cdot GY\bigr),
 \qquad
 \beta=\frac{r^2}{2\sigma^2}.
\]

Here \(Y\sim\mu\) and the independent variable \(G\) is the identity
with probability \(2/k\) and \((2j)\) with probability \(1/k\) for
each \(j=3,\ldots,k\). The orthogonal projection of
\(\mathbb E_G Gp\) onto the supporting plane is \(p/2\).
Indeed the projected adjacent map is \(I-w_jw_j^T\) with
\(w_j=(e_1+e_2)/2-e_j\), and
\(\sum_j w_j=kp_0,\ w_j\cdot p_0=1/2\).
By untied-coordinate symmetry \(\mathbb E Y=m_kp\).
Jensen's inequality gives, on the facet,

\[
 A_r=\frac{a_r}{r\ell}\ \ge\ \sigma^{-1}e^{-X},
 \qquad X=\beta|x|^2+\frac{r m_k}{2\sigma}U .
 \tag{VH5}
\]

In the tied-height model of (VH1), let \(Z_{\mathrm{av}}\) be the mean
of the \(k-2\) remaining coordinates. The projection to the active
space gives \(U=\sqrt{2(k-2)/k}\,(Y-Z_{\mathrm{av}})\).
Integration of the derivative of the tied-height density gives
\(\mathbb E_\mu\lambda(Y)=2\mathbb E_\mu Y/(k-2)\).
The conditional mean of \(Z_{\mathrm{av}}\) is \(-\lambda(Y)\). Hence

\[
 m_k=\sqrt{2(k-2)/k}\,(1+2/(k-2))\mathbb E_\mu Y,
 \qquad m_k^2=(4+o(1))\log k.
 \tag{VH6}
\]

For the asymptotic, put \(z=\sqrt{2\log k}\).
The normalized tied-height density has denominator
\(I_2(k)\asymp\sqrt{\log k}/k^2\) by (VH2).
On \(Y\le(1-\epsilon)z\), the Mills bounds give
\(\Phi(Y)^{k-2}\le \exp(-k^{c_\epsilon})\) for large \(k\) and
each \(0<\epsilon<1\). On \(Y\ge(1+\epsilon)z\), the factor
\(\phi(Y)^2\) alone bounds the unnormalized integral of
\(1+Y^2\) by a polynomial in \(z\) times \(k^{-2(1+\epsilon)^2}\).
After division by \(I_2(k)\), both exterior bounds tend to zero.
Thus \(\mathbb E_\mu Y/z\to1\), as used in (VH6).

For every unit linear functional on the supporting plane, the
Gaussian facet law has centered exponential moment at most
\(e^{t^2/2}\) at argument \(t\). To see this, tilt by the linear
exponential. Completing the square shows that its tilted variance
is a conditional Gaussian variance on a translate of the same
convex cone. The convex Gaussian covariance bound used in
[(S77)](SIMPLEX_INTERFACE_SYMMETRY.md#the-exact-initial-standard-spectral-minimum)
applies also to the translated cone. Integrate the variance bound
twice in the log moment generating function. In addition,
\(|x|^2\) under \(\mu\) is exactly \(\chi^2_{k-2}\), since the cone
restricts only the angular variable. It follows by Hölder and the
chi-square exponential moment that
\(\|X\|_{L^4(\mu)}=O(\log k/k)\) and
\(\|e^X\|_{L^4(\mu)}=O(1)\). For instance
\(\mathbb E e^{8\beta|x|^2}=(1-16\beta)^{-(k-2)/2}=O(1)\),
and the linear exponential is bounded by the preceding centered
estimate and (VH6).

If \(V\) is any unit contrast in \(H^0\), its mean is zero and
\(\mathbb E V^2=\gamma_k\). Its centered moments of every fixed
order, and those of \(U-m_k\), are bounded as above. Thus
\(\mathbb E[V^2|x|^2]=O(k)\) and
\(\mathbb E[V^2 U]=m_k\gamma_k+O(1)\).
Using \(U\ge0\), hence \(X\ge0\), expand \(e^X\) to first order.
The quadratic remainder after multiplication by \(V^2\) has
expectation \(O(\|X\|_4^2\|V^2\|_4\|e^X\|_4)\). With
(VH1)–(VH2), (VH5) and (VH6), this gives

\[
\begin{aligned}
 \mathbb E_\mu[V^2/A_r]
 &\le \gamma_k+\frac{r m_k^2}{2\sigma}\gamma_k
                       +o(\log k/k)\\
 &\le 1-(2+o(1))\log k/k,\qquad
 \mathbb E_\mu[(U-m_k)^2/A_r]\le1+o(1).
\end{aligned}
\tag{VH7}
\]

In the first line the radial correction is \(O(1/k)\) and the
centered linear correction is \(O(\sqrt{\log k}/k)\). For the
second estimate, \(\mathbb E(U-m_k)^2\le1\), and the difference
with the exponential weight is \(O(\|X\|_4\|e^X\|_4)\) using
the bounded fourth moment of \(U-m_k\).

### The interaction beyond degree one

Write \(E_j\) for the nonnegative squared Hermite coefficient forms
in [(S11)](SIMPLEX_INTERFACE_SYMMETRY.md#interaction-sparsity-and-the-sign-sector-hermite-floor),
so \(K_r=\sum_{j\ge0}(-r)^j E_j\).
For every speed on the active facets,

\[
 \frac{r^2 E_2(f)}{D_r(f,f)}=O(r\log^2 k),\qquad
 \left|\sum_{j\ge3}(-r)^j E_j(f)\right|
             \le O(r^2/\ell)\,D_r(f,f).
 \tag{VH8}
\]

For the first estimate, a unit degree-two Gaussian Hermite
polynomial \(P\) restricts to a polynomial on the supporting plane
of Gaussian \(L^2\) norm at most \(1\). Decompose in Hermites of the
normal coordinate, which is set to zero. The normal degree-one term
vanishes, and the degree-two term becomes a constant of magnitude at
most \(1/\sqrt2\) times its coefficient, orthogonal to the
surviving tangential degree-two terms. Gaussian hypercontractivity
from \(L^2\) to \(L^q\), or its immediate consequence for degree-two
Hermites, gives \(L^q\) norm at most \(Cq\) on the plane.
The normalized Gaussian probability of the facet inside that plane
is \(v_k/\phi(0)\asymp\sqrt{\log k}/k^2\).
Hölder with \(q\) of order \(\log k\) therefore bounds the
conditional \(L^4(\mu)\) norm of \(P\) by \(O(\log k)\).
[Nelson's Gaussian hypercontractivity, as stated in Gross, Theorem 1.1](https://doi.org/10.1007/BF02392573).
By (VH5) and the bounded exponential moment,
\(\int_{F_{ij}}P^2 a_r^{-1}d\nu_{ij}
 \le C(v_k/(r\ell))\log^2 k\).
Duality in the degree-two coefficient vector and Cauchy--Schwarz
on each facet give this same bound times its \(D_r\)-energy for the
squared norm of its coefficient vector. The label \(i\) adds at
most \(k-1\) such vectors, so Cauchy--Schwarz in this sum and then
summation over labels give the factor \(2(k-1)\). This proves the
first estimate in (VH8).

For the second estimate fix \(b=1/2\) as a positive comparison
correlation. By (VH5), (VH6) and Young's inequality,
\[
 a_r(x)^{-1/2}\le C(r\ell)^{-1/2}e^{\delta|x|^2/2},
 \qquad \delta=\beta+1/k=O(1/k).
\]
Indeed the extra constant from the linear term of \(X\) is at most
\(\exp[O(k r^2 m_k^2)]\), which is bounded. In the operator on
Lebesgue surface measure the absolute positive-\(b\) kernel,
after multiplying by \(e^{\delta|x|^2/2}\) at both arguments, has
cross coefficient \(b/(1-b^2)\), quadratic coefficient
\((1+b^2)/(4(1-b^2))-\delta/2\), and prefactor
\((2\pi)^{-n/2}(1-b^2)^{-n/2}\), \(n=k-1\).
Choose \(b'>b\) so that

\[
 \frac{b}{b'}(1+(b')^2)=1+b^2-2(1-b^2)\delta,\qquad
 \lambda^2=\frac{b(1-(b')^2)}{b'(1-b^2)} .
\]

For large \(k\) the first equation has a solution \(b'=b+O(1/k)\)
in a fixed compact subinterval of \((0,1)\), since its left side
is strictly decreasing with derivative bounded away from zero near
\(b\). The quadratic and cross coefficients match those of the
positive-\(b'\) kernel at \((\lambda x,\lambda y)\).
Every facet is dilation invariant. Unitary dilation in surface
dimension \(n-1\) therefore multiplies the operator norm by the
prefactor ratio

\[
 \lambda^{-(n-1)}\Big(\frac{1-(b')^2}{1-b^2}\Big)^{n/2}
  =\Big(\frac{b'}b\Big)^{(n-1)/2}
    \Big(\frac{1-(b')^2}{1-b^2}\Big)^{1/2}=O(1).
\]

The absolute row bound in
[(S81)](SIMPLEX_INTERFACE_SYMMETRY.md#a-quantitative-gap-for-all-sufficiently-deep-sectors)
is uniform at \(b'\). Thus the positive-\(b\) interaction after
normalization by \(a_r^{-1/2}\) has norm \(O((r\ell)^{-1})\).
Since every \(E_j\) is nonnegative, the absolute alternating
tail from degree three is bounded by \((r/b)^3\) times the
positive-\(b\) form, which proves (VH8). All comparisons can be
made first for compactly supported speeds, then extended by density
in the weighted spaces. The proof uses the positive comparison
kernel only for an upper bound on Hermite terms.

### The two signs and the low-degree blocks

The antisymmetric constraint is \(Cf=0\), so \(E_0(f)=0\).
In \(D_r-K_r\), degree one contributes \(+rE_1\).
By (VH2), \(\ell\asymp\sqrt{\log k}/k\). Thus (VH8) shows that the
sum of unfavorable terms is \(o(D_r)\) and gives the second part
of (VH4).

For \(D_r+K_r\), \(E_0\) and \(r^2E_2\) have favorable sign.
The adverse tail bound in (VH8) is
\(O(1/(k\sqrt{\log k}))D_r=o((\log k/k)D_r)\).
The only blocks containing degree-one energy are the standard
representation, the rotation hook, and \((k-2,2)\), by (S82).
Their forms are orthogonal since the operators and \(a_r\) are
\(S_k\)-equivariant.

- In the standard scalar model (S49), with \(h\) untied-invariant,
  normalize by its common representation factor and \(r\ell v_k\).
  The local, constant, and degree-one contributions to the form
  are \(\mathbb E_\mu A_r h^2\),
  \(r^{-1}(\mathbb E_\mu h)^2\), and
  \(-\frac12(\mathbb E_\mu Uh)^2\), respectively.
  The last identity follows by differentiating the scalar kernel
  at zero, since its linear term on the tie axis is \(k/2\).
  By \((a+b)^2\le(5/4)a^2+5b^2\), and the centered estimate
  in (VH7), the degree-one term has magnitude at most
  \((5/8+o(1))\mathbb E_\mu A_r h^2+
      (5/2)m_k^2(\mathbb E_\mu h)^2\).
  The constant term absorbs the second term because
  \(1/r\ge k-1\). This leaves a gap bounded below by a positive
  constant.
- On \((k-2,2)\), (S82) gives the rank-one multiplicity
  eigenfield \(P_{H^0}x\) with unweighted interaction eigenvalue
  \((k-1)v_k\gamma_k\). Cauchy--Schwarz with weight \(a_r\)
  bounds the normalized interaction \(rE_1/D_r\) by
  \((1-1/k)\mathbb E_\mu[V^2/A_r]\): the squared length of
  \(P_{H^0}x\) is summed over an orthonormal basis of \(H^0\),
  whose terms give the same integral by untied-label symmetry.
  The first estimate in (VH7) and the adverse tail bound
  give the gap of order \(\log k/k\).
- On the rotation hook, the normalized operator for \(rE_1\)
  has rank at most \(\dim\mathcal R_k\), by (S82) and boundedness
  of the weighted coefficient map. The constant term vanishes
  there because \(C\) maps into the standard representation.
  The remaining interaction has norm \(o(1)\), by (VH8).
  The \(\dim\mathcal R_k\) independent exact rotations already
  have normalized interaction eigenvalue \(-1\) by (NL6).
  If the weighted orthogonal complement had an eigenvalue whose
  modulus exceeds the remainder bound, its eigenspace together
  with the rotation space would contain a nonzero vector in the
  kernel of the rank-\(\dim\mathcal R_k\) part. The normalized
  interaction could then have norm at most the remainder bound
  on that vector, a contradiction. Compactness from (NL5)
  shows the complementary interaction has norm \(o(1)\).

The remaining blocks have no adverse term below degree three.
Combining the bounds proves (VH4).

## Ambient dimension and measurable pairs

The [dimension lift (NL6)–(NL7)](NEGATIVE_SIMPLEX_INTERFACE_LIFT.md#dimension-lift-for-the-paired-form)
preserves the lower bound \(c\log k/k\) in the exact local
weight, after decreasing \(c\) if necessary. The full paired
kernel consists of the simultaneous rotations, including the
rotations between the active and unused coordinates.
For every sufficiently large fixed \(k\) and every fixed
\(0<r\le1/(k-1)\),
[(NL8)](NEGATIVE_SIMPLEX_INTERFACE_LIFT.md#nearby-measurable-pairs-in-every-ambient-dimension)
gives a dimension-uniform strict local comparison for measurable
fractional pairs with matched free masses, with equality on the
common simplex orbit. One may take the two assignments to be equal.
The measurable neighborhood in that implication can depend on \(k,r\).
The [polynomial local radius](NEGATIVE_VOTING_POLYNOMIAL_NEIGHBORHOOD.md)
quantifies this neighborhood and supplies the local step in the
[global same-partition comparison](NEGATIVE_VOTING_INTERVAL.md#exact-gaussian-simplex-minimization).

## A first global transition requires a distinct tied pair

For every sufficiently large fixed \(k\), let \(R_k=1/(k-1)\)
and let \(v(r)\) be the minimum of the matched-mass paired
functional \(P_r\), defined in
[the small negative result](NEGATIVE_SMALL_CORRELATION.md#exact-matched-mass-minimization-near-zero).
The minimum is the same in every finite \(d\ge k-1\).
Consequently, if there is a strict global improvement at some
\(r\le R_k\), the first strict-improvement boundary in \((0,R_k)\)
has an ordinary minimizing pair different from every common
simplex pair. Its value equals the simplex value there. In fact
the simplex is a uniform strict local minimizer along the voting
interval, in the sense that for each fixed sufficiently large
\(k\) a single orbit-distance radius works for all \(d\ge k-1\)
and all \(0<r\le R_k\).

To justify dimension reduction, first fix finite \(d\).
The set of fractional matched-mass pairs is weakly compact in the
Gaussian \(L^2\) product. The operator \(T_{-r}\) is compact there,
so \(P_r\) is jointly weakly continuous. Ordinary pairs have the
same infimum, by independently rounding both assignments on
increasing finite Gaussian subdivisions while preserving all
cell probabilities in every atom. For \(r>0\), reflecting the
second ordinary partition converts the optimization to Heilman's
[Problem 7.1, Lemma 7.3 and Theorem 7.9](https://arxiv.org/html/2209.11216#S7).
The common masses are matched as required by that result.
Theorem 7.9 puts an ordinary minimizing pair in a common
\((k-1)\)-dimensional subspace, up to rotation and null sets.
Extending pairs cylindrically gives the reverse comparison of
minimum values.

In dimension \(k-1\), the minimum \(v(r)\) is continuous for
\(0<r<1\). The reason is that \(\sup_{j\ge0}|r^j-s^j|\to0\)
as \(r\to s<1\), and objective differences are thus small
uniformly on the feasible set. For small \(r>0\) the paired
fixed-cell theorem gives equality of \(v(r)\) with the simplex
value and rigidity. For \(r\) outside that initial interval,
the measurable local comparison from (VH4) covers a correlation
neighborhood of each point of the remaining compact parameter
interval, even at \(R_k\). A finite cover supplies the single
positive orbit-distance radius. In the small initial interval
the global paired inequality already gives strictness off the
common simplex orbit.

Suppose \(r_0\) is a correlation with a strict improvement
and let \(r_*\) be the infimum of such \(r\le r_0\). Continuity
gives no strict improvement at \(r_*\) and \(r_*<r_0\).
Take ordinary minimizing pairs at correlations decreasing to
\(r_*\) from the strict-improvement set, in dimension \(k-1\).
A subsequence converges weakly and its limit minimizes \(P_{r_*}\)
by compactness of \(T_{-r_*}\) and the operator-norm continuity.
If the limit were a common simplex indicator pair, weak convergence
would imply label \(L^1\) convergence. For example,
\(\sum_i\|F_i-1_{T_i}\|_1=2-2\sum_i\langle F_i,1_{T_i}\rangle\)
for any simplex-valued assignment \(F\). The local comparison
would then exclude strict improvements for the subsequence.
Thus the limit is a distinct tied pair.

If needed, obtain an ordinary distinct tied pair as follows.
At fixed common mass vector and fixed second assignment, the
first objective is linear on a weakly compact convex set.
A minimizing extreme point is an ordinary partition. Indeed, if
two labels both carry positive fractional weight on a set of
positive measure, one can exchange a sufficiently small signed
density between those labels, of integral zero, to express the
assignment as a nontrivial midpoint. Use this extremal argument
first on the limit pair's first side and then on its second.
If every ordinary minimizing pair were a common simplex pair
\((A,A)\), the first replacement would be its first assignment
\(A\). With \(A\) fixed, the original limit's second assignment
must then be \(A\), by strict pointwise optimality of the
reflected simplex scores away from the facet planes. The same
argument on the first side would make the original limit
\((A,A)\), a contradiction. This establishes the ordinary
distinct tie. It identifies the first-transition obstruction for the global
matched-mass paired inequality. The same-partition proof uses
the local gap through (PN1), with exponential localization of
same-partition competitors.
