# The adaptive surface variational problem for three Gaussian cells

Fix \(0<r<1\). For an ordinary planar Gaussian partition,
allow the direction in the surface certificate to vary
measurably along each interface. The resulting functional
is lower semicontinuous in Gaussian label distance and
attains its infimum over balanced ordinary partitions.
Its value is strictly above the Gaussian perimeter of
every individual finite-energy partition. Its minimizers
approach the regular simplex orbit as \(r\downarrow0\).
At an actual noise-stability maximizer its square is the
correlation derivative. Thus rigidity of minimizers of
this separate surface problem would rule out right-contact
ties. The simplex minimization claim remains open.

The freedom to vary directions matters for the existence
argument. Each dual test has interface weights equal to
distances between smooth label vector fields.
These weights satisfy the triangle inequality, so a thin
intermediate label cannot reduce that test in the limit.
[Thin-layer obstruction for the fixed-direction energy](THREE_CELL_PAIR_SURFACE_WETTING_OBSTRUCTION.md)

## A dual formula with metric surface weights

Let \(A\) have locally finite perimeter and finite total
Gaussian interface measure. Write \(\nu_{ij}\) for its
Gaussian pair interface measures, \(d_{ij}=e_i-e_j\)
for the corresponding label vectors, and \(N_{ij}\) for
the normal from \(i\) into \(j\). Interfaces are counted
once. For measurable \(\mathbb R^2\)-valued fields
\(u_{ij}\) with \(|u_{ij}|\le1\), put
\[
 \mathcal F_r(A)=\sup_u\|T_{\sqrt r}m^u\|_2,\qquad
            m^u=\sum_{i<j}d_{ij}\otimes u_{ij}\,\nu_{ij}.
\]
The noise operator on a measure uses the Mehler kernel
relative to Gaussian measure. The norm sums the squared
matrix coordinates. Allow the value \(+\infty\) and assign
that value when the partition has infinite Gaussian
perimeter. Up to the square this is the direction
supremum in (WC1).
[Adaptive directions at an actual maximizer](THREE_CELL_WEIGHTED_PERIMETER_CONTACT.md#a-direction-supremum-for-the-correlation-derivative)

The same identification holds at the simplex at every
positive correlation, including outside proved global
intervals: the first extra-coordinate form is nonnegative
by (S37), and its normal direction field attains the
translation equality in the proof of (WC1).
[All-correlation simplex interface gap, S37](SIMPLEX_INTERFACE_SYMMETRY.md#reflected-kernels-close-the-three-cell-spectral-problem)

Take a smooth compactly supported matrix-field test
\(G\) with \(\|G\|_2\le1\), the Hilbert norm summing the
three \(\mathbb R^2\)-valued rows. Put
\(V_i=T_{\sqrt r}G_i\). Then
\[
 \mathcal F_r(A)
   =\sup_{\substack{G\in C_c^\infty\\\|G\|_2\le1}}
                   \sum_{i<j}\int |V_i-V_j|\,d\nu_{ij}.
                                                               \tag{AS1}
\]
For each \(u\), the function \(T_{\sqrt r}m^u\) is
continuous and bounded on compact subsets. This uses the
uniform bound on the Mehler kernel in the source variable
when the observation stays in a compact set. Its extended
\(L^2\) norm is the supremum of its pairings against the
displayed tests, by cutoff and \(L^2\) density on compact
sets. Interchange the two suprema, use Mehler symmetry
and Fubini for the compact test, and choose each
\(u_{ij}\) in the direction \(V_i-V_j\).
This proves (AS1), including the infinite case.

For every fixed \(G\) in (AS1) there is a local calibration
identity:
\[
 \sum_{i<j}\int |V_i-V_j|\,d\nu_{ij}
   =\sup_{\substack{L\in C_c^\infty(\mathbb R^2;\mathbb R^{2\times2})\\
                           \|L(x)\|_{\rm op}\le1}}
       \sum_i\int_{A_i}\operatorname{div}_\gamma(LV_i)\,d\gamma .
                                                               \tag{AS2}
\]
Gaussian integration by parts makes the right integrand
on an interface \((L(V_i-V_j))\cdot N_{ij}\).
Its upper bound is the integrand on the left. To approach
this bound, on the \(i,j\) interface set \(L\) equal to the
rank-one contraction sending the unit vector
\((V_i-V_j)/|V_i-V_j|\) to \(N_{ij}\), and set it to zero
if the difference vanishes. The pair interfaces are mutually
disjoint up to their common surface null sets.
On compact subsets the resulting bounded Borel field can
be approximated in the summed Gaussian surface \(L^1\)
norm by continuous fields taking values in the operator
norm unit ball. For example, use density of continuous
functions for finite Radon measures, project pointwise
onto the closed convex norm ball in the Euclidean matrix
norm, and multiply by a compact cutoff. Mollification
preserves the convex constraint and approximates the
continuous fields uniformly. Exhaust the plane. The
vectors \(V_i\) are bounded and smooth, since \(G\) is
smooth and bounded, so the obtained surface integrals
approach the left side. This proves (AS2).

For fixed \(G,L\) the bulk integrals on the right are
continuous in Gaussian label distance, since their compactly
supported integrands are bounded. Taking the suprema
gives lower semicontinuity whenever the limiting
partition has finite Gaussian perimeter. The perimeter
bound below handles limits outside this class.

## Existence and a strict perimeter bound

Write \(P_\gamma(A)=\sum_{i<j}\nu_{ij}(\mathbb R^2)\).
Choose three constant row vectors \(G_i^0\) forming a
centered regular triangle of side one. Their total
squared norm is one. The directions
\(u_{ij}=G_i^0-G_j^0\) are unit, so the corresponding
smoothed tensor pairs with the constant test \(G^0\) to
give \(P_\gamma(A)\) by Gaussian invariance of the kernel.
Consequently \(\mathcal F_r(A)\ge P_\gamma(A)\).

If \(0<P_\gamma(A)<\infty\) and equality held, the smoothed
tensor for these directions would have equality in
Cauchy–Schwarz with \(G^0\). It would be the constant
field \(P_\gamma(A) G^0\). Injectivity of Mehler
smoothing on finite signed measures would identify
this interface tensor with the nonzero absolutely
continuous measure \(P_\gamma(A)G^0\,d\gamma\), impossible.
Injectivity follows by writing the smoothing as the
Laplace transform of a Gaussian-damped finite measure.
Analytic continuation and Fourier injectivity then recover
each coordinate. Thus
\[
        \mathcal F_r(A)>P_\gamma(A)
                  \quad\text{if }0<P_\gamma(A)<\infty.                 \tag{AS3}
\]
Infinite values satisfy the inequality when the perimeter
is finite.

In particular a minimizing sequence of balanced ordinary
partitions has bounded Gaussian perimeter because the
simplex gives a finite comparison value. On each fixed
disk this bounds the ordinary perimeters, since the
Gaussian density there is bounded below. Local
\(BV\) compactness followed by exhaustion and the
Gaussian tail bound gives a subsequence converging
in Gaussian label distance. The limit is an ordinary
balanced partition of finite Gaussian perimeter.
By (AS1)–(AS2) it minimizes \(\mathcal F_r\).
The same perimeter control and local compactness give
lower semicontinuity also when the proposed limit has
infinite perimeter. In fact any finite liminf would force
finite perimeter for the limit.
[Gaussian perimeter compactness and ordinary limits](FINITE_REFLECTION_TEST_FOR_MAXIMIZERS.md#uniform-gaussian-perimeter-and-a-hermite-tail)

For completeness, if \(\mathcal F_r(A)<\infty\), the
supremum over \(u\) is attained. Indeed the constant
equilateral direction energy is at least
\(\|T_{\sqrt r}\nu_{\rm sum}\|_2^2\), where
\(\nu_{\rm sum}=\sum_{i<j}\nu_{ij}\). To verify this,
write \(q_i\) for the three smoothed nonnegative pair
densities. Their equilateral quadratic energy is
\(\frac12(\sum_i q_i)^2+\frac32\sum_i q_i^2
\ge(\sum_i q_i)^2\) pointwise. Positivity and
dominated convergence now give
\(\|T_{\sqrt r}(1_{|x|>R}\nu_{\rm sum})\|_2\to0\).
This uniformly bounds the tail of the tensor smoothing
up to the fixed incidence norm. On compact source
sets the smoothing operator is Hilbert–Schmidt from
surface \(L^2\) to Gaussian \(L^2\), since
\(\int M_{\sqrt r}(z,x)^2\,d\gamma(z)=M_r(x,x)\)
is bounded there. Weak convergence of a subsequence
of bounded \(u\)'s thus gives strong convergence of
the smoothed tensors, with the norm constraint retained.
This proves attainment.

## Connection with a right endpoint

Let \(T\) be the centered regular simplex and let \(P_*\)
be its Gaussian perimeter. By the Gaussian multi-bubble
theorem \(P_\gamma(A)\ge P_*\) for every balanced ordinary
partition, with equality only for the simplex orbit.
The minimizing value of \(\mathcal F_r\) is strictly
above \(P_*\) by (AS3) and attainment.
If \(A_r\) is any such minimizer, then
\[
 P_*\le P_\gamma(A_r)
       \le \mathcal F_r(A_r)\le\mathcal F_r(T)
         =\sqrt{S'_r(T)}\longrightarrow P_*
                    \quad (r\downarrow0).
\]
Hence perimeter compactness and the equality case of
the Gaussian multi-bubble theorem give convergence of
these minimizers to the simplex orbit in Gaussian label
distance. The final limit also follows from the
degree-one Hermite energy of the simplex.
[Gaussian multi-bubble input and simplex perimeter normalization](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md#external-endpoint-and-application-audit)

At an actual balanced planar maximizer \(A\),
(WC1) gives \(\mathcal F_r(A)^2=S'_r(A)\).
At a tie approached from the right by simplex-optimal
correlations the fixed-partition derivative satisfies
\(S'_r(A)\le S'_r(T)\). Therefore a proof that the
simplex uniquely minimizes \(\mathcal F_r\) at the
remaining correlations would exclude these ties.
Attainment and lower semicontinuity establish the
variational problem for ordinary partitions, without
establishing its minimizer classification.

## Compactly supported interfaces and the minimum profile

The infimum can also be approached by partitions that
are constant outside a disk. In particular the minimum
profile \(r\mapsto\min_A\mathcal F_r(A)\) is continuous
and strictly increasing on \(0<r<1\).

For the approximation fix a balanced \(A\) with finite
\(\mathcal F_r(A)\). Make its exterior beyond \(B_R\)
one label. Away from the cutoff circle this preserves
the old pair interfaces inside and deletes the old ones
outside. Any additional interfaces at the cutoff have
scalar surface measure dominated by the Gaussian area
measure of the circle. Its smoothed \(L^2\) norm tends
to zero: in the squared Mehler integral over the circle,
the product Gaussian surface density times the kernel
is bounded by \(C_r\exp(-R^2/(1+r))\), and the Euclidean
length factor is at most \(4\pi^2R^2\).

Restore the original masses in a fixed compact set.
In each positive-mass label choose two distinct
Lebesgue density points, all of them distinct across
labels. Any surplus label can supply at most two
deficit labels. Transfer the portion of the surplus
label in a small ball about its corresponding density
point. The Gaussian volume of that portion is a
continuous function of the radius, is zero at radius
zero, and is positive at every small positive radius.
As the cutoff errors tend to zero the required ball
radii can therefore tend to zero. They are disjoint
and unaffected by the exterior cutoff for large \(R\).
The extra scalar surface measure inside a repair ball
is bounded by the old total interface measure there
plus the sphere's Gaussian measure; both tend to
zero. For the old measure this follows because the
finite perimeter interface measures in the plane
give zero mass to any point. These repairs change
only such small balls and preserve exact balance.

In comparing the resulting partition with \(A\),
copy any new direction field to the corresponding
unchanged interfaces of \(A\) and use zero elsewhere
on its old interfaces. Differences are supported only
on the added or relabelled interface measures.
Triangle inequality in the smoothed \(L^2\) norm
bounds their cost by a fixed incidence multiple of
the smoothed scalar measures just estimated.
On the fixed compact set for the repairs, the squared
Mehler integrals are bounded by a constant times the
squared mass, since the kernel is bounded there.
Thus the resulting adaptive functional has limsup
at most \(\mathcal F_r(A)\), and its label convergence
to \(A\) and lower semicontinuity give convergence of
the energy as well.

For compact interfaces \(\mathcal F_r(A)\) is
continuous in \(r\): the kernels in the squared norms
converge uniformly on the compact product, with an
error bounded uniformly for all the admissible
directions. The compact-interface approximation
therefore proves upper semicontinuity of the minimum
profile. For lower semicontinuity take optimizing
partitions along a convergent sequence of positive
correlations. Their perimeter bounds are uniform by
comparison with the simplex. Extract a strong label
limit and use joint lower semicontinuity in \(r,A\),
which follows directly from (AS1)–(AS2): for every
fixed compact \(G,L\), the divergence fields in (AS2)
converge uniformly on the support of \(L\) as \(r\)
varies. This proves continuity of the minimum, and
the strong limits just obtained are also minimizing.

Finally, for \(r<s\), choose at any balanced
partition a direction field attaining its finite
\(\mathcal F_r(A)\). Its smoothed tensor is nonconstant
by the same finite-measure injectivity argument
as in (AS3); it also has positive norm. The Mehler
semigroup and strict \(L^2\) contraction on
nonconstant functions give for this field a strictly
larger smoothed norm at \(s\), unless that norm is
infinite, in which case the strict inequality is
immediate. In particular let \(A\) minimize the
finite value at \(s\). Then \(\mathcal F_r(A)\) is
finite and the strict inequality at this \(A\)
shows strict increase of the minimum profile.

## Tagged-interface convergence at simplex energy

At a fixed positive \(r\), suppose \(A_n\) converges in Gaussian
label distance to \(T\) and
\(\mathcal F_r(A_n)\to\mathcal F_r(T)\).
Then each tagged pair surface measure converges locally weakly
to the corresponding simplex pair measure. In particular the
limiting surface cannot include a hidden third-label layer along
a reference pair ray.

Use the fixed Hilbert test
\[
           G_*=\frac{T_{\sqrt r}m_T^{N^T}}{\mathcal F_r(T)}
\]
of norm one. Formula (AS1) also bounds the metric surface
test for this \(G_*\), although it is not compactly supported.
Indeed \(\|T_{\sqrt r}\nu_{\rm sum}\|_2\) is finite whenever
the adaptive energy is finite, by its equilateral test in
the attainment proof. Hence the integrals with \(|G_*|\)
are absolutely convergent by Cauchy–Schwarz and the positive
kernel. Fubini and choice of directions prove the same
dual bound. At \(T\) this test has equality, since
the normal field attains its adaptive norm and pairs
with \(G_*\) to give that norm.
Formula (AS2) still applies locally without any global
bound on \(T_{\sqrt r}G_*\), since the matrix field there
has compact support.

Write \(V_i=T_{\sqrt r}(G_*)_i\). The weighted surface
mass \(\omega_n=\sum_{i<j}|V_i-V_j|\nu_{ij}(A_n)\)
is bounded above in total by \(\mathcal F_r(A_n)\).
The localized form of (AS2) with a nonnegative compact
cutoff gives that every locally weak limit of \(\omega_n\)
dominates \(\omega_T\). This follows by applying the
calibration identity multiplied by that cutoff and passing
to the bulk limit for each test.
The total limiting upper bound is
\(\mathcal F_r(T)=\omega_T(\mathbb R^2)\), so the only
locally weak limit is \(\omega_T\).

Each of the three cyclic pair differences \(V_i-V_j\)
is nonzero, and distinct differences are nonparallel
at every finite point. To see the latter property let
\(u_1,u_2,u_3\) be the cyclic reference normals and
\(H_\ell=T_r\nu_\ell(T)>0\).
Apart from the common positive normalization the
three cyclic differences are
\[
          F_i=2H_i u_i-\sum_{j\ne i}H_j u_j.
\]
They sum to zero. In the basis \(u_1,u_2\),
the determinant of \(F_1,F_2\) is
\(3(H_1H_2+H_2H_3+H_3H_1)>0\), since \(u_3=-u_1-u_2\).
Thus the nonparallelity holds for every pair of differences.

Every limit of the three tagged interface measures is
supported on the simplex rays, because \(\omega_T\)
is supported there and the weights are positive.
On a compact interior patch of a reference ray for
pair \(e=(i,j)\), take the matrix field in (AS2) to
be the rank-one contraction sending the normalized
\(V_i-V_j\) into the fixed reference normal, multiplied
by a nonnegative compact cutoff supported off the
other rays. The limiting volume test saturates
the weighted surface mass on that patch at \(T\).
Therefore, in the approximating sequence, the integral
of the nonnegative difference between the weighted
mass and this calibration tends to zero. For either
other pair \(f\), that difference is bounded below by
its surface weight minus the absolute value of its
projection onto the normalized \(V_i-V_j\). This
quantity has a strict positive lower bound on every
smaller compact patch by the nonparallelity. Thus
the other pair measures tend to zero there. Equality of
the limiting weighted mass then identifies the surviving
measure with \(\nu_e(T)\). There is no mass at the
junction, because \(\omega_T\) has no point mass and
the weights stay positive there. This proves the
local tagged convergence.

The same hypotheses also control tails in the smoothed
surface norm. For fixed \(R>0\), choose the inner directions
on \(B_R\) to pair positively with \(G_*\), so their
smoothed norm is at least \(\omega_n(B_R)\).
Outside \(B_R\) choose any admissible directions.
With both signs of those outer directions the total
squared norm is bounded by \(\mathcal F_r(A_n)^2\).
Averaging the two squared norms cancels the cross term.
Thus the outer squared norm is at most
\(\mathcal F_r(A_n)^2-\omega_n(B_R)^2\).
Using equilateral directions on the tail gives
\[
 \lim_{R\to\infty}\limsup_n
      \|T_{\sqrt r}(\nu_{\rm sum}(A_n)|_{B_R^c})\|_2=0.
                                                               \tag{AS4}
\]
Here \(\omega_n(B_R)\to\omega_T(B_R)\) because the
reference surface has zero mass on the bounding circle.
The smoothed norm bounds the total mass of the positive
tail measure, so all three Gaussian pair perimeters
converge globally as well.

In addition any sequence of energy-attaining directions
at \(A_n\) has a subsequence whose smoothed tensors
converge strongly to a rotated or reflected simplex
normal tensor. In fact the tagged convergence and
bounded directions give locally weak limits
\(u_e\,\nu_e(T)\), with \(|u_e|\le1\).
On each compact disk their smoothed tensors converge
strongly, by the continuous bounded product kernels,
and (AS4) makes the error outside the disk uniformly
small. The limiting norm equals \(\mathcal F_r(T)\).
Equality in the first extra-coordinate form at \(T\)
forces each coordinate of the direction field to be
a translation speed and forces \(|u_e|=1\) on the
reference rays. Hence \(u_e=B N_e^T\) for one planar
matrix \(B\). A symmetric quadratic form that agrees
with the Euclidean form on three distinct unoriented
unit directions has zero difference from it. Applying
this to \(B^{\mathsf T}B\) shows \(B\) is orthogonal.
This gives the stated strong compactness of attaining
tensors at a simplex limit.
