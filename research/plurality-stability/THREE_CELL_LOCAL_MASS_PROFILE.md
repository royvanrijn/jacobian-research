# The local mass profile at the Gaussian tripod

Let \(T\) be the regular centered three-sector partition of the standard
Gaussian plane, and write \(m_*=(1/3,1/3,1/3)\). At every
\(0<r<1\), the largest stability at prescribed masses \(m_*+h\)
among partitions close to the orbit of \(T\) has an exact quadratic
expansion in \(h\). Its coefficient is obtained by minimizing the
simplex interface form over a prescribed mass flux. The result applies
to arbitrary measurable competitors. Translation alone yields a
strictly smaller coefficient.

The expansion separates the tripod branch from other maximizing
partitions when the latter have a nonzero price or a sufficiently large
mass curvature. Near a first balanced tie, this comparison forces
coexistence of separated optimizers at unbalanced quotas. It leaves
open the exclusion of the original balanced tie.

## Exact expansion over measurable competitors

Use the facet convention of
[continuation at positive correlation](POSITIVE_CORRELATION_CONTINUATION.md#the-full-interface-form-at-a-fixed-positive-correlation).
For antisymmetric facet speeds \(f=(f_{ij})\), \(f_{ij}>0\) moves the
interface out of \(T_i\). Put \(d\nu_{ij}=\phi_2\,d\mathcal H^1|_{F_{ij}}\)
and \(\lambda_i^f=\sum_{j\ne i} f_{ij}\nu_{ij}\), a signed measure on
the incident rays. Write
\[
 Cf=(\lambda_i^f(\mathbb R^2))_i,\qquad
 Q_r(f)=
  \sum_{i<j}\int a_{r,ij} f_{ij}^2d\nu_{ij}
   -\sum_i\iint M_r(x,y)\,d\lambda_i^f(x)d\lambda_i^f(y).
\]
Here \(M_r\) is the Mehler kernel and \(a_{r,ij}>0\) is the normal
derivative of the \(i,j\) score contrast on the facet. In the
notation of (B6), \(Q_r=r\mathfrak q_r\) on balanced speeds.
The definition above also allows \(Cf\ne0\).

Fix a compact interval \(J\subset(0,1)\). There is \(\eta>0\) such that
for \(r\in J\), the supremum \(V_r^T(m_*+h)\) over ordinary partitions
of these masses at Gaussian \(L^1\) distance at most \(\eta\) from the
tripod orbit satisfies
\[
 V_r^T(m_*+h)-S_r(T)
     =\Lambda_r|h|^2+o_J(|h|^2),
 \qquad
 \Lambda_r|h|^2=-\min_{Cf=h}Q_r(f).
 \tag{LM1}
\]
The mass plane is equipped with its Euclidean norm as a subspace
of \(\mathbb R^3\). The minimum is attained modulo planar rotations,
and \(\Lambda_r\) is continuous in \(r\).

Here is the sequential expansion underlying the uniform statement.
Let \(r_j\to r\in J\), \(t_j>0\) tend to zero, and
\(h_j=t_j v_j\) with \(v_j\to v\) in the mass plane. Then
\[
 \lim_j \frac{V_{r_j}^T(m_*+h_j)-S_{r_j}(T)}{t_j^2}
              =-\min_{Cf=v} Q_r(f).
 \tag{LM2}
\]
The balanced gap (S37) and (B11) hold uniformly on \(J\) after a
finite covering of \(J\). The two normalizing factors in (B6) are
bounded above and below uniformly on \(J\).
[Three-cell spectral gap](SIMPLEX_INTERFACE_SYMMETRY.md#reflected-kernels-close-the-three-cell-spectral-problem).
[Uniform measurable comparison](POSITIVE_CORRELATION_CONTINUATION.md#an-interface-gap-implies-a-uniform-measurable-local-comparison).

For the upper bound take ordinary \(A_j\) within \(o(t_j^2)\) of this
supremum. A small translation of \(T\) to the required masses gives
\(S_{r_j}(A_j)\ge S_{r_j}(T)-O(t_j^2)\), since the translation
mass derivative is invertible and the unpriced first variation at
\(T\) vanishes. First align \(A_j\) to a closest simplex in Gaussian
\(L^1\). Two compact smooth vector fields with independent mass
flux at \(T\) still have an invertible mass derivative for every
sufficiently close ordinary partition. Their joint flow repairs
the masses of \(A_j\) to balance using parameters \(O(|h_j|)\).
Write \(B_j\) for the corrected partition. The change of variables
in the joint Gaussian integral gives
\[
 |S_{r_j}(B_j)-S_{r_j}(A_j)|
        \le C |h_j|\,d_3(A_j)+C|h_j|^2.
 \tag{LM3}
\]
Indeed the first flow derivatives at \(T\) are zero. Their change
from \(T\) to \(A_j\) is \(O(d_3(A_j))\) by the volume formula (VV3)
before its stationary-price substitution. On the fixed supports of the
fields, the smoothed scores and their first derivatives change by
\(O(d_3(A_j))\), by boundedness of the corresponding Mehler rows as
functions of the source point. The flow second derivatives are
uniformly bounded by the joint density formula.
[Mass transports and the volume formula](THREE_CELL_VOLUME_VALUE_DERIVATIVES.md#ambient-transport-and-uniform-mass-correction).

The flow and its inverse distort Gaussian \(L^1\) distances by bounded
factors. They change any regular three-sector partition by
\(O(|h_j|)\), since the changed region lies within \(O(|h_j|)\)
of its rays on a fixed compact set. Thus
\(d_3(B_j)\le C d_3(A_j)+C|h_j|\) and
\(d_3(A_j)\le C d_3(B_j)+C|h_j|\). Take \(\eta\) small enough that
\(B_j\) lies in the (B11) neighborhood. By (LM3) and (B11),
\[
 c\,d_3(B_j)^2\le C t_j^2+C|h_j|\,d_3(B_j),
\]
so \(d_3(A_j)=O(t_j)\).

Now use the nearest-moment orientation from (B11) and align it to
\(T\). This orientation preserves convergence to the tripod and
makes the moment defect orthogonal to the rotation moment tangent
space. Put \(L_j=L_{r_j}(A_j;T)\), the pointwise score regret, and
\(e_{j,i}=1_{A_{j,i}}-1_{T_i}\). With no balance assumption the exact
stability expansion is
\[
 S_{r_j}(A_j)-S_{r_j}(T)
   =\sum_i\langle e_{j,i},T_{r_j}e_{j,i}\rangle-2L_j.
 \tag{LM4}
\]
In particular the degree-zero term of the interaction is
\(|h_j|^2\).

We first show \(L_j=O(t_j^2)\). Otherwise take a subsequence where
\(L_j/t_j^2\to\infty\), and put \(s_j=\sqrt{L_j/r_j}\). The normal
flux compactness proof of (B13)–(B14) applies to \(e_{j,i}\gamma_2/s_j\):
the saturated margin and (B3) are independent of balance, and the
limiting total flux is zero since \(h_j/s_j\to0\). The moment gauge
passes to the limit by weighted tightness. The resulting speed \(f\)
satisfies
\[
 \sum_{i<j}\int a_{r,ij}f_{ij}^2d\nu_{ij}\le 2r,\qquad
 \lim_j s_j^{-2}\sum_i\langle e_{j,i},T_{r_j}e_{j,i}\rangle
       =\sum_i\iint M_r\,d\lambda_i^f\,d\lambda_i^f .
\]
The mean-zero condition lets one use (B12) on this limit, giving an
upper bound \((1-\epsilon)2r\) for the interaction with some
\(\epsilon>0\). But (LM4) and
\(S_{r_j}(A_j)-S_{r_j}(T)\ge-O(t_j^2)\) give a lower bound \(2r\).

Now normalize the same flux by \(t_j\). The bound \(L_j=O(t_j^2)\)
gives weighted tightness by (B3) and bounded \(L^2\) densities for the
one-sided transfers projected onto their pair planes. After passing
to weak limits on these planes, the signed flux limit is
\((\lambda_i^f)_i\), supported on the actual facets. The supporting
plane is reached because the saturated score margin bounds the
normalized mass a fixed distance from that plane by \(O(t_j)\).
The weak \(L^2\) densities place zero mass on the junction. The limit
satisfies \(Cf=v\). The normal mass-cost estimate in (B13) gives
\[
 \liminf_j 2L_j/t_j^2 \ \ge\
           \sum_{i<j}\int a_{r,ij}f_{ij}^2d\nu_{ij}.
\]
This estimate uses the normal Taylor expansion of the score on the
relative interior of each facet and the fact that the cost of the two
one-sided transfers bounds the squared net flux there. The full
interaction in (LM4) converges after division by \(t_j^2\) to the
Mehler interaction of \(\lambda^f\). On a compact product this follows
by uniform approximation by sums of separated continuous functions.
On the complement, the exponential weight in (B3) controls the kernel
uniformly for \(r_j\in J\). In particular its constant part converges
to \(|v|^2\) and remains in the formula. Together these estimates prove
the upper bound in (LM2).

For the lower bound, take a speed with \(Cf=v\). Approximate it in
the weighted \(L^2\) norm by smooth speeds supported on compact
relative-interior pieces of the facets, using independent mass bumps
to keep the flux exactly \(v\). Change the bump coefficients by \(o(1)\)
to prescribe flux \(v_j\). Place the resulting smooth signed speeds
on normal strips of amplitude \(t_j\), and repair the remaining
\(O(t_j^2)\) mass errors by the same bumps. The unpriced first
variation vanishes at the tripod, so the quadratic coefficient is
\(-Q_r(f)\) for each fixed approximate speed in the limit. This is
the normal-strip expansion of (B9) with the full constant Mehler
term retained. Passing to the speed limit proves the reverse
inequality.

The mass constraint is surjective by translations. The rotation
kernel pairs to zero against every speed by rotation invariance and
unconstrained stationarity of \(T\). On the remaining balanced
subspace the form is coercive, so the minimizer in (LM2) exists
and depends linearly on its mass flux, modulo rotations. Permutation
symmetry of \(T\) makes the optimized quadratic form a scalar on the
mass plane. Its coefficient is continuous in \(r\), by operator-norm
continuity of the form on the common unweighted \(L^2\) facet space
and uniform coercivity on the balanced quotient. The Mehler blocks
have locally uniform Hilbert–Schmidt domination as in (B7). For the
diagonal coefficient, write \(\sigma=\sqrt{1-r^2}\), and put the
outward ray in coordinates where the competing sectors are
\(n\ge0,\ s\ge-n/\sqrt3\) and its \(n\)-reflection. Differentiating the
Gaussian contrast normally at \((0,s)\) gives
\[
 a_r(s)=\frac{2r}{\sigma}
       \int_0^\infty z\phi_1(z)
              \Phi\left(\frac{rs}{\sigma}+\frac z{\sqrt3}\right)\,dz.
\]
It is continuous in \(r\) in the \(L^\infty\) norm on the ray on
every compact \(J\subset(0,1)\). The derivative in \(r\) of the
normal probability is bounded by \(Cs\phi_1(c s)\) for \(s,z\ge0\)
with \(C,c>0\) depending on \(J\). For the minimum in (LM2), one may
write the speed as a fixed translation realizing \(h\) plus a
balanced speed and apply the uniform gap to bound the minimizing
balanced speed. Compactness of \(J\) and (LM2) give the uniform
remainder in (LM1).

The same expansion holds for the cylindrical tripod in every fixed
finite dimension \(d\ge2\). Use the full interface gap (S42) and
the same flux and compact transport proof in that dimension, with
two independent correcting fields supported on regular facet
pieces. The margin, pair-plane projections and weighted kernel
bounds (B2)–(B7) apply to the cylindrical facets in \(\mathbb R^d\).
For the nearest-moment gauge, a regular moment matrix specifies its
orthonormal two-frame by applying a fixed left inverse to the planar
moment matrix. Nearby moment matrices specify nearby frames, and
those can be aligned by orthogonal maps near the identity by
completing the frames locally. The moment tangents of the initial
plane and the rotations mixing it with the unused space give the
full rotation space in (S42). This verifies the same gauge and
coercivity argument of (B11) in the fixed ambient dimension.
Expanding the full interface form in the unused-coordinate
Hermites then gives the same minimum as in the active plane.
Only the zeroth block contributes to the volume flux, and all
other blocks are nonnegative by (S41)–(S43). The planar minimizing
speed can be extended constantly in the unused directions. The
coefficient is therefore independent of dimension. This extension
makes no assertion of dimension-uniform neighborhood or remainder.
[Higher-dimensional interface lift](SIMPLEX_INTERFACE_SYMMETRY.md#every-cell-count-the-active-gap-lifts-to-all-ambient-dimensions).

## The correction to translation curvature

Use the weighted norm
\(\|f\|_a^2=\sum_{i<j}\int a_{r,ij}f_{ij}^2d\nu_{ij}\).
Let \(W\) be the space of balanced speeds orthogonal to rotations
in this norm, and let \(\mathsf L_W\) be the positive definite
operator on \(W\) representing the restriction of \(Q_r\). For a
mass change \(h\), let \(s_h\) be the translation speed with
\(Cs_h=h\), and put \(q=(1-r)/r\). The translation row identity
gives \(Q_r(s_h,g)=-q\langle s_h,g\rangle_a\) for all \(g\).
Minimizing \(Q_r(s_h+g)\) on \(W\) yields the correction formula
\[
 \Lambda_r|h|^2
    =q\|s_h\|_a^2+
            q^2\langle P_Ws_h,\mathsf L_W^{-1}P_Ws_h\rangle_a .
 \tag{LM5}
\]
[Translation row identity](SIMPLEX_INTERFACE_SYMMETRY.md#the-standard-sector-and-an-explicit-full-gap).

On every ray, write \(a=a_r(s)\) with \(s\) arclength from the
junction and \(\nu\) its Gaussian surface measure. In the formulas
below all integrals are on one ray. Translation speeds are the
constant-per-ray speeds of the standard representation. For a
mass vector \(w\) of sum zero, the weighted mass adjoint is
\((C^*w)_{ij}=(w_i-w_j)/a\), and its projection to translations
has the same differences times \(\int1\,d\nu/\int a\,d\nu\).
Thus \(CC^*=3(\int a^{-1}d\nu)I\) on the mass plane, while the
restriction \(C_{\rm tr}\) of \(C\) to translations gives
\(C_{\rm tr}C_{\rm tr}^*=3(\int1\,d\nu)^2/(\int a\,d\nu)I\).
Translations are orthogonal to rotations, so projection onto
\(W\) agrees with projection onto \(\ker C\) on these speeds.
It follows that
\[
 \|P_W s_h\|_a^2=\Theta_r\|s_h\|_a^2,\qquad
 \Theta_r=1-\frac{(\int1\,d\nu)^2}
                     {(\int a\,d\nu)(\int a^{-1}d\nu)}.
 \tag{LM6}
\]
In particular \(\Theta_r>0\). The normal derivative formula above
shows that \(a\) is strictly increasing along the ray. Its endpoint
values are \(r\ell/\sqrt{1-r^2}\) and \(2r\phi_1(0)/\sqrt{1-r^2}\),
where \(\ell=3\phi_1(0)/2\) is the first-moment difference of
the tripod cells, by polar integration over the three sectors.
Thus the lower endpoint is \(3/4\) of the upper.
If \(0<l\le a\le u\), the pointwise inequality \(a+lu/a\le l+u\)
and the arithmetic-geometric mean inequality give
\((\int a)(\int a^{-1})/(\int1)^2\le(l+u)^2/(4lu)\).
For the stated endpoint ratio this implies \(\Theta_r\le1/49\).

Set
\[
 \Lambda_r^{\rm tr}
       =(1-r)\frac{\dot S_r(T)}{W_1(T)},\qquad
 W_1(A)=\sum_i\left|\int_{A_i}x\,d\gamma_2(x)\right|^2.
\]
The bulk translation identity (AM3) gives
\(q\|s_h\|_a^2=\Lambda_r^{\rm tr}|h|^2\), since the correlation-gradient
matrix and the translation centroid Gram matrix of the tripod are
scalar with traces \(\dot S_r(T)\) and \(W_1(T)\), respectively.
[Translation Hessian](GAUSSIAN_ANISOTROPIC_MAXIMIZER_DIMENSION.md#balanced-translations-with-a-matrix-correlation).
Both \(P_Ws_h\) and its inverse image under \(\mathsf L_W\) are
in the standard sector by permutation equivariance. For \(r>1/3\),
(S34) and \(a\le2r\phi_1(0)/\sqrt{1-r^2}\) bound the weighted
standard gap below by \(q(3r-1)/4\). Inserting this bound in
(LM5) gives
\[
 \Lambda_r^{\rm tr}<\Lambda_r
      \le \Lambda_r^{\rm tr}
              \left(1+\frac{4\Theta_r}{3r-1}\right)
      \le \Lambda_r^{\rm tr}
              \left(1+\frac{4}{49(3r-1)}\right).
 \tag{LM7}
\]
The strict lower inequality holds for all \(0<r<1\) by
(LM5)–(LM6) and the all-correlation gap.

## The common leading price near a unique balanced simplex

If the simplex is the unique balanced maximizing orbit at \(r\),
strong maximizing compactness puts all nearby-mass global
maximizers in the orbit neighborhood of (LM1). Thus (LM1)
also expands the global optimal value at that \(r\). In
addition, every maximizer \(A_h\) at \(m_*+h\) has
\[
                 \beta(A_h)=\Lambda_r h+o(|h|).              \tag{LM9}
\]
Indeed compare the value at \(m_*+h\) with that at
\(m_*+h\pm\theta |h|e\), where \(e\) is a unit mass direction
and \(0<\theta<1\). The strong convex price supports (VV8)
give upper and lower bounds for
\((\beta(A_h)-\Lambda_r h)\cdot e\). By (LM1) these bounds
are in magnitude at most \(C\theta |h|+o(|h|)/\theta\),
uniformly in the direction and the maximizer. First take
\(|h|\) to zero for fixed \(\theta\), then let \(\theta\)
decrease to zero. This proves (LM9).
[Global price supports](THREE_CELL_VOLUME_VALUE_DERIVATIVES.md#strong-convexity-across-positive-masses).

## Consequences at a first balanced tie

Suppose at \(r_*\) the tripod and a nonsimplex \(A\) are both
balanced global maximizers, and the balanced tripod is the unique
maximizing orbit immediately below \(r_*\). Put
\[
 d_A=\dot S_{r_*}(A)-\dot S_{r_*}(T)\ \ge0 ,
\]
where the sign follows from simplex optimality to the left.
Let \(\beta(A)\) be the centered prices. Either of the following
conditions forces unbalanced masses tending to \(m_*\) at
correlations increasing to \(r_*\) from below at which there
are two global maximizers uniformly separated in \(L^1\) modulo
rotations and relabelings.

1. The price vector \(\beta(A)\) is nonzero.
2. The pressure-free hypothesis \(\beta(A)=0\) holds and the
   largest translation curvature of \(A\), as a quadratic form
   in its volume changes, exceeds \(\Lambda_{r_*}\).

The second condition has, when \(r_*>1/3\), the sufficient criterion
\[
 \frac{W_1(A)}{W_1(T)}
     <\left(1+\frac{4\Theta_{r_*}}{3r_*-1}\right)^{-1}.
 \tag{LM8}
\]
In fact the largest eigenvalue of the translation curvature of
\(A\) is at least \((1-r_*)\dot S_{r_*}(A)/W_1(A)\). To see
this, compare the quadratic forms \(M_A(v,v)\) and
\(|C_Av|^2\), where \(C_A\) is the derivative of the translation
mass map. Their traces are \(\dot S_{r_*}(A)\) and \(W_1(A)\).
The largest generalized eigenvalue bounds the ratio of these
traces. Equations (LM7) and \(d_A\ge0\) now give (LM8).
[Exact mass translation lower bound](THREE_CELL_VOLUME_VALUE_DERIVATIVES.md#strong-convexity-across-positive-masses).

To prove the coexistence claim, write \(r=r_*-\epsilon\). In the
first case take \(h=C\epsilon v\) with \(\beta(A)\cdot v>0\)
and \(2C\beta(A)\cdot v>d_A\). Translate \(A\) to masses \(m_*+h\).
By the smooth joint volume and stability expansions in the translation
and correlation parameters, it then beats the value of (LM1) for
small \(\epsilon>0\). In the second case take
\(h=C\sqrt{\epsilon}\,v\), where \(v\) is a unit direction of
strictly larger translation curvature. If the excess curvature
is \(\kappa>0\), choose \(C^2\kappa>d_A\). The translated
partition beats (LM1) by \((C^2\kappa-d_A)\epsilon+o(\epsilon)\).
In either case the comparison is between actual ordinary partitions
and the supremum of the full orbit-local class. The Taylor
estimates use positive rank of \(C_A\) and smooth Gaussian
density integrals for the translated fixed partition.

By balanced compactness and (B11), at \((r_*,m_*)\) there are no
global maximizers in some closed annulus between the tripod
orbit and its distant maximizing competitors. Strong compactness
with varying positive correlation and masses preserves this excluded
annulus for all sufficiently close parameter pairs. Indeed the limit
of any sequence of optimizers is ordinary and maximizing, by weak
compactness, strict positivity and the exact mass repairs of (VV1).
Choose the inner radius within the local neighborhood in (LM1).
For each small \(\epsilon>0\), vary the masses from \(m_*\) to
\(m_*+h\) on the line segment at the fixed \(r=r_*-\epsilon\).
At the balanced endpoint all optimizers are in the inner component.
At \(m_*+h\) all are outside the annulus by the strict comparison
with (LM1). The sets of segment parameters admitting an inner and
an outer maximizer are both closed, by strong maximizing
compactness, and their union is the whole segment. Hence they
intersect at an interior parameter. The annulus gives the same
strict separation for every sufficiently small \(\epsilon\).

If all nonsimplex ties at \(r_*\) satisfy \(d_A>0\), then in the
first case the phase coexistence also gives a pressure jump along
\(v\). The set of nonsimplex tied shapes is compact modulo
rotation and relabeling, separated from \(T\) by (B11). Its
minimum derivative gap is therefore \(d_0>0\). Write the coexistence
quota as \(m_*+s_\epsilon v\), with \(0<s_\epsilon<C\epsilon\),
and take a sequence of outer maximizers \(B_\epsilon\) converging
to a nonsimplex tied \(B\) as \(\epsilon\downarrow0\).
Its value at coexistence is \(S_r(T)+O(s_\epsilon^2)\)
by (LM1) for the inner maximizer. Transport \(B_\epsilon\)
to balance at \(r\) and increase the correlation to \(r_*\).
The transported value at \(r_*\) is at most \(S_{r_*}(T)\).
The transport slope (VV3) and uniform correlation derivatives yield
\[
 2s_\epsilon\,\beta(B_\epsilon)\cdot v
        \ge \epsilon d_B+o(\epsilon),
\]
because \(s_\epsilon=O(\epsilon)\). Consequently
\(\liminf 2\beta(B_\epsilon)\cdot v\ge d_0/C\);
the compactness argument applies to every limit subsequence.
The inner maximizers have \(O(s_\epsilon)\) prices by the local
distance estimate and (VV5). Their distinct slopes make the
optimal mass profile nondifferentiable in direction \(v\) at
these unbalanced quotas.
