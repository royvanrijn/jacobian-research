# Terminal-thirteen separation with retained shared rows

For distinct nonnegative integer exponent columns, omit the equal rows.
Write \(\ell,d\) for the smaller exponent and positive gap of a shared
unequal row. Its strength is \(s=d^2/(\ell+d)\), and \(S_0\) is the sum
of these strengths. The two multisets of private gaps are given in
increasing and decreasing orientation. Their total private degree is
the sum of their gaps.

**Theorem.** The factorial moment polynomials at orders twelve and
thirteen have disjoint root arguments under independent positive
rescalings if there is a private row, or if all unequal rows are
shared and some reciprocal \(d/(12\ell+1)\) exceeds \(1/8\).
The two polynomials are therefore coprime over every
characteristic-zero field.

## Retained lists

At orders \(n=12\) and \(N=13\), use the strength cutoffs of the
[fixed-ray proof](MULTIVARIATE_TERMINAL_THIRTEEN_HIGH_CURVATURE.md).
Retain each shared row with \(d/(12\ell+1)>c=1/8\) together with its
orientation. Also retain the private rows in full. The fixed-ray
proof covers the strength at and above each cutoff and every
private degree at least four. Below degree four there is one
private allocation at degree one, three at degree two, and five
at degree three, up to column reversal. The checker uses them all,
including the separate fully shared family. The retained
enumeration uses the closed strength budget, so the cutoff can occur
in both proofs.

For a row above \(c\), the bound \(\ell<(d/c-1)/12\) gives
\(d<B(1+1/(12c))\) in a strength budget \(B\). The integer candidate
rows and their multiplicities are consequently finite. The checker
enumerates their ordered multisets and rejects the empty list only
in the fully shared family, where a retained row is required. When
the two private gap multisets are equal, column reversal preserves
the private case. In that case compare the sorted tuple of retained
rows with the sorted tuple of reversed rows and keep the smaller.
This includes exactly one representative of every reflection orbit.
For unequal private multisets the private orientations are fixed
and the shared orientations stay unrestricted.

## Residual bounds

After retaining the rows, write \(B\) for the remaining strength
budget. Let \(F_{j,k}\) denote the residual factorial coefficient
divided by its endpoint chord and without the binomial coefficient.
Set \(v=(F_{12,1}F_{12,11})^{1/2}\). The Gauss product bounds each
residual row's contribution to \(-\log v\) by
\((d/2)\log(1+11x)\) for \(x=d/(12\ell+1)\). Since
\(d/s<1+1/(12x)\), the monotonicity of
\((1+1/(12x))\log(1+11x)\) gives
\[
 -\log v \le \frac56\log(19/8)B.
\]
For the monotonicity substitute \(z=11x\) and differentiate. The
sign is that of \(z(1+12z/11)-(1+z)\log(1+z)\), which increases
from zero because its derivative is \(24z/11-\log(1+z)>0\).
The residual mean lies between an outward rounding of
\((8/19)^{5B/6}\) and \(1\). On a cell with lower endpoint \(v_-\),
the identity \(v^{-2}=\prod_h(1+11x_h)\) bounds each lower reciprocal
by \((v_-^{-2}-1)/11\). The checker also enumerates admissible
individual rows to refine this cap.

Apply the [endpoint attenuation](MULTIVARIATE_ORDINARY_ADJACENT_COEFFICIENT_ATTENUATION.md#an-all-order-endpoint-inequality)
to the residual shared factors with \(M\ge2\).
Its multiplier is at most \((324/325)^{k^2S_{\mathrm{res}}}\) at
distance \(k\) from either endpoint. By \(\log(325/324)\ge1/325\)
and the exact comparison \(325(5/6)\log(19/8)<235\), this gives
\[
 \begin{split}
 F_{13,k}&\le F_{12,k}v^{k^2/235},\\
 F_{13,k}&\le F_{12,k-1}v^{(13-k)^2/235}
 \qquad(1\le k\le12).
 \end{split}
\]
The Gaussian residual factors are \(v^{e_{12,k}}\) and
\(v^{\rho e_{13,k}}\), where \(e_{j,k}=k(j-k)/(j-1)\) and
\(\rho=-\log (F_{13,1}F_{13,12})^{1/2}/(-\log v)\) when \(v<1\).
For a matched reciprocal \(x'=12x/(13-\gamma x)\),
\(0<\gamma\le1\), there is a lower factor \(1+(144/13)x\) for
\(1+12x'\). The resulting logarithmic quotient is nonincreasing
in \(x\), by concavity of \(y\mapsto\log(1+(144/143)(e^y-1))\).
Thus at cap \(1/8\) the lower allowed quotient is at least
\(\log(31/13)/\log(19/8)>1+1/235\), verified by rational logarithm
enclosures. The reference ratio \(144/143\) is allowed by the
[matched outer-loss bounds](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md).
The endpoint inequalities hold for the Gaussian families because
\[
 (1+1/235)e_{13,k}-e_{12,k}-k^2/235
 =k(k-1)(1/132-13/(12\cdot235))\ge0,
\]
and reflection treats the other endpoint. At \(v=1\) both Gaussian
arrays are the all-ones multiplier.

The independent coefficient bounds and the endpoint inequalities
enclose the deformation from the actual residual to the Gaussian
residual. The checker can also use the
[grouped logarithmic moments](MULTIVARIATE_TERMINAL_TEN_SHARED_LIGHT.md#grouped-moment-and-moving-ray-comparison)
with the [matched tanh transport](MULTIVARIATE_MATCHED_TANH_MOMENTS.md).
At fixed \(v,\rho\), scale skew and mean-defect moments jointly
to zero. The bounds on the mean defects, including their lower
bounds in terms of signed moments, hold throughout the scaling
since \(t^{2l}\le t\) for \(0\le t\le1\) and \(l\ge1\).
The endpoint inequalities persist because they are affine in the
logs of the coefficients and hold at both ends. The refined model
also uses the ordinary row upper bound on the loss ratio and the
[excess-loss transport](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md).
Scaling the moments preserves the transport bounds because the
transport errors scale and the outer losses stay fixed. For the
higher signed moments the refined enclosure centers \(y^{2l-2}\)
at \(\lambda=3/4,2/3,5/8\) for \(l=2,3,4\), respectively, with
\(y=T/T_*\). If \(q=2l-2\) and \(r=1-\lambda\), then
\((q\lambda)^q\lambda\le(q+1)^{q+1}r^q\) holds by direct integer-power
comparisons. The endpoints and the unique interior stationary point
therefore give \(|y^{q+1}-\lambda y|\le r\) on \([0,1]\).
Multiplying by \(T_*^{2l-1}T\), summing over the factors with
their signs, and using the sum of the \(T\)'s bounded by the outer
loss at the corresponding order proves the centered bound.
Afterwards move \(\rho\) to \(144/143\), then \(v\) to its cell midpoint, using
Gaussian arrays. Each retained and private factor is held fixed.
The checker bounds its normalized multiplier by a dyadic interval
using integer power comparisons.

## Exact sectors

The two Gaussian midpoint arrays with their retained factors
have isolating balls for their complex roots. The checker requires
their upper root arguments to have the order pattern \(13,12\)
repeated six times, with strict separation between all cosines.
For each adjacent pair it checks a rational separating cosine in
the gap and verifies strict outward Bernstein projection over all
allowed root radii. The normalization and the outer radius check
use the same annulus as the
[terminal-twelve retained proof](MULTIVARIATE_TERMINAL_TWELVE_RETAINED_PROFILES.md#exact-sector-separation).
Floating-point values merely suggest ray positions. The ray tests
and the comparisons of argument cosines use outward enclosures.

The sector counts persist on the deformation. Joining the last
upper sector to its conjugate shows that the last remaining
order-thirteen root lies on the negative real ray by conjugation.
The order-twelve actual moment is positive on the real line
because it is the independent exponential integral of a nonzero
real binomial raised to an even power. Hence no real-ray collision
is possible. Disjoint root arguments imply coprimality, which
holds over every characteristic-zero field by rational Bézout.


## Direct root disks with a small residual loss

On sufficiently narrow cells near \(v=1\) the checker can instead prove
angular separation directly by
[complex balls and Rouché disks](verify_multivariate_terminal_thirteen_retained_disks.py).
Use \(R=144/143\) and the Gaussian references as before. For each
order the multiplier relative to its reference has logarithm
\(\tau_k+E_k\), where \(\tau_k\) comprises the shift in Gaussian outer
loss and the first four signed skew and first four mean-defect
directions. For a fixed cap \(T_*\) on the tanh parameters at this
order, the magnitude of each moment divided by its order \(p\)
is at most \(T_*^{p-1}W/p\), where \(W\) is the outer loss for that
order. This follows from \(\sum_i T_i\le\sum_i\operatorname{artanh}T_i=W\).
The mean-defect moments are nonnegative. The signed bounds are
symmetric. The remaining log series have absolute geometric bounds
using \(1-T_*^2\), with the skew weights and distances described in
the grouped model. We use \(R\le\rho\) from the ordinary outer-loss
inequality and the usual upper logarithmic quotient bound at the
declared cap.

If \(L_k<1\) bounds \(|\tau_k|+|E_k|\), the error of the quartic
exponential approximation to \(e^{\tau_k+E_k}\) is at most
\[
 \frac{L_k^5/5!+|E_k|}{1-L_k}.
\]
The first term is the exponential tail at \(\tau_k\), and the
second follows from the mean-value theorem and
\(e^{L_k}\le 1/(1-L_k)\). Expand the quartic in the scalar moment
parameters while keeping each complete coefficient direction
together. A group with very small parameter bound can be added to
the independent error after bounding its magnitude at every
coefficient. All computations use outward complex or real balls.

Propose centers from the roots of the midpoint references.
At each proposed center \(z_0\), use the Taylor coefficients
\(a_h\) of the moment polynomials around \(z_0\). Each coefficient
vector is translated before multiplying it by its scalar moment
parameter. For a disk of radius \(r\), the strict inequality
\[
 r\inf |a_1|>\sup|a_0|+\sum_{h=2}^j r^h\sup|a_h|
\]
certifies one root throughout the coefficient enclosure by
Rouché's theorem. For consecutive proposed upper centers \(z_0,z_1\)
the checker also verifies
\[
 \operatorname{Im}(\overline{z_0}z_1)>r(|z_0|+|z_1|+r).
\]
Thus every pair of points from those two disks has strictly ordered
arguments. The upper disks and their conjugates are disjoint.
The remaining odd-order disk is centered on the negative axis
and contains an odd-order real root by conjugation and uniqueness.
The disks account for all roots at both orders. The centers are
proposals and the displayed strict inequalities supply the
certificate.


## Replay

Run the [checker](verify_multivariate_terminal_thirteen_retained_profiles.py)
from the repository root with Python 3, python-flint, and a C++17 compiler.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_thirteen_retained_profiles.py --native all
~~~

The fully shared family has 301 lists, with 713 outer and 63,240
radial cells. The one-private linear family has 143 lists.
Its 480 outer cells use 48,597 radial cells. The refined replay of
its empty shared list accounts for 273 outer and 34,324 radial cells.
The other eight private allocations have 37 lists, with 90 outer
and 8,695 radial cells. An independently written multiplicity
enumeration confirms the list counts. For an unoriented
multiplicity assignment \((m_i)_i\), the number of orientation
assignments is \(\prod_i(m_i+1)\). In a reflection-invariant
private case its reflection fixes one such assignment exactly
when all \(m_i\) are even, so the orbit count is
\((\prod_i(m_i+1)+\mathbf{1}_{\mathrm{all\ even}})/2\).
