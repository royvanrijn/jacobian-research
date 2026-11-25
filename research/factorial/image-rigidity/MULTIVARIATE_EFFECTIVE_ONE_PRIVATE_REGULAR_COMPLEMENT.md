# Effective regular comparison for one private side

Use \(P,p,\sigma,\xi_j\) of the
[uniform one-private regular proof](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md).
Thus \(P\) is total private degree, \(p\) the number of private rows,
and the exact field drift is
\(\xi_N-\xi_n=\Delta\log(N/n)\).

**Theorem.** Given a rational upper total curvature cap \(K\)
and a rational \(\kappa>0\), a terminating finite search gives a
uniform consecutive cutoff for profiles with exactly one private side,
whenever \(P\ge2\) or \(P=1\) and \(\sigma\ge1+\kappa\).
The cutoff allows all physical fields and integer degree differences.
Together with the effective near-critical one-private theorem and
the absolute high-curvature bound, this gives a computable absolute
cutoff for every ordinary profile with exactly one private side.
No cutoff is evaluated here.

## Compact data

There are finitely many private allocations under the cap \(P\le K\).
The [bounded-private annulus](MULTIVARIATE_BOUNDED_PRIVATE_ROOT_ANNULUS.md)
gives a computable real cap on the normalized log radius using
\(V_{\mathrm{sh}}\le K(K+2)\). The endpoint centering error in
the regular proof is bounded by a computable multiple of
\((V_{\mathrm{sh}}+P\log j)/j\) using the explicit
[endpoint factorial estimate](MULTIVARIATE_GROWING_CURVATURE_ENDPOINT_EXPANSION.md#comparison-of-the-orders).
It differs from the field center by at most
\(P+V_{\mathrm{sh}}\) at leading order. For the shared part,
the bound follows by integrating
\(\sum_i d_i(\log(a_i+d_i q)-\log a_i)\) along the real interval,
whose derivative is at most \(V_{\mathrm{sh}}\).
For the private part it follows from \(F(e_i)=e_i\log e_i-e_i\).
Thus a rational \(M\) can be chosen so all possible root fields
satisfy \(|\Re \xi_j|\le M\).

The shared retained and diffuse nets approximate the shared
curvature on fixed compact contour neighborhoods. Integrate
it from zero with the boundary normalization in (3) of the
regular proof to approximate the phase and its derivatives.
The estimates are effective by the
[stationary fiber nets](MULTIVARIATE_EFFECTIVE_STATIONARY_FIBER_LOCALIZATION.md)
and the same integrated error bounds used in the
[one-private interior search](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_INTERIOR.md#compact-phase-family-and-root-bearing-fields).
The private allocation remains fixed. In particular
\(\sigma=B'(0)\) is enclosed effectively.
The limiting compact family includes the purely diffuse and
zero shared parts where the chosen band allows them.

## Interior and negative-boundary witnesses

At a fixed limiting height \(0<\theta=\Im\xi\le\pi\), keep
a small height margin from zero. The
[actual left-integral bounds](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md#compact-phases-and-the-left-integral)
give the nonzero leading endpoint coefficient \(L=1\) for
\(P\ge2\) and \(L=(1-e^\xi)^{-1}\) for \(P=1\).
The relative endpoint error is \(O(1/j)\) on a sufficiently
small field disk. The constants are computable on such a disk.
For \(P\ge2\), choose a rational \(c\in(0,1/\max e_i)\)
on the pinched left line. After the (finite) shift in the
Barnes coordinate from \(-c\) to \(5/2\), the explicit
exponential rate is
\(-\pi P|v|/2+(\pi-\theta)v\).
It gives a uniform integrable exponential majorant on a disk
with positive margin from \(0,2\pi\). The finite private
Gamma factors and the falling Gamma factor bounds give
computable polynomial and constant factors. The residues
at zero, one, and two then have the powers of \(j\) used in
(4) of the regular proof. For \(P=1\) in the stated band
the Euler-product argument for (5) uses \(\theta_\sigma=0\).
Its real-argument comparison and logarithmic error have
computed bounds from \(V_{\mathrm{sh}}\), the integer row
ratios and the height margin. On the pinched line use the
local \(O((1+|v|)^2/j)\) Gamma-ratio bound up to a computable
multiple of \(\log j\) and integrate it against the same
positive exponential decay. Outside this range the
exponential tail is \(O(1/j)\) by increasing the multiple.
These estimates also truncate the original left line at
any sufficiently large fixed \(|\Im q|\), with a computable
exponentially small error. They remain valid on disks
through the negative real parameter.

On the right keep the finite contour of the
[actual right-saddle proof](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md#the-actual-right-saddle-on-compact-field-disks).
It consists of a right crossing, horizontal panels joining
the pinched line, and a local simple saddle segment, with
strict loss on the exterior pieces relative to both the
saddle and the endpoint. At the negative boundary use the
real-gap saddle and its vertical quadratic direction.
Use the two-half-plane domain through the right gap as
in the interior search. The exact right contour is an
open path in this simply connected slit domain with its
two endpoints fixed on the pinched line. A path constructed
at the limiting profile has rational polygon approximation
on its compact strict exterior pieces. The endpoints on
the pinched line move by \(O(1/j)\) with simple short
connectors in the same strict sublevels. The saddle may
move within a Rouché disk with a common rational descent
direction; the Hessian loss and the connectors are
verified by the interval tests of the interior search.
Analyticity of the exact integer-order integrand identifies
these replacement paths with the original one. The
nonzero crossing coefficient at each such field follows
from the regular proof and its continuation from the
right-gap anchor. Only the finite contour and strict
inequalities need be checked.

On these compact paths the integer private exponents
stay in a finite list. The shared Stirling remainders
sum to a computable bound by the \(V_{\rm sh}\) bound.
The logarithmic amplitude and derivatives are uniformly
bounded by the compact path distances to the Gamma cuts
and \(\sum|d_i|/a_i\le V_{\mathrm{sh}}\), as in the
[shared regular expansion](MULTIVARIATE_EFFECTIVE_SHARED_REGULAR_COVER.md#uniform-expansion-and-the-common-root-certificate).
The real right gap uses the ordinary positive-argument
Gamma formula. The Gaussian integral and strict exterior
losses then give the actual expansion
\[
       L(\xi)(1+O(1/j))
          +j^{p/2} A(\xi)e^{jD(\xi)}(1+O(1/j)),
\]
up to the pure \(a\) coefficient. The error bounds are
computable, and \(L,A,D\) are the same functions at the
two orders. Here \(D\) is the actual stationary value
relative to the endpoint, and \(D'=q_s\) with the
continued saddle coordinate \(q_s\). Normalization keeps
the power \(j^{p/2}\) from the \(p\) private Gamma factors.

## Positive boundary witnesses

At height zero the [real endpoint and maximum proof](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md#the-positive-real-boundary-and-the-adjacent-comparison)
gives a unique coexistence between the endpoint and a
strict real interior maximum in the stated regimes.
Otherwise one of them strictly dominates. At a stable
endpoint, retain the initial interval \(k\le\eta j\).
For \(P\ge2\) the coefficient bound
\(C^k(k!)^{P-1}/j^{(P-1)k}\) has successive ratios
less than \(1/2\) on that interval for computably small
\(\eta\) and large \(j\). Its first omitted-from-leading
term is \(O(j^{-1})\). For \(P=1\), stability implies a
strict negative real field, including at coexistence.
Choose \(\eta\) and a field disk so the geometric
majorant in the real endpoint proof has ratio below
one. The initial \(O(\log j)\) terms give
\((1-e^\xi)^{-1}+O(j^{-1})\) by the linked factorial
ratio bound, with the remaining geometric tail of the
same order. Both bounds extend to their nearby complex
fields by modulus estimates.

For a contributing interior maximum, choose a short
real interval, check the nonzero Hessian and a strict
phase loss at the endpoints, and apply the
[positive-boundary lattice expansion](MULTIVARIATE_EFFECTIVE_SHARED_REAL_BOUNDARIES.md#upper-boundary).
Only the relative \(O(1/j)\) bound is needed here, with the
private powers in the exact summand giving the stated
\(j^{p/2}\) after integration. The complement of the
retained maximum and endpoint intervals has strict real
phase loss, which is a finite interval test with the
continuity of \(q\log q\) and the compact phase bounds.
When the endpoint is strictly dominated, the coarse
entropy and factorial bounds on that entire complementary
interval suffice. When the maximum is strictly dominated,
the corresponding strict loss suffices there. Thus the
positive boundary has actual one-term witnesses away
from coexistence and an actual two-term witness at it.

## Linked comparison and termination

At interior height and at the negative boundary,
a two-term modulus balance has at most one field at a
given height by \(\Re D'=\Re q_s>0\) in the regular proof.
At the positive boundary use the unique coexistence.
Consequently for each limit of pairs of fields at equal
heights, either one field has a strict dominance witness
or the two fields agree at a balance point. At that
point the real slope is strictly positive, including
the positive-boundary saddle in \((0,1)\). Use one
convex field chart with \(\Re D'\ge c>0\) for both
fields. Away from balance the strict dominance test
absorbs the power \(j^{p/2}\) after a computable terminal.

At an actual root in the common chart the expansion gives
\[
 j\Re D(\xi_j)=-(p/2)\log j+C(\xi_j)+O(j^{-1}),
 \qquad C=\log|L/A|.
\]
A computable derivative bound \(C_1\) for \(C\) on the
common chart and the logarithmic modulus equations give
\[
       (Nc-C_1)|\xi_N-\xi_n|\le C_2(1+\log n)/n .
\]
This excludes nonzero integer degree difference, since
the exact displacement is then at least \(1/N\).
For \(\Delta=0\) subtracting at identical fields gives
\(\Re D=O(1/n)\), inconsistent with
\(n\Re D=-(p/2)\log n+O(1)\). All thresholds in these
tests use the computed constants on the witness.

The profile bands and the paired closed field rectangle
form a compact family. Enumerate phase enclosures,
rational contour and real-lattice witnesses, and their
strict margins in parallel with refinement. Each limit
of field pairs admits one of the open witnesses just
proved. The retained and diffuse nets and strict
endpoint or contour margins therefore give a finite
subcover and a largest certified terminal. Finally use
the [all-field near-critical exclusion](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_INTERIOR.md#consequence-for-the-whole-subcritical-and-near-critical-margin)
on \(P=1,\ 0\le\sigma\le1+\delta\), the regular theorem
with \(\kappa=\delta/2\), and the
[absolute high-curvature proof](MULTIVARIATE_EFFECTIVE_ABSOLUTE_CURVATURE_BOUND.md)
with \(K=64\).
