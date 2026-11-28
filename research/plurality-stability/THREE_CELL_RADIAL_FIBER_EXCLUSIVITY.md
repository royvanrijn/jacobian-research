# Radial overlap of three-cell maximizers

Fix \(0<\rho<1\). In the planar active space of an ordinary Gaussian
noise-stability maximizer with three prescribed positive masses,
every pair of labels has simultaneous positive angular measure on
a positive measure of circles centered at the Gaussian origin.
The overlap has a uniform positive lower bound on compact positive
correlation intervals and compact sets of positive mass triples.
If the radial images of all three strict-winner sets are intervals,
as when those sets are connected, an open interval of circles meets
all three strict winners.

## Circular rearrangement with one other label per radius

For a partition \(A\), let \(p_i(s)\) be the normalized angular measure
of its \(i\)-slice at radius \(s\). Write \(d\nu(s)=s e^{-s^2/2}ds\)
for the Gaussian radial probability measure. Suppose
\(\min(p_1(s),p_2(s))=0\) for \(\nu\)-almost every \(s\). On each circle,
let \(\ell(s)\) be its sole non-3 label, choosing either value if the
circle is full of label 3. Let \(u_s(\theta)\) be the indicator of the
non-3 slice, so the indicator vector on the circle is
\(e_3+u_s(\theta)(e_{\ell(s)}-e_3)\).

Conditional on the two radii \(s,t>0\), the angular pair measure has
density proportional to
\(\exp(\rho st\cos(\theta-\varphi)/(1-\rho^2))\) with respect to the two
uniform angular measures. Both angular marginals are uniform. The
label dot product is
\[
  1-u_s(\theta)-u_t(\varphi)
   +(1+1_{\{\ell(s)=\ell(t)\}})u_s(\theta)u_t(\varphi).
\]
The first three terms depend in expectation only on the slice masses.
The coefficient of the remaining product is positive. The circular
rearrangement theorem bounds its expectation by that of the two
centered arcs of those masses, with the same angular center. The
theorem applies with the two nonnegative indicator functions, one
nonnegative kernel increasing in the cosine, and the second kernel
constant if using its cyclic formulation.
[Friedberg–Luttinger, *Rearrangement Inequality on a Circular Domain*, Theorem I](https://digital.library.unt.edu/ark:/67531/metadc1025197/m2/1/high_res_d/4446751.pdf)

Replace every non-3 slice by its centered arc about one fixed
direction and keep its assigned label \(\ell(s)\). The radius-dependent
arc length is measurable, and the three masses are preserved. Every
pair of radii improves the conditional agreement weakly, so the new
partition has at least the original stability. Each resulting cell is
invariant under reflection through the axis of the common angular
center. Its Gaussian first moment lies on that axis. At an actual
three-cell global maximizer with positive masses the three first
moments span the active plane. Hence the reflected competitor cannot
be a maximizer, and an actual maximizer cannot have two exclusive
radial labels.
[Full active moment rank](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-active-dimension-is-impossible)

## A stability bound from the radial overlap

Let \(V_\rho(m)\) be the unrestricted planar maximum at the positive
mass triple \(m\), and let \(V_\rho^{\mathrm{ref}}(m)\) be the supremum
over partitions whose cells are invariant under reflection through
one fixed axis through the origin. For an arbitrary partition with
masses \(m\), set
\[
      \delta_{12}(A)=\int\min(p_1(s),p_2(s))\,d\nu(s).
\]
Transfer the smaller of the two slices of labels 1 and 2, with a
fixed tie rule, to label 3 at every radius. The changed label mass is
\(\delta_{12}(A)\), and the resulting two labels are radially exclusive.
Apply the preceding circular rearrangement to obtain a reflected
competitor of at least the same stability. It still has total mass
deficit \(\delta_{12}(A)\) in labels 1 and 2, all held by label 3.
Choose reflection-invariant subsets of label 3 of the missing
masses and transfer them back. This can be done by cutting its
nonatomic radial measure at two quantiles. For any two jointly
Gaussian observations, the agreement indicator changes only if at
least one observation changes labels, so either transfer loses at
most twice its changed mass in stability. The rearrangement itself
loses nothing. Thus
\[
        S_\rho(A)\le V_\rho^{\mathrm{ref}}(m)+4\delta_{12}(A).
                                                        \tag{RF1}
\]
The inequality holds for every pair by relabeling.

The difference \(V_\rho(m)-V_\rho^{\mathrm{ref}}(m)\) is positive and
stays bounded below by a positive number on compact subsets of the
positive parameter range. Otherwise take a sequence of reflected
near-maximizers with \((\rho_n,m_n)\to(\rho,m)\) at an interior point.
The value \(V_\rho(m)\) is continuous in \(\rho\) by operator-norm
continuity of \(T_\rho\) on interior bands, and in \(m\) by direct
mass transfer costing at most \(\|m-m'\|_1\) in stability. Pass to a
weak limit of the indicator vectors in Gaussian \(L^2\). The mass
constraint passes to \(m\), reflection invariance is weakly closed,
and compactness of \(T_\rho\) together with norm continuity gives
stability \(V_\rho(m)\) at the limit. A maximizer in the relaxed class
is ordinary, since a nontrivial two-sided balanced label splitting
would strictly increase the average quadratic objective by positive
definiteness of \(T_\rho\). This contradicts full moment rank of the
resulting reflected maximizer. Therefore (RF1) gives a uniform
positive lower bound on \(\delta_{12}(A)\) for actual maximizers.
[Positive-correlation compactness and relaxation](SMALL_POSITIVE_CORRELATION_INTERVAL.md#a-rigorous-continuation-boundary-in-correlation)

## A common circle for interval radial images

Let \(W_i\) denote the open strict-winner set for label \(i\) under the
adjusted scores. These sets give the cells up to Gaussian null sets.
For almost every radius, a strict winner meets the circle if and
only if its angular slice has positive measure. Here the positive
measure follows from openness of the strict winner, and the converse
uses agreement with the original slice almost everywhere by Fubini.
The positive radii met by \(W_i\) form an open subset of \((0,\infty)\).
If these three radial images are intervals, their intersections are
pairwise nonempty by the overlap result. Three pairwise intersecting
open intervals have a nonempty open common intersection. Every
radius there meets all three strict winners.
[Analytic scores and the strict-winner sets](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)
