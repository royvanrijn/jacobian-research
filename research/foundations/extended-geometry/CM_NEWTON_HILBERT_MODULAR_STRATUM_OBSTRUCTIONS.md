# Hilbert modular geometry of the split one-failure case

Let \(F\) be a totally real cubic field and consider a Hilbert
modular threefold for \(F\) with a fixed neat auxiliary level and
polarization type. At a good rational prime \(p\) split completely
in \(F\), each of the three labelled \(p\)-divisible factors has
height two and dimension one. Its slopes are either \((0,1)\) or
\((1/2,1/2)\). A fixed number-field point with full geometric
endomorphism algebra \(F\) can have exactly one supersingular
factor, as the [example at \(563\)](CM_NEWTON_DISCONNECTED_RM_ORDINARITY.md)
shows. The results below constrain two ways of trying to detect
this reduction by a bounded amount of characteristic-zero geometry.

## The one-failure divisor varies in degree and quaternionic datum

Work over a number field \(K\) containing a normal closure of \(F\)
and, after extension, the field of definition of a chosen geometric
component \(X\). Fix a smooth projective toroidal compactification
\(\bar X\) over a localization of \(\mathcal O_K\), and a
relatively very ample line bundle \(H\) on it. These choices
can be made in characteristic zero and extended after removing
finitely many primes, with the semiabelian scheme and its
\(\mathcal O_F\)-action. Write \(\omega_i\) for the three Hodge
lines extended by this semiabelian scheme. Toroidal
compactification and the extended scheme are the construction of
[Lan–Suh, Proposition 4.2(1)](https://www.kwlan.org/articles/van-tor-aut-gen.pdf#page=22).

At an eligible prime \(v\) of \(K\) over split \(p\), Frobenius
fixes the labels of the three lines. The \(i\)-th partial Hasse
invariant is a section of \(\omega_{i,v}^{\otimes(p-1)}\).
Its reduced zero divisor \(D_{i,v}\) detects the supersingular
factor at label \(i\). At each cusp the invariant has expansion
one, so this divisor lies in the interior of the toroidal model.
The intersections with the other labelled divisors have
codimension two in the threefold when nonempty. Hence the open
stratum with just label \(i\) failing is dense in \(D_{i,v}\).
[Goren, Theorem 2.1(1)–(2) and proof](https://www.math.mcgill.ca/goren/PAPERSpublic/PartialHasseInvariants.pdf#page=8)
· [Andreatta–Goren, Corollary 8.18](https://www.math.mcgill.ca/goren/PAPERSpublic/congruences.pdf#page=35).

**Proposition 1.** Write
\[
  c_i=(c_1(\omega_{i,K})\,c_1(H_K)^2)[\bar X_K].
\]
At all eligible \(v\) the degree of \(D_{i,v}\) with respect to
\(H_v\) is exactly \(c_i(p-1)\). If \(D_{i,v}\) is nonempty at
one such \(v\), then \(c_i>0\). In that case any collection of
effective divisors in the fiber whose union contains the full
one-failure stratum has total degree at least \(c_i(p-1)\).

**Proof.** Intersect the divisor class
\((p-1)c_1(\omega_{i,v})\) with \(c_1(H_v)^2\).
Top intersection numbers of line bundles in the smooth
projective model are invariant under specialization, for example
by the coefficients of the multivariable Hilbert polynomial.
An effective nonzero divisor has positive degree with respect
to an ample line. Every divisor collection containing the open
stratum contains its closure, so at least the fundamental cycle
of \(D_{i,v}\) occurs in the sum of the effective divisor
cycles. \(\square\)

Thus uniform coverage of a full one-failure stratum requires
degree at least linear in \(p\), even if the covering divisors are
chosen afresh in each fiber. A fixed finite list of horizontal
divisors has constant total specialized degree away from
finitely many places. This is a different requirement from
covering just the reduction of an individual section.

The geometry in Tian–Xiao also distinguishes split primes from
the inert examples with projective-line fibers. In their
Theorem 5.2 take the initial quaternion algebra split and
\(\mathtt T=\{i\}\). At a split prime the Frobenius orbit at
\(i\) is a singleton, so it falls under Case \(\beta 2\) of
§5.1. It adds \(i\) and its corresponding finite place
\(\mathfrak p_i\mid p\) to the quaternion ramification set.
The bundle index \(I_{\mathtt T}\) is empty. Thus the closed
stratum in that theorem is an auxiliary quaternionic Shimura
surface in characteristic \(p\), with quaternion algebra
ramified at \(i\) and \(\mathfrak p_i\) and split at the other
places. The quaternion algebra in this description changes
with \(p\). The theorem identifies a vertical stratum and its
prime-to-\(p\) Hecke action. Arithmetic-intersection theorems
for special divisors with a given horizontal defining condition
use different input.
[Tian–Xiao, §5.1 Case \(\beta 2\) and Theorem 5.2](https://arxiv.org/html/1308.0790v2#S5.SS1)
· [Shankar–Tang, Theorem 1](https://arxiv.org/abs/1706.08154).

The fixed example at \(563\) even has geometrically simple
reduction. In its cyclic cubic coefficient field, the
quadratic discriminant \(a^2-4p\) has valuation one at the
selected failing place and zero at both other places. The
three conjugate discriminants thus give independent square
classes and three independent root swaps. The ratio of the
roots of a particular quadratic is nontorsion, because at a
place over \(p\) where that conjugate pair is ordinary its
slopes are zero and one. A ratio between roots of different
quadratics is also nontorsion, as a swap in just one pair
would otherwise contradict the within-pair assertion.
The trace \(a=13e^2+2e-57\) generates the cubic field.
Thus each root and each of its positive powers has six
distinct conjugates. The
Frobenius polynomial remains irreducible after every residue
extension, proving geometric simplicity. Consequently the
geometrically nonsimple reduction detected by Shankar–Tang's
method is a different phenomenon from this one-failure
Newton polygon.

## A fixed-point obstruction using finitely many Hecke translates

**Proposition 2.** Fix a number-field point \(x\in X\), a prime
\(\ell\) away from the tame level, and an integer \(D\).
There are finitely many arithmetic \(\ell\)-power Hecke
translates \(y_1,\ldots,y_m\) of \(x\) and one finite number
field \(K'\) over which they and the isogenies are defined,
with the following property. At every place \(w\) of \(K'\)
outside a finite set, and every \(1\le n\le D\), no nonzero
section of \(H_w^{\otimes n}\), even after extending the
residue field algebraically, vanishes on all the \(y_{j,w}\).
At good split primes with residue characteristic different from
\(\ell\), these points all have the same labelled Newton
pattern as \(x_w\).

**Proof.** The \(\ell\)-power orbit with fixed tame level is
Zariski dense in the complex component. Here is an explicit
density argument. On its analytic uniformization by
\(\mathbb H^3\), take the upper and lower unipotents with
parameters in sufficiently divisible fractional
\(\mathcal O_F[1/\ell]\)-ideals. They preserve all lattice
and level conditions away from \(\ell\). Each such parameter
set is dense in \(\mathbb R^3\), since it includes rescalings
by \(\ell^{-r}\) of a full real lattice for all \(r\geq0\).
The closures contain both upper and lower real unipotent
groups in each factor. These generate
\(\mathrm{SL}_2(\mathbb R)^3\), which acts transitively on
\(\mathbb H^3\). Thus the orbit is analytically dense and
therefore Zariski dense. The rational matrices defining
these correspondences give real-multiplication-compatible
isogenies of \(\ell\)-power degree. The corresponding points
and isogenies are algebraic, since their finite kernels and
level data are defined over finite extensions.

Over \(\overline K\), evaluate the global sections of
\(H^{\otimes n}\) on this orbit. If the evaluations at a
finite set have a nonzero kernel, some section of that
kernel is nonzero at another orbit point by density.
Adding its evaluation reduces the kernel dimension.
After finitely many steps evaluation is injective. Repeat
for the finitely many \(n\leq D\) and take the union of the
chosen translates. They and their isogenies spread out
over one localization of \(\mathcal O_{K'}\). Coherent
base change holds for the finitely many spaces of global
sections after further localization, and an injective
evaluation minor is a unit outside finitely many places.
Thus the claimed injectivity persists, also after residue
field extension. Finally a prime-to-\(p\) isogeny induces
an isomorphism of \(p\)-divisible groups, compatible with the
\(\mathcal O_F\)-idempotents. This gives the Newton statement.
\(\square\)

In particular, if a split reduction of \(x\) has exactly one
nonordinary factor, the same is true for this finite set of
its isogeny relatives. Even a hypersurface of bounded degree
chosen separately in each residue characteristic cannot contain
all of them outside finitely many places. The same applies to a
union of hypersurfaces with bounded total degree by taking the
product of their equations. An approach that uses a
bounded-degree algebraic condition on a single preferred
representative must exploit why this representative is preferred.
The argument uses no density assertion about reductions of
arithmetic points. Characteristic-\(p\) Hecke orbit results
compare points inside a fixed fiber and preserve the labelled
\(p\)-divisible isomorphism class.
[Chai–Oort, Hecke orbit survey, §§2–3](https://www2.math.upenn.edu/~chai/HOC_survey/HOC_v3b.pdf).

## Curve lifts and the remaining arithmetic estimate

For a fixed arithmetic curve through \(x\), the local geometry
also gives a simple-root test. Suppose the tangent of the
curve at \(x\) has nonzero projection on each of the three
Kodaira–Spencer summands. Outside finitely many primes this
persists on reduction. At a split reduction where label \(i\)
fails, Serre–Tate deformations split into the three labelled
one-dimensional deformations. The partial Hasse invariant at
\(i\) has a nonzero first derivative only in the \(i\)-th
direction at its zero, by the local partial Hasse formula.
Its restriction to the curve thus has a simple root at the
fixed arithmetic point on reduction. The condition on the
tangent can always be arranged for a local algebraic curve
through \(x\), by taking a line with all three coordinates
nonzero in an étale coordinate chart. This rules out a
forced multiple root at the one-failure point as a source
of additional divisibility.
[Andreatta–Goren, Corollaries 8.17–8.18](https://www.math.mcgill.ca/goren/PAPERSpublic/congruences.pdf#page=35).

Bogo–Li's Theorems 1 and 2 construct characteristic-zero
twisted modular lifts of powers of partial Hasse invariants
on noncompact genus-zero Kobayashi geodesics with modular
embedding. Their Equation (13) makes the relevant Hasse
polynomial degree grow linearly in \(p-1\) in the split case
for every positive Lyapunov exponent. Equation (10) is a
congruence of formal expansions at a cusp. At a different
fixed arithmetic point with coordinate \(t_0\), all but a
finite number of reductions of \(t_0\) are outside the
residue disk of that cusp. The formal power-series
congruence does not by itself evaluate there. The polynomial
\(\alpha_{\mathfrak p,j}(t_0)\) does evaluate, and the missing
input is a bound on the primes at which one of those
evaluations is zero. The theorem describes the polynomial
and does not give that bound.
[Bogo–Li, Equations (10)–(13), Theorems 1–2](https://arxiv.org/html/2601.16929v1).

There is an algebraicity obstruction to extending the
cusp-expansion argument by a fixed algebraic function.
Suppose Equation (10) uses the same normalized formal
series \(G=F_j^N\) over a fixed number field for infinitely
many degree-one places over split \(p\), and suppose
\(\lambda_j>0\). Then \(G\) cannot be algebraic over the
rational function field of the curve. For if it were algebraic,
remove a fixed finite set of \(t\)-values so that \(G\) is a
unit on the finite cover obtained by adjoining it. One can
use a monic polynomial equation for \(G\) and invert its
constant term and the coefficient denominators to construct
this open set. The same unit statement persists in every
remaining special fiber. At the split degree-one primes
Equation (10) becomes
\(\alpha_{\mathfrak p,j}^{\,l_j^2}=\bar G^{\,1-p}\)
in the function field of the corresponding reduced branch.
This follows from \(j'=j\) and the Frobenius identity for
series over \(\mathbb F_p\). Every branch of the finite cover
has \(G\) invertible above the open set, so \(\alpha_{\mathfrak p,j}\)
can only vanish at reductions of the removed \(t\)-values.
Equation (12) bounds the multiplicity of every root by a
constant depending only on the elliptic orders and the
divisor of \(\varphi_j'\). This gives a uniform degree bound
contrary to Equation (13) and \(\lambda_j>0\).
[Bogo–Li, §1.3](https://arxiv.org/html/2601.16929v1#S1.SS3).

The special control of elliptic points in Bogo–Li's §2.3
covers a different stratum in characteristic zero. On the
Hilbert modular quotient with labelled \(F\)-action, a
nontrivial elliptic stabilizer produces a polarized
automorphism commuting with \(F\) other than \(\pm1\).
Its action on rank-two \(F\)-homology is represented by
a noncentral finite-order matrix of determinant one.
At every real embedding its eigenvalues are nonreal
conjugate roots of unity, so it generates a CM quadratic
extension of \(F\) acting on the fiber. Such a point is of
CM type. The extra control at elliptic points is therefore
outside the independently acting real-cubic case.
