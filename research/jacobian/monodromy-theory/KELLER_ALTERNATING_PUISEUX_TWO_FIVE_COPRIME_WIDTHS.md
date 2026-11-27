# Coprime widths are impossible in the smooth cuspidal alternating profile

**Theorem SixSheetSmoothCuspCoprimeWidthsExclusion.**
Let a plane Keller map of generic sheet degree six have the
alternating profile with no unramified omitted primes. In every
polynomial target frame the two coordinate degrees on the
normalization of the sole ramified curve have a common divisor
greater than one.

**Corollary SixSheetSmoothCuspWidthAtLeastSix.**
Every branch projection width in this profile is composite and
at least six. The ramified curve has ordinary plane degree at least
eight in every polynomial target frame.

The [six-sheet local profiles](KELLER_BOUNDARY_CUSP_PROFILES_AND_OREVKOV_INVARIANT.md#consequences-in-generic-sheet-degree-six)
give the ramified curve \(B\) an affine-line normalization and
exactly one singular analytic branch of type \((2,5)\). Its
source germ has local degree five, so there is no second
analytic branch at its image. At any other self-identification
there are two smooth branches, because three distinct
degree-three germs would exceed the sheet degree. Each local
monodromy image there is abelian, with 3-cycles acting on
disjoint sheet sets. Smooth meridians of \(B\) act by single
3-cycles and the global representation is transitive.

Suppose the two widths are coprime and swap the coordinates
if necessary so their widths \(n,d\) have \(n\) odd.
Both widths are at least five by the
[slice genus bound](KELLER_SIX_SHEET_CUSP_SLICE_GENUS.md#smooth-normalization-with-one-25-affine-cusp).
Choose the affine normalization parameter \(t\) so the cusp
lies at \(t=0\). The proof uses the two-meridian deformation
for the [ordinary cusp profile](KELLER_WIDTH_FIVE_ALTERNATING_ORDINARY_CUSP_EXCLUSION.md),
with the aligned coefficient constraints for the \((2,5)\)
jet verified below.

## Coefficients preserving the \((2,5)\) cusp

For a polynomial pair \(\Phi=(P,Q)\) of maximum degrees \(n,d\),
write \(C_j\in\mathbb C^2\) for its coefficient of \(t^j\).
A cusp of type \((2,5)\) at \(t=0\) is equivalent to
\[
 C_1=0,\quad C_2\ne0,\quad C_3=aC_2,\quad
              C_2\wedge(C_5-2aC_4)\ne0
\]
for a scalar \(a\). In coordinates along and normal to \(C_2\),
the coordinate along \(C_2\) has order two and its third
coefficient is \(a\) times its second. The normal coordinate
has no terms below order four. Subtracting a multiple of the
square of the order-two coordinate cancels the fourth-order
term of the normal coordinate and leaves its order-five term
nonzero exactly under the displayed wedge condition. After
making the order-two coordinate \(s^2\), subtracting the even
part from the normal coordinate and dividing by an analytic
unit in \(s^2\) gives normal coordinate \(s^5\). These analytic
changes can be made along any sufficiently small perturbation
of the displayed jets. In particular they preserve the cusp
ball and its boundary link up to isotopy.

The free parameters are \(a,C_0,C_2,C_4,C_5\) and the permitted
higher coefficients. Imposing the nonvanishing conditions,
including nonzero leading coefficients in both coordinates,
gives an open subset of an affine parameter space. Every
pair there is finite and birational onto its image. Indeed
the parameter ring is finite over \(\mathbb C[P]\), and the
field degree over \(\mathbb C(P,Q)\) divides both \(n\) and
\(d\).

We check that a nonempty Zariski open set of the parameter
space has no other nonimmersed point, no second preimage of
the cusp, and only ordinary nodes away from the cusp. Put
\(\ell_a(t)=t^2+a t^3\). For fixed \(a\), the three
functions \(\ell_a(t),t^4,t^5\) have freely varying
two-vector coefficients. Let \(m\) be the dimension of the
full parameter space, including \(a\). A vanishing derivative
at a nonzero parameter imposes two conditions since the
derivative of \(t^4\) there is nonzero. The equality of its
image with \(\Phi(0)\) likewise imposes two. The incidence
images of these bad conditions thus each have dimension at
most \(m-1\).

For distinct nonzero \(t,u\), equality of images imposes
two conditions since \(t^4-u^4,t^5-u^5\) cannot both be zero.
Consider the determinant of the three linear functionals
\[
          f\longmapsto f(t)-f(u),\qquad
          f\longmapsto f'(t),\qquad f\longmapsto f'(u)
\]
on the ordered basis \((\ell_a,t^4,t^5)\). It is nonzero
as a polynomial in \(a,t,u\), for instance its value is
\(132\) at \(a=0,t=1,u=2\). On its nonzero locus the two
tangent vectors vary freely under the equality of images,
so tangency imposes a third condition. The bad incidence
has dimension at most \(3+(m-1)-3=m-1\).
On the zero locus the variable triple \((a,t,u)\) has dimension
at most two, and equality already gives incidence dimension
at most \(2+(m-1)-2=m-1\).

For three distinct nonzero parameters \(t,u,v\), use the two
evaluation differences at \(u,v\) relative to \(t\), acting
on the same three functions. When the difference map has
rank two, equality of all three images gives four conditions,
for incidence dimension at most \(4+(m-1)-4=m-1\).
The rank-one locus lies over the collinearity equation for
the three points of \((t^4,t^5)\), a proper hypersurface
in \((t,u,v)\). The difference columns of \(t^4,t^5\) have
rank one there, since these functions together distinguish
distinct nonzero parameters. For each such triple, at most
one \(a\) makes the difference column of \(\ell_a\) lie in
that line. Otherwise both the \(t^2\) and \(t^3\) difference
columns would lie in it, making the three points of
\((t^2,t^3,t^4,t^5)\) collinear. This is impossible since
the evaluation matrix of \(t^2,t^3,t^4\) at three distinct
nonzero parameters is a Vandermonde matrix multiplied by
invertible row scalars. Thus the rank-one locus in
\((a,t,u,v)\) has dimension at most two, and its bad
incidence has dimension at most \(2+(m-1)-2=m-1\).
Removing the closures of all these incidence images gives
the asserted open set. It is path connected and dense in
the full affine parameter space.

## Transport and the two-meridian test

The [proper-isotopy argument at infinity](KELLER_WIDTH_FIVE_ALTERNATING_ORDINARY_CUSP_EXCLUSION.md#a-coprime-degree-deformation)
applies unchanged to these coprime degree pairs. For large
\(|P|\) the \(n\) points on a line separate uniformly by
their leading terms \(Q_d(P/P_n)^{d/n}\zeta^d\), with
\(\zeta^n=1\), along any compact path of fixed degrees.
A small perturbation into the open nodal locus carries the
original monodromy to a representation of the perturbed
complement. The
[meridian-linking gluing](KELLER_WIDTH_FIVE_ALTERNATING_ORDINARY_CUSP_EXCLUSION.md#transporting-the-monodromy)
applies at each immersed two-branch ball since its local image
is abelian. The equisingular charts above preserve the
\((2,5)\) cusp ball with its full local representation.
No change is made to the exterior representation up to the
induced isotopy. It remains surjective because one generic
\(P\)-line disjoint from the singular balls originally
generates the global image. Path connectedness of the open
locus therefore gives a transitive representation by
3-cycle meridians for any of its members.

The degenerate pair \(P_*(t)=t^n,\ Q_*(t)=t^2\) belongs
to the affine parameter space with \(a=0\). It may lie on a
leading-coefficient or cusp-jet vanishing locus, which is
permitted for this limit. Choose a pair in the open nodal
locus arbitrarily close to it. The
[bounded braid degeneration](KELLER_WIDTH_FIVE_ALTERNATING_ORDINARY_CUSP_EXCLUSION.md#a-degeneration-with-two-meridians)
then gives a meridian-preserving surjection from the complement
group of \(B_*=\{x^2=y^n\}\) to the chosen complement group.
Here the only test input is that the braids on fixed compact
loops in the \(x\)-plane persist, so the lower degree of
\(Q_*\) does not affect the surjection. Projection of \(B_*\)
to \(y\) has width two. The transitive six-sheet image would
be generated by just two 3-cycles, which is impossible.
This proves the theorem.

For the corollary, if one width were prime, repeatedly shear
the other coordinate by powers of the first to cancel its
leading degrees divisible by that prime. The final degree
cannot be zero, since that would make the branch a
coordinate line and the first width one. The resulting
coprime pair contradicts the theorem, and the existing
lower bound of five then gives width at least six.

In a fixed frame, the ordinary plane degree of \(B\) is the
larger of its two coordinate degrees on the normalization.
A generic affine line cuts \(B\) transversely in smooth
points, and its affine-linear equation pulls back to a
polynomial of that larger degree. Ordinary degree six would
force both coordinate degrees to be six. An invertible
linear target change can cancel one of their leading
coefficients, producing a smaller width, a contradiction.
Ordinary degree seven would give a generic linear
coordinate of prime width seven, also impossible.
