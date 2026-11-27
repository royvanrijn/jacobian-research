# A pure quartic with nonzero binary cubic in a Keller component

**Theorem KellerPureQuarticNonzeroBinaryCubic.** Let \(F\) be a
three-variable complex Keller map of ordinary degree at most five,
and suppose one component \(R\) has \(R_4=z^4\) in a linear
source frame. If \(R_3(x,y,0)\ne0\), then \(F\) is invertible.
Consequently the primitive degree-three branch of the
[rank-two repeated-line reduction](QUINTIC_KELLER_RANK_TWO_REPEATED_LINE_PENCIL.md)
is invertible.

Put \(D(x,y)=R_3(x,y,0)\) and homogenize a fibre \(R=t\) using
a new variable \(w\). Its equation is
\[
       z^4+wR_3(x,y,z)+w^2R_2(x,y,z)+w^3R_1(x,y,z)
                                           +w^4(R_0-t)=0 .
\]
The affine fibre is smooth. At infinity the only possible
singularities are at the zeros of \(D\) on the projective line
\(z=w=0\), because the \(w\) derivative there equals \(D\).
Thus the projective surface has isolated singularities and is
normal. More explicitly, it is reduced and irreducible since a
multiple component would be singular along a surface and
distinct surface components in projective three-space would
intersect along a curve of singular points. The hypersurface
is Cohen–Macaulay and regular in codimension one.

At a simple root of \(D\), in a local projective chart with
coordinate \(u\) along the line, the quadratic jet is \(w\)
times a linear form with nonzero \(u\)-coefficient. It has rank
two. The
[isolated double-point argument](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md#ordinary-double-points-of-the-plane-leading-curve)
shows that the fibre residue form extends holomorphically after
resolution there. This use of the double-point argument requires
only the rank and isolatedness of the surface singularity.
If \(D\) is squarefree, the form therefore extends across the
whole boundary on a smooth resolution. The
[global exactness obstruction](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md#residues-on-a-resolution)
rules this out. The fibre residue is a nonzero exact form on
the affine fibre by the Keller identity, as in
[Wang–Yan, equation (3.28)](https://arxiv.org/pdf/2610.36495v1).

Now \(D\) has a unique repeated root, and a binary linear
change puts that point at \([0:1:0:0]\). In the \(y=1\) chart
its quadratic jet is
\[
                              \eta zw+\theta w^2 .
\]
The \(xw\) coefficient vanishes because the root is repeated.
If \(\eta\ne0\) this again has rank two. Since every other
root is simple, the same global obstruction applies. If
\(\eta=0,\theta\ne0\), the
[quadratic source theorem](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md#the-rank-one-normal-jet-of-a-quartic-component)
gives invertibility. In fact \(R\) then has constant nonzero
quadratic coefficient in \(y\), with the other terms affine
in \(y\).

It remains to consider \(\eta=\theta=0\). Here \(R_4\) has no
\(y\) terms, \(R_3\) has no \(y^3,xy^2,zy^2\) terms, and \(R_2\)
has no \(y^2\) term. Thus \(R=a(x,z)y+b(x,z)\), where
\(\deg a\le2\). If \(a\ne0\) is squarefree, the
[reduced quadratic source theorem](KELLER_QUARTIC_TRIPLE_BOUNDARY_AND_LINEAR_COMPONENTS.md#invertibility-with-a-reduced-quadratic-coefficient)
gives invertibility, including the case of constant \(a\).
If \(a\) is nonzero and nonreduced, it is a scalar multiple of
\(s^2\) for an affine linear \(s=s(x,z)\). In the noninvertible
case the
[parallel coefficient classification](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md)
then gives \(R=\mu+\lambda s+s^2H\) for constants \(\mu,\lambda\)
and a polynomial \(H\) of degree at most two. Let \(s_1\) be
the linear part of \(s\). The identity \(z^4=s_1^2H_2\) forces
\(s_1\) proportional to \(z\). But the cubic part of
\(s^2H\) is divisible by \(s_1\). This gives \(D=0\), contrary
to the hypothesis.

Finally, if \(a=0\), the binary cubic \(D=cx^3\) has \(c\ne0\).
For each fixed \(z\), \(b_x\) is a quadratic in \(x\) with
constant leading coefficient \(3c\). The coefficient of \(x\)
has degree at most one in \(z\), and the constant term has
degree at most two. Both roots, counted with multiplicity,
are \(O(|z|)\) as \(z\) tends to infinity. At either root
\(b_z=4z^3+O(|z|^2)\), since \(b-z^4\) has degree at most
three. The product of the two \(b_z\)-values is polynomial in
\(z\) by symmetry, as the leading coefficient of \(b_x\) is
constant. It has leading term \(16z^6\), so it vanishes at
some finite \(z\). This supplies a critical point of \(R\) and
contradicts the Keller condition.

In the primitive degree-three branch, the
[closed cubic case](QUINTIC_KELLER_RANK_TWO_CLOSED_CUBIC_PENCIL.md)
leaves \(R_4=z^4\), and the
[pure-quartic slice restriction](QUINTIC_KELLER_RANK_TWO_CUBIC_PURE_QUARTIC.md)
gives \(R_3(x,y,0)\ne0\) in any noninvertible case. The theorem
excludes that possibility. The same argument works in
characteristic zero by embedding a finitely generated
coefficient field into \(\mathbb C\) and descending the
unique formal inverse.
