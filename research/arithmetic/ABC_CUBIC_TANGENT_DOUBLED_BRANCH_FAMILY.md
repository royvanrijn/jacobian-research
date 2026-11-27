# The doubled cubic branch family and its missing virtual units

The remaining genus-five correspondence for two rational cubic
tangent parameter maps has an explicit rational model.
It gives simultaneous positive directions and independent
cubic norm classes for infinitely many quadratic fields.
A valuation condition across all specializations prevents
two genuine prime triangles in this family.

## Classification of the doubled branch model

**Theorem DB1.** Suppose two rational cubic weighted tangent
parameter maps have a common nonsquare genus-five quadratic
function field and are not projectively related.
Their elliptic models have rational halves of their deck centers
and have distinct inverse coefficient classes modulo rational cubes.
After translating by those halves and scaling the models to
\[
 E_b:y^2=x^3+b^2,\qquad E_{b^2}:v^2=z^3+b^4,
 \qquad m=b^2,
\]
the maps on their common field, up to an independent sign of each
ordinate and a rational change of source coordinate, are
\[
 \begin{aligned}
 w&=ms^3, & Y^2&=(w^2-1)(w^2+4w-1),\\
 x&=\frac{2ms^2}{w-1}, & y&=\frac{bY}{(w-1)^2},\\
 z&=-\frac{2ms}{w+1}, & v&=\frac{mY}{(w+1)^2}.
 \end{aligned} \tag{1}
\]
Conversely (1) gives cubic parameter maps with these properties
whenever the weighted profiles admit the stated scaling and halves.

By [CP3–CP6](ABC_CUBIC_PARAMETER_TANGENT_OBSTRUCTION.md#the-simple-branch-patterns)
only the doubled irrational pattern remains. In each target,
translate the rational half of the deck center to the origin.
The rational branch value becomes the pole of the elliptic
abscissa. The other three branch values have a common \(S_3\)
splitting field, and cubic Kummer theory in
\(\mathbb Q(\zeta_3)\) makes the two coefficient classes equal
or inverse, as in GT5.

In these abscissas write
\(x=A/B,\ z=C/D\) with pairs of relatively prime binary cubics.
There are three distinct simple poles for each map and the two
pole sets are disjoint. Each pole maps to a finite branch point
in the other abscissa. The six source points giving the doubled
irrational diagonal are the zeros of a separable binary sextic
\(H\). In suitable normalization, equality of quadratic
squareclasses gives
\[
 A^3+b_1^2B^3=D H,\qquad
 C^3+b_2^2D^3=\lambda^{-2}B H,\qquad \lambda\in\mathbb Q^\times .
 \tag{2}
\]
Indeed the binary degree-twelve polynomials given by multiplying
the left sides by \(B,D\), respectively, have exactly the common
twelve simple branch roots. Their ratio is a rational square.

If the coefficient classes are equal, normalize them to
\(b_1=b_2\). The matching of the finite branch values is then
the identity, by uniqueness of the equivariant identification.
Hence \(AD-CB=eH\) for a nonzero rational \(e\).
At a real zero of \(B\), equations (2) give \(e=(D/A)^2>0\).
At a real zero of \(D\), they give
\(e=-B^2/(\lambda^2 C^2)<0\).
Each real binary cubic has a real projective zero.
The encodings by projective cubics and their nonvanishing at
the respective poles justify these divisions also at infinity.
The contradiction excludes equal coefficient classes.

In the inverse case put \(b_1=b,\ b_2=b^2\), so \(m=b^2>0\).
The unique matching of finite branches now sends \(x\) to
\(m/x\), since \(x^3=-m\) gives \((m/x)^3=-m^2\).
Thus \(AC-mBD=eH\), again with \(e\in\mathbb Q^\times\).
It follows from (2) that the image of the parameter line is the
bidegree-\((3,3)\) curve
\[
 (x^3+m)(z^3+m^2)=k(xz-m)^2,\qquad k=(e\lambda)^{-2}>0 .
 \tag{3}
\]
The two cubic parameter fields generate the source field, so
the image is birational to the source and is geometrically
irreducible of arithmetic genus four.

There are three singularities of (3) at the finite branch
matches, where \(x^3=-m,\ z=m/x\).
In coordinates \(x=x_0(1+u),z=z_0(1+v)\) about one such point,
the quadratic term after dividing by \(m^2\) is
\(9muv-k(u+v)^2\).
Its discriminant vanishes only at \(k=9m/4\).
In that case the common tangent line \(u=v\) meets (3)
with order three at the point, by substituting into the
equation. Two smooth branches sharing that tangent would give
order at least four. Both branches in the required doubled
pattern are smooth because both parameter maps are unramified
at the source branch points. Hence \(k=9m/4\) is excluded.
For other values of \(k\) the three fixed singularities are
ordinary nodes.

At infinity the only points of (3) have one abscissa infinite
and the other a simple finite branch root. They are smooth,
as is seen from the derivative in the coordinate of the
finite root. At an affine singular point different from the
three matched branch points, the two partial derivatives
give \(z^3=mx^3\). If either abscissa is zero, both are zero
and (3) gives \(k=m\).
Otherwise write \(U=xz/m,\ V=x^3/m\).
The equations give \(V^2=U^3\) and \(V(U-3)=2U\), so
\((U-1)^2(U-4)=0\).
The value \(U=1\) is one of the three branch matches.
The value \(U=4\) gives \(V=8,\ k=9m\), and supplies three
additional distinct singular points. That many singularities
is impossible for a geometrically irreducible curve of
arithmetic genus four.
At least one additional singularity is needed for a rational
normalization, so \(k=m\) is necessary.

For \(k=m\), the two abscissa functions of (1) give a
rational parametrization of (3), whose inverse is
\[
 s=\frac{z^2+mx}{z(xz-m)} .
\]
For example eliminating \(x\) from these identities leaves
the linear factor \(z(1+ms^3)+2ms\) after removing the factors
at the base points. The two functions are of degree three
and their image therefore fills (3). Substituting them in
\(x^3+m\) and \(z^3+m^2\) gives the square factors and ordinates
in (1). The twelve roots of the new field polynomial are
simple because \(w=ms^3\) is unramified at the four
distinct nonzero finite roots of \((w^2-1)(w^2+4w-1)\).
The resulting quadratic field is of genus five. This also
proves the converse, including nonprojective relation by the
explicit rational inverse.

## Local failure at the coefficient primes

**Theorem DB2.** In (1), for every rational parameter with a
nonsplit quadratic field, let the two cubic norm classes be
represented by
\[
 \alpha=y+b,\qquad \beta=v+m .
\]
If an odd prime \(q\) has \(v_q(b)\ne0\bmod3\), both classes
fail the local virtual-unit condition at the two split primes
over \(q\). If \(v_2(b)\ne0\bmod3\), exactly one of the two
classes fails the condition at the two split primes over \(2\).
In particular at least one fails for every allowable weighted
coefficient \(b\).
If both a dyadic and an odd such prime occur, every nonzero
class in the two-dimensional span fails the virtual-unit test.

The representatives follow from the order-three marked divisor
and the cube-norm comparison of
[GT2](ABC_RECIPROCAL_GENUS_TWO_TANGENT_OBSTRUCTION.md#independent-norm-functions-on-the-reciprocal-cover).
They have norms \(-x^3\) and \(-z^3\).
Let \(a=v_q(b),\ e=v_q(s)\) and \(d=2a+3e\).
When \(3\nmid a\) one has \(d\ne0\) and
\(v_q(w\pm1)=\min(d,0)\). The two relative valuations are
\[
 \begin{aligned}
 v_q(x^3/b^2)&=3v_q(2)+2d-3\min(d,0),\\
 v_q(z^3/b^4)&=3v_q(2)+d-3\min(d,0).
 \end{aligned}
\]
They are positive, and at \(2\) are at least four. Consequently
the normalized square roots \(y/b\) and \(v/m\) lie in
\(\mathbb Q_q\), by the principal unit square test at odd
\(q\) or the \(1\bmod8\) square test at \(2\). The quadratic
field splits there.

For odd \(q\), choose the sign of \(y/b\) congruent to one
modulo \(q\), and separately the sign of \(v/m\) with this
property. The valuations of \(\alpha,\beta\) on those branches
are \(a,2a\), respectively. Their conjugate valuations are
the negatives modulo three by the cube norms.
At \(2\), choose each normalized root congruent to one modulo
four. The valuations on those branches are \(a+1,2a+1\)
and again their conjugate valuations are the negatives modulo
three. For \(a=1,2\bmod3\) exactly one of \(a+1,2a+1\)
vanishes modulo three. At an odd coefficient prime both
entries are nonzero. These two valuation rows have rank
two over \(\mathbb F_3\), regardless of the ordering of the
split primes. This proves the claimed independence and the
absence of any nonzero virtual unit in the span.
The excluded source endpoints \(s=0,\infty\) have split
quadratic algebra and the branch parameters have zero
radicand. Finally \(b\) is noncubic, so some prime valuation
is nonzero modulo three.

## Simultaneous positive independent directions

**Corollary DB3.** The profiles \(k_1=32/343,\ k_2=3087/4913\)
admit two degree-three rational parameter maps with the same
quadratic field in every specialization. On a fixed real
interval the two directions are positive. There are infinitely
many distinct imaginary quadratic fields on that interval,
with two independent norm classes whose span has zero
virtual-unit subspace. The two generic cubic norm classes
are independent even geometrically.

Use \(b=6\) in (1), and the model scalings and halves of GT2.
For example \(1/5<s<1/4\) gives
\(-2+\sqrt5<w<1\), so the common radicand is negative.
Whenever the radicand is negative, the translated abscissa
on each factor is less than \(-|b_i|^{2/3}\), where
\(b_i=b,b^2\). Each marked abscissa comes from the real
point \(D_j-V_i\), using the rational translation \(V_i\),
so it is at least \(-|b_i|^{2/3}\). At abscissa infinity the
parameter is the rational branch point \(V_i\); by GT1 its
three directions have the same nonzero sign up to common
scaling. Along the negative ray up to the indicated limit,
no marked zero is crossed. Thus both weighted directions can
be made positive throughout the interval. This argument is
about the real parameter line even though the field is
imaginary.

At \(2,3\) the coefficient \(b=6\) has valuation one.
By DB2 every rational specialization on the interval has two
independent classes, with zero virtual-unit subspace.
For any fixed squareclass \(d\), those parameters come from
the genus-five twist \(dZ^2=(w^2-1)(w^2+4w-1)\) with \(w=36s^3\).
By [Faltings](sources/faltings-1983-endlichkeitssaetze.pdf) the
twist has finitely many rational points. Hence the fields vary
infinitely.
The generic classes are independent also over the algebraic
closure: if a nonzero linear combination of their divided
divisors were principal geometrically, it would be principal
over \(\mathbb Q\) by descent. The corresponding function
classes would differ from a cube by a rational constant,
which would be a cube by taking the quadratic norm.
That identity would contradict the specialized independence
away from its finite exceptional set.

An exact instance at \(s=1/4\) has field
\(\mathbb Q(\sqrt{-2807})\) and normalized abscissas
\(-72/7,-288/25\). The two profile parameters are
\(-624/1715\) and \(10231200/4721393\).
The classes of
\(6+(30/49)\sqrt{-2807}\) and
\(36+(36/125)\sqrt{-2807}\) have independent valuation
rows at \(2,3\). The
[checker](verify_cubic_tangent_doubled_branch_family.py)
replays the parametrization, the norm identities, the
two exact profile parameters, and the split-prime valuations.

## Consequence for actual prime witnesses

**Corollary DB4.** No two rational parameter maps of positive
degrees at most three into weighted tangent profiles produce
a common nonsquare quadratic field with two independent generic
marked norm classes and simultaneous genuine prime triangles.

If either map has degree one the common curve has genus one.
The maps on Jacobians and the weighted rational Hom calculation
of GT5 force the second map to have square degree, hence also
one in this range. The classes then agree up to inversion.
Degrees two and three together are excluded by CP1.
Two maps of degree two with independent directions cannot both
have positive directions, by
[QP3](ABC_QUADRATIC_PARAMETER_TANGENT_OBSTRUCTION.md#consequence-for-rational-parameter-maps).
For two maps of degree three, CP1–CP6 exclude independent
classes below genus five and outside the inverse-coefficient
doubled pattern or in the projective case.
Theorem DB1 classifies the remaining pattern, and DB2 prevents
the simultaneous virtual-unit condition there.
The cubic norm class of a genuine triangle is a virtual unit by
[CT1](ABC_CUBIC_TANGENT_NORM_RANK_AND_CONDUCTORS.md#the-tangent-identity-and-the-lost-cubic-direction)
and the [weighted lift](ABC_CUBIC_TANGENT_RING_CLASS_LIFTS.md#the-fixed-slope-rank-one-field-and-cofactor).
The unrestricted occupied rank-one norm record count is open.
