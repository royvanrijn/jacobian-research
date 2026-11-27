# Cubic parameter changes and their common tangent covers

For rational cubic changes of weighted tangent parameter, the
common quadratic cover has genus \(1,2,4\) or \(5\). In genus one
or two the marked cubic norm directions are dependent. Genus four
admits only projectively related parameter maps. In genus five a
pair outside this projective case must match the two branch
splitting fields and satisfy a restricted branch incidence pattern.

The parameter maps below are nonconstant rational functions
\(\varphi_i(s)\in\mathbb Q(s)\) on one rational parameter line.
They land in the weighted tangent covers of
[WM1](ABC_WEIGHTED_CUBIC_TANGENT_ELLIPTIC_MARKS.md#the-double-cover-and-its-three-marks).
A common field means equality of the nonsquare quadratic extensions
obtained by pulling back the two Heron quartics. Write \(C\) for the
smooth projective curve of the common field. Its maps to the
elliptic models have degrees \(\deg\varphi_i\) and intertwine their
deck involutions with the hyperelliptic involution of \(C\) when
the genus is at least two.

## Degree-three ramification and low genus

**Theorem CP1.** For a cubic parameter map, the common cover has
genus in \(\{1,2,4,5\}\). Two rational parameter maps of respective
degrees two and three cannot have a common nonsquare quadratic field.
For two cubic maps with a common cover of genus one or two, the two
generic marked cubic norm classes agree up to inversion.

For a weighted model \(E_b:y^2=x^3+b^2\), the four branch points
of the quotient \(P\mapsto T-P\) are the four halves of \(T\).
For an asymmetric profile, \(b=4k(k+1)\) is noncubic by
[FM4](ABC_CUBIC_TANGENT_COLLISION_MOMENT.md#asymmetric-weighted-profiles),
and the Galois action on \(E_b[2]\) is \(S_3\). The affine Galois
action on the four halves lies in \(S_4\) with projection onto
this \(S_3\). Its intersection with the translation kernel
\((\mathbb Z/2)^2\) is zero or all of it, by irreducibility of the
linear action. Thus the group is a point stabilizer \(S_3\) or
all of \(S_4\), giving orbit degrees \(1+3\) or \(4\).
For symmetric profiles the branch polynomial of
[FC1](ABC_CUBIC_TANGENT_FIELD_COLLISIONS.md#the-fixed-elliptic-double-cover)
is irreducible by its reversed Eisenstein model, and its group is
\(D_4\) of order eight by
[FC4](ABC_CUBIC_TANGENT_FIELD_COLLISIONS.md#fixed-projective-relations-and-adjacent-centers).

For a degree-three map, each branch fibre contributes either
three odd-multiplicity preimages if unramified, or one if it
ramifies, by the three partitions of \(3\).
Thus the pulled-back double cover has \(4+2m\) branch points,
where \(m\) is the number of unramified branch fibres.
It has genus \(1+m\). The unramified branch fibres are Galois stable,
so \(m\in\{0,1,3,4\}\).
A quadratic parameter map gives genus at most three. Consequently
a pair of degrees two and three can only have common genus one
or two.

Over \(\mathbb Q\), two weighted elliptic models are isogenous
only if they are isomorphic, and each has endomorphism ring
\(\mathbb Z\), as shown in
[GT5](ABC_RECIPROCAL_GENUS_TWO_TANGENT_OBSTRUCTION.md#classification-with-two-degree-two-projections).
If \(C\) has genus one, its two elliptic projections induce
isogenies between the three Jacobians. After identifying the
elliptic factors, the two pullbacks are rationally proportional.
The ratio of their degrees is therefore a rational square.
This excludes degrees two and three, while two degree-three
pullbacks differ by at most a sign.

Suppose next \(C\) has genus two. Each elliptic pullback is injective
for map degrees two and three because the degree is prime.
Otherwise a nonzero point of that prime order gives a connected
étale cover of the target elliptic curve through which the map
factors over the algebraic closure. A prime degree would then
identify \(C\) with a genus-one curve.
If the two elliptic images in \(\operatorname{Jac}C\) are the
same, the restricted polarization has degree equal to the map
degree, so the degrees are equal and the pullbacks agree after
identifying the targets by an isomorphism. In genus one or two,
the marked cyclic order-three subgroup is unique on each target
by GT5. Equal pullback directions therefore yield the same generic
norm class up to inversion, using the divisor classes and then
the cube norms to remove a rational constant.

For distinct images in genus two, the sum map from the product to
\(\operatorname{Jac}C\) is an isogeny \(\lambda\).
Its pullback polarization has diagonal \((d_1,d_2)\), and
off-diagonal maps vanish unless the targets are isomorphic over
\(\mathbb Q\), in which case identify them and write these entries
as an integer \(n\). Positivity gives \(n^2<d_1d_2\), and
\(\deg\lambda=d_1d_2-n^2\). The kernel lies in the kernel of
the polarization matrix, and meets neither coordinate factor.
The following small torsion facts determine the possibilities.

* No \(E_b\) admits a rational cyclic subgroup of order five.
  Over \(\mathbb Q(\zeta_3)\), a putative invariant line in \(E_b[5]\)
  and its image by \([\zeta_3]\) are distinct invariant lines,
  since \(X^2+X+1\) is irreducible modulo five. Galois preserves
  both lines and commutes with \([\zeta_3]\), so it acts by
  scalars. The Weil pairing then makes the cyclotomic determinant
  a square, impossible since \(\mathbb Q(\zeta_3)\) is disjoint
  from \(\mathbb Q(\zeta_5)\).
* For an asymmetric profile, \(E_b[2]\) is irreducible as a Galois
  module. A rational subgroup of order eight would be cyclic,
  or of type \(\mathbb Z/4\times\mathbb Z/2\). The first gives
  a rational subgroup of order two and the second gives one by
  doubling. Both are impossible.
* For a symmetric profile, scale to \(E_1\). Any rational cyclic
  subgroup of order four would contain a point \(P\) with \(2P=T\)
  and rational \(x(P)\). The two possible abscissas are
  \(-1\pm\sqrt3\) by the duplication formula. Hence any rational
  subgroup of order eight must have type \(\mathbb Z/4\times\mathbb Z/2\).

For \((d_1,d_2)=(2,3)\), the zero off-diagonal case would make the
kernel inject into both \(E_1[2]\) and \(E_2[3]\), impossible for
its order six. If \(n=\pm1\), the kernel has order five and its
projection gives a forbidden rational cyclic subgroup. If \(n=\pm2\)
the kernel has order two, and the second row of the matrix forces
its second coordinate to vanish. This too violates injectivity.
Compare [Shaska, Theorem 1](sources/shaska-genus-two-elliptic-subcovers.pdf)
for the geometric genus-two locus with degrees two and three for
arbitrary elliptic factors. The exclusion here uses the specific
rational torsion of the weighted factors.

For \((d_1,d_2)=(3,3)\) and \(n=\pm2\), the kernel again would give
a rational cyclic subgroup of order five.
For \(n=\pm1\) its order is eight, and it is the graph of
\(P\mapsto-3nP\) on a rational subgroup of order eight of the
common factor. Indeed it lies in the kernel of the matrix, which
gives \((P,-3nP)\) with \([8]P=0\). Only the symmetric subgroup
type \(\mathbb Z/4\times\mathbb Z/2\) is possible, so the graph
lies on the elliptic diagonal \(P\mapsto(P,nP)\).
The polarization on that diagonal has degree eight as a line
bundle. Since \(\lambda\) has degree eight on the diagonal,
its image inherits a principal restriction, which splits off
and contradicts indecomposability of the genus-two Jacobian.
Thus the off-diagonal maps must vanish. The kernel now is a
Galois-equivariant graph identifying the two groups of 3-torsion
points. In particular it identifies their unique rational cyclic
subgroups of order three. The two marked pullback divisor
classes agree up to sign, and the cube-norm argument again
gives equality of norm classes up to inversion.

## The common rational branch point in genus four

**Theorem CP2.** If two rational cubic weighted tangent parameter
maps have a common cover of genus four, the two maps are related by
a rational Möbius transformation on the target parameters. Their
generic marked cubic norm classes agree up to inversion.

Here \(m=3\) in both profiles. Each branch set has a rational point
and an orbit of degree three, and only the rational branch fibre
ramifies. Above that branch point the cubic map has exactly one
odd-multiplicity point, which is rational. A point lying over an
irrational branch value cannot be rational. The two odd branch sets
on the source coincide, so the rational point is common.
Make it \(s=\infty\). The deck centers have rational halves.
After translating to those halves the deck involutions become
elliptic negation. The elliptic abscissas \(x_i\) are rational
Möbius functions of the old target parameters, with poles over the
rational branch point. In the source coordinate they have shape
\[
 x_i(s)=\frac{A_i(s)}{L_i(s)^2},
 \qquad \deg A_i=3,\quad \deg L_i\le1,\quad (A_i,L_i)=1.
 \tag{2}
\]
The denominator is constant if the only pole has order three,
and the square of a linear form if there is an additional double
pole. In translated target coordinates the common field is
given by \(y_i^2=x_i^3+b_i^2\).

The polynomials \(A_i^3+b_i^2L_i^6\) have exactly the same nine
simple roots, coming from the three unramified branch fibres.
Their ratio is a rational square because the fields agree and
the denominators in (2) are sixth powers. By the leading
coefficients this constant square is also a rational cube.
Rescale \(A_2,b_2\) accordingly to get
\[
 A_1^3+b_1^2L_1^6=A_2^3+b_2^2L_2^6
\]
with the leading coefficients of \(A_1,A_2\) equal.
Because the right side of
\(A_1^3-A_2^3=b_2^2L_2^6-b_1^2L_1^6\)
has degree at most six, one has \(A_1-A_2=d\in\mathbb Q\).

If \(d\ne0\), the difference of cubes has strict constant sign
\(\operatorname{sgn}(d)\) on all of \(\mathbb R\).
Two nonproportional nonconstant linear forms on the right give
opposite signs at their two zeros. If just one is nonconstant,
the sign at its zero is opposite to the sign for sufficiently
large \(|s|\). Proportional nonconstant forms have a common
zero, again impossible. If both are constant, the right side
has degree zero instead of six. Thus \(d=0\). The linear
forms are proportional and the abscissas in (2) differ by
rational scaling. The coefficient identity gives exactly the
rational elliptic scaling as well, up to sign of the ordinate.
Consequently the maps differ by a rational elliptic isomorphism
and the prior translations, so they pull back the same marked
cubic norm class up to inversion.

For an exact model in this case take \(Y^2=s^9+36\), with the
degree-three abscissa \(x=s^3\) on \(E_6\).
The two rational halves \((12,-42)\) and \((4,10)\) give profiles
\(k_1=32/343\) and \(k_2=-343/375\), respectively, by the scaling in
[XT1](ABC_WEIGHTED_TANGENT_SHARED_FIELD_TRANSLATIONS.md#translating-the-elliptic-involution).
Their parameters are
\[
 t_1=\frac{-100(s^3-12)}{343(3s^3+13)},\qquad
 t_2=\frac{196(s^3-4)}{375(3s^3+13)}.
\]
At \(s=3\) both sets of directions can be scaled positive.
The common norm class is represented by \(Y+6\); it fails the
virtual-unit condition in that specialization, since the two
valuations above \(3\) are \(1,8\).
The [Sage checker](verify_cubic_parameter_tangent_obstruction.py)
replays the two maps and this specialization.

## The branch incidence restriction in genus five

**Theorem CP3.** Suppose two rational cubic parameter maps have a
common genus-five cover, and neither is a rational Möbius
transformation of the other on their targets. Their target branch
sets have the same splitting field. Both branch sets have type
\(S_4\), \(D_4\), or \(1+3\) with group \(S_3\), of the same type
as each other. Count the \(12\) source branch points by their two
target branch values, obtaining a \(4\times4\) multiplicity matrix.
Up to equivariant identifications of the labels, its possibilities are:

* For \(S_4\), the matrix has ones exactly off diagonal.
* For \(D_4\), the matrix has ones exactly off diagonal in one of the
  two equivariant matchings, or has entries constant on the diagonal
  and opposite pairs with multiplicities \(1,2\) in either order
  and zero on adjacent pairs.
* For \(1+3\) with group \(S_3\), the rational pair has zero and
  the pairs joining a rational and an irrational branch have one.
  The irrational block either has ones exactly off its diagonal
  or has twos exactly on its diagonal.

All branch fibres of each cubic map are unramified. The rows and
columns of the matrix sum to three. A matrix entry cannot be three.
Indeed, that would be a common full fibre above a pair of branch
points. Over the algebraic closure translate each selected branch
point to the origin on its elliptic model, so the new deck
involutions are negation. The new abscissas \(x_1,x_2\) on the
source have the same three simple poles. Then \(x_i^3+b_i^2\)
have the same nine simple zeros and the same poles of order
three. They differ by a nonzero constant.
An identity \(x_1^3-c x_2^3=e\) with \(c,e\ne0\) would give a
nonconstant rational map to a smooth plane cubic. Thus \(e=0\)
and the abscissas are proportional. This implies a Möbius relation
on the old targets, which descends to \(\mathbb Q\) by uniqueness.

Let \(M_i\) be the branch splitting field for target \(i\). In the
joint Galois group the subgroup fixing all the first branch
values induces the normal subgroup
\(\operatorname{Gal}(M_2/M_1\cap M_2)\) of the second group.
Every nontrivial normal subgroup of \(S_4\) or \(D_4\) contains
a fixed-point-free involution on the four branches. For \(S_4\)
it contains the Klein four of double transpositions. For \(D_4\)
it contains the central half-turn. Such a subgroup would pair
the columns within each row, forcing an even row sum.
Thus a transitive branch group forces its splitting field to be
contained in the other one. The \(1+3\) group has order six,
so it cannot be paired with a transitive group of order eight
or twenty-four. When both groups are transitive, their splitting
fields therefore coincide.

In the \(1+3\) case, if the two degree-six splitting fields were
different, each relative subgroup would contain a 3-cycle
acting on its three irrational branches. The matrix would have
constant entries \(D\) on the irrational \(3\times3\) block,
constant \(B,C\) on the two rational-to-irrational blocks,
and \(A\) on the rational pair. The row and column constraints
give \(B=C\), \(A+3B=3\) and \(B+3D=3\), hence \(A=3\).
This contradicts the full-fibre exclusion, proving equality
of the splitting fields also in this case.

When both group actions are \(S_4\), their equivalence follows
because every automorphism of \(S_4\) is inner.
The joint orbits are the diagonal and off diagonal.
Their entry values \(A,B\) satisfy \(A+3B=3\) and \(A<3\),
giving \(A=0,B=1\).
For \(D_4\) both profiles are symmetric, so their branch actions
are identified by the normalization of FC1. The orbit sizes in
a row are \(1,1,2\) for diagonal, opposite and adjacent. The
constraint \(A+B+2C=3\), with \(A,B<3\), gives the stated list.
Finally the \(1+3\) actions agree by innerness of automorphisms
of \(S_3\). The entries joining the rational branch to the
irrational ones are both one and the rational pair has zero,
since the other choice would give a full fibre.
The irrational entries satisfy \(D+2E=2\) for diagonal value
\(D\) and off-diagonal value \(E\). This yields the two options.

## Parity of the product-curve singularities

**Theorem CP4.** In the situation of CP3, give the off-diagonal
homomorphism \(\pi_{1*}\pi_2^*\) the value \(n=0\) for nonisomorphic
elliptic targets and its integer value \(n\) after identifying
isomorphic targets.
In the \(1+3\) case with a doubled irrational diagonal one has
\(n=0\). In every other branch pattern of CP3 the elliptic targets
are rationally isomorphic and \(n=\pm1\).
Neither of the doubled \(D_4\) patterns can occur.
In the other cases with \(n=\pm1\) the two generic cubic norm classes are
independent, even over the algebraic closure.

The two cubic parameter functions generate \(\mathbb Q(s)\)
because the degree of their common subfield divides three, and
degree three would identify their individual parameter fields.
Therefore \(C\) maps birationally to its image \(D\) in the
product of the elliptic curves, since the coordinates of either
elliptic projection then supply the quadratic generator.
The elliptic pullbacks on \(\operatorname{Jac}C\) are injective
by the prime-degree argument in CP1. Their images are distinct:
if the images agreed, the two maps would differ by a rational
elliptic isomorphism and translation and hence give a projective
relation on the parameters. The positivity of the pullback
polarization on their product yields \(n^2<9\).

The divisor \(D\) has degree three on each coordinate elliptic
curve. The off-diagonal map of its polarization relative to the
product polarization has value \(\pm n\), with the sign depending
on the identification with the dual. Indeed the difference of
its restrictions to horizontal fibres over \(Q,Q'\) is the
degree-zero divisor class \(\pi_{1*}\pi_2^*(Q-Q')\) on the first
factor, up to the sign of the translation convention. Thus
\[
 \frac{D^2}{2}=9-n^2.
\]
For distinct isomorphism classes the off-diagonal map is zero
and the divisor class splits. For isomorphic factors one can
check the formula by pulling back under the addition and
subtraction isogeny \((P,Q)\mapsto(P+Q,P-Q)\) of degree four.
The pulled-back polarization is diagonal of degrees
\(6+2n,6-2n\), up to exchanging the factors.
Adjunction on the abelian surface and normalization now give
\[
 \sum_{R\in D(\overline{\mathbb Q})}\delta_R
 = p_a(D)-g(C)=5-n^2.
 \tag{3}
\]

The deck involution on \(C\) acts on \(D\) by the product of
the two target deck involutions. Singularities of \(D\) away
from its fixed corners occur in pairs of equal delta invariant.
At each fixed corner, the branches are exactly the preimages
counted by its multiplicity entry in CP3. To see smoothness
of each individual branch, the two projections on \(C\) are
unramified there: their induced parameter maps are unramified
at the corresponding branch values of the double covers.
In local anti-invariant coordinates \(x,y\) on the two target
elliptic curves, every branch is a graph \(y=h(x)\) with \(h\)
an odd power series. If the corner has two branches, their
intersection multiplicity is the order of the difference of
their odd power series, and is odd. A one-branch corner is
smooth. Thus the parity of (3) is the number of doubled
entries modulo two.

The doubled irrational diagonal has three entries of multiplicity
two. It requires \(n\) even and (3) at least three, hence \(n=0\).
All other patterns have zero or four doubled entries, so \(n\)
is odd and therefore \(n=\pm1\).
For these odd cases the pullback polarization matrix on the two
elliptic factors is \(\left(\begin{smallmatrix}3&n\\n&3\end{smallmatrix}\right)\).
Its determinant is eight. Hence the sum map of the pullbacks
has kernel of order prime to three. The two marked rational
order-three divisor classes therefore remain independent,
which gives geometric independence of the functions modulo cubes.

For the doubled \(D_4\) possibilities, normalize both symmetric
profiles by FC1, so their targets are the same \(E_1\) and their
deck centers are \(T=(-1,0)\).
The doubled entries would occupy either the diagonal or the
opposite matching. Over these four pairs the elliptic projections
satisfy \(P_2=P_1+c\), with one fixed \(c\in\{O,T\}\), because
the opposite parameter action is translation by \(T\).
At all eight distinct Weierstrass points over those pairs, both
\(P_1-P_2=-c\) and \(P_1+P_2=T+c\) are constant.
The elliptic maps \(\pi_1\pm\pi_2\) have degrees \(6\pm2n\),
by pullback and pushforward on Jacobians. One has degree four
since \(n=\pm1\), so it cannot have eight points in one fibre.
This excludes the doubled \(D_4\) patterns.

## Cyclic cubic projections

**Proposition CP5.** Under the assumptions of CP3, the two cubic
parameter maps cannot both be Galois over the algebraic closure.
In fact, if one cubic map to a \(j=0\) elliptic double cover is
geometrically Galois with twelve distinct pulled-back branch points,
the automorphism group of their source branch set is cyclic or
dihedral.

A geometrically Galois cubic map has a cyclic group of Möbius
deck transformations of order three preserving the twelve-point
branch set. Let \(G\) be the full group of Möbius transformations
preserving that set over the algebraic closure. It is finite,
because it acts faithfully on the set.
By the [classification over the algebraic closure](sources/beauville-finite-pgl2-subgroups.pdf),
a group other than the cyclic and dihedral groups is \(A_4,S_4\), or \(A_5\).
In each of the three cases the cyclic group of order three is
contained in a tetrahedral \(A_4\). For \(S_4\) it is in the
alternating subgroup, and in \(A_5\) it fixes two letters,
either of whose stabilizers works.

There is one projective conjugacy class of tetrahedral subgroups,
and their order-three subgroups are conjugate inside them.
For the tetrahedral action the only orbits of length below twelve
have lengths \(4,4,6\). This follows from the fixed axes of its four
order-three subgroups and its three order-two subgroups.
The full twelve-point branch set must therefore be a single free
orbit. In the normalization of
[Gow–McGuire, Section 4, Example 2](sources/gow-mcguire-tetrahedral-invariants.pdf),
the cyclic subgroup acts by \(s\mapsto\zeta_3s\) and the other
generator is \(M(s)=(s+2)/(s-1)\). The degree-twelve invariant
has the expression
\[
 h(s)=s^3+\sum_{i=0}^{2}M(\zeta_3^i s)^3
     =\frac{s^{12}+264s^6+440s^3+24}{(s^3-1)^3}.
\]
The first expression is the trace of \(s^3\) from the cyclic
fixed field to the tetrahedral fixed field. The free orbit is a
finite level \(h(s)=\lambda\), since the pole fibre is not free.
Its image under the cubic map \(s\mapsto s^3\), after writing
\(u=s^3-1\), has quartic polynomial
\[
 u^4+(4-\lambda)u^3+270u^2+972u+729.
\]
The invariant \(I\) of this binary quartic is
\(2916(\lambda+24)\).
A nonsingular quartic would have \(j=0\) only if \(I=0\).
At \(\lambda=-24\) it is \((u+1)(u+9)^3\), so it cannot supply
four distinct branches. This rules out the noncyclic,
nondihedral groups. A cyclic or dihedral group has at most
one subgroup of order three. If both cubic maps were Galois,
they would have the same invariant field and be projectively
related, contrary to CP3.

## The simple branch patterns

**Theorem CP6.** No nonprojective pair in CP3 can have a branch
matrix with twelve distinct off-diagonal matches. Hence the only
remaining genus-five case has two asymmetric profiles with
rational branch points, zero off-diagonal elliptic Hom pairing,
and the doubled irrational diagonal of CP3.
In particular two rational cubic parameter maps to symmetric
profiles cannot give independent generic marked norm classes.

First, consider two equianharmonic branch quartets, meaning their
elliptic double covers have \(j=0\), with a projective matching
between the quartets. The following argument is over the
algebraic closure. Suppose a rational bidegree-\((3,3)\) curve
in the two parameter lines contains the twelve ordered pairs
off the matching and omits all four matched pairs.
Use the projectivity of the matching to identify the two lines.
Normalize the common quartet to the roots \(\alpha_i\) of
\(f(t)=t^4-2t^2-1/3\) in FC1.
This is possible because two \(j=0\) elliptic curves are
geometrically isomorphic, and translation by a half of the
deck center identifies each cover with negation on its curve.

Evaluation at the four branch points identifies each space of
binary cubics in a target variable with its four coordinate
values. Interpolating in both variables gives the equation
\[
 \sum_{i=1}^4 c_i \ell_i(t_1)\ell_i(t_2)=0,\qquad
 \ell_i(t)=f(t)/(t-\alpha_i),\qquad c_i\ne0.
 \tag{4}
\]
In particular the curve is symmetric. Its quotient by interchange
of the coordinates is the irreducible plane cubic \(P=0\) in
\(p=t_1+t_2,q=t_1t_2\) and their projective homogenization.
The original curve double covers this plane cubic, with field
given by adjoining \(\sqrt{p^2-4q}\). It follows that the
normalization of the plane cubic is rational and that the
discriminant conic \(p^2-4q=0\) has exactly two intersections
of odd multiplicity on that normalization.
The cubic has at most one singular point, with at most two
branches there. Thus its restriction \(S(t)=P(2t,t^2)\) on
the diagonal conic has at most two distinct odd-multiplicity
roots as a binary sextic. If all multiplicities are even,
at least one singular point of the plane cubic must lie
on the conic to account for the ramification upstairs.

Write \(S(t)=H(t)^2Q(t)\), where \(H,Q\) are binary quadratics.
Such a factorization exists geometrically when the sextic has
two odd roots, by separating their squarefree part. It also exists
when the sextic is a square, by factoring one linear form from
its cubic square root. In either case \(H\) has no common root
with \(f\), because (4) omits the four matched pairs.

Let \(S_i\) denote the coefficient of \(t^i\) in this sextic.
The four-dimensional space obtained by restricting (4) to the
conic has the three defining relations
\[
 S_0-\tfrac19S_4+\tfrac29S_6=0,\qquad
 S_1+S_3-\tfrac13S_5=0,\qquad
 S_2+\tfrac23S_4-\tfrac73S_6=0.
 \tag{5}
\]
They follow by taking the four trace weights \(1,\alpha,\alpha^2,
\alpha^3\) in \(f(t)/(t-\alpha)\), with \(\alpha\) a root of \(f\).
For \(H(t)=at^2+bt+c\), the determinant of the three linear
equations (5) for the coefficients of \(Q\) is
\[
 -\tfrac13\,(a^2+b^2+2ac-3c^2)\operatorname{Res}(f,H).
 \tag{6}
\]
The resultant is nonzero. Therefore the quadratic \(H\), up to
scale, lies on the conic in (6). Its complete parametrization
is
\[
 H_r(t)=(1-3r^2)t^2+4rt+(1+r^2),\qquad
 H_\infty(t)=1-3t^2.
\]
Along this conic the linear system (5) has rank two and its
kernel is generated by \(Q=(t-r)^2\), with the homogeneous
square at infinity for \(r=\infty\). For example, substituting
the finite parametrization gives the indicated kernel vector
and the \(2\times2\) minors have polynomial gcd one.
The coefficient matrix at infinity has rank two directly.
Consequently \(S\) is a square, so a singular point must lie
on the diagonal conic.

Normalize the plane cubic for finite \(r\) by
\(P_r(2t,t^2)=3H_r(t)^2(t-r)^2\).
Put \(f_+(r)=r^4+2r^2-1/3\).
The same interpolation basis gives the identities
\[
\begin{aligned}
 P_{r,q}(2r,r^2)&=-9f(r)^2,\\
 \operatorname{Res}_t(P_{r,q}(2t,t^2),H_r(t))&=
 11664(r^2-1/3)^2 f(r)^2 f_+(r)^2,\\
 \operatorname{Res}_t(f(t),H_r(t))&=64 f(r)f_+(r).
\end{aligned}
\tag{7}
\]
If \(r^2\ne1/3\), nonvanishing of \(S\) at the four branch points
rules out both \(f(r)=0\) and \(f_+(r)=0\) by the last identity.
The other two identities then show that no point of intersection
with the conic can be singular. There is no intersection at
infinity in this case. For \(r^2=1/3\), the cubic is proportional
to
\[
 3p^2q-3q^2-6q+1
\]
in affine coordinates and its points on the conic are smooth:
at \(t=\pm r\) the \(q\)-derivative is nonzero and the point
at infinity has a nonzero derivative in the homogenizing
coordinate. At \(r=\infty\) the cubic is proportional to
the same expression, so all its points on the conic are smooth.
These facts contradict the required singularity on the diagonal.
The [exact checker](verify_cubic_parameter_off_diagonal_obstruction.py)
computes (5)–(7) directly in the interpolation basis.

It remains to verify the projective matching for the simple
patterns of CP3. By CP4 the two elliptic models are rationally
isomorphic in all these cases. Identify them as one elliptic
curve. The branch points are the four halves of the respective
deck centers. For \(S_4\) and \(1+3\), the equivariant matching of
the two four-point sets is affine for their two-torsion actions,
since every permutation of four points is affine over
\(\mathbb F_2\). Its linear part commutes with the \(S_3\) action
on \(E[2]\) and hence is the identity. It is therefore translation
by a rational point that sends one deck to the other, and gives
a projective parameter matching. For \(D_4\), the matching is
the identity or the opposite map in the symmetric FC1 coordinates.
Both are projective, since the opposite map is \(t\mapsto-t\).
Thus (4) rules out all the simple patterns, leaving only the
doubled irrational pattern stated in the theorem.

For a projectively related pair of cubic maps, the equality of
squareclasses descends to the target parameter field, because the
source extension has odd degree three. The elliptic covers then
identify, and their marked classes agree up to inversion by the
unique rational cyclic order-three subgroup. This also proves
the symmetric conclusion when the maps are projectively related.

## Remaining cubic construction

Independent generic directions from two rational cubic parameter
maps require genus five and the doubled asymmetric matching of
CP3–CP6. The two elliptic targets each have a rational branch
point and their projections have zero off-diagonal pairing.
At least one of the parameter maps must be geometrically non-Galois.
The [doubled branch classification](ABC_CUBIC_TANGENT_DOUBLED_BRANCH_FAMILY.md)
determines that remaining case. Its inverse-coefficient construction
admits two simultaneous positive directions, but its specialized
norms fail a uniform local virtual-unit condition.
