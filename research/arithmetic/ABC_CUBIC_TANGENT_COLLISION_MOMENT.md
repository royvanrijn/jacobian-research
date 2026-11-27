# A power saving for repeated cubic tangent fields

The field collision problem for bounded prime progressions
has a uniform square-root-scale decomposition. For each pair
of small scales, a collision lies on a smooth quartic surface.
The surface bound is uniform in its coefficients. Large
scales are rare by a binary quartic square-divisor estimate.
For an irreducible form, the square divisors also shrink
the boxes in which the determinant method is applied.
This gives a power saving for the full ordered collision
moment in the unweighted family.

## Square-root scales of primitive quartic values

For a fixed integral binary quartic \(F\) with nonzero
discriminant, use the signed factorization
\[
 F(x,y)=d\,s^2,\qquad d\text{ squarefree},\quad s>0,
\]
at primitive nonzero pairs with \(F(x,y)\ne0\).
Represent projective directions by primitive pairs
modulo overall sign, of height \(\max(|x|,|y|)\).
Let \(B\ge2\) bound that height.

**Lemma FM1.** For \(0<\alpha\le1/3\), the
number of these directions with \(s>B^\alpha\)
is \(O_{F,\alpha,\varepsilon}(B^{2-\alpha+\varepsilon})\).

First, a prime \(q>B^\alpha\) with \(q^2\mid F(x,y)\)
contributes \(O_F(B+B^{2-\alpha})\) in the range
\(q\le CB\) for a fixed sufficiently large \(C\).
Outside a fixed finite set it gives at most four
simple projective roots modulo \(q^2\).
Each root gives a two-dimensional integer
lattice of determinant \(q^2\).
For a determinant-\(m^2\) lattice, the number of
ambient primitive integer vectors in the box is
\(O(1+B^2/m^2)\). If the vectors span only a
line there are at most two primitive ones.
Otherwise the second successive minimum is
\(O(B)\), and the standard two-dimensional
successive-minima count gives \(O(B^2/m^2)\).
Summing over \(q\) gives the first assertion.

For \(q>CB\), a fixed factorization of \(F\)
over \(\mathbb Q\) lets us restrict to one
irreducible factor. For large \(q\) the
other factors are coprime to it on
primitive pairs, by their resultants.
Factors of degree at most two have
absolute value less than \(q^2\) on the box
once \(C\) is large, and the value of \(F\)
is nonzero. For degrees three and four,
[Helfgott, Propositions 4.12–4.13](sources/helfgott-square-free-sieve-0309109.pdf)
bounds the total number by
\(O_F(B^{4/3}(\log B)^{A_F})\).

It remains to allow all factors of \(s\)
to be at most \(B^\alpha\).
Multiply prime factors of \(s\), with
multiplicity, up to the first product
\(m>B^\alpha\). Then \(m\le B^{2\alpha}\)
and \(m^2\mid F(x,y)\).
For every \(m\) there are
\(O_F(4^{\omega(m)})\) possible
projective root lattices of determinant
\(m^2\), by Hensel lifting in the two
unit-coordinate charts and the Chinese
remainder theorem. At each of the
finitely many bad primes the number of
projective roots modulo \(q^j\) is
bounded independently of \(j\).
In either chart the polynomial and its
derivative have a fixed nonzero
resultant or Bézout constant. For
large \(j\), the derivative valuation
on a root modulo \(q^j\) is therefore
bounded. Hensel's lemma gives an
actual local root and confines the
class to a bounded number of residues
above that root. The remaining \(j\)
are bounded in terms of \(F\).
Summing the lattice bound over
\(B^\alpha<m\le B^{2\alpha}\)
now gives
\(O_{F,\varepsilon}(B^{2\alpha+\varepsilon}
+B^{2-\alpha+\varepsilon})\).
All three estimates have the asserted
order in the stated range of \(\alpha\).

## Field collisions away from projective symmetries

Let \(\mathcal G_F\) be the group of
rational projective transformations
preserving the four geometric roots
of \(F\). It is finite because three
values determine a Möbius transformation.
Write \(\mathcal D(B)\) for any set
of primitive projective directions
of height at most \(B\) with \(F\ne0\)
and at most one direction from each
\(\mathcal G_F\)-orbit.

**Theorem FM2.** The number of directions
of \(\mathcal D(B)\) whose signed squarefree
value also occurs at a different direction
of \(\mathcal D(B)\) is
\(O_{F,\varepsilon}(B^{47/27+\varepsilon})\).
The number of ordered distinct collisions
with both square-root scales at most
\(B^{7/27}\) has the same bound.
If \(F\) is irreducible, the incidence
bound improves to \(O_{F,\varepsilon}(B^{39/23+\varepsilon})\),
and the improved ordered-pair bound
holds for both scales at most \(B^{7/23}\).

Restrict both scales to \(M=B^\alpha\).
Cancel their gcd and write the remaining
coprime positive scales as \(a,b\).
The collision equation is the quartic surface
\[
 b^2 F(x,y)=a^2 F(x',y')
 \quad\text{in }\mathbb P^3 .
\]
It is smooth over the algebraic closure.
A singularity would require both binary
gradients to vanish, which for such a
form forces both coordinate pairs to be
zero. The four-tuple of the collision
gives a projective point of height at
most \(B\) and determines the two
directions. The signs of the primitive
representatives contribute at most a
constant multiplicity.

No nontrivial collision under consideration
lies on a rational line on this surface.
On a line, each binary coordinate
projection has rank at most two.
If one has rank two, the restricted
binary quartic has four simple roots.
A rank-one projection on the other
side gives zero or a fourth power,
which is impossible. When both have
rank two, the line is the graph of
an invertible rational map preserving
the branch set, and its directions
belong to one orbit. If both projections
have rank at most one, write them as
fixed rational vectors times linear
forms on the line. If both quartic
values are nonzero, their linear
forms must be proportional and cannot
span a line in projective space.
The remaining case is confined to
\(F(x,y)=F(x',y')=0\) or to a zero
coordinate pair. Those points are
outside the counted directions.

As \(a,b\) are coprime, we have
\(a^2\mid F(x,y)\) and \(b^2\mid F(x',y')\).
Use the root lattices of FM1 on both
sides. There are \(O_F(4^{\omega(ab)})\)
lattice pairs since \(a,b\) are coprime.
Let \(\lambda_1\le\lambda_2\) be the
successive minima on a first-side
lattice, and \(\lambda'_1\le\lambda'_2\)
those on a second-side lattice.
The products are comparable to
\(a^2,b^2\), respectively.
Use reduced bases, which give coordinate
bounds \(O(B/\lambda_i)\) and
\(O(B/\lambda'_i)\).
If either second minimum exceeds a
sufficient constant times \(B\), the
corresponding primitive directions
lie in one line, so there is at most
one direction on that side. The other
side then represents a prescribed
nonzero quartic value, giving \(O_F(B)\)
pairs by fixing its first coordinate
and solving a polynomial of degree
at most four in the other. If a
special fiber is constant the
contribution at that coordinate is
also \(O(B)\), and there are only
a bounded number of such coordinates.

Now assume both second minima are
\(O(B)\). Substituting the lattice
bases is a rational projective
change of variables and preserves
smoothness and the line exclusions.
Take size bounds \(B_1\ge B_2\)
and \(B_3\ge B_4\), each at least
one, comparable to the respective
\(B/\lambda_i\) and \(B/\lambda'_i\).
Write \(V=B_1 B_2 B_3 B_4\)
and \(B_*=\max B_i\).
Then \(V\ll B^4/(a b)^2\).
Let \(T_*\) be the largest monomial
size, in the sense of Salberger's
Notation 1.1, on the transformed
quartic. In the first binary block,
the monomial with size \(B_1^4\)
has nonzero coefficient unless the
first basis vector is a root of
\(F\). In that case the simple
root gives the mixed monomial
of size \(B_1^3B_2\).
The second block gives the same
bound with \(B_3,B_4\).
Taking a geometric mean shows
\(T_*\ge V\).

Apply [Salberger, Theorems 3.16 and 0.11](sources/salberger-counting-rational-points-2023.pdf).
The former counts points outside
\(O((V/T_*^{1/4})^{1/4}\log V+1)\)
curves of bounded degree by
\(O_\varepsilon(B^\varepsilon V^{3/8})\),
since \(T_*\ge V\). Each auxiliary
curve of degree at least three
contributes \(O_\varepsilon(B_*^{2/3+\varepsilon})\)
by the latter theorem. The
nonrational lines have at most
one rational point each. Lines
and conics on a smooth quartic
are uniformly bounded in number,
as in the discussion following
Salberger's Theorem 0.5.
Thus, after excluding the rational
lines, the conics and nonrational
lines contribute \(O_\varepsilon(B_*^{1+\varepsilon})\).
This bounds the points for
each lattice pair by
\[
 O_{F,\varepsilon}\!\left(
 B^{3/2+\varepsilon}(ab)^{-3/4}
 + B^{3/4+\varepsilon}(ab)^{-3/8}B_*^{2/3}
 + B_*^{1+\varepsilon}\right).
\]
Use \(B_*\ll B\) and sum the
lattices and \(a,b\le M\).
Including the degenerate lattice
pairs gives
\[
 \ll_{F,\varepsilon} B^\varepsilon
 \left(B^{3/2}M^{1/2}
       +B^{17/12}M^{5/4}+BM^2\right).
\]
With \(\alpha=7/27\) this is
\(O_{F,\varepsilon}(B^{47/27+\varepsilon})\).

For the irreducible improvement,
every nonzero vector \(v\) in
the first lattice has nonzero
\(F(v)\) divisible by \(a^2\).
By quartic homogeneity its length
is \(\gg_F\sqrt a\). Likewise the
other lattice has minimum
\(\gg_F\sqrt b\), so now
\(B_*\ll_F B/\sqrt{\min(a,b)}\).
In the auxiliary-curve term this
replaces \(M^{5/4}\) by \(M^{11/12}\):
\[
 \sum_{a\le b\le M}a^{-3/8-1/3}b^{-3/8}
                     \ll M^{11/12}.
\]
The other terms may be left as
before. Taking \(\alpha=7/23\)
gives the improved ordered-pair bound.

For the incidence estimate, FM1
at the corresponding \(\alpha\)
bounds the number of large-scale
directions by \(O_{F,\varepsilon}(B^{2-\alpha+\varepsilon})\).
A fibre with at least two small-scale
directions is counted by the
small-scale ordered pairs.
A fibre with exactly one small-scale
direction and some large-scale
directions adds at most one more
small-scale incidence per large-scale
direction. This gives both claimed
incidence bounds.

The estimate uses the equality
of squarefree parts even when the
original quartic values are unequal.
For comparison, the essential
representation bound from
[Stewart–Xiao, Lemma 2.4](https://arxiv.org/pdf/1605.03427)
addresses repetitions of the
original integral values.
Uniformity of Salberger's
surface bound is needed here
because the square-root scales
vary in the summation.

## The unweighted second factorial moment

**Theorem FM3.** Let \(N_d(Z)\) count the
genuine progressions \((p-h,p,p+h)\)
with center \(p\le Z\) and tangent
field \(\mathbb Q(\sqrt d)\) as in CT1.
Then
\[
 \sum_d N_d(Z)(N_d(Z)-1)
          \ll_\varepsilon Z^{39/23+\varepsilon}.
\]
In any fixed admissible prime rectangle
of CT4, with either sign of \(R\) and
optionally prime \(h/6\), the number of
distinct fields is the number of
triangles up to an additive
\(O_\varepsilon(Z^{39/23+\varepsilon})\).
The conclusion in that rectangle also
holds with the count restricted to
either sign individually.

Use \(F=-R\).
By [FC4](ABC_CUBIC_TANGENT_FIELD_COLLISIONS.md#fixed-projective-relations-and-adjacent-centers),
the rational projective branch symmetries
are \(t\mapsto t\) and \(t\mapsto-t\).
The ratios \(t=p/h>1\) give a set with
one representative per orbit, and
\(p,h\) are coprime. Irreducibility
follows by Eisenstein at \(3\)
in the reciprocal quartic.
FM2 gives the improved bound for
ordered pairs with both scales at
most \(Z^{7/23}\). By
[FC5](ABC_CUBIC_TANGENT_FIELD_COLLISIONS.md#multiplicity-uniform-in-the-field),
\(\max_d N_d(Z)\ll_\varepsilon Z^\varepsilon\).
Pairs with at least one large-scale
direction are at most twice the count
in FM1 times this maximum. Taking
a smaller \(\varepsilon\) in these
inputs gives the displayed factorial
moment. Removing one triangle from
each nonempty field fibre leaves a
count bounded by the factorial
moment. Restricting to a rectangle
can only decrease the collision count.

For a fixed rectangle with the signs
and prime forms used in
[CT4](ABC_CUBIC_TANGENT_NORM_RANK_AND_CONDUCTORS.md#squarefree-dyadic-fields-and-rank-two-conductors),
the total number is asymptotic to
the positive prime-form constant
times \(Z^2/(\log Z)^a\), where
\(a=3\) or \(4\).
Thus the same asymptotic holds for
the number of distinct fields, with
no condition of quartic squarefreeness.
The field discriminants have absolute
value at most \(C Z^4\).
The count up to that upper limit is
distinct from a prescribed dyadic
window for the squarefree part.

As a bounded test, the
[collision census](verify_cubic_tangent_collision_census.py)
finds \(1365099\) triangles and no
repeated signed squarefree values
for centers at most \(45000\).
This computation checks full
squarefree parts on that interval.
The asymptotic does not settle
whether an unweighted collision
exists at a larger center.

## Asymmetric weighted profiles

**Corollary FM4.** For each fixed asymmetric
weighted tangent profile in WM1, among
primitive prime triangles with all
vertices bounded by \(Z\), at most
\(O_{k,\varepsilon}(Z^{47/27+\varepsilon})\)
are incident to a different triangle on
the same profile with the same quadratic
field. If the weighted quartic is
irreducible the bound improves to
\(O_{k,\varepsilon}(Z^{39/23+\varepsilon})\).
The constant can depend on the
fixed rational changes of coordinates.

The branch points of the weighted
quartic correspond to the halves
of the involution center \(T\) on
\(E_w\), hence form a torsor under
\(E_w[2]\). On an asymmetric profile
\(4w\) is not a rational cube, as
proved in
[WM7](ABC_WEIGHTED_CUBIC_TANGENT_ELLIPTIC_MARKS.md#the-entire-cyclotomic-trace-coset-in-rank-one).
The cubic \(X^3+(4w)^2\) is then
irreducible and has negative
nonsquare discriminant. Its two-torsion
action is all of \(S_3\).
On the four branch points the
Galois group is a subgroup of
\(V_4\rtimes S_3=S_4\) surjecting
onto \(S_3\). Its kernel in the
translation plane is invariant
under \(S_3\) and is therefore
trivial or the whole plane.
The two possibilities are a
point stabilizer \(S_3\) or \(S_4\).
Each has trivial centralizer in
\(S_4\), so the rational projective
branch symmetry group is trivial.
The bounds on the prime vertices
give height \(O_k(Z)\) in the
fixed parameter coordinates.
Each positive primitive triple is
uniquely specified by its rational
direction. Apply FM2 to the associated
integral binary quartic, after
clearing denominators by a square
as necessary. This gives the
claimed incidence bound, including
the fixed profile with the two
explicit cyclotomic triangles.
