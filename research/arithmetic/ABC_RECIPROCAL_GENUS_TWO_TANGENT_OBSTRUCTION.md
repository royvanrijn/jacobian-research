# A reciprocal genus-two cover with incompatible tangent witnesses

The two degree-two projections of a reciprocal
genus-two curve can carry independent generic
cubic norm classes. The same construction imposes
two restrictions on their specializations.
The positive Heron cones of the two projections
have disjoint preimages over the reals.
At any odd prime at which the defining constant
has noncubic valuation, one of the two specialized
norm classes fails the virtual-unit test.

## Which weighted profiles have a rational branch point

Keep \(f_k(t)=-3t^4+4(2k+1)t^3-6k(k+1)t^2+k^2(k+1)^2\)
and the cover of
[WM1](ABC_WEIGHTED_CUBIC_TANGENT_ELLIPTIC_MARKS.md#the-double-cover-and-its-three-marks).
A rational branch point is a rational half
of its deck center \(T\).

**Proposition GT1.** For \(k\ne0,-1\) the
quartic \(f_k\) has a rational root exactly
when
\[
 k=-\frac{2r+1}{r^3(r+2)},\qquad
 r\in\mathbb Q\setminus\{0,-1,-2,-1/2,1\}.
\]
Then a branch parameter is
\(t_0=-(1+r)^2/(r^3(r+2))\).
There are \(O(H^{1/2})\) such
profiles of rational projective
height at most \(H\).

At a finite root of the quartic
the Heron numerator is zero.
None of the three weighted
inputs \(A_i=a_i^2 X_i^3\)
vanishes, since the three
marked fibres have nonzero
square values in WM1.
The identity
\(4A_1 A_3=(A_1+A_3-A_2)^2\)
requires \(A_1,A_3\) to
have the same sign and
squareclass. Writing their
square roots with signs
makes \(A_2\) have the same
sign and squareclass as
well. Hence the \(X_i\)
are all a common nonzero
rational number times
\(z_i^2\), where the \(z_i\)
are nonzero rationals and
their signs may be chosen
to satisfy
\(\sum a_i z_i^3=0\).
The prime line and tangent
conditions also give
\(\sum a_i z_i^2=\sum a_i=0\).
The \(z_i^2\) are distinct,
since equality of any two
would force equality of all
three and give the excluded
point \(h=0\).
The determinant of the rows
\((1),(z_i^2),(z_i^3)\)
is the Vandermonde times
\(z_1z_2+z_1z_3+z_2z_3\).
It therefore forces
that last sum to vanish.
Normalize \(z_1=1,z_2=r\).
Then \(z_3=-r/(1+r)\).
The relations give
\[
 \frac{a_2}{a_1}
    =\frac{1-z_3^2}{z_3^2-z_2^2}
    =-\frac{2r+1}{r^3(r+2)},\qquad
 t_0=\frac{k}{1-z_3^2}.
\]
Conversely the displayed
expressions give distinct
nonzero squares and the
three equations, hence
a root. Its exceptional
\(r\)'s are exactly those
excluded above.
The rational map \(r\mapsto k\)
has degree four and no base
point after homogenization.
Its output height is
comparable to the fourth
power of its input height,
so the count follows.

## Independent norm functions on the reciprocal cover

Put \(E_b:y^2=x^3+b^2\) and take
\(b\in\mathbb Q^\times\).
Consider the genus-two curve
\[
 C_b:\quad Y^2=s^6+b^2
\]
with the two degree-two maps
\[
 \pi_1(s,Y)=(s^2,Y)\in E_b,\qquad
 \pi_2(s,Y)=(b^2/s^2,b^2Y/s^3)\in E_{b^2}.
\]
Both commute with the
hyperelliptic involution
and elliptic negation.

**Proposition GT2.** Suppose the two weighted
profiles have \(4k_i(k_i+1)\) equal to
\(b\) and \(b^2\), respectively, up to
rational cubes. Suppose also that
their deck centers have rational
halves. Then the two projections,
scaled to the profile models and
translated by the rational halves,
give a common hyperelliptic field
parameter with two independent
cubic norm classes. These classes
are represented on \(C_b\) by
\[
 Y+b,\qquad b^2(Y/s^3+1).
 \tag{1}
\]
Each has norm a cube in
\(\mathbb Q(s)\) and divisor
divisible by three.

The two maps are the explicit
reciprocal form of an elliptic
gluing in degree two; compare
[Howe–Leprévost–Poonen, Proposition 4](sources/howe-leprevost-poonen-split-jacobians.pdf).
Here their formulas can be
checked directly.
The target \(E_{b_i}\) of
\(\pi_i\), for \(b_1=b,b_2=b^2\),
is isomorphic to the profile
curve by rational scaling.
Composing with the translation
by half its deck center sends
negation to the required deck
involution. Each corresponding
edge norm class has divisor
class \(-D_i\) after dividing
its divisor by three, as in
[XT1](ABC_WEIGHTED_TANGENT_SHARED_FIELD_TRANSLATIONS.md#translating-the-elliptic-involution).
This agrees with the class
of \(y+b_i\) on \(E_{b_i}\).
Translation is trivial on
degree-zero divisor classes.
Thus the functions differ by
a cube times a rational
constant on each elliptic
factor. Both norms to the
parameter line are cubes,
so the extra constant is a
rational cube. Pullback
gives (1).

The pullback of the order-three
class of \(-D_i\) is nonzero
because pushforward after
pullback is multiplication
by two. The involution
\((s,Y)\mapsto(-s,Y)\) fixes
the first pullback and
negates the second, by the
two map formulas.
Their divisor classes are
therefore independent even
over the algebraic closure.
This implies independence
of their function classes
modulo cubes.

For an exact example take
\(b=6\), \(r_1=-7/2\)
and \(r_2=-17/7\) in GT1.
The profiles are
\(k_1=32/343\) and
\(k_2=3087/4913\).
Their coefficients
\(4k_i(k_i+1)\)
are \(6(20/49)^3\)
and \(36(140/289)^3\),
respectively.
The resulting degree-two
parameter maps are
\[
 t_1(s)=
 \frac{-100(s^2-12)}{343(3s^2+13)},\qquad
 t_2(s)=
 \frac{176400(2s^2-1)}{4913(73s^2+108)} .
\]
For instance at \(s=1\)
the classes (1) are
independent in \(\mathbb Q(\sqrt{37})\).
The second class has
a nontrivial valuation
modulo three.

## Real positivity and the virtual-unit test

**Theorem GT3.** In the reciprocal
construction of GT2, no real
parameter gives both sets of
three positive prime directions,
even before primality is imposed.
In addition, at any odd prime \(q\)
with \(3\nmid v_q(b)\), the
specializations of (1) at rational
finite \(s\ne0\) cannot both have
all their \(\mathfrak q\)-valuations
divisible by three.
Exactly one of the two classes
fails that local condition.

For the real assertion, pull
one weighted parameter map back
to an elliptic factor \(E_{b_i}\).
Let \(V_i\) be the rational
half of its deck center in
these coordinates.
After translation it is a
degree-one rational function
of \(x(Q)\), for \(P=V_i+Q\).
It has no pole on \(x(Q)\ge0\).
Indeed, by the quartic
formula in WM1, the two
points above its parameter
infinity are the scaled
points \(-[\zeta]T\) and
\(-[\zeta^2]T\) on the
original profile.
Subtracting the half
of the center sends them
to \(\pm J V_i\), where
\(J=[\zeta]-[\zeta^2]\).
Writing \(x=x(V_i)\), the
addition law gives
\[
 x(JV_i)=-\frac{x^3+4b_i^2}{3x^2}<0 .
\]
The half is nontorsion
because the symmetric
profiles have no rational
branch point by FC1
and permutation of the
marks. Thus its abscissa
is nonzero.

The three marked roots
in the abscissa \(x(Q)\)
are the three abscissas
of \(D_j-V_i\), with
\(D_j\in\{O,(0,b_i),(0,-b_i)\}\).
For \(V_i=(x,y)\)
these abscissas are
\[
 x,\quad
 \frac{2b_i(b_i+y)}{x^2},\quad
 \frac{2b_i(b_i-y)}{x^2}.
\]
Their product is \(-4b_i^2\).
Exactly one is negative:
this follows from
\(y^2=x^3+b_i^2\)
and the signs of
\(b_i(b_i+y)\) and
\(b_i(b_i-y)\).
They are distinct since
equality up to elliptic
negation of two displayed
points would give
\(2V_i\) torsion.
Any negative abscissa
on \(E_{b_i}(\mathbb R)\)
has magnitude at most
\(|b_i|^{2/3}\).
Hence the larger
positive marked root
\(M_i\) is at least
\(2|b_i|^{2/3}\).

At \(x(Q)=+\infty\)
the profile parameter
is the branch value of
\(V_i\), where the
three direction
coordinates have
the same sign by
the Heron equality
argument in GT1.
On the positive
abscissa ray there
are two distinct
simple marked roots
and no pole.
Crossing either
changes the sign of
one of the three
homogeneous direction
coordinates. After
crossing one or both
the signs are mixed.
Thus all three
directions can be
made positive by
a common scaling
only if \(x(Q)>M_i\).

For the two projections
one has \(x(Q_1)=s^2\)
and \(x(Q_2)=b^2/s^2\).
Their product is \(b^2\).
But \(M_1 M_2\ge4b^2\),
so the two positivity
conditions are incompatible.
At \(s=0,\infty\) the two abscissas are \(0,\infty\).
The factor with abscissa zero still has mixed directions, and the other
is a branch value with Heron value zero.

For the valuation assertion,
put \(a=v_q(b)\) and
\(e=v_q(s)\). Since
\(3\nmid a\), the integers
\(a\) and \(3e\) differ.
As \(q\) is odd,
\(s^6+b^2\) is a
square over \(\mathbb Q_q\):
factor out the dominant
square and take a square
root of the remaining
principal unit.
If \(a<3e\), on the
branch where \(Y/b\)
reduces to \(1\) the first
class of (1) has valuation
\(a\). Its conjugate
valuation is the
negative modulo three
by the cube norm.
The second class has
valuation \(3(a-e)\)
on both branches,
since the \(Y/s^3\)
term dominates \(1\).
If \(a>3e\), the first
class has valuation \(3e\)
and the second has
valuation \(2a\)
on the branch where
\(Y/s^3\) reduces
to \(1\).
This proves the
complementary virtual-unit
condition modulo
three. In particular
the \(b=6\) pair has
this obstruction at \(3\)
for every rational
nonzero finite parameter.

[The Sage checker](verify_reciprocal_genus_two_tangent_obstruction.py)
replays the two rational
halves, their parameter
maps and the finite
norm-class fixture.
The real and local
obstructions apply
throughout the stated
family of reciprocal
degree-two maps.

**Corollary GT4.** In the example \(b=6\)
there are infinitely many distinct real
quadratic fields in which (1) gives
two independent norm-cube classes.
At most one dimension of their
span consists of virtual units.

Adjoin \(\zeta_3\) and cube roots
of both functions (1) to the generic
quadratic field over \(\mathbb Q(s)\).
GT2 gives degree nine for the Kummer
extension after the quadratic
and cyclotomic adjunctions.
The full extension is Galois over \(\mathbb Q(s)\), since conjugating
\(Y\) inverts both classes modulo cubes and the cube roots of unity
are included.
[Hilbert specialization](sources/zywina-quantitative-hit-2010.pdf)
preserves the full specialization
degree over \(\mathbb Q\) for
infinitely many rational parameters.
For each fixed signed squarefree
\(d\), its corresponding
parameters come from the smooth
genus-two curve \(dZ^2=s^6+36\).
By [Faltings, Satz 7](sources/faltings-1983-endlichkeitssaetze.pdf)
there are only finitely
many such parameters in
one field. The fields
therefore vary infinitely.
The valuation homomorphism
modulo three from their
two-dimensional span is
nonzero at \(3\) by GT3.
Its kernel contains the
virtual-unit subspace.

## Classification with two degree-two projections

**Theorem GT5.** Let \(C/\mathbb Q\) be a smooth geometrically integral
genus-two curve. Suppose \(\pi_i:C\to E_{B_i}\), for \(i=1,2\), are
degree-two maps to the weighted tangent curves of WM1, with
\(B_i=4k_i(k_i+1)\), intertwining the hyperelliptic involution with
their weighted deck involutions. If the two pullbacks have the same
elliptic image in \(\operatorname{Jac}C\), the two generic norm classes
agree up to inversion. If the images are distinct, both profiles are asymmetric, \(B_1B_2\) is a rational cube,
and \(C\) with its two maps is the reciprocal construction of GT2
up to rational scaling of the elliptic models and translation by the
two rational halves. In the distinct-image case no real point of \(C\) projects to
two positive tangent directions.

Each \(\pi_i^*\) is injective. Its kernel is killed by two from
\(\pi_{i*}\pi_i^*=[2]\). A nonzero point in the kernel over the
algebraic closure would give an étale double cover of \(E_{B_i}\)
through which \(\pi_i\) factors, by triviality of the pulled-back
2-torsion line bundle. Since \(\pi_i\) has degree two, this would
identify \(C\) with a genus-one curve.

For the same-image case, the injective pullbacks give a rational
isomorphism of the elliptic curves and the pushforwards agree after
that identification. The original maps differ by the isomorphism and
an elliptic translation, so the pullbacks of the order-three divisor
classes \(\pm D_i\) agree up to sign. Indeed, the only rational
cyclic subgroup of order three is the marked subgroup, as shown
below. As in GT2, functions with this property agree up to inversion,
a cube and a rational constant. Their norms to the hyperelliptic base
are cubes, so the constant is a cube as well.

For distinct images, write \(\lambda=(\pi_1^*,\pi_2^*)\) for the sum map from the product
to \(\operatorname{Jac}C\). The distinct images make this an isogeny.
With the canonical principal polarizations, its pullback polarization
has diagonal entries \([2]\) and off-diagonal entries \(\pi_{i*}\pi_j^*\).

Two curves \(E_b\) and \(E_c\) of the indicated form are
\(\mathbb Q\)-isogenous only if \(c/b\) is a rational cube, in which
case they are \(\mathbb Q\)-isomorphic. Indeed, a scaling using a
cube root of \(c/b\) identifies the curves over the algebraic closure,
and the endomorphism algebra after this identification is
\(\mathbb Q(\zeta_3)\). Over \(\mathbb Q(\zeta_3)\) the
endomorphisms are all defined, while Galois acts on the scaling by the
cubic Kummer character of \(c/b\). A nontrivial character acts on the
space of homomorphisms by multiplication by a nontrivial cube root
of unity and has no nonzero fixed homomorphism. A rational number
that becomes a cube over the quadratic cyclotomic field was already
a cube by taking norms. Over \(\mathbb Q\) the endomorphism ring of
each factor is \(\mathbb Z\), since conjugation interchanges the two
nontrivial automorphisms of order three.

If the factors are not isogenous over \(\mathbb Q\), the off-diagonal
entries vanish. Otherwise identify both factors with one \(E_b\).
The polarization matrix is then
\(\left(\begin{smallmatrix}2&n\\n&2\end{smallmatrix}\right)\) for an
integer \(n\), and positivity gives \(n\in\{0,1,-1\}\).
Suppose \(n=\pm1\). Then \(\deg\lambda=3\), and its kernel is an
order-three subgroup of the graph \(\{(P,nP):P\in E_b[3]\}\).
Indeed the polarization matrix kills the kernel, and its own kernel
is that graph. The only rational cyclic subgroup of order three on
\(E_b\) is \(\langle D\rangle\), for \(D=(0,b)\).
The division polynomial is \(3x(x^3+4b^2)\), and a rational root in
its second factor would require \(2b\) to be a rational cube.
For the weighted coefficient \(b=4k(k+1)\), this would require
\(k(k+1)\) to be a cube. Writing \(k=a/c\) in lowest terms would
make the pairwise coprime nonzero integers \(a,c,a+c\) all signed
cubes, contrary to [Fermat's theorem for exponent three](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#witness-selection).

Now the elliptic graph \(P\mapsto(P,[n\zeta_3]P)\) contains
\(\ker\lambda=\langle(D,nD)\rangle\), since \([\zeta_3]D=D\).
The polarization restricts to degree three on this graph, as
\(4+\zeta_3+\zeta_3^2=3\).
Its image under \(\lambda\) thus inherits a polarization of degree one
as a line bundle, because the graph maps to its image with degree three.
An elliptic subvariety with principal restricted polarization splits
off as a principally polarized factor by the adjoint of its inclusion.
This contradicts indecomposability of the genus-two Jacobian, whose
theta divisor is geometrically irreducible. It follows that the
polarization matrix of \(\lambda\) is diagonal in every case.

Consequently the kernel of \(\lambda\) has order four and lies in the
product of the two groups of 2-torsion points. By injectivity of each
pullback, the kernel has
trivial intersection with both factors and is the graph of an
isomorphism of the two Galois modules of 2-torsion. This is the
degree-two gluing of
[Howe–Leprévost–Poonen, Proposition 3](sources/howe-leprevost-poonen-split-jacobians.pdf).

For an asymmetric profile \(B_i\) is noncubic by the coefficient
calculation in
[FM4](ABC_CUBIC_TANGENT_COLLISION_MOMENT.md#asymmetric-weighted-profiles).
Its 2-torsion has Galois image \(S_3\), by the roots of \(x^3+B_i^2\).
For a symmetric profile the coefficient is a rational cube and the
image has order two. Hence a symmetric and an asymmetric profile
cannot be paired. For two asymmetric profiles, the graph, if it exists,
is unique because the centralizer of \(S_3\) on the three nonzero
points is trivial. Existence gives equality of the splitting fields
\(\mathbb Q(\zeta_3,\sqrt[3]{B_i})\).
By cubic Kummer theory, the two coefficient classes are equal or inverse,
first over \(\mathbb Q(\zeta_3)\) and then over \(\mathbb Q\) by norms.
In the equal case the graph comes from the rational isomorphism of the
factors. It gives the decomposable polarized quotient since
\((P,Q)\mapsto(P+Q,P-Q)\) has the same graph kernel after that
identification and pulls back the product polarization to the diagonal
matrix with entries \([2]\). Thus only the inverse case can occur.

In the inverse case, scale the factors to \(E_b,E_{b^2}\) with
\(b=B_1\). The reciprocal \(C_b\) has those factors with the same
2-torsion graph, because its two projections are independent and their
off-diagonal entries vanish. The quotient isomorphism preserves the
principal polarizations, since the pullbacks on the product agree.
[Strong Torelli for hyperelliptic curves](sources/serre-strong-torelli-2001.pdf)
gives an isomorphism \(C\simeq C_b\) over \(\mathbb Q\) that
induces this polarized isomorphism. The two projections then induce
the same pushforwards on Jacobians as the reciprocal projections by
duality. Maps to an elliptic curve inducing the same pushforward
differ by a constant, here a rational point. Intertwining the deck
involution makes twice each translation equal to the associated
deck center. GT3 therefore applies.

Finally, two symmetric factors scale to \(E_1\times E_1\).
The centralizer of the order-two Galois image on the three nonzero
2-torsion points consists of the identity and the transposition.
The identity gives the decomposable quotient. The transposition gives
the reciprocal \(C_1\), because its two projections are orthogonal by
their opposite signs under \((s,Y)\mapsto(-s,Y)\), and its Jacobian is
indecomposable. The same Torelli and translation argument would give
rational halves of the symmetric deck centers. These centers have
order two or six and \(E_1(\mathbb Q)\) has order six, as in WM2.
The halves do not exist. This rules out the symmetric case.
