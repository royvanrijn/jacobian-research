# Rank two in both quadratic mirror class groups from marked Mestre twists

Two independent Mordell–Weil points need not supply two
independent virtual units. The following specialization of
Mestre's rank-two construction controls the two Kummer
functions at every finite prime. Two split residue tests
give independence, while both functions are cubes at \(3\).
This gives nontrivial rank two in each mirror class group
and on a fixed elliptic curve's quadratic twist.

## Two marked functions on a rank-two twist

**Theorem MT1.** Let
\[
 E:\ y^2=x^3+(3x+17)^2 ,
\]
and use the short coordinates \(Z=x+3\), in which
\(y^2=Z^3+75Z+37\).
For a rational parameter \(u\ne0\), write
\[
 A=u^4+u^2+1,\quad B=u^2+1,\quad
 Z_1=-\frac{37A}{75u^2B},\quad Z_2=u^2Z_1,\quad x_i=Z_i-3.
\]
Put
\[
 P(u)=-75\cdot37\,B\bigl(37^2 A^3+75^3 u^4 B^2\bigr).
 \tag{1}
\]
The points with abscissas \(x_1,x_2\) on \(E\), and
ordinates \(y_1=Y/(75^2u^3B^2)\) and \(y_2=u^3y_1\)
on \(Y^2=P(u)\), give independent generic cubic norm
classes
\[
 a_i=y_i+3x_i+17,\qquad \operatorname N(a_i)=-x_i^3.
 \tag{2}
\]
The generic twist by \(P(u)\) has Mordell–Weil rank exactly two.

The coordinates are the standard two-abscissa identity in
[Mestre's construction as given by Rubin–Silverberg, Theorem 3.7](sources/rubin-silverberg-rank-frequencies-quadratic-twists.pdf).
In fact, if \(f(Z)=Z^3+75Z+37\), then
\(f(u^2Z_1)=u^6 f(Z_1)\) follows on substituting the
displayed \(Z_1\). Since \(A-u^2B=1\),
clearing the square denominator gives (1).
The polynomial is squarefree of degree fourteen, as is
already seen by its reduction at \(7\), so this is a
hyperelliptic curve of genus six.
Each abscissa map has degree four and thus each elliptic
projection has degree four.
The function \(y+3x+17\) has divisor three times the
difference of the marked point \((0,-17)\) and infinity.
Each pullback of its divided divisor is nonzero by
degree-four pushforward. The involution
\((u,Y)\mapsto(-u,Y)\) negates the first elliptic point and
fixes the second. Their order-three divided divisor classes
therefore lie in opposite eigenspaces and are independent
even over the algebraic closure. The same eigenspace argument
on rational Hom groups makes the two nonconstant elliptic
maps independent.

For the upper rank bound, on the Jacobian of (1) the
Frobenius polynomial at the good prime \(7\) is
\[
 (T^2+T+7)^2 (T^4+T^3-2T^2+7T+49)^2 .
\]
The first quadratic is the Frobenius polynomial of \(E\)
and occurs exactly twice. Thus the rational isogeny
multiplicity of \(E\) is at most two.
The \(j\)-invariant \(108000000/63869\) is nonintegral,
so the endomorphism ring over the algebraic closure is
\(\mathbb Z\). Rational maps to the constant elliptic curve
modulo constants inject into the rational Hom group from
the Jacobian. This gives the upper bound for the
anti-invariant rank, and hence exact rank two for the twist.

## The common local conditions and the mirror classes

**Theorem MT2.** Suppose \(u\) satisfies
\[
 u\in\mathbb Z_2^\times,\qquad
 u\equiv1\bmod27,\qquad u\equiv2\bmod17,\qquad
 u\equiv2\bmod7,\qquad u\equiv13\bmod19,
 \tag{3}
\]
where each congruence includes integrality at its prime.
Then \(K=\mathbb Q(\sqrt{P(u)})\) is imaginary and
both classes in (2) are virtual units, independent in
\(K^\times/K^{\times3}\), and cubes at every prime over \(3\).
Both \(K\) and the real mirror
\(K'=\mathbb Q(\sqrt{-3P(u)})\) have class-group \(3\)-rank
at least two. The quadratic twist of \(E\) has rational
rank at least two, witnessed by the displayed points.

The polynomial in (1) is negative on the real line.
At \(2\), for \(u\) odd one has \(v_2(B)=1\) and
\(A\) odd. Both \(Z_i\), and hence both \(x_i\), have
valuation \(-1\) at the rational prime. In each local
extension the pole of \(y_i\) has order divisible by
three and dominates \(3x_i+17\), giving the virtual-unit
condition. The same pole argument applies at any prime.
Where \(x_i\) is integral at an odd prime not over \(17\),
a positive valuation of \(a_i\) requires \(x_i\) to
vanish modulo that prime by the cube norm.
Then \(3x_i+17\) is a unit, the other norm conjugate is
a unit, and the valuation is a multiple of three.
At \(17\) the prescribed residue gives \(x_1=11,\ x_2=2\)
modulo \(17\); both norm valuations are zero.
This proves the virtual-unit condition at all primes.

At \(3\) both \(x_i\) are congruent modulo nine to their
common value \(-187/50\) at \(u=1\).
Indeed \(Z_i\) has a coefficient of valuation \(-1\)
at \(3\), and its other denominator factors are units.
Changing \(u\) modulo \(27\) thus preserves \(x_i\bmod9\).
The ordinate square is \(-4\bmod9\) and
\(3x_i+17\equiv2\bmod9\).
The completion of \(K\) at \(3\) is the unramified
quadratic extension, which contains \(i\) with \(i^2=-1\).
We have \(a_i\equiv2\pm2i\bmod9\).
The eighth power of this unit is \(1\bmod9\).
In this unramified quadratic field, any such unit is
a cube: cubing is bijective on its Teichmüller group of
order eight and carries \(1+3\mathcal O\) onto
\(1+9\mathcal O\). The latter follows from the expansion
\((1+3r)^3\) and Hensel's lemma over its ring of integers.

At \(7\) and \(19\), the following exact residues in one
split branch give independent cube characters:
\[
 \begin{array}{c|r|r|r|r|r}
 q&u&y_1&a_1&a_2&
 (a_1^{(q-1)/3},a_2^{(q-1)/3})\\ \hline
 7&2&3&4&4&(2,2)\\
 19&13&3&16&15&(7,11).
 \end{array}
\]
The first row gives two equal nontrivial characters
and the second gives two inverse nontrivial characters.
Changing the prime over \(q\) inverts the whole row,
by the cube norms. These rows are independent over
\(\mathbb F_3\), so the classes in (2) are independent
in every specialization satisfying (3).

The field is not \(\mathbb Q(\sqrt{-3})\) because it is
unramified at \(3\), so its units are cubes and the
two virtual classes give independent ideal classes.
For the mirror use the eigenspace argument of
[MC1](ABC_CUBIC_TANGENT_MIRROR_CLASS_FAMILY.md#the-local-cube-and-the-real-cyclic-cubic-extension)
on both classes at once.
They remain independent after adjoining \(\zeta_3\)
by the degree-two norm, and their two cubic radicals
give an elementary abelian degree-nine extension there.
It is unramified away from \(3\) by the virtual valuations
and locally trivial at \(3\).
Conjugation fixing the mirror acts trivially on its
Kummer Galois group because the norms of both classes
are cubes. Quotienting by the subgroup of order two gives
an elementary abelian unramified degree-nine extension
of \(K'\). In particular each of its four cyclic cubic
subextensions supplies a distinct totally real rational
cubic field of fundamental discriminant
\(\operatorname{disc}(K')\).

For the rank assertion, the two points are anti-invariant,
and the function in (2) is the ordinary
[Cohen–Pazuki isogeny evaluation](sources/cohen-pazuki-2009-three-isogeny.pdf),
up to inversion, after passing to the twist coordinates
\(x'=dx,\ y'=d\sqrt d\,y\) for the squareclass \(d\).
The division polynomial is
\(3x(x^3+12x^2+204x+1156)\) with irreducible cubic factor.
Thus the only points of order three in a quadratic
extension lie in the rational marked line, and none of
its nonzero points is anti-invariant.
The anti-invariant group modulo multiplication by three
has dimension equal to its Mordell–Weil rank. Its
images in the Kummer group include the two independent
classes just computed, so its rank is at least two.

## Quantitative distinct fields

**Corollary MT3.** At least
\[
 X^{1/7}\exp\!\left(-O\!\left(\frac{\log X}{\log\log X}\right)\right)
\]
distinct imaginary quadratic fields of absolute discriminant
at most \(X\) satisfy MT2. Their real mirrors are distinct
with discriminants of absolute value at most \(3X\).
Thus each of those positive fundamental discriminants
supports at least four totally real cubic fields.
The count gives the indicated rank-two twists and does
not assume a squarefree-value theorem for (1).

The congruences (3) admit one rational residue class with
denominator coprime to \(2\cdot3\cdot7\cdot17\cdot19\).
On a fixed interval such as \(1/2<u<2\) this gives
\(\gg H^2\) primitive rational parameters of height at
most \(H\). Since (1) has degree fourteen, the field
discriminants are \(O(H^{14})\).
The fixed elliptic curve \(E\) and the first abscissa of
degree four give the same uniform bound as the
[height-packing argument of DM2](MORDELL_CUBIC_DOUBLED_SIMULTANEOUS_TWISTS.md#counting-distinct-simultaneous-twists).
Its cubic algebra for ordinary two-descent is fixed across
all quadratic twists, so the rank is \(O_E(1+\omega(d))\).
This yields \(\exp(O(\log H/\log\log H))\) parameters per
field. Take \(H\asymp X^{1/14}\) with the constant chosen
to keep the discriminants below \(X\). The mirror
discriminants have size \(3|D_K|\), because \(3\) is
unramified in the imaginary field.

The [Sage checker](verify_marked_mestre_mirror_rank_two.py)
replays the twist identity, the generic-rank Frobenius
factor, and all three local residue computations. It also
replays the congruences at the rational parameter \(271/433\).
The construction concerns virtual units and actual
Mordell–Weil points of a fixed curve. It makes no assertion
about prime support for the original weighted tangent
witnesses.
