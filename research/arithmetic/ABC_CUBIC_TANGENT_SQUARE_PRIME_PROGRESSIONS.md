# The cubic tangent square and its single coherent class direction

The cubic endpoint of the higher-prime tangent-square
construction has genus one. Its two consecutive coherent
classes coincide, but one nontrivial marked class survives.
The quartic radicand and a uniform height count then give
a larger family of genuine reflected cubic fields. The same
witnesses give positive rank on a fixed Mordell twist.

## The exact class relation and the elliptic map

**Proposition CTQ1.** Put
\[
 F(t)=-3t^4+6t^2+1,\qquad C:y^2=F(t),
\]
and \(A=(t-1)^3,\ B=4t^3,\ C_0=(t+1)^3\). The two cubic
norm functions of the tangent-square construction are
\[
 \alpha_1=(A+B-C_0+2y)/4,\qquad
 \alpha_2=(B+C_0-A+2y)/4 .
\]
They have the same nonzero class in
\(\overline{\mathbb Q}(C)^\times/
 \overline{\mathbb Q}(C)^{\times3}\).
The following identities hold over \(\mathbb Q\):
\[
 \frac{\alpha_1}{\alpha_2}
 =\left(\frac{y-1-t^2}{2t(t+1)}\right)^3,\qquad
 1+V=\alpha_1
       \left(-\frac{y+3t^2+2t+1}{4t^2}\right)^3,
 \tag{1}
\]
where
\[
 U=\frac{F(t)}{4t^2},\qquad
 V=\frac{y(-3t^2-1/t^2)}{8t}.
\]
The map \((U-1,V):C\to E_1\), \(E_1:Y^2=X^3+1\), has
degree four and is anti-invariant under the hyperelliptic
involution. The generic quadratic twist of \(E_1\) by
\(F(t)\) has rank one.

The tangent-square polynomial of
[HP5](ABC_COHERENT_HIGHER_PRIME_PROGRESSION_FIELDS.md#removing-a-tangent-square-from-the-field-parameter)
at the cubic exponent is \(3p^4-6p^2h^2-h^4\).
The genus-one model is the symmetric instance of
[WM1](ABC_WEIGHTED_CUBIC_TANGENT_ELLIPTIC_MARKS.md#the-double-cover-and-its-three-marks).
The quartic \(F\) has distinct roots and the map is the
even-quartic formula (5) of HG1. Equivalently one can
substitute to verify \(V^2=(U-1)^3+1\) and (1).
The degree of \(U-1\) as a rational function on the
parameter line is four, and the two hyperelliptic covers
give the stated degree on \(C\).
The first cubic function has divided divisor with zeros
at \((0,1),(1,2)\) and poles at the two infinities, all
simple in the divided divisor. If it were principal,
its function would be regular on the affine curve and
have at most simple poles at both infinities.
A function \(P(t)+yQ(t)\) with these pole orders has
\(Q=0\) and \(\deg P\le1\). This cannot vanish on
only one sheet at \(t=0\).
Thus the class is nonzero. The nonconstant anti-invariant
map contributes the single rational Hom direction from
the genus-one Jacobian, which is isogenous to \(E_1\).
This proves the exact generic rank.

## Prime witnesses and reflection at three

**Theorem CTQ2.** There are genuine coherent normalized
\((1,1,2)\) progression triangles
\((p_1,p_2,p_3)=(p-h,p,p+h)\) with rational prime vertices
for either sign of the quadratic field \(K=\mathbb Q(\sqrt D)\).
On the congruences specified below, their two consecutive
norm classes coincide and are noncubic virtual units.
The \(D\)-twist of \(E_1\) has positive rank and conductor
\(4D^2\), where \(D\) is the fundamental discriminant of \(K\).
Their mirror has a cyclic unramified cubic extension
split above three. Its associated rational cubic field
has fundamental discriminant \(-D/3\) and has \(3\)
completely split. In the imaginary base field case the
three oriented prime ideal classes form a full coset
of one order-three subgroup.

Take \(p\) odd, \(p\equiv1\bmod3\),
\(h\equiv0\bmod18\), and \(p\equiv5,\ h\equiv1\bmod19\).
One can use any fixed compact interval of \(t=p/h>1\)
inside either sign region of \(F(t)\).
The prime forms are admissible after these congruences.

Use \(A=p_1^3,\ B=4p^3,\ C_0=p_3^3\) and
\(R=3p^4-6p^2h^2-h^4\), so the radicand is \(-R\).
Write \(E=h^2R\) and
\(u_1=(A+B-C_0)/2,\ u_2=(B+C_0-A)/2\).
The elements
\[
 \alpha_{12}=(u_1+\sqrt{-E})/2,\quad
 \alpha_{23}=(u_2+\sqrt{-E})/2,\quad
 \alpha_{13}=(B-A-C_0)/2+\sqrt{-E}
\]
have the three prime-cube norms and satisfy
\(\alpha_{12}\alpha_{23}=p^3\alpha_{13}\).
They are integral by their trace and norm.
The prime-unit calculation of HP5 applies at exponent
three: the outer constant \(4-2^3=-4\) is nonzero at
the odd vertex primes, while at the middle prime the
outer powers are opposites. Thus these are genuine
one-sided oriented prime-cube ideals, with the common
scale times \((1,1,2)\) for their determinant coefficients.

The ratio identity in (1) holds also for the two
homogeneous generators, since both specialize with the
same factor \(h^3\). Modulo \(19\), \(t=5\) and the
branch \(y=9\) gives \(\alpha_1=6\), whose sixth power
is \(11\ne1\). Thus the common class is noncubic.
For the matching imaginary orientation the ideal
classes satisfy
\[
 [\mathfrak p_2]-[\mathfrak p_1]
       =[\mathfrak p_3]-[\mathfrak p_2].
\]
The common difference is nonzero of order three.
Indeed in an imaginary field all units are cubes
apart from the possible \(\mu_3\) direction in
\(\mathbb Q(\sqrt{-3})\). Even those units become
cubes modulo \(19\), which is \(1\bmod9\),
so the noncube residue excludes a unit class.

The local argument of PA1 works at \(\ell=3\).
The two generators divided by \(p^3\) differ from
one by terms of \(3\)-adic valuation greater than
two. Thus both are local cubes. Here
\(-R\) has squareclass \(-3\) at \(3\), since
\(-R/3\equiv-p^4\equiv-1\bmod3\).
Consequently \(K\) ramifies at \(3\), its mirror
\(K'=\mathbb Q(\sqrt{-3D})\) splits there, and
\(D_{K'}=-D/3\).
A noncubic virtual unit with rational cube norm
and local cubes at \(3\) gives an unramified cyclic
cubic over \(K'\), split at \(3\), by the descent
of MC1 or MT2. The cubic trace polynomial is
\[
 T^3-3p_1p_2T-u_1.
 \tag{2}
\]
Its normal closure is that reflected extension
and its field discriminant is \(D_{K'}\), by the
unramified relative conductor formula.
The case of a real base field needs only the
virtual-unit noncube test. That test survives even if
the ordinary ideal class of the virtual unit is trivial.

The second relation in (1) gives noncubic Kummer
evaluation \(1+V\) on the anti-invariant point.
The cubic factor of the \(3\)-division polynomial of
\(E_1\) is \(X^3+4\), irreducible. Thus no quadratic
anti-invariant \(3\)-torsion can account for this
nonzero evaluation, and the twist has positive rank.
The residue \(-R\equiv5\bmod8\) gives \(D\equiv5\bmod8\).
Hence \(2\) is nonsplit unramified, with the fixed
local squareclass \(5\), and \(D\) is odd.
The original conductor is \(36\). The quadratic
twists by \(5\) and \(-3\) retain its local exponents
two at \(2,3\), respectively. At every prime
\(q\ge5\) dividing \(D\), the base curve has good
reduction and the twist has conductor exponent two.
This gives \(36(D/3)^2=4D^2\).

## Counting with the genus-one multiplicity bound

**Corollary CTQ3.** For each base sign the family of
CTQ2 has at least
\(\gg X^{1/2}/(\log X)^B\) distinct base fields of
absolute discriminant at most \(X\), for a fixed
constant \(B>0\). The same lower count holds for
their reflected cubic fields, with absolute
fundamental discriminant at most \(X/3\).
The positive-rank twist conductor count is
\(\gg Y^{1/4}/(\log Y)^B\).

Fix a positive rectangle \(p,h\asymp Z\) on one
side of the quartic sign change, at \(p/h>1\)
and away from the zeros of \(F\).
The prime-form theorem used in HP5 gives
\(\gg Z^2/(\log Z)^3\) progression parameters on
the fixed congruence slice. The radicand \(R\) is
a fixed separable binary quartic of size \(O(Z^4)\).
The prime-factor truncation of
[PP3](ABC_COHERENT_PRIME_PROGRESSION_NORMS.md#quantitative-count-of-genuine-fields)
applies with four in place of six for its
root bound. Outside \(O(Z^2/(\log Z)^4)\)
primitive pairs, the squarefree part has
\(\omega(D)=O(\log\log Z)\).

For fixed \(D\), apply the
[height-packing argument of DM2](MORDELL_CUBIC_DOUBLED_SIMULTANEOUS_TWISTS.md#counting-distinct-simultaneous-twists)
to the fourth-degree map of CTQ1 and the fixed curve
\(E_1\). Its fixed étale cubic algebra for \(2\)-descent
can have more than one field factor. The same argument
bounds rank by \(O(1+\omega(D))\), since the unit and
class-group contributions from that algebra are fixed.
The \(X\)-coordinate of the anti-invariant point is a
rational function of \(p/h\) of degree four, so its
logarithmic height is \(O(\log Z)\) and each abscissa
has at most four parameters.
The uniform packing estimate yields
\[
 \ll \left(1+C\sqrt{\frac{\log(2Z)}{\log(2|D|)}}\right)^{
 O(1+\omega(D))}
\]
parameters over \(D\), at most
\(\exp(O(\log Z/\log\log Z))\) uniformly in \(D\).
Together the fields with \(|D|\le Z\) therefore
account for \(Z^{1+o(1)}\) parameters.
Above this range and on the prime-factor truncation,
the displayed multiplicity is a fixed power of
\(\log Z\). Dividing the count by it and taking
\(Z\asymp X^{1/4}\) with a sufficiently small constant
gives the claim.
Both maps from \(D\), to \(-D/3\) and to \(4D^2\)
on a fixed sign, are injective.

## Prescribing the two outer decompositions

**Corollary CTQ4.** The two outer vertices of the
reflected cubic field may independently be prescribed
split or inert. The middle vertex always splits.
The lower bounds of CTQ3 hold with either sign and
either fixed pair of outer prescriptions.

The local normalization of PA5 holds at exponent
three for the dependent pair. At both outer vertices
the common marked direction has cube residue that of
\(1/2\), since \(2/(4-2^3)=-1/2\) has the same cube
class. The middle residue is trivial.
All three primes are \(1\bmod3\) on CTQ2.
Apply the fixed Artin restrictions for
\(\mathbb Q(\zeta_3,\sqrt[3]{2})\) at the two outer
linear forms, as in
[Frei–Tinková, Theorem 1.4](https://arxiv.org/pdf/2408.17441).
The identity and the conjugate nonidentity rotations
restrict trivially on the common cyclotomic quadratic
subfield. The existing prime form congruences are
therefore compatible and all local factors remain
positive. The genus-one multiplicity argument of
CTQ3 is unchanged.

The [exact checker](verify_cubic_tangent_square_prime_progressions.py)
replays the two cube identities and the elliptic map,
and checks the rank-detecting residue, fundamental
reflected cubic discriminants, local prime types and
conductors for both signs. For example the imaginary
base fixture \(p=2551,h=1008\) has
\(D=-86341391658723\) and cubic polynomial
\(T^3-11808579T-12498602366\), of field
discriminant \(28780463886241\).
