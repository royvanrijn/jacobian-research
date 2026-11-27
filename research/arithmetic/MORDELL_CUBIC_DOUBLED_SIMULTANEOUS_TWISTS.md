# Simultaneous negative quadratic twists from the doubled cubic family

The doubled branch model gives a quantitative family of negative
quadratic twists of two fixed \(j=0\) curves. Both twists have
positive rank, even \(3\)-Selmer dimension, and a prescribed
nontrivial Kummer evaluation at a coefficient prime.
Their fields have the split-prime class constraint of
[DB5](ABC_CUBIC_TANGENT_DOUBLED_BRANCH_FAMILY.md#cubes-of-split-prime-ideal-classes).

Use \(b=6\), \(m=36\), \(w=36s^3\), and the two projections
of [DB1](ABC_CUBIC_TANGENT_DOUBLED_BRANCH_FAMILY.md#classification-of-the-doubled-branch-model).
Write
\[
 \Delta(s)=(w^2-1)(w^2+4w-1),\qquad E_b:y^2=x^3+b^2 .
\]
For rational \(s\in(1/5,1/4)\), let \(d\) be the signed squarefree
part of \(\Delta(s)\).
The quadratic twists will be written
\[
 E_b^{(d)}:y^2=x^3+b^2d^3
 \qquad (b=6,36).
\]
For comparison, [Wong, Theorem 1](sources/wong-quadratic-twists-of-pairs.pdf)
gives infinitely many simultaneous positive-rank quadratic twists
for pairs of elliptic curves not both of \(j\)-invariant zero.
His Theorem 2 applies to \(j=0\) after admitting a common sextic
twist of the original curves. The present family keeps the two
specific \(j=0\) curves fixed and additionally imposes the
coefficient-prime evaluations and splitting of their quadratic fields.

## The marked points and conductor pair

**Theorem DM1.** Every \(d\) from that interval is negative,
with \(d\equiv1\bmod24\). Both curves \(E_b^{(d)}\) have a point
of infinite order that is not divisible by three, and both have
even \(3\)-Selmer dimension at least two.
Their conductors are \(972d^2\) and \(243d^2\), respectively.

The field is imaginary by DB3 and splits at \(2,3\) by DB2,
which gives the congruence. The two points given by (1) of
DB1 are anti-invariant in the quadratic field and therefore
give rational points on \(E_b^{(d)}\). On \(E_b\) over the
quadratic field the \(3\)-isogeny Kummer evaluation is
\(y+b\), up to inversion, by
[Bhargava–Elkies–Shnidman, Remark 25](sources/bhargava-elkies-shnidman-mordell-isogeny-selmer.pdf).
By DB2 both evaluations are nontrivial at \(3\).
Thus the points are not divisible by three.
The division polynomial \(3x(x^3+4b^2)\) shows that the only
points of order three defined over a quadratic field are the
rational points \((0,\pm b)\): \(2b\) is noncubic for \(b=6,36\).
None of these nonzero points is anti-invariant.
Any torsion point in the anti-invariant subgroup would therefore
have order prime to three and have trivial \(3\)-isogeny evaluation.
This proves that the displayed points have infinite order.

The base conductors are \(972=2^2 3^5\) and \(243=3^5\),
as in the exact checker below. The local twists at \(2,3\)
are isomorphic to the base models. At each \(q\ge5\) dividing
\(d\), the original model has good reduction and its twist has
type \(I_0^*\), since the short Weierstrass constant has
valuation three. Its conductor exponent is two, while all
other primes of this size have good reduction.
This gives the two conductors.

For the rational isogeny from the constant term \(b^2d^3\)
to \(-27b^2d^3\), use the local Selmer ratios of
[Bhargava–Elkies–Shnidman, Propositions 34 and 39](sources/bhargava-elkies-shnidman-mordell-isogeny-selmer.pdf).
At \(q\ge5\) dividing \(d\) the coefficient valuation is
three, so the ratio is one. At good primes the ratio also
equals one. At \(2\) the character is split, and the valuations
of four times the coefficient are \(4,6\) for \(b=6,36\);
the two ratios are \(1/3,1\).
At \(3\) the valuations of the coefficients are \(2,4\),
respectively, and their unit parts are \(1\bmod3\).
The ratios are again \(1/3,1\).
The real ratio is one because the coefficient is negative.
Thus the global ratios are \(1/9,1\).
There is no rational three-torsion on either twist by the
anti-invariance argument. Proposition 49 of the same source
therefore gives an even finite \(3\)-Selmer dimension for both
curves. It is at least two by the exhibited points.

## Counting distinct simultaneous twists

**Corollary DM2.** There are at least
\[
 X^{1/6}\exp\!\left(-O\!\left(\frac{\log X}{\log\log X}\right)\right)
\]
distinct negative \(d\) of absolute value at most \(X\) with
all the conclusions of DM1 and DB5.
Equivalently the number of distinct conductor pairs with
both conductors at most \(Y\) is at least
\(Y^{1/12}\exp(-O(\log Y/\log\log Y))\).
The estimate needs no squarefree-value assertion for the
degree-twelve polynomial.

Write \(s=a/c\) in lowest terms with \(c>0\), and use height
\(\max(|a|,c)\le H\).
The fixed interval contains \(\gg H^2\) rational parameters
of that height by counting primitive lattice vectors.
Their field discriminants have size \(O(H^{12})\):
clear denominators from the degree-twelve polynomial
\(\Delta(s)\) using \(c^{12}\), which is a square.
The abscissa \(x=2ms^2/(w-1)\) on \(E_6\) has logarithmic
height at most \(3\log H+O(1)\).
For the height count we use only this first factor.

For a fixed field \(K=\mathbb Q(\sqrt d)\), take the
anti-invariant Mordell–Weil group on the fixed \(E_6(K)\),
which is identified with \(E_6^{(d)}(\mathbb Q)\).
Canonical height on the fixed model satisfies
\(2\widehat h(P)=h(x(P))+O(1)\) as in
[FC5](ABC_CUBIC_TANGENT_FIELD_COLLISIONS.md#multiplicity-uniform-in-the-field).
For an anti-invariant point \(Q\) of infinite order with
\(x(Q)=u/v\in\mathbb Q\) in lowest terms, the integer
\(v(u^3+36v^3)\) has squareclass \(d\).
It follows that
\(\widehat h(Q)\ge \tfrac18\log |d|-O(1)\).
The finitely many small fields can be absorbed by taking a
smaller positive constant, so distinct lattice vectors are
separated by at least \(c_0\sqrt{\log(2|d|)}\).
Torsion of the fixed curve in extensions of degree at most
two has bounded size, by the same height comparison and
Northcott finiteness. The degree-three abscissa map then
bounds the number of parameters per lattice vector by a
fixed constant.

Put \(r_d=\operatorname{rank}E_6^{(d)}(\mathbb Q)\).
The étale algebra in the ordinary \(2\)-descent is the
fixed cubic field \(L=\mathbb Q(\sqrt[3]{-36})\), because
the cubic \(x^3+36d^3\) has a root \(d\sqrt[3]{-36}\).
The descent injection and prime support of
[Stoll, Corollary 4.7 and Lemma 4.9](sources/stoll-two-descent-hyperelliptic-jacobians.pdf)
bound \(r_d\) by \(O(1+\omega(d))\).
Indeed the valuation coordinates are supported on at most
three primes of \(L\) over each prime dividing \(6d\);
the everywhere-even classes are bounded by the units and
the two-torsion of the class group of the fixed \(L\).
The descent hypotheses hold for the separable monic cubic
and its rational point at infinity.

Packing in the \(r_d\)-dimensional height lattice gives,
for one field, at most a constant times
\[
 \left(1+C\sqrt{\frac{\log(2H)}{\log(2|d|)}}\right)^{r_d}
\]
parameters. As in FC5, this is
\(\exp(O(\log H/\log\log H))\) uniformly, since
\(\log(2|d|)\ll\log H\). For completeness, when
\(\log(2|d|)\ge\sqrt{\log(2H)}\) use
\(\omega(d)\ll\log(2|d|)/\log\log(4|d|)\) and
\(\log(1+u)\le u\). In the complementary range the
logarithm of the bound is
\(O(\sqrt{\log H}\log\log H)\).
Divide the parameter count by this multiplicity estimate
and put \(H\asymp X^{1/12}\) with a small fixed constant.
The conductor assertion follows since each conductor is
a fixed multiple of \(d^2\).

## The generic ranks

**Proposition DM3.** Each of the two quadratic twist curves
over \(\mathbb Q(s)\) by \(\Delta(s)\) has generic
Mordell–Weil rank exactly one.

The two explicit maps give nonconstant anti-invariant points
on the constant curves over the function field of
\(C:Y^2=\Delta(s)\). For the upper bound use the
identification of rational maps from \(C\) to a constant
elliptic curve with a constant translation and a homomorphism
from its Jacobian. Here \(C\) has rational points above \(s=0\),
so the identification holds over \(\mathbb Q\).
The hyperelliptic involution acts as \(-1\) on the Jacobian.
Thus the rank of the anti-invariant group is bounded by the
rank of the rational Hom group from \(\operatorname{Jac}C\)
to the corresponding target.

At the good prime \(13\), the characteristic polynomial of
Frobenius on \(\operatorname{Jac}C\) is
\[
 (T^2-5T+13)(T^2-2T+13)(T^2+7T+13)
 (T^4-4T^3+3T^2-52T+169).
\]
The two target Frobenius polynomials are \(T^2-5T+13\)
and \(T^2+7T+13\), both of multiplicity one.
Any rational isogeny decomposition specializes at \(13\),
so each target contributes with multiplicity at most one
over \(\mathbb Q\). The rational endomorphism ring of each
target is \(\mathbb Z\). This proves the upper bounds and
hence both generic ranks are exactly one.
The [Sage checker](verify_mordell_cubic_doubled_simultaneous_twists.py)
replays the base conductor and local isogeny computations
and the finite-field characteristic polynomial.
