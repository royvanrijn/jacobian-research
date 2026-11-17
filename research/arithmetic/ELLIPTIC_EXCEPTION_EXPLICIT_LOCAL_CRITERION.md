# An explicit all-primes local criterion on the existing elliptic slice

On the fixed elliptic exception, every varying tame ramification test can
be read from one rational function. Outside the three primes \(2,3,5\),
Hasse failure is exactly a divisibility-by-three condition at primes
\(p\equiv2\pmod3\). This is an equivalence, not an infinitude theorem.

## 1. The fixed slice and the criterion

Retain the admissible matching-resolvent slice
\((\pi,\ell)=(128/25,9/10)\) from the
[elliptic arithmetic note](ELLIPTIC_EXCEPTION_ARITHMETIC.md). Write
\[
 E_0:\quad y^2=x^3-10656x^2+27091730x,\qquad
 x=X+3552=\frac{6400}{3}v-432.
\]
Its quadratic factor and cubic factor are
\[
 q(W)=W^2+\frac9{10}W+v,\qquad
 h(W)=W^3-\frac{12}{5}W^2+(s_0-v)W+t(v),
\]
where
\[
 s_0=\frac{1816347}{819200},\qquad
 t(v)=\frac{10}{9}v^2-\frac{10s_0}{9}v+\frac{125}{1024}.
\]
The admissible open has
\(\operatorname{Disc}(q)\operatorname{Disc}(h)\operatorname{Res}(q,h)\ne0\).
In this note a **full field fibre** means an admissible point for which
\(h\) is irreducible over \(\mathbb Q\). The previous arithmetic theorem
makes \(q\) irreducible at every admissible point. Thus the complete fibre
is \(K\times L\), where \(L=\mathbb Q[W]/h\) is a non-Galois cubic field
and \(K\) is its quadratic resolvent. Such a fibre has no rational point.

Put
\[
 z=v-\frac{243483}{819200}
   =\frac{3F(x)}{819200},\qquad F(x)=128x-25865,
 \qquad \lambda=\frac{10}{9},\quad \kappa=\frac{42649}{16384}.
                                                               \tag{1.1}
\]
Admissibility implies \(z\ne0\). The exact translation \(W=4/5+w\) gives
\[
 h(4/5+w)=g(w)=w^3-zw+z(\lambda z-\kappa).             \tag{1.2}
\]
The integer \(42649\) is prime and is \(1\pmod3\).
All valuations below are normalized by \(v_p(p)=1\).

**Theorem ELV1 (exact local criterion).** A full field fibre on this
slice is Hasse-failing if and only if:

1. at each of \(2,3,5\), at least one of \(q,h\) has a root in
   \(\mathbb Q_p\); and
2. for every prime \(p\equiv2\pmod3\), \(p\ge11\),
\[
 \boxed{v_p(F(x))\in3\mathbb Z.}                       \tag{1.3}
\]

There is no exception at \(1889\), even though that prime occurs in
the discriminant of the polynomial defining the branch values.
The three remaining small-prime tests are now explicit: use
[ETC1](ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md)
at three and [ETF2--ETF3](ELLIPTIC_EXCEPTION_TWO_FIVE_ADIC_CLASSIFICATION.md)
at two and five. They retain every rational point on the same slice;
they do not require a new parametrized family or a maximal-order
computation for each cubic field.

**Proof at a prime \(p\equiv2\pmod3,\ p\ge11\).** Both \(\lambda\)
and \(\kappa\) are units: their exceptional primes are \(2,3,5,42649\),
and \(42649\equiv1\pmod3\). Let \(k=v_p(z)\). Since \(3/819200\) is a
unit, \(k=v_p(F(x))\).

If \(k>0\), then \(c=\lambda z-\kappa\) is a unit. The lower Newton
polygon of \(g\), in ascending powers of \(w\), has its single edge
from \((0,k)\) to \((3,0)\); the point \((1,k)\) lies strictly above it.
Every root has valuation \(k/3\). If \(3\nmid k\), \(h\) has no
\(\mathbb Q_p\)-root. Moreover
\[
 \operatorname{Disc}(g)=z^2\{4z-27(\lambda z-\kappa)^2\}             \tag{1.4}
\]
has squareclass \(-3\): the bracket is \(-27c^2\) times a unit
congruent to one modulo \(p\). Such a unit is a local square by Hensel's
lemma. Matching discriminant squareclasses therefore identify the
local quadratic factor with \(\mathbb Q_p(\sqrt{-3})\). This is
nonsplit when \(p\equiv2\pmod3\), so neither factor supplies a point.
If \(k=3m>0\), write \(z=p^{3m}u\) and \(w=p^m U\). Then
\[
 p^{-3m}g(p^mU)=U^3-p^m uU+u(\lambda p^{3m}u-\kappa).             \tag{1.5}
\]
Its reduction is \(U^3-\bar u\bar\kappa\). Cubing is a bijection of
\(\mathbb F_p^*\), giving a nonzero simple root, which Hensel lifts.

If \(k<0\), one has \(v_p(c)=k\). The Newton polygon is the single
edge from \((0,2k)\) to \((3,0)\), and all roots have valuation \(2k/3\).
For \(3\nmid k\) there is no cubic root. In (1.4), factor the bracket
as
\[
 -27\lambda^2z^2
 \left\{\left(1-\frac{\kappa}{\lambda z}\right)^2
                 -\frac4{27\lambda^2 z}\right\}.
\]
The expression in braces is one modulo \(p\), so the discriminant
again has squareclass \(-3\) and the quadratic supplies no point.
If \(k=-3m<0\), put \(z=p^{-3m}u\), \(w=p^{-2m}U\). This time
\[
 p^{6m}g(p^{-2m}U)
   =U^3-p^m uU+\lambda u^2-\kappa p^{3m}u.             \tag{1.6}
\]
Its reduction \(U^3+\bar\lambda\bar u^2\) again has a nonzero simple
root.

If \(k=0\), \(g\) is integral and its linear coefficient is a unit.
It cannot reduce to \((w-r)^3\): its zero quadratic coefficient
would force \(r=0\), contradicting that unit linear coefficient.
A singular cubic over \(\mathbb F_p\) without a triple root has a
simple root. A separable reducible cubic also has a simple root.
Both cases give a local root of \(h\). In the remaining case the
residue cubic is irreducible and separable. Frobenius acts on its
three roots by a three-cycle, an even permutation, so its unit
discriminant is a square in \(\mathbb F_p\). Hensel lifts this
square, and matching squareclasses make the quadratic split.
Thus the product has a point whenever \(k=0\).

These arguments prove that the product has a \(\mathbb Q_p\)-point
exactly when \(3\mid k\). They do not assume good reduction of the
branch-value polynomial, and apply at \(1889\) as well.

At every \(p\equiv1\pmod3\) the product is locally soluble by
[RDF1--RDF3](RANK_FIVE_DISCRIMINANT_CLASSIFICATION.md): a non-total
cubic local algebra supplies either a root or a split resolvent,
whereas a totally ramified cubic has resolvent squareclass \(-3\),
which splits at these primes. The real cubic has a real root.
The remaining finite places are exactly \(2,3,5\), giving both
directions of the theorem. \(\square\)

The Newton and Hensel inputs are Milne, *Algebraic Number Theory*,
version 3.00, Theorem 7.32 and Proposition 7.44, printed pp.114 and 118--119
(PDF pp.122 and 126--127):
[author-hosted primary text](https://www.jmilne.org/math/CourseNotes/ANT300.pdf),
[retained copy](sources/milne-ant300.pdf).
The polygon above uses ascending coefficient indices; Milne states
the equivalent convention with descending coefficients. The finite-field
discriminant assertion follows by applying Frobenius to the
Vandermonde square root: it fixes that square root exactly for an
even permutation.

## 2. The maximal-order conductor away from four primes

Write \(D_L=d_K f^2\), where \(D_L\) is the cubic **field**
discriminant and \(d_K\) the quadratic fundamental discriminant.
This \(f\) is the invariant of RDF1, not the squarepart of an
arbitrary defining-polynomial discriminant.

**Theorem ELV2 (tame conductor support).** For every full field fibre
on this slice and every prime \(p\notin\{2,3,5,42649\}\),
\[
 \boxed{p\mid f\quad\Longleftrightarrow\quad
                    v_p(F(x))\notin3\mathbb Z.}        \tag{2.1}
\]
When \(p\mid f\) here, \(v_p(f)=1\).

**Proof.** Now \(p\) may be either one or two modulo three.
The same coefficient valuations and the normalizations (1.5)--(1.6)
apply.

For \(k\ne0\) not divisible by three, every root of \(g\) has
nonintegral valuation with denominator three. Consequently there is
no root in \(\mathbb Q_p\), so the cubic is irreducible there. The
valuation of a root forces its degree-three extension to have
ramification index three. It is a totally ramified cubic field.
To pass to its maximal-order discriminant, choose a local
uniformizer. Its Eisenstein minimal polynomial generates the maximal
local ring, and its discriminant has valuation two for \(p\ne3\);
hence \(v_p(f)=1\), exactly as proved in RDF3. This does not use the
discriminant of a nonmaximal scaled presentation.

For \(k\ne0\) divisible by three, the normalized reductions are
\(U^3+\text{unit}\). They are separable for \(p\ne3\), whether or not
they have a residue root. Hensel factorization lifts each irreducible
residue factor to the corresponding unramified extension. The local
cubic algebra is therefore a product of unramified fields, so \(p\)
does not divide \(f\).

For \(k=0\), suppose instead that the local cubic were a totally
ramified field. Its integral generator \(w\) and all its conjugates
would reduce to one element \(r\in\mathbb F_p\), making the reduction
\((w-r)^3\). The zero quadratic coefficient and unit linear
coefficient rule this out as above. Thus the cubic is not totally
ramified. RDF3 again says \(p\nmid f\). This proves (2.1) and its
exponent assertion. \(\square\)

The maximal-order uniformizer statement is Milne, Proposition 7.55
and Remark 7.56, printed pp.122--123 (PDF pp.130--131); its application
and discriminant calculation are retained in RDF3. Hensel
factorization is Milne, Theorem 7.33, printed p.114 (PDF p.122).
These are imported local-field results, not conclusions of the
finite controls below.

For the certified
[elliptic Hasse example](ELLIPTIC_EXCEPTION_HASSE_EXAMPLE.md),
\[
 x=\frac{173387072000}{9601^2},\qquad
 F(x)=\frac{5\cdot13\cdot1321\cdot230703199}{9601^2}.               \tag{2.2}
\]
Outside the exceptional set, (2.1) recovers exactly
\(13,1321,9601,230703199\), agreeing with the independently certified
[maximal-order conductor](ELLIPTIC_EXCEPTION_FIELD_DISCRIMINANT.md).
The exceptional prime five illustrates why (2.1) is not asserted
there.

## 3. An integer form of the sieve

**Corollary ELV3 (numerator and denominator).** Write a finite rational
point in reduced form \(x=a/c^2\), with \(a\in\mathbb Z\),
\(c\in\mathbb Z_{>0}\), and \(\gcd(a,c)=1\). Put
\[
 N=128a-25865c^2.
\]
Then condition (1.3) is equivalent to
\[
 \boxed{v_p(N)\in3\mathbb Z,\qquad v_p(c)\in3\mathbb Z
       \quad(p\equiv2\pmod3,\ p\ge11).}               \tag{3.1}
\]
The numerator is nonzero at every admissible point.

**Proof.** The denominator of an elliptic \(x\)-coordinate in this
integral monic model is a square: at a negative valuation, the term
\(x^3\) uniquely has smallest valuation in the cubic, so \(3v_p(x)\)
must be even. Also \(\gcd(N,c)\) divides \(128\).
Thus, for the odd primes in (3.1), the two terms in
\[
 v_p(F(x))=v_p(N)-2v_p(c)
\]
cannot both be positive. Since two is invertible modulo three,
their difference is divisible by three exactly when each term is.
\(\square\)

## 4. Fixed cubeclasses cannot give infinitude

**Theorem ELV4 (genus-four obstruction).** For each fixed
\(r\in\mathbb Q^*\), there are only finitely many rational points of
\(E_0\) with
\[
 F(x)\in r\mathbb Q^{*3}.                             \tag{4.1}
\]
Consequently finitely many cubeclasses of \(F(x)\), or a fixed finite
set of primes supporting those classes, give only finitely many
admissible targets. Any infinite Hasse subfamily on this slice must
use infinitely many cubeclasses and cannot keep their prime supports
inside a fixed finite set.

**Proof.** Impose \(F(x)=rt^3\). Substitution gives the curve
\[
 y^2=f_E\!\left(\frac{25865+rt^3}{128}\right),\qquad
 f_E(x)=x^3-10656x^2+27091730x.                        \tag{4.2}
\]
The polynomial \(f_E\) has three distinct roots. None equals
\(x_0=25865/128\), since
\[
 f_E(x_0)=\frac{25865}{128}
                \left(\frac{639735}{128}\right)^2
           =\frac{10585532408369625}{2097152}\ne0.      \tag{4.3}
\]
For \(r\ne0\), each root of \(f_E\) therefore has three distinct
nonzero preimages under \(t\mapsto(25865+rt^3)/128\).
Thus the right side of (4.2) is squarefree of degree nine. Its smooth
projective hyperelliptic completion is geometrically integral and
has genus \((9-1)/2=4\). Every rational elliptic point satisfying
(4.1) lifts to a rational point of (4.2). Faltings's theorem makes
that set finite.

Take a finite union for finitely many cubeclasses. A rational
cubeclass supported on a fixed finite set \(S\) has a representative
\(\prod_{p\in S}p^{e_p}\), \(e_p\in\{0,1,2\}\); signs are cubes.
There are at most \(3^{|S|}\) such classes, proving the last claims.
\(\square\)

For full field fibres this also rules out keeping the ramified-prime
support of the conductor \(f\) fixed. Indeed ELV2 would then force the
cubeclass of \(F(x)\) to be supported on that finite set together with
\(\{2,3,5,42649\}\). ELV4 gives only finitely many such points.
Thus any infinite Hasse family on this slice must introduce infinitely
many conductor primes; fixing one ray modulus cannot supply it.

The exact finiteness input is Faltings, *Endlichkeitssatze fur
abelsche Varietaten uber Zahlkorpern*, Inventiones Mathematicae
**73** (1983), Satz 7, printed p.365:
[primary scan](https://math.uchicago.edu/~drinfeld/Deligne%27s_conjecture_Manin_conf/Faltings_argument/Faltings.pdf),
[retained copy](sources/faltings-1983-endlichkeitssaetze.pdf).
It is the same imported theorem used in the elliptic arithmetic note.
No effective list of the points in (4.1), or uniform bound in \(r\),
is asserted.

## 5. What remains and what is checked

**Corollary ELV5 (a complete integer test for rational field fibres).**
Write a full field point in reduced integral coordinates
\[
 x=a/c^2,\qquad T=b/c^3,\qquad c>0,\qquad
 \gcd(a,c)=\gcd(b,c)=1,
 \qquad N=128a-25865c^2.
\]
Its fibre is Hasse-failing if and only if all the following hold:

1. `v_p(c)` is divisible by three for every prime `p=2 mod3`,
   including two and five.
2. `v_p(N)` is divisible by three for every prime `p=2 mod3`,
   `p>=11`.
3. If `3` does not divide `c`, then either `a=2 mod3`, or
   `a=1 mod3` and `27` divides `b`. If `3` divides `c`, put
   `n=v_3(c)`; then `3|n` and
   \[
    -(a-3552c^2)(c/3^n)b^{-1}=1\text{ or }-1\pmod9.
   \]
   Here `b` is a unit modulo nine, so its inverse is defined.

**Proof.** The integral monic elliptic model gives the displayed
coordinate form by the same negative-valuation argument as ELV3.
For either `p=2` or `p=5`, ETF2--ETF3 say exactly that a denominator
formal depth must be divisible by three; the proven rational
two-adic image supplies the hypotheses at two. Thus these two
conditions, together with ELV3, are precisely items 1--2.

The finite saturation certificate ETC3a determines the rational
three-adic closure as `E(Q_3)[2]+E_2(3)`. In its finite residue discs,
ETC1 consequently permits every point with `a=2 mod3`, permits
`a=1 mod3` exactly when `v_3(T)>=3`, and forbids `a=0 mod3`.
If `3|c`, the point is in the identity disc, of depth `n=v_3(c)`.
The exact formula
\[
 t=-X/T=-(a-3552c^2)c/b
\]
turns ETC1's cube condition on `t` into item 3. ELV1 now proves
both directions. All tests involve integer prime valuations and
congruences; no cubic-field maximal-order computation is required.
\(\square\)

The finite-place sieve in
[ELS1](ELLIPTIC_EXCEPTION_LOCAL_SIEVE.md) supplies infinitely many
points passing any fixed finite collection of local tests. ELV1
identifies the remaining simultaneous condition as the rational
valuation sieve (1.3), and ELV4 rules out fixing finitely many
cubeclasses to obtain that infinitude. Neither result controls all
new primes as the elliptic point varies. Infinitely many everywhere
locally soluble fibres on this elliptic slice remain unproved.

Run from the repository root:

```sh
python3 research/arithmetic/verify_elliptic_exception_explicit_local_criterion.py --write-receipt
```

The [checker](verify_elliptic_exception_explicit_local_criterion.py)
uses exact standard-library rational and polynomial arithmetic.
It verifies the translation, discriminant and coordinate identities,
the exceptional constants, the seed factorization, and the squarefree
degree-nine calculation. The bounded controls cover all primes
\(7\le p\le101\), units \(1\le u\le\min(p-1,8)\), and nonzero
\(-6\le v_p(z)\le6\); all nonzero unit residues for those primes
and for \(1889\); and the integer equivalence for
\(-40\le a\le40\), \(1\le c\le16\), \(\gcd(a,c)=1\), at all primes
\(11\le p\le59\), \(p\equiv2\pmod3\). Four fixed cubeclass
representatives control the degree-nine calculation, while the written
preimage argument proves squarefreeness for every nonzero class.
These are controls of the unrestricted proofs, not finite proofs of
all local fields or rational elliptic points. The
[receipt](elliptic_exception_explicit_local_criterion_checks.json)
records source hashes and coverage. The local-field theorems, RDF
classification and Faltings finiteness remain explicit proof inputs.
No formal verification, independent whole-proof replay or external
review is claimed.
