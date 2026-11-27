# Prime progressions with the same cubic tangent field

The quadratic field of a cubic tangent triangle is the
field of a point on a fixed elliptic double cover. The
three prime vertices mark the fibres over \(0,1,-1\).
Multiplication preserves or merges those fibres. Their
prime factorizations force the point to be primitive
in a trace coset. Two distinct genuine progressions in
one field therefore require a companion quadratic twist
of rank at least two.

## The fixed elliptic double cover

Write \(R(p,h)=3p^4-6p^2h^2-h^4\) and let
\((p-h,p,p+h)\) be distinct primes greater than \(3\), with
\(0<h<p\). By [CT1](ABC_CUBIC_TANGENT_NORM_RANK_AND_CONDUCTORS.md#the-tangent-identity-and-the-lost-cubic-direction)
the quadratic field is \(K=\mathbb Q(\sqrt d)\), where
\(d\) is the signed squarefree part of \(-R(p,h)\).
Put
\[
 f(t)=1+6t^2-3t^4,\qquad
 E:y^2=x^3+1,\qquad T=(-1,0),\qquad
 \tau(P)=\frac{x(P)+1}{y(P)} .
\]
Let \(\sigma\) be the nontrivial involution of \(K\).

**Proposition FC1.** The smooth double cover
\(U^2=f(t)\) identifies over \(\mathbb Q\) with \(E\),
with \(t=\tau(P)\) and deck involution \(P\mapsto T-P\).
The admissible rational parameters with field \(K\) are
the finite \(\tau\)-values of the trace coset
\[
 S(K)=\{P\in E(K): P+\sigma P=T\}.
\]
For fields of genuine tangent progressions,
\(K\ne\mathbb Q(\sqrt{-3})\). Every finite
\(t\) from \(S(K)\) in such a field has a reduced
denominator divisible by \(6\). The anti-invariant
subgroup \(A(K)=\ker(1+\sigma)\) identifies with the
rational points of \(E^{(d)}:y^2=x^3+d^3\).

For \(t\ne0\) the identification and its inverse are
\[
 x=\frac{U+1+t^2}{2t^2},\qquad y=\frac{x+1}t,\qquad
 U=2t^2 x-1-t^2.
\]
In particular
\(\tau^2=(x+1)/(x^2-x+1)\).
The fibre at zero consists of \(O,T\), so its sum is
\(T\) and the deck involution is \(T-P\).
The rational parameters with \(f(t)\) in the nontrivial
squareclass \(d\) give the indicated trace. There
are no rational roots of \(f\), since \(-f(1/t)t^4\)
is Eisenstein at \(3\) up to sign. The two poles of
\(\tau\) are the other points of order two and are
defined over \(\mathbb Q(\sqrt{-3})\).
The twist identification on an affine anti-invariant
point is \((x,y)\mapsto(dx,d\sqrt d\,y)\).

For \(t\in\mathbb Z_3\), \(f(t)\equiv1\bmod3\) is a
local square. For \(t\in2\mathbb Z_2\), \(f(t)\) is
\(1\bmod8\). For \(t\in\mathbb Z_2^\times\), write
\(t^2=1+8j\). Then \(f(t)=4(1-48j^2)\) is again a
local square. If \(v_2(t)<0\) or \(v_3(t)<0\), the
factorization
\[
 f(t)=-3t^4\left(1-\frac2{t^2}-\frac1{3t^4}\right)
\]
shows that its local squareclass is \(-3\) at the
respective prime. The bracket is \(1\bmod8\) at \(2\)
and \(1\bmod3\) at \(3\). In a genuine progression
\(6\mid h\) and \(p\) is a unit at \(2,3\).
Consequently those local classes are \(-3\), and
another parameter of the same field must have
negative valuation at both primes.

For the exceptional field \(d=-3\), the twist
\(y^2=x^3-27\) has rank zero by two-isogeny descent.
After translating \(x\) by \(3\), the coefficients
\((a,b)\) in \(y^2=x(x^2+ax+b)\) are \((9,27)\)
and the isogenous coefficients are \((-18,-27)\).
The two squareclass images both have size two.
The first is \(\{1,3\}\) since the quadratic factor
is everywhere positive. On the second curve the
two-torsion gives class \(-3\). The candidate class
\(-1\) would require coprime integers \(u,v\) with
\[
 n^2=-u^4-18u^2v^2+27v^4 .
\]
Reduction modulo \(3\) forces \(3\mid u\), after
which the right side has valuation \(3\). This rules
out \(-1\) and therefore \(3\). The rank formula
\(2^{r+2}=|\alpha(E)|\,|\alpha(E')|\) gives \(r=0\).
See [Dixit–Murty–Pathak, Section 2](https://www.imsc.res.in/~anupdixit/rankproblem.pdf)
for the two-isogeny formula and its quartics.

In any quadratic extension the reductions of \(E\)
at \(5,7\) have residue degrees at most two, with
point counts \(6,36\) and \(12,48\), respectively.
They exclude all torsion of prime order at least
five. The order-three division polynomial is
\(3x(x^3+4)\), so the only order-three points over
a quadratic field are rational and hence cannot
lie in \(A(K)\). Its two-torsion is rational since
negation fixes it, and a point of order four in
\(A(K)\) would halve \(T\). Such a point would be
a branch point of \(\tau\) with a degree-four
parameter. Thus \(A(K)_{\rm tors}=\{O,T\}\).
For \(K=\mathbb Q(\sqrt{-3})\) the two nonrational
points of order two already have trace \(T\).
As \(A(K)=\{O,T\}\) in this case, they exhaust the
trace coset and both give \(t=\infty\). This
excludes the exceptional field for finite witnesses.

## The three marked fibres force primitivity

Put \(H(K)=\{P\in E(K):P+\sigma P\in\{O,T\}\}\).
When the trace coset is nonempty \(A(K)\) has index two
in \(H(K)\), so the rank is that of \(E^{(d)}(\mathbb Q)\).

**Theorem FC2.** If \(P\in S(K)\) gives a genuine
prime progression, its image in \(H(K)/H(K)_{\rm tors}\)
is primitive. Points of \(S(K)\) giving two different
positive ratios \(p/h>1\) of genuine progressions
have independent images over \(\mathbb Q\) in this
free group. In particular rank one on \(E^{(d)}\)
allows at most one such progression in the field.
Each progression gives a nontorsion point on
\(E^{(d)}\) by \(P\mapsto2P-T\) and the twist map.

First \(H(K)_{\rm tors}=\{O,T\}\) for \(K\) of FC1,
since there are no further points of order two
over \(K\), no halves of \(T\) in a quadratic
field, and the odd torsion of \(H(K)\) has trace
zero and therefore lies in \(A(K)\). In particular
the free image of \(P\) is nonzero. If it is
divisible by \(n\), then \(P=[n]G\) up to a translate
by \(T\), with \(n\) odd and \(G\in S(K)\).
Translating by \(T\) or negating changes \(\tau\)
by sign, so we may take
\(\tau(G)=a/b>0\) and
\(|\tau([n]G)|=p/h\), with \(b>0\) and
\(\gcd(a,b)=1\). FC1 gives \(6\mid b\). Thus
\(|a-b|,a,a+b\) are nonzero, pairwise coprime,
and prime to \(6\).

The fibres of \(\tau\) at \(0,1,-1\) on the
good model at every prime \(q>3\) consist,
respectively, of
\[
 \{O,T\},\qquad \{D,T-D\},\qquad \{-D,T+D\},
 \qquad D=(0,1),\quad [3]D=O .
\]
Indeed the nonzero fibres displayed have affine
\(x\)-coordinates \(0,2\). The map \(\tau\)
extends to the good integral model. At \(O\) use
\([X+Z:Y]\) on the plane cubic, and near \(T\) use
\(y/(x^2-x+1)\). It is étale at the six displayed
sections in residue characteristic \(q>3\), since
the values \(f(0),f(1),f(-1)\) are \(1,4,4\).
Hence if \(q^j\) divides one of the three
source markings, \(G\) is congruent modulo
\(q^j\) to its section in the corresponding
fibre, by étaleness. Multiplication by \(n\)
preserves the congruence on the good model.
At the target the denominator \(h\) is a unit
at this \(q\), since the reduced value of
\(\tau([n]G)\) is finite.

For \(3\mid n\) all three marked fibres map
to the zero fibre. At least two source
markings exceed \(1\) in absolute value since
\(b\ge6\). Their distinct prime divisors
would both divide the prime \(p\).
For \(3\nmid n\) the zero fibre is preserved
and the two outer fibres are preserved or
exchanged by a common permutation, allowing
also the sign of the target. Each nonunit
source marking must therefore equal the
corresponding target prime in absolute
value, since all its powers divide that
prime. If \(a=1\), the two outer target
primes would be \(b-1,b+1\) with midpoint
\(b\), divisible by \(6\). If \(a>1\) and
\(|a-b|=1\), then \(a=p\) and the remaining
outer prime \(2p\pm1\) must be less than
\(2p\). The only possibility would yield
the other target vertex \(1\).
All three source markings therefore exceed
\(1\) in absolute value. The identities
\(a=p\) and \(a+b=p+h\) give \(b=h\), since
\(a+b>a\). The source ratio is the same
as the positive target ratio, so
\([n]G=\pm G\) modulo \(T\), using the
two-point fibre and its negative. Nontorsion
implies \(n=\pm1\).

To prove independence for two progressions,
a common rational line would contain a
primitive lattice vector on which both
images are \(\pm1\) multiples. Their
points would then agree up to negation
and translation by \(T\), giving the
same positive ratio. Reduction to
primitivity uses all three markings
even though the center already
detects divisibility along a progression
that is known to be prime.

The divisibility mechanism is related to
[Everest–Ingram–Mahé–Stevens, Theorems 1.4–1.5](https://arxiv.org/pdf/0712.2696).
Their bounds concern prime terms of elliptic
divisibility sequences under a magnifying
isogeny. Here the constraint concerns the
three numerator markings of a torsion-shifted
degree-two quotient. The three target primes
force primitivity exactly, without a height
threshold or a primitive-divisor estimate.

## Field-preserving maps of the parameter line

**Corollary FC3.** The nonconstant \(F\in\mathbb Q(t)\)
with \(f(F(t))/f(t)\in\mathbb Q(t)^{\times2}\)
are exactly the maps induced under \(\tau\) by
odd multiplications on \(E\), with an optional sign.
Their degrees are \(n^2\) for odd \(n\).
No such map sends a genuine progression to a
different genuine progression.

A square ratio lifts \(F\) to a rational self-map
of \(U^2=f(t)\) by multiplying \(U\) by its square
root. It is a morphism on the smooth projective
curve. Every nonconstant \(\mathbb Q\)-morphism
\(E\to E\) is \([n]+B\) for \(B\in E(\mathbb Q)\).
Indeed after subtracting the image of \(O\)
it is an endomorphism. The complex analytic
multiplier of a \(\mathbb Q\)-endomorphism is
rational, from its action on a rational
differential. A rational multiplier preserving
the period lattice is an integer. Commutation
with \(P\mapsto T-P\) gives
\(2B=(1-n)T\). There is no rational half of
\(T\), so \(n\) is odd and \(B=O\) or \(T\).
Conversely those maps descend and commute
with the deck involution, so the new \(U\)
is \(U\) times a rational function of \(t\).
The degree of \([n]\) is \(n^2\), hence so
is the induced degree. The prime exclusion
follows by FC2, since one trace-coset
point is a multiple of the other up to
\(T\). These are elliptic quotient maps
of the type described in
[Milnor, Definition 2.1](https://people.math.harvard.edu/~demarco/Math295/Milnor_Lattes.pdf).

In the prime fixture \((p,h)=(17,6)\) with
\(d=-186843\), the positive ratios for \([3]P\)
and \([5]P\) are
\[
 \frac{121647372209}{82553876370},\qquad
 \frac{6481558012053794672265292343537}
      {1613150365317644553963830999454}.
\]
The first center is divisible by \(11\cdot17\cdot23\).
In the second ratio the three markings inherit
\(23,17,11\), respectively, and the center is
again composite. Both have the same quadratic
field, positive markings and denominator
divisible by \(6\). The exact arithmetic is
reproducible with
[the Sage checker](verify_cubic_tangent_field_collisions.py).

In fact each field containing a genuine
progression has infinitely many relaxed
positive ratios \(t>1\) given by odd
multiples of that point. For a real field,
the point is nontorsion on the connected
compact real Lie group \(E(\mathbb R)\).
Even multiples are dense in its circle,
so the odd multiples return infinitely
often to an open arc at \(P\) with
\(\tau>1\). For an imaginary field, use
the compact subgroup of \(E(\mathbb C)\)
defined by \(Q+\overline Q=O\), with
one-dimensional identity component.
The point \(2P\) lies in this subgroup.
A nonzero multiple of \(2P\) in the
identity component is nontorsion and
generates a dense subgroup there.
Adding its multiples to \(P\) gives
the required returns in the trace
coset. Each finite \(\tau\)-fibre has
only two points. FC2 allows only the
original prime progression on this
cyclic line.

## Fixed projective relations and adjacent centers

**Proposition FC4.** Fix a rational Möbius
transformation. It relates only finitely
many pairs of distinct genuine progressions
in the same quadratic field. In particular
for each fixed rational \(s\ne0\) there are
only finitely many same-field pairs with
the same step \(h\) and centers \(p,p+sh\).
Such a pair is impossible whenever
\[
 p\ge \max(4h,\,2|s|h,\,10|s|h^3).
\]

The quartic branch set of \(\tau\) has Galois
group \(D_4\) of order eight. Indeed its
roots are \(\pm\sqrt{1+2/\sqrt3}\) and
\(\pm\sqrt{1-2/\sqrt3}\). The first
positive root generates a real field
of degree four by irreducibility. The
imaginary root generates a further
quadratic extension since its square
already lies in the real field. The
group permutes the opposite pairs, so
this is the full dihedral group. A rational
Möbius map commuting with Galois on
these roots and preserving the branch
set must induce an element of the
centralizer of \(D_4\) in \(S_4\).
This centralizer has order two, generated
by simultaneous sign reversal. Three
values determine a Möbius map, so the
only possibilities preserving the set
are \(t\) and \(-t\). They cannot relate
different ratios \(>1\).
All other maps have disjoint branch
sets for \(f(t)\) and \(f(F(t))\), since
their intersection is a Galois-invariant
subset of the four conjugates.
The equality of field squareclasses
then cuts out the double cover with
equation \(Y^2=f(t)f(F(t))\) and eight
branch points, hence genus three.
[Faltings, Satz 7](sources/faltings-1983-endlichkeitssaetze.pdf)
gives finitely many rational parameters.
The fractions \(p/h\) for genuine
progressions are reduced, so each
rational parameter determines at most
one positive \((p,h)\). This finiteness
holds for each fixed transformation.

For the elementary exclusion, take
\(A=p^2-h^2\), \(B=(p+sh)^2-h^2\) and
\(N=3AB-4h^4\). Then
\[
 R(p,h)R(p+sh,h)=N^2-12h^4(A-B)^2 .
\]
Under the stated bounds,
\(A\ge15p^2/16\), \(B\ge3p^2/16\), and
\(N>p^4/2\). The subtracted term is
positive and at most \(75s^2h^6p^2\),
which is at most \(3p^4/4\).
The product lies strictly between
\((N-1)^2\) and \(N^2\). For a prime
pair the quantities are integers.
The squareclass equality would make
the integer product a square.
If both signed quartic values are
already squarefree, any same-step
collision is also excluded by strict
increase of \(R(p,h)\) for \(p>h\).

## Multiplicity uniform in the field

**Theorem FC5.** Let \(N_d(Z)\) count the
genuine progressions with center at most \(Z\)
in a given signed squarefree field \(d\), and
let \(r_d=\operatorname{rank}E^{(d)}(\mathbb Q)\).
For an absolute constant \(C\) one has
\[
 N_d(Z)\le
 \left(1+C\sqrt{\frac{\log(2Z)}{\log(2|d|)}}\right)^{r_d}.
\]
In particular \(N_d(Z)\le
\exp(O(\log Z/\log\log Z))\) uniformly over
the represented fields. For each sign of \(d\)
this gives
\[
 X^{1/2}\exp\left(-O\!\left(\frac{\log X}{\log\log X}\right)\right)
\]
distinct genuine tangent fields of absolute
discriminant at most \(X\), with no sieve on
the quartic square factor. The statement
also holds when \(h/6\) must be prime.

Use the canonical height normalized by
\(\widehat h(P)=h(x(P))/2+O_E(1)\), with
\(h\) the logarithmic absolute projective
height. The difference is uniformly bounded
on the fixed \(E\), by
[Silverman, Theorem VIII.9.3](https://www.math.ens.psl.eu/~obenoist/refs/Silverman.pdf).
The degree-two morphism in \(x\) for
\(\tau^2\) gives
\(h(x(P))=h(\tau(P))+O(1)=\log p+O(1)\).
In the height space of \(H(K)\) modulo
torsion, each counted progression thus
gives a vector of squared length at most
\(\tfrac12\log Z+O(1)\). Two vectors from
distinct progressions have nontorsion
difference in \(A(K)\).

For a nontorsion \(Q\in A(K)\), write
\(x(Q)=u/v\) in lowest terms. The
integer \(v(u^3+v^3)\) has squareclass \(d\),
so \(|d|\le 2\max(|u|,|v|)^4\).
It follows that
\(\widehat h(Q)\ge\tfrac18\log|d|-O(1)\).
This gives a lower bound
\(c\log(2|d|)\) with \(c>0\) for nonzero
height differences. The adjustment for
bounded \(|d|\) uses the positive minimum
of each of the finitely many
Mordell–Weil lattices. Packing disjoint
balls of half that separation in the
height space, contained in one ball
with radius enlarged by the same amount,
gives the first bound.

The translated twist has \((a,b)=(-3d,3d^2)\)
and its two-isogenous partner has
\(b'=-3d^2\). Both descent images are
supported on the signed squarefree
divisors of \(3d\), giving \(r_d\le2\omega(3d)\).
Put \(L=\log(2Z)\) and \(J=\log(2|d|)\).
Since \(|d|\ll Z^4\), \(J\ll L\).
When \(J\ge\sqrt L\), the bound
\(\omega(3d)\ll J/\log J\) and
\(J\log(1+C\sqrt{L/J})\ll\sqrt{JL}\)
give \(\log N_d(Z)\ll L/\log L\).
For \(J<\sqrt L\) one has the stronger
\(\log N_d(Z)\ll\sqrt L\log L\).
The implied constants are uniform.

Finally the prime-forms step of
[CT4](ABC_CUBIC_TANGENT_NORM_RANK_AND_CONDUCTORS.md#squarefree-dyadic-fields-and-rank-two-conductors)
gives \(\gg Z^2/(\log Z)^s\) ordered pairs
of either sign, for \(s=3\) or \(4\),
before its squarefree sieve. Divide by
the uniform field multiplicity and take
\(Z\) a suitable constant times \(X^{1/4}\).
The discriminant has absolute value at
most \(X\) since \(|d|\ll Z^4\) and
\(d\equiv1\bmod4\).
The conclusion has a weaker count than
the squarefree-family result CT4,
with the new control applying to all
square factors. It also supplies an
explicit rank input for collision
moments. A stronger sum over fields
would require bounds on the twist ranks
or the three marked numerators on the
higher-rank trace cosets.
