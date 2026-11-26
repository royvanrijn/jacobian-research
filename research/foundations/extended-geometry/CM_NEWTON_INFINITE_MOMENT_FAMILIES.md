# Infinite sphere-face exclusions from Frobenius coefficients

There are two distinct ways to lift a supporting-face obstruction to
unbounded degrees. A power of the angular face and a matching change of
height preserve the original connected cyclic cover. When that cover
has a nonvanishing Frobenius scalar on the required primes by the CM
inequality, this gives a characteristic-zero result with arbitrary
strictly lower-weight terms. Alternatively an iterated Frobenius at
one fixed prime can certify new connected covers of unbounded order.
In that second construction the observed Frobenius line can be
invertible even when the primitive Jacobian factor is nonordinary.

## Isolating the coefficient with a prime power

Use the sphere normal form and literal Laurent chart of the
[face criterion](CM_NEWTON_POLYGONS_AND_MOMENT_OBSTRUCTIONS.md).
Suppose a phase \(f\) has highest literal weight \(j+\rho i=n/q\),
with angular polynomial \(z^{-b}P(z^q)\) on the equality line.
Assume \(\rho>1\), \(a,b>0\), \(q\deg P=a+b\), and \(P(0)\ne0\). Let
\(n\mid(p^r-1)\) at an odd prime where the coefficients are integral,
and write \(k=(p^r-1)/n\).
At moment order \(qk\) the angular constant term has height at most
\(p^r-1\). Integration divides an even-height coefficient of
height \(j\) by \(j+1\). For \(0\le j<p^r-1\), the number \(j+1\) is
divisible by at most \(p^{r-1}\). Therefore multiplication of the
moment by \(p^r\) and reduction gives
\[
 p^r\mathcal A(f^{qk})\equiv
 [z^{bk}]P(z)^{qk}\pmod p.
\]
The endpoints of the weight line need no assumption on smaller
terms. All products of the top height use only equality nodes.

Assume the matching cover \(w^n=F(z)=z^{n-b}P(z)^q\) is connected,
with reduced nonzero exponents at both endpoints and the roots
of \(P\). At good reduction the coherent character-one line has
Čech generator \(w/z\). The \(r\)-fold iterate of Frobenius maps
this to \(w F(z)^k/z^{p^r}\); projection to the generator gives
the same coefficient. The hypothesis \(p^r\equiv1\bmod n\)
suffices for that iterate to preserve the line. Intermediate
Frobenius maps can permute characters. The calculation uses the
standard cyclic-cover eigensheaves and Hasse–Witt coefficient
formula. [Bouw, §§4–5](https://doi.org/10.1023/A:1017513122376)

## A uniform infinite family from each arithmetic face

Suppose a base face with data \((n,q,\rho,a,b,P)\) as above is
regular on the sphere and its connected cover has nonzero
character-one Frobenius scalar at every good absolute degree-one
prime outside a finite set after finite extension. Assume this
arithmetic statement holds for every number-field specialization
of the stated branch configuration. The simultaneous CM theorem
gives such a statement for a single mixed primitive pair with
multiplicities \((1,h-1)\) and the other pairs definite. It also
gives it for rank two with at most two mixed pairs.
[CM theorem](../../../unpublished/cm-newton-polygon-constraints/main.tex)

For any integer \(L\ge1\), let a new phase have highest
weight \(Ln/q\) for the same \(\rho\), with angular polynomial
\[
 z^{-Lb}P(z^q)^L.
\]
Its equality monomials are regular liftable, since their
angle-height nodes are sums of \(L\) base nodes. There is no
restriction on the smaller weights. After adjoining
\(\zeta_{nL}\), absolute degree-one primes satisfy
\(p=nLk+1\). At moment order \(qk\) the coefficient isolated by
the prime is
\[
 [z^{Lbk}] P(z)^{Lqk}
 =
 [z^{p-1}]\bigl(z^{n-b}P(z)^q\bigr)^{Lk},
\]
which is precisely the scalar of the original connected cover
at that prime. The finite extension for its ordinarity and the
cyclotomic extension can be taken together, because the arithmetic
conclusion applies at all sufficiently good absolute degree-one
primes of the chosen extension and persists over further extensions.
This excludes all-order moment vanishing. For coefficients in an
arbitrary characteristic-zero field, specialize a finitely generated
algebra to a number field after inverting the nonzero coefficients,
root differences, and any denominators needed to preserve the face
and branch configuration. The rational moment identities would
persist under this specialization.

The leading negative endpoint and the bound on the normal degree
give the exact least sphere degree \(Ld\) whenever the base lift has
degree \(d\) and the smaller-weight terms have normal degree at
most \(Ld\). The base examples in degrees six and eight thus give
infinite families at degrees \(6L\) and \(8L\). The symmetric
order-nine base of degree 13 likewise gives an infinite family
at degrees \(13L\). In this propagation the connected cover
remains the base cover. A formal equation of degree \(Ln\) with
right side \(F(z)^L\) is disconnected for \(L>1\).

## Selective certification by a split trace norm

The rank-two norm test can force a particular line at a split prime
while permitting, and even verifying, another nonordinary pair.
Suppose the relative rank is two, the field degree is \(2e\), and
\(b\) of the conjugate pairs are balanced. If the relative
Frobenius trace norm is nonzero and \(p>4^e\), the local bound
\(v_p(N\tau)\ge e-b+2B\) for \(B\) nonordinary pairs, together
with \(|N\tau|\le4^e p^e\), gives \(B\le\lfloor b/2\rfloor\).
At a good absolute degree-one prime in the simultaneous CM theorem
one can instead use the strict bound \(B\le\lfloor(b-1)/2\rfloor\).
In particular for three balanced pairs both bounds allow only one
bad pair. If a different pair is known to be bad, the selected
line must be invertible. In the local large-prime version, any
ordinary pair also certifies \(N\tau\ne0\).

For an explicit example take
\[
 p=73,\qquad P(z)=(z-1)^2(z-5)^3,\qquad
 w^9=z^5P(z).
\]
The primitive order-nine factor has three balanced pairs.
At this split prime the character-two Hasse–Witt scalar is zero
and the one in character four is \(22\). Thus one pair is bad,
another is ordinary, and the Weil comparison at \(73>64\)
forces the remaining pair, including character one, to be ordinary.
The exact scalar in character one is \(11\), but its computation
is unnecessary for the nonvanishing argument. The bad pair
\((2,7)\) makes the primitive sixfold itself nonordinary.

For any \(L\ge1\) coprime to 73, choose \(r\) with
\(73^r\equiv1\bmod 9L\). A phase with highest weight \(j+2i=9L\)
and angular polynomial \(z^{-4L}P(z)^L\) there, with arbitrary
strictly lower terms integral at 73, has a nonzero moment of order
\(m=(73^r-1)/(9L)\). The prime-power residue at the highest
angular-constant height is the \(r\)-fold scalar on the same
connected split order-nine cover, equal to \(11^r\). The least
degree of the lift is \(13L\). This result is local at the stated
prime. It leaves the characteristic-zero families with arbitrary
denominators to the uniform arithmetic argument when its hypotheses
hold.

## An infinite selective family with connected covers

Fix \(p=37\) and
\[
 P(z)=(z-1)^2(z-3)^3,\qquad N_r=\frac{37^r-1}{4},
 \qquad r\ge2 .
\]
Consider any rational phase whose literal support has greatest
weight \(j+2i=N_r\), with angular polynomial \(z^{-4}P(z)\)
on that line. More generally its angular polynomial can have
coefficients integral at 37 and reduce to this polynomial.
All strictly lower weights are allowed with coefficients integral
at 37. The same proof works over a number field when integrality
and the specified reductions hold at a prime above 37.

The fourth moment has angular-constant height at most
\(4N_r=37^r-1\). Prime-power isolation gives
\[
 37^r\mathcal A(f^4)\equiv
 [z^{16}](z-1)^8(z-3)^{12}
 \equiv32\pmod {37}.
\]
So the fourth moment is nonzero. The lift of the most negative
angle is a nonzero multiple of \(y^4t^{N_r}\).
The least sphere degree is \(N_r+4\) when all smaller-weight
normal terms have degree at most \(N_r+4\).

For the displayed \(P\) the connected cover has genuinely increasing
order. Its good reduction at 37 has equation
\[
 w^{N_r}=z^{N_r-4}P(z)
\]
and branch exponents \(N_r-4,2,3,N_r-1\) at \(0,1,3,\infty\).
Their sum is \(2N_r\); they are nonzero and reduced, and their
greatest common divisor with \(N_r\) is one. The primitive factor
has relative rank two over \(\mathbf Q(\zeta_{N_r})\). At
37 its character-one \(r\)-fold Frobenius scalar is \(32\).

This primitive factor is nonordinary at the certifying prime.
Writing \(x=i/N_r\) for a primitive character \(i\), its coherent
dimension is
\[
 d_i=1+\lfloor4x\rfloor-\lfloor2x\rfloor-\lfloor3x\rfloor.
\]
Indeed the fractional-part formula has exponents
\((-4,2,3,-1)\) modulo \(N_r\) and none of the four fractional
parts vanish at a unit character. Thus \(d_i=2\) for
\(1/4<x<1/3\), and \(d_i=1\) for \(0<x<1/4\).
Take \(i=N_r/3-1\), a unit as \(9\mid N_r\) and
\(\gcd(N_r/3-1,N_r)=1\). Its ratio lies strictly between
\(1/4\) and \(1/3\). On the same Frobenius orbit,
\[
 37^{r-1}i\equiv N_r/3-37^{r-1}\pmod{N_r}.
\]
This representative has ratio
\(1/3-4/(37-37^{-(r-1)})\) strictly between zero and \(1/4\).
Hence a Frobenius orbit in primitive coherent cohomology contains
components of dimensions two and one. Frobenius cannot be invertible
on that primitive factor. The isolated moment only uses the
character-one line, on which the \(r\)-fold map is invertible.
This family is a local Frobenius-coefficient application and uses
no assertion of ordinarity or prime density.

## What the split CM inequalities alone imply for the selected line

The raw weight-face construction supplies one coherent line,
with mixed signature \((1,h-1)\) on the selected pair. At rank two,
making just that pair nonordinary is numerically compatible with the
simultaneous bound when \(b\ge3\). In rank \(h\ge3\), set the
slopes on its dimension-one side to
\(s=(0^{h-2},1/2,1/2)\) and use complementary slopes for the dual
component. Its contribution to the rewritten bound is
\[
 D_j(s)=\left\lfloor\sum_{k=h-j+1}^h s_k\right\rfloor
       -\left\lceil\sum_{k=1}^j s_k\right\rceil .
\]
This is nonnegative for every \(1\le j<h\) and is zero for \(j=1\).
For example \(D_{h-j}=D_j\), so it suffices to inspect
\(j\le h/2\), where the lower sum is zero and the upper sum is
nonnegative. If another mixed pair stays ordinary it contributes at
least one for every \(j\). Thus this numerical obstruction remains
even if several singleton pairs are selected, provided some mixed
pair lies outside the selection.

A deck normalizer transports the selected zero by multiplication
on the primitive character group. Its action on conjugate pairs
has equal-size orbits, because
\((\mathbf Z/n)^\times/\{\pm1\}\) acts regularly by multiplication.
The balanced pairs form a union of these orbits. At rank two,
when the orbit size violates the bad-pair bound, every balanced
orbit is forced ordinary. Thus the uniform orbit argument for
primitive cyclic characters gives full ordinarity of that factor
whenever it forces the selected scalar nonzero.

The failure of a generic implication is visible even at a split
prime satisfying the strict relative trace-norm bound. For the
order-nine family \(w^9=z^5(z-1)^2(z-\lambda)^3\), take
\(p=73\) and \(\lambda=6\). The primitive coherent
scalars for characters \(1,2,4,5,7,8\) are
\[
 (0,22,30,30,22,0)\pmod{73}.
\]
The selected pair \((1,8)\) is the only nonordinary pair.
At its two split places the relative trace valuations are equal.
Indeed the dual Frobenius spectrum is obtained by replacing each
\(\lambda_i\) by \(p/\lambda_i\), and for rank two this sends the
trace \(\tau\) to \(p\tau/\det(\pi)\). The determinant has
valuation one in the selected pair. The two valuations are at least
one by its two slopes \(1/2\) and integrality. The other four place
traces are units. If the two bad valuations were at least two,
the nonzero \(E\)-norm would be divisible by \(73^4\), violating
the Weil bound \(64\cdot73^3<73^4\). Hence the exact norm valuation
is two, strictly below \(e=3\), and the selected line still vanishes.
This is an obstruction to a pointwise argument from the split
inequality alone. It places no restriction on the possible
density of this pattern for the fixed curve.

## Linked Frobenius failures on full sphere leading faces

Let \(n>2\) and let \(h(U,V)\) be a binary form of degree \(2n\)
over a characteristic-zero field. Its distinct geometric roots
have nonzero multiplicities \(a_i<n\), with
\(\gcd(n,a_1,\ldots,a_s)=1\). Use the normalized connected cover
\(W^n=h(z,1)\) and its primitive Jacobian factor. Suppose that deck
normalizers transport the character-one line to one character in
every mixed primitive conjugate pair. Then for every \(L\ge1\), a
sphere phase with leading binary form \(\gamma h^L\),
\(\gamma\ne0\), has infinitely many nonzero power moments. Every
lower sphere-degree term is allowed.

In the convention where \(W\) has character one, the coherent
dimension in a primitive character \(i\) is
\(-1+\sum_j\{ia_j/n\}\). The dimensions of a conjugate pair sum to
\(s-2\). In particular character one has dimension one. The primitive
factor carries \(\mathbb Q(\zeta_n)\) with Rosati acting by
conjugation, since the deck action preserves the Jacobian
polarization. The dimension formula is Serre dual to the
holomorphic differential formula of
[Moonen, Lemma 2.7 and equation (3.2.2)](https://arxiv.org/pdf/1006.3400).

Work first over a number field. Include the roots of unity of order
\(nL\), the deck normalizers, and the extension of the
[simultaneous CM bound](../../../unpublished/cm-newton-polygon-constraints/main.tex)
in one finite base extension. For a good absolute degree-one prime
of characteristic \(p\) in this extension put \(m=(p-1)/(nL)\)
and \(k=Lm\). The
[full sphere leading-form extraction](SPHERICAL_LEADING_BINARY_FORM_MOMENT_OBSTRUCTIONS.md#3-a-leading-coefficient-detected-at-good-primes)
at weight \(i+j\) gives
\[
       p\,\Lambda(f^m)\equiv
           \gamma^m [z^{p-1}]h(z,1)^k\pmod{\mathfrak p}.
\]
The only interval denominator divisible by \(p\) is at height \(p-1\),
and that height uses only the leading face. This applies to every
lower sphere-degree perturbation, regardless of the heights of its
individual Laurent monomials. The conductor of the normalized
connected cover stays \(n\) for all \(L\).

On its character-one eigensheaf, \(W\) and \(z^{-2}W\) generate over
the finite and infinite charts, respectively. All finite branch
exponents are below \(n\), and the exponent at infinity is
\(2n-\deg h(z,1)\) in \([0,n-1]\). The Čech class of \(W/z\) spans
coherent \(H^1\) in character one. Raising it to the \(p\)-th power
and projecting back gives the scalar \([z^{p-1}]h(z,1)^k\).
Coherent \(H^1\) of the curve and of its Jacobian agree equivariantly,
so this scalar observes the corresponding primitive factor.

If the scalar were zero, the deck normalizers, defined over the
residue field at these primes, would transport it to a zero on a
one-dimensional coherent component of every mixed pair, by the
[normalizer transport argument](CM_NEWTON_POLYGONS_AND_MOMENT_OBSTRUCTIONS.md#transport-of-a-vanishing-line-by-a-deck-normalizer).
Each such pair would be nonordinary. Its dimension-one block has
slopes summing to one. A unit slope would make all the other slopes
zero and the pair ordinary. Thus the smallest slope of the dual
block is strictly positive, contributing at least one to the
rounded \(j=1\) sum. Each definite pair contributes one as well.
For \(e=\varphi(n)/2\) the sum is at least \(e\), in contradiction
with the CM bound \(e-1\). Hence the coefficient is nonzero at every
such prime outside the finite set excluding denominators, bad
reduction and the vanishing of \(\gamma\). These primes have
unbounded rational characteristics, so the certifying orders \(m\)
are unbounded. Over a general characteristic-zero field, assumed
eventual vanishing would survive specialization of the finitely
generated coefficient algebra to a number field after inverting
\(\gamma\), the branch differences and denominators. One may
include the branch coordinates and normalizer data in that algebra
to preserve their defining identities.

### The symmetric degree-nine phase

Put
\[
 H_9=t(-2y)^3(x-2t-2y)^2(x-t-2y)^3
 \quad\text{on}\quad t^2+2xy=1.
\]
Substitution of \((U^2,-V^2/2,UV)\) gives the leading binary form
\(h_9=UV^7(U-V)^4(U^2-UV+V^2)^3\). In particular any phase
\(f=\gamma H_9^L+R\) with \(\deg_{\mathrm{sphere}}R<9L\)
has exact sphere degree \(9L\). On \(W^9=F(z)\) write
\(F(z)=z(z-1)^4Q(z)^3\) and \(Q(z)=z^2-z+1\). The transformation
\(S(z)=1/(1-z)\) cycles \(0,1,\infty\) and fixes both roots of \(Q\).
The assignments \(1,4,7\) on the cycle and \(3,3\) on the fixed
points are multiplied by \(4\) modulo \(9\). Explicitly
\[
  \frac{F(S(z))}{F(z)^4}=\frac{1}{((1-z)^3Q(z))^9},
  \qquad
  (z,W)\longmapsto\left(S(z),\frac{W^4}{(1-z)^3Q(z)}\right).
\]
The lift transports the character-one line to characters \(4\)
and \(7\). The coherent dimensions in \(1,4,7\) are all one,
and in \(8,5,2\) they are all two. These three mixed \((1,2)\)
pairs meet the criterion above and prove infinite nonvanishing.

### A movable conductor-seven family

Let \(\zeta\) be a primitive cube root of unity and take \(a,b\ne0\)
with \(a^3\ne b^3\). The binary form
\[
 h_7(U,V)=\prod_{j=0}^2
     \bigl((U-\zeta^j aV)(U-\zeta^j bV)\bigr)^{2^j}
\]
has two disjoint three-point orbits of branch exponents \(1,2,4\)
under \(S(z)=\zeta z\). Set \(B_j(z)=(z-\zeta^ja)(z-\zeta^jb)\)
and \(F(z)=h_7(z,1)\). The identity
\(F(\zeta z)/F(z)^2=(\zeta^2/B_2(z))^7\) gives the lift
\((z,W)\mapsto(\zeta z,\zeta^2 W^2/B_2(z))\). The order-seven
curve has three primitive mixed pairs of type \((1,3)\). Their
singleton characters are \(1,2,4\), one orbit under the lift.
Consequently every sphere phase with leading binary form
\(\gamma h_7^L\), for \(\gamma\ne0\) and \(L\ge1\), has infinitely
many nonzero power moments, with arbitrary lower sphere-degree terms.
Every such leading form has a sphere lift by the degree-\(7L\)
normal basis.

### Two single-orbit examples

Let \(\zeta_k\) be a primitive \(k\)-th root of unity. In each row
take \(h_n(U,V)=\prod_{j=0}^{k-1}(U-\zeta_k^j V)^{a_j}\)
with the indicated ordered list of exponents.

| \(n\) | \(k\) | Deck multiplier \(u\) | \((a_0,\ldots,a_{k-1})\) | Mixed pair types |
| --- | --- | --- | --- | --- |
| 11 | 5 | 3 | \(1,3,9,5,4\) | Five of type \((1,2)\) |
| 21 | 6 | 2 | \(1,2,4,8,16,11\) | Six of type \((1,3)\) |

In both rows the exponents sum to \(2n\), and the exponent \(1\)
ensures connectedness. Rotation by \(\zeta_k\) multiplies
the exponents by \(u\) modulo \(n\), including across the last
entry. The \(k\) powers of \(u\) give distinct characters modulo
conjugation and cover all primitive pairs. The orbit starts
from coherent dimension one, so each paired dimension is
\(k-3\), as asserted in the table. For both rows this proves
infinite nonvanishing for every sphere phase with leading binary
form \(\gamma h_n^L\) with \(\gamma\ne0,\ L\ge1\), with arbitrary
lower sphere-degree terms.

### Balanced homogenization and a finite search criterion

For any of these families and every \(D\ge\deg_{\mathrm{sphere}}f\),
let \(H_D(f)\) denote the balanced homogenization of \(f\) using the
two-pair Hopf radius. Exact contraction gives
\(\mathcal E_2(H_D(f)^m)=(mD+1)!\Lambda(f^m)\). Since the two-pair
Image space is the kernel of \(\mathcal E_2\), each lift lies
outside its radical. This uses infinitely many nonzero pure
moments of each individual phase.
[Hopf normalization](HOPF_AUTOMORPHISM_ORBIT_SPECTRA.md#5-moment-and-original-image-interpretation)
and [Image-kernel identity](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md).

For fixed \(n\), a search with one nontrivial Möbius normalizer can
enumerate branch data directly. A finite branch set of size at least
three forces any Möbius symmetry of it to have finite order \(k\).
After a geometric conjugation the symmetry is
\(z\mapsto\zeta_k z\), with two fixed points and length-\(k\)
orbits elsewhere. At a fixed point a reduced exponent is either
zero or satisfies \((u-1)a=0\bmod n\), where \(u\) is the deck
multiplier. A full orbit has successive exponents
\(a,ua,\ldots,u^{k-1}a\bmod n\), all taken as nonzero reduced
representatives. Connectedness gives \(u^k=1\bmod n\).
The total exponent is \(2n\), so there are at most \(2n\) branch
points and \(k\le2n\). Thus the exponent data can be enumerated by
finite integer lists for each \(n\). Conversely, each such list
with total \(2n\) and gcd one with \(n\) is geometrically realizable
by disjoint orbits of the rotation, with the fixed-point exponents
as specified.

For each list compute the fractional-part dimensions and the
orbit of character one. Write \(\rho\) for the relative height,
\(m_i\) for the smaller dimension of pair \(i\), and \(c\) for
the number of pairs reached by the orbit. For the information
that those \(c\) pairs are nonordinary, the simultaneously least
excess at index \(j\) is
\[
    c\bigl(\mathbf 1_{j=1}+\mathbf 1_{j=\rho-1}\bigr).
\]
Indeed every reached pair has \(m_i=1\). Its largest nonordinary
endpoint multiplicities are \(z=\rho-2,\ o=0\) on the
dimension-one side. The
[least-polygon theorem](../../../unpublished/cm-newton-polygon-constraints/numerical.tex)
gives the two indicated contributions, and unreached pairs may
be ordinary. Reject the assumed scalar vanishing if this excess
exceeds \(\sum_i\min(j,\rho-j,m_i)-1\) for any \(1\le j<\rho\).
The all-definite case has ordinary reduction directly. Passing
the inequalities establishes numerical compatibility of the
listed failures. It gives no vanishing assertion for any actual
cover. Normalizers with several generators can be tested by the
same inequalities after their branch actions have been verified
and their joint character orbit computed.

## Exact replay

The [standard-library checker](verify_cm_newton_infinite_moment_families.py)
recomputes the order-nine scalar rows. It tests the first three
connected orders \(342,12663,468540\) by rational multiplication of
the normalized lifts with smaller terms added, and verifies the
different coherent dimensions in the indicated Frobenius orbit.
The internal replay passes and gives residue 32 in each case.
The symbolic congruence and orbit calculation establish the
unbounded statements.

    python3 -I verify_cm_newton_infinite_moment_families.py
