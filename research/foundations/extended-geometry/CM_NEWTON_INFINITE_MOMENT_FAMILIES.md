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
