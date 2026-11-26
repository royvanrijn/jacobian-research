# Exact complementary products determine three primitive rows at subpower cost

For a positive primitive row \(A+B=C\), use the two fractions
\(A/C\) and \(B/C=1-A/C\). The products of both fractions over
three rows define a genus-one curve. Cycling the rows is a
translation of exact order three. The resulting uniform
bounded-height count gives small triple-product fibres
for the exact records and also for their full binary digit
records at a fixed power of the radical modulus.

## Completion from two exact products

Write \(H(t)\) for the height of a rational number in lowest
terms. Fix \(K>0\). For \(u,v>0\) rational and \(X\ge2\), consider
ordered triples \((a,b,c)\) of rational numbers in \((0,1)\)
of height at most \(X^K\), satisfying
\[
                   abc=u,\qquad(1-a)(1-b)(1-c)=v .
                                                               \tag{EC1}
\]
The count is \(O_{K,\varepsilon}(X^\varepsilon)\), uniformly in
\(u,v\) for any \(\varepsilon>0\). A nonempty fibre has
\(H(u),H(v)\le X^{3K}\).

The affine curve in coordinates \(a,b\) is
\[
                 (1-a)(1-b)(ab-u)-vab=0.
                                                               \tag{EC2}
\]
On its open set where \(a,b,c=u/(ab)\) are finite and different
from 0 and 1, it is isomorphic to the two equations in (EC1)
with no positivity restriction. The gradients of the two
products there have coordinates \(u/a,u/b,u/c\) and
\(-v/(1-a),-v/(1-b),-v/(1-c)\). A singular point therefore
has \(a=b=c=z\), requiring \(z^3=u\) and \((1-z)^3=v\).
For positive real \(u,v\), such a \(z\) must be real in \((0,1)\).
In fact, if \(z\) is nonreal, the two nonreal summands \(z,1-z\)
are cube roots of positive numbers with opposite imaginary
parts. Their real parts are both negative, so they cannot add
to 1. In the real case \(u^{1/3}+v^{1/3}=1\). Hölder gives
\[
 \left(\prod_i x_i\right)^{1/3}+
 \left(\prod_i (1-x_i)\right)^{1/3}
          \le \left(\prod_i (x_i+1-x_i)\right)^{1/3}=1
 \qquad ((x_1,x_2,x_3)=(a,b,c)).
\]
Equality for \(0<x_i<1\) forces the \(x_i\) equal. Thus in the
case of an open-set singularity there is at most one admissible
triple.

In the remaining case the closure in
\(\mathbb P^1\times\mathbb P^1\) is smooth at the boundary as
well. Its homogeneous equation in \([a:p],[b:q]\) is
\[
            (p-a)(q-b)(ab-upq)-vabpq=0.
\]
The boundary consists of
\[
          (0,1),\ (1,0),\ (0,\infty),\ (\infty,0),\
          (1,\infty),\ (\infty,1).
\]
In particular \((\infty,\infty)\) does not occur because
the left side there equals 1 with \(a=b=1,p=q=0\).
At \((0,1)\) the derivative of (EC2) in \(b\) is \(u\).
In the chart near \(b=\infty\) the equation is
\((1-a)(q-1)(a-uq)-vaq=0\), whose derivative in \(a\)
at \((a,q)=(0,0),(1,0)\) is respectively \(-1,1\).
Symmetry checks the other three points. Smoothness implies
geometric irreducibility here. Two distinct components of
bidegrees \((d,e),(d',e')\) meet if \(de'+ed'>0\); disjoint
components would all be fibres of the same ruling, contrary
to total bidegree \((2,2)\). By adjunction the curve has genus
one, and its listed boundary points are rational.

Cycling the three letters extends to an automorphism
\(\sigma(a,b)=(b,u/(ab))\) of this smooth projective curve.
It has order three. Its only possible fixed points on the
open set would lie on the excluded diagonal. At the boundary
it cycles \((0,1),(1,\infty),(\infty,0)\) and cycles
\((1,0),(0,\infty),(\infty,1)\). Indeed the first two
steps of the first cycle and the steps from \((1,0)\) and
\((\infty,1)\) follow by taking limits of \(u/(ab)\), and
\(\sigma^3=1\) determines the others. Choosing a rational
boundary point as origin makes the curve an elliptic curve.
Write \(\sigma(P)=h(P)+Q\), with \(h\) an origin-preserving
automorphism. If \(h\ne 1\), the nonzero endomorphism
\(1-h\) is surjective over \(\overline{\mathbb Q}\), and a
point mapping to \(Q\) would be fixed. Consequently
\(\sigma\) is translation by a nonzero rational point of
exact order three.

### A model with polynomial height

For a point of (EC1), the variable \(Y=a(1-a)(b-c)\) obeys
\[
       Y^2=[a(1-a)+u(1-a)-va]^2-4ua(1-a)^2
           =u^2+f_1a+f_2a^2+f_3a^3+a^4.
\]
All \(f_i\) have polynomial rational height. The quadratic
discriminant identity comes from \(bc=u/a\) and
\(b+c=1+u/a-v/(1-a)\). Conversely \(a,Y\) recover that
sum and \(b-c\), hence both \(b,c\) at every admissible
point. Put \(r=f_1/(2u)\) and
\(\xi=(Y+u+ra)/a^2\). Then \(Y=\xi a^2-u-ra\), so
\[
 (\xi^2-1)a^2-(2r\xi+f_3)a+r^2-2u\xi-f_2=0.
\]
Taking the discriminant in \(a\) gives the model
\[
 Z^2=(2r\xi+f_3)^2-4(\xi^2-1)(r^2-2u\xi-f_2),
 \qquad Z=2(\xi^2-1)a-(2r\xi+f_3).
\]
Its cubic leading coefficient is \(8u\ne0\). The changes are
birational, with only boundedly many exceptional points
when \(\xi^2=1\), since the original curve is irreducible.
In particular the cubic model is nonsingular. The change
\((\xi,Z)\mapsto (8u\xi,8uZ)\) makes the cubic monic,
of polynomial coefficient height, and the first coordinate
of every image of an admissible point has polynomial height.
Its elliptic curve has a rational point of order three.
An isomorphism taking the chosen origin to the
Weierstrass origin is obtained by rational translation if
necessary, so existence of the torsion transfers. Apply
Dujella's Corollary 1.3 with a height threshold a suitably
large power of \(X\), meeting both the coefficient and
point height requirements. Its bound is \(O_{K,\varepsilon}(X^\varepsilon)\).
This proves (EC1).
[Dujella, Corollary 1.3](https://arxiv.org/html/2312.03655v1#S1)

## Uniform product fibres and moments

Let \(\mathcal F\) be any family of \(T\) distinct positive
primitive rows \(A+B=C\le X^K\), with no radical assumptions.
Set
\[
                        f(A,B,C)=(A/C,B/C)\in\mathbb Q_{>0}^2 .
\]
Use coordinatewise multiplication. A fixed two-row product
has at most two ordered representations. Indeed for the
\(A/C\) fractions \(x_1,x_2\), the products \(u,v\) determine
\(x_1x_2=u\) and \(x_1+x_2=1+u-v\). A fixed three-row
product has subpower multiplicity by (EC1). For \(m\ge3\),
an \(m\)-row product has multiplicity \(O(T^{m-3}X^\varepsilon)\)
by fixing all but three rows and dividing out their contributions.
This last assertion also covers any fixed-weight product with
three equal nonzero integer weights. Fix the other rows and
take the unique positive root for that common weight in
each component, when it exists.

In particular the equality count between two ordered triples
of product records is \(O(T^3X^\varepsilon)\), and the
threefold product image has size at least
\(T^3X^{-\varepsilon}\) up to a constant. Equivalently, if
one restricts the two coordinate triple products separately,
the product of the sizes of their two dictionaries is at least
\(T^3X^{-\varepsilon}\). This gives a two-coordinate exact
energy estimate independent of support-overlap assumptions.
It does not bound either coordinate energy separately.

## Full digit records at radical precision

Write \(n=X^L\) and suppose each row in \(\mathcal F\) also
satisfies both pairwise radical caps
\[
       \operatorname{rad}(AC),\operatorname{rad}(BC)
                                               \le n^\varrho.
\]
For a prescribed prime pool \(\mathcal P\) of polynomial
prime product, use the binary slot
\[
            M_j(W)=\prod_{p\notin\mathcal P}p^{b_j(v_p(W))},
\]
where \(b_j(e)\in\{0,1\}\) is the \(j\)-th binary digit of \(e\).
This is the quotient \(x_{2^j}(W)/x_{2^{j+1}}(W)^2\)
of projected powers outside the same pool. Choose
\(J\) so that \(2^{J+1}>K\log_2 X\) and \(J=O_K(\log\log X)\).
For a fixed \(m\ge2\) take a modulus \(D>4n^{m\varrho}\)
whose primes belong to \(\mathcal P\). The \(2(J+1)\)-component
record of a row is the pair of residues
\[
                      (M_j(A)/M_j(C),\,M_j(B)/M_j(C))
                               \pmod D \qquad (0\le j\le J).
\]
It is permitted to take a suitable power of two as \(D\) and
\(\mathcal P=\{2\}\). All ratios in the record are units.

At each slot, each component of a prescribed \(m\)-row
product is the residue of a positive rational whose reduced
numerator times denominator is at most \(n^{m\varrho}\).
Record its numerator and denominator dyadic intervals.
There are \(O(\log^2 X)\) interval choices per component,
so the cost over all slots is \(X^{o_X(1)}\). Within a fixed
nonempty pair of intervals the product of the upper endpoints
is at most \(4n^{m\varrho}<D\). Two fractions in those
intervals with the same unit residue have determinant zero
and coincide. The reconstructed exact rational components
at every slot, raised to powers \(2^j\) and multiplied,
give \(\prod_i A_i/C_i,\prod_i B_i/C_i\) outside the fixed
pool. The missing factors on that pool have polynomial
height and cost only \(X^{o_X(1)}\) by fixed-support
enumeration. Thus a fixed modular product gives only
\(X^{o_X(1)}\) choices of exact two-coordinate products.
[Binary digit reconstruction at radical precision](ABC_MODULAR_DIGIT_SUPPORT_AND_ATOMIC_COMPLETION.md#binary-differences-at-the-radical-scale)
· [Uniform fixed-support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

Applying (EC1) after reconstructing the products shows
modular threefold fibres of size \(O(X^\varepsilon)\)
at this precision. More generally \(m\)-fold fibres for
\(m\ge3\) have size \(O(T^{m-3}X^\varepsilon)\).
The exact addition in this argument requires only
\(A/C+B/C=1\). It does not require a shape, an anchor
row or an affine coefficient dictionary. This gives
stronger fibre bounds at radical precision than the
higher positive-product estimates based on a bounded
transcript at an arbitrarily small modulus. Those
small-modulus estimates use a shorter record and
remain applicable there.
[Higher positive products at small modulus (PP3)](ABC_POSITIVE_PROJECTION_TRANSCRIPT_PRODUCTS.md#fixed-positive-products-of-more-rows)
