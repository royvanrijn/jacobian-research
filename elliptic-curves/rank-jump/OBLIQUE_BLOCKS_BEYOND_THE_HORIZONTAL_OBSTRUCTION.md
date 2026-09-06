# A split residual cubic gives two directions beyond the horizontal obstruction

A common quadratic cover can produce two independent rational directions
even when **every rational horizontal pair is excluded**. The following
small construction closes that broader mechanism:

\[
 \boxed{E_t:\ y^2=(x+1)(x-2)(x-\tfrac12)+t x^2,
 \qquad t=n^2.}
\]

The arithmetic generic rank is exactly **1 over Q(t)** and exactly
**3 over Q(n)**. The fixed n=5 fibre has three certified independent
points and an odd inert valuation v₁₇(−A)=1 in its short model.
Thus it fails the [horizontal norm test](HORIZONTAL_BLOCK_OBSTRUCTION_ON_A1.md)
while supporting this independent two-direction block.

This is a constructive mechanism control, not a new high-rank search or
a claim that a production family has this representation.

## The single simultaneous-solubility condition

Set

\[
 g(x)=(x+1)(x-2)(x-\tfrac12)
      =x^3-\tfrac32x^2-\tfrac32x+1.
\]

The point R=(0,1) is rational on E_t before the base change. For each
root a of g,

\[
 F_t(a)=t a^2.
\]

The three values are t,4t,t/4: they have one nontrivial squareclass [t]
over Q(t), despite being unequal. The rational double cover n²=t gives

\[
 P=(-1,-n),\qquad Q=(2,2n),\qquad T=(\tfrac12,\tfrac n2).
\]

All three points lie on the **oblique** line y=nx, so

\[
 P+Q+T=O.
\]

There are at most two new directions from this triple. The incidence
and independence proofs below show that this bound is attained.

The condition t=n² supplies actual rational points simultaneously;
it is a solubility condition, not an estimate from Selmer parity or
point-search visibility. The cover n²=t is rational, so its parameter
is explicit. No class-group or descent calculation is needed to make
these three sections rational.

## What splits, and what does not

Over Q(n), use the shear Y=y−nx. The equation becomes

\[
 \boxed{Y^2+2nxY=g(x).}
\]

At Y=0 the residual cubic splits into its three fixed rational roots.
The square condition is exactly what makes this shear and these points
rational over the base. The plane cubic defining E remains smooth and
irreducible at the generic point.

This split residual cubic is **not the two-division cubic**. The
completed ordinate is 2Y+2nx, which is nonzero at those three roots
when n≠0. In particular, they are not newly rational two-torsion points.

Both the parent and base-changed two-division cubics have Galois group
S₃. The discriminants below are nonsquare; irreducibility follows from
the monic anchor specialization modulo 7. Thus this two-direction rank
increase does not require a drop of the two-division Galois group.

This also illustrates why raw-y patterns can mislead. The three points
have raw Y=0 after the shear, but their completed ordinates are −2n,4n,n.
Their absolute values are distinct for n≠0, consistently with the
horizontal obstruction.

## Exact geometric and arithmetic generic ranks

The original Weierstrass equation is

\[
 y^2=x^3+(t-\tfrac32)x^2-\tfrac32x+1.
\]

Its invariants are

\[
 c_4=16t^2-48t+108,
 \qquad\Delta=-64t^3+324t^2-972t+729.
\]

The discriminant is squarefree and coprime to c₄. There are three
finite I₁ fibres. At infinity, the regular chart x=v⁻²X,y=v⁻³Y
has orders (v(c₄),v(Δ))=(2,9), hence an I₃* fibre. The surface has
Euler number 12 and is rational. Its reducible-fibre root lattice is
D₇, so Shioda–Tate gives

\[
 \operatorname{rank}E_t(\overline{\mathbf Q}(t))=10-2-7=1.
\]

After t=n², the discriminant has six simple finite roots, all type I₁;
Δ(0)=729 ensures that ramification at n=0 is over a smooth fibre.
At infinity the invariant orders are (0,6), giving I₆. The new surface
is again rational, with root lattice A₅. Thus

\[
 \operatorname{rank}E_{n^2}(\overline{\mathbf Q}(n))=10-2-5=3.
\]

These are applications of the standard rational elliptic surface and
Shioda–Tate formulas in
[Schütt–Shioda, §§6 and 8](https://arxiv.org/pdf/0907.0298).

The three rational sections R,P,Q specialize to three independent
classes at the fixed anchor below. Hence they are generically
independent and attain the geometric upper bound. R alone attains the
parent upper bound. This proves

\[
 \boxed{\operatorname{rank}E_t(\mathbf Q(t))=1,
 \qquad\operatorname{rank}E_{n^2}(\mathbf Q(n))=3.}
\]

The deck involution n↦−n fixes R and negates P,Q. Accordingly the
two new directions lie together in the anti-invariant character. There
is no remaining geometric direction whose field of definition is unknown.

## An independent block on a horizontally obstructed fibre

At n=5, t=25, the equation is

\[
 y^2=x^3+\tfrac{47}{2}x^2-\tfrac32x+1.
\]

The integral transport x_int=4x,y_int=8y gives

\[
 y_{\rm int}^2=x_{\rm int}^3+94x_{\rm int}^2-24x_{\rm int}+64
\]

and the prescribed points

\[
 R=(0,8),\qquad P=(-4,-40),\qquad Q=(8,80).
\]

Its cubic is irreducible modulo 7, so there is no rational two-torsion.
Good-prime split-root Kummer fingerprints, together with the real
signature, have rank three. For example the finite columns at 5 are
(6,5,5), and at 29 are (0,5,6); their concatenation already has rank
three. This is exact subgroup independence, not a numerical height
estimate. The whole curve rank is not determined here.

Completing the original equation to short form gives

\[
 -A=\frac{2227}{12}=\frac{17\cdot131}{12}.
\]

Because 17≡2 mod 3, its odd valuation prevents every distinct rational
horizontal pair on this same curve. The obstruction also rules out
horizontal rational-function pairs on either parameter line by the
previous specialization lemma. The oblique block therefore cannot be
turned into a horizontal one by a rational Weierstrass change.

At n=0, by contrast, the three new displayed points become two-torsion
and contribute no new free direction. Generic independence is not a
blanket specialization theorem. Once the double cover is adopted as
the base, R,P,Q form its generic rank-three subgroup and must not be
counted again as an additional jump in that presentation.

## A guard against dependent shared-cover identities

Shared squareclasses can also produce arbitrarily many **dependent**
representations. On any short curve F(X)=X³+AX+B, duplication gives

\[
 R_2(X)=\frac{X^4-2AX^2-8BX+A^2}{4F(X)},
\]
\[
 H_2(X)=\frac{X^6+5AX^4+20BX^3-5A^2X^2-4ABX-8B^2-A^3}{8F(X)^2},
\]

with the universal identity F(R₂(X))=F(X)H₂(X)². It uses precisely
the same quadratic cover as the original point, but represents 2P,
so it adds no Mordell–Weil direction. The test suite verifies this
identity coefficientwise and derives H₂ from the doubling slope.

More generally, suppose E/K has no rational two-torsion and
End_K(E)=Z. A nonconstant universal identity

\[
 F(R(X))=F(X)H(X)^2,\qquad R,H\in K(X),
\]

defines the odd map (X,Y)↦(R(X),H(X)Y) from E to itself. It extends
over the smooth projective curve. Oddness forces the image of O to
be rational two-torsion, hence O. An origin-preserving map is an
endomorphism, so it is [k], and deg R=k². See
[Milne, *Elliptic Curves*, Proposition 1.5](https://www.jmilne.org/math/Books/EC2.pdf)
for the group-homomorphism statement. Such a universal one-variable
identity cannot supply an independent second direction under these
hypotheses.

Our construction is instead an identity at particular abscissae over
the parameter field. For example F_t(2)=4F_t(−1); the candidate
universal substitution R(X)=−2X,H=2 fails away from the prescribed
roots, since

\[
 F_t(-2X)-4F_t(X)=-3(X+1)(2X-1)^2.
\]

The separate independence certificate is essential. Merely finding
two formulas using the same radical would not distinguish this block
from the multiplication examples.

## Experiment and current mechanism ranking

The [frozen protocol](OBLIQUE_SPLIT_CUBIC_PROTOCOL.json) permits one
30-second computation, one fixed anchor n=5 and fingerprint primes
at most 199. No parameter sweep, point search, class group or descent
ran. The [certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_oblique_split_cubic_v1.json)
retains the symbolic invariants, fibre data, point transports and all
fingerprints. Five tests independently check the fingerprints, the
oblique-line identities, the local norm obstruction, the universal
duplication identity and the branch-point torsion collapse.

```sh
sage -python elliptic-curves/rank-jump/oblique_split_cubic.py check
sage -python -m unittest discover -s elliptic-curves/rank-jump -p test_oblique_split_cubic.py
```

1. **Constructively supported:** a split residual cubic after subtracting
   d(t)h_t(x)² can give several rational multisections on one cover.
   A square specialization of d can make an independent block rational.
   This example proves two directions while preserving S₃ and failing
   the horizontal norm gate.
2. **Insufficient alone:** residual roots give points but have a line
   relation; shared radicals can also come from multiplication maps.
   Independence modulo the correct generic subgroup remains mandatory.
3. **Still missing for production:** a representation of an existing
   A1/MW16 or R17 family with a nontrivial d(t), constructed without
   exceptional points, and enough independent multisections to explain
   a large quotient. This control supplies no such representation.
4. **Potential use for Agent 1:** a proved shared-cover condition could
   be a simultaneous-solubility selector. A residual cubic fitted to
   known points at a single fibre would be retrospective and cannot
   enter that selector. No active search policy was changed.
