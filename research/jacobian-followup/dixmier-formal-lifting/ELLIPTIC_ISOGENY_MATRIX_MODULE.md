# A finite matrix differential model of the elliptic three-isogeny

14 November 2026. **QIM1.** On every smooth generic elliptic energy
fiber identified in [QEI1](ELLIPTIC_ISOGENY_AND_QUANTUM_ENERGY.md),
the actual three-isogeny admits a free rank-three differential
module with an explicit basis and full algebraic right action.
Its entire source curve differential-operator algebra embeds in
three-by-three matrices over the target curve operator algebra.
There is also an explicit faithful six-by-six polynomial matrix
differential realization, with spatial generator degree at most
four and no denominators in the position coordinate.

This closes the frame and right-action questions for this fixed
curve model, rather than inferring them from formal freeness.
It does not construct the four Weyl images of Long. The curve
parameter is central here; Long's quantum energy is noncentral.
Moreover the actual source coordinate algebra cannot remain
finite over the target across the zero-energy divisor. Both
boundaries are proved below.

The basic isogeny is standard: see
[Pizzo--Pomerance--Voight, equation (2.7)](https://math.dartmouth.edu/~carlp/3isog-070319.pdf).
QEI1 supplies the identification with the retained Long map.
The scalar obstruction uses
[Bell--Ingalls--Ramkumar, Theorem 1.1](https://arxiv.org/html/1411.3574):
an embedding between differential-operator division fields on
curves cannot lower the genus. No novelty priority is asserted.

## A global primitive on the smooth affine isogeny

Let k be a characteristic-zero field, fix b in k nonzero, and put
\[
 A=k[u,v]/(v^2-u^3+27b),\qquad
 B=k[x,x^{-1},y]/(y^2-x^3-b).
                                                               \tag{1}
\]
The negative standard three-isogeny induces A -> B by
\[
             u=(x^3+4b)/x^2,\qquad
             v=-y(x^3-8b)/x^3.                                \tag{2}
\]
The source is the elliptic curve with the three kernel points
removed; the target is the elliptic curve with its origin
removed. In characteristic zero the restricted isogeny is
finite etale of degree three.

A particularly useful primitive is w=y/x. Direct substitution
gives
\[
 \begin{split}
       &4w^3-3uw+v=0,\\
       &x=(4w^2-u)/3,\qquad y=(2uw-v)/3=xw,\\
       &x^3-ux^2+4b=0,\qquad
                         x^{-1}=(ux-x^2)/(4b).
 \end{split}                                                    \tag{3}
\]
Consequently
\[
 B=A[w]/(w^3-3uw/4+v/4)
               =A\oplus Aw\oplus Aw^2.                         \tag{4}
\]
For a direct algebraic proof, the formulas (2) define the map
from the displayed quotient to B. Conversely (3) defines
x,y,x^{-1} in that quotient and verifies y^2=x^3+b and (2).
The maps are inverse. The cubic is monic, so the three stated
elements are a free basis in every degree.

Its discriminant is
\[
 \operatorname{disc}_w(w^3-3uw/4+v/4)
    ={27\over16}(u^3-v^2)={729\over16}b,                       \tag{5}
\]
a unit. Equivalently its derivative is (9/4)x, a unit by (3).
This directly proves etaleness; no search over bases, poles,
operator orders, or coefficient supports is involved.

All formulas hold over the parameter ring Q[b,b^{-1}], not
merely for numerical b. The restriction b!=0 is essential.

## The full curve differential action is finite and free

On A let
\[
          \delta(u)=v,\qquad \delta(v)=3u^2/2.
\]
It is a nowhere-zero regular vector field. On B it lifts to
\[
        \widetilde\delta(w)=u/2-w^2,\qquad
        \widetilde\delta(x)=-y,\qquad
        \widetilde\delta(y)=-3x^2/2.                            \tag{6}
\]
Thus it is the negative of the usual invariant source vector
field. These identities preserve both defining equations and
the cubic (3), and agree with delta on A.

The invariant vector fields generate the respective tangent
modules: their two coordinate components do not vanish
simultaneously on a smooth fiber. For a smooth affine curve
in characteristic zero, functions and derivations generate
its differential-operator ring. In this case the rings are
explicit Ore extensions
\[
 {\cal D}_h(A)=A[q;h\delta],\qquad
 {\cal D}_h(B)=B[z;h\widetilde\delta],\qquad q\longmapsto z.
                                                               \tag{7}
\]
At h=1 these are the ordinary curve differential-operator
rings; the polynomial h family records the order filtration.
Localization of differential operators is compatible with
the punctures; see
[Stacks, Lemma 10.133.12](https://stacks.math.columbia.edu/tag/0G36).

The same basis (1,w,w^2) makes D_h(B) free of rank three as
a right D_h(A)-module. Write each Ore polynomial with its
coefficient on the left, expand that coefficient in the A
basis, and collect powers of z=q. Uniqueness follows by
comparing the top Ore degree and the free A coefficients.
The left version follows by reversing the Ore normal order.

This is an actual algebraic module with the complete right
operator action. It does not rely on a formal module
identification, stable freeness, or a bounded-charge assumption.

## Explicit three-by-three matrices

Use columns for the basis (1,w,w^2). Left multiplication by w
and z is represented by
\[
 W=\begin{pmatrix}
 0&0&-v/4\\1&0&3u/4\\0&1&0
 \end{pmatrix},\qquad
 Z=qI_3+hC,\quad
 C=\begin{pmatrix}
 0&u/2&v/2\\0&0&-u/2\\0&-1&0
 \end{pmatrix}.                                                \tag{8}
\]
For example delta(w^2)=-uw/2+v/2, explaining the last column
of C. The actual commutator is
\[
 [Z,W]=h\bigl(\delta(W)+CW-WC\bigr)
                           =h(uI_3/2-W^2).
\]
Define
\[
 X=(4W^2-uI_3)/3,\quad Y=(2uW-vI_3)/3,\quad
                         X^{-1}=(uX-X^2)/(4b).                \tag{9}
\]
Then
\[
 \begin{gathered}
  Y^2=X^3+bI_3,\qquad XX^{-1}=I_3,\\
  [Z,X]=-hY,\qquad [Z,Y]=-3hX^2/2,\\
  uX^2=X^3+4bI_3,\qquad
                         vX^3=-Y(X^3-8bI_3).
 \end{gathered}                                                \tag{10}
\]
These relations include the full source curve algebra and
its invariant derivation, rather than only a commuting
spectral pair.

Left multiplication on the free right module is faithful:
an operator annihilating every element annihilates 1.
Thus (8)--(10) define an injective homomorphism
\[
                       {\cal D}_h(B)
                              \longrightarrow
                         \operatorname{Mat}_3({\cal D}_h(A)).
                                                               \tag{11}
\]
At nonzero h, the quotient division extension has left and
right degree three by the coefficient-field argument in
[QDI1](QUANTUM_DIVISION_AND_COMPLETION.md). This finiteness
concerns elliptic curve operator rings, not a Weyl algebra
finite over the image of its own endomorphism.

## Polynomial six-by-six Weyl matrices

The target coordinate ring A is free over k[u] with basis
(1,v). Let p be the ordinary Weyl momentum, [p,u]=h, and set
\[
 V_2=\begin{pmatrix}0&u^3-27b\\1&0\end{pmatrix},\qquad
 Q_2=V_2p+h\begin{pmatrix}0&3u^2/2\\0&0\end{pmatrix}.            \tag{12}
\]
Multiplication by v is V_2, and h delta is Q_2. They satisfy
\[
 V_2^2=(u^3-27b)I_2,\quad
 [Q_2,uI_2]=hV_2,\quad [Q_2,V_2]=3hu^2I_2/2.
\]
These are polynomial matrix differential operators, even at
the branch points v=0. Inverting v to use the generator
partial_u instead would create artificial position poles;
the invariant derivation avoids them.

This representation is faithful in every operator degree.
For a nonzero Ore polynomial in q with top coefficient a
and degree j, the top p coefficient is M_a V_2^j. Over k(u),
multiplication by any nonzero function a on the integral
curve is invertible, and V_2 is invertible as well. The top
coefficient cannot vanish.

Apply this two-by-two representation to every entry in (8)
and (9). It gives
\[
                {\cal D}_h(B)\longrightarrow
                       \operatorname{Mat}_6(A_{1,h}(k)).
                                                               \tag{13}
\]
For generators (w,x,y,x^{-1},z), the respective spatial degrees
in u,p are (3,3,4,4,4). The bound four ignores the scalar
parameters b and h. Every entry is polynomial in u,p,h with
coefficients in Q[b,b^{-1}]; no spatial completion is used.

Equations (11) and (13) are finite representations of a curve
differential algebra. They are not endomorphisms of a Weyl
algebra or of a full matrix Weyl algebra. Their existence
therefore proves no Dixmier counterexample.

## An actual Weyl pair and an obstruction in every matrix size

The source curve algebra contains the explicit Weyl pair
\[
                Q=-{2\over3}x^{-2}z,\qquad [Q,y]=h.
\]
Indeed [z,y]=-3hx^2/2, and x commutes with y. Its images in
(13) are finite polynomial six-by-six Weyl matrices. This
has spatial degrees (5,4) for (Q,y) and supplies an actual
full rank-one Weyl action inside the
curve model, not just formal differential symbols. It is
still not an endomorphism of the full matrix Weyl algebra:
no compatible matrix units making such an endomorphism are
supplied.

The same identity gives a precise obstruction to a naive
central-energy promotion, for every matrix size. Consider
the total relative curve Ore algebra
R_0=B_0[z;h delta], with B_0=k[b,x,x^{-1},y]/(y^2-x^3-b),
before inverting b, where delta(x)=-y,
delta(y)=-3x^2/2 and delta(b)=0. Let P be any positive-order
finite scalar ordinary differential operator over a
characteristic-zero differential function field, with h a
nonzero differential constant. There is no unital
homomorphism
\[
               R_0\longrightarrow\operatorname{Mat}_m(D)
\]
sending b to P times the identity and h to its scalar value.

Every image entry would commute with P because b is central
in R_0. The scalar finite-operator centralizer of a
positive-order operator is commutative by the Schur theorem;
see [Previato--Rueda--Zurro, Corollary 2.4](https://arxiv.org/html/1902.01361).
All matrices would therefore have entries in one commutative
ring. Taking the matrix trace of [Q,y]=hI_m gives 0=mh,
a contradiction. No faithfulness assumption or order bound
on the matrix entries is required.

The trace argument is general: any algebra containing a Weyl pair
and a central parameter has no finite-matrix differential realization
that sends that parameter to a scalar finite operator outside the
differential constants. For an order-zero coefficient f with
delta(f)!=0, its finite-operator centralizer is the coefficient field:
the top commutator coefficient of an order-j operator is
j*a_j*delta(f), so a commuting operator must have j=0.
This covers the remaining nonscalar order-zero case as well.

This excludes the direct central-parameter substitution in
any finite matrix rank. It does not exclude a construction
that changes the relative algebra, the image of the energy,
or the Hamiltonian reduction itself.

## Why the scalar version of this model is excluded

Over algebraically closed characteristic-zero constants
and nonzero h, any injective scalar rational differential
realization of D_h(B) would extend to an embedding of its
Ore quotient division ring into D_1(k), the operator division
ring of a rational curve. B is an open subset of a smooth
projective genus-one curve, so its quotient operator field
has genus-one coefficient curve.

Bell--Ingalls--Ramkumar's cited theorem would require
1<=0, a contradiction. Thus the full elliptic curve operator
algebra cannot be carried faithfully by scalar rational
differential operators, regardless of their finite orders
or coefficient poles. Matrix representations evade this
specific scalar obstruction, as (12)--(13) exhibit.

This does not exclude a scalar Weyl lift of Long: that image
is a Weyl algebra, not this elliptic curve differential
algebra. No embedding of the latter inside a hypothetical
Long image has been established. Nor does this argument
exclude a quotient or Hamiltonian reduction of the matrix
model; it excludes faithful scalarization of the whole
curve operator algebra.

## The actual source algebra is not finite across zero energy

It is tempting to use the free basis (4) uniformly across b=0.
Here is the intrinsic obstruction to doing so for the actual
source ring, rather than just for one chosen basis.

Consider the total families
\[
 A_0=k[b,u,v]/(v^2-u^3+27b)\simeq k[u,v],\qquad
 B_0=k[b,x,x^{-1},y]/(y^2-x^3-b)\simeq k[x,x^{-1},y].
\]
The formulas (2) still define A_0 -> B_0. Its fraction-field
degree is three. In the finite monic algebra A_0[w], the
multiplication determinant of x=(4w^2-u)/3 is
\[
                            \operatorname{Norm}(x)=-4b.       \tag{14}
\]
Hence Norm(x^{-1})=-1/(4b). Since A_0 is a normal polynomial
ring, x^{-1} cannot be integral over A_0: its norm would
otherwise be integral and lie in A_0. Consequently B_0 is
not finite over A_0.

Equivalently,
\[
 B_0=A_0[w,1/x],\qquad A_0[w]\ne B_0.
\]
The free monic algebra does extend over b=0, but it loses the
required unit x. No change of module basis makes the actual
B_0 finite. At b=0, the equation for x is
\[
                             x^2(x-u)=0;
\]
the two x=0 branches are precisely omitted from the source.
On the remaining generic branch the map is u=x, v=-y, of
degree one. Thus the degree-three isogeny family itself has
a three-to-one versus one-to-one boundary degeneration.
This uses the energy parameter b, not the quantum parameter
h, and is not a counterexample to a Long-specific h-index
comparison.

For Long, b=16r^2 lambda on the classical elliptic fiber.
The construction is therefore valid on the generic smooth
energy locus, not across lambda=0 or r=0. Global admission
cannot be inferred from its fiberwise freeness.

## The noncentral-energy bridge remains open

There is a second, logically separate issue. In (1)--(13),
b is a central scalar and x,y are commuting curve functions.
The operator z is an independent cotangent direction.
Thus D(B) quantizes the cotangent bundle of a curve; it
does not quantize the original two-dimensional family of
energy levels with its symplectic bracket.

In a canonical Weyl pair U,V with [U,V]=c nonzero central,
the cubic energy P=V^2+4U^3 satisfies
\[
                       [P,U]=-2cV,\qquad [P,V]=12cU^2.
\]
It is noncentral. Substituting it for b in a central-parameter
matrix construction does not preserve the defining
commutation relations. Fixed-energy geometry is not a
quantum coefficient curve without an additional theorem.

A genuine bridge would have to construct a compatible
deformation or Hamiltonian reduction with the correct
noncentral energy, retain all four Long images, and pass
global polynomial admission at the exceptional divisors.
It cannot simply embed the whole elliptic operator algebra
into scalar operators, by the genus obstruction.

The present result supplies the explicit finite frame and
right action on smooth elliptic fibers. It also locates two
remaining obligations: noncentral-energy compatibility and
global extension over the singular-energy locus. It proves
neither a common bound for a full Long lift nor an intrinsic
obstruction to every such lift.

## Exact controls and proof boundary

[The checker](check_isogeny_module.py) verifies the monic
cubic, constant discriminant, inverse coordinate formulas,
the full three-by-three differential action, its six-by-six
polynomial realization and degree bounds, and the norm and
boundary identities. It binds QEI1's retained Long-isogeny
identification by hashes. There is no numerical sampling
of energies or increasing-order search.

Freeness and faithfulness in every operator degree,
nonfiniteness of the actual total source algebra, and the
scalar exclusion use written arguments and the cited genus
theorem. No independent theorem replay, formal verification,
external review, novelty priority, bounded full Long lift,
all-representative nonlifting theorem or DC2 solution is
claimed.
