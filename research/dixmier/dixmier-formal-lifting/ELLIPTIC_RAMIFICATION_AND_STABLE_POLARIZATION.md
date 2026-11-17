# An elliptic obstruction to every stable cotangent polarization of Long

14 November 2026. **QHR1.** Let Phi be Long's full polynomial Poisson
endomorphism, with images (R,T,D,S). For H in its image, consider
local nilpotence of the Hamiltonian derivation {H,-} on the entire
ambient polynomial algebra. Then
\[
 \{H\in\operatorname{im}\Phi:\{H,-\}\text{ is locally nilpotent}\}
                                  =K[R].                       \tag{1}
\]
No assumption that H commutes with R is needed. This strengthens
[QPA1's fiber-image statement](RAMIFICATION_AND_POSITION_AVOIDANCE.md).

The result is stable. After adjoining k identity canonical pairs,
every such H belongs to K[R, spectator variables]. The maximum
number of algebraically independent, pairwise Poisson-commuting,
ambient locally nilpotent Hamiltonians in the stabilized image is
exactly k+1, while the symplectic rank is k+2. Consequently no
polynomial symplectic source and target automorphisms, in any degrees,
can put Long or any of these stabilizations into a form with all
position images positional. In particular they cannot be cotangent
lifts in any such polynomial canonical coordinates.

This is unconditional and does not assume prescribed arithmetic
support or a quantization of the canonical automorphisms. It also
excludes a positional form after regular formal quantum source and
target automorphisms with regular inverses, by classical reduction.
It does not exclude arbitrary high-order quantum lifts, parameter-pole
changes without invertible classical reduction, or other target
changes that do not seek a full position-preserving form.

The central calculation is an exact elliptic differential identity.
It proves an all-degree Poisson centralizer theorem and turns a
nonproperness component into an invariant under canonical coordinates.

## The discriminant as a nonproperness component

Use source output coordinates (r,t,d,s), with
{d,r}={s,t}=1, and let Phi send them to (R,T,D,S).
All geometric statements can be proved after algebraic closure of K;
the resulting polynomial conclusions descend to every characteristic-
zero field.

[QPA1](RAMIFICATION_AND_POSITION_AVOIDANCE.md) gives, for r nonzero,
the actual fiber map
\[
 s=2v+4a-3ra^2,\qquad
 t=(ra^3-a^2-va)/2,\qquad
 f(X)=rX^3-2X^2+sX+4t,\quad f'(a)=2v.                           \tag{2}
\]
The target discriminant is
\[
 \Delta=4s^2-4rs^3+128t-144rst-432r^2t^2.                        \tag{3}
\]
Here r=R, rather than R/2.

The d coordinate is free in this chart. To check this without assuming
a symplectic source chart, use the ordinary polynomial coordinates
(x,Y,beta,D) in
[Long, Proposition 4.1 and Section 6](https://arxiv.org/html/2608.23777v1).
They have the explicit inverse
q=Y+x*beta/3,
p=3q^2*(beta+9q^2)+2(1-3xq)D0,
z=(1+3xq)*(beta+9q^2)/2-3x^2*D0,
where D0=D-H(x,Y,beta) and H is the retained polynomial correction.
The determinant of the displayed momentum change is -1.
Since H is independent of D0, shifting D0 to D preserves invertibility.

On x nonzero the further ordinary coordinate change is
\[
 x=v^{-1},\qquad Y=a-v,\qquad
            \beta=5v^2-3va-rv^3,\qquad D=d.                     \tag{4}
\]
Its inverse is v=1/x, a=v+Y and
r=2x-3x^2Y-x^3*beta. Thus over output r nonzero every preimage has
x nonzero and the whole map is (2), with r,d unchanged and v invertible.
These ordinary coordinates are used only for properness; no canonical
property is attributed to (4).

Complete the chart by allowing v=0. The map in (2) then is finite of
degree three, because a satisfies the monic cubic f/r and v is its
displayed polynomial derivative. The removed divisor v=0 maps onto
Delta=0. Its parameterization is
s=4a-3ra^2, t=(ra^3-a^2)/2.
Off Delta the original chart is the whole finite cover, hence proper.
On Delta, fix r,a,d and let v approach zero: the target extends while
x=1/v has a pole. The valuative criterion therefore gives nonproperness.
Consequently
\[
       N(\Phi)\cap\{r\ne0\}=\{\Delta=0,\ r\ne0\},                \tag{5}
\]
where N is the set of target points over which the point map is not
proper on any neighborhood.

The set N is closed: compactify the source in projective space,
take the closure of its graph, and project its boundary to the target.
It is proper because a generically finite polynomial map is finite
over some nonempty open set, by clearing the denominators of integral
equations for its finitely many source coordinate generators.
Therefore the closure of the divisor in (5) is an irreducible
component of N.

For irreducibility and later use, set
\[
       X=3rs-4,\qquad Y=108r^2t+18rs-16.
\]
Over K(r) this is an invertible affine coordinate change and
\[
                    Y^2+4X^3=-27r^2\Delta.                    \tag{6}
\]
The cusp polynomial Y^2+4X^3 is irreducible. The coefficient 128
of the monomial t makes Delta primitive over K[r], so it is
irreducible in K[r,s,t] as well. Thus Delta=0
is the required component globally. No assertion about other possible
components of N is needed.

Adjoining identity spectator pairs takes this component to its product
with the spectator affine space. The same argument applies verbatim.

## An elliptic differential with no algebraic primitive

Put
\[
 V(g)=g_s\Delta_t-g_t\Delta_s,\qquad
 A=\frac{4rst-32t-s^2}{20r^2}.
\]
A direct identity in K(r)[s,t] is
\[
                      \Delta_r=V(A)+\frac{8\Delta}{5r}.         \tag{7}
\]
It is essential that the final term is nonzero. Work on the generic
level curve Delta=lambda over K(r,lambda), allowing any extra
spectator variables in the constant field. Equation (6) becomes
\[
                       Y^2=-4X^3-27r^2\lambda.                 \tag{8}
\]
Its smooth projective completion has genus one: projection to X
is a double cover branched at the three distinct roots and infinity.
The differential
\[
               \omega=\frac{ds}{\Delta_t}
                              =-\frac{dX}{24rY}                \tag{9}
\]
is nonzero and holomorphic on that completion. At a finite branch
point the local coordinate Y removes the apparent pole; at infinity,
X has pole order two and Y pole order three, so dX/Y is regular
and nonzero.

Since omega(V)=1, equation (7) is the identity
\[
                     \Delta_r\omega
                       =dA+\frac{8\lambda}{5r}\omega.           \tag{10}
\]
It proves that Delta_r has no rational primitive for V on the
generic curve. Otherwise subtracting A would give a rational
primitive of a nonzero holomorphic differential. Any pole of a
rational function creates a pole of its differential in characteristic
zero; with none, the function is constant on the projective curve.

The same obstruction survives every finite algebraic cover of the
generic curve: pullback of omega remains nonzero and holomorphic.
Thus no algebraic primitive exists either. This is a direct
holomorphic-differential argument, without estimating candidate
primitive degrees or extrapolating from samples.

These elliptic curves have j=0. Their family is isotrivial after
an algebraic scaling; nonisotriviality is not used. The nonzero
differential class in (10), rather than variation of j, is the
actual obstruction.

## The full Poisson centralizer of Delta

Let w denote any number of independent spectator coordinates.
Write k=K(r,w), temporarily treating them as constants for V.
The full Hamiltonian commutator with Delta on k[s,t,d] is
\[
                     \{F,\Delta\}=V(F)+\Delta_r F_d.           \tag{11}
\]
Although the spectators have their own Poisson brackets, Delta is
independent of them, so they contribute no terms to (11).

The constants of V on k(s,t) are exactly k(Delta). Indeed on the
generic geometrically integral curve (8) a rational function with
zero differential is algebraic over k(Delta) in characteristic zero.
Geometric integrality makes that field relatively algebraically
closed in the curve's function field: an additional finite constant
extension would split it after algebraic closure. Intersecting with
k[s,t] gives k[Delta]. For this last assertion, a reduced rational
function in Delta has coprime numerator and denominator by the
univariate Bezout identity; a nonconstant polynomial in Delta is
not a unit in k[s,t].

Now let F=sum_{j=0}^m a_j(s,t)d^j commute with Delta.
Its highest d coefficient has V(a_m)=0, so a_m belongs to
k[Delta]. If m>0, the next coefficient in (11) gives
\[
                       V(a_{m-1})=-m a_m\Delta_r.              \tag{12}
\]
On the generic level, divide by the nonzero scalar m*a_m(lambda).
This is the forbidden rational primitive from (10).
Therefore m=0, and we have the all-degree centralizer
\[
           \operatorname{Cent}_{\rm Pois}
                   (\Delta,\ k[s,t,d])=k[\Delta].              \tag{13}
\]
No restriction on d degree, s,t degree or spectator number enters
the proof.

## Locally nilpotent Hamiltonians preserving the component

Suppose F is a polynomial Hamiltonian in the full source algebra,
{F,Delta}=0, and {F,-} is locally nilpotent.
Equation (13) makes it independent of d and, rationally, a function
of r,Delta,w. In particular it fixes r and Delta.

Every additive-group orbit with generic fixed values r and Delta
projects, via s,t, to the affine elliptic curve (8). That projection
must be constant: a nonconstant map from the affine line to the
projective completion would extend to a map from the projective line
to a genus-one curve, contradicting
[Riemann--Hurwitz](https://stacks.math.columbia.edu/tag/0C1B).
The coordinates of the orbit are polynomials in its group parameter,
so the projected map is defined everywhere on the affine line.

Consequently the Hamiltonian kills s and t. In the expression from
(13) this forces its Delta derivative to vanish; hence
\[
                         F\in K[r,w].                         \tag{14}
\]
To justify using generic values, if either derivative were nonzero,
it would be nonzero on an open set meeting r*Delta nonzero. An orbit
from such a point would contradict the preceding argument.
Because the original F is polynomial, the rational independence of
s,t,d in (14) leaves an actual polynomial in r,w.

For no spectators this says F belongs to K[r], and every such
polynomial does generate a locally nilpotent Hamiltonian, translating
d. With spectators, (14) is a necessary condition; not every
Hamiltonian in those variables is asserted locally nilpotent.

## Transfer from nonproperness to the whole image

Let H=Phi(F) generate a locally nilpotent Hamiltonian on the entire
ambient polynomial algebra, with any identity spectators added.
Poisson compatibility gives
{H,Phi(G)}=Phi({F,G}). Injectivity implies that {F,-} is locally
nilpotent as well. Their additive-group actions make the point map
equivariant. Properness is invariant under algebraic isomorphisms,
so the source-output action preserves N(Phi).

An additive group is connected and therefore fixes each irreducible
component of this closed set. It preserves Delta=0.
For a locally nilpotent derivation delta on a domain, preservation
of a principal prime ideal forces delta(Delta)=0. Here is an
elementary proof: the nilpotence degree of a product is the sum
of those of its factors, by the exponential polynomial in the
group parameter. If delta(Delta)=B*Delta were nonzero, its degree
would be both one less than that of Delta and at least that degree,
a contradiction.

We may therefore apply (14) to F. It follows that
\[
 \{H\in\operatorname{im}(\Phi\times{\rm id}):
          \{H,-\}\text{ locally nilpotent on the ambient ring}\}
                   \ \subseteq\ K[R,w].                       \tag{15}
\]
When there are no spectators, the reverse inclusion holds because
R is positional, proving (1). This conclusion no longer assumes
H commutes with R; it proves that restriction for every such H.

## A stable obstruction under all polynomial canonical coordinates

For a polynomial Poisson image B in symplectic rank n, let nu(B)
be the largest number of algebraically independent, pairwise
Poisson-commuting elements of B whose Hamiltonians are locally
nilpotent on the entire ambient polynomial ring.
A polynomial symplectic target automorphism preserves nu by
conjugating Hamiltonian derivations. A polynomial source
automorphism leaves the image algebra unchanged.

For Long plus k identity pairs, (15) places all these elements
inside K[R,w_1,...,w_{2k}]. Its Poisson form has rank 2k and
one-dimensional radical, generated by dR, over its fraction field.
Differentials of algebraically independent commuting functions
span an isotropic subspace. A skew form of rank 2k with radical
dimension one has no isotropic subspace of dimension greater than
k+1. The functions R and the k spectator positions attain this
bound and have locally nilpotent Hamiltonians. Hence
\[
                   \nu(\operatorname{im}
                        (\Phi\times{\rm id}_{2k}))=k+1.         \tag{16}
\]

A position-preserving polynomial Poisson endomorphism in rank
k+2 has at least k+2 such elements: its independent position
images. A positional Hamiltonian is locally nilpotent because
it lowers momentum degree. Thus (16), and its invariance under
source and target symplectic automorphisms, prove:
\[
 \boxed{\text{Long and every identity stabilization have no
 polynomial canonical position-preserving form.}}              \tag{17}
\]
In particular no cotangent form is possible, even after stabilization,
arbitrary polynomial canonical degrees and arbitrary source/target
canonical word lengths. This concerns a single position-preserving
normal form, not every composition of several noninvertible
cotangent maps interleaved with canonical transformations.
The latter class has a separate conditional arithmetic exclusion
in [QFL3](ARITHMETIC_COTANGENT_EXCLUSION.md).

The distinction from ordinary polynomial coordinates is sharp. After
one identity pair, Long is ordinarily polynomially equivalent to
G times the three-dimensional identity, by its global ordinary
coordinates (x,Y,beta,D).
The cotangent lift of the three-variable Keller core G is also
ordinarily equivalent to that same map: use the polynomial fiber
change pi -> J_G(u)^{-T}pi, whose inverse is J_G(u)^T pi.
Thus these two rank-three maps are ordinarily equivalent, while
(16) gives values two and three for nu and forbids a polynomial
symplectic equivalence. Ordinary equivalence does not transport
this Hamiltonian condition or a Weyl realization.

The conclusion is stronger than position avoidance in one chosen
target chart, but it remains a classical geometric restriction.
It is not a theorem that arbitrary Weyl endomorphisms must admit
such a position-preserving form.

## Consequence for formal quantum coordinate changes

Consider any full formal lift with classical reduction Long plus
identity pairs. Let source and target quantum automorphisms and
their inverses be regular in h, coefficientwise polynomial, with
polynomial symplectic automorphisms as classical reductions.
They may have unbounded spatial degree across h orders.
If their transformed lift had all position images positional,
its reduction would violate (17).

Thus no such regular two-sided formal normalization restores a
position-preserving or cotangent form. This does not require
keeping R exact, or even keeping its quantum image positional.
It applies to the full formal gauge class under the stated
regularity hypothesis. Parameter-pole changes without invertible
classical reduction are not covered.

Nothing here implies the nonexistence of a bounded full Weyl lift.
Such a lift may have genuinely nonpositional, high-order images.
Local nilpotence of source coordinate Hamiltonians on the image
does not imply local nilpotence on the ambient target. That missing
implication is neither assumed nor proved.

## Verification and reusable method

[The checker](check_elliptic_ramification.py) verifies (6)--(7)
exactly, the smooth Weierstrass discriminant, the source-chart
inverse and branch formulas, and the arbitrary-d-power identity
underlying (12). It checks explicit spectator controls, both
source-coordinate Poisson brackets and the retained input hashes.
These are symbolic identities. The nonproperness component,
absence of algebraic primitives, complete centralizer, local
nilpotence classification and stable invariant are written proofs,
not finite searches. See the
[control receipt](elliptic_ramification_controls.json).

The reusable principle is: identify an invariant nonproperness
component, compute its Hamiltonian centralizer using the
cohomology of its generic level curve, and bound the number of
ambient additive Hamiltonian directions in the image. In this
case one explicit holomorphic elliptic class excludes all positive
conjugate degree at once and supplies an obstruction under every
polynomial canonical polarization.

No independent theorem replay, formal verification, external
review, novelty priority, general finite-cotangent-word exclusion,
parameter-pole normalization theorem, bounded full Long lift,
all-representative quantum nonlifting theorem or DC2 solution is claimed.
