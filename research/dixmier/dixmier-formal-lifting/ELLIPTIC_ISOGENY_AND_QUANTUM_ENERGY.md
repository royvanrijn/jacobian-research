# The elliptic three-isogeny and the forced quantum energy jump

14 November 2026. **QEI1.** Long's generic elliptic fiber map is,
after explicit scalings and a sign, the standard degree-three
isogeny between the j=0 curves y^2=x^3+b and y^2=x^3-27b.
This identifies a concrete connection to elliptic-isogeny research.

Its exact quantum counterpart has a rigid cubic energy. For every
finite scalar Weyl pair X,Y with [X,Y]=c nonzero constant, set
P=Y^2+4X^3. If N is the differential order of P, then
\[
             \operatorname{ord}X\le2(N-1),\qquad
             \operatorname{ord}Y\le3(N-1).                     \tag{1}
\]
Hence a spatially bounded full Long lift with exact undeformed R
cannot keep its quantum discriminant at scalar order two. Its
classical discriminant is quadratic, but its quantum image must
have order at least three, regardless of every other image order.

There are also exact source-algebra results. The finite cubic
oscillator has polynomial centralizer, no nonzero polynomial in
it is a commutator with it, and its energy-preserving finite plane
endomorphism is the identity. The full localized source quantum
discriminant has an explicit centralizer and first-order normalizer.
These statements use all-degree filtered arguments. They do not
claim that a centralizer inside a quantum image is the entire
ambient centralizer, or that bounding one formal energy image
automatically bounds all images.

## The generic fiber is exactly a three-isogeny

Use r=R and the retained classical fiber coordinates of
[QPA1](RAMIFICATION_AND_POSITION_AVOIDANCE.md):
\[
 s=2v+4a-3ra^2,\quad t=(ra^3-a^2-va)/2,\quad
 \lambda=4v^2((3ra-2)^2-8rv).
\]
Define
\[
 \begin{split}
 x&=8rv,& y&=8rv(3ra-2),& b&=16r^2\lambda,\\
 x'&=16-12rs,& y'&=432r^2t+72rs-64 .
 \end{split}
\]
On a generic fixed level lambda these satisfy
\[
              y^2=x^3+b,\qquad (y')^2=(x')^3-27b,
\]
and substitution gives
\[
 x'=\frac{x^3+4b}{x^2},\qquad
 y'=-\frac{x^3-8b}{x^3}\,y.                                   \tag{2}
\]
This is the negative of the usual normalized three-isogeny:
see [Pizzo--Pomerance--Voight, equation (2.7)](https://math.dartmouth.edu/~carlp/3isog-070319.pdf).
Our assertion is the explicit identification of Long's retained
map with these formulas; no novelty priority is claimed.

For completeness, the formulas extend to a morphism between the
smooth projective generic curves and send the point at infinity
to infinity. The x' function has poles of order two at infinity
and at the two points x=0. Since the target x coordinate has a
double pole at infinity, the map has degree three. A morphism
of elliptic curves taking origins to origins is an isogeny.
Equivalently its geometric kernel is the infinity point and the
two order-three points (0,+sqrt(b)), (0,-sqrt(b)).
The pullback of the invariant differential is
\[
                             \frac{dx'}{y'}=-\frac{dx}{y}.      \tag{3}
\]
Indeed, at fixed b, dx'/dx=1-8b/x^3, and this factor cancels in
(2). These are generic smooth levels; the singular level
lambda=0 is not being treated as an elliptic curve.

This isogeny identification does not provide a finite operator
quantization. The rest of this note gives a restriction on such
a quantization that applies to arbitrary finite correction orders.

## An exact quantum elliptic normalizer

In the source Weyl algebra use [d,r]=[s,t]=h, with the other
generator commutators zero. Work over K(h)(r), or over a
characteristic-zero scalar Laurent-series field with r adjoined;
all plane expressions have finite PBW support. Let
\[
 \begin{split}
 \widehat\Delta
  &=4s^2-4rs^3+128t-72r(st+ts)-432r^2t^2,\\
 \widehat A
  &=\frac{2r(st+ts)-32t-s^2}{20r^2},\qquad
 D_*=d-\widehat A,\\
 X&=3rs-4,\qquad Y=108r^2t+18rs-16 .
 \end{split}
\]
These are symmetric orderings of QHR1's polynomials. Direct
commutator identities are
\[
 \begin{gathered}
 [X,Y]=324hr^3=:c,\qquad
 \widehat\Delta=-\frac{Y^2+4X^3}{27r^2},\\
 [D_*,r]=h,\quad
 [D_*,X]=\frac{6h}{5r}X,\quad
 [D_*,Y]=\frac{9h}{5r}Y,\quad
 [D_*,\widehat\Delta]=\frac{8h}{5r}\widehat\Delta .              \tag{4}
 \end{gathered}
\]
There are no omitted higher quantum terms. In the (s,t) plane,
A has degree at most two, so its normalized Weyl commutator
with any polynomial is exactly its Poisson bracket. The r
derivative comes from the generator d, which is linear.
The identity from
[QHR1](ELLIPTIC_RAMIFICATION_AND_STABLE_POLARIZATION.md)
therefore lifts exactly.

On the finite cover r=zeta^5, put
\[
        X_0=\zeta^{-6}X,\quad Y_0=\zeta^{-9}Y,\quad
                         \partial_\zeta=5\zeta^4D_* .
\]
Then
\[
 [X_0,Y_0]=324h,\qquad
 [\partial_\zeta,\zeta]=h,\qquad
 [\partial_\zeta,X_0]=[\partial_\zeta,Y_0]=0,\quad
 \widehat\Delta=-\zeta^8(Y_0^2+4X_0^3)/27.                     \tag{5}
\]
This is an exact ramified separation of the source connection,
not a polynomial canonical change on the original affine space.
It uses r denominators and an algebraic cover. It does not
supply the bounded full target images required by QFL2.

## A polynomial cubic oscillator in every finite degree

Let C be any characteristic-zero field, let c in C be nonzero,
and put W=C<X,Y>/([X,Y]-c), L=Y^2+4X^3.
Use weights wt(X)=2, wt(Y)=3. Products add weights, and
commutators lower weight by at least five. The graded Poisson
bracket is {X,Y}=c.

For ell=Y^2+4X^3 in the commutative graded plane, its Poisson
centralizer is C[ell]. One way to see this is to work on the
generic smooth elliptic level ell=lambda. A rational constant
of the nonzero Hamiltonian field is algebraic over C(lambda);
geometric integrality leaves only that field. The univariate
Bezout argument then gives the polynomial intersection C[ell].

Moreover
\[
 \{\{f,\ell\}:f\in C[X,Y]\}\cap C[\ell]=\{0\}.                  \tag{6}
\]
If the bracket were g(ell) nonzero, on the generic elliptic
curve df would be g(lambda) times the nonzero holomorphic
differential dual to that Hamiltonian field. A rational function
cannot have a nonzero holomorphic differential as its derivative:
a pole would persist in its derivative, and a pole-free
function on a projective curve is constant.

The same statements hold exactly in W:
\[
 \operatorname{Cent}_W(L)=C[L],\qquad
                 \{[F,L]:F\in W\}\cap C[L]=\{0\}.              \tag{7}
\]
Here is the finite all-degree argument. If the highest weighted
symbol of F commutes with ell, it is a scalar multiple of a
power of ell. Subtract the corresponding power of L and lower
weight. Repeat, finitely, until either F is a polynomial in L,
or its leading symbol f has nonzero bracket with ell.
In the latter case the leading commutator symbol is {f,ell}.
It cannot vanish or be a nonzero polynomial in ell by (6).
This proves both assertions without a bound on F's weight.

The result applies with C containing h and r and
c=324hr^3. It concerns finite operators over the scalar field,
not arbitrary formal spatial completions.

## The full finite source centralizer and normalizer

Let B be the finite source algebra with rational r coefficients,
written as C<X,Y>[D_*] with the relations in (4), where
C=K(h)(r) or K((h))(r). Nonzero polynomials in r are inverted. The d to D_* change is finite and invertible.
Then
\[
                  \operatorname{Cent}_B(\widehat\Delta)
                                      =C[\widehat\Delta].      \tag{8}
\]
For F=sum a_j D_*^j commuting with Delta, its leading
D_* coefficient commutes with Delta, so belongs to C[Delta]
by (7). If its D_* degree m is positive, the next coefficient is
\[
 [a_{m-1},\widehat\Delta]
       =-\frac{8mh}{5r}a_m\widehat\Delta .
\]
The right side is a nonzero polynomial in Delta. This contradicts
the second assertion of (7), so m=0, proving (8).

The same argument gives the complete first normalizer:
\[
 \begin{split}
 \{F\in B:[F,\widehat\Delta]\in C[\widehat\Delta]\}
     &=C[\widehat\Delta]D_*+C[\widehat\Delta],\\
 \{[F,\widehat\Delta]:F\in B\}\cap C[\widehat\Delta]
     &=\widehat\Delta\,C[\widehat\Delta].                       \tag{9}
 \end{split}
\]
A D_* degree at least two would give the same forbidden next
coefficient. At degree one, (7) forces the remaining plane
coefficient to commute; (4) gives exactly the multiples of Delta.
In particular no finite operator in this source localization is
a constant-commutator mate of Delta.

Additional spectator Weyl pairs can be included tensorially:
expand in their PBW monomials and apply (7) to every coefficient.
The centralizer in (8) acquires precisely the spectator factor,
and the normalizer has the same first-order D_* form.

Under an injective bounded lift these identities identify the
centralizer and normalizer *inside its localized image*. They
do not identify the whole ambient centralizer. Such a transfer
would require a separate theorem.

## Cubic energy bounds the orders of any finite pair

The following identities require only [U,V]=c, with c commuting
with both U and V. For P=V^2+4U^3,
\[
 [P,U]=-2cV,\qquad [P,V]=12cU^2,\qquad
                    [P,[P,U]]=-24c^2U^2.                     \tag{10}
\]
Suppose a degree function on a domain is additive on products,
subadditive on sums, and commutators lower degree by at least
e>0. Equations (10) imply
\[
 \deg U\le2(\deg P-e-\deg c),\qquad
 \deg V\le3(\deg P-e-\deg c).                                 \tag{11}
\]
Indeed the last identity gives
2deg(c)+2deg(U)<=2deg(P)+deg(U)-2e. The first identity
then gives the bound for V. These are necessary bounds for
finite elements; no existence of a pair is inferred.

For finite scalar differential operators, e=1 and c has
differential order zero, proving (1).
For the weighted source plane, e=5 and deg(L)=6. Therefore
any finite endomorphism of W fixing L sends X to degree at
most two and Y to degree at most three:
\[
                 U=aX+a_0,\qquad V=bY+b_1X+b_0.
\]
The commutator gives ab=1. Equality V^2+4U^3=L gives
b^2=a^3=1, hence a=b=1. The XY and Y coefficients force
b_1=b_0=0, and the X^2 coefficient forces a_0=0. Thus
\[
       \text{every finite plane endomorphism fixing L
                         is the identity.}                   \tag{12}
\]
No inverse was assumed.

Likewise a finite full source Weyl endomorphism fixing both
the source coordinate r and Delta fixes s,t by (12), after
localizing at r. Its last image is d+f(r), by the joint
centralizer of r,s,t. It is automatically an automorphism.
Fixing this source coordinate r is different from a proposed
Long lift assigning r to R; this stabilizer theorem is not
a nonlifting theorem.

## The quantum discriminant cannot remain quadratic for Long

Assume a full Long lift has one common spatial degree bound
and keeps R exact. Its fiber images S_h,T_h are finite scalar
operators over F=K(r)((h))(v), by QGR1's centralizer and
localization. Form the actual source expressions
\[
 U_h=3rS_h-4,\qquad
 V_h=108r^2T_h+18rS_h-16,\qquad
 P_h=V_h^2+4U_h^3=-27r^2\Phi_h(\widehat\Delta).                 \tag{13}
\]
They satisfy [U_h,V_h]=324hr^3.

If P_h had scalar order N<=2, (1) would give orders at most
two and three for U_h,V_h; if N<=1 it would already force
two commuting order-zero operators. But
[QPA1](RAMIFICATION_AND_POSITION_AVOIDANCE.md) excludes every
finite scalar Weyl pair in the entire completed fiber image
with either order at most two. Consequently
\[
 \boxed{\operatorname{ord}\Phi_h(\widehat\Delta)\ge3
        \quad\text{for every bounded exact-R Long lift.}}      \tag{14}
\]
The classical reduction of this operator has order exactly two,
since the displayed lambda in (2) is quadratic in a.
Thus even cancellations between arbitrarily high-order S_h and
T_h cannot leave a quadratic finite quantum discriminant.

This excludes quantizing the elliptic isogeny while preserving
its second-order energy, regardless of partner orders or
coefficient degrees. It does not exclude energy order three
or higher, a nonpositional deformation of the first image,
or an entirely different bounded full carrier.

## Why a formal energy bound still is not enough

The finiteness assumption in (11) cannot be discarded.
In W[[tau]], with tau an independent formal parameter, the
Hamiltonian flow exp(tau*ad(L)/c) fixes L and all Weyl relations.
Its images of X,Y have coefficients of unbounded spatial weight.
On the graded elliptic level, the n-th Hamiltonian derivative
of X has pole order 2+n at infinity: the dual differential is
holomorphic and nonzero there. Thus it never vanishes and has
weight 2+n. These are also the leading symbols of the
corresponding iterated normalized commutators.

The same example can be made h-adically small by using
exp(h^N*ad(L)/h) in the h-dependent Weyl algebra, for N>=1.
At the h order of the n-th classical iterate its highest
weight cannot cancel against the lower iterates. It is an
actual full formal canonical change, but not a finite plane
endomorphism. Equation (12) gives an independent reason
that a nontrivial finite alternative cannot represent that
same exact energy-preserving change.

The identity is still a finite map with this same energy. The
control refutes automatic boundedness of a chosen formal pair;
it does not refute an existence theorem for a different bounded
pair sharing its energy. Such an effectivity theorem for a
Long-compatible energy remains open. Here the established result
is the restriction on every candidate finite realization.

## Verification and next gate

[The checker](check_quantum_elliptic_energy.py) verifies the
three-isogeny identities with the correct sign, both elliptic
equations, the exact Ore commutators in (4), oscillator identities
(10), low-weight stabilizer coefficients and finite formal-flow
controls. It binds the retained QPA1 and QHR1 inputs by hashes.
See the [receipt](quantum_elliptic_energy_controls.json).

The all-degree centralizer, normalizer, degree bounds, finite
stabilizer and energy-order obstruction are written proofs.
Finite controls do not certify the absence of all finite
quantum lifts, and no such conclusion is claimed.

The next constructive route must allow higher-order quantum
energy and preserve all four images and all six relations.
The isogeny formula supplies geometry and exact normalizer
constraints, not a finite quantization. No independent theorem
replay, formal verification, external review, novelty priority,
bounded full Long lift, all-representative nonlifting theorem
or DC2 solution is claimed.
