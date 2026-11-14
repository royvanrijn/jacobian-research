# Quartic contact exclusion of both cusp classes without source finality

14 November 2026. PF2R7Q1 is an exact algebraic obstruction: neither
retained PF2R7CUSP1 polynomial pair admits a conic tangent at the four
generic points of any quartic polynomial parameter map. PF2R7EX2 combines
it with the cubic-cover results to exclude both nonimmersed normal-r7
classes under the entire single ramified affine-prime hypothesis,
**without any conclusion of Borisov's Theorem 3.8 or its finality step**.

The earlier boundary theorem, its type1 residue and finite cyclic-quotient
coexceptional results, and minimal-resolution positivity remain explicit
written inputs. No formalization, external review or novelty priority is
claimed. The immersed multiple-root and general quartic family, additional
ramification patterns, normal r5, nonnormal slices and JC2 remain open.

## 1. The geometric case now to be excluded

Assume a complex plane Keller map has geometric degree six, and the
entire affine finite-normalization ramification divisor is one prime E,
of generic (e,f)=(2,1), with branch image one of the two exact retained
nonimmersed polynomial degree-(6,10), gap-21 classes. All unramified
omitted primes are retained.

The retained classification uses affine parameter changes and affine
target translations, scalings and shears q->q+C*p. These preserve conics,
tangency, the target line at infinity and quartic polynomial base maps.
Excluding each exact representative therefore excludes its whole stated
coefficient-normalization class.

PF2R7INF1 supplies the rank-eight rational pencil F=K+H on the minimal
smooth source Y, H.F=2, and its three possible horizontal-infinity
budgets. PF2R7V1 excludes the varying-direction parabola budget.
PF2R7G1 forces a multiple root of p' in the reciprocal budget, whereas
both exact classes have literal squarefreeness certificates. Thus only
the infinity pair (1,5) with fibre degrees (2,0) can remain, and

    deg(E^nu -> pencil base)=4.                       (1)

The general pencil fibre maps isomorphically onto a smooth target conic:
its two points on the index-one infinity curve have distinct images on
the target line, so its image cannot be a line. Adding the pencil
parameter to C(x,y) therefore generates the whole degree-six source
field. A primitive conic equation has degree six in that parameter, by
Gauss irreducibility as proved in PF2R7V1.

## 2. The quartic parameter on E must be a polynomial after inversion

Every exceptional curve over the common boundary point A is vertical:
the index-one infinity curve alone accounts for H.F=1*2=2. The affine
normalization is smooth by PF2R7INF1, so all exceptional curves lie over
A, and that connected exceptional fibre has a single pencil value z0.
The infinity point of E and the unique point of the index-one curve
mapping to a=[0:1:0] meet this fibre. Consequently their pencil values
agree. Choose a base coordinate z with z0=0.

On the index-one infinity curve, g is an isomorphism onto target infinity,
and the pencil has degree two. Its incidence with target directions is
therefore the graph of a degree-two morphism P1 -> P1. The quadratic
homogeneous part of the primitive conic equation is

    h(z)*(A2(x,y)-z*B2(x,y)),

where A2,B2 are coprime binary quadratics, A2(0,1)=0 and B2(0,1)!=0.
Possible special vertical factors are retained in the polynomial h.
The full equation is

    h(z)*(A2-z*B2)+l(z)*x+m(z)*y+n(z)=0,
    deg(h)<=5, deg(l),deg(m),deg(n)<=6.                (2)

These bounds follow from the actual parameter degree six of the primitive
equation, not an assumed bounded family. The coefficient of y^2 in the
quadratic part is a nonzero multiple of z*h(z), proving the bound on h.

Let k be the order of z(t) at t=infinity on E^nu; by (1), 1<=k<=4.
Suppose k<4. Since p,q have pole orders six and ten, A2(p,q) has pole
order at most16, whereas z*B2(p,q) has pole order20-k>16. Let h0,l0,m0
be the nonnegative vanishing orders of h,l,m at z=0. The leading orders
of the three terms in (2), apart from regular n, are

    k*(h0+1)-20, k*l0-6, k*m0-10.                    (3)

Tangency of the conic to the branch image gives the same three orders
plus six in F_x+s*F_y=0, with s=q'/p', and their leading coefficient
ratios to incidence are respectively 10/3,1,5/3. These are distinct.
The first order in (3) is negative, since h0<=5 and k<=3. Therefore the
minimum is negative; regular n cannot cancel it. A single leading term
is impossible, and any two cannot cancel both equations with distinct
ratios. All three must have the same order. This requires

    k*(h0+1-l0)=14, k*(h0+1-m0)=10.

Hence k divides both14 and10, so k is one or two. But h0+1<=6, while
the first equation requires h0+1-l0>=7 even for k=2. This contradiction
excludes every k<4.

Thus z has its unique zero at infinity, of order four. Since z(t) has
degree four, z(t)=constant/Z(t) with Z a quartic polynomial. Scaling the
base and subtracting its constant term puts

    Z(t)=t^4+u*t^3+v*t^2+w*t.                         (4)

No nonlinear target transformation has been used. The argument in this
section applies to every gap-21 image satisfying the quartic escape,
including immersed images, not only the two cusp examples.

## 3. Four contacts force a rank drop modulo a doubled quartic

Use a separate transcendental base symbol Z and let

    W(t)=t^4+u*t^3+v*t^2+w*t-Z.

It is separable over C(Z). A conic with coefficients in C(Z) through the
four conjugate points (p(t),q(t)) and tangent there gives

    W(t)^2 divides c0+c1*p+c2*q+c3*p^2+c4*p*q+c5*q^2
    in C(Z)[t],                                      (5)

for a nonzero coefficient vector. Indeed incidence gives divisibility
by W; tangency gives the same for the t derivative, which is equivalent
to a double factor because W is separable. This is a necessary condition;
no sufficiency or global Keller-source reconstruction is inferred.

Reduce the six columns (1,p,q,p^2,p*q,q^2) modulo W^2 in the ordered basis
1,t,...,t^7. The first and p columns pivot in rows zero and six because
p is monic of degree six. After those pivots, retain rows 1,2,3,4,5,7
and columns q,p^2,p*q,q^2 to obtain a 6-by-4 matrix M. Thus (5) requires
every four-row determinant of M to vanish identically as a polynomial
in Z. Determinant row labels below are the zero-based indices in this
six-row matrix, not powers of t.

## 4. Exact coefficient certificates exclude every quartic (4)

For either exact cusp class, the determinant on rows (0,2,3,4) satisfies

    [Z^10] det M_(0,2,3,4)=-80*u^5.

Hence u=0. At u=0, the determinant on rows (0,1,2,3) has coefficient

    class I:  [Z^10] = -4*(90*v+31)*(9*w+4)^2/243,
    class II: [Z^10] = -4*(90*v+157)*(45*w+8)^2/6075.

Over C this leaves just two branches in each class:

| Class | Fixed v branch | Fixed w branch |
|---|---|---|
| I | v=-31/90 | w=-4/9 |
| II | v=-157/90 | w=-8/45 |

On each branch, two further coefficients of that same determinant form
a unit ideal in the one remaining parameter. Their powers of Z are
respectively (9,8),(9,7),(9,8),(9,8) in the table's class/branch order.
[certificates.json](certificates.json) freezes both coefficient arrays
and their rational Bezout multipliers. Literal multiplication gives

    U(parameter)*f(parameter)+V(parameter)*g(parameter)=1.

The checker reconstructs the selected determinants and coefficients from
the original exact p,q and W^2, then separately replays each unit identity
using dense Fraction polynomial multiplication. Thus there are no complex
parameters u,v,w satisfying (5) in either class. No numerical roots,
finite coefficient slice, modular lift or fresh polynomial gcd is trusted
as the exclusion certificate. This proves the algebraic claim PF2R7Q1.

## 5. Scoped geometric conclusion and its precise inputs

Section 2 forces (4) in the only remaining infinity escape; Sections 3--4
exclude it for both cusp classes. Combined with PF2R7INF1, PF2R7V1 and
PF2R7G1, this proves PF2R7EX2: neither nonimmersed coefficient class can
occur in a degree-six Keller map under the entire unique ramified affine
prime (e,f)=(2,1) hypothesis.

Unlike PF2R7EX1, this deduction does not import source finality or the
singularity-at-infinity conclusion of Borisov's Theorem 3.8. The common
boundary point, type1 residue degree one, Theorem 3.6's finite
coexceptional-curve requirement, minimal-resolution positivity and the
exact two-class classification remain explicitly imported. The older
source-relative proof and its certificates are preserved as a separate
route; they are not retroactively described as independently reviewed.

The entire unique ramified-prime hypothesis is still not derived for
every normal-r7 terminal passport. The immersed quartic case has (4) as
a new necessary condition but has not been excluded by the two exact
class checks. The special immersed reciprocal case with a multiple
root of p' also remains. No whole normal-r7, degree-frontier or JC2
closure is claimed.

## 6. Replay boundary

[verify.py](verify.py) recomputes the selected exact minor coefficients,
the four univariate unit certificates and the bounded k<4 order
obstruction. A conic curve control with q=p^2 confirms the contact-rank
test can accept a genuine quadratic relation; it is explicitly outside
the birational degree-(6,10) hypotheses. Mutating a unit certificate
fails the separate Fraction replay. [receipt.json](receipt.json) binds
the local sources and original class data.

Irreducible conic/field identification, the infinity incidence graph,
contact geometry and polynomial inversion of the totally ramified map
are written proof obligations. The algebraic checker does not formalize
them or independently reprove the earlier surface theorems.
