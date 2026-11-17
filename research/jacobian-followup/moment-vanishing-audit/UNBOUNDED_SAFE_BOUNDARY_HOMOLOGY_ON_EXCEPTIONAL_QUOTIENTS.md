# Unbounded safe boundary homology on the exceptional quotients

16 November 2026. Internal written proof; independent replay, formal
verification and external review are not claimed. Review remains deferred.
Both COMPLETE exceptional Mathieu--Zhao cases remain **OPEN**.

The [cycle/linking theorem](HAAR_CYCLE_VARIATION_AND_FIBRE_BOUNDARY_LINKING.md)
gives a sufficient condition using semisimple fibre monodromy at 1 and
vanishing invariant boundary image. The second hypothesis CANNOT hold for
every pure phase on either exceptional quotient. In fact there are pure
safe phases with identity fibre AND complement monodromy, yet arbitrarily
large invariant boundary images. The actual variation cycle is zero.

This closes a proposed universal hypothesis, not a new sector of the MZ
problem. Single-weight safety was already known. The new result is the
unbounded topology of actual invariant safe phases and the obstruction to
using boundary-space dimension in place of the actual Haar-cycle coupling.
See the [checker](verify_haar_boundary_excess.py),
[receipt](haar_boundary_excess_verification.json) and
[source audit](haar_boundary_excess_source_audit.json).

## 1. An all-degree existence theorem

Let Gamma=2O or 2I, G=Gamma/{+-I}, and N=|G|. Put

    (N,e)=(24,12) for 2O, (60,30) for 2I.                (1)

Here e is the exponent of the projective rotation group. For EVERY
integer k divisible by e with k>=2N, there is a nonempty Zariski-open
set of Gamma-invariant bihomogeneous forms

    H(x,y;z,w), bidegree (k+1,k-1),                      (2)

whose restrictions P to xw-yz=1 have these properties on X=SL2/Gamma:

* P has matrix-filtration degree exactly 2k and left torus weight 2.
* P is pure and safe for EVERY original fixed marker Q, with the
  sufficient cutoff 2m>deg(Q), independent of k.
* Over C* the fibre AND complement families are algebraically trivial.
  Thus T_F=T_U=I and the actual Haar-cycle variation d_P is zero.
* A smooth projective compactification of F_1={P=1} has a boundary
  component of genus

      g=k(k-2)/N.

  The nearby-fibre intersection radical therefore satisfies

      dim_Q Rad(B_F)^T_F >= 2k(k-2)/N.                  (3)

In particular the invariant boundary image is unbounded on the PURE SAFE
locus of each WHOLE exceptional algebra. The explicit threshold k>=2N
is a conservative interpolation bound, not a computed least degree or
a degree search. We do not output expanded coefficients of a generic H.
The existence statement, including nonempty openness, is proved below.

## 2. Invariant first jets and a smooth curve

Let Y_0=P1 times P1, with the diagonal projective G action. The line bundle
O(k+1,k-1) has the natural binary-group linearization. Its total degree
2k is even, so the central involution acts trivially on sections and the
linearization descends to G. Write V for its invariant section space.

First record an elementary interpolation bound. If a,b>=2N-1, then
H^0(Y_0,O(a,b)) surjects onto first jets at any orbit of at most N points.
For each target point p and every other orbit point q, choose a linear
form in one factor which vanishes at q and not at p. The product of their
squares has bidegree at most (2(N-1),2(N-1)), is a unit at p, and vanishes
to order two at all other points. Multiply by sections of the remaining
positive bidegrees to prescribe an arbitrary value and two first
derivatives at p. Sum over the orbit. Averaging over G shows that V
surjects onto every stabilizer-compatible first-jet datum at p.

For k>=2N, both degrees in (2) meet this bound. The nonfree locus of
G on Y_0 is finite. A nontrivial element of PGL2 fixes exactly two points,
and a finite point stabilizer is cyclic: after fixing infinity its finite
affine action has a common second fixed point, and is conjugate to a
finite multiplicative subgroup.

Let S_off be the nonfree points off the diagonal. They are the ordered
pairs of the two fixed points of a nontrivial cyclic stabilizer. If a
binary lift of a generator has eigenvalues lambda,lambda^-1, local
coordinates xi,eta at such a pair transform by lambda^-2,lambda^2.
The value multiplier of (2) is lambda^2, so an invariant section must
vanish at this point. However its xi first jet is allowed, because
lambda^2*lambda^-2=1. Interpolation therefore permits a nonzero first
derivative. A generic invariant section is smooth at every point of S_off.

At a nonfree diagonal point the value multiplier is lambda^(2k)=1,
since e divides k. A generic invariant section avoids these finitely
many points. Elsewhere the stabilizer is trivial and invariant first
jets are arbitrary. A singular zero imposes three independent linear
conditions (value and two derivatives), while the point varies in a
surface. Their incidence therefore has codimension at least one in V.
Tangency to the diagonal imposes two independent conditions with a point
varying on a curve, and has the same codimension bound.

Consequently a nonempty Zariski-open subset of V defines a smooth curve
C_H, transverse to the diagonal, avoiding its nonfree points and passing
smoothly through every point of S_off. This is a first-jet incidence
proof, not a claim based on sampled coefficients. The curve is irreducible:
if a smooth divisor of positive bidegrees were disconnected, its distinct
components would have zero pairwise intersection. Effective bidegrees
then force all components to be fibres of the same ruling, contradicting
positivity of both total degrees.

The genus calculation is elementary. For a smooth curve of bidegree
(a,b), the divisor exact sequence and
chi(O(i,j))=(i+1)(j+1) give genus (a-1)(b-1). Thus

    g(C_H)=k(k-2).                                       (4)

The construction does not make C_H contain the diagonal. Its defining
H is therefore not divisible by xw-yz, which also proves the asserted
matrix-filtration degree 2k after imposing xw-yz=1.

## 3. The exact quotient genus

The G action on C_H is faithful: any nonidentity element has only finitely
many fixed points on Y_0 and cannot fix this curve pointwise. Its only
ramification points lie in S_off. These points are in bijection with
the nonfree points u of P1, by mapping u to (u,v_u), where v_u is the
other fixed point of its cyclic stabilizer. Stabilizer orders agree.

Riemann--Hurwitz for P1->P1/G=P1 therefore gives the total ramification
on C_H as

    R=sum_(p in S_off)(|G_p|-1)=2N-2.                   (5)

Applying Riemann--Hurwitz to C_H->C_H/G gives

    2k(k-2)-2=N(2g(C_H/G)-2)+(2N-2),
    g(C_H/G)=k(k-2)/N.                                 (6)

This is an equality for EVERY k in the theorem. The quotient curve is
smooth as an abstract projective curve; its strict transform in a surface
resolution has the same genus. The characteristic-zero formula used
here is [Stacks Project, Riemann--Hurwitz](https://stacks.math.columbia.edu/tag/0C1B).

## 4. The nearby fibre is the complement of that curve

The map taking the two row lines identifies the left-torus quotient

    T\SL2 = (P1 times P1) minus diagonal.                 (7)

Distinct row lines admit determinant-one representatives, unique up to
the left torus. Quotienting on the right gives the surface

    B=((P1 times P1) minus diagonal)/G.

All left weights in C[X] are even. The effective left torus T/{+-I}
is an algebraic G_m, with weight half the original weight. The phase P
in (2) has effective weight ONE. Hence

    X minus {P=0} = F_1 times G_m,
    (y,a) -> a*y,  x -> (P(x)^-1*x,P(x)).               (8)

The quotient of the left side by G_m identifies F_1 with B minus the
image of C_H. All nonfree points of the G action on the surface lie on
that curve; thus this complement is smooth, as F_1 is. The same action
trivializes the entire fibre family F_s and complement family U_s over
C*: scale the X coordinate by s. This proves T_F=T_U=I and d_P=0,
without an abstract choice of a different lift.

Let Z=(P1 times P1)/G. It is a rational projective surface. Indeed the
first projection gives Z a fibration over P1/G=P1 with generic fibre a
form of P1. The diagonal gives it a rational section. A smooth genus-zero
curve with a rational point is P1 over the base field, so the function
field is C(t)(v). See
[Stacks Project, characterization of P1](https://stacks.math.columbia.edu/tag/0C6U).

Resolve Z and the boundary consisting of the diagonal image and curve
image. All quotient singularities have cyclic stabilizers and lie in
the removed boundary; resolving and then blowing up boundary intersections
does not change F_1. This produces a smooth projective rational surface
Y with an SNC boundary D and F_1=Y minus D. The strict transform C of
C_H/G is a boundary component of genus g from (6). In particular

    H_3(Y;Q)=0.                                         (9)

## 5. The genus produces independent boundary tori in the radical

Here is a direct topological proof of the lower bound (3). Let N(D) be
a small plumbing neighbourhood of the SNC boundary, let L be its boundary,
and K the complementary compact core of F_1. Choose 2g standard homology
loops on C away from the finitely many intersections with other boundary
components. Their inverse images in the normal circle bundle are tori
t_1,...,t_(2g) in L.

These tori are independent in H_2(L;Q). To see this, take loops on C
dual under its symplectic intersection pairing and lift each to a closed
loop in the circle bundle. A circle bundle over a loop has a section.
Their intersections with the t_i give the nonsingular surface intersection
matrix. All loops avoid the plumbing points, so the calculation occurs
in the ordinary circle-bundle piece of L.

Each t_i maps to zero in H_2(N(D);Q): it bounds the normal disc bundle
over its base loop. Mayer--Vietoris for Y=K union N(D), together with (9),
makes

    H_2(L;Q) -> H_2(K;Q) direct_sum H_2(N(D);Q)           (10)

injective. Therefore the images of all 2g tori remain independent in
H_2(K;Q)=H_2(F_1;Q). They come from the boundary, so belong to the
intersection radical by Poincare--Lefschetz duality and the pair exact
sequence, as in the retained linking theorem. This proves

    dim Rad(B_F)>=2g.

Since T_F=I, this is already the INVARIANT radical. This proof uses neither
a moment computation nor an identification of a generic boundary space
with the original unit's V-filtration defect.

## 6. Purity, safety and the implication that must remain open

Under diag(t,t^-1), the phase P has weight 2. Haar kills every nonzero
weight, so h(P^m)=0 for all m>=1. An original marker of matrix degree q
has weights in [-q,q]. All terms of QP^m have positive weight when 2m>q,
proving the uniform marked cutoff. This reuses the completed weight theorem;
it is not new evidence for its safety conclusion.

Combining this with (3) gives actual safe phases with arbitrarily large
boundary homology on BOTH exceptional quotients. For instance the first
degrees covered by the conservative construction are k=48 on 2O and
k=120 on 2I, with boundary-radical lower bounds 184 and 472 respectively.
These numbers are formula evaluations, not expanded examples or least
thresholds. The theorem holds for all the indicated multiples of k.

The full MZ target remains vanishing of the ACTUAL d_(P,1) for EVERY pure
phase. Proving Rad(B_F)^T_F=0 on the entire pure locus is now ruled out.
One must control how the Haar lift couples to the boundary space, rather
than bound or eliminate that space. The sufficient theorem from the
linking note remains correct, but is not an exhaustive characterization.
No right-2O/2I character argument or boundary-rank test supplies the missing
cycle vanishing. Both complete cases stay open, with review deferred and
no new degree, support, density or prefix campaign.

The constructed phases lie in the already closed one-weight sector. This
does not disprove a boundary-vanishing statement restricted to some further
unresolved candidate locus. No such restricted theorem is established here.
The barrier concerns a uniform argument on the entire pure locus, and shows
why the actual-cycle target must be retained when applying it.
