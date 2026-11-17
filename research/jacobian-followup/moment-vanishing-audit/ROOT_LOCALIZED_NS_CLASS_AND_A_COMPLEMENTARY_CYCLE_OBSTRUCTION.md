# A root-localized NS class and a complementary cycle obstruction

17 November 2026. Written geometric obstruction to a proposed shortcut.
Both FULL exceptional Mathieu--Zhao cases remain **OPEN**; review remains
deferred. This control is an ACTUAL original PGL2 phase, but it is neither
pure nor right-O/I invariant. It does not refute the desired exceptional
pure-to-safe theorem or a theorem with genuine volume-regular hypotheses.

The [balanced ruling construction](RULING_CLASS_BALANCING_AND_THE_ORIGINAL_HAAR_SECTION.md)
localizes an UNBALANCED original NS representative above a chosen left line.
That support alone cannot localize the intrinsic corrected cycle. The example
below places the entire unbalanced degree on two rational tails, while the
entire nonzero component remains on the complementary four-cycle. Its
balancing divisor is explicit. Thus one cannot discard the compensation
and contract the complementary subgraph on NS-support grounds alone.

## 1. Use the actual original phase and move the chosen left line

On PGL2 take

    f=xw/Delta, Delta=xw-yz.

On SL2 this is the actual regular phase xw. For s!=0,1 its compact fibre is
the smooth quadric

    Y_s={(1-s)xw+s yz=0}.

Its original affine fibre is (C*)^2. The boundary in Y_s is the square of
four lines

    T={x=y=0}, B={z=w=0},
    R_0={x=z=0}, R_infty={y=w=0}.

The two row-zero lines have the same ruling class U; the column-zero
lines have ruling class V. Thus U^2=V^2=0, U.V=1. The original ruling
bundle A_Q=O_Q(1,0) has degrees (0,0,1,1) on this ordered square.
The actual boundary Picard section, after subtracting the balancing
class U, has holonomy s/(s-1), as computed from its node gluing in
the preceding note. Its original ordinary Cauchy period is

    C_1(s)=integral_0^1 dt/(s-t)=log(s/(s-1)).         (1)

Fix ANY kappa in C* and choose the left line

    ell_kappa={x=kappa z, y=kappa w} subset Q.

This line is distinct from T and B. Substitution into the equation of Y_s
gives kappa*z*w=0, independently of s. Hence ell_kappa meets Y_s transversely
at the two fixed points

    P_0=[0:kappa:0:1] on R_0,
    P_infty=[kappa:0:1:0] on R_infty.

Each is a smooth point of one boundary line, away from every square node.
Blowing up ell_kappa in P3 therefore restricts to

    Y'_s=Bl_(P_0,P_infty) Y_s.

This follows either in transverse local coordinates or by restricting the
centre ideal: on Y_s it is the ideal of these reduced smooth points.
The new boundary is the same square, with one exceptional rational tail
E_0 at R_0 and another E_infty at R_infty. Resolving the remaining pencil
indeterminacy preserves this generic fibre or makes further boundary
subdivisions; it does not remove its cycle.

## 2. Compute both the raw support and its compensating divisor

Use the basis U,V,E_0,E_infty of NS(Y'_s)_Q. The boundary component classes are

    T=U, B=U,
    R'_0=V-E_0, R'_infty=V-E_infty,
    E_0, E_infty.

The E_i have self-intersection -1 and are orthogonal to U,V and to each other.
In that order the surface boundary intersection matrix is

    S = [0 0  1  1  0  0]
        [0 0  1  1  0  0]
        [1 1 -1  0  1  0]
        [1 1  0 -1  0  1]
        [0 0  1  0 -1  0]
        [0 0  0  1  0 -1].                          (2)

It has rank 4. The pulled-back original ruling bundle has degree vector

    a=(0,0,1,1,0,0),
    beta=U, S*(1,0,0,0,0,0)^t=a.                    (3)

Let E denote the ambient exceptional divisor over ell_kappa. Its restriction
to Y'_s is E_0+E_infty. The root-localized representative

    A_root=A tensor O(-E)|D

is trivial on the strict transform Q'. Its new generic boundary degrees are

    a_root=(0,0,0,0,1,1),
    beta_root=U-E_0-E_infty,
    S*(1,0,0,0,-1,-1)^t=a_root.                      (4)

Thus the UNBALANCED ruling representative has degree supported entirely
on the two tails. The balancing class includes the square component T.
It cannot be confined to those tails: a divisor supported on E_0,E_infty
which has degree 1 on each tail also has nonzero degree on its attached
column line. Equations (2)--(4) make the compensation necessary.

## 3. The intrinsic component stays entirely on the complementary cycle

Let b:C'_s->C_s be the boundary contraction of the two tails. Pullback gives
an isomorphism on Pic^0: each tail is a rational curve attached at a single
point, so it contributes neither a Jacobian nor a graph cycle. These are
the ACTUAL fibre curves, not an ambient divisor graph.

Write L_old=A|C_s tensor O_Y(-T)|C_s. On the new boundary the balanced bundle is

    L_root=A_root|C'_s tensor O_Y'(-T+E_0+E_infty)|C'_s
          =b^*L_old.                                (5)

The equality is an equality of the actual line bundles, including node
gluing: the two ambient exceptional terms cancel. Therefore the unchanged
square cycle has

    holonomy s/(s-1), normalized local valuation 1 at s=0.  (6)

All generic normalizations are rational, so there is no abelian correction
which could hide this value. In graph terms the square and its two tails
have six vertices, six edges and first Betti number 1. The root-complement
square itself carries that whole H_1; the two tails carry zero H_1.

If one incorrectly contracted the root-complement square on the grounds
that A_root is trivial there, the resulting graph would have three vertices
and two edges, a tree. Equation (6) proves that such a contraction fails the
exact bounded-descent hypothesis. Its obstruction on the proposed contracted
subgraph is precisely 1, and cannot be removed by a bounded vertex gauge.

This holds for EVERY kappa in C*, not a selected numerical root or a scan
of degrees. The general logical failure is

    support of unbalanced original NS representative above root
       DOES NOT IMPLY
    corrected component supported on cycles above root.    (7)

The [bounded graph descent theorem](BOUNDED_GRAPH_DESCENT_AND_GOOD_NODE_EXTENSION.md)
requires vanishing of the actual cycle restriction before contraction;
(6) is an explicit original geometric obstruction that its input can consume.

## 4. The fixed complementary model has an actual nonfree-node obstruction

The four square lines and their nodes are constant as generic curves over
the phase parameter. Nevertheless their degree-zero line bundle varies.
Choose a cycle gauge putting three gluing scalars equal to 1 and the fourth
equal to g(s)=s/(s-1). At s=0 its remaining factor is s times an R-unit,
R=C[[s]]. On the fixed constant nodal square, after absorbing that unit,
the natural lattice at the last node is

    M={(u(x),v(y)):u(0)=s v(0)}
      over R[[x,y]]/(xy).                            (8)

It has generators (x,0),(s,1), both necessary at the closed node. Its closed
node fibre has dimension 2 and it is not locally free of rank 1 there.
The gluing-index cokernel is R/(s), of length 1. Thus an actual original
phase, rather than an arbitrarily assigned Picard point, realizes the
nonfree-node obstruction on this proposed complementary model.

The full proper phase resolution and its marked semistable models retain
extra special curves. The fixed complementary square model need not map
to that full boundary at s=0. One cannot use its generic constancy to claim
the actual raw bundle extends invertibly on this smaller model. On the full
model its exceptional multidegrees retain (6), in accordance with the
preceding raw-degree theorem.

## 5. Exact limit of the obstruction

This square is NOT volume-good. In z=1 coordinates, put y=xw(s-1)/s.
The ORIGINAL volume gives, with the fixed conventional sign,

    Omega=dx wedge dy wedge dw/Delta^2
         =-dlog x wedge dlog w wedge ds.              (9)

Hence the fibre two-form has logarithmic poles on the square boundary.
Blowing up either smooth point P_i on one logarithmic boundary component
adds a tail on which the fibre form is regular: locally
(dn/n) wedge dt pulls back under n=u,t-t_0=u v to du wedge dv.
Further resolving the phase pencil gives the corresponding horizontal
boundary components with that same regularity. Raw NS support and volume
pole support therefore must not be identified in this example.

Also h(f^m)=1/(m+1), m>=0, so f is not pure. Its top left tuple has two
simple roots rather than a unique strict dominant root. No exceptional
invariant pure phase or fixed original bad marker is constructed. In
particular this note does NOT disprove that purity together with the
exceptional capacity and volume-regularity inputs could force zero
corrected periods on genuinely good components.

The required exceptional step remains to use THOSE additional inputs to
control the compensating divisor in (4)'s general analogue. Triviality of
the unbalanced ruling bundle on a region, a tree for its visible support,
or constant generic boundary curves is insufficient. All original volume,
marker, branch, capacity, automorphism, primitive and family-transfer goals
and results are preserved. Both full cases remain OPEN; review is deferred.
