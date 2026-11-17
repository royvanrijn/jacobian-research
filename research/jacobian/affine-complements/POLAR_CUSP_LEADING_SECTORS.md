# Generic polar budgets and marked cusp sectors narrow pole(6,9)

14 November 2026. PF2RULED19 adds two actual-finiteness constraints to
[PF2RULED18's cover-period test](CUBIC_COVER_BOUNDARY_PERIOD.md).
Every root of a common leading section in L^d must have multiplicity
strictly greater than d/2. For gcd3 this forbids every simple common
nonic root, at the contact point as well as away from it.

For pole(6,9), contact multiplicities1 and8 are therefore excluded.
The contact polar budget also excludes the entire ordinary72,108 case.
Together with the previous contact2 exclusion, the remaining necessary
source groups are now:

| Nonic contact multiplicity k | Ordinary source degrees(P,Q) |
| --- | --- |
| 0 | (96,144) |
| 3,4,5,6,7,9 | (80,128) |
| 7,9 | (48,128) |

The subsequent [PF2RULED20 generic quartic polar gap](QUARTIC_GENERIC_POLAR_GAP.md)
excludes the elliptic222/contact3 profile and the noncontact partitions
(5,2,2),(3,2,2,2). The earlier five-partition and polar17/29 handoffs
below are historical necessary stages. Current noncontact types are
(7,2),(5,4),(4,3,2); the sixteen marked two-root leading candidates
and this note's general polar/cusp identities remain valid.

In the noncontact branch, only five noncube root partitions remain.
Its two-root sectors have exactly sixteen marked coefficient points
passing the cover period and the highest four cusp-jet equations. These
are leading candidates, not finite maps. Three/four-root positions and
the full Jacobian, cusp and finiteness equations remain open. No whole
height-nine pair or actual compatible-completion theorem is proved.

## A generic-polar budget at every boundary point

Let P have exact B-pole m, and let x be any root of its leading symbol,
with multiplicity l. At x choose regular coordinates w,z with B=(z=0)
and x=(0,0), and write

 A=z^m*P, F=A-c*z^m,
 A(0,w)=w^l*unit.

Choose a generic c, simultaneously outside the finitely many local
polar exceptions, fibre critical values and any value making E a
component. The generic P fibre is understood with all components.
The polar A_w is not divisible by z and

 I(z,A_w)=l-1.

On each normalized branch of A_w, genericity gives
ord(A-c*z^m)=min(ord A,m*ord z). Counting the polar branch
multiplicities therefore proves

 I(F,A_w)<=m*I(z,A_w)=m(l-1).                           (1)

This elementary pencil intersection argument includes repeated polar
branches. It needs no Newton-face simplicity or chosen Puiseux type.

Write the actual volume as e(w,z)*z^-3 dz wedge dw, with e=0 the
local part of E (or a unit away from q). Along a normalized branch of
F=0, the Jacobian identity gives

 dQ/dt=-e*z^(m-3)*(dz/dt)/A_w.

If z has order e_i and Q has pole order nu_i, comparison of orders gives

 nu_i=ord_i(A_w)-(m-2)*e_i-ord_i(e).

Finiteness is essential here: every missing point on the normalized
generic fibre is a pole of Q. Indeed the affine fibre is finite over
the Q-line, so a place above a finite Q value has its center in that
affine fibre. Over all components, all positive pole orders sum to4.
Each root x of the leading P-symbol supplies at least one missing
branch, since F is not divisible by z and meets z=0 at x.

Summing over the branches at x yields the exact identity and bound

 sum_x nu_i = I(F,A_w)-(m-2)*l-I(F,E)
            <=2l-m-I(F,E).                             (2)

The left side is a positive integer. In particular2l>m at every root,
including q. Also the entire leading symbol has support at most four,
because each distinct boundary point contributes at least one Q pole.
Both statements hold at every pole height under the actual finite-map
and pole-three-volume hypotheses. Meromorphic positive controls without
finiteness need not satisfy them.

If a common section h in L^d has root multiplicity k_i, the component
symbol h^(m/d) has l=m*k_i/d. Thus(2) forces

 2*k_i>d.                                              (3)

For d3, every nonic root has multiplicity>=2. For d2 this also forbids
a simple common sextic root at q, even though a pre-finiteness rational
period could allow it. The condition is a finite-map requirement, not
a correction to that period classification.

The checker retains balanced local polar controls attaining(1), and
the strict control A=w^6+z^5 with m6, for which the polar intersection
is25<30. The upper bound is not silently used as equality. The
all-height generic-polar and properness arguments above are written
proofs, not conclusions extrapolated from those controls.

## The contact correction and the72,108 exclusion

On this compatible completion E.B=4 and E meets B only at q. Put
v=D on E minus q, the fixed actual omitted-prime parameter. If
P0=P|E has degree d0 in v, its closure fibre has class mB and

 I_q(F,E)=4m-d0.                                       (4)

One may see this directly by restricting z^m(P-c) to E: z has
order four at q and P0 has a pole of order d0 there. Generic c also
handles d0=0. Equation(4) does not change the marked source parameter.

In the former pole(6,9)72,108 branch, PF2RULED18 gives nonic contact3,
so l=6, m=6, and d0=18. Hence I_q(F,E)=6 and(2) gives

 0<sum_q nu_i<=12-6-6=0,

an impossibility. This excludes that whole numerical branch. It does
not exclude the different72,108 resonances of the gcd-one height-nine
pairs or that ordinary-degree exception at arbitrary pole heights.

Nonic contact1 is a forbidden simple root at q. Contact8 leaves an
away simple root, also forbidden. The contact2 obstruction can also
be seen directly here: l=4 and I_q(F,E)>=4, so(2) has upper bound-2.
The smaller intersection bound follows by substituting z=w^4*unit
in A: its boundary term has order4 and every remaining term has
order at least4. Thus the three source groups in the first table are
the current necessary pole(6,9) handoff.

## The remaining root types and an exact four-place profile

The noncontact common nonic cannot be a cube, by the earlier rational
period and actual cusp obstruction. Combining(3) with its total degree9
leaves exactly

 (7,2), (5,4), (5,2,2), (4,3,2), (3,2,2,2).             (5)

Every such noncube must still pass the exact cubic-cover polynomial
test. No three/four-root position classification or realization is
inferred merely from(5).

At nonic contact3, its finite-root degree is6. The surviving possibilities
are a single root6 (the cube branch), two roots(4,2), or three roots
(2,2,2) with the equianharmonic root-position condition of PF2RULED18.
The other degree-six period configurations contain a forbidden simple
root. At contact6 only a single root3 survives; contacts4,5,7 have the
single finite roots5,4,2 respectively. Contact9 is the all-contact cube.

The three-(2,2,2) contact3 case has four distinct boundary roots including
q. Therefore the generic P fibre must have exactly four missing points,
one over each root, and all four Q poles must be simple. This is a
necessary profile, not an exclusion. In its source80,128 branch,
d0=20 and I_q(F,E)=4. Equation(2)'s exact polar identity forces

 I(F,A_w)=17 at each of the three away roots (l=4),
 I_q(F,A_w)=29 at q (l=6).

Possible further local obstructions to that profile remain to be proved;
no Puiseux classification or global Q is supplied here.

## Marked high cusp jets in the noncontact branch

Now keep the actual base frame from PF2WRIGHT1: r=s/t, q=(r=0),
g0=beta*t^2. Set w=1/r. This is fixed, so the marked negative-section
point is w=0. Write the monic common nonic as

 H(w)=w^9+a*w^8+b*w^7+c*w^6+d*w^5+...,
 h(r)=r^9*H(1/r).

No affine normalization of its roots is substituted for an automorphism
of the actual surface. Scale the target components so their symbols
are h^2,h^3. The actual cusp identity has derivative gcd b_R=v-v*,
polynomial multiplier lambda and lambda(v*)!=0, with
b_R*F1=lambda*F0' for every target-algebra element.

For noncontact symbols, PF2RULED14 gives

 lambda=(3/4)*b_R+delta, delta!=0,
 J(F)=F1-(3/4)*F0', b_R*J(F)=delta*F0'.

Here v is the omitted-prime coordinate, different from the base w.
For every section F in V_m, m6 or9, with leading binary coefficients
c0,...,c4, its normalized top four restriction coefficients are

 [c0,c1,c2-m*c0,c3-m*c1],

and its normalized top four corrected-jet coefficients are

 [-c1/4,-(m*c0+c2)/2,-(m*c1+3*c3)/4,
  -c4+m*(m-1)*c0/2].

The normalization divides by(-1)^m; it cancels from the target identity.
The checker verifies these linear identities on all70 and145 Cox
sections, respectively. Lower-pole terms cannot affect these four
positions, and no lift or marked root position is omitted.

Substitute the binary coefficients of h^2 and h^3 and compare the
highest four coefficients in b_R*J(F)=delta*F0'. They force

 delta=-a/48, a!=0,
 v*=(24*b+72-11*a^2)/(12*a),
 E=a^2*b+11*a^2+3*a*c-4*b^2-24*b-36=0,
 F=a^3*b+16*a^3+5*a^2*c-4*a*b^2-30*a*b
   +4*a*d-66*a-6*b*c-18*c=0.                           (6)

After the first two substitutions, the pole-six residuals are
-E/(2a),-F/(2a); the pole-nine residuals are
-3E/(4a),-3(F+aE)/(4a). Thus(6) exactly characterizes this highest
four-coefficient block for both components. It is necessary, but not
sufficient, for the full cusp-jet identity or a map.

## Exactly sixteen marked candidates in the two-root sector

For each of the two root types(7,2),(5,4), write

 H=(w-x)^k*(w-y)^(9-k), k=2 or4.

The root positions x,y are actual fixed marked coordinates, not roots
normalized to0,1. These two-root noncube periods pass at arbitrary
distinct positions by PF2RULED18. Insert their top coefficients in(6).
The exact resultant in x factors into a nonzero constant times

 y*q2(y)*P8(y),

with the following rational univariate representations:

| k | q2(y), forbidden zero-cusp-unit branch | P8(y) |
| --- | --- | --- |
| 2 | 21y^2-4 | 1225y^8-1568y^6+5600y^4-13440y^2+8192 |
| 4 | 15y^2-8 | 25y^8+2800y^6-4480y^4+1920y^2-256 |

For k2,

 x=-y*(37975*y^6+4592*y^4+217504*y^2-253824)/123904.

For k4,

 x=-y*(65*y^6+7180*y^4-22864*y^2+15936)/5888.

The checker retains the full resultant, the linear subresultant and
its invertible leading coefficient modulo P8. At y=0 an exact Bezout
identity excludes solutions. Over q2=0 the exact Groebner basis contains
k*x+(9-k)*y=-a, so every solution has the forbidden a=0. Thus every
allowed point must satisfy P8=0 and the displayed x formula.

Conversely, substitution into E,F reduces to zero modulo P8. Both P8
polynomials are squarefree of degree8. Exact Bezout certificates make
a, x-y, x and y units modulo P8, so every one of the eight roots in
each sector has distinct nonzero positions and nonzero cusp unit.
The two sectors have different multiplicities and are disjoint. This
proves exactly sixteen marked coefficient points satisfying the period
and highest four cusp-jet equations in this two-root sector. It does
not count surface-isomorphism classes, and none is asserted to extend
to full P,Q, zero fibre residues or a finite map.

## Assurance and remaining gates

The [checker](verify_polar_cusp_sectors.py),
[certificates](polar_cusp_sector_certificates.json) and
[receipt](polar_cusp_sector_receipt.json) retain the complete finite
jet basis, every elimination and unit certificate, balanced/strict
polar controls and the remaining numerical frontier. The all-height
polar intersection, fibre properness and contact-intersection arguments
are written proofs. Shared exact polynomial utilities are used; no
formal verification, independent full reproof, conditional2026 upgrade,
external review or novelty priority is claimed.

The first open pairs remain(5,9),(6,9),(7,9),(8,9). The sixteen marked
two-root points, the three/four-root noncontact sectors, the retained
contact cases and four-place profile still require actual global
Jacobian, full cusp and finiteness analysis. Actual compatible embedding
of an arbitrary quartic normalization, higher heights, all-height
pole-four coordinates and both full quartic packets remain open.
The original full goal remains active.

The elliptic polar17/29 realization gate is closed by PF2RULED20.
The next surviving local profiles require polar33/25/18 in the
noncontact432 case, or33/18 away and29 at q in the contact42 case.
These higher-polar profiles and the sixteen two-root leading points
remain open; no whole pole69 pair or actual compatibility is settled.
