# Odd-degree images force affine intersections and singleton packets

14 November 2026. PF2OMIT15 adds an intersection obligation for the
unknown images in the four residual rows of PF2OMIT13. It uses the
actual gap-21 Laurent expansion, rather than assuming that two images
with the same infinity point must meet in the affine plane.

## 1. Intersection theorem

Let C have a birational polynomial normalization (p(t),q(t)) of degrees
(6,10), normalized as in [PF2R7FAM1](../jc2-r7-image-family/PROOF.md).
In the formal coordinate s=p(t)^(1/6), with s=t+O(t^-1), its expansion
has leading term s^10, no odd exponent above -11, and a nonzero
coefficient gamma at s^-11. This is the meaning of exact gap21 in
[PF2D6O1](../../foundations/plane-jc/F2_DEGREE_6_10_POLYNOMIAL_GAP.md).
The sign of gamma depends on the obstruction convention; only gamma!=0
is used.

Let B be any other reduced irreducible polynomial image with birational
normalization (a(v),b(v)) of degrees (3k,5k), where k is positive and
odd. This applies to the second ramified image and also to a distinct
unramified omitted image. Normalize a to be monic by scaling v, and
write c!=0 for the leading coefficient of b. Let H_B(U,V) be the reduced
implicit equation, monic of V-degree3k. Then

    R_CB(t)=H_B(p(t),q(t)) != 0,
    9k <= deg R_CB <= 30k,
    20k <= I_infinity(C,B) <= 41k.                 (O1)

In particular C and B intersect in the affine plane. Their total affine
intersection multiplicity is deg R_CB, hence at least9k. This is not a
lower bound for the number of distinct intersection values.

Here is a formal norm proof, valid for every odd k, without a bounded
coefficient search. Put s=z^k and w=z^2. The equation

    a(v)=s^6=w^(3k)

has3k distinct formal roots v_j in C((w^-1)), with leading terms
zeta_j*w, zeta_j^(3k)=1. This follows from the formal implicit function
theorem after dividing v by w; the leading polynomial has simple roots.
Thus every b(v_j) has only even powers of z and leading coefficient
c*zeta_j^(5k) at exponent10k. The values zeta_j^(5k) comprise the three
cube roots of unity, each with multiplicity k.

The expansion q(t(s)) at s=z^k has first odd z-exponent -11k, with
coefficient gamma. Since k is odd, no b(v_j) can cancel that coefficient.
Every factor q(t(z^k))-b(v_j) therefore has highest exponent at least
-11k. If c^3!=1, all3k factors have highest exponent10k. Otherwise
exactly2k factors have highest exponent10k and the remaining k have
highest exponent between -11k and10k. The product of these factors is
H_B(s^6,q(t(s))). Pole orders add, so its highest z-exponent is between

    (2k)*(10k)+k*(-11k)=9k^2 and (3k)*(10k)=30k^2.

Primitivity of the parametrization is essential: the norm is the reduced
equation, rather than a power of it. The polynomial R_CB(t), evaluated
at t=t(s), has highest s-exponent deg R_CB because t(s)~s. Its highest
z-exponent is consequently k*deg R_CB, proving the degree bounds.
The product cannot vanish since each factor has the uncancellable term.
Also B and C have different degrees when k is odd.

The normalization maps extend to P1. Each has its sole infinity point
at [0:1:0], and deg C=10, deg B=5k. Pulling the degree5k homogeneous
equation of B back to C's normalization gives total degree50k, with
finite zero divisor R_CB. This proves the infinity bounds and identifies
the finite multiplicities with local intersection multiplicities,
including when the images have affine singularities.

For k=1 the possible degree list tightens to

    deg R_CB in {9,10,12,14,16,18,20,22,24,26,28,30}. (O2)

Indeed after leading cancellation the one matching factor has highest
s-exponent in {8,6,4,2,0,-2,-4,-6,-8,-10,-11}; the other two contribute20.
Without leading cancellation the degree is30. This is a list of necessary
possibilities, not a realization assertion for every entry.

For even k the parity argument fails: -11k is an even z-exponent and
can occur in b(v_j). No corresponding even-k affine-intersection lower
bound is asserted here.

## 2. Which intersections are necessarily singleton values

Apply (O1) to the actual B in the four residual rows. The complete
[length ledger](RESIDUE_LEDGER.md) and the
[nonimmersed correction](NONIMMERSED_LEDGER.md) give:

| Row/grouping | Retained length at C/B intersection | Forced branch counts |
|---|---|---|
| A, H=C | 6-3c-2b | c=b=1, retained1 |
| A, H=B | 6-2c-3b | c=b=1, retained1 |
| C, index-three S | 6-2c-3b | c=b=1, retained1 |

At an intersection c,b>=1. Nonnegative length forces the displayed
unique solution in each row. Both image branches are smooth there:
each ramified image is immersive in these rows, and each has only one
normalization parameter at the value. Transversality does not follow.

Thus if k is odd, each of these three groupings has at least one
singleton source value, and the sum of tangency/intersection multiplicities
at all of its C/B singleton values lies in[9k,30k]. The omission packets
are distinct even when their images agree:

| Row/grouping | Original omitted local lengths | Radial omitted indices |
|---|---|---|
| A, H=C or H=B | 2+2+1 | 2+2+1 |
| C | 2+3 | 2+3 |

Each packet has one additional retained reduced point and all residue
degrees one. The length-one omitted point in row A is still a separate
point of the original normalization. Dropping it would incorrectly
give retained length2. PF2OMIT10's A2 or C inverse bounds apply at every
one of these values, not merely at an unspecified possible singleton.

In rows B,N, a smooth C/B intersection can instead have retained length2.
The ordinary cusp correction can make an intersection a singleton but
its occurrence is not forced by (O1). Row A with three distinct images
also requires the independent H branch count. Consequently this theorem
does not force singleton values in those other rows.

For source degrees(75,125), the [pole budget](GENERIC_LINE_ACCOUNTING.md)
limits the forced-singleton odd choices to k=1,3,...,19 for H=C,
k=1,3,...,11 for H=B, and k=1,3,...,21 for row C. These are conditional
degree labels, not constructed covers. The intersection theorem itself
has no source-degree cap.

## 3. Exact controls and the remaining global calculation

The [checker](verify_odd_intersections.py) constructs reduced implicit
norms and intersection polynomials for explicit target curves. With
the L,d=0 first image, the controls B=(v^3,v^5) and (v^3,2v^5) give
degrees24 and30. The birational odd k=3 control (v^9+v,v^15) gives74;
doubling its second coordinate gives90. These are target curves only.

A closer immersed quintic control is defined over a^3=-7/75 by

    a_B(v)=v^3+a*v-2/15,
    b_B(v)=v^5+(5a/3)*v^3+(5a^2/9)*v.

The exact resultant composed with the same C has degree12 and its
derivative resultant is a unit in that cubic algebra. The replay uses
coefficient reductions, not approximate complex roots. This control
does not assert that B avoids every forbidden C value or extends to
any finite normalization with source A2. The bound9k is not asserted
sharp by these examples.

Once an actual B is provided, compute R_CB by the multiplication norm
of b(v)-q(t) in C[t,v]/(a(v)-p(t)). Factor its nonzero one-variable
polynomial, compute the normalization parameters and branch counts
over each value, and apply the full omitted-prime length formula.
Multiple roots record tangency, not additional omitted points. This
gives a finite exact interface for an actual missing exceptional locus.

No actual B,H equations, infinity divisors, global cover or source
identification are supplied here. Existence of a singleton does not
prove properness. A target translation can move such a value to zero,
but the singleton-zero properness question remains open.
