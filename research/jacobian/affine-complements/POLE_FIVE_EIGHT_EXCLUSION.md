# Height-eight coprime exclusions: (5,8) closes, (7,8) retains contact

14 November 2026. PF2RULED12 excludes every finite quartic map with
coordinate poles(5,8), and also every noncontact-symbol(7,8) case,
under the full compatible F1/contact-four hypothesis. The first open
reduced height-eight pairs are now(6,8),(7,8). The latter can occur
only in the retained ordinary-degree72,108 contact alternative.

The proof combines [PF2RULED11](BOUNDARY_LEADING_PERIOD.md) with an
intrinsic first-Newton-edge bound and a generic-fibre pole budget.
Repeated face roots and arbitrary regular lift terms are retained.
The contact base uses only the pinned2022 degree alternative; no gcd16
shortcut or conditional2026 degree upgrade is used. Neither remaining
pair, the higher-pole problem nor actual compatible embedding is settled.

## Complete initial reduction

Let m=5 or7 and n=8. The leading differential equation gives symbols
tau^m,tau^8, with tau a nonzero cubic in L=O_B(B)=O_P1(3), after
nonzero target scalings. The leading boundary period restricts tau
to a triple root, or a mixed divisor containing the contact point q.
If tau vanishes at q, ordinary degrees are at most64,112 for m=5,
and96,112 for m=7. The first case satisfies neither the2022 alternative
of maximum degree>=125 nor its72,108 exception, because64<72.
For m=7 only(72,108), in this order, remains. We do not exclude it.

It remains to exclude the noncontact branch for each m. The period
then forces tau to have a sole triple root p distinct from q. No
position relative to the marked point is discarded: p may be the
marked point or the point at infinity in an initial chart. Choose
regular local coordinates w,z at p, with B=(z=0), so tau has order
three and the prescribed volume is

 Omega=U*z^-3 dz wedge dw, U a unit at p.                 (1)

All arguments below are intrinsic at this away-contact point.

## A complete first-Newton-edge bound

For generic c0 put F=z^m(P-c0). Its expansion has the form

 F=w^(3m)*unit+sum_{i=1}^{m-1} A_i(w)*z^i+z^m*C(w,z,c0).

Let h_i=ord_w(A_i), allowing infinity, and consider the largest first
Newton slope

 lambda=max_i (3m-h_i)/i.                               (2)

All coefficients are regular. Terms of z-degree at least m have slope
at most three. Thus a value lambda>3 is attained at an index i<m.
We prove

 lambda<=3+5/(m+6).                                     (3)

If instead lambda=b/a exceeds this bound, take positive coprime
weights a,b, and c=b-3a>0. Choosing the largest slope in(2) retains
the w^(3m) term and puts every other term above or on the same face.
The weighted initial of P is

 P_in=t^(-m*c)*H(Z), H=(A+sum_{i=1}^{m-1} B_i Z^i)/Z^m,

under w=t^a,z=t^b Z, where A!=0 and the numerator is nonconstant.
This uses the full Newton diagram; selecting a face while allowing a
lower term to dominate would not justify the argument.

The Q symbol tau^8 gives its weighted pole K>=8c. Its initial G(Z)
is a polynomial in Z^-1 with pole at most eight, by Q's actual normal
pole and regular local coefficients. The prescribed volume has value
a-2b=-5a-2c. Since (m+6)c>5a, this is greater than -(m+8)c.
The initial wedge must vanish, with ODE

 m*c*H*G'-K*H'*G=0.

Its pole at Z=0 forces K=8c, hence G^m=constant*H^8. Since gcd(m,8)=1,
each zero of H must have multiplicity divisible by m. A nonconstant
numerator of degree at most m-1 cannot satisfy that, even when its
roots coalesce. This proves(3), without a global grading assumption.

The complete integer order bounds and possible slopes strictly above
three are:

| m | bounds on h_1,...,h_(m-1) | possible lambda>3 after(3) |
| --- | --- | --- |
| 5 | 12,9,5,2 | 10/3,13/4 |
| 7 | 18,15,11,8,5,1 | 10/3,13/4,16/5,19/6 |

Each bound is ceil(3m-(3+5/(m+6))*i). Enumerating the integers between
that bound and3(m-i)-1 gives exactly the slope list. This is a finite
exhaustion of the possible first-edge slopes, not a numerical search
for functions. Regular higher terms cannot give another slope>3.

## Every surviving high edge exceeds the quartic pole budget

Let lambda=b/a in lowest terms be one of the listed slopes. Its face
has a nonzero constant and at least one other term, so it has a
nonzero complex root. On some normalized branch over that root,

 ord(w)=a*e, ord(z)=b*e, e>=1.

Since every F term has weighted value at least3ma, differentiation
gives

 ord(F_z)>=(3ma-b)*e.                                   (4)

A simple face root gives equality with e=1 in the binomial cases;
a repeated root may require further ramification and raises the
F_z order. It never lowers the bound. In particular the possible
coalescent m=7 face A+B*Z^3+C*Z^6 at slope10/3 is included.
Newton-Puiseux normalization provides branches for repeated roots
as well as simple ones. On a generic reduced P fibre F_z is nonzero
on such a branch, whose w projection is nonconstant.

The fibre differential from(1) is

 dQ/dt=unit*z^(m-3)*w'(t)/F_z.                           (5)

It follows from(4) that Q has a pole of order at least

 e*L, L=(3m-1)*a-(m-2)*b.                               (6)

Indeed the derivative order is at most -e*L-1. The complete table is

| m | lambda | a,b | F_z lower bound/e | Q pole lower bound/e |
| --- | --- | --- | --- | --- |
| 5 | 10/3 | 3,10 | 35 | 12 |
| 5 | 13/4 | 4,13 | 47 | 17 |
| 7 | 10/3 | 3,10 | 53 | 10 |
| 7 | 13/4 | 4,13 | 71 | 15 |
| 7 | 16/5 | 5,16 | 89 | 20 |
| 7 | 19/6 | 6,19 | 107 | 25 |

Every entry exceeds degree four at a single missing point. A generic
P fibre of a finite quartic map is finite of total degree four over
the Q-line; after normalization and completion every missing point
is a pole of Q and all pole orders sum to four, over all components
if necessary. This contradiction excludes every lambda>3.

Thus, for both m=5 and m=7, every actual coefficient satisfies

 h_i>=3(m-i), 1<=i<m.                                   (7)

For m=5 the two binomial controls give exact poles12 and17. For m=7
only lower bounds are needed; no face-root simplicity or genericity
of the face polynomial is assumed.

## The remaining coefficient cone makes P a polynomial in one lift

Lift tau to T in V1. We prove by descending pole degree that
P is an actual polynomial in T. Initially P-T^m lies in V_(m-1).
Suppose

 P=sum_{j=l+1}^m a_j*T^j+R_l, R_l in V_l.

The coefficient of normal pole l in each displayed T^j has order at
least3l at p: write T=tau/z+regular and expand; every contribution
at that normal pole contains at least l factors of tau. Equation(7)
gives the same order for the coefficient of P. Consequently the
symbol of R_l has order at least3l. That symbol is a section of L^l,
of total degree3l, so it is a scalar multiple of tau^l. Subtract that
multiple of T^l and continue. At l=0 the remainder is constant on the
projective completion. Therefore

 P=f(T), deg(f)=m.

The [polar-centralizer lemma](LOW_COORDINATE_POLE.md) now excludes any
independent companion with the stated nonzero volume. At each positive
pole j of Q, the leading equation forces its symbol to be a scalar
times tau^j. Subtracting the corresponding actual power of T preserves
dP wedge dQ because P is a polynomial in T, and strictly lowers j.
The process ends at a constant and contradicts the nonzero volume.
These auxiliary subtractions are not claimed to be target automorphisms.

This covers every lower coefficient and all regular lifts, including
the boundary marked point. Contact cases were handled before applying
the unit-volume calculation. Hence(5,8) is excluded in full, and the
noncontact-symbol(7,8) case is also excluded in full.

## Current remaining gate

The first open reduced height-eight pairs are(6,8),(7,8). For(6,8),
the PF2RULED11 table still allows the noncontact common sextic with
one root and ordinary degrees(96,128), or its retained contact72,108
alternatives. For(7,8), only its contact72,108 alternative remains.
Neither remaining pair is closed in full. The pole-four triple and
contact-mixed layers, arbitrary higher poles and actual ruled
compatibility remain open.

## A source-version shortcut that was not adopted

A secondary abstract for arXiv:1708.09367v1 advertises a gcd lower
bound of25. The
[revised primary version v2](https://arxiv.org/html/1708.09367v2)
explains that the available minor intersection formula is an inequality
and that the desired infinite-family exclusions do not follow. The
[retained source review](sources/gcd-boundary-review.json) pins that
version boundary. No gcd>=25 theorem is imported from the superseded
abstract, and gcd>=16 alone does not exclude equality16. No global
gcd16 exclusion is adopted here. The proof above uses neither shortcut.

## Reproducible evidence and assurance

The [checker](verify_pole58.py) controls the weighted ODE, non-fifth-power
multiplicities, complete integer first-edge bounds, exact depression,
the complete slope tables for both m5 and m7, supplemental symbolic
Newton controls, coefficient-cone induction data and pole budgets. Its [receipt](pole58_receipt.json) separates
these exact controls from the written all-family Newton coverage,
section-degree argument, local normalization and finiteness proof.
It also retains the primary source-version guard without adopting it
as a theorem. No formal verification, independent full reproof,
external review or novelty priority is claimed.

The goal remains active: this closes one pair, not the remaining
finite-map realizability problem, actual ruled compatibility, either
full quartic packet or JC2.

The subsequent [PF2RULED14 global cusp-jet theorem](GLOBAL_CUSP_LEADING_OBSTRUCTION.md)
excludes every one-root noncontact symbol, including all noncontact(6,8)
profiles. The current height-eight pairs(6,8),(7,8) are both contact
cases with ordinary source degrees72,108; the older tables retain their
source-specific necessary inputs.

The subsequent [PF2RULED16 full height-eight exclusion](POLE_HEIGHT_EIGHT_EXCLUSION.md)
closes both remaining pairs(6,8),(7,8), including every contact stratum.
The earlier height-eight handoffs in this note are historical. The current
minimum height is9; first reduced pairs are(5,9),(6,9),(7,9),(8,9).
