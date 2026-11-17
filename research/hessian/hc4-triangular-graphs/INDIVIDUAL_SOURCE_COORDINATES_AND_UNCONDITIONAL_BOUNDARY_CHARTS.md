# Individual source coordinates and unconditional boundary charts

Every added unramified prime of the actual canonical plane Keller
completion has an individual polynomial source coordinate regular on
its whole component open. Consequently its complete algebraic boundary
is a smooth affine line and the whole normalized plane root charts
exist. The coordinate can depend on the component. Simultaneous
alignment is not a conclusion here.

The proof combines an elementary intersection-sign argument with a
specific part of Borisov's minimal-coordinate valuation analysis.
It does not infer a coordinate from negative definiteness alone,
does not adopt Xie's full valuation-intersection theorem, and does not
use a bound on any polynomial degree.

## 1. Valuation criterion and actual consequence

Let R=C[x,y], L=C(x,y), and let v be a primitive divisorial valuation
of L centered outside the marked affine plane X=Spec R. Normalize its
value group to Z. Suppose
```
ord_v(dx wedge dy)=0,
trdeg_C(R intersect O_v)=2.                              (1)
```
Then there are polynomial source coordinates (a,b) of R and lambda in C
such that
```
v(b-lambda)>0.                                          (2)
```
The coordinates are obtained after a polynomial automorphism of the
GIVEN source; they are not Darboux coordinates of an abstract rational
field.

For an actual complex polynomial Keller map K:X->Y, let Z be the full
finite normalization and U its entire canonical etale locus. The
[canonical boundary theorem](CANONICAL_ETALE_COMPLETION_BOUNDARY_TOPOLOGY_AND_FREE_DECK_ACTION.md)
gives smooth affine component opens
```
M_i=X union D_i subset Z,
```
where the whole D_i is present, its normalization is A1, and it is
homeomorphic to C. The actual unit volume on U gives the first clause
of(1). The inclusion
```
O(U) subset R intersect O_(v_i)
```
and Frac O(U)=L give the second clause. Each v_i is centered at infinity:
otherwise its center in the separated M_i would also lie in the open X,
contradicting uniqueness of the valuation center.

Thus(2) applies to every D_i. The coordinate b_i is regular on all of M_i.
Indeed at source height-one points it is a polynomial, and at the only
remaining prime divisor D_i its value is nonnegative by(2).
Normality of M_i then extends it across every closed point.

The [conditional whole-chart theorem](COMMON_COORDINATE_FORCES_SMOOTH_BOUNDARY_AND_WHOLE_ROOT_CHARTS.md)
now applies without an additional coordinate assumption for THIS
component. Hence every original added unramified D_i is algebraically
A1 and smooth, and its whole normalized plane charts have
```
b_i=lambda_i+z^r_i,
a_i=z^(1-r_i)*t+h_i(z),   h_i finite Laurent.
```
The actual source volume has a global regular primitive on M_i.
After a polynomial t-shear the principal part is nonzero, with no
z^-r_i coefficient and with largest pole p_i>=r_i+1. These conclusions
are inherited from that proved whole-chart theorem, rather than
assumed in proving(2).

This settles the INDIVIDUAL coordinate and component-chart obligation.
It does not prove that b_i and b_j are affine functions of the same
source coordinate when i!=j. No common-pencil, global LND or inverse
conclusion is asserted.

## 2. Determinant convention

Extract v by point blowups at infinity from some polynomial-coordinate
P2 compactification of X, and write
```
E_0,...,E_N
```
for the whole boundary, including the strict transform of the initial
line at infinity. All these classes form a basis of the Picard group.
Let I=(E_i dot E_j) be their intersection matrix. Its signature is (1,N),
and det I=(-1)^N, since the total-transform basis has diagonal
(1,-1,...,-1).

The divisor dual to E_v is the boundary-supported class C_v with
```
C_v dot E_j = Kronecker_delta_(v,j).
```
Define
```
d(v)=det(-I with row and column v deleted).
```
The cofactor identity gives exactly
```
C_v^2=(I^-1)_(v,v)=d(v).                                (3)
```
This is the unnormalized determinant convention, not the determinant
of a normalized valuation matrix.

Let kappa(v) be its augmented canonical label. Since the rational form
dx wedge dy represents the canonical divisor and the boundary is
reduced,
```
kappa(v)=ord_v(dx wedge dy)+1=1.                         (4)
```
These labels are unchanged by further extractions and by polynomial
source automorphisms. A source automorphism changes the volume form
only by a nonzero constant.

## 3. Richness forces d(v)<0: a direct surface proof

Choose algebraically independent p,q in R intersect O_v. The source
polynomial map (p,q):X->A2 is dominant and generically finite.
Resolve its projective rational extension and extract v on a common
smooth surface T obtained from the source P2 by point blowups at infinity.
This is possible without modifying X, since p,q are regular there.
Write
```
phi:T->P2_target,
H=phi^*(target line at infinity).
```
The pullback satisfies
```
H^2=[L:C(p,q)]>0.                                       (5)
```
It is supported on the boundary of X. It has coefficient zero on E_v,
because p and q have nonnegative value at v and its generic image lies
in the affine target chart.

Let V be the real span of ALL boundary classes except E_v.
Then H belongs to V and is positive by(5). The whole boundary
intersection form has signature (1,N). A subspace containing a positive
vector is nondegenerate with exactly one positive direction: split off
that vector; its orthogonal complement in the whole space is negative
definite. Therefore V has signature (1,N-1), so its minus-intersection
determinant is strictly negative. This proves
```
d(v)<0.                                                 (6)
```
No converse criterion is needed. In particular this proves the needed
direction of the singleton valuation-intersection statement directly;
it is not an adoption of Xie's theorem from its displayed statement.

The same argument can place several zero-coefficient boundary divisors
in the orthogonal complement of one positive H, but no simultaneous
coordinate conclusion is inferred from it in this note.

## 4. The precise minimal-coordinate input from Borisov

We use the following consequence of the proof of Theorem4.3 in
Borisov's primary author version of *On two invariants of divisorial
valuations at infinity*, printed pp12--15, especially Lemmas4.3--4.4.
Choose a polynomial-coordinate P2 compactification minimizing the
finite ancestor count of v. Its Case2, where the second blowup is at
the intersection of the first exceptional curve with the line at
infinity, has no final label
```
kappa=1, d<0.
```
In that case the first positive-label vertex has label1 and determinant
c>=0. Further vertices retaining label1 arise only by blowing the
adjacent 0--1 edge, with determinant
```
(m+1)^2*c+m^2+m >=0,   m>=0.
```
Once a label at least2 is reached, later ancestors cannot return to1.
The alternative free centers in the preceding negative-label chain
contradict the chosen minimal compactification, by the source
automorphisms in that proof. Thus (4),(6) force its Case1:
the first TWO centers are free.

This is a specialized imported proof consequence, not an inference from
the theorem's bounded-family conclusion alone. The primary minimality
argument, diagrams and non-positive/positive transitions were inspected.
The exact source locators and adopted classical blowdown foundations
are recorded in Section7. No effective bound on ancestors is needed.

## 5. Two free infinity centers give an actual polynomial coordinate

We now translate that graph conclusion into the needed regular source
coordinate. Let (a,b) be affine coordinates for the minimizing source
P2, with the first center at the a-infinity point [1:0:0].
Near that point take
```
u=1/a,   w=b/a.
```
After the first blowup, away from the strict line at infinity the chart
is
```
u,   b=w/u.
```
The first exceptional curve is u=0, with b its finite coordinate.
A second FREE center is therefore
```
u=0, b=lambda
```
for some lambda in C. The function b-lambda is regular in a neighbourhood
of that center and vanishes there. On the second exceptional curve its
value is positive; for example the free chart is
```
a=s^-1,   b=lambda+s*t.
```
It is regular on a neighbourhood of the ENTIRE second exceptional curve,
including its intersection with the first exceptional curve.

Every later point in the minimal extraction of v lies above this second
center. Pullback preserves regularity of b-lambda there and its positive
exceptional order. Consequently
```
v(b-lambda)>0,
```
including when later centers are intersections of boundary components.
This proves(2) in the actual polynomial source ring.

The normality argument in Section1 then gives global regularity on
M_i, with no assertion that an arbitrary rational function regular at
one boundary prime is already regular on the whole canonical U.

## 6. Individual coordinates are not a common coordinate theorem

The distinction is necessary even with zero volume orders. The retained
[contracted-axes control](VALUATION_MATRIX_AND_CONTRACTED_SOURCE_AXES_CONTROL.md)
has valuations
```
v1(x)=-1, v1(y)=2,
v2(x)=2,  v2(y)=-1,
```
both with volume order0, negative determinants, and individual regular
source coordinates y and x. Their simultaneous nonnegative ring
C[xy,x^2*y,x*y^2] contains no nonconstant source coordinate.
That control contracts source axes and is not an actual Keller source.
It does not contradict the present individual theorem.

For actual components, proving simultaneous alignment must use further
source geometry beyond the individual criterion. This note neither
assumes nor supplies such an alignment. It also does not identify a
regular source coordinate with an integral target coordinate, descend
the source translation to the target, or prove a Keller inverse.

## 7. Primary inspection and exact verification scope

The essential primary is Alexander Borisov,
*On two invariants of divisorial valuations at infinity*, author PDF,
metadata13 February2013; published J. Algebraic Combinatorics39(2014),
691--710, DOI10.1007/s10801-013-0462-9.
The retained exact cache is
[sources/borisov-divisorial-valuations-2013.pdf](sources/borisov-divisorial-valuations-2013.pdf),
from the
[author's primary URL](https://people.math.binghamton.edu/borisov/documents/papers/divisorialvaluations.pdf).

The complete Theorem4.3 proof on printed pp12--15 was read; those are
PDF zero-index pages11--14. All four pages were also rendered and
visually inspected, including the Case1/Case2 diagrams, the minimality
moves, the first-positive-label discussion and Lemma4.4's determinant
formula. Definitions and Proposition2.1 on printed pp3--6, and the
determinant recurrence Lemma3.7 on pp9--10, were checked for conventions.
The earlier end of Theorem4.2, printed pp11--12, was read for the
replacement-P2 blowdown used by minimality. The usual contraction
and smooth rational-surface facts in that primary argument are adopted
classical inputs; no independent reproof of every such fact is claimed.

Xie's primary definition and Theorem1.2 were read for comparison,
but its full potential-theoretic proof is NOT adopted as a dependency.
Section3 replaces the exact implication needed here by a surface proof.
Resolution of a rational plane map by point blowups, the Picard
intersection form of a blown-up P2 and the projection formula in(5)
are the classical surface foundations used there.

The checker gives exact cofactor/sign-convention controls, blowup
determinant and canonical-label recurrences, the nonnegative label1
formula, local free-center coordinate identities, and controls of the
individual-versus-simultaneous distinction. These finite controls do not
prove the unrestricted minimality classification or replace its primary
proof. The source audit and receipt pin the primary and all retained
proof inputs. Internal checking is not formal verification, external
human review or an HC4/JC2 solution.
