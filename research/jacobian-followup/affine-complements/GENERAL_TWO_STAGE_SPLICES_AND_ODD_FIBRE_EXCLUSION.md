# General two-stage splices and the odd-fibre exclusion

15 November 2026. PF2RULED80 removes the literal-binomial hypothesis
from the ZERO-normal-degree two-stage exclusion and closes an entire
new negative-normal-degree sector: BOTH last cabling numbers odd.
Every remaining genuine two-stage quartic packet must have

    eta=-2h<0, h>=1,
    k=r*(pq-p-q)+1-2h>=2, gcd(r,k)=1, rk EVEN.             (1)

Here p,q>=2 are coprime and (r,k) are the ACTUAL ordinary winding
and cabling coefficient of its reduced infinity splice. Arbitrary
lower-weight terms in the key are allowed. No global coordinate
normalization of those terms, or equivalence of affine complements,
is assumed. The remaining negative-even sector, longer target
sequences, higher degrees and the full goal stay OPEN.

The source criteria remain degree-free under PF77-79's explicit full
smooth/exact/whole-line hypotheses. The application is to PF46/50's
actual quartic packet with exactly one fold R and one unramified D.
There is no arbitrary higher-degree packet classification.

## 1. Actual splice data already determines the peripheral model

Work with an irreducible polynomially parametrized curve with a genuine
reduced TWO-stage infinity splice, first torus pair(p,q), p,q>=2,
and final pair(k,r), k,r>=2. The coprimalities are
gcd(p,q)=gcd(k,r)=1. Neumann--Norbury Section4, Figure1, gives the
positive ordinary cable

    K=C_(r,k)(T(p,q)), k<rpqr.                            (2)

The virtual companion-to-final-knot linking number is k. Hence this
is the integer coefficient relative to the ZERO-linking preferred
longitude of the companion, not a framing chosen from a key's name.
The semigroup condition k in <p,q> is necessary for actual plane-curve
infinity data. It is an additional unrestricted integer condition,
not replaced by a modular coprimality test.

Theorem5.1 on complete PDF9-11 distinguishes actual key polynomials
with lower-weight corrections from literal binomials. Those corrections
can change the final integer k relative to the uncorrected binomial
defect. Once the ACTUAL splice is known, they do not add another
parameter to the oriented knot or its preferred peripheral pair.
No connectedness theorem about polynomial moduli is needed here.

On the efficient transverse completion used in PF59, Section5.2's
canonical recurrence (complete PDF12-13) gives

    d1=p+q, Delta2=k-pqr,
    d2=r*d1+Delta2=k-r*(pq-p-q),
    eta=d2-1=k-r*c-1, c=pq-p-q.                          (3)

This computation depends on the splice, not on a Wronskian identity
for X^q-Y^p. In the actual quartic, PF59 gives eta<=0.

If eta=0, (2)-(3) give EXACTLY the same marked cable model and integer
parameters as PF63/79. The earlier argument after that model is fixed
uses the full source charts, the marked branched end and its actual
sheet fillings; it never uses coefficients of X,Y or the binomial key.
An oriented ambient knot isotopy preserves positive meridian and
zero-linking longitude. It transports every representation and every
true meridian-cycle filling, including the fixed sheets and R mark.
Thus PF79's end obstruction applies to EVERY general two-stage eta0
packet as well. This is transfer of the MARKED END MODEL, not a
global polynomial automorphism or affine-cover construction.
Consequently eta=0 is now closed without a literal-binomial hypothesis.

For specificity, the dependency audit is: PF60's binomial Wronskian
was used to FORCE eta0; here eta0 is the case assumption. PF63's
finite equations then use only the knot peripherals. The companion
covers, their Seifert models and all meridian fillings in PF65-76 are
covers of those standard knot pieces, independent of the affine key.
PF77-79 compare them with PF50's whole source charts. The actual packet
hypotheses providing those charts and the S4 action remain in force.
None of these steps requires producing a binomial representative of
the general affine curve, or matching its finite singularities.

## 2. General cable peripherals, with k independent of r*c+1

For arbitrary coprime p,q,r,k in (2), choose integers

    q*u0+p*v0=1, r*u-k*v=1.

The companion and cable groups with their actual peripherals are

    <a,b,z | a^p=b^q, z^r=M^k*L^r>,
    M=a^u0*b^v0, L=a^p*M^(-pq), H=a^p,
    F0=M^k*L^r=z^r,
    mu=M^u*L^v*z^(-v), lambda=mu^(-rk)*F0.              (4)

Clay--Watson Section2, complete PDF6, supplies the root amalgam and
these meridian/longitude words. Their winding p and coefficient q
are our r,k; their root t is z. The positive Bezout solution can be
replaced by any integral solution: u->u+kt,v->v+rt changes mu by
F0^t*z^(-rt)=1, using that F0 commutes with the boundary torus and z.

Let tau be the image of mu in the actual quartic S4 monodromy.
PF60 supplies tau a transposition and lambda=1. Therefore

    image(F0)=tau^(rk).                                  (5)

When rk is EVEN the cable fibre image is identity, as in PF63.
When rk is ODD it is a TRANSPOSITION. Deleting that second case
would incorrectly generalize the old eta0 enumeration.

## 3. An exact finite-period problem for every two-stage key

All relations and words in (4)-(5) can be evaluated in S4 using only
(p,q,r,k) modulo12. For both Bezout pairs use ANY modular solution.
Actual coprimality lets a modular solution lift to an exact one after
adding multiples of12: for example r*u-k*v=1+12T can be corrected
by choosing r*m-k*n=-T. Two modular solutions differ by (k*t,r*t)
modulo12, and (4) makes their images identical. The companion pair
has PF63's analogous independence argument.

Thus the general, necessary and sufficient PERIPHERAL S4 test is
the following finite problem, with no degree cutoff:

    a^p=b^q, z^r=M^k*L^r,
    mu=M^u*L^v*z^(-v) a transposition,
    mu^(-rk)*M^k*L^r=1, <a,b,z>=S4.                     (6)

There are96 admissible ordered(p,q) residue pairs and96 admissible
ordered(r,k) pairs, hence9216 cells. These encode only NECESSARY
modular coprimality. Actual integer coprimality, k<rpqr, the semigroup
condition and eta inequalities must still be checked for actual keys.

The complete certificate has1344 nonempty cells,66048 representations
and2752 simultaneous-conjugacy classes, each orbit of size24.
Exactly256 classes have rk odd. The other2496 have rk even;
448 of those have a cycle in the lifted-piece graph and128 tree
classes have positive companion base genus. PF62/69 exclude them
for the actual source at every degree under their stated hypotheses.
The remaining1920 even-product planar tree TYPES are a necessary
finite frontier, not realizations or a count of polynomial curves.

All1152 old eta0 cells, including their596 classes and tree counts,
are recovered exactly by substituting k=r*c+1. This comparison uses
the correct partial root conjugation between the old meridian convention
and (4); it does not identify unrelated generator tuples.

## 4. Odd product forces a single companion sector

Assume r,k odd. Equation(5) gives z^r=tau. In S4 an odd power can
be a transposition only when z IS that transposition. Hence z=tau.
Put beta=M^u*L^v, so mu=beta*z^(-v) gives beta=tau^(v+1).
The boundary-basis inverse of (4) is

    M=F0^(-v)*beta^r, L=F0^u*beta^(-k).

Since r,k odd and r*u-k*v=1, these identities give

    M=tau, L=1, <a,b>=S4, H=a^p=b^q=1.                 (7)

Thus pq must be even. Write e for its even exponent and o for its
odd exponent. The odd generator is a 3-cycle. A double transposition
on the even side generates at most A4; a transposition joining the
fixed sheet of a 3-cycle gives a 4-cycle meridian, contrary to (7).
Therefore the even generator is a FOUR-cycle, and

    e=4B, o=3t, gcd(4B,3t)=1.                            (8)

This is a transitive S4 companion and C2 cable image, but its generic
cable fibre has degree TWO rather than PF72's even-product degree1.
The complete odd part of (6) verifies precisely8 companion residue
types times32 admissible odd(r,k) pairs:256 classes and6144
representations, all with tree incidence and planar bases.

## 5. The odd lifted fibre, actual fillings and zero R mark

Retain all fixed sheets. The main two-sheet cable piece has generic
fibre degree2 and base degree1. Its old cone has order r, still odd.
The branching along the ordinary target fibre contributes an exceptional
order2. Here is its local framing check. On the ramified boundary,
the kernel torus has basis(mu^2,F0*mu^-1); the primitive lifted regular
fibre is

    F=F0^2=mu^2+2*(F0*mu^-1) in additive notation.

Filling mu^2 therefore makes R an order2 fibre, with

    Delta(meridian_R,F)=2.                               (9)

The filled main cable P is an ACTUAL disk Seifert piece with cones(r,2),
and drilling R gives the annulus piece with its ONE genuine cone r.
Its H1 is free rank1: the two cone orders are coprime. On its outer
boundary use m=M^2,L. The branched map P to the ACTUAL companion
solid torus V sends m to0 and L to a generator. Hence its free H1 map
is onto and, since it has no torsion, is an isomorphism. Consequently

    H1(P)=Z*L, [m]=0, [F]=2r*L, [K_R]=r*L.              (10)

The last equality uses F=2K_R and torsion-freeness, not division in
an arbitrary finite group. Both other cable pieces have degree1 and
give their actual fixed-sheet meridian fillings.

After those caps, the companion E has base disk with cones(B,t,3t),
with orders1 omitted. Its regular fibre H is H=(eo/2)m+L. PF72's
closed companion computation W=E(m), independent of the cable, gives

    H1(W)=Z/t, [H]=0.

Explicitly its fibre degree is1, base degree4 and Euler number-4/(eo);
lcm of cone orders is eo/4, so PF72's exact fibre-order formula gives
order1. The cone presentation gives the displayed cyclic group.
Gluing (10) simply performs that m filling, with NO extra cable torsion.
It follows, integrally, that the full target end satisfies

    H1(end)=Z/t, [K_R]=r[H]=0.                           (11)

PF67's common-source nonzero R self-linking excludes EVERY common
source in this odd sector. Its hypotheses concern the actual source
charts and do not require a literal binomial target key.

## 6. Every nonsolid odd-sector companion is excluded

If E is nonsolid, its disk(B,t,3t) has at least two genuine cones,
negative orbifold Euler number and incompressible boundary. The root
annulus-cone-r piece also does. In the actual basis(m,L), their fibres
are H=(eo/2)m+L and F=k*m+2r*L, with intersection

    Delta(H,F)=r*eo-k>0.                                (12)

Thus there is exactly ONE canonical JSJ edge. PF79's terminal/Klein
and full peripheral-recognition argument applies. The source root
cannot be rho=1 (distance3), nor rho>=3 (regular R, distance1), since
the target distance is2 by (9). A corner merger has two root cones,
whereas the target has one. Only the rho=2 collar is possible.

Write its initial D data p_D=Q*a,sigma=l*a, Q=r, 1<=l<Q,
gcd(l,Q)=1. Its marked distance gives

    epsilon=Q-3l=+2 OR -2.                              (13)

If B,t>1, the target terminal disk has THREE genuine cones, violating
PF79's two-cone budget immediately. Only two nonsolid families remain.
There must be exactly two D nodes; its terminal orders are(a,kappa).
If l=1 its initial contact is integral, otherwise the first gcd drops
from l*a to a. In BOTH cases, with last drop j_*=l*a-1+kappa,
all supports allowed by these data have

    S_D=(l-1)*a*(Q*a)+(a-1)*[(2Q-l)*a+1-kappa],
    n_D=[(l+1)Q-l]*a^2+(1-kappa)*a+kappa,
    det Q_src=3n_D-(Q*a)^2
       =C*a^2+3a+3kappa*(1-a),
    C=(2-epsilon)Q+epsilon in {2,4Q-2}.                  (14)

Higher noncharacteristic coefficients and terms below the final gcd
drop cannot alter these contact sums. A failed gcd or support inequality
already removes a proposed datum; none is assumed realizable.

For B=1,t>=3, the target terminal orders are(t,3t), and the target
closed end has order t. If(a,kappa)=(t,3t), (14) gives

    epsilon=2: det=t(12-7t)<0;
    epsilon=-2: det=t[(4Q-11)t+12]>t, Q>=7.

If(a,kappa)=(3t,t), it gives respectively

    det=3t(3t+4)>t;
    det=3t[(12Q-9)t+4]>t.

These violate source negative definiteness or (11).

For t=1,B>=2, the target terminal orders are(B,3) and its end has
order1. Since 3 does not divide B, and C=2 mod3 for BOTH signs in
(13), formula(14) gives det=2 mod3 when(a,kappa)=(B,3), or det=0
mod3 when(a,kappa)=(3,B). Neither can be1. This excludes EVERY
nonsolid-companion complementary source in the odd sector.

## 7. The sole solid companion, with both cap signs retained

Only B=t=1 remains, first pair(e,o)=(4,3). E is an ACTUAL solid
torus with cone3. As in PF79, its H=6m+L, H1(E(m))=0 and H=3v
in a solid-torus generator v give [m]=+/-v. Its meridian is therefore

    eta_E=3m+L OR 9m+L.                                  (15)

Both signs are kept; no orientation or section twist is guessed.
Capping the root annulus with fibre F=k*m+2rL gives the SINGLE target
disk with cones(r,J), where

    J=6r-k OR 18r-k.                                    (16)

Indeed PF59 and (3) give k<=5r+1; oddness and gcd(r,k)=1 strengthen
this to k<=5r-2. Thus both J are positive and J>=r+2. This target
has negative orbifold Euler number, centreZ, and R-meridian distance2.
PF79's whole-Seifert source classification forces the rho=2 single D
disk: the alternative full merger has THREE genuine cones.
Its orders require

    p_D=J, kappa=r, sigma=J+1-r,
    |p_D-3sigma|=|3r-2J-3|=2.

Hence J=(3r-1)/2 or(3r-5)/2. The exact source determinant is

    3n_D-p_D^2=2J^2+3(1-r)J+3r
               =6r-1 OR 5,                             (17)

neither of which is1, the target order in (11). If the proposed
characteristic gcd fails it was impossible already. This closes both
solid-cap signs and completes the ENTIRE odd-product two-stage
quartic sector in BOTH source pencils, at all target/source orders.

The target exteriors used in Sections6-7 are irreducible,
boundary-irreducible and aspherical. Their genuine edge, or an essential
vertical annulus in the single disk piece, gives sufficient largeness.
The full peripheral Waldhausen realization used in PF77-79 therefore
preserves the R meridian and its distances. No homology-only recognition
or deletion of D from the full finite map is used.

## 8. Exact remaining frontier

c=pq-p-q is odd. In (3), odd rk gives odd eta and even rk gives even
eta, by gcd(r,k)=1. PF59 excludes eta>0, Section1 excludes eta=0
for general two-stage keys, and Sections4-7 exclude odd rk. Thus (1)
is necessary for EVERY surviving genuine two-stage quartic packet.

The remaining finite periodic test can equivalently use
(p,q,r) mod12 and h mod6, substituting k=r*c+1-2h. It has6144
admissible cells before actual inequalities and source gates; its1920
planar tree representation TYPES are retained in the complete census.
They are not declared affine extensions or all excluded. The eta0
closures must not be promoted to negative eta merely because some
residue data coincide. Actual k and eta remain unrestricted integers.

The [checker](verify_general_two_stage.py),
[receipt](general_two_stage_receipt.json) and
[manifest](sources/general-two-stage/manifest.json) seal the unbounded
peripheral reduction, complete9216-cell certificate, old-slice comparison,
odd classification, generic fibre-orbit genera, exact marked contact and
cap identities and ordinary full branched-cover controls. New primary
audits inspected complete Neumann--Norbury PDF6/7/9/10/11/12/13 and
Clay--Watson PDF6. No global-moduli connectedness, Seifert surgery
orderability result or polynomial key normalization is imported.

The written splice-transfer, local odd-fibre, contact and cap arguments
prove the universal claims, together with explicitly imported topology.
The finite group certificate is exhaustive AFTER the written period
proof; the bounded polynomial/source controls are not a universal proof.
No full two-stage/quartic/higher-degree exclusion, finite Keller map,
formal verification, independent full reproof, external review or
novelty/priority claim is asserted. The full goal remains active.
