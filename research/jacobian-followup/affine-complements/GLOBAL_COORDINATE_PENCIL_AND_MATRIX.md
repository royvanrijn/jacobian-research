# The unramified primes share one coordinate pencil

15 November 2026. PF2RULED50 gives an all-degree compatibility theorem,
using Xie's intersection-of-valuation-rings theorem instead of further
principal-part searches.

**For every actual plane Keller finite normalization, all unramified
omitted primes have the SAME intrinsic regular source-coordinate pencil.**
For every index-two omitted prime and every unramified omitted prime,
their intrinsic coordinate pencils are either the same or complementary:
there is no third, noncomplementary coordinate configuration. Other primes
and arbitrary finite field degrees are allowed and retained.
The common unramified coordinate can still have poles on ramified primes;
this is not an assertion that a source coordinate is regular on all of N.

In the remaining quartic packet this removes the noncomplementary-coordinate
gap of PF2RULED48/49. The common-ruling and remaining complementary sectors
are still open; the result does not exclude the whole packet.

## 1. A necessary matrix condition, not an affineness certificate

Write U=A2 and let R,D be any two omitted prime valuations of a normal
affine N containing U. Every global regular function of N is a polynomial
on U and has nonnegative order at R,D. Therefore

    O(N) subset B_(R,D)=C[U] intersect O_(v_R) intersect O_(v_D),
    trdeg_C Frac(B_(R,D))=2.                          (1)

No affineness or finiteness assertion is made after deleting other primes.

Choose any common smooth projective completion exposing both valuations.
For every boundary prime E let checkE be the divisor supported on the
whole boundary with checkE.E'=delta_(E,E'). Xie, *Intersections of
valuation rings in k[x,y]*, Theorem1.2, states that the algebra in (1)
has transcendence degree two iff the matrix

    M_(R,D) = ((checkR^2,checkR.checkD),
               (checkR.checkD,checkD^2))               (2)

is negative definite. We import precisely that theorem. In particular
(1) forces negative definiteness of every pair submatrix. Passing this
test is NOT an affineness certificate: PF2RULED48's simple/simple p=q=2
pair passes the matrix test but has a complete curve and is excluded.

Dual pairings are unchanged under further boundary blowups, because the
pullback of an old dual divisor is orthogonal to every new exceptional
curve. Thus a P1 times P1 completion can be used for computing (2), with
a common refinement dominating P2 if desired. We use the UNNORMALIZED
primitive integer valuations throughout; no skewness scaling is omitted.

## 2. The intrinsic coordinate pencils and a cross pairing

For each unramified prime the smooth affine component open U union D
has the complete classification of PF2RULED46. For each index-two prime
the component open is smooth by PF2OMIT7, affine by the retained BGV
component-open theorem, and its source volume is exact. PF2RULED45
therefore supplies a regular source coordinate on that component open.
These uses are componentwise and do not delete primes from the full map.

PF2RULED38 says that a regular source coordinate at a given omitted prime
is unique up to a nonzero affine change. It consequently specifies an
intrinsic pencil. Its boundary value is constant; shifting that individual
coordinate to zero defines its positive order. Its complementary source
coordinate has an intrinsic positive pole height q or p.

Let f and g denote distinct such pencils, represented by coordinates.
Put

    m=(generic f-fibre).(generic g-fibre) in A2
      =[C(x,y):C(f,g)]>=1.                            (3)

This is an intersection number of SOURCE coordinate fibres, not the
Keller field degree. In coordinates (x,f), a coordinate g not affine
in f has a constant nonzero leading x coefficient and degree_x(g)=m.
Indeed the arbitrary-length Jung--van der Kulk alternating normal form
used in PF2RULED38 starts with (a*x+P(f),b*f+c). After a nonlinear
swap step the x degrees are distinct, and the larger degree's leading
coefficient is constant. Further nonlinear steps preserve this assertion;
the final affine combination either takes the larger degree, the smaller
degree or the original f pencil. This proves the claim for every word.
Since v_R(f)>0 and v_R(x)=-q, no lower x-degree term can cancel it:

    v_R(g)=-m*q.                                     (4)

The cross pairing is exactly

    checkR.checkD=m*q*p.                             (5)

Here is a proof retaining the actual pencils. Choose a source coordinate
h complementary to g and use P1_h times P1_g. All centres extracting R
are over g=infinity, by (4). Extract R first. A generic curve g=beta
has projective closure meeting the horizontal section E_h={h=infinity}
once, away from all these centres. Pairing its principal divisor with
checkR shows that the coefficient of E_h in checkR is -v_R(g-beta)=m*q.

The D extraction is over E_h at the finite value g|D. Its centre is
disjoint from the entire R cluster. The total pullback of E_h has
coefficient p=v_D(1/h) at D, and every other old boundary divisor has
coefficient zero there. Hence pushing checkD back to the R-extracted
model gives p times the divisor dual to E_h. Pairing with checkR gives
p times its E_h coefficient, proving (5). The D zero order, all lower
Laurent coefficients and all lengths of the relative automorphism word
are retained in this argument.

When m=1, (f,g) is an actual source automorphism pair: monicity and
degree_x(g)=1 give g=a*x+P(f), a!=0. When the pencils agree, m in
(3) is not defined and formula (5) must NOT be used. Their centres may
overlap, and their cross pairing can be negative. The actual common-ruling
controls in PF2RULED49 are preserved.

## 3. The diagonal entry from the full source chart

Use a primitive entire chart

    f=z^r, x=h(z)+z^(e-r)*t, h=sum c_j*z^-j, max j=q,
    gcd(r,support h)=1.                              (6)

For e=2,r=1 this means x=f*t+A+h(f), with A retained. For r>=2 we
use precisely the whole charts of PF2RULED43/46/48, with negative h
and e-r<=0. Define

    j(a)=max{j:c_j!=0, r does not divide a*j},
    S=sum_(a=1)^(r-1)(2q-j(a)), n=S+2q-r+e.           (7)

At r=1 the sum is empty. Then

    checkR^2=-n.                                    (8)

We give the intersection computation, rather than infer (8) from end
homology. Take two generic constants t0,t1 in (6). The irreducible full
norm curve P_t0 has x-degree r, includes x^r*f^q, and has total degree
q+r in the source
coordinates (x,f). Its germ at x=infinity,f=0 is one primitive curvette
of R; its other infinity germ has a different point on the source P2.
Its order after substituting the second whole chart is

    v_R(P_t0)=S+2q-r*q+e-r=n-r*q.                     (9)

This is the full conjugate-difference derivative identity, including
every contact j(a); the first chart-varying power is e-r. Homogenizing
P_t0 at that infinity point multiplies by x^-(q+r). Its local
intersection with the second curvette is therefore

    I_local=(n-r*q)+q*(q+r)=q^2+n.                    (10)

For completeness, the elementary blowup intersection calculation gives
checkR^2=q^2-I_local. In the orthogonal total-transform basis of a
successive point-blowup resolution, checkR has line coefficient q
and exceptional coefficients the successive curvette multiplicities.
Its square is q^2 minus the sum of their squares. The local Noether
intersection calculation for two generic curvettes is exactly that sum;
they meet different points after the last extraction. This proves (8).
The r=1,e=2 constant A changes no order or local intersection in (9),(10).

Apply the same formula to D, with order s, pole p, index1 and sum S_D:

    checkD^2=-n_D, n_D=S_D+2p-s+1.                   (11)

Combining (2),(5),(8),(11) gives the exact necessary condition for
DISTINCT pencils:

    n_R*n_D > m^2*q^2*p^2.                           (12)

The same test applies to any pair of charts (6) for which they are
established, with their own indices. It does not assert existence of a
regular coordinate or smoothness at arbitrary higher-index primes.

## 4. Exact unramified volume and index two make m small

For an unramified prime, exactness removes c_s and the entire support
begins at s+1. Thus p>=s+1 and every j(a)>=s+1. Consequently

    0<n_D<=2s*p-s^2-s+2<=p^2,                        (13)
    p^2-(2s*p-s^2-s+2)=(p-s)^2+s-2>=0.

The final inequality follows from p>=s+1, including s=1. Equality is
possible at s=1,p=2. For TWO unramified primes with distinct pencils,
(12) and (13) would give

    p1^2*p2^2>=n1*n2>m^2*p1^2*p2^2>=p1^2*p2^2,

a contradiction. Hence ALL unramified omitted primes share one pencil.
Their boundary constants need not be equal; this is a pencil theorem,
not a claim that all primes belong to the same special fibre.

For index-two R with r=1, exactness gives q>=2 and n_R=2q+1<3q^2.
For r>=2 the support begins at r-1, q>=r-1, and

    0<n_R<=2r*q-r^2+r+1<=3q^2.                      (14)

Writing q=r-1+t with t>=0, the last difference is exactly

    (2r-1)*(r-2)+(4r-6)*t+3t^2>=0.

This retains the sharp borderline r=2,q=1,n_R=3. Equations (12)-(14)
give m^2<3, hence m=1 for distinct pencils. Thus an index-two prime's
coordinate pencil is either the unramified pencil itself or complementary
to it. No assumption on pole heights, characteristic defects, finite field
degree or target normal positivity is used for this dichotomy.

## 5. A uniform pole constraint in the complementary sector

At m=1, (12),(13) also require n_R>q^2. For r=1, this gives q=2.
For r>=2 there is a sharper conclusion than the coarse estimate (14).
Put g=gcd(r,q). Precisely r-g nontrivial characters have j(a)=q;
the remaining g-1 have j(a)>=r-1. Thus

    S_R<=(r+g-2)*q-(g-1)*(r-1),
    n_R<=(r+g)*q-g*r+g+1,
    q^2-n_R>=(q-r)*(q-g)-g-1.                       (16)

If q>=r+1, set a=q-r>=1. Since g divides a, g<=a, and the last
expression is at least a*(r-1)-1>=0. This contradicts n_R>q^2.
The value q=r is impossible by the exact residue c_r=0. Hence

    r=1: q=2;   r>=2: q=r-1.                        (17)

In the latter row the entire principal part has just the one term
c_(r-1)*z^(-(r-1)); there is no higher or multicharacteristic Laurent
sector left on R. The r=1 row likewise retains only c_2*f^-2 and the
simultaneous constant A. This is a uniform geometric collapse of the
complementary sector, without restrictions on D order or pole and with
every source automorphism accounted for. The exact product (12) remains
an additional necessary condition.

With r=1,q=2, D order s<=2 is entirely excluded after PF2RULED48/49.
For s=1, (12) forces p=2, already excluded. For s=2, (13) gives
n_D<=4p-4; 5*n_D>4p^2 and p>=3 force p=3. Then gcd(p,s)=1,
Delta_D=3 and q=2>=s, exactly the PF2RULED49 exclusion. Thus a surviving
complementary simple-R quartic case must have q=2 and s>=3, as well as
(12). No exclusion of those higher D orders is asserted here.

If D has simple zero and r>=2, (17) gives n_R=r^2-r+1 and n_D=2p.
For r=2, (12) leaves only p=2,3,4,5. For r=3 it leaves p=2,3,
and PF2RULED48 excludes p=3; for every r>=4 it leaves only p=2.
Thus simple D and r>=3 require p=2, while the full problem still includes
higher D order and common rulings. These are necessary conditions, not
finite Keller constructions.

## 6. Sources, controls and the remaining finite-map problem

The [manifest](sources/coordinate-pencil-matrix/manifest.json) pins Xie's
Theorem1.2 and its completed proof, using the dual-divisor formulation
on PDF3; PDF32 containing the proof's final step was also visually checked.
The preliminary Section2 sketch expressly contains unproved steps and
is NOT imported as a proof. The whole-chart and coordinate-normal-form
inputs remain the sealed PF2RULED38/43/45/46/48/49 results.

The [checker](verify_coordinate_pencil_matrix.py) writes
[the receipt](coordinate_pencil_matrix_receipt.json). It replays actual
point-blowup matrices in complementary, one-shear and longer relative
coordinate systems, including shared initial centres, and verifies both
diagonal entries, the cross pairing and canonical orders. Symbolic
nonnegative identities prove the uniform bounds; finite controls are
support for the written universal argument. Negative controls retain the
failure of the unramified bound without exactness, common-pencilled affine
models, and the fact that a passing matrix is not an affineness certificate.

The quartic gate now has two coordinate configurations, common and
complementary. PF2RULED48/49 and (12) constrain the latter. The actual
common-ruling affine controls still have no finite Keller map. No Ga
action on the target, monogenic normal algebra, target ruling, external
review, independent full reproof or full-goal completion is claimed.
