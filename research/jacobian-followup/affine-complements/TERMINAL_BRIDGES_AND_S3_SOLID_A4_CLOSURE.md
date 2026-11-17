# Terminal bridges close S3 and still-solid A4 histories

15 November 2026. PF2RULED95 excludes ALL remaining incoming S3
complementary trajectories, at every solid-history depth and later
tail length. It also excludes ALL A4 trajectories still solid at
first simple arrival, through every tail length. Together with the
earlier common-pencil theorems, these sectors are closed in BOTH
pencils. The remaining A4 complementary cases have an EARLIER
nonsolid reset, with terminal pairs(2,8) or(4,10). A first cyclic3
loop after such a reset, when present, now has a finite winding list.
This is not a full quartic or research-goal closure.

## 1. A degree-free terminal bridge law

Retain the full smooth exact two-line complementary source and actual
unramified D chart of PF77/79/89. Suppose its canonical D path has
at least THREE genuine nodes. Its last winding is t>=2, free-tail
order kappa>=2, and penultimate winding w>=2. Write j for the
penultimate gcd-dropping support index, i for the preceding one
(or p when the preceding node is the first node). Then

    t divides sigma,p,i,j,
    t*w divides sigma,p,i,
    j_* = sigma-1+kappa, gcd(t,j_*)=1.

The two fibre-intersection numbers adjacent to the penultimate node,
toward the marked root and toward the terminal node respectively, are

    D_out=(i-j)/t, D_in=j-j_*>0.                       (1)

Its order-w cone has normalized coefficient

    +/- inverse(D_out) modulo w.                      (2)

These formulas retain the cone COEFFICIENT, not just its order.
Here is their all-depth local derivation. For u=1/f the successive
characteristic exponents are p,2p-j_1,2p-j_2,... . Earlier reciprocal
products have exponents divisible by the preceding characteristic
gcd, so the first nondivisible term at each drop cannot cancel.
Successive differences are i-j, including p-j at the first such edge.
Resolve the preceding truncation and subtract its integral terms.
The next local projection has order g_before and transverse first
order equal to that characteristic difference; dividing by g_after
gives its coprime Newton pair. The determinant-one Euclidean changes
transport adjacent primitive fibres with intersection difference/g_after.
In particular g_after=t at the penultimate drop and1 at the last,
giving(1). The same local pair at the penultimate node is
(w,D_out), modulo adding multiples of w to the transverse order.
Its regular-fibre slope is D_out*mu+w*ell. A primitive section
therefore gives(2). Earlier depth and lower noncharacteristic terms
do not alter these differences or residues.

This is the local Euclidean argument extending PF85/89's terminal
and first-edge calculations. Neumann--Wightwick math/9911104v1,
Sections2.1.5/2.2, independently checks the characteristic-difference
and edge-determinant convention; its signed rooted splice determinant
is not silently equated with our positive torus intersection. The
positive conversion and normalized residues are proved above and
checked on complete integral plumbing graphs. Retained Hatcher and
Neumann--Raymond primitive-section conventions remain imported.

Strip the LAST winding t exactly as in PF89. Set

    sigma'=sigma/t, p'=p/t,
    kappa'=j/t+1-sigma',
    delta'=(delta+h*(t-1)*(t+kappa))/t^2,              (3)

where the full source determinant is delta=h*n_D-q^2*p^2.
The stripped object is a primitive LOCAL chart, not an assumed
affine source. Since w divides sigma',p',i/t, formulas(1)--(3)
give

    D_out=1-kappa' mod w,
    D_in=t*kappa'-(t+kappa-1),
    delta'=h*kappa' mod w.                            (4)

In particular D_in=1-kappa mod t. The last congruence is the
actual terminal modulus for the stripped chart, which has at least
two nodes. The penultimate cone and its inner bridge therefore obey
the UNIVERSAL residue condition: if its target coefficient is
+/-inverse(gamma) modulo w, then

    t*gamma = +/- (D_in+kappa-1) mod w.                (5)

Indeed gamma=+/-D_out and t*D_out=1-kappa-D_in modulo w.
Allowing both signs is a weaker necessary condition, so a failure
excludes either oriented matching convention.

More generally, if target peripheral data give

    D_out=+/- f*D_in mod w,                            (6)

for an integer f, then(4) forces w to divide AT LEAST ONE of

    C_plus =(f*t+1)*delta'
                     -h*(f*(t+kappa-1)+1),
    C_minus=(f*t-1)*delta'
                     -h*(f*(t+kappa-1)-1).            (7)

This follows by multiplying the two signed equations in(6) by h
and substituting h*kappa'=delta' modulo w. No inversion of h is
used. If a constant is zero, that branch gives no finite bound;
no uniform divisor bound is claimed then. This source/marked-slope
lemma has no cover-degree, characteristic-depth or support cutoff.
Actual primitive coefficients, index order and canonical pieces are
required. It is not a classification of all higher-degree packets.

## 2. Every nonempty S3 solid history is excluded

PF90 leaves only solid trefoil-origin incoming S3 histories.
PF94 proves A=W^2, forbids a factor3 in EVERY solid winding and
later-tail winding, and keeps the terminal ordinary trefoil canonical,
including its product-root first-tail merger qualification. Any
remaining complementary trajectory has a nonempty later simple tail.
The source has rho=2,h=3,q=1; its terminal trefoil assignment MUST
be(winding,free)=(3,2). The alternative(2,3) cannot break the last
residual gcd2. The closed end is cyclic with exactly one factor3.

If the ordinary companion has any solid loop, consider its FIRST
loop next to that terminal trefoil. Its winding w>=3 is odd, with
3 not dividing w, and actual coefficient

    k0=2w+epsilon, epsilon=+1 or-1.

It has cone coefficient +/-inverse(k0) modulo w. Its inner fibre
intersection with the ordinary trefoil is

    D_in=6w-k0=4w-epsilon>0.                          (8)

There is an outer neighbour, either an arrival root or a later tail;
PF94 excludes the possible product mergers. Thus Section1 applies
with t=3,kappa=2,gamma=k0. Formula(5), reduced modulo w, gives

    3*epsilon=+/- (1-epsilon) mod w.

For epsilon=+1, this forces w|3, contrary to3 not dividing w and
w>=3. For epsilon=-1, it forces w|5 or w|1, hence only w=5.
But(4) also requires D_in=2 mod3; at w=5,epsilon=-1,
D_in=21=0 mod3. Contradiction. Therefore EVERY nonempty ordinary
solid S3 history is excluded, independently of later history depth,
arrival coefficient or tail length. No prefix affine source is used.

## 3. The depth-zero S3 cases are also excluded

Only the initial ordinary trefoil remains: W=A=1,c=1,B=6.
PF92/94 leave either r=3 or r=3s>=9 with absolute cap order1.
The actual closed cyclic group and R mark are still PF94's

    H1=Z/abs(D0), D0=k/2-s, [R]=s*k*T.

The source terminal winding3 and free2 give, from(3),

    delta=abs(D0)=9*delta'-30=6 mod9.                 (9)

This exact terminal-strip congruence is stronger than merely one
factor3. A nonempty tail ensures at least two D nodes, so stripping
the last winding is valid even if it leaves a single node.

For r=3, the genuine infinity bound gives 2<=k<18, k even and
gcd(3,k)=1. Thus D0=k/2-1 is between0 and7. D0=0 is already
excluded by positive first Betti number; (9) forces EXACTLY

    D0=6, k=14.

Its drilled root has cone abs(3-k/2)=4 and is not a product.
The canonical inner edge to the trefoil has intersection18-k=4,
contradicting(4)'s required2 modulo3. The product candidates k=4,8
instead have end orders1,3 and already fail(9). This finite base
deduction follows from the exact infinity bound after all histories
were excluded; it is not a bounded-depth search.

For r=3s>=9 the retained unit-cap reduction gives

    k=6s+2*zeta, zeta=+1 or-1,
    delta=2s+zeta, 3 not dividing s.

The root is an annulus with its one genuine cone s and an outer
tail neighbour. Its normalized cone coefficient is +/-inverse(2)
modulo s: PF81's cone coefficient is -v with r*u-k*v=1,
so -2*zeta*v=1 mod s. The actual inner fibre intersection is

    D_in=6r-k=12s-2*zeta.

Its residue modulo3 must be2 by(4), forcing zeta=-1.
Then(5), with t=3,kappa=2,gamma=2, requires

    6=+/-3 mod s, hence s|3 or s|9.

Since s>=3 and3 does not divide s, this is impossible. Thus ALL
depth-zero remaining S3 trajectories are excluded too.

Together with PF90/91 and PF93, ALL incoming S3 trajectories are
now excluded in BOTH source pencils, at arbitrary companion depth,
arrival parameters and later finite tail length, in the retained
genuine positive iterated-torus quartic class.

## 4. Still-solid A4 arrivals close through every tail

PF89 reduces a still-solid first simple arrival to terminal
(winding,free)=(2,4), source h=3,delta=2, and at least TWO later
simple tails. Empty and one-tail arrivals are already excluded.
Every remaining simple tail has odd winding, odd k, and winding
prime to3, by PF83's zero-R full-source lattice.

Take the FIRST simple tail next to that terminal arrival disk.
Its winding w>=3 is an internal D winding, with an outer tail
neighbour. Its cone coefficient is +/-inverse(2k) modulo w.
The ordinary arrival fibre parameter B0 is even, so PF83's actual
primitive-fibre determinant gives

    D_in=w*B0-k>0, 2k=-2D_in mod w.

Matching its cone coefficient with(2) requires(6) with f=2.
For terminal t=2,kappa=4, (3) gives delta'=5 exactly. Formula(7)
has the COMPLETE constants

    C_plus=-8, C_minus=-12.

Hence w divides8 or12. No odd w>=3 prime to3 can do so.
Contradiction. This closes EVERY still-solid A4 arrival with two or
more simple tails, at arbitrary solid-history depth and parameters.
The proof needs no separate enumeration of PF94's all-plus cases
or last-minus5 exception: ALL are excluded together.
With PF89's empty/one-tail and earlier common exclusions, this
ENTIRE A4 origin-to-still-solid-arrival sector is closed in BOTH
pencils at every finite length.

## 5. Earlier A4 resets have bounded first post-reset winding

Only PF89's earlier nonsolid-reset A4 complementary families remain,
with source terminal(t,kappa)=(2,8) or(4,10), h=3,delta=2.
PF94 forces an ALL-PLUS solid history before that reset, exact
A=W^2=1 mod24, and its four actual reset assignments. All those
conditions remain. There may be no cyclic3 loop after reset.

If such a loop IS present, consider the FIRST one immediately next
to the terminal reset disk. Let its winding be w>=2; w,k are units
modulo3. Its cone coefficient is +/-inverse(3k) modulo w. On the
index3 companion boundary, primitive fibres are

    H_old=B0*m+3L, H_new=k*m+3wL,
    D_in=3*(w*B0-k)>0.

Here B0 is prime to3 since both reset parameters are units modulo3;
the actual fibre degree is3. Thus 3k=-D_in modulo w and cone
matching requires(6) with f=1. Both source terminal pairs strip
to delta'=8. Formula(7) now yields the COMPLETE restrictions

| Source terminal(t,kappa) | Constants | First post-reset winding w |
| --- | --- | --- |
| (2,8) | -6,-16 | 2,4,8,16 |
| (4,10) | -2,-12 | 2,4 |

The lists use w>=2 and3 not dividing w. In particular an ODD first
post-reset cyclic3 winding is impossible. Subsequent nonsolid cyclic3
loops and later simple tails have not been declared closed or
silently bounded. If there is NO post-reset loop, the old one-simple-
tail exclusion remains, while two or more simple tails remain open.

## 6. Exact controls and the current frontier

The [checker](verify_terminal_bridge_closure.py),
[receipt](terminal_bridge_closure_receipt.json), and
[manifest](sources/terminal-bridge-closure/manifest.json) retain the
symbolic local identities, signed cone comparisons, complete divisor
lists, exact actual-cover controls and full integral source graphs.
The new [polynomial cluster helper](exact_poly_cluster.py) performs
the SAME two-curvette blowup algorithm as the retained symbolic
helper using exact numerator/denominator polynomials over QQ,
cancelling polynomial gcds. It uses no truncated series or numeric
jets. Small full-row comparisons with the prior helper and complete
three/four-node graphs check the local bridge and cone-cofactor
formulas. A larger completed four-node diagnostic is retained under
tmp/wright95; it is not rerun in the routine checker or used to infer
the universal theorem. All control blowup guards are finite and are
not asserted to prove general termination or the unbounded claims.

The current genuine iterated-torus quartic frontier is ONLY the
complementary earlier-reset A4 families (2,8)/(4,10), with every
PF89/94 restriction and the new first-post-reset winding lists.
PF93 closes the full common class. Initial full-S4, cyclic incoming,
incoming S4, S3 and still-solid A4 sectors are all closed by the
complete PF82 first-arrival graph and these proofs. KIG28 is not a
bridge identifying intermediate A4 with full generic monodromy.

The retained whole smooth source, entire two-line complementary
charts, actual preferred peripherals, genuine positive stages and
simple-meridian full-final-S4 hypotheses remain. Degenerate stages,
higher-degree packets, the remaining reset A4 families, full quartic
and the full research goal remain OPEN. Deep topology and peripheral
rigidity are imported, not independently reproved. No finite Keller
map, formal verification, external review or novelty claim is made.
