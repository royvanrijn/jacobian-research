# Cyclic caps, genus deficits and prime-mark history reductions

15 November 2026. PF2RULED94 gives a degree-free square recurrence
for actual cyclic solid caps and a common genus-deficit argument for
the S3 and A4 histories. It restricts the complementary A4 frontier:
every earlier nonsolid reset has an ALL-PLUS solid history; a still-solid
arrival can have a minus history only at its LAST solid winding5,
followed by the smaller arrival parameter. The later-tail families
are not declared closed. A second degree-free cyclic-lattice lemma
forces exact3-primary restrictions on every remaining complementary
S3 trajectory. PF93's common-class closure is unchanged.

## 1. The actual cyclic cap and its square recurrence

Consider a positive(r,k)-cable pattern, r,k>=2, gcd(r,k)=1, in its
standard connected cyclic n-cover. Require n>=2 and BOTH r,k units
modulo n. Keep actual outer/inner meridian/preferred-longitude bases

    (m_o,l_o)=(lift(M^n),lift(L)),
    (m_i,l_i)=(lift(mu^n),lift(lambda)).

The lifted cable exterior has annulus base and ONE cone r. Indeed
the downstairs regular fibre has cyclic image M^k, of exact order n;
the generic fibre degree is n, the base degree is1, and the cone r
is preserved because gcd(n,r)=1. The primitive OUTER regular fibre is

    H=k*m_o+n*r*l_o.                                  (1)

The coefficients are coprime. The annulus-one-cone presentation has
free H1 of rank2, since its cone relation is primitive. Thus PF90's
unfilled integral transfer has no added summand here and gives

    H1=Z*m_i plus Z*l_o,
    m_o=r*m_i, l_i=r*l_o.                             (2)

Attach an ACTUAL solid companion with primitive cap A*m_o+b*l_o,
b>0, gcd(A,b)=1. Its meridian-fibre distance is

    J=abs(n*A*r-b*k).                                 (3)

The capped piece is disk(r,J), with order1 absorbed. In particular
J=1 makes it an ACTUAL solid torus, by its Seifert model. This is
not a recognition from cyclic homology. At J=1, gcd(r,b)=1, and
the relation A*r*m_i+b*l_o=0 in(2) has primitive coefficients.
Writing its generator as t gives

    [m_i]=b*t, [l_o]=-A*r*t, [l_i]=-A*r^2*t.

The new primitive meridian is therefore

    A_new*m_i+b*l_i, A_new=r^2*A.                     (4)

This is an all-degree statement about the specified cyclic restriction
and actual solid cap. Other boundary pieces must be genuinely capped
as specified before using it. It is not a classification of arbitrary
higher-degree representations, nonsolid companions, or k nonunit n.
It preserves the coefficient b, without replacing a nonunimodular
cap by coefficient1. The sign of A is fixed by the actual preferred
longitude. We use A>0 in the genus statements below.

## 2. A unified genus deficit and a last-minus obstruction

Suppose this solid-cap iteration is part of a genuine positive
ordinary iterated torus knot, with canonical last fibre parameter B,
c=2g(J)-1, and strict infinity inequality k<r*B. At a unit cap write

    b*k=n*A*r+epsilon, epsilon in{+1,-1}.

By(4), the ordinary cabling genus identity, and B_new=r*k,

    A_new=r^2*A,
    c_new=r*c+(r-1)*k,
    D=n*A-b*c,
    D_new=r*D-(r-1)*epsilon,
    B_new=(n*A_new+epsilon*r)/b.                      (5)

If D>=1, then D_new>=1. Moreover

    k/r-c=(D+epsilon/r)/b>0.                          (6)

Hence a positive L-space ordinary companion remains a positive
L-space knot at EVERY such solid continuation. This uses the
retained primary Hedden/Hom cabling theorem (Hom0912.4046v1,
Introduction Theorem); it is not independently reproved here.

After a MINUS loop of winding w, B=(n*A-w)/b. Any putative next
solid loop of winding r and sign epsilon' would require

    k_next-r*B=(epsilon'+r*w)/b>0,                    (7)

contradicting the strict infinity inequality. Thus a minus sign can
occur ONLY at the LAST solid loop, at every depth, with no parameter
cutoff. Starting at D=1, all-plus histories keep D=1; a last minus
of winding w has D=2w-1. The accumulated cap is A=A0*W^2,
where W is the product of the solid windings. No proper prefix is
assumed to be an affine source, and no prefix eta bound is imposed.

PF91's solid S3 histories are n=2,b=1,A0=1,D0=1.
The A4 histories below are the same argument with n=3,b=2.
The formulas do not assert that all cyclic-cover histories start
with such a solid cap or have positive deficit.

## 3. Exact A4 normalization at every solid depth

By PF87/89, any remaining complementary A4 trajectory originates
at the trefoil T(2,3), not a later A4 birth or non-trefoil initial
torus. Its actual index3 companion, after the true fixed-sheet
meridian cap, is PF87's disk(1,1,1,3), hence a solid torus with
one genuine cone3. In its actual basis
m=lift(M^3), L=lift(lambda), its primitive cap is

    m+2L.                                            (8)

To fix the coefficient and its sign, the actual trefoil cover
presentation a^2=b^3, mu=a*b^(-1), lambda=a^2*mu^(-6), after
the true fixed meridian fill, has H1=Z*t with m=2t, L=-t;
its kernel is exactly(1,2). All initial A4 epimorphisms form ONE
simultaneous S4-conjugacy orbit: a is one of the three double
transpositions, b one of the eight3-cycles; the24 generating pairs
are exactly the24 conjugates of any one pair. The checker exhausts
these finite possibilities, including the exchanged torus exponents.
This is a finite base-case certificate, not a sampled inference
about arbitrary histories. It also retains the unfilled
integer rows, the fixed-sheet fill and these primitive images.
The Seifert model supplies solidity, rather than those homology data.

PF89 already excludes every cyclic3 loop with3|k by its unmarked
homology summand and zero-mark source constraint. All remaining
loops have r,k units modulo3; their active restriction is the
standard cyclic3 pattern. The degree1 fixed-sheet piece performs
the genuine fixed meridian cap, as in PF87. Sections1--2 therefore
apply. A solid loop satisfies

    k=(3*A*r+epsilon)/2, A_new=r^2*A.                 (9)

In particular r is odd and prime to3, hence r>=5. At ALL solid
depths, before the first nonsolid outcome,

    A=W^2, W odd and prime to3, A=1 mod24,
    D=3A-2c>=1,
    all-plus: D=1, c=(3A-1)/2,
    last-minus w: D=2w-1, B=(3A-w)/2.                (10)

Every ordinary companion is a positive L-space knot. A minus loop
is last among the solid loops. Subsequent nonsolid loops are not
subject to the solid-cap recurrence and are not silently omitted.

## 4. Still-solid first simple arrival: only one minus possibility

PF89 leaves a still-solid A4 arrival only with r_f=2,k_f=3K,
K odd, terminal pair(winding,free)=(2,4), and at least TWO later
simple tails. Its cap order4 gives

    K=A+2*zeta, zeta in{+1,-1}.

Its ACTUAL arrival degree is, by PF91's all-depth genus bridge,

    eta_arrival=k_f-2*c-1=D+6*zeta-1.                (11)

For an all-plus history this is +6 or -6. The smaller parameter
requires A>2, so the empty history cannot use it.

For a last-minus history of winding w, the strict infinity inequality
becomes

    3A+6*zeta<2B=3A-w, hence 6*zeta<-w.

The plus choice is impossible. The minus choice requires w<6;
by(9), w is odd, at least5, and prime to3. Therefore EXACTLY

    w=5, K=A-2, D=9, eta_arrival=2                  (12)

is the only possible minus history. In particular25 divides A.
These are necessary data, not a realization; even this case still
needs at least two later simple tails and every PF89 restriction.
There is NO eta0 still-solid arrival in the surviving terminal class.

## 5. Earlier nonsolid reset: every pre-reset solid loop is plus

For a first nonsolid cyclic3 cap, let r be its winding and
J=abs(3Ar-2k). PF89's terminal comparison leaves only ordered source
pairs(2,8) or(4,10). Consequently the ACTUAL unordered reset orders
are{r,J}={2,8} or{4,10}. Both possible assignments are retained.

If the preceding solid history ended with a minus winding w>=5,
write2k=3Ar+zeta*J. The infinity inequality would require

    zeta*J<-r*w.

The positive choice is impossible, and the negative one requires
J/r>w. But across the COMPLETE four assignments, J/r is at most4,
less than w. Contradiction. Thus the ENTIRE solid history before
any surviving nonsolid reset is ALL-PLUS. In particular A=W^2=1 mod24
and c=(3A-1)/2, without a depth bound or a coefficient search.

The possible reset parameters and reset degrees now have the exact
following form. The sign in each column is the sign of zeta.

| r | J | k, positive sign | eta_reset | k, negative sign | eta_reset |
| --- | --- | --- | --- | --- | --- |
| 2 | 8 | 3A+4 | 4 | 3A-4 | -4 |
| 8 | 2 | 12A+1 | 4 | 12A-1 | 2 |
| 4 | 10 | 6A+5 | 6 | 6A-5 | -4 |
| 10 | 4 | 15A+2 | 6 | 15A-2 | 2 |

Indeed eta_reset=k-r*c-1=(r+zeta*J)/2-1. Positivity,
coprimality and the actual infinity bound remain required; the
table does not assert that every row/sign is realized. These are
degrees at a CYCLIC3 RESET, not necessarily at first simple arrival
or the final stage. No final-stage eta bound is imposed on them.
All later nonsolid cyclic3 loops, the first simple arrival and the
later simple tails retain PF89's order and length restrictions.

## 6. A degree-free prime-mark lemma for the full source lattice

Consider the full two-line source presentation with marked first
covector R, cross q*p, and positive determinant

    Q=[[-h,q*p],[q*p,-n]], delta=h*n-q^2*p^2>0,
    gcd(h,q)=1.                                      (13)

Suppose coker Q is CYCLIC. For a prime ell dividing BOTH h and p,
Smith reduction forces ell not dividing n: otherwise every matrix
entry is divisible by ell and the first Smith factor is nontrivial.
Killing R leaves a cyclic group of order gcd(q*p,n), prime to ell.
Thus R generates the ENTIRE ell-primary subgroup of coker Q.

Write a=v_ell(h), t=v_ell(p). The two summands in delta have
valuations a and2t. If a!=2t, then

    v_ell(delta)=min(a,2t).                           (14)

If a=2t, cancellation can occur; no equality is asserted there.
In particular, if ell divides h EXACTLY once and divides p, then
ell divides delta EXACTLY once and R generates that order-ell
subgroup. This is a general integral lattice/marked-group statement,
not a rational linking approximation. Cyclicity is essential.

## 7. Every remaining S3 survivor has an exact single factor3

Retain PF92's remaining complementary S3 histories, all with solid
incoming cap A=W^2, an ordinary positive iterated-torus companion
whose terminal disk is the trefoil, and at least one later simple
tail. Even-k tails are impossible by PF83's canonical-root budget,
with the product-root qualification checked in the next paragraph.
Write T for the
product of their odd windings, and r=3s for the first simple arrival.
The complete actual closed marked homology is

    H1=Z/abs(D0), D0=k/2-A*s !=0,
    [R]=s*k*T times its generator.                    (15)

This is PF73's actual solid-cap calculation, retained through the
odd-k tails by PF82. The coefficient sign does not affect the prime
argument. D0=0 is already excluded by positive first Betti number.
The later tail's normal/fibre distance2 forces the source rho=2,
so h=3,q=1 in(13).

Here is the product-root qualification. At r=3,abs(J0)=1 the seed
E_R is the ordinary E(J), and its normal slope is t*n+L1 with
t=2A +/-1; only an all-plus history can have such a unit cap.
For a nonempty solid history, its last ordinary fibre has normal
distance B-t=w_last minus/plus1, an EVEN integer. An odd-k first
tail could merge fibres only if that distance equals2*w_tail,
hence w_last=2*w_tail +/-1. The merged annulus then has two
NONCONSECUTIVE cones(w_last,w_tail), so fails the source corner
budget if final, or the internal budget if another tail follows.
An even-k tail has old-boundary fibre distance w_tail ODD, so
cannot merge with this even distance. For the empty ordinary
history the seed is a trefoil; any merger adds at least one genuine
cone to its terminal disk(2,3), forbidden with or without later tails.
Thus all possible mergers in the product case are excluded. All
surviving first-tail edges are genuine; subsequent edges persist
by PF83's actual determinant. In particular the terminal trefoil
used below remains canonical, and all surviving tails have odd k.

There are at least two D nodes. The target terminal trefoil orders
are(2,3). Source assignment(winding2,free3) is impossible: the last
index sigma-1+3 is even, while the residual gcd2 must be broken.
Thus the source terminal winding is3 and its free order is2,
as in PF91. Its pole p is the product of ALL D windings, so3|p.
Equation(15) makes the end homology cyclic. Section6 now forces

    v_3(abs(D0))=1, and R generates its3-primary part. (16)

This excludes BOTH 3 not dividing D0 and9 dividing D0, regardless
of tail length, without a bounded determinant search.

For r=3, D0=k/2-A. If3|W, then3|A and gcd(3,k)=1 makes
3 not divide D0, contradiction. Hence3 not divide W. For r>=9,
PF92 already forces an all-plus history and k=2Ar+2*zeta,
zeta=+1 or-1. Then D0=2As+zeta; again3|W is impossible,
and(16) also requires3 not divide s. In ALL remaining cases,
gcd(3,s*k)=1. The mark in(15) generates the3-primary part only
if3 not divide T. Consequently

    EVERY solid-companion winding and EVERY later simple-tail
    winding is prime to3; at r>=9, also3 not divide s. (17)

Equivalently every such odd winding is1 or5 modulo6. The exact
square cap also satisfies A=1 modulo24. The exact
arrival restrictions are

    r=3: k/2-A=3 or6 mod9;
    r>=9: 2As+zeta=3 or6 mod9,
           k=2Ar+2*zeta, all-plus history.             (18)

The canonical source pole can also be read exactly from its rooted
path: p=3*W*C*T, where C=abs(3A-k/2) at r=3 (omit C=1),
or C=s at r>=9's unit cap. Each edge persists by PF83; the product
root contributes no extra winding. Equations(16)--(17) imply
3 not dividing C and v_3(p)=1. This exact pole formula is an
additional necessary match, not a realization of its chart.

The lattice argument uses only the full marked source group,
canonical winding product and actual target mark, not its Floer
homology. It does not extend the common-pencil LO theorem to a
complementary source or assume that the latter is an L-space.

## 8. What changes, what remains open

The old arbitrary odd solid-cap coefficient has become the precise
positive square A=W^2=1 mod24, with its ordinary genus and last fibre
parameter known at every depth. This closes ALL minus solid histories
before the earlier-reset complementary A4 families, and ALL minus
still-solid arrivals except(12). It provides a general cap/genus
mechanism rather than a longer enumeration of local coefficients.

The remaining complementary A4 families are:

* Still solid until first simple arrival: terminal(2,4), at least
  two simple tails; all-plus history and K=A +/-2, or the last-minus5
  exception(12).
* Earlier nonsolid reset: terminal(2,8) or(4,10), all-plus pre-reset
  history, one of the four reset assignments above, then the retained
  possible nonsolid cyclic3 loops and later simple tails. PF89's
  single-tail/no-intervening-loop exclusion and final winding5 mod6
  remain in force.

Complementary S3 survivors retain PF92's arrivalr3 or all-plus r>=9,
k=2Ar +/-2 with later tails, NOW with(16)--(18): exactly one factor3
in the end order,3-primary generating R mark, and no factor3 in any
solid-loop or later-tail winding. PF93's all-length COMMON genuine
class is closed. No complementary full closure, higher-degree packet
classification, degenerate-stage exclusion, finite Keller map, full
quartic exclusion or completion of the research goal is claimed.
KIG28's full generic A4 Keller action is not identified with these
intermediate/peripheral states.

The [checker](verify_cyclic_solid_cap_genus.py),
[receipt](cyclic_solid_cap_genus_receipt.json), and
[manifest](sources/cyclic-solid-cap-genus/manifest.json) retain
symbolic all-degree identities, actual cyclic cap controls in several
degrees, unexpanded composed A4 histories, the exact initial cap,
full arrival/reset/tail controls and cyclic-lattice prime countercontrols.
The written arguments establish
universality, not the finite controls. The retained topology, genus
and primary Floer inputs are imported; no formal verification,
external review, independent full reproof or novelty claim is made.
