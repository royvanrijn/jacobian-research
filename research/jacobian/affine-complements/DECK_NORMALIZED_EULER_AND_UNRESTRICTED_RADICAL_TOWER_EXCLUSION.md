# Deck-normalized Euler bounds remove the radical-tower inertia restriction

PF2RULED133. This extends the credited knot/Wright/HC4 coupling in
[PF132](LOW_EULER_CYCLICITY_AND_DECORATED_RADICAL_TOWER_EXCLUSION.md).
The ENTIRE natural radical-tower action H_q, q>=2, is excluded for an
actual plane Keller map, now with ARBITRARY generic inertia. The old
decorated proof remains correct at its sealed scope. The new argument
uses its low-Euler surface lemma on a different open quotient; it does
not claim that arbitrary inertia has zero residual ramification.

More generally, every noncyclic subgroup of the generic deck group has
fixed-field degree at least FIVE over the plane target. This applies in
all original degrees and without degree minimality. The full finite-target,
trace-energy and EMPTY ORIGINAL fibre gates remain OPEN.

## 1. Use the actual deck quotient even when the relative map ramifies

Fix an actual polynomial complex plane Keller map, with field extension
L/F of degree d>1, full finite normalization N, and ENTIRE etale completion
U=U_N. Keep the whole original A2 and ALL omitted primes. Write b0 for
the unramified omitted primes, A for the total moved branch support,
and Eta for the nonnegative total finite excess. Retain PF129:

    chi(U)=1+b0=d-A-Eta.

Let G be ANY subgroup of Aut_F(L), with order g; it need not be the full
generic deck group or normal in it. HC4MYLA0/1 gives a regular FREE action
of G on U, so

    T=U/G is smooth affine,   pi1(T)=G,
    c=chi(T)=(d-A-Eta)/g is a positive integer.

T has negative logarithmic Kodaira dimension, because it is dominated
by the whole original A2. Its etale map to the plane target supplies
the regular exact nowhere-zero form dP wedge dQ. Thus PF132 applies
to T, and

    G NONCYCLIC ==> c>=3.                              (1)

Let M be the genuine normalization in E=L^G, with b=[E:F]=d/g, and
let V be the ENTIRE etale locus of M. T is an open subset of V; in
general it need not equal V. Apply PF130 to N->M->A2, retaining every
residual critical curve and its finite excess:

    g*chi(V)-chi(U)=Delta+Eta_res >=0,
    c=chi(V)-(Delta+Eta_res)/g <=chi(V).                (2)

Consequently

    chi(V)<=2 ==> G cyclic,                            (3)

even if N->M ramifies over V. This strengthens the earlier use that
required a finite etale relative cover of ALL V. The finite etale cover
used HERE is U->T; it is not incorrectly identified with h^-1(V)->V.
Nonzero residual defect lowers c, so it cannot rescue a noncyclic G
when chi(V)<=2.

This also applies to any given genuine factor with L/E Galois: its
relative automorphisms form just such a subgroup G. Without this Galois
premise, a noncyclic relative monodromy group is not a deck subgroup and
is not excluded by (3).

## 2. An exact integer budget in deck orbits

On the d generic sheets, G centralizes the transitive original monodromy
and is SEMIREGULAR: a commuting permutation fixing one letter fixes all
letters, so is the identity. Its orbits all have size g and form the
genuine b-block quotient corresponding to E=L^G.

For every original branch meridian sigma_C, its fixed set is G-stable.
Hence its support w_C is divisible by g. Put

    s_C=w_C/g,    S=sum_C s_C=A/g.

The quotient action moves a_C blocks, and s_C>=a_C. This includes
invisible inertia a_C=0 and movement inside fixed quotient blocks.
Since b and c are integers, the original Euler identity also proves
that J=Eta/g is a NONNEGATIVE INTEGER, although no individual local
excess or branchwise excess is asserted divisible by g. We have exactly

    b=c+S+J.                                         (4)

For noncyclic G, (1) implies b>=3. In particular b>1. Its connected
finite intermediate M->A2 must have a branch component: otherwise
purity and simple connectedness of A2 would make b=1. Some a_C>=2,
because a nonidentity permutation cannot move just one letter. Thus
S>=2, and (4) gives

    G NONCYCLIC ==> b>=5,   |G|<=d/5.                 (5)

This is a fixed-field INDEX bound, not a bound on the order of the
full monodromy. It excludes every relative noncyclic Galois group whose
fixed-field degree is at most4, without assumptions on its inertia,
residue degrees, smoothness of the full N, or degree minimality.
In particular every actual plane Keller map of generic degree d<20
has CYCLIC generic deck group, since a noncyclic finite group has order
at least4. It does not assert inversion for those degrees.

There is a useful support refinement. Let mu be the minimum support of
a nonidentity element of the ACTUAL original permutation monodromy.
At least one branch exists, so A>=mu, and

    G NONCYCLIC ==> 3*g <= d-mu.                      (6)

Thus a half-support-minimum action requires d>=6*g. This criterion can
be tested without enumerating its possible inertia classes.

## 3. The entire natural H_q family is excluded

Retain the knot construction's exact permutation action

    H_q=Aff(W,R):S4 on W x R,
    W=F2^2, R=F2^q, d=4*2^q.

The group of affine FUNCTIONS W->R acts by fibre translations. PF132's
centralizer proof, independent of any branch decoration, gives

    Aut_F(L)=C_Sym(W x R)(H_q)=R,
    g=2^q,    [L^R:F]=4.

For q>=2 this deck group is noncyclic, contradicting (5). Therefore
NO actual plane Keller map has THIS NATURAL H_q monodromy action in
ANY generic inertia classes. All degrees16,32,64,... are excluded.

One can also check (6) directly: any nonidentity base permutation moves
at least two of the four whole R-blocks. If the base permutation is
identity, a nonzero affine function W->R vanishes on at most two points
(three distinct points of F2^2 determine it). Thus mu=d/2. The resulting
bound g<=d/6 contradicts g=d/4. Both arguments keep the original source
and all omitted primes; neither silently removes a residual branch.

Order-four lifts and vertical radical permutations remain perfectly valid
elements of H_q and can have positive residual defect. They are NOT
rejected as abstract knot-cover permutations. The different quotient T
and its Euler bound now exclude their actual plane-surface realizations
in the full natural monodromy family. The q=1 case has cyclic deck C2
and is not excluded. Other permutation actions of the same abstract
group are outside this conclusion. The knot representations remain valid.

## 4. Exact first possible noncyclic-deck profile

If equality b=5 occurs in (5), write (4) as

    0=(c-3)+(S-2)+J.

All terms are nonnegative, so c=3, S=2 and Eta=0. At least one visible
branch has a_C>=2; every other nontrivial original branch would cost
at least one further s_C. Hence there is exactly ONE original branch,
its quotient inertia is a transposition, w_C=2g and a_C=2, with Delta=0.
Transitivity and meridian normal generation make the degree-five quotient
the natural S5 action. PF130 then gives

    chi(V)=3, Eta_M=0, U=h^-1(V), U->V finite etale,
    b0=3g-1.

The retained zero-excess geometry makes the intermediate smooth with one
index-two critical A1 and immersive image normalization. Original critical
residues are1 and original N is smooth as well, by the same PF129 criterion.
These are NECESSARY conditions, not an actual cover construction or an
exclusion of this index-five profile. They show exactly where the general
argument stops, without proposing a new parameter campaign.

For binary primitive intermediates, PF124 gives chi(V)=r*(2-kappa),
b=4r, with kappa=0 or1. Thus any Galois relative group above degree4,
and above a degree8 one-cusp intermediate, is cyclic by (3), EVEN WITH
residual ramification. The remaining degree8 immersive profile is not
excluded by this numerical test.

## 5. Provenance and assurance

The low-Euler lemma is credited to the coordinator coupling packet and
HC4MYML0, as adopted with its classical inputs by PF132. The present
extension uses only the retained HC4 generic-deck geometry and PF129/130
Euler identities. It imports no new surface/group classification.
PF132 and its exact primary pins are retained unchanged; no predecessor
mathematical checker is rerun. Small exact controls check the integer
ledger, semiregular support divisibility, residual-defect direction,
the same valid nondecorated permutations and the index-five boundary.
They support the written unbounded proof, not a finite-sample theorem.

A literature comparison with Blanc's *Symplectic birational transformations
of the plane*, Osaka J.Math.50(2013), Section1.1, concerns preservation of
dx wedge dy/(xy). It supplies no classification for the ACTUAL exact
polynomial volume used here. No finite-Cremona linearization assertion
or consequence is imported from that comparison.

Full Wright/HC4/JC2, natural A/S or full-AGL compatibility, cyclic generic
deck cases, index-five noncyclic profiles and the trace-energy/EMPTY
ORIGINAL fibre gates remain OPEN. Formal verification, independent replay,
external review, novelty and full-goal completion remain unclaimed.
