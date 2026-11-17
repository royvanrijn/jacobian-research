# Low Euler cyclicity excludes the specified knot radical towers

PF2RULED132. Integrated from the coordinator's
[written coupling proof, Sections1--2](../cross-programme-connections/INWARD_COUPLED_RESTRICTIONS_2026_11_16.tex),
with its [scope handoff](../cross-programme-connections/inward-coupled-restrictions-handoff-2026-11-16.json).
The coupling argument is credited to that packet; this note audits and
adopts its Wright/HC4 deduction. Its Mathieu/Factorial sections are NOT
adopted by this claim. No independent review or novelty is asserted.

The conclusion is an INFINITE FAMILY EXCLUSION: natural monodromy
H_q=Aff(F2^2,F2^q):S4 is impossible for an actual plane Keller map when
q>=2 and EVERY generic branch permutation is one of the specified lifted
transpositions below. Degrees are16,32,64,..., with no upper bound or
degree-minimality assumption. Other inertia types for the same abstract
group, q=1, and the full finite-target problem remain OPEN.

## 1. The low-Euler exact-volume surface lemma

Let X be a smooth connected AFFINE complex surface satisfying ALL of:

* negative logarithmic Kodaira dimension;
* finite topological fundamental group;
* Euler characteristic chi(X)<=2;
* a regular nowhere-zero EXACT algebraic two-form.

Then pi1(X) is CYCLIC. Euler characteristic alone does not give this
conclusion. Nor does trivial canonical bundle replace exactness.

We adopt the classical existence of an A1-fibration rho:X->B from
Dubouloz--Kishimoto, arXiv:1408.1328v1, Section1.2 and Theorem2,
printed/PDFpage3, recalling Miyanishi--Sugie. Over C its finite constant
field extension is trivial. The handoff's Section1.1.2 locator corresponds
to different version numbering; the pinned v1 locator is Section1.2.
The whole relevant page and theorem proof were read, not the original
Miyanishi--Sugie proof independently reproduced.

Retain HC4MYML0's precisely sourced structure package: fibres are disjoint
unions of A1s with positive multiplicities, and irreducible reduced fibres
give a Zariski A1-bundle. For a fibre with multiplicities m_ij put
m_i=gcd_j(m_ij). The fundamental-group presentation is the base group
with the extra relations gamma_i^m_i=1. Here is the needed topological
argument: off the finitely many special fibres the A1-bundle has
contractible fibre. Restoring each smooth divisor component kills its
meridian. At a general point its local base parameter is a unit times
u^m_ij, so that meridian projects to gamma_i^m_ij. Their common normal
closure imposes exactly the gcd relation. This also matches the retained
Miyanishi2005 report, Lemma2.3, printed1111/PDF39.

The presentation surjects onto pi1(B). Finiteness leaves B=A1 or P1.
For A1 it is a free product of cyclic groups, hence finite only with
at most one nontrivial factor. This case is already cyclic.

For P1, Euler integration, ignoring scheme multiplicities, gives

    chi(X)=2+sum_i(s_i-1),

where s_i>=1 is the number of reduced components of the fibre. Thus
chi(X)<=2 forces chi(X)=2 and EVERY fibre irreducible. Its fundamental
group is the orientable sphere-orbifold group with cone orders m_i.

Suppose this finite group is noncyclic. Adopt the classical good/bad
orientable two-orbifold classification in Thurston, *The Geometry and
Topology of Three-Manifolds*, Chapter13, Theorem13.3.6 and its adjacent
table, printed312/PDF18. The relevant complete page was read. A finite
noncyclic sphere-orbifold group is a good spherical triangle group;
there is a finite holomorphic cover P1->P1 with EXACT local ramification
orders m_i. The bad teardrop and unequal spindle orbifolds do NOT have
that cover; their groups are already cyclic and were not put in this case.
The spherical triangle signatures are (2,2,n),(2,3,3),(2,3,4),(2,3,5).
This is an explicit adopted classical classification, not a new proof.

Normalize the corresponding base change to obtain X'->X. It is connected:
the generic field is C(B)(z), so extending the base field preserves its
geometric integrality. It is finite. Away from special fibres it is
etale. At a special fibre m_i*C_i the base ramification order is exactly
m_i; the local normalized root extension has ramification index
m_i/gcd(m_i,m_i)=1. Thus it is unramified in codimension one. Purity over
the smooth X makes it finite ETALE everywhere, so X' is smooth affine.
This purity step retains possible codimension-two points; a generic local
equation alone is not a proof of global etaleness.

The induced ruling X'->P1 has REDUCED fibres. They MAY SPLIT into several
A1 components. The exact nonvanishing form pulls back to one of the same
kind on X'. Therefore the FULL component-selection form of HC4MYML0
applies: choose a multiplicity-one component in every fibre, obtaining
an open A1-bundle over P1, which exact volume in an affine surface forbids.
This contradiction proves cyclicity. No irreducibility of pulled-back
fibres, nonexistent cover of a bad orbifold, or affine-base ruling was
assumed.

## 2. The knot tower, with its branch decoration retained

Credit the [knot construction](../meridional-rank/LIFT_DIMENSION_AND_RADICAL_TOWERS.tex),
Section2. Put W=F2^2, R=F2^q and

    H_q=Aff(W,R):S4,   |H_q|=24*2^(3q),   D=4*2^q.

The normal additive subgroup consists of AFFINE FUNCTIONS W->R, not
just constant translations. For an edge e={a,b} of the four points of W,
write s_e for its transposition and delta_e for the indicator of {a,b}.
The specified generic branch permutations are EXACTLY

    T_(e,r)(x,z)=(s_e(x), z+delta_e(x)*r),   r in R.

Each is an involution with support D/2 and fixes EVERY sheet in EACH
fixed base block. This last condition is essential to the proof.

The full permutation centralizer of the transitive H_q action is R,
acting by CONSTANT radical translations. Indeed a commuting permutation
must commute with those translations, hence is (x,z)->(f(x),z+c(x)).
Commuting with every affine-function translation forces f(x)=x, since
such functions separate base points. Commuting with the base S4 makes
c constant. For q=0 the natural S4 centralizer is trivial. The constants
are central in H_q and act regularly on each of its four radical blocks.
This is a deck group, not the full monodromy H_q.

## 3. Wright's support budget and the entire quartic intermediate

Assume an ACTUAL plane Keller map has this natural monodromy and the
specified branch permutations. Retain its whole original A2, full finite
normalization N, and ALL omitted primes. PF129 gives

    sum_C supp(sigma_C)+Eta <= D-1.

Every summand is D/2. Two branch components are impossible, while zero
components would give a connected finite etale cover of A2 and degree1.
There is therefore ONE original branch component.

The four radical blocks define a GENUINE degree-four intermediate
N->M->A2, relative degree t=2^q, with quotient group S4. No intermediate
is substituted by a plane. Its meridian is a transposition of support2.
For the ENTIRE etale locus V of M, PF129 gives exactly

    1 <= chi(V)=4-2-Eta_M <= 2.

The retained KU/PF128 construction makes V smooth affine, with finite
pi1, negative log Kodaira dimension (domination by the original plane),
and exact nonvanishing volume dP wedge dQ=d(P*dQ). These verify every
hypothesis of Section1.

There is NO residual ramification over V. The specified inertia fixes
every original sheet in each fixed quotient block; thus its residual
support defect is0. There is no invisible branch image either, because
every permitted generic inertia projects to a nonidentity transposition.
PF128's exact criterion, or PF130's relative Euler formulation, gives

    Z=h^-1(V) -> V finite etale,   Z=U_N.

The whole original A2 is open in Z. Its complement consists of disjoint
unramified A1 divisors; moving loops off them makes pi1(Z)=1. The constant
translations R give t distinct relative deck automorphisms, so their
degree-t cover is Galois and universal. Equivalently, its invariant field
is exactly the four-block intermediate field. Hence

    pi1(V)=R=(C2)^q.

Section1 forces this group cyclic. For q>=2 this is a contradiction,
excluding the ENTIRE DECORATED family in degrees16,32,64,..., without
minimal-counterexample induction or a degree cutoff.

## 4. General consequences and exact boundaries

Every genuine intermediate etale locus V with chi(V)<=2 has cyclic
fundamental group, by the same surface lemma and retained KU/PF128
hypotheses. If a relative cover over that ENTIRE V is finite etale,
its relative deck group is cyclic. Residual ramification cannot be silently
discarded in applying this statement.

The same lemma strengthens HC4MYLA1's generic deck cost. For ANY subgroup
G of the actual generic deck group, U_N/G is smooth affine, has pi1=G,
negative log Kodaira dimension, and exact nonvanishing volume descended
from the plane target. Thus

    G noncyclic ==> c=chi(U_N/G)>=3,
    b0+1=|G|*c ==> b0>=3*|G|-1.

Here b0 counts ALL unramified omitted primes; ramified ones are separate.
This improves the old 2*|G|-1 bound using the NEW low-Euler lemma. The
quotient need not be the ENTIRE etale locus of its finite intermediate;
the surface lemma still applies to it. No preservation of the original
plane chart or permutation of its boundary basis is assumed.

The tower exclusion RETAINS its branch-permutation hypothesis. Other
inertia classes can retain residual ramification, so the same abstract
H_q alone is not excluded here. The q=1 case has cyclic deck C2 and is
not excluded by this argument. Knot representations remain valid: knot
complements lack these additional affine-surface constraints. PF131's
proper-linear primitive quotient result is complementary: the present
H_q action has invariant radical blocks and is imprimitive for q>0.

The coupling packet and its exact controls are pinned without rerunning
their unrelated Mathieu calculations. The new checker verifies only the
relevant source seals, orbifold/Euler arithmetic, decorated permutations,+centralizers and residual-defect countercontrols. Written proofs plus
the named imports establish the unbounded conclusions. The classification
inputs and HC4MYML0 are credited, not re-proved by finite checks.

The full Wright/HC4 finite-target, trace-energy and EMPTY ORIGINAL fibre
gates remain OPEN. Formal verification, independent replay, external
review, historical priority and full-goal completion remain unclaimed.
