# Complete ruling charts and the Euler-two multiplicity pair

HC4MYMO0. 16 November 2026. Written all-degree deductions, conditional
on the stated ruling multiplicities. HC4/JC2 and the original full
trace-pairing/marked-energy implication remain OPEN. No independent
replay, formal verification, external review or novelty is asserted.

## 1. A complete ruling with at most one bad fibre supplies a plane

Let X be a smooth connected complex AFFINE surface with an affine
A1-fibration rho:X->P1. Call a fibre bad when NONE of its components
has multiplicity one. Assume at most ONE fibre is bad. Then X contains
an OPEN A2.

Here and below reuse HC4MYML0's local structure input: reduced fibre
components are disjoint A1s, there are finitely many degenerate fibres,
and an affine A1-fibration with every fibre one reduced A1 is a Zariski
A1-bundle. Its precise source and proof-sketch limits are retained;
this is not a new adoption or independent proof of that package.

Delete the entire bad fibre, or any fibre if there is none. The remaining
base is A1. In every remaining degenerate fibre choose ONE component
of multiplicity one, and delete its other components. These are finite
unions of closed divisors. Since X is smooth they are effective Cartier
divisors, whose complements are affine open immersions. Thus the
selected map to A1 is affine. It is flat over the smooth base curve
(the source is integral and the map is dominant), and ALL its fibres
are one reduced A1. The adopted structure theorem makes it a Zariski
A1-bundle. Such a bundle is trivial because Pic(A1)=H1(A1,O)=0.
The resulting open is A2.

This lemma does NOT require an exact two-form. If X does carry an exact
nowhere-zero two-form, HC4MYML0 excludes zero bad fibres, but allows
exactly one. The retained higher-power Danielewski controls already
have exactly one bad complete fibre outside a whole plane chart.
Neither an additive action on all X nor X=A2 follows from this lemma.

## 2. Euler two compresses the affine-base test to one pair

Let V be a smooth connected complex AFFINE surface with negative log
Kodaira dimension, finite fundamental group, chi(V)=2, and an exact
nowhere-zero algebraic two-form. These hypotheses hold for a genuine
entire intermediate etale locus of Euler two in the retained plane
Keller setting. PF2RULED132 supplies an A1-ruling and cyclic pi1(V).
Every chosen ruling has base A1 or P1, because its fundamental group
surjects onto that of its smooth base curve.

First suppose rho:V->A1. Euler integration gives

    2=1+sum_i(s_i-1),

where s_i is the number of reduced A1 components. Exactly ONE fibre
has two components, and all other fibres are irreducible. Write its
multiplicities as mu,nu, and put

    g0=gcd(mu,nu),  (a,b)=(mu/g0,nu/g0).

The pair (a,b) is coprime. Every irreducible fibre automatically has
a component equal to its whole-fibre gcd. Thus PF2RULED130 Section3's
component criterion for THIS ruling is equivalent to

    min(a,b)=1  <=>  mu divides nu OR nu divides mu.             (1)

When (1) holds, the universal finite etale cover W->V contains an
OPEN A2. Explicitly, finite pi1 permits at most one whole multiple
fibre; its gcd m gives pi1(V)=C_m. The normalized degree-m base root
cover used in PF130 is the universal cover. At its branched fibre a
component of multiplicity m*c lifts to m disjoint A1s of multiplicity
c. Other special fibres are reproduced over unramified base points.
Condition (1) gives a unit-multiplicity component in every W fibre,
and PF130's selection over A1 supplies the plane.

When (1) fails, a,b>=2 and gcd(a,b)=1. This is a failure of this
specified ruling construction, NOT a proof that W has no plane chart
or that every other ruling fails. Exactness does not make (1)
automatic: PF130 Section5 already supplies the smooth exact-volume
affine surface x^2*z=(y^2-x)*((y-1)^3-x), with special multiplicities
(2,3), gcd one and trivial pi1. Its prior certificate is reused without
reconstruction or replay. No plane chart or etale plane target was
claimed for that control.

## 3. The same divisibility test for a complete-base ruling

Now suppose rho:V->P1 under the SAME hypotheses. Euler integration
gives 2=2+sum_i(s_i-1), so EVERY fibre is irreducible. The retained
orbifold presentation and classification from PF132 imply that finite
cyclic pi1 allows at most TWO genuine multiple fibres. Three cone
points with finite group give one of the noncyclic spherical triangle
groups; four or more do not give a finite group. This reuses PF132's
classical input, rather than adopting another classification.

There is at least one multiple fibre: otherwise HC4MYML0 would give
the forbidden open A1-bundle over P1. Represent a single multiple
fibre by orders (m1,m2)=(1,m2), m2>1, placing the ordinary marked
fibre at0 and the multiple fibre at infinity. For two multiple fibres
put their orders m1,m2>1 at0,infinity. Set

    g=gcd(m1,m2),  (a,b)=(m1/g,m2/g).

The presentation gives pi1(V)=C_g. Normalize the base change
t=u^g, P1_u->P1_t, branched at0,infinity. For g=1 this is the
identity. The extension is connected of degree g since the generic
field of V is C(t)(z). At either marked fibre the local normalized
root extension has ramification index g/gcd(g,m_i)=1. Away from them
the base change is etale. Consequently the normalization W->V is
unramified in codimension one and is finite ETALE by purity over
smooth V, including at codimension-two points. It has degree |pi1(V)|
and hence is the universal cover.

The induced ruling W->P1_u has ordinary reduced A1 fibres. Each
marked original A1 lifts to g disjoint A1s, because a finite etale
cover of A1 is trivial. Their multiplicities are respectively m1/g
and m2/g: pull back the fibre divisor using t=u^g. Thus the pullback
fibres may SPLIT; they are not asserted to remain irreducible.

Equal orders m1=m2>1 would make both normalized orders one. All
fibres would then have a unit component, contrary to HC4MYML0 applied
to the pulled-back exact form on W. Equal orders are therefore
impossible here. If one original order divides the other, equivalently

    min(a,b)=1,                                                  (2)

W has at most ONE bad fibre. Section1 gives an OPEN A2 in W. This
includes the single-multiple-fibre case g=1, where V itself has the
plane chart. If (2) fails, the normalized pair is again coprime with
a,b>=2. The construction then supplies no plane chart; it does NOT
rule out charts obtained otherwise. No automatic divisibility of
the two orders is proved.

The cover and chart assertions use only the same normal-root,
purity, finite-etale comparison and bundle inputs already retained
by PF130/PF132. The above paragraphs check their hypotheses anew.

## 4. Minimal generic degree forces the universal cyclic relative cover

Return to an ACTUAL complex polynomial plane Keller map, its full
finite normalization N, and a genuine intermediate N->M->Y=A2.
Write d=bt, b=deg(M/Y)>1 and t=deg(N/M). Let V be the ENTIRE
etale locus of M->Y, not a selected smaller open or an assumed plane.
Suppose chi(V)=2, and some A1-ruling passes (1) or (2). Then its
universal finite etale cover W->V, of degree m=|pi1(V)|, has a plane
chart. PF130 Section2 gives an ACTUAL polynomial plane Keller map

    A2 -> W -> V -> Y,  generic degree b*m,  m divides t.          (3)

The open chart need NOT surject or be almost surjective onto V.
The degree in (3) is nevertheless b*m because an open immersion
does not change the function field. The target composition is etale,
so its Jacobian determinant is a nonzero constant on the source A2.

For an original counterexample of MINIMAL GENERIC DEGREE, m<t is
impossible: (3) would be a smaller counterexample, of degree b*m>1.
PF130's finite lift and normality argument then gives

    m=t,  h^-1(V)=W=U_N,
    h|U_N universal finite etale with group C_t,
    Delta=Eta_res=0,  b0=2*t-1.                                  (4)

ALL omitted primes and residual ramification were retained before
this deduction. Equation (4) proves their relevant residual defect
vanishes under the chart premise; it was not discarded in advance.
This argument also works for t=1 and gives no contradiction there.

Thus, if an Euler-two intermediate of a minimal counterexample does
NOT reduce to (4), EVERY one of its A1-rulings must fail the comparable
pair test: either an affine-base two-component fibre or a complete-base
two-multiple-fibre profile has a normalized coprime pair a,b>=2.
This is a necessary obstruction to this construction, not an exhaustive
classification of surfaces without plane charts.

## 5. Full deck cyclicity after compression, through Euler five

There is also a general strengthening of PF130 Section4. Whenever
that plane-chart compression gives a relative subgroup C_t of the
FULL generic rational deck group D, let c=chi(V) and |D|=t*h.
HC4MYLA0's free-deck Euler formula gives h|c and

    chi(U_N/D)=c/h.

If D is noncyclic, h>=2 because its relative subgroup has order t
and is cyclic. PF132's strengthened noncyclic cost gives c/h>=3.
Therefore

    D noncyclic ==> c>=6.                                       (5)

For c<=5 the FULL D is cyclic, improving PF130's older cutoff c<=3.
This requires the preceding cyclic relative compression; it is not
a theorem about arbitrary intermediates of small Euler characteristic.

In the Euler-two case (4), |D| is t or2*t. Since |D| divides d=bt,
an odd quotient degree b forces |D|=t. A proper primitive block
quotient also gives D=C_t by PF130/PF128's retained normalizer
argument. Without those additional conditions the order2*t case
is retained; relative and full deck groups are not identified.
No generic deck action on the original plane chart is assumed.

## 6. Primary-source audit and the remaining obligation

Miyanishi's retained *Affine pseudo-planes and affine pseudo-coverings*,
Oberwolfach Report19/2005, printed1110--1113/PDF38--41, was read in
full for these definitions and statements. PDF39 and PDF40 give the
orbifold presentation and cyclic-fibre-space context; complete PDF40
was visually inspected. An affine pseudo-covering is etale AND almost
surjective (its missed locus has codimension at least two).

Our chart A2 subset W -> V can omit a WHOLE fibre. Its generic field
extension is Galois of degree g in Section3, but it need not be an
affine pseudo-covering. This avoids a false conflict with report
Lemma2.6(2). Lemma2.7 states a DIFFERENT, almost-surjective nonGalois
plane cover of degree m2 when m1 divides m2 strictly; we do not use
or newly adopt that statement, whose proof is not supplied there.
The primary report supplies context and a terminology check, not
an independent proof of the new chart lemma.

The concurrently registered PF133/PF134 were also read in full before
registration here. They prove general deck-index/group-container and
noncyclic universal-ruling restrictions, including a normalized mixed
pair at Euler three. They are credited as separate results; the new
positive complete-chart lemma and cyclic Euler-two divisibility test
here neither replace nor re-adopt those theorems. No predecessor proof
or checker is modified or executed.

The Cambridge primary abstract for Gurjar--Masuda--Miyanishi--Russell,
*Affine Lines on Affine Surfaces and the Makar--Limanov Invariant*,
DOI10.4153/CJM-2008-005-8, was screened. Its ML0 hypotheses supply no
automatic divisibility or additive action for the remaining rigid/ML1
branch. Its full paper was not read or adopted. The retained Q-factorial
and equivariant quotient papers were likewise screened without a new
import. There are ZERO new external theorem imports or primary caches.

The written proofs are unbounded. A NEW preregistered checker supplies
only exact arithmetic, transcription, source-seal and boundary controls;
it is not a formal proof or a finite census offered as a theorem.
No registered mathematical checker is rerun.

The unresolved normalized pairs include (2,3), with no upper bound on
either entry. More fundamentally, the cyclic alternative (4) still
requires an actual finite-target or full trace/marked-energy argument.
This result does not supply the missing target coframe, rule out every
complete ruling, force an additive action, or prove HC4/JC2. The
published degree-four PLANE exclusion remains retained, without
treating a quartic intermediate as a plane. Review/outreach is DEFERRED.
