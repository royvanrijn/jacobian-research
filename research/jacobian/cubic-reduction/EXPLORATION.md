# Retained cubic-reduction discovery attempts

These are discovery experiments, not optimum certificates.

1. **Balanced monomial splitting.**
   [search.py](search.py), [log](search.log) and
   [checkpoint](search-checkpoint.json) retain a beam of width 40 through
   dimension 19. Each move adds carriers for a factorization of one
   maximum-degree monomial (one carrier for a square), applies a tame
   source/target gadget and clears standalone carrier monomials by target
   row operations. This run reached degree four, not three. Beam pruning,
   restricting moves to maximum-degree monomials, and the dimension cap
   leave other reductions untested. It is not a negative dimension theorem.

2. **Single coordinate-graph slices of the external 19-variable map.**
   [slice_external.py](slice_external.py) and
   [single-slices.json](single-slices.json) test the sixteen displayed
   carrier equations individually. Each straightforward slice raises
   degree above three. This does not exclude linear combinations of
   carriers. Those combinations produce the verified 17-variable map.

3. **Quadratic target corrections after selected single slices of the
   17-variable map.**
   [repair_slices.py](repair_slices.py),
   [log](repair-slices.log) and
   [case data](quadratic-slice-repair.json) test the eight quadratic
   carrier slices that leave the remaining carrier outputs cubic.
   None of their top three outputs can all be made cubic by quadratic
   polynomials in the remaining carrier outputs in this test.
   Other slices, simultaneous slices, higher-degree target corrections,
   general directions and other constructions remain outside the result.
   These experimental linear systems are not used as premises of the
   dimension-17 existence theorem or as a global lower bound.

The successful step was to inspect repeated carrier directions before
performing more product splitting. The external top coordinates use
\(w_{11},w_{12}\) only through \(6w_{11}-w_{12}\), and
\(w_{14},w_{15}\) only through \(3w_{14}-w_{15}\). The two complementary
directions are absent from every retained output, allowing exact
destabilization.

To improve the cubic dimension upper bound, change the carrier realization
or use transformations beyond these tested slices. Repeating the failed
beam with the same parameters or a single displayed coordinate slice
does not test that larger question.


## Historical reused-output construction

The [successor proof](PROOF13.md) records a beam search that can reuse
existing output coordinates in elementary target products. Ten new
carriers suffice, producing a degree-three map in13 variables with47
nonlinear terms. One exact transverse shear reduces the count to42.
Independent complete coefficient and determinant replays verify existence.

The wider order experiment and selected graph slices do not find a
smaller witness within their stated finite search. Their scopes are
recorded in the successor proof; they do not establish dimension13 or
support42 as minima.

## Corrected comparison and actual sparsity progress

The13-variable42-term map, like the17-variable58-term map, is dominated by
BCRDB1's pre-existing12-variable40-term cubic. Their exact verifications
did not establish a better size bound. See [the correction](BASELINE_CORRECTION.md).

The [current36-term integral map](PROOF12.md) starts from that established
twelve-variable baseline. Rational elementary-line descent performs two
shears; diagonal conjugation clears denominators. Direct independent
polynomial checks verify the40-term input and36-term output.

The completed graph tests on both13-variable frames remain historical
local misses: [47-term frame](quadratic-slice-repair13.json) and
[42-term frame](quadratic-slice-repair42.json). They do not motivate
repeating those dominated routes in place of the twelve-variable problem.


## Nine fixed graph slices: arbitrary target degree is now closed

[KCUB4](GRAPH_DELETION_OBSTRUCTION.md) independently checks finite dual
functionals that obstruct all polynomial target correction degrees for
one raw component in each of the nine fixed KCUB3 graph slices. Raising
the raw target degree alone is therefore not a next search.

Small initial quotients and the first capped adaptive run were
inconclusive in seven cases; the wider adaptive run terminates with all
nine explicit separators. The final result is certificate-based, not
an inference from those caps. A subsequent bounded cubic-pullback target
relation search finds no new first-stage direction involving the
obstructed output, but is retained only as exploratory evidence.

The separate shared-cube finite-space experiment tests a larger set of
polarization forms than the earlier termwise greedy choices. Its status
is recorded in power-search/checkpoint.json; no dimension improvement
is accepted without an exact reconstruction, determinant and collision.

## Shared-cube basis exchange: verified89-variable successor

The fixed pool has178 forms and cube-span rank121. The60-second SMT
support89 query terminated UNKNOWN by timeout; the1280-step plateau
basis search stayed at90. The successful search allowed temporary
support growth to94 and found an exact89-form identity, now converted
to an integral cubic-linear map with an independently checked
determinant argument and rational collision. See [the proof](POWER89.md).

Neither the finite pool nor the discovered89-form support is proved
optimal. Continue from the saved89-form identity rather than repeating
the older90-form baseline. The saved matrix and exact verification
do not require rerunning the time-limited discovery.

## Warm starts and weighted factor sums

The saved89 identity remains best after16,980 exchanges in the178-form
pool (60 seconds) and1,966 exchanges in an858-form weighted pool
(90 seconds). Both pools have cube-span rank121. These terminal bounded
misses do not establish any optimality. See the scoped
[derivative-space bound and next gate](FIXED_TENSOR_BOUND.md), which
independently proves34 as a lower bound for this fixed padded tensor
only. Its first Koszul diagnostic also gives34, not a stronger bound.

The36-term base's50-restart termwise decomposition gives91 forms.
The exact companion bridge shows that this base is in the previous
tensor's linear orbit, so it already inherits89 forms and the same
minimum. Its lower monomial count is not a new companion mechanism.
See [the bound note](FIXED_TENSOR_BOUND.md) for the explicit change
and the retained primary identity checks.

## Unrestricted single-cube replacement and changed search spaces

[The single-cube obstruction](SINGLE_CUBE_OBSTRUCTION.md) independently
checks86 explicit unit identities across15 source-support cases. No
shorter representation can retain a subset of the saved89 forms and
introduce just one new form, even with arbitrary coefficients over a
characteristic-zero extension. This supersedes the pair-line restriction,
not the open question of simultaneous replacements or smaller maps.

The quadratic-pencil pool has256 forms and rank154;10,554 bounded
exchanges stay at89. The nonmonotone ordinary-cubic beam stays at36
after69,747 line values; primary coordinate translation lines also
stay at36. Their exact caps and replay inputs are linked in the proof
note. These misses are not global optimum certificates.

## Fixed-span optimality and the two-new-form test

[The fixed-span theorem](SPAN_OPTIMALITY.md) classifies all cube points
in the saved89-dimensional scalar span and proves that every
representation staying there needs89 forms. Its factor trees and
annihilator premises are independently checked. This does not prove a
global minimum for the tensor or for Keller maps.

A possible two-new-form improvement must therefore use proportional
output columns and cancel two nonzero scalar residues outside the old
span. The930-form union pool has1774 such admissible pairs; exact
independent replay cancels at most two original forms, hence gives no
shorter witness. Unrestricted forms and larger simultaneous changes
remain open. The65 isolated points, eight binary families and exact
next gate are recorded in the canonical proof and receipts.

## Arbitrary two-form closure and raw13 coefficient rigidity

[The complete two-new-form obstruction](TWO_CUBE_OBSTRUCTION.md) uses
a mixed-second-derivative lemma to cover arbitrary source supports,
then independently checks607 unit identities across120 clique pairs.
It strengthens the finite930-form result. Three-form pool and tangent
searches are retained as scoped misses, not an arbitrary three-form
exclusion.

[The raw13-pattern classification](RAW13_PATTERN_RIGIDITY.md) records
the pre-existing three-variable degree-seven map's13 raw nonlinear
occurrences and proves rigidity in its exact coefficient pattern and
specified collision family. Seven initial boundary timeouts are not
negative certificates; diagonal normalization and15 independently
replayed identities complete that family by a different method. This
is not a global raw term bound or a normalized-axis witness.

## Weight-compatible additions beyond the raw13 pattern

[The19-atom tail classification](WEIGHTED_TAIL_CLASSIFICATION.md) allows
all compatible degree-seven tail atoms for weights(1,-1,-2), with the
first output fixed. Three exact branch identities plus an elementary
polynomial-degree argument force all six additions to zero.

[The marked degree theorem](WEIGHTED_MARKED_DEGREE.md) allows new
first-output terms through degree six. A cleared parameter certificate
forces theta^7=0; the certified zero slice requires a degree-seven atom.
Thus degree seven is sharp only in this marked weighted class.
The original23-atom degree-seven raw12-term query retains UNKNOWN.
The subsequent KCUB13 proof below settles its mathematical class; global
minimum questions and other collision/weight classes remain open.


## Full23-atom marked class closed (KCUB13)

The previously proposed degree7 first-output parameter-field gate is now
complete. The full field gives a unit but its denominator has an unwieldy
additional factor; generic success alone was not used as an exclusion.
A structural boundary argument replaces that route. For
f=1-3u-w+R*u^3+P*u^2+Q*u*w+T*w^2, R=P-4Q+16T, a Keller completion
requires f=0 smooth. IfR*T!=0, its smooth compactification has one point
at infinity and nowhere-zero differential omega. The identity
d(h/g^2)=-omega/g^3 forces finite pole degree2m but a zero at infinity
of order3m+1, impossible. Singular boundaries are excluded directly by
the determinant. ThusR*T=0.

Four cleared parameter identities then reduce the residual branches to
P=Q=T=0, with a fifth exact unit at the sole nonzero parabolic exception
T=7/4. The unaugmented cubicgraph and parabolic field queries timed out
at20 seconds; necessary polynomial-boundary Wronskian degree constraints
made both branches tractable. One polynomial renderer failed after linear
elimination introduced1/T; its error and source snapshot are retained,
and the corrected renderer clears the numerator. No timeout became UNSAT.
The primary full-field certificate is retained but is not part of the proof.

The independent Fraction checker reconstructs all23 source atoms, the full
universal3D determinant and linear parameterization, both boundary charts,
Wronskians and necessary equations, and multiplies all1834 certificate
monomials (maximum multiplier degree15). The curve/pole-order argument is
written proof, not formal machine verification. See
[the canonical proof](WEIGHTED_FULL_DEGREE7_CLASSIFICATION.md) and
[validation](weighted-newfirst7/validation.json).

The unique marked weighted degree<=7 map is the established raw13-term
representative; all ten extra atoms vanish. This is a matching raw13-term
minimum only in this class. The old30-second real query remains UNKNOWN
as a solver artifact, while its mathematical class is settled byKCUB13.
The next raw sparsity gate must change weights, first-output coefficients
or collision marks, for example freebeta. Ordinary cubic12/36,
cubic-linear89 and normalized3D[10,48] remain the separate current bounds;
the full global simplicity goal is still open.


## Collision-free weighted sector closed (KCUB14)

Freeing beta was superseded by a stronger experiment: remove every
collision equation and keepR,P,Q,T independent in the full23-atom
weighted degree7 determinant model. The smooth-cubic pole lemma still
forcesR*T=0. Generic graph/conic/cubicgraph lifts clear toQ^16,
T^2*(P*T-Q^2/4)^7 andP*R^5. Direct polynomial boundary Wronskians
supply complete linear reductions. Full-ring saturation units then
close the nonzero purecubic, quadraticgraph and nonparallel parabola
branches, with no omitted generic-denominator roots.

The parallel-line branch has f=phi(L),L=3u+w,phi=1-L+T*L^2.
Its seven linear boundary equations makeg=A(L)u+B(L) and
h=C(L)u^2+D(L)u+E(L), degA<=2,degC<=1. The determinant first forces
C=c*A^2. A target shear removesC, and the next coefficient forces
A=alpha*phi,D=-3. The remaining equation integrates to
phi*(alpha*phi*E-D*B)=L, impossible. This written algebra proof closes
that branch. Generic parallel postprocessing was explicitly stopped
after its complete Singular export became redundant. Direct and
simplified full-ring lifts and a combined liftstd attempt remain
30-second timeouts; none is used as a negative certificate.

For fixedf=1-3u-w the full determinant-only tail solves before/after
linear reduction each time out20 seconds. Splitting its boundary
G0=A*u+B atA1 orB0 instead gives two exact unit identities. Off those
exceptions beta=B/(A-1) is a nonzero fixed point and automatically
supplies precisely the collision covered byKCUB11. Thus the normalized
map is unique. No collision has been assumed or numerically sampled.
The raw13-support-only primary ideal is also zero-dimensional with the
known vector; it is retained but not used as a new witness or proof.

The independent Fraction audit passes all eight used certificates,
30337 multiplier monomials and maximum multiplier degree15. It
reconstructs the full source-variable determinant without collision
conditions and certifies complete boundary row reductions with only
explicitly invertible parameter monomials. Its first fixed pivot order
could not certify one parabola row; selecting any remaining allowed
monomial pivot resolves that audit issue without inverting a new
parameter polynomial. Initial replay logs are retained. The written
curve and parallel-line degree/divisibility arguments are not formalized.

[KCUB14 proof](WEIGHTED_UNMARKED_DEGREE7_CLASSIFICATION.md) and
[validation](weighted-unmarked7/validation.json). Diagonal conjugacy
removes the two fixed first coefficients whenever both are nonzero,
preserving all13 terms. Consequently changing collision marks or varying
these coefficients among nonzero values cannot improve raw sparsity in
this degree7 weight sector. The next gate is a zero firsta^2b or a^3c
coefficient, different weights/coordinates, or another construction.
Global simplicity minima remain open; ordinary cubic12/36,
cubic-linear89 and normalized3D[10,48] bounds have not changed.


## Full weighted degree-seven sector closed (KCUB15)

The coefficient-zero gate is complete. All three normalization cases
(mu,nu)=(0,-1),(-3,0),(0,0) retain21 free coefficient parameters and no
collision conditions. The smooth-cubic boundary argument extends to
arbitrarymu,nu by the completed square. Nine cleared field identities
and four units close the residual parameter branches; they have12342
independently replayed multiplier monomials, maximum multiplier degree15.
The mu-zero case has no Keller completion. The nu-zero case reduces to
f=1-3u; its direct tail solve times out20 seconds, but a high-w split
uses a13-monomial unit for the sole remaining nonzero quadratic tail.
A normalized auxiliary high-w calculation is retained but not needed.
The generic nu-zero parabola lift has an extraneous factor; its reduced
full-ring saturation unit covers every nonzeroT instead of sampling roots.

Two written univariate boundary lemmas replace further generic solves.
Their global linear/quadratic transverse tails give ACprime-2AprimeC=0;
a target shear removes the quadratic term. The remaining coefficient
relation and constant determinant force a nonconstantphi to divide its
normal variable despitephi(0)=1. Both complete source-degree/linear
reductions and the coefficient identities are checked independently.
The initial abstract coefficient replay rejected a wrong sign in a term
that vanishes after the shear; its corrected full identity passes, and
the original log remains retained. It does not change any solver status.

The mu=nu0 remainder is f=1. Then quotient(g,h) is a plane Keller map
of ordinary degree<=4. Nakai-Baba1977, Theorem3(1)-(2), printedp408 is an
explicit primary import, not a reproduced theorem. The original PDF and
checksum/page record are retained, and the complete theorem page was
visually inspected with bundled Poppler. The plane inverse(U,V) lifts
polynomially to(a,a^-1U(aB,a^2C),a^-2V(aB,a^2C)): zero constant terms and
no first-target linear term inV exclude negative powers. This separates
automorphism families from counterexamples in the determinant-only model.
A literature degree-four lead was checked against Gao's actual text:
it denotes geometric fibre degree, not an improved ordinary degree.

[Complete theorem](WEIGHTED_COMPLETE_DEGREE7_CLASSIFICATION.md) and
[validation](weighted-zero7/validation.json). Together withKCUB14, this
proves the full weight(1,-1,-2) sector through degree7 has exactly13 terms
in every counterexample, with weighted minimum ordinary degree7. The
13-term sparsity minimum has the stated degree cutoff. No global minimum
or smaller construction is claimed. Next inspect prior work and change
weights, for example(1,-1,-3), seek fewer terms above degree7, change
polynomial coordinates or return to dimension/special-form construction.
Ordinary cubic12/36,cubic-linear89 and normalized3D[10,48] remain separate
unchanged bounds. The global simplicity objective remains active.


## Sparse new-weight search and a complete13-slot exclusion (KCUB16)

The new bounded models use weights(1,-1,-3),(2,-1,-3),(1,-2,-3),
ordinary degree<=7, identity linear part and every compatible atom.
Their nonlinear slot counts are20,16,13 and free parameter counts19,16,13.
All unrestricted coefficient bases time out20 seconds; all raw<=12,
nonlinear-first real queries returnUNKNOWN30 seconds. The source and
receipts are retained, and no process is left running.

A general one-positive-weight lemma makes explicit collision variables
unnecessary for these construction searches: nonlinearF1=a*f gives a
zero offa0, and the negative-weight restriction ata0 is a triangular
onto map. This is a sufficient collision gate in any dimension/degree;
it does not say every counterexample must have nonlinearF1. The linear
first-component sector must be handled separately in a full classification.
The independent audit reconstructs all three source determinants, their
complete linear eliminations and all negative-plane restrictions.

The13-slot(1,-2,-3) model was prioritized because deleting its highest
first-output atom already meets the12-term target. Diagonal normalization
and the degree-bounded(2,3)-Wronskian close every first-coefficient case.
Two cleared identities leaveP0 and the rational exception-810/2021;
four separate unit identities cover the required zero/root slices.
A written quadratic univariate primitive-divisibility argument and the
explicit plane degree<=3 consequence of the pinned Nakai-Baba theorem
finish the remaining cases. The weighted plane inverse has no negative
a powers. Every map in this complete degree7 weight sector is therefore
an automorphism, not a smaller counterexample.

[Proof](WEIGHTED_123_DEGREE7_EXCLUSION.md),
[validation](weighted-sparse7/1_-2_-3/validation.json). Six independently
multiplied identities have9017 multiplier monomials,maxdegree13.
The balanced boundary controlG=t^2+1,H=t^3/3+t/2 has(2,3)-Wronskian1,
so the old(1,2) degree<=1,<=2 lemma must not be imported into this model.
All exceptional roots are resolved exactly. The written gate, degree and
integral proofs and imported plane theorem are not formalized.

Only(1,-2,-3) is closed by the later proof. For the two open weights,
next change to complete diagonal first-coefficient cases or structural
boundary/quotient equations before another bounded calculation. The mixed
weight invariant ring may have an algebraic relation and diagonal normalizations
may require algebraic extensions. Do not re-run the unchanged20/30 caps.
No construction bound improves: raw13,ordinary cubic12/36,cubic-linear89
and normalized3D[10,48] remain separate. Global minima remain open.


## Mixed-weight paired-infinity obstruction closes(2,-1,-3) (KCUB17)

Changed first-coefficient field models and eight fixed low-support/rank
controls all give primary unit ideals quickly. The controls are construction
attempts, not an exhaustive parameter proof. The generic full-field lift
has an additional factor and is retained but unused. The useful route
instead uses the nonzero-a covera=s^2,b=u/s,c=v/s^3. The quotient firstf is
even and the tails odd; the original source determinant independently
checks the1/2 and3/2 coefficients in the quotient identity.

For nonzeroP*T, completing the square gives a smooth even quartic with
two infinity points. The sign involution exchanges them, forces equal
pole orderm of oddg at both, and makesrho=h/g^3 invariant. The differential
identityd(rho)=-2*omega/g^4 forces all2m finiteg zeros simple andrho pole
degree6m. Its two infinity values agree, so subtracting that value gives
zero degree at least8m+2, impossible. This written curve/divisor proof
forcesP*T0; it is not an argument from generic coefficient sampling.

Two cleared identities then force the conic discriminant zero onP0 or
cross coefficientnu0 onT0. Only69 multiplier monomials,maxdegree7, are
used. Complete global transverse degree constraints and parity turn the
residual cases into a written univariate polynomial-operator contradiction.
All source parity/degree conditions and linear parameterizations, the
completed quartic and coefficient identities are independently checked.
The constant-first family has a generic plane inverse overK(a) by the
explicit Nakai-Baba prime/degree4 import up through degree5. Its unimodular
linear matrix overK[a] makes the unique formal inverse have polynomial
coefficient functions, so the generic inverse is globally polynomial.

[Canonical proof](WEIGHTED_213_DEGREE7_EXCLUSION.md) and
[validation](weighted-213/validation.json). The full16-slot degree7 weight
class consists of automorphisms. The original20-second basis timeout and
30-second realUNKNOWN remain unchanged as experiments. Only the separate
(1,-1,-3) class is still open among the initial new-weight models. Next use
its normalized first-output parameters and structural boundary equations,
not its unchanged dense queries. No construction bound improves; the full
global simplicity goal remains active. Curve theory, the parity argument
and the plane theorem are written/imported rather than formalized.


## Complete20-slot (1,-1,-3) exclusion closes the last initial model (KCUB18)

The preceding open-weight statements describe their earlier checkpoints.
The new determinant-only model restores all20 original coefficients before
first-parameter normalization. Its generic cleared identity forces the first
cross coefficientQ=0 exactly (targetQ^13,229 multiplier monomials). The original
unrestricted coefficient/real support queries remain timeout20/UNKNOWN30.
Changed qzero/both/mu-zero/both-linear field queries time out20 as recorded;
these failures are not mathematical exclusions.

When the first quotient depends onw, normalize its coefficient to-1 and
mu to-1or0. The boundaryg hprime-3gprime h=1 gives degG<=2,H<=6 for cubic
graphs, orG<=1,H<=3 for quadratic/linear graphs. Append all boundary and
native affine-in-tail equations; the complete rank6 reductions leave9 tails
and onlyR/P monomial denominators with separately covered zero strata.
Five full-ring unit identities close all but the nonzero-mu quadratic graph.
Its reduced, parameter-field and factor-cleaned queries still time out20.
They are preserved, with their actual input and status.

The useful change is a written highest-transverse-degree case split. In
v=f coordinates the equation is(G+vGv)Hu-Gu(3H+vHv)=1. Quadratic tails
force either a contradiction or an affineG, excluded by a degree-free
polynomial eigenvector argument at the source. The remaining affineG1
case reduces to five polynomials inT,P,beta; their unit identity has only
eight multiplier monomials,maxdegree2. The other nu0 branches are excluded
by squarefree boundary roots, polynomial degrees and antiderivative
divisibility. The constant-first family uses the precisely attributed
Nakai-Baba prime/degree4 theorem through plane degree5, followed by the
unimodular formal-inverse argument overK[a].

[Canonical proof](WEIGHTED_113_DEGREE7_EXCLUSION.md),
[validation](weighted-113/validation.json). Independent Fraction replay
reconstructs all20 atoms, the full determinant, all6 complete linear
parameterizations, source bounds, coefficient identities, final core source
conditions and all7 exact identities:5728 multiplier monomials,maxdegree13.
Written case coverage/divisibility and the imported plane theorem are not
formalized or reproved by the arithmetic replay.

All three additional initial weight classes are now automorphic through
degree7, while the separate(1,-1,-2) class retains its established13-term
counterexample minimum through that degree. Do not repeat unchanged closed
models. Next change weights, degree, polynomial coordinates or construction;
ordinary cubic below12/36, cubic-linear below89 and normalized3D below48
remain independent open avenues. No construction bound improves and the
full global simplicity objective remains active.


## Cubic carrier realization: exact four-parameter family (KCUB19)

The new construction attempt returned to KCUB3 instead of extending the
closed three-dimensional weight sectors. All792 elementary monomial
quadratic conjugacy lines have only zero cubic parameter. A direct written
leading-degree argument strengthens this: every arbitrary nonlinear
one-axis source shear leaves an unchanged output of degree>3, even after
any correction confined to that same target row. This is not a bound on
simultaneous nonlinear transformations or other frames.

The displayed carriers have12 nonlinear atoms (the initial preparation
count13 was corrected before solving). Varying all12 coefficients gives
14 high equations and a two-parameter diagonal family, always36 terms.
Adding all10 missing weight-compatible base-variable monomials through
degree3 gives22 atoms and forces every new coefficient zero. The next
model adds all weight-compatible quadratic monomials in preceding carriers,
for47 total coefficients. The source change is constructed by exact
triangular recursionC_old(E)=C_new; thus determinant and collision are
inherited fromKCUB3, and only cubic-degree conditions need solving.

The47-variable primary basis is four-dimensional. Rather than treating
that basis as proof, the exporter saves43 necessary ideal memberships with
123 multiplier monomials,maxdegree3. Independent Fraction arithmetic
reconstructs all47 atoms, the complete source change/composition,93 high
coefficient equations,96 low slots, all memberships and the complete family.
Product identities forcep,q!=0 before division; the remaining parameters
s,t are free. Substitution verifies every high coefficient vanishes.
Generic support45 factors into33 units,4s,5t,one(3p+2s),two(pq+3t), giving

```
36+4[s!=0]+5[t!=0]-[s=-3p/2]-2[t=-pq/3].
```

The exact minimum36 requires s=t0; the entire equality family is checked
to be diagonal conjugates of the old map. No37/38-term member occurs in
this fixed realization. [Canonical proof](CARRIER_REALIZATION_LIMIT.md)
and [validation](carrier-reparam-wide/validation.json). The bound does
not include arbitrary subsequent coordinate changes, different triangular
orders or carrier-dependent cubic monomials beyond the stated support.

Nine changed raw graph restrictions atp=q1 all give primary unit ideals,
without arbitrary target repairs or an independent all-degree graph proof.
Six changed(s,t) frames followed by at most6 greedy rational linear passes
start at39,39,42,40,41,45 and finish39,36,39,40,36,40. None improves36.
These completed finite experiments are retained separately from the theorem.

Two exporter/audit details were caught and fixed before a pass: Singular
matrix-to-scalar comparison does not test that every matrix entry is zero,
so the exporter checks each entry of its residual row individually. The
shared polynomial substitution helper handles nonnegative powers; applying
it to Laurent coefficients discarded negative p/q exponents and failed the
independent diagonal-conjugacy comparison. The final checker specializes
s=t0 by exact coefficient filtering, preserves all p/q exponents, and asserts
nonnegative powers before its polynomial substitutions. The shared helper
and older theorem sources were not changed.

No positive bound improves. Next change carrier-dependent cubic atoms,
triangular order/graph or coordinated source/target transformations, or use
another construction. Do not repeat a higher-degree one-axis shear or the
closed47-coefficient realization as if it were a new search. Ordinary
cubic12/36, cubic-linear89, raw3D13 and normalized3D[10,48] stay distinct;
the full global simplicity objective remains active.


The subsequent read-only substitution-domain audit covered exactly seven
current weighted/carrier checkers, intercepting their output writes. It
found the same helper-domain crossing in the previous KCUB18 core source
normalization: its beta-denominator coefficients were incorrectly passed
as polynomial input. `source_pair_substitution` now changes onlyu,v and
retains every parameter exponent. The entire KCUB18 checker passes with
those rational source conditions preserved, and the new domain audit
passes all seven modules. The first failing audit and old KCUB18 receipt
remain in `substitution-domain-before-fix.json` and
`weighted-113/validation-before-laurent-fix.json`. Only KCUB18's intended
receipt/source seal was refreshed; the other older receipts were not
written. The shared helper was unchanged. This is a checker correction,
not a new mathematical exclusion or a change in construction bounds.


## Complete weighted cubic carriers and companion rank (KCUB20)

The88 missing carrier-dependent cubic atoms extend KCUB19's47-atom model
to all135 degree2/3 weighted triangular atoms in the original order. Exact
composition gives540 high equations and184 low slots. The primary basis
has dimension8. A proposed explicit family satisfies all high equations;
a127-target exact lift (410 multiplier monomials,maxdegree5), independent
composition/membership replay and written unit-justified elimination prove
that it is the complete family. The coefficient ideal is the two-unit,
six-free-variable Laurent polynomial ideal displayed in the canonical note.

After diagonal normalizationp=q1, all64 generic nonlinear coefficients
are affine in six parameters. The primary QF_LRA support<=35 query is
UNSAT in0.001seconds. Its solver status is not a proof premise. The
independently reconstructed factor inventory has31 units and16 affine
factor groups. A written indicator inequality proves sharp36 support,
only the old diagonal family at equality, and no37/38 terms.

Weight-block coefficient minors and a three-coefficient Bezout identity
prove nonlinear-output rank12 everywhere. Thus no affine frame supplies
an affine-linear component for deletion. This does not close nonlinear
graph deletion or global cubic dimension. Cubic-output rank has minimum5,
with the exact locus s=u=v=w=z0, t arbitrary. The entire locus is checked
linearly conjugate to KCUB3, including the t-dependent quadratic terms.
Consequently the standard n+r+1 companion has minimum18 here and gives
no new tensor at its minimum. A new tensor from this family has rank>=6,
hence at least19 companion variables in that recipe. It could still have
smaller shared-cube complexity; no global89-variable lower bound is made.

[Canonical proof](FULL_CUBIC_CARRIER_CLASSIFICATION.md) and
[validation](carrier-reparam-full/validation.json). The independent checker
uses guarded polynomial substitutions and explicit Laurent coefficient
units; it does not call the primary CAS/SMT tools. Mathematical elimination,
indicator counting and companion covariance remain written, not formalized.

Ten changed triangular orders were also prepared. The complete weighted
quadratic/cubic support is enumerated separately for each. All nine primary
identity-carrier charts and four additional same-weight carrier-combination
charts are units in every order; last-pivot/singleton charts are exactly
the identity queries. Two orders already have unit cubic-completion ideals.
The ten-case set comprises five adjacent swaps plus moving the zero-weight
carrier4 last, moving12 first, moving positive-weight6 last, reversing all
carriers, and putting factors6,7,8 first. Their model/checkpoint inputs and
actual primary statuses are retained in `carrier-order-search/`; this does
not exhaust all orders or all12-output combinations and is not an
independently certified global exclusion.

A joint rational orthogonal line preserving2*x2^2-x11*x7 andx2+x10 has
60 generic nonlinear occurrences. Its exact rational cancellation values
0,1/3,3/8,3/4,3,6 give36,54,59,59,59,59. This changed simultaneous frame
therefore gives no improvement on that line, without an arbitrary-frame
or algebraic-parameter optimum claim.

No positive construction bound changes. The original-order degree2/3
carrier realization is now complete. Next change base/frame/order or
construction, or build and decompose a genuinely new rank>=6 companion
from the family. In the fixed old tensor, at least three arbitrary new
forms are still needed to evade the KCUB9 exchange obstruction; a changed
companion is a different problem. Raw3D13, ordinary cubic12/36,
cubic-linear89 and normalized3D[10,48] remain separate current metrics.
The global simplicity objective remains active.


## Changed19-variable companion pencils (KCUB21)

The next construction attempt left the fixed old tensor. Two source-pencil
pairs preserve cubic degree by exact quartic cancellation:
E_s changesx5 bys*b*c/2 andU bys*b*R; E_w changesx8 byw*c*R andx12 byw*R^2.
Each is two commuting triangular quadratics, so determinant and collision
are inherited fromKCUB3. Their nonzero cubic-output rank6 gives a standard
19-variable companion. This is not a new size record over the known18D
homogeneous or89D cubic-linear examples.

The old89 cube identity was transported exactly and extended by polarization
cubes for the changed terms. The s-pool has277 forms/rank165; mixed product
bases enlarge it to427/rank171. The w-pool has271/rank163. Exact basis
exchanges optimize a common nonzero rational parameter by carrying both
affine coefficient matrices. Forty-second runs change104->98 (5313 steps),
98->98 (2040 steps), and101->98 (5821 steps). The best parameters are
s=-6 andw=-3/2. Every saved primary identity is checked exactly. The
mandatory pure-tau padding is retained throughout; this is part of the
construction class, not a claim about all decompositions.

Two separate LRA/cardinality queries in the277-form s-pool at bounds88 and97
both returnUNKNOWN after20 seconds. These are not UNSAT and close no pool.
The later SAT numeric guard fix did not rerun them; their executed source
is retained as `smt_companion_pencil19_before_sat_guard.py`, with its exact
hash named by the checkpoint. Generic/pencil rank and bounded scores are
not lower-bound proofs.

A separate exact derivative rank certificate uses only rational Laurent-unit
pivots in the parameter:37 pivots per pencil leave a zero residual, so rank
is37 for every nonzero parameter, with no hidden exceptional roots. At zero
an independent rational rank is35. The derivative span lower-bounds the
number of scalar cubes, giving the conditional bound37. Independent
Fraction replay reconstructs each full source composition and companion,
checks every pivot and cube coefficient, matrix rank19, pure-tau padding,
uniform100/101-cube affine identities and complete rational base/companion
collision points. Each pencil has a98-cube identity at its best parameter;
the least count over its nonzero parameters therefore lies in[37,98].
Neither98 is claimed minimal, and no dominated98-variable matrix was built
or advertised as an improvement.

[Canonical note](CHANGED_COMPANION_PENCILS.md),
[validation](companion-pencils19-validation.json). The next changed gate is
the coupledw,z rank-six component. Normalize its extra cubic output toa*c*R,
use embeddingw*e10+z*e11, and beginz=-4/3 withfreew so that the cubic R-part
has the repeated factor(3b-2S)^2/3. This is a proposed changed model, not an
executed search or a theorem about lower cube counts. Do not rerun the
unchanged three pools or their UNKNOWN queries. All global simplicity
metrics remain open at the existing upper bounds.


## Complete signed orthogonal block (KCUB24)

The [new subgroup theorem](ORTHOGONAL_BLOCK_MINIMUM.md) extends the old
one-parameter line to both signs of the entire orthogonal block, with the
open and boundary charts included. Independent Fraction reconstruction
proves minimum36 and next count47; no smaller positive map is found.
Further affine sparsification must leave this block group.

The first29-atom combined base-axis calculation in `base-quadratic-conjugacy/`
was redundant with KCUB19's all-degree one-axis theorem. Its original
Singular renderer error, source snapshot and corrected primary bases are
retained, without a new mathematical claim. Read the full canonical prior
scope before preparing a model, not only the short next-step description.


## Coupled affine centres and stable conjugacy obstruction (KCUB25)

The97 new coupled carrier-zero graph centres in `joint-translation12/`
all produce denser inverse-derivative-normalized cubics, with best145.
All counts and the four retained full maps/collisions are independently
reconstructed. This is a finite miss; arbitrary translations or subsequent
linear frames are not excluded. Do not repeat the same grid.

The separate [formal invariant proof](STABLE_CONJUGACY_OBSTRUCTION.md)
closes every stable polynomial-conjugacy attempt to remove an identity
factor from the saved12-variable map, in arbitrary degree. It uses a
rank11 quadratic minor and a cubic functional that kills the entire
quadratic Lie image. Unequal source/target changes are outside this
obstruction and remain a valid dimension-reduction route. No smaller
positive map or global cubic dimension lower bound is claimed.


## Target rigidity through degree5 and uncapped single-coordinate changes (KCUB26)

The [new target theorem](TARGET_CUBIC_RIGIDITY.md) certifies full rank6175
for all nonlinear target monomials throughdegree5 whose pullbacks would
need to remain cubic. It applies after arbitrary affine source/target
frames. The degree4/5 extension completes in its60-second primary budget;
there are no unresolved modular dependencies. An independent row-minor
replay verifies the exact characteristic-zero exclusion.

A separate polynomial inverse-Jacobian calculation gives twelve rank443
operators on the455-dimensional cubic source space. Their kernels are
exactly the affine spans of the other eleven outputs. Thus even rational
single-coordinate target changes of arbitrary degree cannot preserve
cubic degree in the displayed axes. No smaller positive map was found.
Future work must use the explicitly remaining coordinated changes,
temporary degree growth, higher target degree after a changed output frame,
other graph directions or another construction. This is not a universal
closure of all nonlinear left-right changes.


## Uniform weighted graph completion and partial carrier frontiers (15 November2026)

The [weighted graph proof](WEIGHTED_GRAPH_COMPLETION.md) replaces finite
parameter samples by two universal polynomial dual identities. The
[primary derivation search](graph_pencil_derivations12.py) retained its
order5 false positive, failed bounded dual supports, order7 models and the
preliminary parameter factor t-3/8. The [adaptive support search](graph_pencil_adaptive12.py)
adds the omitted-root certificate; the independent checker verifies a unit
pairing valid at every parameter. All nine same-weight projective charts
are now excluded for the specified one-stage polynomial repair. This does
not close cross-weight/nonlinear graphs or ordered stages.

The [complete subset experiment](carrier_subset_frontier12.py) constructs
all512 literal carrier-zero restrictions. Its independent sequential
replay and full determinant checks establish the
[finite-family frontier](CARRIER_SUBSET_FRONTIER.md), including intermediate
7-variable degree4/23-occurrence,5-variable degree5/17-occurrence and
4-variable degree6/15-occurrence maps. These are explicit witnesses;
no global-record or optimum claim is made. The
[three-baseline primary comparison](carrier-subset-comparison/checkpoint.json)
is separately retained without independent replay, and does not test the
40-term BCRDB1 conjugate or arbitrary carrier combinations. No repeated
unchanged subset census is needed.


## Complete linear carrier-subspace charts (15 November2026)

The weighted-only first screen used84 charts for6-variable quartic maps,
126 for7-variable quartics, nine for4-variable sextics and36 for5-variable
quintics. Every excluded weighted chart had a fixed high coefficient;
the surviving weighted parameters left the retained maps unchanged.
These checkpoints are preserved underweighted-subspaces-nXX-dY.

The successor [full unweighted producer](linear_carrier_subspaces12.py)
allows every lower free carrier coefficient in each highest-pivot relation.
Its ten completed checkpoints cover all512 chart families in dimensions
3..12. Most degree obstructions remain single fixed coefficients; two
6-variable and seven7-variable charts need constant combinations of
several coefficient equations. The original parameter-equations-unresolved
statuses remain in the primary screens; their actual resolution is the
[exact certificate file](linear-carrier-subspaces12/certificates.json).
No timeout, generic Groebner basis or parameter sample is a premise.

The [independent replay](verify_linear_carrier_subspaces12.py) verifies
all508 unit exclusions, the three proper surviving normal forms and both
quartic branches, all generic source embeddings/collisions and all ten
old upper witnesses. Counting fixed nonlinear coefficient slots plus one
exceptional two-coefficient unit proves sharp13 occurrences throughout
the entire raw linear-carrier graph class. The
[canonical theorem](LINEAR_CARRIER_SUBSPACE_LIMITS.md) records the fixed
graph parametrization and output-frame boundaries. This is not a global
sparsity or dimension result and does not close subsequent target repair.


### Completed output-frame audit after KCUB29

The candidate simultaneous row changes2-3*row7-12*row10,
3-4*row9-20*row8 and9+3*row8 give a verified integral cubic with31
nonlinear and17 linear terms. All48 total terms remain. The
[canonical proof](OUTPUT_FRAME_TRADEOFF.md) gives the exact six-point
output-only Pareto frontier; no identity-linear-part improvement is claimed.
The [independent checker](verify_output31.py) reconstructs all coefficients,
the full symbolic determinant and collision, all42 coefficient columns and
support-bound premises. A source-normalized identity-linear-part conjugate
has74 nonlinear terms, so this candidate does not improve36 under that
normalization. All data and intermediate stages are retained in the
[full witness](output-frame31/witness.json). General source/output changes
and output frames of other graph maps remain outside this calculation.


### Completed five-parameter joint source/output family

The proposed unitriangular source matrix on(2,7),(2,10),(3,9),(3,8),(9,8)
is now completely reconstructed. The convenient parameterization
S_(3,8)=z+w*t makes its inverse source3=x3-w*x9-z*x8 polynomially simple
and covers every original parameter. The
[full model](joint-unipotent5/model.json) has90 nonlinear coefficient
slots and19 factors; the source-only form has68 nonlinear source columns.

The [independent checker](verify_joint_unipotent5.py) verifies every
coefficient/factor and all68 source-vector quotient units, all generic
collisions, complete branch indicator identities and two matching controls'
full symbolic determinants. The [canonical theorem](JOINT_UNIPOTENT_LIMITS.md)
proves36 only at the identity normalized frame and next39; any changed
source needs at least34 terms even with arbitrary output mixing, while31
remains at the unchanged source. These are exact family limits, not global
or total-term optima. No sampled grid or generic-only result is used.
The old36 and74 points remain controls, not new progress claims.

The proposed weight(2,-1,-1) gate is now completed and extended to every
equal-negative ratio by [KCUB32](EQUAL_NEGATIVE_ALL_RATIOS.md). It excludes
all unequal primitive ratios fully throughdegree7 and imports only the
strict111 branch; no whole unrestricted profile or positive bound changes.


## Complete primitive equal-negative ratios (KCUB32)

The first211-only screen has ten correct primary normal charts. The next
[17-chart first pass](equal-negative-smallweights/README.md) contains three
off-locus quartic(1,2) templates; those models and certificates are retained
unused. The correct relation Q^2*S^5=constant*U^4 reduces that branch to
two line charts. The v2 producer asserts the leading coefficient is zero
on each template, and the independent checker reconstructs all valuation
possibilities and full determinants before multiplying the certificates.

The complete [new proof](EQUAL_NEGATIVE_ALL_RATIOS.md) uses the primitive
weight lattice to reduce all unequal ratios throughdegree7 to ten models.
It verifies147 atoms,30 generic identities and14 unit identities for145
generators with376 multiplier terms. No bounded parameter sample is a
proof premise. Ten non-affine triangular inverse controls pass. The full111
outside-strict sector and its older timeouts stay unresolved, and the
strict partition stays62 excluded/12 remaining.

The proposed weight(2,-1,-2) gate is now completed by KCUB33 below.


## Full primitive opposite-weight classification (KCUB33)

Prior-scope preflight found no existing closure of the complete212 family.
The boundary b=0 of its full20-atom model removes all u-only terms and
leaves11 coefficients. The discovery calculation first produced five
primary saturated unit certificates across212 and312. All inputs, scripts
and logs remain in `opposite-weight-reflection/`; the early212-only screen
is also retained. The producer's empty-support Python integer/substitution
error was fixed before the completed fourteen-model checkpoint.

A shorter proof supersedes those large units as a proof method, without
claiming their independent replay. The [three universal coefficients](opposite_weight_structure.py)
force every mixed uv coefficient to zero by the highest active support
identity, with seven complete nonzero support cases and symbolic m,q.
The residual Jacobian factors into univariate polynomials; its leading
integer coefficient forces all three multiplicative factors to one.

The [canonical theorem](OPPOSITE_WEIGHT_CLASSIFICATION.md) classifies every
primitive weight(m,-q,-m),m>=2 throughdegree7 as an explicit triangular
shear after diagonal normalization. The [independent checker](verify_opposite_weight_reflection.py)
reconstructs169 atoms,43 reduced atoms, fourteen complete boundary and
actual determinant dictionaries, three symbolic identities, seven support
masks, nine degree cases and fourteen two-sided inverse controls. A written
module degree bound handles all larger weight sums, so the result is not
a finite ratio scan. No external plane theorem is imported.

The m=1 exception remains outside the argument: its second-output module
has an additional generator, and after permutation it includes the known
weight112 counterexample. The proposed ratios(1,-1,-r),r=4..7 are now preflighted and completed
by KCUB34 below. No positive bound or full invariant-profile count changes;
all global minima remain open.


## Exceptional positive-unit/opposite-unit weights (KCUB34)

The previous goal turn made concrete progress as KCUB33; there was no
live solver to await. Full prior-scope preflight read KCUB15-18 and KCUB33,
the current index/registry and retained source/archive text. No r4..7
closure was found. The [new complete model producer](exceptional_weight1.py)
then enumerated15,15,11,11 atoms for weights(1,-1,-r),r4..7, with72
native determinant coefficient equations. The11-slot cases were direct
raw13 improvement targets, not support-pruned templates.

The [structural calculation](exceptional_weight1_structure.py) replaces
full coefficient solving by a univariate root/degree argument and a graph
Wronskian degree split. Only r5 retains a nonaffine graph resonance. A
preliminary nine-equation primary Groebner basis was[1], retained as
unused discovery data. Three constant pivots instead reduce this branch
to two univariate polynomials with a five-monomial rational unit. Two
more coefficient identities exclude the constant-first quadratic tails.

The initial independent replay stopped at a decoder/export mismatch:
pair-format primary polynomial rows did not match the existing named-field
decoder. Both producers were changed to the repository format and every
model regenerated before the complete passing replay. No failed run or
primary solver result is promoted to an exclusion.

The [accepted theorem](EXCEPTIONAL_UNIT_WEIGHT_CLASSIFICATION.md) now
classifies the entire r>=4 family through7 as compositions of two elementary
shears, including exact parameter product restrictions and all coefficient
boundaries. The [independent checker](verify_exceptional_weight1.py) verifies
52 source atoms,72 generators, three univariate and four graph identities,
28 vector-field eigenvalues, all three pivots and the five-term unit, both
quadratic-tail identities and eight generic two-sided inverse charts.
Written complete modules cover all r>=8, so this is not a bounded ratio
search. No external plane theorem is imported.

No smaller positive map is found and no unrestricted profile is closed.
The r2/r3 prior proofs and their original UNKNOWN/timeouts remain unchanged;
full11177* still needs a changed method. The proposed mixed weights(2,-1,-4),(2,-1,-5) are now preflighted and
completed, extended to every r>=4 by KCUB35 below.


## Complete mixed positive-weight2 family (KCUB35)

The previous turn made concrete progress as KCUB34. Current scope preflight
read the full213 proof and the33/34 boundaries, then scanned the registry
and retained source/archive for214/215. After the first structural basis
inspection, preflight was extended to all r>=4 before any certificates.
No earlier closure was found. The [full producer](mixed_weight2.py) retains
ten complete models r4..13:94 atoms,142 native coefficient equations.
The eleven-slot r5/r7 maps were actual raw13 improvement targets, not
support-pruned or collision-axis-normalized models.

Even r has complete polynomial invariant modules. The b=0 boundary kills
all V-only terms; three diagonal coefficients with seven nonzero supports
remove the remaining mixed terms of r4. Positive univariate leading
coefficients finish every even ratio. Odd r5/r7 retain one mixed-first
coefficient nu. Two complete nu1 native models each have a full-ring unit
with nine multiplier terms,maxdegree1, obtained by the bounded
[primary certificate producer](certify_mixed_weight2.py) in less than a
second. Their actual exact multipliers, not primary statuses, are replayed.

The [independent checker](verify_mixed_weight2.py) verifies source descent
of the half-cover, all ten actual source determinants and complete odd
supports/parities, five even boundary factorizations, three diagonal
identities/seven supports, eighteen leading cases, both units and all34
generators, three formal-r odd univariate identities and the constant-first
top coefficient, and eleven complete generic inverse charts. There is no
new timeout or pending process. Old r1/r2/r3 UNKNOWN/timeouts are untouched.

The [canonical theorem](MIXED_WEIGHT2_CLASSIFICATION.md) gives exact
one/two-shear normal forms for all r>=4 through7, with an infinite written
module degree bound rather than a finite ratio scan. No external plane
theorem is imported. No smaller positive example or whole-profile closure
is found; all global simplicity minima remain open.

The proposed complete resonance-plane census is now completed by KCUB36
below, with a universal rank-one lemma and a finite fourteen-type frontier.
No new positive bound or global optimum follows.


## Complete native weight-ray frontier (KCUB36)

The previous turn made progress as KCUB35. New preflight found the archived
degree-four arrangement note/script, which already supplies the general
pair-kernel idea. No old radical or unregistered classification was imported.
A new written all-degree/all-dimension rank-one lemma gives identity or one
monomial shear, leaving rank-two native relation kernels for degree7.

The [primary census](weight_resonance7.py) produces348 labelled atoms,
182 relations,154 primitive relation planes and896 weight rays. Exact
acyclic/rank gates and prior full weighted imports leave463 bases. The
[bounded coefficient producer](certify_weight_resonance7.py) constructs451
complete mixed native ideals and exports coefficient-power identities for
all positive-output nonlinear coefficients in447 models. The four remaining
models retain full native five-second timeouts, scripts and logs, unchanged.
No support pruning, coefficient grid or generic parameter-field root is used.

The first independent replay stopped before certificate multiplication:
a generic block determinant variable shadowed the census dictionary. After
renaming those symbols, the [full checker](verify_weight_resonance7.py)
passes in a few seconds. It reconstructs pair kernels by rational row
reduction and native ideals by derivative-atom Leibniz expansion, then
multiplies594 identities with153878 multiplier monomials,maxdegree15 and
targetpower16. It also verifies99 generic rank-one shear inverses and343
leading cases. Two zero-weight cases close by polynomial units.

The [canonical theorem](WEIGHT_RESONANCE_CENSUS7.md) reduces every displayed
diagonal-Gm degree7 raw<=12 collision to exactly fourteen retained weights:
four mixed signatures,full111 and nine opposite-sign zero-weight types.
Known112 remains the exact13-term control. This is a necessary weighted
frontier, not a global lower bound or a new positive witness. Every older
UNKNOWN/timeout, including broader111, is untouched.

The complete12-slot(1,-2,-4),(1,-2,-5) models are now closed by KCUB37
below. Their original native timeouts remain unchanged. The14/15-slot
(4,-1,-2),(3,-1,-2) models,full111,the nine zero-weight cases,source-frame
changes and unweighted constructions remain routes toward the global goal.


## Both complete twelve-slot frontier types (KCUB37)

The previous turn made concrete progress as KCUB36. Current preflight read
its full scope, native timeout inputs and the earlier123 boundary method;
no124/125 closure was found. Neither old native ideal was rerun or given
more time. The [changed model producer](weighted124125.py) retains both
full twelve-slot maps,24 atoms and30 native coefficient equations.

At the nonlinear-first graph boundary, integer ratio4/2 allows subtraction
of a leading G0^2 power, ruling out the apparently balanced quadratic/quartic
branch. Both graph cases then have affine G0. The exact transverse resonance
3deg(H1)=r+1 immediately excludes r4; r5 reduces through a determinant-2
constant solve to two one-variable polynomials with a four-monomial unit.
Univariate-first cases close by roots, field units and a weight-compatible
target-square shear; no old degree bound is reused after the shear.

The [independent checker](verify_weighted124125.py) passes both full source
and graph reconstructions, all universal coefficient/power identities,
the tiny unit, and five complete generic inverse charts. The remaining
linear-first maps have four explicit automorphism families. One r5 family
preserves kappa*b^2-a*c and inverts by lambda sign reversal; no nonlinear
first component or collision survives. No external theorem is imported
by these two complete closures.

The [canonical result](WEIGHTED_124_125_EXCLUSION.md) tightens the displayed
diagonal-weight raw<=12 frontier to twelve types. The four original native
timeout records remain byte-for-byte terminal evidence; two weights are now
mathematically settled and the312/412 types remain unresolved. Their next
query requires new complete structural reductions, not a longer unchanged
solve. All positive bounds, unrestricted62/12 profiles and global minima
remain unchanged/open.


## Last mixed signatures and a sharp restricted minimum (KCUB38)

The previous goal turn made progress as KCUB37. Scope preflight found no
312/412 closure. The full new models have29 atoms and52 equations. The
412 b0 unit boundary and three active-support diagonal identities give
only triangular shears. The312 L0 stratum closes by a boundary unit and
a formal coefficient product derivative in independent U,V.

For312 L!=0,the changed producer normalizes L1 and compiles all four first
coefficient strata with constant pivots. Three units finish immediately;
mu1 leaves seven variables/twenty equations and a retained10-second timeout.
No unchanged native ideal was rerun. Its actual b0 determinant gives a
complete product locus with four branches,including T=1/3. Each new branch
has an exact unit. All original four native five-second timeouts and the
changed mu1 timeout remain terminal solver history,not relabelled successes.

The [independent checker](verify_weighted3412.py) reconstructs both original
sources,all eight reductions and48 pivots,and multiplies seven units for123
generators with1096 multiplier terms,maxdegree7. It directly reconstructs
the existing13-term map with full determinant1,degree7,weights and the two
integer collision points. The [new theorem](DISTINCT_NONZERO_WEIGHT_MINIMUM.md)
therefore gives matching raw13 throughdegree7 and degree7 in the distinct
nonzero diagonal-weight class. No new witness,global minimum or source-frame
sparsity bound is asserted.

The current broader weighted raw<=12 frontier has ten types:full111 and
nine opposite-sign zero-weight signatures. A proposed next route is their
complete degree-seven t0 unit reduction,only after reading the retained
degree-four quotient work. General source frames,unweighted constructions
and higher-degree support remain alternatives. All global simplicity goals
stay open.


## KCUB39: all-degree exactly-one-zero closure (15 November 2026)

The saved next gate was executed:current scopes/catalogue,zero-weight history
and the full strict-only KCUB23 proof were read before computing. Archived
degree-four quotient/radical work was provenance,not a degree-seven closure.
The t0 product gives the complete pure-x boundary. A prospective one-layer
degree-seven argument was superseded by the stronger attained-support-slope
filtration proof:positive slope gives a positive rank-one principal map with
constant Jacobian,contradicting its leading coefficient. The R0 tail is a
univariate unit. Same-sign complete modules have elementary or matrix inverses.

The unbounded proof is written in ZERO_WEIGHT_CLASSIFICATION.md. The producer
and independent Fraction checker cover all nine complete degree-seven models
(240 atoms,438 coefficient generators,162 pure-x removals,78 remaining atoms,
twelve principal faces),27 rank-one and16 tail leading controls,nine generic
inverses and a full4x4 identity stabilization of the established raw13 control.
This gives sharp collision dimension4 for exactly one zero diagonal weight.
No new witness or solver query was made;all native and111 timeouts and unused
primary certificates are preserved. The prior goal turn was progress KCUB38;
this turn is progress KCUB39,not goal completion. The degree7 diagonal-frame
raw<=12 frontier is now only full111;the strict lower-flag linear-Gm subcase
is closed using KCUB23 only in its stated strict scope. No unrestricted whole
profile or global positive bound changes. Next full111 work needs new structure,
or leave this counted diagonal frame/construction.


## KCUB40: full111 closure by curve/Euler-Wronskian lifts (15 November 2026)

Previous goal turn:progress KCUB39. Current preflight read root instructions,
current registry,index,memory,fullKCUB23 proof,all unresolved chart scopes and
UAAFFINER1's prior Riemann-Hurwitz Laurent-support conclusion. No unchanged
solver query was run. The full111 leading products have coprime degrees7
and4/5/6,forcing f3 to be a cube. This gave a complete curve split.

A formal quadratic-algebra place excludes the Z2 branch. For the punctured
line,a new elementary proof of the needed Laurent-support lemma uses the
constant Euler Wronskian,coprimality,D2a=Ha and highest/lowest exponents. Native
pole bounds leave k1/k2 lifts. Three direct symbolic coefficient identities
close those and the graph branch;vertical roots reduce to a polynomial unit.
The independent Fraction checker reconstructs33 native atoms,the full source
determinant,21 post-leading coefficients and every finite algebraic control.
The universal coverage and pole/Laurent/lift proofs are written,not inferred
from bounded tests or externally reviewed. Producer completed in about6sec;
there are zero new solver queries or timeouts. Every old input/script/log/
checkpoint is pinned and unchanged;the large old primary unit remains unused.

This is progress KCUB40:the degree7 diagonal-frame raw<=12 frontier is empty.
Only112 collides and has exactly13 terms;this proves matching raw13/degree7
in3D and all-degree dimension3 for every nontrivial displayed diagonal weight.
It does not close arbitrary source-frame raw sparsification or the global
simplicity objective. Next audit the112 general source/output frame orbit
structurally before any search,or leave this torus/construction class. Existing
positive bounds and62/12 unrestricted whole profiles are unchanged.


## KCUB41: complete affine orbit of raw13 (15 November 2026)

Previous goal turn:progress KCUB40. Current preflight inspected instructions,
registry,memory,index,the raw13 source and existing coefficient/frame results.
The apparent general-source escape has a simpler invariant:the sum of
homogeneous output-span ranks is13 at the origin and is invariant under
all linear source/output changes. Full translated jets extend this to every
source centre,with exact sums13 on the c-axis and14 off it. Three subsecond
two-variable SymPy rank-ideal calculations returned units during discovery;
all are retained unused. Two selected-minor units with31 multiplier terms
now independently prove the remaining normalized-centre ranks. No grid or
old solver query was rerun. Producer/checker call no solver.

Equality forces one top monomial,so UFD makes the source matrix monomial.
Nonzero c-axis translations have18 source columns;at the origin13 pure-axis
columns force monomial output matrices. Thus the complete affine left/right
orbit has sharp13,and only original-centre monomial frames attain it. The
independent checker uses full six-variable translation rather than producer
binomial entries and checks all348 coefficients,both units,all rank branches
and the old determinant/collision. Universal invariance and equality proofs
are written,not formalized or externally reviewed. Combined KCUB40 closes
affine-coordinate sparsification of every classified degree7 torus collision.

This is progress KCUB41,not a smaller positive witness or global optimum.
Next leave this affine3D orbit,use genuinely nonlinear coordinated changes
or change construction. The12-variable cubic affine problem remains separate.
All positive bounds and unrestricted62/12 whole profiles are unchanged.


## KCUB42: complete low-source-degree target function field (15 November 2026)

Previous goal turn:progress KCUB41. Current preflight read instructions,
registry,memory,index,KCUB26 scope and the raw13 source;the all-output target
problem was not closed by the12-variable one-axis derivation theorem. A
naive leading-degree shortcut failed because the cusp relation produces
S(H),T(H) of degree8. These exact expressions are retained in discovery.
SAGBI documentation was inspected,but no completion query was started.

The changed method derives the full generic cubic field and explicit source
recovery. This makes arbitrary rational-target membership finite on source
degree<=7 and<=8. Exact good-point minors certify ranks116/120 and158/165,
with independent companion/adjugate reconstruction and column elimination.
The kernels are precisely the four affine functions and seven functions
including H1^2,S(H),T(H),respectively. There is no target-degree or denominator
cap and no candidate coefficient grid;only fixed normal-form functionals
are specialized. A first8-checker power table stopped at7 and was corrected
to the actual chosen source degree before successful replay;no certificate
was accepted from that failed run.

The complete9-parameter degree8 target model has27 determinant equations.
The primary Singular basis and lift jobs both completed below their five-
second budgets;all logs/checkpoints remain. Seven exact memberships have
17 multiplier terms,maxdegree1,and leave only the generic quadratic shear
Y+(0,lambda,mu)*Y1^2 with a checked inverse. Nonzero members have degree8
and sharp fixed-source nonlinear count19,attained by a directly checked
old-map target shear. This is not a smaller witness.

This turn is progress KCUB42. All target-only degree<=7 routes are closed
even after affine frame changes and temporary degree growth. Next audit
affine source/output frames of the explicit degree8 two-parameter family;
the fixed-source19 count does not settle it. Nonlinear source coordination
and changed/unweighted constructions remain open. No global positive bound
or unrestricted62/12 profile count changes.


## KCUB43: all affine frames of the first nonlinear target layer (15 November 2026)

Previous goal turn:progress KCUB42. Current preflight read registry,memory,
index,KCUB41/42 and their next gate. A centre-dependent determinant-one
output shear replaces P^2 by(P-P0)^2,avoiding spurious translation terms.
Five formal c2 coefficients and complete alpha0 formulas prove exact ranks
for every centre and L/M branch,using KCUB41 base ranks. Rank sums14/15/16
already forbid improving13. Coupled top/cubic rigidity strengthens the
nonlinear lower bound to16:nonmonomial squares cannot have two terms;the
top layers force source axes;three cubic occurrences force the original
monomial frame,which has19 source columns. All <=15 branches are exhausted.

The producer/checker reconstruct483 complete formal shear coefficients,
all base entries and square rows,and both equality models,without a solver
or coefficient grid. All zero strata remain. The matching old target-shear
control has degree8/19 terms,full determinant1 and the saved collision.
The nonlinear-family minimum is only within[16,19],not proved sharp.
The complete degree<=8 rational-target plus affine-source/output class has
sharp13 onlyat the old zero-shear equality frames. No smaller witness or
unrestricted profile/global positive bound change follows.

This is progress KCUB43. Do not run another affine search in that family.
Next proposed gate:a genuine3D positive-axis quadratic source shear with
full rational-target repair;KCUB19/26 scopes were read and concern the
separate12-variable cubic,so they do not close it. Derive the new field
operator and every parameter branch before any universal exclusion.


## KCUB44: uniform nonlinear positive-source obstruction (15 November 2026)

Previous goal turn:progress KCUB43. Preflight inspected current registry,
memory,index,full KCUB42 field machinery and KCUB19/43 scopes. Seven cheap
fixed-parameter rank pilots suggested low-source7 kernel2 for nonzero
quadratics,but their saved discovery statuses are not universal exclusions.
No symbolic parameter-minor search was inferred from those values.

The changed proof uses the source torus:parameter weights3/4/5 contract
every nonzero tuple to zero. A new full base-field9 certificate has rank
213/220 at107 good target points;independent companion/adjugate reconstruction
and column elimination give determinant116435059 mod1000000007. The exact
kernel is the old7 functions,all of degree<=8. Formal inverse-source
coefficients stay in the original field. Taylor orders d,2d,3d have five
nonzero heads with distinct a exponents5/2/3/1/0,forcing every weighted
initial field member into span(1,H1). The finite filtered-space dimension
then proves transformed low-source8 kernel2 for every nonzero coefficient
triple,including all mixed and degenerate quadratics. Zero q is separate.

No rational target repair can be Keller throughdegree8,with no target
degree/denominator bound. Generic source maps have row degrees7/9/10;the
q=b2 control has37 nonlinear terms,full determinant1 and transported
collision. Its degree10 is only an upper against lower9;degree9 in the
transformed field remains unknown. The original-field gap9 also extends
the old original-map target/affine closure through endpoint degree9.
Universal descent/induction/initial-space proofs are written,not formal
or externally reviewed;all finite premises are independently checked.

This is progress KCUB44,not a smaller global witness or optimum. Next
preflight the negative-axis quadratic source shear with coefficient weights
-3/0/3,including its1+3s=0 resonance and old degree9/19 control. The positive
contraction argument does not transfer. Whole-profile62/12 and every
headline positive bound remain unchanged.
