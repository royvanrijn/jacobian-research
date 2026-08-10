# Closing attacks for the plane one-pair and degree frontiers

This note ranks the next attacks by whether success changes the rigorous
degree frontier. It separates theorem-bearing reductions from useful but
nondecisive reinterpretations.

## One-pair closure ledger — 2026-07-23

The boundary-only program needs one correction.  A weighted SNC tree that
passes the intrinsic `A^2` gate describes a possible completion of the
*source*.  It does not yet describe a possible Keller map.  The `(72,108)`
package demonstrates the distinction: its complete source boundary is a
valid `A^2` boundary and passes every present numerical gate, even though
the coefficient equations are inconsistent.

The object to classify must therefore be the map-decorated package

\[
 {\cal B}=(\Gamma_X,Q_X,k,p;\ \Gamma_Y,\varphi;\ e,f;\ V_D,\Sigma_D).
\]

Here `Gamma_X,Q_X,k` are the source boundary and canonical data, `p=f^*L`
is the pole vector, `Gamma_Y` is a resolved target boundary, `varphi` is the
map of dual graphs, `e,f` are normal and residue degrees, and `V_D` is the
three-section linear series induced on every noncontracted boundary
component, with vanishing data `Sigma_D` at nodes and critical points.
The current compiler supplies the first four entries and part of `e,f`; the
next attacks construct and test the missing entries.

Every proposed finite-normalization reconstruction now has two preliminary
global gates before a finite-support residue or a large coefficient search is
meaningful:

1. saturate the different away from the declared omitted boundary and reject
   any residual height-one ramification;
2. when the normalized open contains a certified `G_m^2` core, compute the
   valuations of both character generators on **all** codimension-one fill
   primes and require the resulting matrix to be unimodular.

The first gate detects the companion tame divisor in the canonical odd
wild-boundary rows.  The second detects `Cl=Z/(N-1)` in the balanced rows
after their different has vanished.  The enlarged coefficient audit also
detects a free unit and `Z/(N+1)` at gluing exponent zero, a wild index-`N`
branch at every exponent `a>=N`, and an extra different divisor for every
factor of a general `C(P)` away from `P=0`.  For arbitrary `C(P,Q)` the same
support gate leaves only `cP^aQ^b`; the base-change Euler/transfer theorem and
the uniform full-source calculation
`Cl(U_(2,c))=Cl(D(xu))=(Z/2)^2` now close every `b>0` row.
The same proof gives an immediate atlas sieve for future architectures: a
principal reduced packet of fibres with multiplicities `m_i` and generic
unit rank one forces `(direct_sum Z/m_i)/<diagonal>`, so any prime shared by
two multiplicities rejects affine-plane reconstruction before a full class
group computation.
For multi-retained balanced covers there is now an earlier support gate:
`A*H_T-A'*H=P^(N-1)*Q*(A-T*A')` forces
`A=a0+T*B(T^p)`.  Consequently only cover degrees `1 mod p` enter the source
search.  The normalization chart supplies the stronger global gate
`#(C_A-E_A)(F_q)=q^2+(n_q(A)-1)q` and
`chi_c(C_A-E_A)=deg(A)` after the roots split.  It excludes every nonlinear
retained polynomial, including the former degree `7,10,13,11,15`
support-only rows in characteristics `3,5,7`.
The complete degree-seven scan makes the mechanism exact: four coefficients
fail over `F_3`, while `A=T^4+T+1` and `A=T^4+2T+1` match `A^2` over
`F_3` and `F_9` but have `810!=27^2` open points over `F_27`.  No unit or
class-group calculation is needed to reject them.  The next search must
change the retained-boundary or gluing architecture.
This is the general
[normal-core Smith criterion](BOUNDARY_LATTICE_PREFILTER.md#dual-torus-core-localization),
not the local two-by-two Smith profile of the logarithmic differential in
Attack E.  The same gate computes `ord([L1])=N-1`, rather than only the
abstract cokernel.  Passing both remains necessary, not sufficient unless an
extended affine toric action is separately certified: affineness, coordinate
reconstruction, and the coefficient equations are still later tests.

There is also a new exact reduction.  In terminal Case 2 of the archived
`(72,108)` pair, the a priori residue-cover degrees were `1,2,4`.  For
degrees `2` and `4`, reconstructing a completely general monic polynomial
right component from the leading coefficients and dividing both residue
coordinates by its powers gives exact unit ideals with 9 and 12 generators.
This uses `(J3),(J2)` and the part of `(J1)` that determines `G`, but no
`(J1)` compatibility equation and no `(J0)`.  It first reduces Case 2 to
cover degree one and image degree twelve.

The remaining row is now also empty.  The forced degree-twelve vertex is
`G_12 != 0`.  Localizing at it with `w*G_12-1`, the seven residual `(J1)`
compatibility equations--each an eight-term cubic--generate the unit ideal
over the exact degree-35 first-block field.  No `(J0)` equation is used:

\[
 \boxed{\text{terminal Case 2 is excluded at the }(J1)
        \text{ endpoint-compatibility stage}.}
\]

The certificates are
[`cas/audit_case2_residue_strata.py`](cas/audit_case2_residue_strata.py)
and
[`cas/case2_infinity_resolution.py`](cas/case2_infinity_resolution.py).

## Poisson-square primary filtration — associated-prime gap completed

The three-band coefficient scheme is not Cohen--Macaulay.  The earlier
normalized free-resolution attack was aimed at proving the opposite and is
now retired.  On the exact `d0=1` chart, separator saturations of the three
minimal components fail to reconstruct the coefficient ideal.  A coefficient
filtration replaces the intractable full decomposition:

\[
 I_0\subset I_1=I_0+(d_3)\subset
 I_2=I_0+(d_3,d_2).
\]

The exact primary counts are:

| ideal | primary radicals | dimensions |
| --- | --- | --- |
| `I0:d3` | `S,S_C,S_A,K_C,K_A` | `2,2,2,1,1` |
| `I1:d2` | `T,S,S_C,S_A` | `3,2,2,2` |
| `I2` | `T,C0,A0` | `3,3,3` |

Here `K_C=S∩S_C` and `K_A=S∩S_A` are the two irreducible
core/intersection curves.  The two multiplication exact sequences prove
that the normalized algebra has exactly eight associated primes.  Restoring
the `G_m` factor gives three dimension-four minimal primes, three
dimension-three embedded primes, and two dimension-two embedded primes.
The executable certificate is
[`cas/poisson_square_normalized_defect.sing`](cas/poisson_square_normalized_defect.sing).

This closes the associated-support gap and changes the next attacks:

1. **Normal-module fibers — first milestone completed.**  The two cyclic
   presentations have Fitting ideals `I0:d3` and `I1:d2`.  Exact transverse
   Hilbert vectors are now known on every associated stratum.  On the `d3`
   layer they are `(1,2)`, `(1,2,2)`, `(1,2,3,1)`,
   `(1,3,5,7,6,3)`, and `(1,4,8,10,6,1)` on
   `S,S_C,S_A,K_C,K_A`; on the `d2` layer they are `(1)`,
   `(1,3,3,1)`, `(1,2,1)`, `(1,1)` on `T,S,S_C,S_A`.
   These vectors have now been replayed over the rational function field of
   every branch, so they are generic rather than sample-only.  The next
   presentation step is also complete on the `d2` surfaces: `S` needs three
   quadrics and one cubic and has socle dimension two, `S_C` has an exact
   four-relation socle-one presentation, and `S_A` is the dual numbers.  The
   `d3` presentations are now finite as well: generator/relation counts are
   `2/3`, `2/3`, `2/4`, `3/9`, `4/18` on
   `S,S_C,S_A,K_C,K_A`.  Their standard monomial bases determine exact
   multiplication tables.  The remaining milestone is to simplify the last
   three presentations to coordinate-invariant normal forms.
2. **Filtered lower-band action — implemented at reduced-support level.**
   When a Newton architecture produces the
   same top three bands, reduce every proposed lower bracket layer on the
   eight associated-prime modules before forming its coefficient ideal.
   A candidate dies as soon as one lower layer is a unit on a required
   dense chart.  The new search-facing filter reports `preserved`, `cut`, or
   `eliminated` after localized exact substitution.  The next step is to
   compute the action on the nilpotent primary fibers, not only their reduced
   supports.
3. **One-band enlargement.**  Add one Laurent monomial orbit at a time and
   recompute the `d_top,d_next` colon ladder.  The milestone is a finite
   adjacency graph recording which support additions preserve the tangent
   component and which merge or remove the five embedded strata.
4. **Local conductor comparison.**  Complete the five colon-primary pieces
   along `S,S_C,S_A,K_C,K_A` and compare their conductors with the
   log-boundary differents.  A mismatch eliminates a Newton chain before
   its lower coefficient equations are expanded.

These attacks are smaller than a full primary decomposition and interact
directly with the log-boundary compiler: the colon filtration supplies the
scheme modules on which boundary valuations and lower bands must act.

## New reduction: the repeated-tail row is a triple-root problem

The apparent source conflict around the 2017 row

\[
 (8,40)\longrightarrow(8,28)\longrightarrow(11/4,7),
 \qquad (m,n)=(3,2),
\]

is not a contradiction. The 2017 complete-chain program emits a necessary
combinatorial over-approximation. In its source, the filter coming from the
lower-side prohibition of corners
\(\wp(n',n'-1)\) is commented out. The 2016 paper says in a remark that the
row leads to the forbidden corner `(8,4)`, but does not prove that transition.
The 2022 paper proves the analogous transition for the companion row starting
at `(8,32)`, where the residual vertical factor has degree one.

The part of that proof depending only on the common tail does extend. The
edge from `(8,28)` to `(11/4,7)` has weight `(4,-1)` and forces

\[
 \operatorname{en}(F_1)=(6,21)=\frac34(8,28).
\]

Thus \(q_1=4\). The divisibility theorem gives \(4=q_1\mid d_0\), while
\(d_0\leq\gcd(8,28)=4\), so \(d_0=4\). Consequently

\[
 \ell_{1,0}(P)=R^{4m},\qquad
 R=\kappa x^2y^7p(y),\qquad \deg p=3,
\]

with nonzero constant and leading coefficients. Translating a nonzero root
of multiplicity \(r\) to the origin produces the normalized corner
`(8,4r)`. Any simple root therefore gives the forbidden corner `(8,4)`.
The partitions `(2,1)` and `(1,1,1)` are excluded. The source conflict is
reduced to the single branch

\[
 \boxed{R=\kappa x^2y^7(y-\lambda)^3.}
\]

The remaining triple root translates the vertical edge to

\[
 (8,40)\longrightarrow(8,12).
\]

The complete-chain length bound for this edge is three.  An exact enumeration
using a *superset* of the possible-last-lower-corner list has open-chain
counts

```text
1 -> 6 -> 3 -> 0
```

and no final corner.  The same implementation recovers the published final
corner `(11/4,7)` from the companion `(8,32)->(8,28)` input.  Thus the
triple-root branch has no complete-chain escape and the repeated-tail row is
excluded.  The other five raw `(96,144)` rows are unaffected.

The factor audit is
[`cas/frontier_96_144_source_audit.py`](cas/frontier_96_144_source_audit.py);
the no-escape certificate is
[`cas/complete_chain_no_escape.py`](cas/complete_chain_no_escape.py).

## Attack A — triple-root continuation — completed

This attack meets its kill criterion at the complete-chain stage, before any
approximate-root or bracket-band expansion.  It removes one of the six raw
`(96,144)` chains.

## Attack B — F2 polyhedral coefficient route — B0 completed, route parked

This was the original route to an F2 exclusion.  It remains mathematically
valid, but it is no longer the programme-wide critical path: the Kummer-orbit
and terminal-residue theorems bypassed the missing lower masks for the
selected principal chain and produced the exact target row `(e,f)=(1,6)`.
The coefficient route is retained here with its present status.

**Milestone B0 (support envelope) — completed.**  Pull each original polynomial monomial
through `x=X^5`, `y -> y+lambda/X`.  A descendant has

\[
 (i,j,k)\longmapsto(5i-k,j-k),\qquad 0\leq k\leq j,
\]

so `a-b=5i-j` is invariant along the translation string.  Combine this with
the exact standard-pair Newton polygon and the two forced translated edges
to emit a finite set of candidate lattice points through the entire bracket
range.  Every point must carry either an original polygon inequality or a
forced-edge provenance tag.  The older `(50,75)` values `gamma=2,3` cannot
be used: their preliminary reduction is explicitly unproved and concerns
the `(2,3)` member.

The exact degree/terminal-halfspace envelope is now emitted by
[`cas/classify_f2_75_125_layers.py`](cas/classify_f2_75_125_layers.py).  Its
upper window contains all 35 zero layers, 665 band-pair incidences, and 978
jet-reduced linear parameters.  The complete record extends through layer
`-200` and has `240` zero layers, `13,741` band pairs, and `2,418`
jet-reduced parameters.  This meets B0 as a certified over-envelope; it is
not an exhaustive list of B1 polygon masks.

If the direct coefficient route is resumed, introduce binary variables for
its lattice points below the common-power band and impose, in order:

1. convexity and the two forced edges;
2. endpoint nonvanishing;
3. Minkowski compatibility of the `(3,5)` leading powers;
4. vanishing of bracket layers `39` down to `5`; and
5. the unique monomial \(t^4z^4\) on layer `4`; and
6. vanishing of every lower layer `3` down to `-200`.

At each layer, use support incidence before coefficient equations: a uniquely
represented bracket monomial is forbidden, while a required right-hand-side
monomial must have at least one representation.

**Milestone B1.** A finite list of support masks, with a machine-checkable
proof that every omitted lattice point is incompatible with one of the six
conditions.

**Milestone B2.** Split only those masks by coefficient cancellation and
derive the actual gamma branches.

**Kill criterion.** No support mask survives, or every survivor reaches a
weighted-Wronskian block with a nonzero de Rham class.

This reverses the expensive order “guess gamma, then expand 35 layers”:
support incidence should eliminate most branches before any field extension.

## Attack C — selected F2 branch scale and terminal row — completed locally

The Kummer-character transfer proves that a nonzero fifth-root fiber is one
orbit, excludes the zero-root strata, and makes simple cofactor roots
spectators.  The selected squared factor has source ray `(12,-17)` and target
ray `(5,2)`.  Its residue map is the parameter-free degree-six cover

\[
h(s)=\frac{125s(s+1)^5}{(9s^2+15s+5)^3},
\]

with `(e,f)=(1,6)`, passport `(5,1)|(3,3)|(3,1,1,1)`, residue-different
packet `(4,2,2,2)`, and geometric monodromy `A_6`.  Its natural action is
four-transitive and has trivial target-fixed deck group.  This completes the
positive local-scale objective for the selected chain; it does not make the
full B1 coefficient masks exhaustive.

The rescaled cover is Belyi and its regular `A_6` closure has triangle
signature `(5,3,3)` and genus `25`.  More immediately for the boundary graph,
the target toric nodes have three preimages in the source-divisor interior.
They force three attachment points carrying different contributions
`(4,2,2)`; the remaining contribution `2` is at the source endpoint over the
smooth third branch value.

The resulting global consequences are already fixed: geometric degree is at
least six, or at least twelve for two distinct packets over the same target
divisor; the infinity-centered row receives no affine-sheet increment;
purity forces a separate affine ramification row; and global geometric
monodromy has `A_6` as a nonabelian simple composition factor.  The live
attack is now to glue the one-packet and two double-packet attachment cases,
classify spectators, and run the class-group, unit, canonical, different, and
meridian ledgers.
<!-- status-consumer: PF2GC1 33dbc5ff48b5d064 -->

## Attack D0 — non-birational Case-2 residue strata — completed

The old Case-2 cover rows `delta=2,4` are empty by the exact decomposition
certificate above, and the last `delta=1` row is empty by the `(J1)` endpoint
certificate.  This retires the former six-row valuation table: only the
three Case-1 rows remain.

## Attack D1 — run the same decomposition sieve on Case 1

**Input.**  The exact alternate-chart residue
`[1:P2(0,r/u):Q2(0,r/u)]` on each of the two sign branches.

**Band prerequisite -- completed.**  The archived certificate records `P`
only through `z^-5` and `Q` only through `z^-4`, but the exact continuation
in [CASE1_FULL_BAND_CONTINUATION.md](CASE1_FULL_BAND_CONTINUATION.md) solves
bracket layers `-4,...,-11` and reconstructs the full tails through `z^-8`
and `z^-12`.  All eight systems have full column rank and add 66 necessary
compatibility equations without new moduli.  The next input-construction step
is now the certified alternate-chart substitution and restriction to `X=0`;
zero-filling is neither needed nor permitted.

**Calculation.**  Reconstruct the general monic right component of degree
two and four from the top residue coefficients, reduce both coordinates
modulo its powers, and test the remainder ideals before importing the final
branch unit identities.

**Kill criterion.**  A unit ideal removes that cover row.  A nonempty ideal
must emit its dimension, generic residue field, and a sample point; merely
timing out is not a result.

**Why it matters.**  This is now the only live case/cover sieve, and the
Case-2 implementation is reusable.  The search has narrowed from an
unspecified Newton-tail derivation to a deterministic transport followed by
two already-implemented polynomial remainder tests.

## Attack D2 — Case-2 target resolution — closed before implicitization

The degree-twelve endpoint open is empty after the seven `(J1)`
compatibility equations, so the proposed implicitization and degree-29
harmonic-cover calculation have no surviving input.

For comparison, before imposing compatibility, the infinity chart has
orders `(4,12)`.  Cancelling the common cubic tangent gives the
translation-invariant characteristic numerator

\[
 K_{13}=2C_8G_{11}-3C_7G_{12}.
\]

On `G_12*K_13 != 0` the generic infinity branch is the `(4,13)` cusp.  Its
seven exceptional rays are
`(1,1),(1,2),(1,3),(4,13),(3,10),(2,7),(1,4)`, with self-intersections
`-2,-2,-5,-1,-2,-2,-2`.  This supplies the target graph that D2 would have
used, but the exact endpoint unit ideal proves that no compatible
degree-twelve stratum occurs.

## Attack D3 — boundary linear series and the Pluecker budget — retired for Case 2

For a noncontracted rational boundary component `D`, the three target
sections restrict to a basepoint-free linear series

\[
 V_D\subset H^0(D,\mathcal O_D((Qp)_D)).
\]

Its vanishing sequences at boundary nodes, target contacts, and critical
points must satisfy the Pluecker formula.  This data is invisible to
`Q,k,p` and to the normalization different.

On the smallest Case-2 gcd stratum the residue net is a `g^2_12`.  The
forced origin orders give vanishing sequence `(0,2,4)`, of weight `3`.
Homogenizing coordinate degrees `(0,8,12)` gives sequence `(0,4,12)` at
infinity, of weight `13`.  Since a `g^2_12` has total ramification weight

\[
 3(12-2)=30,
\]

every solution must place exactly `14` further units.  Equivalently,

```text
C'*G''-C''*G' = t^3*W_14(t),  deg(W_14)=14.
```

**Calculation.**  Factor or subresultant-stratify `W_14` together with the
target singularity resolution from D2.  Allocate every root either to a
singular branch, a flex, or a boundary node, and impose the corresponding
vanishing sequence on the harmonic graph cover.

For higher gcd degree use the exact identities

```text
C'*G''-C''*G' = H^2*(c*g'-c'*g),
K^2*(C'*G''-C''*G') = -2*H^3*(A*g-c*E).
```

The seven a priori gcd degrees `1,...,7` leave relative-Wronskian degrees
`15,13,11,9,7,5,3`.  Attack D4 now excludes the last two rows exactly, so
the remaining pre-compatibility list is `deg(H)=1,...,5`, with degrees
`15,13,11,9,7`.

**Kill criterion.**  On the degree-one gcd row, the forced node/contact
weights exceed `14`, or every partition of `14` violates local Hurwitz.
On a higher row, use its displayed relative-Wronskian degree instead.

This remains a useful geometric description of the pre-compatibility
family, but it is no longer a closure attack: D2's endpoint certificate
already makes the entire Case-2 family empty.

## Attack D4 — Case-2 gcd strata as compact coefficient certificates

For the already excluded `(72,108)` case, impose

\[
 H=\gcd(C',G'),\quad C'=Hc,\quad G'=Hg,\quad
 B=Kc,\quad F=Kg
\]

directly on the exact `(J3),(J2)` solution. Compute certificates separately
for the exact gcd degrees.  The maximal row is now complete.  When
`deg(H)=7`, `C'` divides `G'`; the three degree-`0,1,2` coefficients of
`remainder(G',C')` together with the coefficient of `t^19` in `(J0)`
generate the unit ideal over the exact degree-35 field.  The input hash is
pinned in
[`cas/audit_case2_maximal_gcd.py`](cas/audit_case2_maximal_gcd.py), and no
residual `(J1)` compatibility equation is used.

The degree-six row is now complete too.  Write `C'=H*(t+v)` with
`deg(H)=6,H(0)=0`.  The five equations `C'(0),H(0)`, the last two
coefficients of `remainder(G',H)`, and `(J0)_{19}` generate the unit ideal.
The exact input is pinned in
[`cas/audit_case2_gcd6.py`](cas/audit_case2_gcd6.py).

**Next calculation.**  Compute certificates for the five
pre-compatibility degrees `deg(H)=1,...,5`; in the first stratum impose the
forced origin orders `(1,2,3,3)`.

**Success criterion.** Replace the archived four-residual unit ideal by
smaller stratum-aware identities whose factors have geometric meaning.
On `deg(H)=1`, use the exact factor `t^3 W_14` and stratify the remaining
fourteen roots instead of treating all coefficients symmetrically.

**Priority.** Lower. This improves the intrinsic explanation but cannot
raise the degree bound, because `(72,108)` is already exactly excluded.

## Attack E — logarithmic second-Chern defect

The
[`log-conductor degree-shift theorem`](LOG_CONDUCTOR_DEGREE_SHIFT.md)
now proves that this is the first intrinsic invariant not already erased by
the complete determinant ledger.  It exhibits two integrable local
Jacobian matrices with determinant `uv` and identical generic branch Smith
profiles whose cokernels are respectively `R/(uv)` and
`R/(u) direct-sum R/(v)`.  Their nodal `Fitt_1` ideals are `R` and `(u,v)`.
It also proves that a conductor mismatch belongs naturally to `H_Z^1`, so
the `H_Z^0` purity theorem cannot replace this calculation.

<!-- status-consumer: LCDS1 5b4d92acd50d6c41 -->

The first terminal-node audit is now complete.  The three interior F2
attachments have local form `(pi,xi)=(tau*w^-2*unit,w^e*unit)` with
`e=5,3,3`.  Two source blowups are required; both nodes in each resulting
chain have invertible tame exponent matrices.  Together with `s=0`, all four
terminal preimages of target nodes have zero logarithmic cokernel.  At the
smooth-target endpoint `s=infinity`, boundary support forces the cyclic model
`R/(w^3)`, whose reduced support is smooth and whose normalization defect is
also zero.  The carrier continuation now closes the next packet as well: the
squarefree spectators, the double-row fivefold attachment, and all aligned
principal-arm nodes have tame exponent determinants `1`, `5`, and `3`.
Their common fan refinements strengthen the source lower bounds to `27/48`
components but again have zero normalization defect.  Attack E must therefore
move to the upstream carrier-extraction chain, the outgoing terminal tail,
the affine purity row, uncompiled global centers, or a positive global
localized-second-Chern remainder.  Assigning a positive localized defect to
any marked terminal, carrier-local, arm, or spectator slot is no longer
viable.

<!-- status-consumer: PF2LNP1 e4f0f231bf7494d5 -->

<!-- status-consumer: PF2CLP1 41625dd5d3f8f898 -->

The determinant of the logarithmic differential gives the current
`R_log=k+1+2p` gate, but a rank-two bundle map has a second invariant.  For
an SNC completion of `A^2`,

\[
 c_2(\Omega_X^1(\log D))=e(\mathbb A^2)=1,
\]

while the integrated class
`c_2(f^* Omega^1_{P^2}(log L))` equals `deg(f)`.  Together with the already
known first Chern classes, these determine the codimension-two Chern
character of the boundary torsion cokernel.

**Calculation.**  In every node chart record the full two-by-two logarithmic
matrix up to invertible row and column operations, `Fitt_0`, `Fitt_1`, the
two generic branch elementary divisors, and the normalization-defect length.
Compute its localized second Chern class in two ways: globally from `Q,k,p`
and locally from these profiles.  A determinant-only or generic-profile-only
record is now known to be insufficient.

**Kill criterion.**  The required residual zero-cycle has negative length
or a length smaller than the forced node contributions.

**Why it matters.**  This is the next intrinsic invariant after the
determinant.  It uses map data but avoids the full coefficient ideal.

## Attack F — bounded valuation-budget falsification

Enumerate the three surviving Case-1 cover rows:

```text
Case 1: 26, 23, 17
```

This is deliberately weaker than D2: attach additional divisorial valuations
through admissible blowups of the audited source tree and rerun the intrinsic
`A^2`, pole, ramification, and Hurwitz gates after every attachment.

**Kill criterion.** No completion realizes a remainder.

**Stop criterion.** One explicit completion for every row proves that this
numerical route is realizability-neutral; then it should be retired rather
than refined.

## Execution order

For the one-pair closure program:

1. Derive the omitted Case-1 bands or prove the D1 truncation lemma.
2. Run D1 on the three Case-1 cover degrees `1,2,4`.
3. Apply E on every Case-1 graph survivor.
4. Continue D4 only as an optional geometric compression of the already
   excluded Case-2 certificate, and only when it yields a smaller
   certificate.
5. Run F once on the Case-1 rows as a falsification control; retire it if
   witnesses survive.

For movement of the numerical degree frontier:

1. Globally attach the certified F2 `(1,6)` terminal row and the simple
   spectator orbits to the original completions.
2. Separate the double-root same-target row (`d>=12`) from the distinct-target
   row (`d>=6`), place the purity-forced affine ramification and the residue
   different across the three node attachments and the endpoint-over-smooth
   incidence, and run the global ledgers.
3. Resume B1 support masks and lower-layer coefficient/de Rham propagation
   only as an optional independent exclusion route; B0 is already complete.

Every stage has a finite artifact: a unit certificate, a target resolution,
a list of harmonic graph covers, a Pluecker partition, or a support mask.
This prevents “more boundary geometry” from becoming an unbounded search.

<!-- status-consumer: C1FBC1 0f14ef01fff25097 -->
