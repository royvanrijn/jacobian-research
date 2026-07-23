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

## Attack B — polyhedral support completion for F2 `(75,125)`

The current F2 record has three proof obligations, but they should not be
attacked as one theorem.  A paper audit adds a prerequisite: the complete-
chain theorem fixes the two displayed edges but does not bound the entire
lower support after the Puiseux translation.  Therefore a binary-support
enumeration is finite only after the following envelope lemma is proved.

**Milestone B0 (support envelope).**  Pull each original polynomial monomial
through `x=X^5`, `y -> y+lambda/X`.  A descendant has

\[
 (i,j,k)\longmapsto(5i-k,j-k),\qquad 0\leq k\leq j,
\]

so `a-b=5i-j` is invariant along the translation string.  Combine this with
the exact standard-pair Newton polygon and the two forced translated edges
to emit a finite set of candidate lattice points through bracket layer
four.  Every point must carry either an original polygon inequality or a
forced-edge provenance tag.  The older `(50,75)` values `gamma=2,3` cannot
be used: their preliminary reduction is explicitly unproved and concerns
the `(2,3)` member.

**B0 success criterion.**  A machine-checkable finite support envelope whose
complement is excluded by cited polygon inequalities.  Until this exists,
B1 must not silently choose a bounding box.

After B0, introduce binary variables for its lattice points below the
common-power band and impose, in order:

1. convexity and the two forced edges;
2. endpoint nonvanishing;
3. Minkowski compatibility of the `(3,5)` leading powers;
4. vanishing of bracket layers `39` down to `5`; and
5. the unique monomial \(t^4z^4\) on layer `4`.

At each layer, use support incidence before coefficient equations: a uniquely
represented bracket monomial is forbidden, while a required right-hand-side
monomial must have at least one representation.

**Milestone B1.** A finite list of support masks, with a machine-checkable
proof that every omitted lattice point is incompatible with one of the five
conditions.

**Milestone B2.** Split only those masks by coefficient cancellation and
derive the actual gamma branches.

**Kill criterion.** No support mask survives, or every survivor reaches a
weighted-Wronskian block with a nonzero de Rham class.

This reverses the expensive order “guess gamma, then expand 35 layers”:
support incidence should eliminate most branches before any field extension.

## Attack C — nested branch scale at the F2 resonance faces

The branch-scale fan work shows that equal radial values do not determine the
stable target at a multiple leading collision.  Attack A removes the
triple-root repeated-tail face.  The remaining live target is a repeated root
of the F2 common-power polynomial \(H\).

For that target, compute the first translated residue polynomial and its
cross-ratio or jet after fixing the radial scale. Feed that second scale into
the existing recursive resonance atlas and then into the log-boundary
compiler.

**Kill criterion.** The nested scale forces a source basepoint or target
component whose pole, ramification, or intrinsic \(A^2\) boundary data is
incompatible.

**Stop criterion.** If the nested parameter is genuinely free and all
intrinsic gates pass, do not add more blowups: return to the coefficient or
support attack with the free residue recorded as a parameter.

## Attack D0 — non-birational Case-2 residue strata — completed

The old Case-2 cover rows `delta=2,4` are empty by the exact decomposition
certificate above, and the last `delta=1` row is empty by the `(J1)` endpoint
certificate.  This retires the former six-row valuation table: only the
three Case-1 rows remain.

## Attack D1 — run the same decomposition sieve on Case 1

**Input.**  The exact alternate-chart residue
`[1:P2(0,r/u):Q2(0,r/u)]` on each of the two sign branches.

**Prerequisite gap.**  This input is not present in the archived
certificate.  Case 1 records `P` only through `z^-5` and `Q` only through
`z^-4`; its full Newton polygons continue through `z^-8` and `z^-12`.
Those omitted bands contribute to the alternate residue.  The next valid
step is therefore either to derive the missing bands from the Newton descent,
or to prove a truncation lemma showing that the degree-two and degree-four
composition remainders use only archived coefficients.  Assigning zero to
the omitted bands is forbidden.

**Calculation.**  Reconstruct the general monic right component of degree
two and four from the top residue coefficients, reduce both coordinates
modulo its powers, and test the remainder ideals before importing the final
branch unit identities.

**Kill criterion.**  A unit ideal removes that cover row.  A nonempty ideal
must emit its dimension, generic residue field, and a sample point; merely
timing out is not a result.

**Why it matters.**  Once the prerequisite is met, this is now the only live
case/cover sieve, and the Case-2 implementation is reusable.  Until then it
is blocked by missing mathematical input, not by CAS cost.

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

**Calculation.**  Record the local two-by-two Smith profile of the log
differential at each boundary generic point and node.  Compute its localized
second Chern class in two ways: globally from `Q,k,p` and locally from the
profiles.

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

1. B0: prove and emit the finite translated support envelope.
2. B1: support-only F2 enumeration inside that certified envelope.
3. C only on resonance faces actually surviving B1.

Every stage has a finite artifact: a unit certificate, a target resolution,
a list of harmonic graph covers, a Pluecker partition, or a support mask.
This prevents “more boundary geometry” from becoming an unbounded search.
