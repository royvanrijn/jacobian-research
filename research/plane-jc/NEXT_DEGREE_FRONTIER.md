# Candidate frontier from 125 through 150

> **Bounded computational regression.**  The full complete-chain algorithm is
> published as pseudocode in Guccione--Guccione--Horruitiner--Valqui,
> arXiv:1708.07936.  No author-supplied implementation or source repository was
> located.  The local script therefore does not claim an independent
> reimplementation: it encodes every relevant row of the paper's accepted
> Section 7 tables, recomputes \((m(a+b),n(a+b))\), checks the printed maximum
> and coprimality, quotients by coordinate reversal, and asserts the exact
> distinct output.  This is a deterministic table regression.

Run:

```bash
python3 plane-jc/cas/frontier_125_150.py
```

It returns 13 normalized unordered pairs:

\[
\begin{gathered}
(75,125),(84,126),(96,128),(88,132),(90,135),\\
(56,140),(84,140),(96,144),(108,144),\\
(42,147),(63,147),(98,147),(100,150).
\end{gathered}
\]

In particular, 125 is an admissible **maximum** in the enumeration, through
\((75,125)\); it is not a claim that a Keller counterexample exists.

## Compiler handoff and loose ends — 2026-07-22

The reusable middle end is now implemented and tested:

```text
published exhaustive Laurent normal form
    -> Laurent lattice supports and bracket layers
    -> weighted-Wronskian block
    -> superelliptic character eigenspace
    -> canonical de Rham obstruction coordinates
    -> optional exact local-system reuse certificate
```

For `(72,108)`, one typed normal-form certificate expands the published chain
to both Laurent polygon cases.  Both compile to the same genus-three first
block, and the six compatibility equations are certified coordinates in the
six-dimensional compact de Rham eigenspace.  The implementation and current
boundary are documented in
[`LAURENT_BAND_FRONTEND.md`](LAURENT_BAND_FRONTEND.md) and
[`NEWTON_DERHAM_COMPILER.md`](NEWTON_DERHAM_COMPILER.md).

Open work, in dependency order:

1. **Exploit the completed `(72,108)` log replay.**  The primary proof yields the
   smooth-branch rays `(2,1),(3,1),(4,1)`, while the Laurent transformations
   have longer adapted base ideals `(t,x^4),(t,x^6),(t,x^8)`.  Their residual
   source chains compile in `4,6,8` blowups, and additive composition proves
   that all three cases share the order-four eight-blowup graph.  Its
   `F_4` transition and affine-plane fill give a unimodular `10 x 10` source
   boundary passing adjunction.  Alternative residue labels are encoded
   symbolically, and the unselected order-three factor is proved to separate
   from the common order-four center and the filled divisor.  Target-infinity
   pole orders are exact on all eight exceptionals.  The full common-graph
   pole vector has no dicritical component, so at least one additional
   global cluster is forced.  The first weighted-Wronskian equation now
   selects a base-order-ten crossing blowup at `E3 intersect E4` followed by
   ten simple children, and excludes the formerly numerical smooth-`E3`
   candidate.  In terminal Case 2, two further target-infinity blowups have
   poles `12,0` and produce a degree-twelve dicritical; the combined intrinsic
   package passes with remaining self-intersection `29`.  On the final
   plane-return edge, the top bracket forces `A=a*r^2,C=c*r^3` for a quartic
   `r`.  The edge alone permits five homogeneous root partitions, and every
   corresponding intrinsic model passes.  The primary split-factor formula
   is sharper: it forces `r=(s-beta)^4`.  Selecting the other order-three
   factor gives the exact adapted transverse parameter and resolves the
   second blowup with poles `12,0` and degree-twelve residue.  Thus Case 1
   selects the same 23-component numerical package as Case 2.  The
   chain-to-boundary translation is complete, but the resulting package is
   allowed by every current intrinsic gate;
   see [`FRONTIER_LOG_SCALE_AUDIT.md`](FRONTIER_LOG_SCALE_AUDIT.md).
   The transformed pair has `[P,Q]=X^2`, so its actual ramification is
   `K+3H+div(X^2)`, not the boundary-supported `K+3H` representative.
   The corrected dicritical normal indices are `3` in Case 1 and `5` in
   Case 2; the total ramification intersection remains `35`.  The exact
   degree-twelve residue has a priori normalization-cover degree `1,2,4`,
   with image degrees `12,6,3`.  Case 1 retains all three rows, with
   degree-29 valuation contributions `3,6,12`.  Its formerly omitted bands
   `P:z^-6,...,z^-8` and `Q:z^-5,...,z^-12` are now reconstructed by the
   exact eight-layer continuation in
   [`CASE1_FULL_BAND_CONTINUATION.md`](CASE1_FULL_BAND_CONTINUATION.md); the
   remaining sieve input is the transported exact residue coefficient vector.
   In Case 2, general
   polynomial-right-component remainder ideals exclude cover degrees `2`
   and `4` already from `(J3),(J2)` and the determined part of `(J1)`.
   The remaining degree-twelve row is also empty: after localization at its
   forced endpoint `G_12 != 0`, the seven residual `(J1)` compatibility
   cubics generate the unit ideal, without `(J0)`.  Thus the former
   contribution `5` and remainder `24` are no longer live rows.  The earlier
   bottom-equation reduction, writing
   `H=gcd(C',G')`, `C'=H*c`, `G'=H*g` reduces `(J0),(J1)` to
   `B=K*c,F=K*g` and
   `2H(Ag-cE)+K^2(cg'-c'g)=0`, with `deg(K)<=deg(H)+1` and `K=0` or `t|K`.
   Its initial coefficients force `t|H`, eliminating gcd degree zero.  The
   maximal gcd degree `7` is also excluded: three coefficients of
   `remainder(G',C')` plus the terminal `(J0)` coefficient generate the unit
   ideal over the exact degree-35 field, with no residual `(J1)`
   compatibility equation.  The same extreme-stratum method excludes gcd
   degree `6` using a reconstructed linear cofactor, two terminal
   `G' mod H` coefficients, and `(J0)_{19}`.  Its remaining exact gcd
   degrees `1,...,5` are now retired as pre-compatibility strata.  The
   degree-one row is already reduced to the single origin pattern
   `ord(B),ord(E),ord(G'),ord(F)=(1,2,3,3)`; its only other valuation pattern
   is excluded by `(J1)`.  After those coefficient checks comes an accounting
   of the remaining divisorial valuation degrees.
2. **Glue the explicit F2 target row globally.**  The B0 replay in
   [`F2_75_125_DERIVATION.md`](F2_75_125_DERIVATION.md) remains an
   over-envelope rather than an exclusion.  The corrected exact tangent keeps
   every P-band direction (first-five kernel dimensions `6,6,7,7,10`); the
   layer-35 extra is a commuting `C0^4`, not `lambda*C0^(-1)`, whose formal
   resonance occurs only at layer 10 and is not a source kernel.  Nonlinear
   rows recover a source root through descent 7 and isolate the first local
   defect as `27*y^2-9*y+1=0` at descent 8.  Its fixed Kummer supports are
   excluded by the Q-band-one endpoint; only the nonzero double-root stratum
   of `R` remains, and it passes the first local target jet.  The complete B0
   recurrence now includes all `240` zero layers through `-200`; the fifth
   multiple is a layer-zero lower-tail Fitting row, not `E5=0`.  The complete
   target and layer-zero cokernels are now explicit rank-`14` modules, but
   old B0 generators span both before their earlier triangular equations are
   imposed.  The sparse lowest-`u` edge witness completes only as an infinite
   formal shear.  Its all-`r` polynomial second-order repair never terminates
   quadratically, and exact `r=3` unit ideals exclude cubic and quartic
   termination.  The exact `v^5` Kummer packet is nevertheless surjective
   over both descent-eight branches, and its edge Bezout recursion remains
   surjective through `v^10`.  Uniformly in `r`, its `18*r-1` minor pivots
   have strict-interior polynomial source lifts (`53` at `r=3`), so the
   current two-edge Newton data cannot remove them.  Exact CRT quotients the
   controllable `w=0` block at `r=3`, leaving a
   rank-`24` global Hermite module over the rank-two candidate algebra.  The
   fixed `w=1` block is now exactly eliminated by one normalized identity and
   a `10 x 10` determinant `75000`, leaving thirteen cokernel coordinates.
   The substitution has been carried into the pre-lower-tail system: it has
   `1,061` active source coordinates and degree at most eight, not thirteen
   scalar variables.  Tridiagonal unit minors remove `134` endpoint-disjoint
   new-Q coordinates on layers `39..29`, leaving `927`; the next live
   calculation is the coupled Schur/Fitting block beginning at descent `12`
   (layer `28`) and ending at descent `37`.  First-defect
   <!-- status-consumer: PF2ER1 64378dad616fc3f2 -->
   spacings `9..90` still require the same full
   target/tail descent if it is eliminated.  The retained
   contact audit in
   [`F2_BOUNDARY_HANDOFF.md`](F2_BOUNDARY_HANDOFF.md) proves that raw contact
   multiplicities cannot be promoted to boundary rows.  The later
   [`Kummer-orbit theorem`](F2_KUMMER_ORBIT_TRANSFER.md) supplies the missing
   rigidity: zero-root strata are impossible, simple cofactor roots are
   spectators, and the row has one principal chain or two copies on the
   nonzero double-root stratum.  The
   [`terminal residue theorem`](F2_TERMINAL_RESIDUE_COVER.md) computes source
   ray `(12,-17)`, target ray `(5,2)`, transverse index one, residue degree
   six, passport `(5,1)|(3,3)|(3,1,1,1)`, and geometric `A_6` monodromy.  The
   action is four-transitive and primitive, its target-fixed deck group is
   trivial, the
   cover has no `2`-by-`3` factorization, the transverse different is zero,
   and the same residue map occurs on every admissible principal stratum.
   The residue-different packet is `(4,2,2,2)`.  Three of those contributions,
   `(4,2,2)`, occur at three interior preimages of the target toric nodes and
   force boundary-attachment points; the last `2` occurs at the source
   endpoint over the smooth third branch value.  The rescaled cover is Belyi
   and its regular `A_6` closure has genus `25`.  The target valuation equality
   forces geometric degree at least six.  The global target coordinates have
   orders `(5,2)`, a uniformizer of order one, and a nonconstant residue
   coordinate, so both double-root packets necessarily dominate the unique
   extracted target divisor and force degree at least twelve.  Since this
   divisor is centered at infinity,
   no affine-sheet `+1` applies.  Purity instead forces a separate affine
   ramification row, and the global geometric monodromy has `A_6` as a simple
   composition factor.  The live global split is therefore one squarefree
   packet or two identical double-root packets over the same target component.  The
   [`global attachment compiler`](F2_75_125_GLOBAL_ATTACHMENT_COMPILER.md)
   now tracks the translated chart to the original nonmonomial valuation
   `nu(x),nu(y_old),nu(x*y_old^5-c)=(-25,5,12)`: six blowups extract the
   shared carrier `(-5,1)`, and six further blowups at each marked carrier
   point extract a principal arm.  It also completes the four-blowup target
   boundary with weights `(-2,-2,-1,-3,-2)`.  Each terminal source component
   has five forced boundary neighbors.  Exact transverse orders show that
   every interior neighbor needs two source blowups and is log-etale after
   resolution; the remaining smooth endpoint is the cyclic thickening
   `R/(w^3)` with zero normalization defect.  The one/two-packet principal
   source terminal-resolved skeletons have
   `(components,leaves)=(19,6)/(31,10)`, while the two live
   normalization equations are `d=6+rho_T` and `d=12+rho_T`.  Candidate mode
   rejects the former distinct-target row and runs the class-group, unit,
   canonical, purity, spectator, finite-normalization, and meridian gates once
   the missing global geometry is supplied.  The
   [`carrier Wronskian classifier`](F2_75_125_CARRIER_WRONSKIAN_CLASSIFIER.md)
   now extracts the deeper target ray `(5,36)`.  It forces the unique
   squarefree value `R=(v^2-3v+3)/25`, whose carrier map is the cyclic cubic
   `1+1/(v-1)^3`, so both simple spectators are unramified.  On the double
   row it forces `rho^2-3rho+1=0`; the degree-six carrier map is exactly the
   terminal Belyi map after a linear source change.  The
   [`carrier log-node theorem`](F2_CARRIER_LOG_NODE_PROFILE.md) now resolves
   the two squarefree spectator branches, the double-row fivefold point, and
   the common fan on every principal arm.  All marked exponent determinants
   are `1`, `3`, or `5`; the refined lower bounds are
   `(components,leaves)=(27,8)/(48,11)` before global purity.  The
   [`upstream extraction theorem`](F2_UPSTREAM_CARRIER_EXTRACTION_PROFILE.md)
   then proves that the carrier-zero ladder is unimodular and forces
   `coker(theta_log)=R/(W^3U^18)` at the extraction root, with branchwise
   quotient length `54`.  The cyclic blowup theorem combines this raw length
   with component self-intersections to give the stable Cartier charge `27`.
   The kernel-line theorem identifies the actual cyclic contribution as
   `deg(K_root)+27`, reducing the old twist ambiguity to the zero/pole divisor
   of one logarithmic kernel direction.  Since this packet is contracted, the
   kernel is the tautological line of a map to `P^1`; its nonnegative degree
   `e_root` changes the contribution to `27-e_root<=27`.  Full logarithmic
   divisibility of the fixed tangential coordinate proves `e_root=0`, so the
   cyclic root contribution is exactly `27`.  The outgoing-terminal theorem
   then closes the remaining terminal fan as unimodular log-etale.  The
   affine-purity frontier proves that a new component with `e>1` over an
   affine nonproperness curve is unavoidable, strengthening the component
   floors to `28/49`.  It bounds the geometric degrees by `6..9375` and
   `12..9375` and the target parametrization degree by `124`, but its coarse
   ledger survives at every remaining degree.  The target-curve theorem
   sharpens the curve input to 24 singular normalization charts `(3k,5k)`,
   `1<=k<=24`, with curve degree `5k<=120` and a proper
   divided-difference collision/critical ideal.  The live geometric work is
   now to compile the solved `k=1` collision quartic and its four generic
   nodal conductor points against the source pullback, extend the
   stratification to `k=2..24`, and compile possible cancellation centers.
   The `k=1` puncture is transverse on `(5,2)`; `lambda=125/729` decides
   whether it follows the special carrier centers, but the resolved terminal
   neighborhood supplies no affine-divisor slot and does not determine `e`.
   Exact complement monodromy now rules out a one-row completion on every
   immersed, distinct-image `k=1` collision partition and on the generic
   one- and two-cusp faces.  Each complement group is `Z`, so the fixed
   affine-sheet remainder forces a second ramified affine component.  The
   conditional component floors on those seven strata are `29/50`, and
   their Chern compiler must contain two affine packets.  The first escape
   is `E_6+A_1`: its noncyclic group admits a transitive degree-six
   fixed-sheet action, so its multiplicity-three boundary packet is now the
   sharp `k=1` target.  Severe image mergers and the other 23 charts retain
   the unconditional `28/49` floors.
   Locate the row at another unresolved boundary locus, then prove the global
   localized-Chern/descent identity.  The independent coefficient
   route still must impose the exact carrier points on the lower Laurent system.
   The degree pair `(75,125)` remains unexcluded.
   <!-- status-consumer: LCCT1 2fc6ecea7a7c8b49 -->
   <!-- status-consumer: LKGD1 8a357250b5005186 -->
   <!-- status-consumer: LTKT1 32ac27318f16c20c -->
   <!-- status-consumer: PF2OTT1 af25012e34020e11 -->
   <!-- status-consumer: PF2APF1 192055eb737d3140 -->
   <!-- status-consumer: PF2ATC1 9ab722c45c586b73 -->
   <!-- status-consumer: PF2K1C1 358a6ba820e8b2f1 -->
   <!-- status-consumer: PF2PPA1 b24c4d80c2f8230e -->
   <!-- status-consumer: PF2K1M1 fafcbb3c2e6ceb2b -->
   At the concentrated `E_8` endpoint, exact degree-six enumeration leaves
   one fixed-sheet permutation class, the exceptional `A_5` action.  Its
   peripheral action forces two distinct `(2,1)` source rows and excludes a
   single `(2,2)` row.  With total negativity `N=n_1+n_2`, maximal carrier
   contact must pass four scale-free equations and forces
   `4N>=36+s_X`; the minimal `N=2` packet instead leaves by `b<=4`, before
   the first transport-independent carrier equation at `b>=5`.  The next E8
   input must therefore fix a target transport globally or bound `N`.
   Exhausting the `A_5` coset actions gives the broader fixed-sheet list
   `d=6,10,15,30`, with `r=2,4,6,14` distinct `(2,1)` rows.  The uniform
   squarefree doubled residual `7d-62+4N-4r(b-6)-s_X` excludes only the
   minimal maximal-contact degree-six and degree-ten cases; degrees 15 and
   30 remain.  Thus the immediate global target is to bound `N` or couple
   these rows to the already compiled terminal/carrier filtration.
   <!-- status-consumer: PF2K1E8M1 bbb282c6bcfa62fc -->
   The universal simple-inertia quotient has order `240`, and its full
   subgroup atlas gives `13` fixed-sheet actions in degrees
   `6,10,12,15,20,24,30,40,60,120`.  Its central order-four action produces
   `(2,f)` rows with `f=1,2,4`; only the minimal maximal-contact degrees
   `6,10,12` fail the squarefree Chern ledger.  This reduces the remaining
   simple-inertia endpoint to a finite stable-charge/contact problem, while
   inertia greater than two remains a separate branch.
   <!-- status-consumer: PF2K1E8O1 4251750ed4e43c89 -->
   The complete-chain identity further removes raw negativity and contact:
   after the full determinant cycle is subtracted, only `u-1` point units
   remain, while the E8 cusp costs at least `2R`.  All simple-inertia rows
   fail.  The immediate endpoint is therefore the sign of the unresolved
   normalization/`Fitt_1` correction, with required deficits explicitly
   ranging from `3` to `113`.
   <!-- status-consumer: PF2K1CB1 5cc386dba344a867 -->
   Cyclic-submodule positivity proves that every isolated `Fitt_1` correction
   is an effective quotient, so the negative escape cannot occur.  This
   excludes the complete one-component simple-inertia E8 atlas.  Continue
   with E8 inertia greater than two, multiple affine components, and then
   the remaining `k=1`/higher-`k` cusp strata.
   <!-- status-consumer: LCSP1 8658eebeb1d65671 -->
   The
   [`common-power carrier theorem`](COMMON_POWER_CARRIER_WRONSKIAN.md)
   now makes this reusable: for every primitive degree-`k` common edge with
   coprime powers it forces the first nonshear descent, turns fixed-carrier
   reconstruction into one unique projective linear kernel, and assigns an
   explicit three-point passport to every root-multiplicity partition.  This
   supplies a dessin-first entry point for later frontier rows; the resonant
   `k=2` and imprimitive partitions remain separate cases.
   The next lower-system handoff is also exact: the
   [`carrier specialization compiler`](F2_75_125_CARRIER_SPECIALIZATIONS.md)
   routes the rational squarefree carrier outside the descent-eight
   double-root component and specializes every exposed linear map for the two
   double carriers over `Q(sqrt(5))`.  It pins the factor-quotient block of
   `53` coordinates and distinguishes it from the full `347`-coordinate
   Laurent cokernel.  The downstream
   [`nonlinear forcing compiler`](F2_75_125_NONLINEAR_FORCING.md) now carries
   all ten endpoint circuits and the upper tangent solutions into those
   rows, appends the final `7+6` functionals and five incidence rows, and
   routes the squarefree carrier through all later spacings `9..90`.  The
   remaining tasks are a localized ideal test for the double component and
   separate target-and-tail compilers for the squarefree ledger.
   <!-- status-consumer: PF2CS1 666da98d2d24669e -->
   <!-- status-consumer: PF2NF1 cfd1da5136c0b6d0 -->
   <!-- status-consumer: PCW1 94b10929118f151d -->
   <!-- status-consumer: PF2CW1 a7774b0fa736b64c -->
   <!-- status-consumer: PF2GA1 57dea3062b1147fb -->
   <!-- status-consumer: PF2LNP1 e4f0f231bf7494d5 -->
   <!-- status-consumer: PF2CLP1 41625dd5d3f8f898 -->
   <!-- status-consumer: PF2UCE1 7f15bc756cc73fff -->
   <!-- status-consumer: LCBBC1 b3eb4679f781c55f -->
   <!-- status-consumer: PF2GC1 6ba3fd9eb6a0bcdf -->
3. **Retire the excluded `(96,144)` repeated-tail branch.**  The source
   statements are now reconciled in
   [`FRONTIER_CLOSING_ATTACKS.md`](FRONTIER_CLOSING_ATTACKS.md).  The 2017
   table is a necessary over-approximation whose source leaves the
   lower-corner filter commented out; the 2016 full transition is stated only
   in a remark.  Reusing the proved common-tail part of the 2022 `(8,32)`
   calculation gives `q1=d0=4` and reduces the `(8,40)` vertical factor to a
   cubic.  Every cubic root partition containing a simple root yields the
   forbidden corner `(8,4)`.  The remaining factor
   `R=kappa*x^2*y^7*(y-lambda)^3` translates the initial edge to
   `(8,40)->(8,12)`.  Its maximum complete-chain length is three; even with a
   permissive superset of possible lower corners, the open-chain counts are
   `1,6,3,0` and no final corner occurs.  This repeated-tail row is excluded
   before Laurent compilation.  The other five `(96,144)` rows remain live.
4. **Test reuse only after a new band derivation exists.**  Derive the new
   curve polynomial `A`, give an explicit base/coordinate map, pass the exact
   cyclic-curve identity check, and then compare the forcing differential,
   primitive support, residual scaling, and gauge/exact class.  Equal coarse
   fingerprints alone are not evidence of reuse.
5. **Optional engineering simplification.**  Emit one machine-readable report
   per normal-form branch containing provenance, polygons, bands, bracket
   layers, Wronskian data, de Rham certificate, and any reuse certificate.

Scope guard: this handoff records a compiler programme, not a new family
exclusion or bound.  The externally published bound remains 108; the local
125 statement retains its explicit dependence on the published normal-form
reduction and local `(72,108)` certificate.

## Ranked worklist

“System” below means the coefficient system still to be derived; none is
silently treated as solved.

| Pair | gcd / ratio | first corner and chain data | Remaining obstruction | Expected difficulty | Does the \((72,108)\) proof apply unchanged? |
| --- | --- | --- | --- | --- | --- |
| \((75,125)\) | 25; \(3:5\) | family F2, \(A_0=(5,20)\), final \((7/5,2)\); certified terminal row `(e,f)=(1,6)` | eliminate the earliest movable-double-`R` full-tail Fitting branch, then the later-spacing regimes; in parallel globally test its two-principal-packet gluing | high: one exact earliest coefficient stratum plus a finite spacing ledger | no: different ratio, corner, bands, and approximate-root exponents |
| \((84,126)\) | 42; \(2:3\) | \(A_0=(7,35)\to(19/7,5)\), and a second realization \((12,30)\to(16/3,10)\to(11/6,3)\) | two chain realizations must both be excluded | high | no: the ratio agrees but neither chain is \((8,28)\to(11/4,7)\) |
| \((96,128)\) | 32; \(3:4\) | F24, \((8,24)\to(14/4,6)\to(5/4,0)\to(19/8,3)\) | longest family-derived corner chain at this maximum | very high | no |
| \((88,132)\) | 44; \(2:3\) | \((11,33)\to(19/4,8)\) | new \(2:3\) coefficient system | high | no |
| \((90,135)\) | 45; \(2:3\) | \((9,36)\to(17/9,4)\) with \((m,n)=(3,2),(2,3)\); \((9,36)\to(9,24)\to(11/3,8)\); and \((12,33)\to(11/3,8)\), both with \((2,3)\) | four published realizations; their Laurent polygons and coefficient systems remain to be derived | very high | no |
| \((56,140)\) | 28; \(2:5\) | F11, \((7,21)\to(13/7,3)\) | different approximate-root multiplicities | high | no |
| \((84,140)\) | 28; \(3:5\) | F9, \((7,21)\to(11/7,2)\) | family companion to the previously excluded \((56,84)\) | medium-to-high; prior \((7,21)\) work may help | no, but older \((7,21)\) reductions may be reusable |
| \((96,144)\) | 48; \(2:3\) | \((8,40)\to(8,28)\to(11/4,7)\), \((m,n)=(3,2)\); four \((2,3)\) chains from \((12,36)\), through respectively \((12,33)\to(11/3,8)\), \((9,24)\to(11/3,8)\), \((21/4,9)\to(19/4,8)\), and \((21/4,9)\to(12/4,5)\); and \((12,36)\to(12,30)\to(16/3,10)\to(11/6,3)\), \((m,n)=(3,2)\) | The repeated-tail row is excluded: simple-root partitions give `(8,4)`, and the triple-root translation `(8,40)->(8,12)` has no complete chain.  The other five rows still have no Laurent systems. | five high-difficulty rows remain | no |
| \((108,144)\) | 36; \(3:4\) | \((8,28)\to(7/4,3)\) | same first corner as the audited case but a different final corner and exponent ratio | high | no |
| \((42,147)\) | 21; \(2:7\) | F7, \((6,15)\to(7/3,4)\) | large degree ratio and approximate-root exponent | high | no |
| \((63,147)\) | 21; \(3:7\) | F8, \((6,15)\to(8/3,5)\) | companion chain to F7 | high | no |
| \((98,147)\) | 49; \(2:3\) | \((7,42)\to(13/7,6)\), both orientations | new \(2:3\) chain | high | no |
| \((100,150)\) | 50; \(2:3\) | \((10,40)\to(16/5,6)\to(23/10,3)\) and \((10,40)\to(18/5,8)\to(8/5,3)\), both with \((m,n)=(3,2)\) | two distinct Laurent/residual systems expected but not yet constructed | very high | no |

The table records every one of the 24 published chain realizations in this
window.  A “remaining system” is deliberately recorded as *not yet constructed*
when the cited paper supplies corner data rather than explicit Laurent
coefficient equations; the number of coefficient systems need not equal the
number of chains until the polygonal branching is derived.

## Does the method eliminate a family?

No family theorem follows from the reproduction.  The pair
\((72,108)=(2d,3d)\) with \(d=36\) is controlled not merely by its gcd or
ratio, but by the specific admissible tail

\[
 (8,28)\longrightarrow(11/4,7),
\]

the resulting two Laurent polygons, their band widths, and the degree-35
first-block field.  The frontier contains many other \(2:3\) pairs, and their
corner chains differ.  Even \((108,144)\), which shares \(A_0=(8,28)\), uses
ratio \(3:4\) and final corner \((7/4,3)\), so the coefficient system changes.

A defensible parameterized compiler conjecture is:

> A repeated admissible tail should permit reuse of its valuation/band descent
> after the earlier chain has been separately normalized.

The raw \((96,144)\) row through
\((8,40)\to(8,28)\to(11/4,7)\) initially looked like the first test.  A
source-level rereading changes that status.  In the 2016 lower-side paper,
the remark following Proposition 3.29 says that
\(B_0=(8,28),B_1=(8,40)\) leads to the impossible last lower corner
\((8,4)\) and can be discarded.  The 2017 Section 7 table nevertheless lists
the length-two complete chain.  The formal proposition proves the
impossibility of the relevant class of last lower corners; the application to
this particular `B_0,B_1` pair is stated in the remark as straightforward.

The reconciliation is now complete.  The common tail forces a cubic residual
factor.  Simple-root partitions give `(8,4)`; the sole triple-root partition
gives the translated edge `(8,40)->(8,12)`, whose permissive complete-chain
enumeration has counts `1,6,3,0` and no final corner.  The source-aware
compiler records the row as pre-excluded and refuses to manufacture a Laurent
block.  See [`NEWTON_DERHAM_COMPILER.md`](NEWTON_DERHAM_COMPILER.md) and
[`cas/complete_chain_no_escape.py`](cas/complete_chain_no_escape.py).

<!-- status-consumer: C1FBC1 0f14ef01fff25097 -->
