# Plane Keller search policy

This policy applies only to searches for two-dimensional Keller maps.

1. Do not spend coefficient-solving time on arbitrary counterexample searches
   with \(\max(\deg P,\deg Q)<125\).  The validated published reduction plus
   the locally replayed \((72,108)\) certificate excludes that range.
2. Low-degree work remains useful for automorphisms, regression examples,
   implementation tests, and validation of individual Newton/valuation
   lemmas.  Such runs must not be advertised as counterexample searches.
3. Every proposed plane counterexample must first be reduced to a standard
   \((m,n)\)-pair and pass the possible-corner/admissible-complete-chain filter.
4. Record orientation, first corner, all successor/final corners, gcd, degree
   ratio, and every vertex nonvanishing condition before constructing a
   coefficient ideal.
5. Before Gröbner elimination, supply exhaustive positive local branch
   scales (including nested resonance charts), compile them with the
   [log-boundary compiler](LOG_BOUNDARY_COMPILER.md), and run the chart-aware
   [boundary-lattice prefilter](BOUNDARY_LATTICE_PREFILTER.md).  Use every
   boundary prime, not only dicritical components; distinguish primitive
   divisor classes from pullback multiplicities; and declare whether the
   chart is `A^2` or a Laurent chart.  Inspect the compiled semigroup
   conductors, differents, and residue degrees.  If the candidate has a
   certified balanced squarefree retained-polynomial presentation with one
   omitted fierce boundary, supply the `retained_root_euler` block and reject
   `deg(A)>1` before normalization, Smith, or coefficient searches.  Do not
   apply this shortcut to nonsquarefree collisions or undeclared boundary
   presentations.  On a complete `A^2`
   resolution, also run the
   [intrinsic adjunction/Noether gate](INTRINSIC_A2_BOUNDARY_GATE.md), supply
   the global target pole vector, and require effective ordinary and log
   ramification plus an intrinsic dicritical.  A nonproper candidate must
   have canonical free depth at least three.  Corners alone are not a
   compiler input.  After a Puiseux/Kummer substitution with transformed
   bracket \(X^r\), record every band exponent modulo \(r+1\).  Apply a
   constant-Jacobian sparse theorem only to the trivial-character sector
   after proving descent through \(u=X^{r+1}/(r+1)\); a small raw Laurent
   support is not itself a Keller support.
6. Treat the source boundary tree as only the first layer of the obstruction.
   A closure claim must also resolve the target nonproper curve, give the
   harmonic map of source and target dual graphs with normal/residue degrees,
   and record the three-section linear series on every noncontracted source
   component.  Bare partitions of a field-degree remainder are not boundary
   packages.
7. If a Poisson-square truncation has the three-layer support
   `P=X^3A+X^2B`, `Q=X^2C+XD` with degree bounds `(3,4;2,3)`, replace its
   reduced coefficient system by the three classified components from
   [POISSON_SQUARE_RIGIDITY.md](POISSON_SQUARE_RIGIDITY.md): the forced
   tangent-pencil closure, `C=0`, and `A=0`.  Retain nonreduced structure
   explicitly when the downstream calculation is scheme-theoretic.  The
   exact associated-prime filtration has three minimal components, the
   three intersection branches, and two deeper core/intersection branches.
   Run proposed lower equations through
   [`cas/poisson_square_filtered_modules.py`](cas/poisson_square_filtered_modules.py)
   before constructing a new global coefficient ideal.  It reports whether
   each dense associated chart is preserved, cut, or eliminated and carries
   the certified `d3,d2` transverse Hilbert vectors needed by later
   scheme-theoretic checks.
8. Before a residue coefficient elimination, run the general polynomial
   right-component remainder sieve for every cover degree dividing the gcd
   of its coordinate degrees.  Do not replace a general right component by
   a parity or fixed-critical-point ansatz.  The sieve requires the complete
   residue coefficient vector.  If an archived Newton calculation omits lower
   bands, compile the branchwise ledger from
   [`CONDUCTOR_JET_TRUNCATION.md`](CONDUCTOR_JET_TRUNCATION.md).  Prefer its
   dependency-sensitive form: name every input and matching/residue output,
   certify a complete expression tree, propagate derivative, pole, and other
   contact losses only along paths that actually use that input, and derive
   available normal-jet orders from a certified valuation frontier of omitted
   Newton exponents and a normal-valuation vector.  On Laurent cones, require
   the frontier certificate to bound the entire omitted support, not merely
   its coordinatewise antichain.  Use the scalar branch maximum only as a
   backward-compatible fallback.  The archive is sufficient only when
   every displayed dependency reaches `conductor + path loss`.  Otherwise
   recover precisely the failing input band and exact deficit; never fill
   omitted bands with guessed zeros.
9. Any localization in a coefficient solve must ship with a complementary-
   strata audit.  A basis containing `1` is not an adequate artifact without
   the input generators, field, order, saturation factors, and an explicit
   identity or independently checkable resultant chain.
10. In parallel with degree-ordered coefficient searches, construct the
    [canonical finite-normalization package](FINITE_NORMALIZATION_PROGRAM.md).
    Its underlying module is automatically finite free in dimension two, so
    do not introduce a closed-point flatness gate.  Instead record every
    missing-boundary prime, its target curve, ramification and residue degrees,
    the affine sheets over the same curve, the boundary class-group basis, and
    the resolved projective log ledger.  Compile its coarse row through
    `cas/finite_normalization_signatures.py` before treating it as a new
    boundary type; use the resulting Pareto coordinates only as a bounded
    classification device, not as an existence claim.
    Before assembling its conductor quotient or local-cohomology residue, run
    the shared retained-root Euler gate whenever its hypotheses are certified;
    a finite-support residue cannot repair the global term
    `(deg(A)-1)[A^1]`.
    Do not pursue a purely local exclusion of singular unibranch packets:
    [`UNIBRANCH_SPECTATOR_COUNTERMODELS.md`](UNIBRANCH_SPECTATOR_COUNTERMODELS.md)
    gives smooth integral finite-free models in every rank at least four,
    with an étale spectator and saturated Euler budget.  The usable
    obstruction must retain the global distinguished `A2` open, its trivial
    unit group, the free boundary class group, and connected monodromy.
11. Make the full nodal logarithmic matrix the primary Plane-JC theorem
    target.  The logarithmic cotangent cokernel is pure Cohen--Macaulay, but a
    conductor mismatch lives naturally in degree-one local cohomology and
    cannot be lifted faithfully into its degree-zero torsion.  The complete
    determinant divisor and generic branch Smith profiles are insufficient:
    compile `Fitt_1`, the normalization defect, and the localized second
    Chern length at every node.  Type-I bracket nonvanishing alone is not a
    finite-support residue; the `(75,125)` terminal bracket normalizes to
    one.  Use
    [`LOG_CONDUCTOR_DEGREE_SHIFT.md`](LOG_CONDUCTOR_DEGREE_SHIFT.md) and
    [`F2_LOG_NODE_PROFILE.md`](F2_LOG_NODE_PROFILE.md) for the terminal
    profiles.  The
    [`F2 carrier profile`](F2_CARRIER_LOG_NODE_PROFILE.md) closes the marked
    carrier, aligned principal-arm, and spectator nodes as tame log-etale.
    The subsequent
    [`upstream extraction profile`](F2_UPSTREAM_CARRIER_EXTRACTION_PROFILE.md)
    proves that its carrier-zero ladder is also unimodular but forces
    `coker(theta_log)=R/(W^3*U^18)` at the extraction root, with length-`54`
    branchwise matching quotient.  Treat this as a degree-one class: the live
    target is a global localized-Chern/descent identity.  Raw node lengths are
    not invariant; use the
    [`cyclic boundary charge`](LOG_CYCLIC_BOUNDARY_BLOWUP_CONSERVATION.md)
    `D^2/2`, which has stable F2 root value `27`, and separately compile the
    kernel line from
    [`LOG_CYCLIC_COKERNEL_TWIST.md`](LOG_CYCLIC_COKERNEL_TWIST.md).  The actual
    contracted-packet theorem
    [`LOG_KERNEL_GAUSS_DEGREE.md`](LOG_KERNEL_GAUSS_DEGREE.md) makes this
    `27-e_root<=27`, where `e_root` is the nonnegative degree of the
    logarithmic kernel-direction map to `P^1`.  The
    [`tangential-kernel theorem`](LOG_TANGENTIAL_KERNEL_TRIVIALIZATION.md)
    proves `e_root=0`: the fixed covector `dz` pulls back into the full ideal
    of `3E+18L`.  Retain the exact cyclic contribution `27`; do not reopen an
    independent root twist or nilpotent-jet variable.  The
    [`outgoing-tail theorem`](F2_OUTGOING_TERMINAL_TAIL.md) also closes the
    remaining terminal rays as unimodular log-etale.  Compile noncyclic
    corrections and possible cancellation centers.  The
    [`affine-purity frontier`](F2_AFFINE_PURITY_FRONTIER.md) already proves
    that a new affine-branch component is necessary, raising the source
    floors to `28/49`, but also proves that generic ledgers survive for every
    `d` in `6..9375` or `12..9375`.  Do not rerun a free row enumeration;
    recover the actual nonproperness curve, factor its pullback, and locate
    the resulting proximity chain.
    The
    [`target-curve atlas`](F2_AFFINE_TARGET_CURVE_ATLAS.md) reduces the first
    operation to the 24 charts `(deg p,deg q)=(3k,5k)`, `1<=k<=24`, on the
    proper divided-difference collision/critical locus.  Start there; reject
    any unit collision ideal before implicitization.
    On `k=1`, use the exact quartic in
    [`F2_AFFINE_TARGET_K1_COLLISION.md`](F2_AFFINE_TARGET_K1_COLLISION.md):
    compile its four generic nodal roots and all discriminant, diagonal,
    tangent, and merged-image degenerations before attempting a generic
    two-variable elimination.
    Compare the leading residue with `125/729` to determine whether the
    target curve follows the special carrier cluster.  Do not use the
    terminal index `3` as the affine-row index: the terminal neighborhood is
    already resolved and admits no divisor dominating the affine curve.
    Locate the affine divisor at a different unresolved source-boundary locus.
    Use
    [`UNIVERSAL_COMPLETE_CHAIN_BOUNDARY_SATURATION.md`](UNIVERSAL_COMPLETE_CHAIN_BOUNDARY_SATURATION.md)
    as the claim boundary.
    <!-- status-consumer: UCBS1 824720a8f727bdf8 -->
    <!-- status-consumer: LCDS1 5b4d92acd50d6c41 -->
    <!-- status-consumer: PF2LNP1 e4f0f231bf7494d5 -->
    <!-- status-consumer: PF2CLP1 41625dd5d3f8f898 -->
    <!-- status-consumer: PF2UCE1 7f15bc756cc73fff -->
    <!-- status-consumer: LCBBC1 b3eb4679f781c55f -->
    <!-- status-consumer: LCCT1 2fc6ecea7a7c8b49 -->
    <!-- status-consumer: LKGD1 8a357250b5005186 -->
    <!-- status-consumer: LTKT1 32ac27318f16c20c -->
    <!-- status-consumer: PF2OTT1 af25012e34020e11 -->
    <!-- status-consumer: PF2APF1 192055eb737d3140 -->
    <!-- status-consumer: PF2ATC1 9ab722c45c586b73 -->
    <!-- status-consumer: PF2K1C1 358a6ba820e8b2f1 -->
    <!-- status-consumer: PF2PPA1 b24c4d80c2f8230e -->
12. Rank residual degree-frontier regression work by the tables in
    [FRONTIER_CLOSING_ATTACKS.md](FRONTIER_CLOSING_ATTACKS.md) and
   [NEXT_DEGREE_FRONTIER.md](NEXT_DEGREE_FRONTIER.md).
   The pair \((75,125)\) is the first numerical maximum.  Its selected F2
   chain now has a certified `(e,f)=(1,6)` terminal row, so its live priority
   is the finite global attachment ledger: one packet or two packets on the
   same target divisor, together with spectators, the target equation and
   full pullback factorization of the purity-forced affine branch component,
   three forced target-node attachment points, and the
   source-endpoint-over-smooth branch incidence.
   Exhaustive lower Laurent masks are an optional independent coefficient
   route, not a prerequisite for that gluing.  Multiple chains at
   \((84,126)\), \((90,135)\), and \((96,144)\) may still offer reusable
   structural tests.
    <!-- status-consumer: PF2GC1 6ba3fd9eb6a0bcdf -->
    <!-- status-consumer: PWB7 19f4f4ffc96227a3 -->
    <!-- status-consumer: CJT1 afb70f90ff10f3d7 -->

The frontier is a lower bound, not an attainability prediction.  Search
documentation must keep JC(2) separate from the repository's JC(3)
construction program.
