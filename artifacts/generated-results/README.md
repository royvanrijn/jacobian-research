# Generated artifacts

This directory contains reproducible outputs from the retained
three-dimensional work.  Only outputs with a named generator and verification
path below are treated as reference artifacts.

- `two_pair_sic_bidegree33_rank_two_interior_cyclic_split.json` records
  the exact characteristic-zero critical-algebra decomposition
  \(18=14+2+2\) at the first integral rank-two point.  The period function
  \(P\) is primitive on the interior and both endpoint summands, giving the
  cyclic interior basis \(1,P,\ldots,P^{13}\).  The checker also proves
  that the saturated toric logarithmic critical ideal equals this interior
  ideal, so its exact Picard--Fuchs target rank is \(14\).  Exact
  endpoint-idempotent
  period values are nevertheless nonzero, proving that an \(m\)-dependent
  divergence reduction is required.  It also certifies the first gradient
  lift of the degree-14 eliminant and records the nonzero endpoint
  restrictions of its 6791-term \(t\)-certificate.  Their ordinary normal
  forms occupy all \(14+2+2\) critical-algebra coordinates.  Regenerate it with
  `.venv/bin/python scripts/verify_two_pair_sic_bidegree33_rank_two_interior_cyclic_split.py`.
  Its whole-file SHA-256 is
  `c3cbe744b784e10b06e565276569eb9736209b068d6d41a4d3e5218513e088e9`.
  The canonical source is
  `extended-geometry/TWO_PAIR_SIC_BIDEGREE33_RANK_TWO_HOLONOMIC_PROBE.md`.
- `two_pair_sic_bidegree33_rank_two_dmodule_picard_fuchs_research.json`
  records the exact fixed-fiber rational \(D\)-module seed for the
  generating integrand \(u^2/(u^3-zQ)\).  Macaulay2 computes 34
  first-order annihilators of \((u^3-zQ)^{-1}\) and 76 annihilators of
  the specific numerator \(u^2\); both ideals are certified holonomic of
  rank one.  The stored run deliberately stops before pushforward, so it
  is an all-order integrand certificate rather than a Picard--Fuchs
  operator for the interval period.  Regenerate it with
  `.venv/bin/python scripts/research_two_pair_sic_bidegree33_rank_two_dmodule_picard_fuchs.py --annihilator-only`.
  Its whole-file SHA-256 is
  `6d570888a35b178adc4f5b67dc666094f55da374d90ba5fd607af56c0016605f`.
  The canonical source is
  `extended-geometry/TWO_PAIR_SIC_BIDEGREE33_RANK_TWO_HOLONOMIC_PROBE.md`.
- `two_pair_sic_bidegree33_rank_two_ore_gcd.json` records six exact
  modular shift-Ore Euclidean calculations comparing the sampled
  order-\(18\) and order-\(27\) recurrence operators.  Their common right
  factor has order \(14\), primitive coefficient degree \(58\), and is
  verified directly on 487 moment rows at each sample.  Its order matches
  the interior length in the exact \(2+2+14\) logarithmic-Jacobian
  decomposition, but no universal Picard--Fuchs identification is
  claimed.  Regenerate it with
  `.venv/bin/python scripts/verify_two_pair_sic_bidegree33_rank_two_ore_gcd.py`.
  Its whole-file SHA-256 is
  `0fb908186a2ac48a240a0f0c6a8928cea28a2ee6d36836caa80345fdc9f9f310`.
  The canonical source is
  `extended-geometry/TWO_PAIR_SIC_BIDEGREE33_RANK_TWO_HOLONOMIC_PROBE.md`.
- `two_pair_sic_bidegree44_rank_two_swap_slice.json` records the exact
  exclusion on
  \(F_P=\xi_1^4P(z_1,z_2)-\xi_2^4P(z_2,z_1)\).  The first five even
  moments have a length-twelve zero scheme with radical
  \(P=(z_2\pm z_1)^4\), where coefficient rank drops to one.  It includes
  the seven-prime lex reconstruction, exact rational tail reductions
  through order sixteen, the complete good-prime projective-boundary
  unit certificates, and the sharper three-moment parity-even exclusion.
  Regenerate it with
  `.venv/bin/python scripts/verify_two_pair_sic_bidegree44_rank_two_swap_slice.py`.
  Its whole-file SHA-256 is
  `8a953186d90d6b95622d60fa6a77b485ac88a37c18ff19cbf4eeff5100bddc35`.
  The canonical source is
  `extended-geometry/TWO_PAIR_SIC_BIDEGREE44_RANK_FRONTIER.md`.
- `two_pair_sic_frobenius_curvature.json` records the primitive
  first-order recurrences for nine SIC2C4 radial propagations, direct
  good-prime recurrence-curvature norms, differential \(p\)-curvature of
  the normalized beta period at every odd prime through \(101\), and the
  bounded prime-power valuation/re-entry correlation audit.  The exact
  good-prime recurrence curvature is \(d^d(m^p-m)^d\); the differential
  computation is bounded evidence, not an all-prime proof.  Regenerate it
  with
  `.venv/bin/python scripts/research_two_pair_sic_frobenius_curvature.py`.
  Its whole-file SHA-256 is
  `2a1ed9dff308448f6e1cb856e06a661637c3cefbde9d0870c4e4f5d585a06682`.
  The written analysis is
  `extended-geometry/TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md`.
- `keller_tschirnhaus_descent_567.json` records the exact ranks five, six,
  and seven comparison between the primitive coordinates `r` and `r+r^2`
  on the same split finite-etale algebra.  It pins the inverse interpolation,
  nonzero projective residuals, compiled universal-map seed parameters,
  stable boundary fingerprints, and relative/promoted coordinate degrees.
  It is generated and replayed by
  `scripts/verify_keller_tschirnhaus_descent_567.py`; pass `--write` only
  when intentionally refreshing the pinned certificate.  Its whole-file
  SHA-256 is
  `4177843a58dc77ce89aee5e2e5a2c781afc6d29ffb8bb17b013fd1bf0150ffda`.
- `degree_five_cubic_h7_unit_certificate.sing` is the exact
  characteristic-zero ten-variable consistency ideal for the full reduced
  cubic \(\mathbb A^{27}\) fifth-order lift component.  Its 401 generators
  contain 27 nonzero constants.  The selected \(X^{18}\) residual has an
  explicit inverse in
  `Q[a]/(94*a^3+335*a^2+400*a+160)`, and the program verifies both that
  direct identity and a one-generator degree-zero Singular lift of `1`.
  It is generated by
  `scripts/analyze_degree_five_cubic_fifth_order.py --exact-cubic
  --seventh-component-elimination --seventh-component-program-output
  artifacts/generated-results/degree_five_cubic_h7_unit_certificate.sing`
  and replayed by
  `scripts/verify_degree_five_cubic_h7_unit_certificate.py`.  Its
  whole-file SHA-256 is
  `86eeadee714614dba8794eb392d087e3bcedcb53ce51517b0d83acda8200e980`.
- `two_pair_sic_bidegree33_corrected_boundary_deepest.json` records the
  content-preserving branch exports through \(\mu_{10}\), the known rank
  \(6,6,5,5,5\) fiber shapes on the generic, \(L\), \(Q\), \(J\), and
  \(L=Q\) strata, the exact rank-fifteen
  \((\mu_3,\mu_4,\mu_5)\) algebra on \(L=Q=0\), primitive exact
  fraction-free normal forms for \(\mu_6,\mu_7\) there, the imported exact
  five-coordinate \(\mu_6,\mu_7\) normal forms on \(J=0\), and the exact
  unit ideals through \(\mu_{10}\) on all five \(t_0=0\) strata.  Those
  saturated strata partition the adapted \((L,Q)\)-plane, giving a complete
  branchwise exclusion of the \(t_0=0\) divisor.  On \(t_0\ne0\), it also
  records the exact \(t_0=1,\ u=s_0^{-1}\) presentation in which \(\mu_3\)
  is a fiber-independent base equation, \((\mu_4,\mu_5)\) gives a
  rank-six algebra in \(s_6,s_5\), and \(\mu_6\) reduces to all six
  standard monomials; its leading coefficients factor as
  \(K,H,Q_*KJ_*H\).  It also records
  \(K=4A_*-15Q_*\), \(H=4J_*-15A_*Q_*\), the changed rank-six basis on
  the generic \(K=0\) divisor, and a second exact rank-six basis on the
  reduced \(K=H=0\) linear slice.  A dense rational parametrization of
  \(H=0\) gives a third exact rank-six basis with leading ideal
  \((s_6^2,s_5^3)\); its omitted point belongs to \(J=0\).  It is
  generated by
  `scripts/verify_two_pair_sic_bidegree33_corrected_boundary.py`.  The
  \(t_0\)-open common-root equations and exceptional-factor radicals remain
  open; the artifact proves neither a semistable point nor a full
  characteristic-zero radical equality.  Its whole-file SHA-256 is
  `96881a3428d1b0a1cd279208ab13f9f235b41cc371acc5137606f502a3d0e19d`.
- `two_pair_sic_bidegree33_t0_open_fixed_fiber.json` records the exact
  rational base (5.12u) on the \(t_0\)-open rank-six chart.  It verifies
  \(\mu_3=0\), \(Q_*J_*KH\ne0\), quotient length six for
  \((\mu_4,\mu_5)\subset\mathbb Q[s_6,s_5]\), and the unit ideal after
  adjoining \(\mu_6\).  Thus the first multiplication norm is not
  identically zero on the local base component; the artifact does not
  compute its exceptional divisor.  It is generated by
  `scripts/verify_two_pair_sic_bidegree33_corrected_boundary.py
  --branch generic --prime 0 --t0-open-fixed-fiber`.  Its whole-file
  SHA-256 is
  `70e1fb7f2d33b3f471de75fbfacf37715594cec67bce1708d0d77035788eebf3`.
- `two_pair_sic_bidegree33_t0_open_curve_norm.json` records the exact
  rational \(\mu_3=0\) curve (5.12v), the primitive degree-\(198\)
  irreducible numerator and degree-\(144\) denominator of
  \(\det M_{\mu_6}\), and the next coefficient of
  \(\det(M_{\mu_6}+zM_{\mu_7})\).  That coefficient has numerator and
  denominator degrees \(209,153\), and its numerator is coprime to the
  \(\mu_6\) norm numerator.  Hence the degree-\(198\) exceptional divisor
  has no common \(\mu_6,\mu_7\) root on the curve's border open.  The
  denominator factors have degrees \(2,3,4\): exact calculations on the
  \(Q_*=0\) and \(J_*=0\) number-field fibers give length five and the
  unit ideal through \(\mu_7\), while the cubic factor is the curve pole
  and has coprime numerator.  Thus every defined point of the rational
  curve is excluded.  The artifact is generated by
  `scripts/verify_two_pair_sic_bidegree33_corrected_boundary.py
  --branch generic --prime 0 --t0-open-curve-norm`.  Its whole-file
  SHA-256 is
  `6ea0bdbdc489530a9d74479966feab4b8e2b066cd6990f3fb2a494792a9c73f9`.
- `two_pair_sic_bidegree33_t0_fitting_degree_scout.json` aggregates
  fifteen bounded directional reconstructions of
  \(\det M_{\mu_6}\) and the linear coefficient of
  \(\det(M_{\mu_6}+zM_{\mu_7})\) modulo the quadratic \(\mu_3(s_3)\).
  Across two base points and primes \(1019,2039\), 6750 paired-root
  samples fit the line functions and 750 held-out pairs verify them.
  Every direction matches the denominator models
  \(a_2^{41}Q_*^3J_*^3\) and \(a_2^{42}Q_*^4J_*^4\).  The artifact is a
  modular degree scout, not a multivariate Fitting certificate.  It is
  replayed by
  `scripts/verify_two_pair_sic_bidegree33_t0_fitting_degree_scout.py`.
  Its whole-file SHA-256 is
  `e6aec153eee2da72bae801efdb85da3ffb2b6c689d86a3eaba304f294cb2e2c6`.
- `two_pair_sic_bidegree33_t0_pencil_random_scout.json` aggregates
  forty-four deterministic bounded random shards of the complete
  determinant pencil \(\det(M_{\mu_6}+zM_{\mu_7})\) at primes
  \(43,47,59,71\).  The shards test 19800 paired bases, hence 39600 roots
  of the quadratic \(\mu_3(s_3)\).  Twenty points make all seven pencil
  coefficients vanish; direct Gröbner replay gives a reduced length-one
  common quotient through \(\mu_7\) at each.  Nineteen use the leading
  \(M_6\) pivot chart and one uses a second pivot; \(M_8\) restores full
  block rank at all twenty.  The aggregate also includes direct
  specialization-safe \(Q,J,K,H\) scouts at \(p=43\), 900 roots per
  divisor.  \(Q,J\) each contain one sampled reduced through-\(\mu_7\)
  point and both are excluded by \(\mu_8\); \(K,H\) contain none in the
  sample.  This is bounded modular evidence, not a global or divisor
  common-root exclusion.  The summary and shard hashes are checked by
  `scripts/verify_two_pair_sic_bidegree33_t0_pencil_random_scout.py`.
  Its whole-file SHA-256 is
  `31def59d50b689efb7001aea5b8731fc72a079d377bd25f775eda363a0cbee68`.
- `two_pair_sic_bidegree33_t0_strata_rank_continuation.json` aggregates
  the guarded lower-stratum and rank-complement continuation at \(p=43\).
  Thirty-one direct-stratum artifacts contain seventeen sampled reduced
  common roots through \(\mu_7\), all excluded by direct corrected
  \(\mu_8\) evaluation.  Twelve rank shards test 6300 further
  \(\mu_3\)-roots; four miss both selected pivots but have full joint
  \(M_6,M_7\) rank, and no rank-at-most-four point occurs.  Nine
  `liftstd` calculations extract the exact modular leading-coefficient
  borders and check them against every retained fiber.  On \(Q,J,JK\)
  their sampled zeros are exactly the quotient-length drops; the enlarged
  \(Q,J\) clouds have no extra polynomial relation through degree four.
  Five factored border resultants give the residual degree-\(20\) and
  degree-\(24\) components in (5.12aa), and the linear
  pseudo-remainder \(A s_3+B\) is coprime to every residual factor.
  The \(Q\)-row is also promoted over \(\mathbb Q\): its irreducible
  degree-\(36\) border has 588 terms, and its exact resultant factors as
  \(c\,u^{20}J_Q^4R_{20}^2\), with irreducible degree-\(20\), 200-term
  residual factor and a coprime dense linear pivot.  The other projected
  rows remain finite-field calculations, and no residual component is
  yet excluded through later moments.  The summary and all constituent
  hashes are checked by
  `scripts/verify_two_pair_sic_bidegree33_t0_strata_rank_continuation.py`.
  Its whole-file SHA-256 is
  `c8c35d8ff3749ccb2c866dbc0f4319bf1b475f66e063509eb6f91a573422f0af`.
- `two_pair_sic_bidegree33_boundary_generic_quotient.json` records the
  characteristic-zero generic and \(L,Q,J\) divisor quotient certificates,
  including the exact quadratic-extension \(\mu_6,\mu_7\) normal forms on
  \(J=0\), together with the characteristic-\(47,101\) replays.  It is
  generated by
  `scripts/verify_two_pair_sic_bidegree33_boundary_generic_quotient.py`.
  Its whole-file SHA-256 is
  `3b833ef7ad49555990a62c4ddfc1c7008f79377cd4df628e102425f99f908b48`.
- `backward_cubic_reduction_calibration.json` records the exact fixed-level
  restriction `G20 -> M19`, the stable factorization
  `M19=A_B o (F13 x I_6) o S_gamma`, survival of the rational collision pair,
  the separate direct-cubic and cubic-homogeneous terminal keys, a
  non-coordinate fixed-covector determinant regression, and a generic
  companion determinant identity.  It also proves that every nonzero
  homogenizing slice is a scaled copy of the base map, that the zero slice
  is triangular and injective, and transports the MacFarlane collision
  exactly to `tau=2`.  It is generated by
  `scripts/verify_backward_cubic_reduction.py`; no twelve-variable map is
  claimed.
- `backward_cubic_current_applications.json` applies those keys to thirteen
  retained Pareto representatives from ten active restricted-minima
  archives.  It exactly reconstructs the best retained direct source and
  the best retained homogeneous quotient and checks their surviving
  collision pairs.  It also records the then-current external-certificate
  consequences `n_cub<=20` and `n_HN,4<=40`; the later coordinate-pair
  reduction improves these to `n_cub<=19` and `n_HN,4<=38`.  The archive was selected under different
  objectives, so this is not a lower bound.  It is generated by
  `scripts/verify_backward_cubic_current_applications.py`.
- `macfarlane_g20_dimension_reduction_audit.json` records the exact
  MacFarlane `F13/G20` collisions, cubic-output rank six, zero constant input
  kernels, full nonlinear output spans, and the linear Keller-hyperplane
  obstruction.  It is generated by
  `scripts/audit_macfarlane_g20_dimension_reduction.py`.
- `macfarlane_f13_low_degree_invariants.json` records the exact torus grading
  and good-prime full-column-rank certificates proving that the pullback-fixed
  polynomials through degree three are only constants.  It is generated by
  `scripts/audit_macfarlane_f13_low_degree_invariants.py`.
- `macfarlane_f12_coordinate_pair_reduction.json` records the exact
  source-coordinate/target-square reduction from 13 to 12 variables, its
  direct determinant and rational collision, cubic-output rank six, and the
  resulting 19-variable cubic-homogeneous parent.
- `k12_coordinate_pair_frontier.json` classifies every linear target
  coordinate whose pullback is a polynomial graph coordinate and records
  exact unit-ideal obstructions to a raw degree-three restriction.  For the
  literal triangular coordinates it also excludes one target-shear stage of
  degree at most three, and excludes the two closest deletions through target
  degree four.  It is generated by
  `scripts/audit_k12_coordinate_pair_frontier.py`.
- `k12_parameterized_completion_frontier.json` records fixed-minor and
  unit-ideal certificates over all six quadratic graph-coordinate families
  of the twelve-variable map. It excludes quadratic target completion for
  every family and cubic target completion for all five single-defect
  families. It is generated by
  `scripts/audit_k12_parameterized_completion.py`.
- `k12_z8_cubic_completion_frontier.json` records the sparse minor-first
  closure of the remaining multi-defect quadratic graph family at target
  degree three. Three exact determinant opens generate the unit ideal, and
  every augmented determinant is `9/7` times its column determinant. It is
  generated by `scripts/audit_k12_z8_cubic_completion.py`.
- `k12_single_defect_quartic_completion_frontier.json` records constant
  full-column and augmented minors excluding target completion through
  degree four on the five single-defect quadratic graph families. It is
  generated by `scripts/audit_k12_single_defect_quartic_completion.py`.
- `hvc38_cross_construction_frontier.json` compares the public
  eleven-variable/rank-seven and local twelve-variable/rank-six routes to
  dimension 38.  It records seven inconsistent quadratic completed-pivot
  systems on the public lift, a quadratic-completion obstruction for the
  nonlinear local `z8` pivot, and a fourteen-dimensional coordinated
  degree-preserving local source-shear kernel whose exact rank-at-most-five
  Schur system is inconsistent.  It is generated by
  `scripts/audit_hvc38_cross_construction_frontier.py`.
- `hvc38_gap_closure.json` excludes the public nonlinear `d` pivot and the
  local nonlinear `z8` pivot through target degree eight.  It also records
  the 36-dimensional exact high-degree kernel of a 932-column quadratic
  source-target system, its seventeen integrable triangular directions,
  and the unit Gröbner basis excluding cubic rank five on their combined
  degree-three locus.  It is generated by
  `scripts/audit_hvc38_gap_closure.py`.
- `hvc38_maximal_block_closure.json` enumerates all six maximal jointly
  affine source blocks of `K12`, records the exact rational lifts of their
  complete quadratic source-target high-degree kernels, and gives fixed
  Schur rank-six witnesses.  It integrates every kernel direction into a
  full triangular family and records Singular unit ideals for the
  degree-three equations plus pinned rank-drop minor packets.  It is
  generated by `scripts/audit_hvc38_maximal_block_closure.py`.
- `automatic_missing_invariants_d3_d6.json` records the refined
  apolar-even/odd invariant dimensions through polynomial degree six for
  `V_d=End(Sym^d)`, `d=3,4,5,6`, together with the dimensions left after
  quotienting by moment monomials.  It is generated by
  `scripts/research_completed_moment_algebra.py`.  The same run stores exact
  modular moment-Jacobian ranks, bounded Hilbert necessary tests, and the
  Casimir values on the propagated all-order witnesses.  It also records
  the exact radial Casimir multipliers, the closed odd-cubic dimension
  formula, a beta-sum regression through `m=32` for the all-power
  Casimir-ladder proof, and sparse exact projections of `F_4^m` through
  `m=12`.  The theorem gives first nonzero quadratic
  `q_(2 ceil(m/2))` for every positive `m`; the bounded projections also
  store the stronger exact torus-phase supports.  It proves
  that the first missing space is the even quadratic completion of dimension
  `d-1`.  The full phase-support formula remains conjectural, and the scan
  does not classify every semistable component or prove the candidate
  augmented systems finite.
- `two_pair_sic_bidegree33_casimir_fiber.json` joins the corrected cubic
  moment degrees to the completed invariant algebra.  It records exact
  Hilbert-coefficient tests and good-prime rank-thirteen certificates for
  the corrected moment-only system, the smaller one-`q_2` system, and
  mixed `q_2`/`q_4` systems of the same total invariant degree.  Its
  complete weight-14 evaluation matrices prove that `mu_14` is independent
  from the pure Casimir span modulo the lower-moment span, both with `q_2`
  alone and with `q_2,q_4`.  The optional null-quadratic run records the
  exact seven-column linear normal-symbol matrix, its generic rank three,
  and the irreducible common cubic divisor `P` of its nonzero maximal
  minors.  Modulo `32003`, the rank-drop ideal after removing `P` has
  dimension three and multiplicity four in the five allowed coordinates.
  At allowed-coordinate base `(20,27,36,47,60)`, the full normal
  restriction of `mu_2,...,mu_12` is zero-dimensional of quotient length
  `195`; its recorded coordinate-power memberships certify the finite
  modular fiber.  Good reduction yields characteristic-zero transverse
  isolation at that base, hence on a nonempty open synchronized locus.
  The exceptional run decomposes the quotient-minor support into two
  quadratic-field components and one lower rational locus, proves all
  three are disjoint from `P=0`, and records good reductions of exact
  algebraic points on all four exceptional strata.  The normal fibers
  on `P=0` and the two top-dimensional residual components have length
  `195`; the lower residual locus has length `197`.  All are
  zero-dimensional, proving transverse isolation on a nonempty
  characteristic-zero open subset of each stratum.
  The default generator does not run the older direct boundary
  standard-basis probes.  These results do not prove global nullcone
  equality, cover every proper closed subset inside the exceptional
  strata, handle the `F_2=0` chart, or establish a smaller minimal moment
  cutoff.
- `lr_rees_sagbi_module_computation.json` records the finite target-invariant
  SAGBI basis, target-field modules and initial lifts in weights
  `p=+-1,+-2`, the saturated normal quotient, and the complete quadratic
  matrices for the only surviving weights `p=1,2`.  It is generated by
  `scripts/compute_lr_rees_sagbi_modules.py`.  Exact module reduction leaves
  the new `p=2` class
  `II_(F,2,-2)(partial_A,A^2 partial_A)` with remainder
  `-987/395*e_C` modulo the full `p=1` image.  The dependency-free
  `scripts/audit_lr_rees_sagbi_module_certificate.py` independently checks
  the separating covector `(0,-144/79,1)` at `(u,gamma)=(1/6,0)`.  The
  new quotient has exact annihilator `(gamma,6*u-1)` and length one.  Its
  whole-file SHA-256 is
  `c2f2700b051afd539d42b1fda50ff057f66a1f61297b5d49b0d064edf42732a8`.
- `lr_rooted_tree_normal_classes.json` records the exact constant-direction
  pre-Lie rooted-tree compiler for the \(F_2\) target lift and the family
  `tau_2=B(C)`, `tau_3=A(C(C))`, `tau_(n+2)=B(C(tau_n))`.  Its specialized
  `3 x 3` transfer matrix has a positive-coefficient Cayley--Hamilton
  recurrence proving a nonzero third saturated normal residue in every order
  `n>=2`; the stored expansions through order twelve are regression data, not
  the basis of the all-order claim.  The dependency-free
  `scripts/audit_lr_rooted_tree_normal_classes.py` replays the rational matrix
  identity and sign induction.  This certifies individual tree classes, not
  noncancellation in the full BCH/LR forcing sum.  Its whole-file SHA-256 is
  `1672290d4d326af04dd85248a0afa346b66d6b2e0d807a3cf2d2d2fa28936859`.
- `lr_mixed_bch_classes.json` records the exact balanced mixed-BCH
  continuation with `X=N*(x,0,-3z)` and the commuting lifted directions
  `D_B,D_C`.  It stores the complete two-component boundary-face operator and
  the triangular recurrence
  `c_(k+1)=-73440*(k+3)*(2k+7)*c_k`.  Together with the nonzero even Bernoulli
  coefficients, this proves a nonzero third saturated normal residue in the
  linear-in-`X`, bidegree-`(k,k)` BCH coefficient for every odd order
  `2k+1`.  The certificate also records the exact status boundary: target
  amplitudes multiply the class by `s^k*t^k`, so it vanishes on `s*t=0` and
  is not a universal lower-jet obstruction by itself.  Generate it with
  `scripts/compile_lr_mixed_bch_classes.py` and replay its rational
  coefficients with `scripts/audit_lr_mixed_bch_classes.py`.  Its whole-file
  SHA-256 is
  `877560240a592dcf6813948bd71ce95597e43d32d3b8b98bd2c1d9298ffe5740`.
- `degree_four_tau_even_parameters.json` records the exact modular
  Jacobian certificate for twenty-two apolar-even trace parameters of
  degrees `1,2^4,3^9,4^8`, together with the equal generic cotangent span
  of the first twenty-two moments.  It is generated by
  `scripts/verify_degree_four_tau_even_parameters.py`.
- `degree_four_moment_field_bounded_relations.json` records the bounded
  full-column-rank search excluding
  `Q(mu)+c_234^2*P(mu)=0` through invariant weight sixteen.  It is
  generated by `scripts/research_degree_four_moment_field.py`; this is not
  a proof of fixed-field equality.
- `degree_four_diagonal_moment_field.json` records the characteristic-zero
  parameter quotient of length 120 and the exact two-point
  reversal-related first-six-moment fiber on the diagonal quartic slice.
  It is generated by
  `scripts/verify_degree_four_diagonal_moment_field.py` and proves generic
  degree two on the raw parameter field; Weyl reversal identifies the two
  points in the invariant quotient.
- `degree_four_single_phase_moment_fields.json` records exact
  characteristic-zero parameter quotients of length 360 and two-point
  reversal-related first-seven-moment fibers on all ten coordinate
  single-phase parameter spaces.  It is generated by
  `scripts/verify_degree_four_single_phase_moment_fields.py`.  It proves
  generic degree two on the raw parameter field for a nonempty
  Zariski-open family of direction pairs in each nonzero phase.  Only the
  four cross-direction coordinate pairs in phases one and two are
  certified to move off the apolar-fixed quotient locus.
- `completed_moment_single_phase_fields.json` records exact finite-origin
  and two-point-fiber certificates on one single-phase quotient slice in
  every nonzero phase for `d=3,5`.  It is generated by
  `scripts/verify_completed_moment_single_phase_fields.py`.  Finite-field
  standard bases at `32003`, weighted-projective properness, and Nakayama
  prove characteristic-zero integrality and exact degree two on each raw
  slice.  The cross-direction tests in quintic phases one and two have
  nonzero `c_234` and therefore detect two distinct invariant-quotient
  orientations; the other slices are fixed-locus controls.
- `degree_four_q2_augmented_nullcone_local.json` records the first exact
  normal-jet calculation for the global quartic `q2`-augmented origin
  problem.  It is generated by
  `scripts/research_degree_four_q2_augmented_nullcone.py`.  On one
  normalized nonzero-`Sym^2` synchronized branch point, four linear
  pivots reduce the twelve forbidden directions to eight; the quadratic
  and cubic jet ideals then have dimensions six and four over `F_32003`.
  This is a bounded local frontier, not a global nullcone certificate.
- `degree_four_q2_cubic_decomposition.json` decomposes that cubic
  normal-jet support.  It is generated by
  `scripts/research_degree_four_q2_cubic_decomposition.py`.  The artifact
  records the exact radical and binary-cubic factorization on `x6=0`,
  the degree-nine generic fiber of the off-`x6` saturation, and the
  collapse of both dominant cubic sheets to one three-plane under the
  quartic moment jet.  The off-axis quartic calculation remains a
  timeout, so this is still a bounded finite-field normal-slice result.
- `degree_four_phase_one_chart_modular.json` records the exact
  \(F_{101}\) reduced four-point first-ten-moment fiber on the
  eight-dimensional two-direction phase-one chart, agreement through
  moment eleven, the odd-cubic sign split, and an explicit
  `SL_2` orbit witness pairing the extra branch with the apolar reversal.
  It also verifies the exact rational branch with \(u=5/3,w=6\), odd
  cubic values \(1728,-1728\), and its characteristic-zero orbit matrix.
  The Jacobian of the first eight moments is nonzero at all four
  reconstructed rational points, proving that those branches are reduced
  and isolated.
  It is generated by `scripts/research_degree_four_phase_one_chart.py`.
  Exclusion of additional characteristic-zero components remains open,
  so this is not a characteristic-zero generic-degree theorem.
- `arithmetic_keller_quintic.json` is the portable proof object generated
  from `arithmetic/specifications/ramified_quintic.json`.  It records local
  stability radii, coefficient CRT residues, real isolation, irreducibility
  and local-action witnesses, the selected translation and target, the
  inverse identity, and the expanded determinant-one map hash.  It is
  generated by `scripts/compile_arithmetic_keller_certificate.py` and replayed
  independently by the standard-library Python and PARI/GP verifiers.  Its
  current whole-file SHA-256 is
  `2ce5ec3dcccb09d355bae2c57004ec981cdd255b1c0906971632b309975cc4c2`;
  the generated Lean specialization has whole-file SHA-256
  `3015fa3b156058a82f1916b6abbd214c898f9ca50e0e20f67e9d649e2ce9590a`.
- `arithmetic_keller_quintic_stable_m2.json` is generated from the same
  ramified-quintic specification with `--stable-parameter 2`.  It retains
  the local algebra, inverse polynomial, and target while recording Fitting
  Newton area `13`.  Both independent verifiers reconstruct its shifted map.
  Its whole-file SHA-256 is
  `b8f1f851884bd351f00d906b85f749489468a1650fd5ccd4a1ef66f6190b9284`;
  its generated Lean specialization has whole-file SHA-256
  `7d25390a400e49816754e2f8136a5ee32b486a4f88fb64f0c8f77aa82582d4ba`.
- `arithmetic_keller_cubic_stable_n7.json` is generated from
  `arithmetic/specifications/connected_cubic_stable_n7.json`.  It records the
  connected cubic `T^3-T-1`, exponent `n=7`, and boundary-component count
  `7`; both independent verifiers reconstruct its cubic lift.  Its whole-file
  SHA-256 is
  `f58d0f463f6e1d8d655f03027a32fdffd41ad2aae5e7593a042bad4772e24f43`;
  its generated Lean specialization has whole-file SHA-256
  `c26ebd434fa214f6c4cd3e67efe405345b5f98549afa3d872501f4403789d142`.
- `fixed_quintic_certificate_ledger.json` is the generated ten-row ledger for
  the split-seed quintic arithmetic zoo.  It records every target, real-root
  count, Galois or factorization witness, the seven modulo-`7` partitions,
  and the coefficient Jacobian.  It is generated and checked by
  `scripts/verify_fixed_quintic_certificate_ledger.py`.
- `universal_quintic_target_line_search.json` records the bounded incidence
  search through the exceptional targets of projective height at most `30`
  and rational line parameters of height at most `120`.  It is generated by
  `scripts/search_universal_quintic_target_lines.py`; the negative result is
  bounded computation, and PARI is used only to classify screen survivors.
- `hc4_finite_field_sparse_search.json` records the collision-normalized
  degree-`5` through degree-`8` search over `F_11` and `F_13`.  It exhausts
  45,181,194 potentials supported on at most two vectors of the full linear
  collision-kernel basis and finds no exact constant-Hessian candidate.  It
  is generated by `scripts/search_hc4_finite_field_potentials.py`; the result
  is a bounded experiment, not an unrestricted `HC_4` theorem.
- `hc4_finite_field_dense_support_search.json` records 192 exact modular
  coefficient-ideal calculations: 96 deterministic supports of sizes
  `6,8,10,12`, each over `F_11` and `F_13`.  All reduce to the unit ideal,
  with no timeout.  It is generated by
  `scripts/search_hc4_finite_field_dense_supports.py`; the supports are
  sampled rather than exhaustive.
- `hc4_finite_field_axis_support_search.json` records 64 further exact
  modular coefficient-ideal calculations.  Its 32 supports preferentially
  include the directions capable of changing both forced determinant defects
  on the normalized collision axis.  Every ideal over `F_11` and `F_13` is
  again a unit ideal.
- `hc4_finite_field_cone_bridge_search.json` records 256 exact modular
  coefficient-ideal calculations on 128 further supports.  Each top
  homogeneous correction omits `x2` or `x3`, making its Hessian determinant
  identically zero, while lower-degree monomials bridge through the omitted
  variable.  Every full ideal over `F_11` and `F_13` is nevertheless a unit
  ideal.
- `hc4_finite_field_oblique_cone_bridge_search.json` records 288 exact
  modular coefficient-ideal calculations on 144 non-coordinate families.
  They use `u=x2+lambda*x3` for `lambda=-1,1,2`, put the top correction in
  three variables, and add lower complementary bridges.  Every ideal over
  `F_11` and `F_13` is a unit ideal.
- `hc4_fitting_denominator_extraction.json` records the canonical
  degree-three cube-torsion presentation, its sign-character block sizes,
  complete parameter-plane scans for four coefficient-monomial orbits over
  `F_11,F_13,F_17,F_19`, and exact rational specializations at the radial
  and reconstructed mixed point `(-5/3,-1/6)`.  At the mixed point the
  specialized coefficient quotient has dimension 60, exactly the three
  `x_i^2*x_j^2` cubes survive, and all fourth powers vanish.  This is an
  exact nilpotence-jump fiber calculation; the full integral zeroth Fitting
  ideal and associated primes remain open after 900-second timeouts.
- `hc4_canonical_signed_quadratic_cubic_words.json` records the complete
  fixed-order signed census of \(T_{H_2}\circ T_{H_1}\), with quadratic
  \(H_1\), cubic \(H_2\), and all flow and mixed-line signs in
  `{-1,1}`.  It classifies the 648 noncommuting words by their factored
  Poisson bracket before expansion.  Exact support leaves 216 words with
  two coordinate affine pivots; every reduced Hessian pencil is generically
  rank four, every parent Hessian is nonconstant, and all 34992 descended
  repair determinants are nonconstant.  It is generated by
  `scripts/search_hc4_mixed_quadratic_words.py`.  Its whole-file SHA-256 is
  `0baf8823d8ed090356d756b0ba9689b91402ce5fc6c6ea9044de40077d2b1f3d`.
- `hc4_symbolic_quadratic_cubic_words.json` records the coefficient-uniform
  parent-Hessian obstruction for all 54 noncommuting shared-dual
  support/sign patterns \(T_{bL_2^3}\circ T_{aL_1^2}\).  On the open locus
  `a*b != 0`, exact determinant differences give 14 localized monomial
  certificates and 40 saturated unit standard bases over `Q`; there is no
  parent-constant specialization.  It is generated by
  `scripts/verify_hc4_symbolic_quadratic_cubic_words.py` using Singular.
  Its whole-file SHA-256 is
  `c7ebbc6345d3037a9fdbecb080ca152ed7390cd6db7fcc5dd36dd7b8ffde082d`.
- `hc4_symbolic_cubic_quadratic_words.json` records the reverse-order
  coefficient-uniform family on the same 54 noncommuting shared-dual
  support/sign patterns.  Every pattern has an exact parent-preserving
  coefficient line: `a=+/-1/2` in the 48 one-sided cases and
  `a=+/-1/4` in the six reciprocal cases, with arbitrary `b != 0`.
  All 102 coordinate two-pivot reductions have generic rank four for every
  such `b`.  It is generated by
  `scripts/verify_hc4_symbolic_quadratic_cubic_words.py` using Singular.
  Its whole-file SHA-256 is
  `1c33fdbcb2296efe04df6e6e86d79bd407793a1adc5142f89e11912b419a36d9`.
- `hc4_mixed_quadratic_cubic_commutators.json` records the 162 distinct
  noncommuting unit quadratic--cubic group commutators.  Exact modular
  Hessian-chain-rule witnesses show that all 162 transformed parent
  determinants are nonconstant, across both one-sided bracket-incidence
  types and the reciprocal type.  It is generated by
  `scripts/search_hc4_mixed_commutator_words.py`.  Its whole-file SHA-256 is
  `95cd7757483cc71e97c8ed8925a0bce9e2d351794b2366a8ee83dab41f1ab359`.
- `fixed_quintic_hasse_curve_search.json` records the exact exclusion of
  every affine-linear base curve through the `Q(sqrt(-31))` Hasse point on
  the fixed-discriminant common-resolvent double cover, the exact exclusion
  of every degree-at-most-two curve on all three coordinate-fixed slices,
  and the bounded six-coefficient general quadratic search.  It is generated
  and checked by `scripts/search_fixed_quintic_hasse_rational_curves.py`;
  the unrestricted quadratic conclusion remains bounded.
- `cubic_homogeneous_counterexample.json` is generated by
  `scripts/cubic_homogeneous_reduction.py`.
- `cubic_linear_counterexample.json` is generated by
  `scripts/cubic_linear_reduction.py`.
- `long_bcw_79_counterexample.json` is the provenance-faithful conservative
  route generated by `scripts/verify_long_bcw_79_route.py`.
- `shared_bcw_33_counterexample.json` is the repository shared-factor
  optimization generated by `scripts/verify_shared_bcw_33_route.py`.
- `rank_compressed_bcw_24_counterexample.json` is the rank-compressed
  homogenization of that 16-variable map, generated by
  `scripts/verify_rank_compressed_bcw_24_route.py` and independently replayed
  by `scripts/audit_rank_compressed_bcw_24_independent.py`.
- `constant_kernel_bcw_22_counterexample.json` is the quotient of the
  24-variable map by the two-dimensional constant kernel of its homogeneous
  Jacobian, generated by `scripts/verify_constant_kernel_bcw_22_route.py`;
  `scripts/audit_bcw_22_linear_quotients.py` rules out any further
  collision-preserving linear quotient of this map.
- `index_reduced_bcw_22_counterexample.json` is the polynomial-circuit BCW
  witness of exact Jacobian nilpotency index 18, generated by
  `scripts/verify_index_reduced_bcw_22_route.py` and checked by the
  dependency-free `scripts/audit_index_reduced_bcw_22_independent.py`.
- `rank_reduced_bcw_24_counterexample.json` is the two-atom circuit witness
  of exact generic Jacobian rank 17 and exact index 18, generated by
  `scripts/verify_rank_reduced_bcw_24_route.py` and checked by
  `scripts/audit_rank_reduced_bcw_24_independent.py`.
- `hessian_rank_reduced_bcw_22_counterexample.json` is the two-atom
  22-variable cubic source whose 44-variable cotangent quartic has exact
  generic Hessian rank 37.  It is generated by
  `scripts/verify_hessian_rank_reduced_bcw_22_route.py` and checked by
  `scripts/audit_hessian_rank_reduced_bcw_22_independent.py`.
- `hessian_rank_35_identity_slice_counterexample.json` specializes the
  identity coordinate of that source to one.  It records a 21-variable
  nonhomogeneous nilpotent-Jacobian collision of exact rank 17 and its
  42-variable degree-2/3/4 HN potential of exact generic Hessian rank 35.
  It is generated by
  `scripts/verify_hessian_rank_35_identity_slice.py` and independently
  replayed by
  `scripts/audit_hessian_rank_35_identity_slice_independent.py`.  It does
  not change the homogeneous quartic rank or dimension endpoints.
- `hessian_rank_34_double_identity_slice_counterexample.json` uses the exact
  output relation `K_9-3*K_1-K_6=0` in the rank-35 slice.  Restricting its
  zero identity hyperplane gives a 20-variable nilpotent-Jacobian collision
  of exact rank 17 and a 40-variable degree-2/3/4 HN potential of exact
  generic Hessian rank 34 and kernel excess zero.  It is generated by
  `scripts/verify_hessian_rank_34_double_identity_slice.py` and independently
  replayed by
  `scripts/audit_hessian_rank_34_double_identity_slice_independent.py`.
- `identity_slice_hessian_rank_search.json` rescored all 140 terminals of the
  frozen width-64 circuit census after identity slicing.  The unique best
  modular profile is `(Hessian rank, Jacobian rank, excess, dimension) =
  (35,17,1,21)`.  `identity_slice_local_perturbation_search.json` continues
  64 neutral low-degree perturbations around that terminal and scores 414
  terminals, again finding no rank below 35.  Both are bounded diagnostics,
  not lower bounds; they are generated by the correspondingly named
  `scripts/search_identity_slice_*.py` programs.
- `restricted_minima_frontier.json` records the exact frozen-search intervals for
  the cubic rank/index and quartic HN rank/dimension minima, including the
  then-current pinned-external-certificate bounds `n_cub<=20` and
  `n_HN,4<=40` and the separate internal replay provenance.  The F12
  coordinate-pair artifact supersedes those two ambient endpoints; the
  frozen search ranks and indices are unchanged.  This artifact is generated by
  `scripts/verify_restricted_minima_frontier.py`.
- `cotangent_kernel_excess_frontier.json` verifies the cotangent block-rank
  decomposition at twelve deterministic points and records the excess
  profiles `(2,2,4,1)` of the four certified cubic sources, the exact
  excess-one value of the rank-37 witness, and the excess histogram of the
  140-terminal combined search.  It is generated by
  `scripts/analyze_cotangent_kernel_excess.py`.
- `index_three_inverse_degree_model.json` records an exact triangular
  cubic-homogeneous map with Jacobian nilpotency index three and two-sided
  inverse of degree nine.  It is generated by
  `scripts/verify_index_three_inverse_model.py`; it is a lower calibration,
  not a sharp full-class bound.
- `index_three_degree_bound_counterexample.json` verifies van den Essen's
  dimension-five cubic-homogeneous automorphism with generic Jacobian rank
  three, weak nilpotency index three, and inverse degree thirteen.  Its
  nonzero degree-eleven term realizes `Omega_11` on the full coefficient
  variety `(JH)^3=0`.  It is generated by
  `scripts/verify_index_three_degree_bound_counterexample.py`.
- `index_three_inverse_tree_obstruction.json` is the exact rational
  rooted-tree expansion through inverse degree eleven.  Euler and
  differential/context consequences of `(JH)^3=0` leave a three-dimensional
  quotient and one explicit normal form for `K_11`; van den Essen's exact
  tensor evaluates this form nontrivially, proving it cannot be killed by the
  full coefficient ideal.  It is generated by
  `scripts/derive_index_three_tree_obstruction.py`.
- `index_three_rank_normal_form_exclusion.json` proves that index three
  forces `lambda=0` in a binary-cubic four-variable normal-form family.
  The resulting locus has Jacobian rank
  at most two and an explicit inverse of degree at most five; it is generated
  by `scripts/verify_index_three_rank_normal_form.py`.
- `two_real_gmc_frontier.json` records the quadratic, two-weight, and
  affine-circular-source positive theorems for the unresolved GMC(2)
  frontier, the exact 27-support/72-chart cubic three-weight exclusion, the
  reduction of 33 cubic four-weight supports to four supports and 24 charts,
  and three exact bounded-support exclusions; it is generated by
  `scripts/verify_two_real_gmc_frontier.py`.
- `two_real_gmc_symmetric_chart.json` records three good-prime
  quotient-algebra certificates covering all four charts of the symmetric
  four-weight support.  Each order-eight quotient has dimension 84 and each
  tenth-moment multiplication matrix has modular rank 84 at a certified good
  prime; swapping the circular coordinates supplies the fourth chart.  It is
  generated by `scripts/verify_two_real_gmc_symmetric_chart.py` and lowers
  the live four-weight frontier from 24 charts on four supports to 20 charts
  on three supports.
- `two_real_gmc_remaining_four_weight.json` records seven exact rational
  unit-ideal certificates covering the last three supports and 20 charts.
  Centering identifies the two zero-weight coefficient choices, and
  circular-coordinate reflection supplies the reflected charts.  Moments
  through order six suffice in every representative.  It is generated by
  `scripts/verify_two_real_gmc_remaining_four_weight.py`; combined with the
  symmetric-support artifact, it closes all 121 four-weight cubic charts.
- `cubic_gaussian_null_cone_closure.json` records 31 exact
  characteristic-zero reduced bases `[1]`: a finite audit of all eleven
  two-weight presentations, the last three five-weight chart orbits, all
  fourteen six-weight presentations, and the three reflection orbits on the
  full seven-weight support.  Moments through order eight close the five- and
  six-weight strata; moments nine and ten close the full support.  It is
  generated by `scripts/verify_cubic_gaussian_null_cone_closure.py` and proves
  the radical moment-ten null-cone containment for cubics.
- `two_real_gmc_three_weight_low_degree.json` records the 31 exact rational
  unit-ideal certificates for the three-level support `{-1,0,1}` through
  total degree six: 6 charts in degree four, 10 in degree five, and 15 in
  degree six.  It also records the Bessel--factorial moment identity
  independently checked against direct Wick contraction.  It is generated by
  `scripts/verify_two_real_gmc_three_weight_low_degree.py`.
- `minimal_counterexample_scoreboard.json` records the separated ambient
  dimension, rank, nilpotency-index, and degree intervals.  Its dimension-20
  cubic-homogeneous and dimension-40 homogeneous-quartic entries depend on
  the pinned external MacFarlane determinant certificate; it is generated by
  `scripts/verify_minimal_counterexample_scoreboard.py`.
- `dvorsky_gvc5_counterexample.json` records the exact five-variable
  unrestricted constant-coefficient GVC witness and its five-pair SIC lift.
  It is generated by the dependency-free
  `scripts/audit_dvorsky_gvc5_counterexample.py`.
- `binary_degree_five_gvc_face_search.json` records the exhaustive modular
  screen of the two residual cubic-leading weighted faces at
  \(p=101,103,107\), totaling 6,696,142 raw triples, and every projective
  top form at every squarefree quartic cross-ratio over seven further
  primes, totaling 2,082,612 tuples.  Only the two predicted one-sided
  lines and the four root fifth powers survive.  Regenerate it with
  `scripts/search_binary_degree_five_gvc_faces_mod_p.py`.  This is a bounded
  experiment; the characteristic-zero proof is in
  `extended-geometry/BINARY_DEGREE_FIVE_GVC_FRONTIER.md`.
- `binary_degree_five_gvc_second_moments.json` records the normalized first
  equation and complete second-moment jet layers for all eight
  cubic-leading quintic normal forms.  Regenerate it with
  `scripts/explore_binary_degree_five_gvc_frontier.py`; add
  `--triangular-components` for a heuristic branch dump without replacing
  the default artifact schema.
- `binary_repeated_quartic_gvc_jet_search.json` records bounded modular
  samples of the first migrating repeated-root quartic jets
  \((\Lambda_4+\Lambda_5,P_5+P_4)\) and the conditioned defect-two
  \((\Lambda_4+\Lambda_5+\Lambda_6,P_5+P_4+P_3)\) slice.  The primes are
  larger than every input degree in the recorded moment windows.  Regenerate
  it with `scripts/search_binary_repeated_quartic_gvc_jets_mod_p.py`.  It is
  an experiment, not an exhaustive search, proof, or counterexample.  Its
  whole-file SHA-256 is
  `299297fcf936c9041cb2da7ab4f7d124271504e09af02066f92d0cfc4e930f14`.
- `three_pair_image_mathieu_counterexample.json` records the eight-term
  three-pair SIC witness with multiplier `g=z`, its bidegrees, and the exact
  all-order contraction identities.  It is generated by the dependency-free
  `scripts/verify_three_pair_image_mathieu_counterexample.py`.
- `factorial_moment_witnesses.json` records the Dvorsky--Long
  torus-diagonal factorial translation, its order-two nonmultiplicativity
  certificate, exact cyclotomic linear searches, and witness-derived
  quartics with sharp finite zero prefixes.  It is generated by the
  dependency-free `scripts/verify_factorial_moment_witnesses.py`.
- `sparse_factorial_moment_frontier.json` records the exact exhaustion of
  3,276 three-term binary supports and 4,950 involution-paired four-term
  supports, together with projective Gröbner certificates for the sharp
  homogeneous binary cutoffs in degrees one through four.  It is generated
  by `scripts/verify_sparse_factorial_moment_frontier.py`.
- `restricted_bcw_circuit_search.json` is the bounded modular discovery
  record for the circuit-level search.  Its profiles remain diagnostic; the
  separately frozen index-reduced artifact is the certified result.
- `restricted_bcw_circuit_search_v2.json` and
  `restricted_bcw_circuit_search_v2_w64.json` are the five-atom Pareto
  searches.  The width-64 record retains the modular Hessian-rank-37
  candidate subsequently promoted by the exact generator and audit above.
- `restricted_bcw_circuit_search_xxs_w32_d28.json`,
  `restricted_bcw_circuit_search_yvyb_structural_w22.json`, and
  `restricted_bcw_circuit_search_ayb_yvyb_w36.json` retain the bounded
  negative tests of the alternative linear--cubic, balanced, and
  gate-sharing circuit factorizations.  None beats a certified endpoint.
- `restricted_bcw_circuit_search_xvvz_v2h_structural_w25.json` and
  `restricted_bcw_circuit_search_xvvz_v2h_mixed_w25.json` retain the two
  bounded beams that share the exposed `v*z` gate across target components.
  Their best diagnostic vector is `(18,18,41,23)`.
- `rank37_gate_perturbation_search.json` continues eight coefficients of the
  no-new-gate `lambda*(xy)*(xy^2)` shear from the frozen rank-37 circuit.
  Every family leader reproduces `(18,18,37,22,35)`.
- `restricted_bcw_circuit_search_aspert_m12.json` restarts the exceptional
  `-12*(xy)*s` shear after it cancels a frozen cleanup step.  Its 33
  terminals have best vector `(18,18,41,23)`.
- `restricted_bcw_circuit_search_xxs_rank_hybrid_w24.json`,
  `restricted_bcw_circuit_search_xxs_v2r_w16.json`, and
  `restricted_bcw_circuit_search_xxs_y2vb_w16.json` test whether the
  zero-excess `xxs` gate composes with the known rank-reducing atoms.  Their
  best vectors are respectively `(22,19,45,28)`, `(21,19,44,27)`, and
  `(19,18,43,26)`; none improves an endpoint.
- `restricted_bcw_circuit_search_all_w64.json` is the combined nine-atom,
  width-64, 24-depth search over 140 terminal maps.  It finds no profile
  below Hessian rank 37 and rediscovers the certified `qb+x2s` route.
- `weighted_seed_scan.json` is an untracked exploratory output generated by
  archived `archive/tooling/scan_weighted_seeds.py`.
- `branch_scale_fan_degree6.json` records the exact degree-six `(2,2,2)`
  critical-value valuation vectors, the six unimodular radial braid-fan
  cones, and the triple-resonance refinement generated and checked by
  `scripts/verify_branch_scale_fan.py`.
- `branch_target_graph_degree6.json` records the six-line target arrangement,
  the four-point Kapranov blowup of `P^2`, its ten-curve Petersen boundary,
  and the four reduced triple-Maxwell pullback branches checked by
  `scripts/verify_degree_six_branch_target_graph.py`.
- `admissible_equal_scale_degree6.json` records the degree-six central
  component, its three quadratic tails, the two admissible target bubbles,
  and the four-branch normalization of the common index-two node deformation
  checked by `scripts/verify_degree_six_admissible_equal_scale.py`.
- `admissible_radial_atlas_degree6.json` records all thirteen ordered radial
  scale types for three quadratic clusters, their degree-six bubble
  decompositions, node-index partitions, and Kummer normalization counts
  checked by `scripts/verify_degree_six_admissible_radial_atlas.py`.
- `admissible_maxwell_atlas_degree6.json` records the three pairwise and one
  triple Maxwell source-node charts, their Kummer branches, and all
  intersections with the radial boundary checked by
  `scripts/verify_degree_six_admissible_maxwell_atlas.py`.
- `central_hurwitz_selection_degree6.json` records the two ambient central
  Hurwitz classes and the exact square-cubic invariant proving that labelled
  source roots select the polynomial class, checked by
  `scripts/verify_degree_six_central_hurwitz_selection.py`.
- `labelled_node_saturation.json` records the general phase quotient,
  label-preserving inertia test, permutation equivariance, anchored and
  unanchored cyclic-tail families, and the corrected subgroup-quotient
  degree checked by `scripts/verify_labelled_node_saturation.py`.
- `branch_wonderful_pullback.json` records the universal
  `Mbar_0,n` boundary building set, maximal nested-set counts, permutation
  action, normalized pullback rule, and the common degree-five/degree-six
  recovery checked by `scripts/verify_branch_wonderful_pullback.py`.
- `source_vertex_rigidity.json` records the exhaustive genus-zero
  two-fiber reconstruction, divisor-permutation invariance, third-flag
  normalization, and the central, cluster-tail, and cyclic-tail degree-six
  instances checked by `scripts/verify_source_vertex_rigidity.py`.
- `general_radial_source_atlas.json` records the all-multiplicity
  connector/local-tail/identity rule, ordered scale-type counts, target
  bubble and node checks, lcm saturation, label equivariance, the
  full-chain radial inertia product on all 48,580 bounded types, and the
  unequal `(2,3,4)` order-24 example checked by
  `scripts/verify_general_radial_source_atlas.py`.
- `polynomial_monodromy_forests.json` records every reduced polynomial
  transposition factorization through degree six, all collision subforests,
  nested refinements, componentwise Riemann--Hurwitz checks, and the common
  Maxwell/caustic node rule checked by
  `scripts/verify_polynomial_monodromy_forests.py`.
- `monodromy_inertia_characters.json` records polynomial-tree and cyclic
  connector deck centralizers, all anchored and unanchored simple-resonance
  node inertias through degree six, and the common radial/Maxwell/caustic
  character rule checked by `scripts/verify_monodromy_inertia_characters.py`.
- `recursive_resonance_atlas.json` records the framed nested residue screens,
  their explicit affine gauge action, exact one- and two-step contraction
  maps, normalized initial-form flag equations including nonfactorized
  smoothing families, arbitrary tame
  full-centralizer extraction with concrete inertia subgroups, the
  unequal-multiplicity radial product formula on 76 bounded charts, every
  interval-nested degree-six
  matching-tree chart, all thirteen equal-multiplicity radial types, a
  surviving unrigidified ghost example, and the order-four pair--triple
  Maxwell inertia checked by
  `scripts/verify_recursive_resonance_atlas.py`.
- `stack_inertia_degree6.json` separates normalization branches from actual
  label-preserving inertia on all thirteen radial types and the
  pairwise/triple Maxwell bubbles, checked by
  `scripts/verify_degree_six_stack_inertia.py`.
- `stacky_fan_descent_degree6.json` records the four Maxwell root divisors,
  their pair--triple face intersections, four radial `S_3` quotient orbit
  types, the codimension-two inertia audit, and the separate local
  three-pair wreath-product descent checked by
  `scripts/verify_degree_six_stacky_fan_descent.py`.
- `universal_cubic_quartic_kernel_saturation_frontier.json` records 28
  exact full-support quartic-kernel lines across the seven squarefree cubic
  symbols and one exact smooth coordinate ten-space.  It also serializes
  the ten components of the universal 24-parameter tensor and their
  deterministic hash.  Over the full parameter ring it certifies that all
  parameter-dependent cotangent entries begin in collision degree three
  and splits six unit pivots, reducing the raw 12-by-31 presentation to a
  cokernel-equivalent 6-by-25 matrix.  Every recorded family has saturated
  relative cotangent presentation and flat length-six `Ext^2`; the artifact
  itself stops before universal saturation.
- `universal_cubic_cotangent_saturation.json` records the exact
  three-variable formal-gauge matrices `C`, `G`, the cubic-modulus tensor
  `eta`, and the explicit lift `L` with
  `G*L=[x*eta,y*eta,z*eta]`.  Exact module reduction proves
  `ker(C)=im(G)+A*eta` and `ker(C)/im(G)=Q` in collision degree three.
  A direct dual-number expansion derives all nine columns of `G` from the
  determinant-twisted finite gauge action.  The artifact also stores an
  explicit 9-by-24 linear-polynomial matrix `Q` satisfying
  `G*Q=[psi_1,...,psi_24]`, together with all exact matrices rather than
  only their hash.
  It is checked by
  `scripts/verify_universal_cubic_cotangent_saturation.py` and proves the
  smooth-symbol 24-parameter cotangent saturation theorem without
  computing the full universal saturation.
- `cubic_formal_gauge_cokernel_atlas.json` records the exact compatibility
  matrix and determinant-twisted gauge matrices for all ten ternary-cubic
  symbols.  Singular computes the Hilbert numerator of
  `ker(C)/im(G_h)` in each row.  The artifact proves that smooth is the
  unique symbol whose gauge cokernel vanishes above collision degree
  three, and records exact quartic nongauge dimensions
  `0; 2,4,4,6,6,8; 11,16,24`.  It also records the exact principal
  annihilators of all ten cokernels; the non-squarefree rows are faithful
  of generic ranks `1,2,4`.  Nongauge tensor moduli are not asserted to be
  cotangent torsion.  It is checked by
  `scripts/verify_cubic_formal_gauge_cokernel_atlas.py`.
- `nodal_cubic_formal_slice.json` records the cyclic identification
  `ker(C)/im(G_nodal)=Q[y,z](-3)`, an explicit lift of `x*eta`, and the
  complete 9-by-24 quartic gauge-lift matrix.  The first two fixed quartic
  directions give a complementary two-plane; the other 22 are gauge.
  The full-support sum/alternating-sum plane has slice matrix
  `[[1,-1],[1,1]]`, of determinant two.  Exact Singular replays show that
  both transverse planes have saturated cotangent presentation and the
  central length-six Ext block.  For the deterministic row-reduced lift,
  the artifact also stores the complete three-component degree-five
  quadratic curvature, which vanishes on the coordinate slice and is
  nonzero on the dense slice.  It then stores the five-dimensional
  quartic-lift kernel, its rank-four action on the six slice--gauge
  curvature coordinates, the two quotient cross forms, and the three
  lift-independent pure-gauge quadrics.  Their reduced zero scheme is two
  rational planes; the unreduced ideal has one embedded quadratic socle
  generator at the origin.  On both reduced planes, the artifact records
  an exact quadratic correction of the degree-five term and the resulting
  degree-six classes `27/8*(q*y+p*z)^3*eta` and
  `27/8*(q*y-p*z)^3*eta`.  It also checks that the 15-dimensional
  correction ambiguity acts trivially on the degree-six quotient.  It is
  checked by
  `scripts/verify_nodal_cubic_formal_slice.py` and does not claim an
  all-order nodal normal form, quartic-lift independence at degree six, or
  continuation of the embedded socle.
- `universal_cubic_filtered_syzygy_frontier.json` records the smooth
  central minimal cotangent resolution of ranks `7 -> 13 -> 6`, the exact
  regular boundary form `x+y+z`, the collision-degree comparison of all
  150 entries in the reduced universal matrix, and the twelve nonzero
  remainders obtained by applying the unchanged central syzygies
  universally.  It is checked by
  `scripts/verify_universal_cubic_filtered_syzygy_frontier.py` and
  explicitly leaves corrected Rees strictness and universal saturation
  open.
- `ritt_cellular_prototype_completion.json` records the complete labelled
  degree-42 and degree-30 factor/move/braid-cell diagrams, scalar incidence
  matrices, totalized coefficient complexes, every certified filtration
  cohomology row, uniform prototype `H2=0`, and the exact order-three
  sector--spectator obstruction consumed from HRCELL2--HRCELL5.  It is
  checked by `scripts/verify_ritt_cellular_prototype_completion.py`.
- `hessian_ritt_cotangent_descent.json` records exact face-bar-to-cellular
  subdivision comparisons, the genuine permutohedron top-cell boundary,
  all six degree-42 Hessian vertex and move tangent ranks, and the common
  Dickson-plus-excess tangent space.  It also records the uniform
  `(5,6,6,7)` first conormal flag for the three labelled degree-42 sectors.
  It is checked by `scripts/verify_hessian_ritt_cotangent_descent.py`.
- `degree42_ritt_rotated_conormal_jet_237.json` and
  `degree42_ritt_rotated_conormal_jet_327.json` record the exact fourth
  maximal-adic jets of the cut-`14` and cut-`21` sectors.  Their spectator
  dimensions agree as `(1,3,6)`, while their sector dimensions are
  `(1,4,9)` and `(1,4,10)`.  They are generated by the two specialized
  commands for
  `scripts/explore_degree42_ritt_rotated_conormal_flags.py`; neither
  artifact claims a completed-local-ring result.
- `degree42_ritt_rotated_source_ideals_237.json.gz` is the compressed
  one-time cache of the 76 thick, 76 thin, and 102 boundary residual
  equations on the normalized `237` chart.  It is consumed by
  `scripts/verify_degree42_ritt_cut14_postnikov_overlap.py`.
- `degree42_ritt_cut14_postnikov_overlap.json` records the exact Nakayama
  containment proving completed thin-path/boundary equality, the cutoff-4
  sector source, and the cutoff-5 zero quadratic overlap for the composite
  omission `14`.  The equal overlap colengths are `34`; hence the completed
  first-Postnikov conormal sequence is short exact on this half-braid.
- `degree42_ritt_rotated_source_ideals_327.json.gz` is the corresponding
  compressed cache of the 76 thick, 76 thin, and 96 boundary residual
  equations on the normalized `327` chart.
- `degree42_ritt_cut21_postnikov_overlap.json` records the exact Nakayama
  containment and cutoff-4/cutoff-5 calculations for composite omission
  `21`.  Its quadratic numerator has 61 generators and its equal overlap
  colengths are again `34`, proving completed first-Postnikov exactness.
- `degree42_ritt_rotated_tensor_matrices_237_b4_n1.txt.gz` caches the exact
  Singular bases and action reductions for the cut-`14` base-order-four
  tensor presentation.
- `degree42_ritt_cut14_tensor_split_q4.json` records its exact
  `9 -> 13 -> 4` module extension, an explicit section, and equal
  coboundary/augmented ranks `32=32`.  It proves splitting through order
  four.
- `degree42_ritt_rotated_tensor_matrices_237_b7_n1.txt.gz` and
  `degree42_ritt_rotated_tensor_matrices_327_b7_n1.txt.gz` cache the exact
  order-seven tensor presentations and action reductions for both rotated
  sectors.
- `degree42_ritt_inverse_limit_sections_q5_q7.json` constructs orders five
  and six as quotients of those order-seven modules.  It records compatible
  sections with dimensions `12 -> 17 -> 5`, `15 -> 21 -> 6`, and
  `18 -> 25 -> 7` in both sectors, as well as the two-dimensional
  cokernels of the section-difference restriction maps.
- `degree42_ritt_completed_presentations_237.txt.gz` and
  `degree42_ritt_completed_presentations_327.txt.gz` cache the exact
  two-variable module presentations after killing the seven graph-normal
  directions.
- `degree42_ritt_completed_splits.json` records the polynomial completed
  sections
  `e4+(-3*(1+tau)^2+2*zeta)*e6` and
  `e4+(-4*(1+tau)^3+8*(1+tau)*zeta)*e7`.  Exact module reduction proves
  that both rotated first-conormal extensions split over
  `Q[[tau,zeta]]`, so their completed extension and inverse-limit torsor
  obstruction classes vanish.  Full braid restriction coherence is not
  claimed.
- `boundary_obstruction_theory.json` records exact module saturation and
  boundary-torsion controls, a finite-jet tower with surjective transitions
  but unbounded annihilation exponents, node and cusp conductor kernels and
  their rank-three free tensor blocks, and strict versus non-strict filtered
  lifting controls.  It is checked by
  `scripts/verify_boundary_obstruction_theory.py`.
- `conductor_first_one_chart_obstruction.json` records explicit nodal and
  cuspidal finite marked-root algebras with descended discriminants, the
  complete conductor divisibility equations, exact low-order unit-ideal
  regressions, and the stable unit-group obstruction for the separated
  one-conductor one-chart ansatz.  It is checked by
  `scripts/verify_conductor_first_one_chart_obstruction.py`.
- `conductor_first_foundational_cusp_keller.json` begins with the cusp
  conductor `Q[u^2,u^3] subset Q[u]`, reconstructs the foundational cubic
  tangent incidence and its marked-root discriminant, and verifies the
  reconstruction pole, exact denominator cancellation, distributed
  determinant ledger, Jacobian `-1/2`, and three-point collision.  It is
  checked by `scripts/verify_conductor_first_foundational_cusp_keller.py`.
- `support_saturation_universal_cubic_symbols.json` records exact module
  saturations, associated primes, and regular boundary elements for all ten
  homogeneous ternary-cubic symbol strata.  It does not claim arbitrary
  nonhomogeneous universal cotangent saturation.
- `support_saturation_cubic_annihilator_frontier.json` imports the proved
  formal-gauge cokernel atlas into the support-saturation workflow.  It
  closes the smooth-symbol quartic saturation search, queues the six
  singular squarefree quotients by annihilator type
  `(x),(x^2),(yz),(y^3),(xyz),(x^3)`, and records the
  generically-étale/Keller gate for double-line, triple-line, and zero
  symbols.  It is a routing certificate, not a new saturation theorem.
- `support_saturation_degree42_ritt_fiber_mod32003.json` records an exact
  characteristic-zero nonzero support class on the full reduced degree-42
  core at `(e1,e2,t)=(1,2,3)`, together with the full order-six and
  order-seven support saturations over `GF(32003)`.  Their local cohomology
  has boundary exponents one and two, and the transition is surjective.  The
  stored order-six representative is not itself torsion at order seven but
  has a different torsion lift.  The artifact does not contain the full
  untruncated saturation or associated-prime list and does not lift the
  finite-jet result to characteristic zero.
- `degree42_c6_macaulay_certificate.json` closes the fixed order-six
  characteristic-zero gap on the same specialized fiber.  It stores the
  normalized 16-generator sparse core, rational Macaulay syzygies for
  `w0*c6` and `w2*c6`, a nine-term finite-support functional taking `c6` to
  one and annihilating the whole ideal, and transcripts for eight 31-bit
  prime block-Wiedemann runs for each linear system.  The certificate proves
  that `c6` is a nonzero `(w0,w2)`-torsion class modulo `(u,v)^6`; it does not
  compute the full rational saturation or the order-seven module.  Regenerate
  it with `scripts/compile_degree42_c6_macaulay.py` and replay it independently
  with `scripts/verify_degree42_c6_macaulay.py`.
- `support_saturation_plane_jc_boundary_layer.json` records the local
  cohomology and associated primes of the normalized cyclic `d3`
  multiplication-kernel layer in the Poisson-square coefficient scheme.
  It is not the still-undefined Case-1 conductor/residue matching module.
  The three Singular support-saturation artifacts are generated by
  `scripts/compile_support_saturation_cases.py` and checked structurally by
  `scripts/verify_support_saturation_compiler.py`; the degree-42 Macaulay
  artifact has the separate generator and independent checker named above.
- `hc4_projective_polar_atlas.json` records all positive, degree-bounded,
  log-concave projective-degree lists for degree-two through degree-four
  maps `P4 --> P4` with top degree two or three, together with the
  corresponding Segre-degree signatures.  It also records the
  zero-dimensional lengths, smooth-lci-curve numerical filter,
  graph/full-polar calibrations, aggregate cotangent and Meng--Yang
  controls, and the exact consequences of Wang's theorem and `HC4CQ1`.
  For the 626 quartic-gradient rows it records the leading-quintic Hessian
  rank coverage matrix: rank zero and the aligned rank-three branch are
  closed, while ranks one and two require kernel synchronization and the
  nonaligned rank-three branch is an explicit ternary quintic--cubic Schur
  divisibility problem.  The degree-nine-versus-degree-eight adjugate
  argument closes that branch whenever the ternary-quintic Hessian
  determinant is squarefree; an exact squarefree witness proves this is
  the generic case.  Only its nonsquarefree discriminant locus remains.
  The Hessian matrix alone does not assign individual signatures, but the
  attached top-gradient/Rees sieve partitions the rows in the other
  direction: generic smooth essential Hessian ranks one, two, and three
  feed only codimensions two, three, and four, while singular top quintics
  feed the stated lower-codimension columns.  The smooth rank-three vertex
  colength theorem then excludes the two codimension-four
  affine-degree-two/three rows, leaving 318 and 306 numerical signatures.
  The codimension-three gradient sieve further closes the squarefree
  binary-Hessian open when \(h_4|_K\ne0\), forces \(\sigma_3=16\) on its
  nonsquarefree remainder, and records the Schur-cubic incidence at
  ordinary rank-three singularities; it does not delete a row
  unconditionally.
  Of the
  codimension-three rows, 18 for each affine degree pass
  the additional smooth-integral-curve genus and Castelnuovo test.  The
  rows are a numerical pre-Keller atlas, not existence results or a
  classification by Hilbert polynomial.
  The attached `HC4PPG9` packet records a further conditional split of the
  codimension-two repeated-binary-root rows: on the active-unit stratum a
  partition with \(q\) distinct roots has \(\sigma_2=5-q\), and the
  generic double-root packet retains 51 and 50 rows.
  Regenerate it with
  `.venv/bin/python scripts/verify_hc4_projective_polar_atlas.py`; run the
  independent calibration with
  `M2 --script scripts/verify_projective_polar_calibrations.m2`.  The
  canonical source is `HC4_PROJECTIVE_POLAR_GEOMETRY.md`.
- `projective_gradient_segre_registry.json` records the dimension-free
  projective-degree/Segre transform convention and attaches it to the
  repository's explicit triangular, cotangent, stabilization, Meng--Yang,
  Schur-descended, and restricted-minima HN families.  Full vectors,
  top-degree-only controls, and open computations are separate fields.
  Regenerate it with
  `.venv/bin/python scripts/verify_projective_gradient_segre_machinery.py`
  and independently certify the cotangent and stabilization multidegrees
  with
  `M2 --script scripts/verify_projective_gradient_segre_families.m2`.
  The canonical source is `PROJECTIVE_GRADIENT_SEGRE_MACHINERY.md`.
- `projective_gradient_normal_slices.json` records theorem `PGS2` for every
  regression tuple \(2\le n\le10\), \(2\le m\le7\), and \(1\le r<n\).
  It gives the smooth-essential Jacobian Hilbert function and length
  \(m^r\), the compactifier-truncated active length \(m^{r+1}\), the
  filtered missing-generator bound, and the exact unit-penultimate Segre
  law.  The range is a regression ledger, not a dimension bound.  Generate
  it with
  `.venv/bin/python scripts/verify_projective_gradient_normal_slices.py`;
  replay exact complete-intersection calibrations with
  `M2 --script scripts/verify_projective_gradient_normal_slices.m2`.  Its
  whole-file SHA-256 is
  `5853c8fa609879663b31f680591a5e612ab944b1637902de5dcd115c9400837b`.
  The canonical source is `PROJECTIVE_GRADIENT_SEGRE_MACHINERY.md`.
- `projective_gradient_singular_slices.json` records theorem `PGS3`: an
  essential singular component \(C\) of dimension \(s\) and degree \(d\)
  gives the infinity component
  \(\operatorname{Join}(\mathbf P^{n-r-1},C)\) of codimension \(r-s\),
  while its transverse multiplicity depends on the Jacobian length
  \(\mu_C\) and the finite \(K[[X_0]]\)-module profile
  \((\rho_C;a_{C,j})\).  It also records the exact truncated length
  \(m\rho_C+\sum_j\min(m,a_{C,j})\) and the binary-quintic profiles of
  lengths \(8,3,2\).  Generate it with
  `.venv/bin/python scripts/verify_projective_gradient_singular_slices.py`;
  independently replay the three local lengths with
  `M2 --script scripts/verify_projective_gradient_singular_slices.m2`.
  Its whole-file SHA-256 is
  `c6971874b5359e4aed11a8918328804f9ffdd6e67811f49c9ff79b2a8c5d7b72`.
  The canonical source is `PROJECTIVE_GRADIENT_SEGRE_MACHINERY.md`.
- `hc4_quintic_infinity_rees_strata.json` records the universal
  56-coefficient top quintic, its Euler/Hessian/curl/Koszul and
  midpoint-parity identities, the generic smooth essential
  Hessian-rank-one/two/three Rees models, their pure-top
  complete-intersection Segre calibrations, and the exact support-codimension
  intersection with the 626-row atlas.  It separates those pure-top vectors
  from the unknown lower-layer normal-cone multiplicities.  Regenerate it
  with
  `.venv/bin/python scripts/analyze_hc4_quintic_infinity_rees.py`;
  independently certify the linear-type and projective-degree claims with
  `M2 --script scripts/verify_hc4_quintic_infinity_rees_strata.m2`.
  The recorded independent replay uses Macaulay2 1.22 with `Cremona` and
  `ReesAlgebra` over \(\mathbb Q\).
  The canonical source is `HC4_PROJECTIVE_POLAR_GEOMETRY.md`.
- `hc4_rank3_vertex_colength.json` records the exact filtered-length
  obstruction on the smooth essential rank-three quintic vertex.  The
  active ternary Jacobian complete intersection has Hilbert function
  `(1,3,6,10,12,12,10,6,3,1)` and length \(64\); after truncation at
  \(\epsilon^4\), a nonaligned missing gradient component generates an
  ideal of length at least six.  This excludes the two codimension-four
  affine-degree-two/three atlas rows.  Generate it with
  `.venv/bin/python scripts/verify_hc4_rank3_vertex_colength.py`; replay the
  exact Fermat and deformed calibrations with
  `M2 --script scripts/verify_hc4_rank3_vertex_colength.m2`.  Its whole-file
  SHA-256 is
  `c610f57af67061d0b4eb9523cb018569a7e8220a51dbd2350b71eb7007bfe473`.
  The canonical source is `HC4_PROJECTIVE_POLAR_GEOMETRY.md`.
- `hc4_codim3_gradient_strata.json` records the exact conditional sieve on
  the two codimension-three top-gradient packets.  For smooth essential
  rank two, a nonzero \(h_4|_K\) synchronizes a constant kernel direction;
  the squarefree binary-Hessian branch reaches `HC4CD5`, while its
  nonsquarefree remainder has forced generic multiplicity
  \(\sigma_3=16\).  For essential rank three, the Schur cubic vanishes at
  every isolated singular point where the top Hessian has rank two.
  Generate it with
  `.venv/bin/python scripts/verify_hc4_codim3_gradient_strata.py`; replay
  the radical powers and transverse lengths with
  `M2 --script scripts/verify_hc4_codim3_gradient_strata.m2`.  Its
  whole-file SHA-256 is
  `8759875cf431d18f35321631984d9120c72a2335dcae31d107fa191ae539e5a3`.
  The canonical source is `HC4_PROJECTIVE_POLAR_GEOMETRY.md`.
- `hc4_binary_root_partition_segre.json` records theorem `HC4PPG9`, the
  first direct singular-top application of `PGS3`.  At a binary-quintic
  root of multiplicity \(e\), the transverse Jacobian length is \(e-1\).
  On the active-unit lower-layer stratum this is its exact contribution to
  \(\sigma_2\), so a root partition with \(q\) distinct roots forces
  \(\sigma_2=5-q\).  The generic \(2+1+1+1\) packet retains 51 and 50
  atlas rows.  Generate it with
  `.venv/bin/python scripts/verify_hc4_binary_root_partition_segre.py`;
  replay multiplicities \(2,3,4\) with
  `M2 --script scripts/verify_hc4_binary_root_partition_segre.m2`.  Its
  whole-file SHA-256 is
  `09f6a57c735b2751d0f890b8cd216822001bae875fd9a5156e2a27550f8e71ad`.
  The higher-torsion failure locus remains open.  The canonical source is
  `HC4_PROJECTIVE_POLAR_GEOMETRY.md`.
- `quartic_coefficient_kuranishi_mod32003.json` records the modular
  36-variable affine-normal quadratic slice and coordinate-axis jet screen at
  the integer-root quartic map.  The adjacent `.sing` and `.m2` files are
  reproducible primary-decomposition inputs, not completed decomposition
  certificates.  Generate them with
  `scripts/research_quartic_coefficient_kuranishi.py`; certify the full
  characteristic-zero quadratic rank independently with
  `scripts/verify_quartic_full_box_kuranishi.py`.
- `quartic_generic_component_mod32003.json` records the 22-variable normal
  quadratic slice at a rational generic point of the explicit
  27-dimensional quartic reduced family.  Its coordinate-axis cubic and
  greedy higher-jet screens are modular discovery data; the greedy
  order-four failures are not obstruction certificates.  The adjacent
  `.sing` and `.m2` files are decomposition inputs, and no completed radical
  or associated-prime calculation is claimed.  Generate them with
  `scripts/research_quartic_generic_component.py`.
- `quartic_generic_component_order3_mod32003.sing` contains the optional full
  canonical cubic homogeneous Kuranishi layer: 305 independent cubics
  together with the 22 quadrics.  It is a reproducible Singular input, not
  the completed local ideal or a standard-basis certificate.  Its capped
  12 GB standard-basis run ended without a basis.
- `filtered_source_tangent_profiles_mod32003.json` records the source-only
  first-order degree-filtration breakpoints for `F_4,F_5,F_6`, including the
  visible seed subspaces.  It is a good-prime filtered computation, not an
  optimal two-sided or characteristic-zero contact theorem.  Generate it with
  `scripts/research_filtered_source_tangent_profile.py`.

Run `make verify-normal-forms` for the original normal-form artifacts and
`make verify-external-consequences` for all BCW route artifacts.
