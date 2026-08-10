# Exact replay commands

The normalized sparse-support certificate is:

```bash
.venv/bin/python plane-jc/cas/verify_sparse_support_exclusions.py
```

It replays the arbitrary-degree singleton-versus-five classification's
nonempty Gröbner charts and their two-sided shear inverses.  For the
balanced `2+2` class it exhausts all 256 divergence/determinant presence
masks and the 20 residual collision partitions with exact integer
arithmetic, proving that every arbitrary-degree support has a singleton
Keller coefficient and hence an explicit Rabinowitsch unit identity.  The
former `14,653,584`-row census through degree twelve is retained as an
independent regression.  The same replay closes support five: `1+4` forces
the cubic shear `Q=b(x+a*y^m)^3`, while all `2+3` and transposed charts are
unit ideals after an exact `2048 -> 321 -> 98 -> unsat` collision sieve.
It then closes support six in arbitrary degree: `1+5` is the quartic shear
`Q=b(x+a*y^m)^4`, `2+4` and `4+2` are unit ideals, and the only `3+3`
collision support is the directional quadratic shear on
`{x^2,x*y,y^2}`.  A separate `5,290,000`-support census through degree six
finds that same unique `3+3` collision support.  The replay checks the
resulting affine-normalized support lower bound seven.
The pinned artifact is
`artifacts/generated-results/jc2_sparse_support_exclusions.json`.
Intentional regeneration requires `--refresh`.  The theorem and its
coordinate-dependent claim boundary are in
[`../CERTIFIED_SPARSE_SUPPORT_EXCLUSIONS.md`](../CERTIFIED_SPARSE_SUPPORT_EXCLUSIONS.md).

The attempted Newton/boundary bridge is replayed by:

```bash
.venv/bin/python plane-jc/cas/verify_affine_support_newton_bridge.py
.venv/bin/python plane-jc/cas/classify_f2_75_125_layers.py
.venv/bin/python plane-jc/cas/audit_f2_75_125_boundary_handoff.py
.venv/bin/python plane-jc/cas/test_f2_75_125_frontend.py
.venv/bin/python plane-jc/cas/generate_f2_modified_system.py --include-equations --output artifacts/generated-results/jc2_f2_modified_laurent_family.json
.venv/bin/python plane-jc/cas/verify_f2_modified_chart_bridge.py
.venv/bin/python plane-jc/cas/verify_f2_kummer_orbit_transfer.py
.venv/bin/python plane-jc/cas/verify_f2_terminal_residue_cover.py
.venv/bin/python plane-jc/cas/verify_f2_a6_simple_spectator_gluing.py
.venv/bin/python plane-jc/cas/verify_f2_75_125_global_attachment.py
.venv/bin/python plane-jc/cas/test_f2_75_125_global_attachment.py
.venv/bin/python scripts/verify_f2_carrier_log_node_profiles.py
.venv/bin/python plane-jc/cas/verify_f2_75_125_carrier_wronskian.py
.venv/bin/python plane-jc/cas/test_f2_75_125_carrier_wronskian.py
.venv/bin/python plane-jc/cas/verify_f2_75_125_carrier_specializations.py
.venv/bin/python plane-jc/cas/test_f2_75_125_carrier_specializations.py
.venv/bin/python plane-jc/cas/compile_f2_75_125_nonlinear_forcing.py
.venv/bin/python plane-jc/cas/test_f2_75_125_nonlinear_forcing.py
.venv/bin/python plane-jc/cas/test_sparse_circuit_modp.py
.venv/bin/python plane-jc/cas/probe_f2_75_125_nonlinear_modular.py
.venv/bin/python plane-jc/cas/compile_f2_75_125_tangent_obstruction.py
.venv/bin/python plane-jc/cas/compile_f2_75_125_formal_homotopy.py
.venv/bin/python plane-jc/cas/compile_f2_75_125_formal_homotopy.py --regular-gauge --artifact artifacts/generated-results/jc2_f2_75_125_formal_homotopy_regular_gauge.json
.venv/bin/python plane-jc/cas/compile_f2_75_125_formal_homotopy.py --prime 61 --rho 19 --y 19 --maximum-order 8 --artifact artifacts/generated-results/jc2_f2_75_125_formal_homotopy_mod61.json
.venv/bin/python plane-jc/cas/verify_common_power_carrier_wronskian.py
.venv/bin/python plane-jc/cas/test_common_power_carrier_wronskian.py
```

The first command certifies that fixed geometric degree, trivial
nonproperness, and three Newton vertices do not upper-bound
affine-normalized support, using dense triangular automorphisms.  It also
checks the Kummer descent gate for monomial-Jacobian blocks.  The second
command corrects the F2 chart to `[t,z]=-z`, proves the finite
degree/edge-halfspace band envelope, and classifies both the 35 upper zero
layers and the complete lower tail through layer `-200`.  The full record has
`2,418` jet-reduced parameters, `240` zero layers, and `1,327,026` exact
compressed generators.  On the earliest movable-double-root branch it also
derives the exact fixed-endpoint dependency cone and eliminates all eleven
`w=1` coordinates: one is normalized and the other ten have determinant
`75000`, leaving thirteen global Hermite coordinates.  Its pinned
artifact is
[`../../artifacts/generated-results/jc2_f2_75_125_character_layers.json`](../../artifacts/generated-results/jc2_f2_75_125_character_layers.json).
It is an exact B0 envelope, not an exhaustive polygon normal form or an F2
exclusion.
The third command classifies the quadratic common edge into its four
exhaustive contact partitions and audits the handoff to the
finite-normalization/log-boundary programmes.  It proves that contact orders
do not determine toroidal scales and that even the strongest naive
contact-to-ramification surrogate survives the finite-flat packet budget.
The retained contact-only warning is
[`../F2_BOUNDARY_HANDOFF.md`](../F2_BOUNDARY_HANDOFF.md).

The final attachment compiler tracks the translated terminal normal to the
original nonmonomial orders `(-25,5,12)`, constructs the six-blowup carrier,
the six-blowup principal arms and the target extraction fan, orients all five
endpoint/interior attachment slots, and uses the exact transverse order
`-2` to resolve every interior slot by two source blowups.  This gives the
corrected 19/31-component principal lower-bound trees and emits the two live global
normalization cases plus an audited rejection of the distinct-target double
case by target-valuation uniqueness.  Its optional `--candidate FILE.json` mode
rejects incomplete geometry and checks declared source class/unit/canonical,
purity, spectator, and meridian data.  Passing remains a necessary-gate
result, not a constructed Keller map; see
[`../F2_75_125_GLOBAL_ATTACHMENT_COMPILER.md`](../F2_75_125_GLOBAL_ATTACHMENT_COMPILER.md).
<!-- status-consumer: PF2GA1 57dea3062b1147fb -->
<!-- status-consumer: PF2LNP1 e4f0f231bf7494d5 -->
The carrier-profile command continues the exact log-matrix audit through the
marked carrier points and principal arms.  It constructs the common regular
fans, verifies local determinants `1`, `3`, and `5`, and strengthens the
source-boundary lower bounds to `27/48` components.  It does not compile the
upstream carrier-extraction chain, outgoing terminal tail, affine purity row,
or other global centers; see
[`../F2_CARRIER_LOG_NODE_PROFILE.md`](../F2_CARRIER_LOG_NODE_PROFILE.md).
<!-- status-consumer: PF2CLP1 41625dd5d3f8f898 -->
The carrier-Wronskian commands continue from that skeleton at the generic
carrier.  They remove the seven rational target shears before descent 36,
solve the resulting low-degree Wronskian, extract target ray `(5,36)`, and
reduce the two live cofactor strata to one rational squarefree point and two
`Q(sqrt(5))` double-root points.  Their residue maps are respectively a
cyclic cubic and the already certified terminal Belyi map; see
[`../F2_75_125_CARRIER_WRONSKIAN_CLASSIFIER.md`](../F2_75_125_CARRIER_WRONSKIAN_CLASSIFIER.md).
<!-- status-consumer: PF2CW1 a7774b0fa736b64c -->
The carrier-specialization commands perform the next exact handoff.  They
prove that the squarefree carrier is outside the descent-eight double-root
component, and specialize every exposed zero-row, target, and Hermite linear
map for the two double carriers over `QQ(rho)`.  The result pins the 53 raw
successive forcing cokernels through layer `3`, the final `7+6` global
cokernel, and the quartic compositum needed by the nonlinear forcing; see
[`../F2_75_125_CARRIER_SPECIALIZATIONS.md`](../F2_75_125_CARRIER_SPECIALIZATIONS.md).
<!-- status-consumer: PF2CS1 666da98d2d24669e -->
<!-- status-consumer: PF2NF1 cfd1da5136c0b6d0 -->
The nonlinear-forcing compiler presents all `366` geometric equations as an
exact degree-at-most-seven arithmetic circuit.  The separate modular probe
adds the `a!=0` localization equation over `F_31`, differentiates the circuit
sparsely, and identifies a consistent `169`-coordinate spacing-four
staircase tangent chart of rank `57` and affine dimension `112`.  Exact
univariate interpolation proves that the pinned particular tangent line and
six selected coordinate lines have unit raw common gcd.  After corrected
descending-column back-substitution, their projected remainders vanish.
The complete `112`-coordinate census and eight dense mixed lines are also
identically zero in the `153`-dimensional full-Jacobian cokernel.  A
fixed-Jacobian homotopy then lifts through order `16` over `F_31` and order
`8` over `F_61`; a seven-coordinate regular gauge removes a repeated
`(1-lambda)^(-2)` pivot mode.  Neither order-16 truncation specializes to a
modular point at `lambda=1`.  This is a localized good-reduction and formal
deformation experiment, not an F2 exclusion; see
[`../F2_75_125_NONLINEAR_FORCING.md`](../F2_75_125_NONLINEAR_FORCING.md).
The final command pair proves the reusable primitive common-power theorem.
For arbitrary coprime top powers `(m,n)` and degree-`k` carrier polynomial it
forces descent `k*(m+n-1)+1`, reduces the first Jacobian row to the fixed-`c`
linear kernel `k*c*D'-(k-1)*c'*D-kappa*c=0`, and computes the complete
three-point passport from the root-multiplicity partition.  It explicitly
defers the resonant `k=2` and imprimitive-multiplicity loci; see
[`../COMMON_POWER_CARRIER_WRONSKIAN.md`](../COMMON_POWER_CARRIER_WRONSKIAN.md).
<!-- status-consumer: PCW1 94b10929118f151d -->
The fourth command checks the forced chain, terminal normalization, and live
character profile `P={1,4}`, `Q={0,1,3}` modulo five.  The coarse-bridge
obstruction is in
[`../AFFINE_SUPPORT_NEWTON_BRIDGE.md`](../AFFINE_SUPPORT_NEWTON_BRIDGE.md).
The modified-system generator reproduces the published `r=2` systems, emits
the conditional 14- and 22-function `r=3` windows, eliminates the power rows
in polynomial coordinates, and presents the residue by an Artinian Fitting
ideal.  It also proves the universal endpoint-binomial section, parametrizes
the `d=2` residue as a rational fourfold, and certifies a smooth formal
coefficient-torus branch plus cubic-invariant reduction for `d=3`.  Its
uniform congruence-support gate excludes every `d=2,3` congruence section
under the certified `X^4` weight.  In the surviving `d=3,h=2` section it
eliminates the antiderivative constant and isolates the exact finite
candidate `B=1+y^2*t^6+y^3*t^9/3`, `lambda=0`, `F_-6=5*y^7/81`.  It does
not prove the common-power ansatz, `d=2,3`, or the lower Laurent-`y` ledger.
The modified-chart checker independently derives `gamma=2`, the monomial
chart, and every possible nonnegative-`xi` support from the corner chain.  It
then retains the binomial-jet relations: at `r=3` the P/Q source images have
ranks `74/83` and `196/215`, the formal terminal point violates an explicit
top-band relation, and the complete projected top-band ideal is `(1)`.  The
P-only gap algebra has length `27`; the first Q gap has nonzero
resultant/Fitting determinant in it.  This kills all branches of the literal
polynomial projection and proves that naive negative-tail deletion cannot be
the missing modified-Laurent theorem; it does not exclude the full Laurent
F2 row.
<!-- status-consumer: PF2MCB1 6ff13314e0090f52 -->
The optional `--extended-r5` flag requires Singular and proves exactly over
`QQ` that the `r=5` P-only projected top-gap ideal is also `(1)`.
The terminal checker first proves the
all-`r` degree-`2r` passport and geometric `A_(2r)` theorem.  The Kummer and
terminal commands then certify the degree-six specialization, including
primitive
four-transitive geometric `A_6` monodromy, trivial deck group,
indecomposability, zero transverse different, residue-different packet
`(4,2,2,2)`, and parameter-free residue data.  The terminal checker also
verifies discriminant `5^17*r^4*(729*r-125)^2`, arithmetic `S_6` over `Q(r)`,
the `(5,3,3)` genus-25 regular `A_6` closure, the three interior preimages of
the target toric nodes, the geometric-degree floor six, and the same-target
double-packet floor twelve.  It explicitly records that the target center is
at infinity, so no affine-sheet increment applies.  Their output reopens F2
only at the global gluing stage and does not exclude `(75,125)`.  The final
command exhausts a stated two-transposition spectator model; all six gluing
classes survive and generate `S_7`.  It also proves that the terminal
five-cycle normalizer is the generic Kummer group `AGL(1,5)`, with common
quadratic character `Q(sqrt(5))`, and audits the stronger fivefold model.
With a rational connected source boundary, that model forces degree eleven
and leaves one inertia-supported unoriented `S_11` class.  Both models remain
conditional, so neither supplies an F2
contradiction.  The
consolidated theorem and claim boundary are in
[`../F2_MODIFIED_LAURENT_FAMILY.md`](../F2_MODIFIED_LAURENT_FAMILY.md).
<!-- status-consumer: PF2GC1 33dbc5ff48b5d064 -->

The log-boundary and Poisson-square prefilters are:

```bash
.venv/bin/python plane-jc/cas/test_intrinsic_a2_boundary.py
.venv/bin/python plane-jc/cas/test_plane_boundary_exclusion.py
.venv/bin/python plane-jc/cas/test_finite_normalization_signatures.py
.venv/bin/python plane-jc/cas/verify_unibranch_spectator_models.py
.venv/bin/python plane-jc/cas/test_log_boundary_compiler.py
.venv/bin/python plane-jc/cas/test_poisson_square_rigidity.py
.venv/bin/python plane-jc/cas/log_boundary_compiler.py --frontier-72-108
```

The first reconstructs the canonical and log-canonical classes of a complete
`A2` boundary from its intersection matrix, checks `K_X^2+rho(X)=10`, and
audits target pole vectors against nefness, degree, ramification, and the
intrinsic dicritical condition.  It also verifies the sharp canonical
free-depth-three obstruction.  The second checks the residue-immersion
Riemann--Hurwitz budgets through degree eight, now for every ordered
positive puncture profile rather than only one or two punctures, and applies
the resulting one-puncture obstruction to the first numerical degree-six
package.  It also
checks that the primitive minimal-sheet formula `d=e+1` is incompatible
with the `2e` fiber length forced by a conductor identification.  These are
structural boundary regressions, not low-degree counterexample searches.
The same regression now exercises a typed finite-normalization gate.  It
checks `d=e*f+a`, refuses to infer target transfer or exhaustive affine-sheet
data from a source tree, excludes the sheet-deficient range `a<e`, and
records why the `(72,108)` ledgers `29=3+26` and `29=5+24` do not yield a
conductor contradiction.
The third command enumerates the complete coarse finite-normalization
signature atlas at fixed geometric degree: boundary rows `(e,f)`, positive
affine-sheet partitions, puncture counts, and the residual
Riemann--Hurwitz cost `f+s-2`.  It checks the forced ramified degree-three
`2+1` row, the collapse `f=s=1` under residue immersion, and Pareto
antichains through degree eight.  The same regression audits the
componentwise and global residual-different identity
`available neighbor ramification = companion intersection + 2f-2`,
detecting exposed, exactly paid, and overdrawn boundary leaves.
The intrinsic boundary regression now feeds its complete intersection
matrix, pole vector, and reconstructed ramification vector into the same
audit.  For a dicritical over a target curve of degree `c`, it recovers
`f=(Qp)_i/c` and `e=r_i+1`, rejects failure of the divisibility condition,
and infers the forced companion-sheet intersection
`M.E=available-2(f-1)`.  A negative value is a graph-level exclusion; known
companion geometry can then be checked for exact equality.  The
first free-depth-three package consequently forces companion intersection
two on its degree-one dicritical.
The same intrinsic module now contracts every `H`-null boundary curve by an
exact Schur complement, producing the Mumford intersection form on the
normal finite Stein model.  Its projection audit keeps both adjunction
corrections visible: the surface different created by contracted chains and
the normalization conductor `(c-1)(c-2)` of a rational plane image curve.
The free-depth-three example gives `E^2=-1/3`, surface different `2/3`, and
corrected companion intersection `4/3`.  On the two terminal 23-component
graphs the formal Keller-class audit contracts nine curves and gives
`E^2=33/8`; every candidate image-degree row balances after its target
conductor cost is included.  Because those terminal Newton pairs have
bracket `X^2`, these last values are comparison ledgers, not a Keller-map
exclusion.
[`target_conductor_atlas.py`](target_conductor_atlas.py) performs the next
bounded reduction.  For a rational plane curve it distributes
`delta=(c-1)(c-2)/2` among singular points, enumerates every branch count
allowed by `delta_q >= binomial(r_q,2)`, and records necessary branchwise
conductor weights.  Exhaustive degrees three through five show that the
unique Pareto-minimal profile concentrates the entire conductor at one
unibranch point.  The explicit curve `y^(c-1)z=x^c` proves this face exists
in every degree, so total conductor cannot by itself force a collision of
distinct normalization points.  The direct minimal-face regression scales
through degree 125.  Applying the already-proved residue-immersion gate
removes every singular unibranch packet.  The unique remaining Pareto
minimum is one two-branch singularity, with complexity `(1,2,1,2)`;
arbitrary conductor remains possible through branch tangency, but the two
normalization points now activate the finite-flat packet inequality
`d>=2e`, equivalently `a>=e` for one residue-degree-one boundary row.
The unibranch-spectator replay tests the proposed direct local exclusion.
For every `n>=3` it verifies the finite-free rank-`n+1` map
`(T,u)->(u,T^(n+1)-T^n+u*T)`, whose singular unibranch fiber has a
length-`n` boundary point and a reduced étale spectator.  The packet
saturates Orevkov's Euler budget.  Its exact failure is global: deleting the
ramification curve gives `A1 x G_m`, with a nonconstant unit, rather than
the Keller open `A2`.  See
[`../UNIBRANCH_SPECTATOR_COUNTERMODELS.md`](../UNIBRANCH_SPECTATOR_COUNTERMODELS.md).
The next command turns certified monomial branch
scales into regular toroidal
blowups, a proximity graph, complete boundary and intersection matrices,
valuation/different/conductor labels, and chart-aware Smith invariants.  It
extracts the local `(2,1),(3,1),(4,1)` rays proved by the `(72,108)` Laurent
case tree and separately audits the longer map-base ideals
`(t,x^4),(t,x^6),(t,x^8)`, compiling their isolated nested source chains of
lengths `4,6,8`.  Additive composition proves that cases `a,b,c` share the
same eight-blowup graph.  `F_4` absorbs the final involution and records its
swap of the two base divisors; filling pre-transition `Xinf`, which is
post-transition `X0`, gives a unimodular `10 x 10` affine-plane boundary
passing the intrinsic adjunction/Noether audit.  The unselected order-three
factor is a unit at the common order-four center and avoids the filled
post-transition `X0` divisor.  Target-infinity pole orders are exact on all eight
exceptionals.  Together with the original-boundary orders they give the full
common-graph pole vector `(1,24,1,1,1,12,9,6,3,2)`.  Its intrinsic audit has
no dicritical component, proving that at least one additional global cluster
is required.  The first weighted-Wronskian equation forces the actual source
cluster at `E3 intersect E4`: one pole-three exceptional and ten simple
pole-two children.  The low monomial excludes the numerical smooth-`E3`
candidate.  Terminal Case 2 then forces a `12,0` target cluster and a
degree-twelve dicritical; the combined 23-component boundary passes all
intrinsic gates with remaining self-intersection `29`.  The final upper edge
forces `A=a*r^2,C=c*r^3` with quartic `r`.  Its five multiplicity partitions
form an exhaustive edge-only comparison family.  The compiler emits every
regular fan, full boundary and intersection matrix, target pole and
ramification vector, normalization different, and source/image conductor.
All five packages pass with remaining self-intersection `29`.  The primary
split-factor formula then forces `r=(s-beta)^4`; the alternate legal factor
supplies an exact transverse chart, so Terminal Case 1 selects and resolves
the same 23-component `(4)` package as Terminal Case 2.  The legacy aggregate
record is rejected only because the generic IR cannot serialize the
nonmonomial first-block cluster.  Because the final pair has `[P,Q]=X^2`,
the compiler now corrects the boundary-supported `K+3H` representative by
`div(X^2)`: the actual dicritical normal indices are `3` in Case 1 and `5`
in Case 2, while the total ramification intersection remains `35`.  The
a priori exact residue cover degrees are `1,2,4`.
[`audit_case2_residue_strata.py`](audit_case2_residue_strata.py) excludes
the Case-2 degree-two and degree-four polynomial-composition strata exactly,
using neither `J0` nor the residual `J1` compatibility equations.  Thus
only the degree-twelve row remains at that stage.
[`case2_infinity_resolution.py`](case2_infinity_resolution.py) then
localizes the seven residual `J1` compatibility cubics at the forced
endpoint `G_12 != 0`; the resulting exact ideal is the unit ideal, without
`J0`.  Terminal Case 2 is therefore excluded.  Before compatibility, its
generic infinity branch has characteristic `(4,13)` and a seven-ray regular
toric resolution, also recorded by that audit.  The compiler also rewrites
the Case-2 bottom equations as
`B=K*c,F=K*g`,
`2*H*(A*g-c*E)+K^2*(c*g'-c'*g)=0`, where
`H=gcd(C',G')`.  The first coefficients force `t|H`, so the degree-zero
gcd stratum is already excluded.  If `deg(H)=1`, the only surviving origin
orders are `ord(B),ord(E),ord(G'),ord(F)=(1,2,3,3)`.
At the opposite end,
[`audit_case2_maximal_gcd.py`](audit_case2_maximal_gcd.py) excludes
`deg(H)=7` exactly: three low coefficients of `remainder(G',C')` and the
terminal `t^19` coefficient of `J0` generate the unit ideal, without any
residual `J1` compatibility equation.
[`audit_case2_gcd6.py`](audit_case2_gcd6.py) similarly excludes
`deg(H)=6`: write `C'=H*(t+v)` and use `C'(0),H(0)`, the last two
coefficients of `G' mod H`, and `J0` at `t^19`.  The standalone exact audits
are:

```bash
make verify-plane-case2-residue-strata
make verify-plane-case2-j1-endpoint
make verify-plane-case2-maximal-gcd
make verify-plane-case2-gcd6
```

The third proves that the entire
geometric reduced weighted-tangent three-layer support box has exactly three
components: the cubic tangent-pencil closure and the explicit `C=0` and
`A=0` components.  Exact transverse families and tangent ranks show that
their generic multiplicities are respectively `2,3,1`.  Their pairwise
reduced intersections have branch counts `(2,2,1)`: a common lower-Wronskian
core plus one constant-`D` tangent/`C=0` branch and one constant-`B,C`
tangent/`A=0` branch.  Dense-chart scheme tangent dimensions are `8,7,6`;
chosen tangent-kernel slices on the two extra branches both have length five
but Hilbert/socle data `(1,4;4)` and `(1,3,1;1)`.  The exact saturation
identity `I:d0^infinity=I` proves that no associated prime is supported on
`d0=0`; a `G_m` normalization reduces the primary problem to `d0=1`.
An exact `d3,d2` colon filtration then gives the complete associated-prime
set: three minimal components, the three intersection surfaces, and the two
core/intersection curves.  Their normalized dimensions are `(3,3,3)`,
`(2,2,2)`, and `(1,1)`.

An optional exact Singular audit independently computes the radical and its
three minimal primes, and checks that every pairwise and triple reduced
component intersection has dimension three.  A second fast audit verifies
the principal-chart saturation and reconstruction identities:

```bash
Singular -q plane-jc/cas/poisson_square_radical.sing
Singular -q plane-jc/cas/poisson_square_primary_charts.sing
Singular -q plane-jc/cas/poisson_square_separator_primary.sing
Singular -q plane-jc/cas/poisson_square_normalized_defect.sing
.venv/bin/python plane-jc/cas/test_poisson_square_filtered_modules.py
```

The third command proves that generic primary closures leave a genuine
separator defect and prints exact torsion witnesses.  The fourth decomposes
`I:d3`, `(I+(d3)):d2`, and `I+(d3,d2)` and certifies exactly eight associated
primes.  It also resolves the second gluing layer into two primary curve
components and computes exact transverse Hilbert vectors on all associated
strata.  The final command checks the reusable lower-band filter, including
localized `preserved/cut/eliminated` decisions.

The reusable superelliptic leading-block reducer is documented in
[`../SUPERELLIPTIC_DERHAM_ENGINE.md`](../SUPERELLIPTIC_DERHAM_ENGINE.md).
It includes exact Gauss--Manin matrices and cyclic-vector scalar
Picard--Fuchs extraction.  Its fast regression is:

```bash
python3 plane-jc/cas/test_superelliptic_derham.py
```

The source-aware compiler IR and the explicit `(72,108)` tail-basis
certificate are documented in
[`../NEWTON_DERHAM_COMPILER.md`](../NEWTON_DERHAM_COMPILER.md).  Run:

```bash
python3 plane-jc/cas/test_newton_derham_compiler.py
```

This also verifies that the source-excluded repeated-tail `(96,144)` row and
the incomplete `(75,125)` **coefficient record** are rejected rather than
assigned invented Laurent bands.  The latter route-specific refusal is
compatible with the separately certified F2 terminal target row.

The exact forced F2 `j=1` skeleton and its machine-readable residual
obligations are tested separately:

```bash
.venv/bin/python plane-jc/cas/classify_f2_75_125_layers.py
python3 plane-jc/cas/test_f2_75_125_frontend.py
python3 plane-jc/cas/f2_75_125_frontend.py
```

The first command replays the exact full B0 envelope and the corrected
top-tangent profile: the first-five kernels are `6,6,7,7,10`, while the
formal `C0^(-1)` resonance is at layer 10 and is not a source-band kernel.
The same replay uses the next nonlinear rows to force source-root continuation
through descent 7 and isolate `27*y^2-9*y+1=0` at descent 8.  The Q-band-one
normalization excludes its fixed Kummer supports; at this earliest spacing
only the nonzero double root of `R` remains, and it passes an exact local target-jet interpolation.
The descent-40 fifth multiple is recorded as a lower-tail Fitting row, not as
the invalid primitive equation `E5=0`.  The complete target cokernel is
split into twelve local jets and two triangular residues, while layer zero
becomes a rank-14 quotient of an exact length-15 Artinian algebra.  The
lowest-`u` target edge has an exact Bezout witness.  Raw old-band generators
span both cokernels.  The edge witness is also tested uniformly in `r`: it
completes as an exact formal shear, but the shear has an unavoidable infinite
tail.  The exact polynomial second-order repair never terminates
quadratically, and the `r=3` cubic and quartic termination ideals are units.
The first Kummer-return band `v^5` is now reconstructed from the exact source
jets.  Its fifth-binomial correction cancels the former terminal conflict;
the remaining `8 x 10` map has unit minor `3^5*5^16*e^13` after base change
to the rank-two descent-eight algebra.  Both branches survive, and a second
unit Bezout minor proves that the standalone edge recursion remains
surjective through `v^10`.  The replay now gives the two all-`r` unit-minor
formulas and traces their `18*r-1` pivots to exact original-polynomial source
combinations strictly inside both certified supporting edges (`53` at
`r=3`).  Its confluent-CRT determinants quotient the
triangular `w=0` control block, leaving a rank-`24` global Hermite module over
the rank-two candidate algebra.  The fixed `w=1` block is then one normalized
identity plus a `10 x 10` affine-linear block of determinant `75000`.
Exact `w=0`-preserving followers eliminate those ten variables, leaving
thirteen cokernel coordinates.  The carried reduction is implemented by
`reduce_f2_75_125_endpoint_system.py`: its ten solutions have `1,489`
straight-line terms, the resulting bracket circuits have degree at most
eight, and `1,061` active source coordinates remain.  The endpoint-disjoint
new-Q operators on layers `39..29` are tridiagonal with unit minors; they
eliminate `134` coordinates and leave `219` upper Fitting slots.  The exact
coupling boundary is descent `12` (layer `28`), where P3/Q13 first enters the
endpoint solution.  Later
<!-- status-consumer: PF2ER1 64378dad616fc3f2 -->
first-defect spacings `9..90` remain explicitly enumerated.
The final command emits JSON.  Its `frontend_complete` field is intentionally
false until the B0 over-envelope has been cut to the actual lower Newton row
or eliminated directly.

The source reconciliation for the repeated-tail `(96,144)` row is:

```bash
python3 plane-jc/cas/test_frontier_96_144_source_audit.py
python3 plane-jc/cas/frontier_96_144_source_audit.py
```

It proves `q1=d0=4`, reduces the vertical residual factor to a cubic, removes
the two root partitions containing a simple root, and passes the remaining
triple-root factor to:

```bash
python3 plane-jc/cas/test_complete_chain_no_escape.py
python3 plane-jc/cas/complete_chain_no_escape.py
```

The latter reproduces the published companion final corner and proves that
the triple-root edge has open-chain counts `1,6,3,0`, hence no complete-chain
escape.

To compile both audited Proposition 4.3 polygons through lattice supports,
the Laurent chart, all upper bracket layers, and the genus-three first block,
run:

```bash
python3 plane-jc/cas/test_laurent_band_frontend.py
```

The exact 90 MB Zenodo attachment and its extracted source snapshot are pinned
in the repository at
`plane-jc/external/zenodo-21479814/jc2-72-108-exact-certificates-v1.0.1.zip`
and `plane-jc/external/zenodo-21479814/bilLkarkariy-jc2-72-108-exact-certificates-d9ea4fd/`.
First check the attachment:

```bash
md5 plane-jc/external/zenodo-21479814/jc2-72-108-exact-certificates-v1.0.1.zip
shasum -a 256 plane-jc/external/zenodo-21479814/jc2-72-108-exact-certificates-v1.0.1.zip
```

Expected values:

```text
MD5    91255150c8c689b26dc6fb61f9d80aec
SHA256 f7f0876de12d35badbed2be6a773d4a9dada50aff778c080126aab541deefcde
```

After extracting the outer archive and its
`release_bundle/jc2_72_108_exact_replay_v1.0.1.zip`, create a Python
environment from the archived `requirements.txt` and run from
`release_bundle/exact_replay`:

```bash
PYTHON=/absolute/path/to/venv/bin/python ./verify_all.sh
```

Expected final marker:

```text
JC2_72_108_EXACT_REPLAY_PASS
```

The extracted archive is immutable provenance.  A local portable CPU adapter
checks both archived hash manifests, removes two duplicate calculations, and
runs the independent branches with at most four concurrent processes:

```bash
/absolute/path/to/venv/bin/python \
  plane-jc/cas/verify_72_108_exact_fast.py --jobs 4
```

Equivalently, when `PYTHON` names that FLINT-enabled environment:

```bash
make verify-plane-72-108-exact-fast
```

It uses the archived equations and verifiers without changing the pinned
release.  Its hard-certificate adapter parses large decimal integers directly
with GMP and performs the formal branch transport after one, rather than two,
replays of the 89 MB identity.  `--jobs 1` retains the same optimized
calculation graph without concurrency.

Engineering benchmark on 2026-07-30: the archived serial driver took
`64.02 s`, while this adapter took `22.51 s` with `--jobs 4` on an AMD Ryzen
9 9950X3D under WSL2, using Python 3.12.3, gmpy2 2.3.1,
python-flint 0.9.0, and SymPy 1.14.0.  This is a runtime measurement, not a
new certificate or mathematical result.  Apple M2 timing remains to be
measured with the same command.

Optional accelerators remain a possible future backend for bounded modular
searches, batched finite-field calculations, or candidate discovery.  They
are deliberately not dependencies of the exact replay: any accelerated
experiment must emit a reproducible artifact with a portable CPU verifier.
No CUDA- or Metal-specific target is currently maintained.

For the independent hard-certificate check, run from the repository root:

```bash
/absolute/path/to/venv/bin/python \
  plane-jc/cas/verify_h_certificate_independent.py \
  /absolute/path/to/exact_replay/hne0_polred.pkl \
  /absolute/path/to/exact_replay/hard/h_certificate_exact.txt
```

This checker shares only the serialized four-polynomial input and the text
certificate with the primary computation.  It does not import the generating
code or a CAS.

For the primary Case-2 unit ideal, the following independent checker replaces
the standard-basis output by a projective Cramer/resultant proof:

```bash
.venv/bin/python plane-jc/cas/verify_case2_resultant_proof.py
```

The support shape reduces the four residuals to degree-eight univariate
eliminants on two charts.  Extended-gcd identities at the good place
`p=101,u=55` certify that their characteristic-zero resultants are nonzero;
the same checker excludes the singular Cramer branches and the origin.

The larger characteristic-zero Bézout identity remains as an independent
fallback:

```bash
.venv/bin/python plane-jc/cas/verify_case2_syzygy_independent.py
```

It pins the serialized certificate, verifies the stored degree-35 field
polynomial is irreducible, and checks
`1=T_0 R_0+T_1 R_1+T_7 R_7+T_9 R_9` coefficient by coefficient without
importing the equation generator, `exact_core`, or Singular.  The mathematical
argument and Macaulay dimensions are recorded in
[CASE2_EXPLICIT_SYZYGY_PROOF.md](../CASE2_EXPLICIT_SYZYGY_PROOF.md).

To verify that the hard membership identity and the specialized \(h=0\)
identity compose to a single unit certificate for the seven pre-division
Case-1 equations, run with the archive's FLINT-enabled Python environment:

```bash
/absolute/path/to/archive/venv/bin/python \
  plane-jc/cas/verify_case1_unit_composition.py \
  /absolute/path/to/exact_replay
```

This checker does not expand the combined multipliers.  It verifies their
factored straight-line representation: all six elimination lifts, the
invertible degree-five descent, and the lift of the \(h=0\) identity to
`1 = sum(A_i G_i) + B h`.

## Critical pinned inputs

| Artifact in `exact_replay/` | SHA-256 |
| --- | --- |
| `firstblock_Q_exact.sing` | `90fb933bf4ae75accddecb69993957db8e289b44c1ae2285dc91b9130f30c062` |
| `firstblock_Q_exact.out` | `2e965d03b39d87531228943cd634de4438d3311fb7a0660dd6dc43a768ae05cb` |
| `case2_compact4_exact.sing` | `38694ed8e9e3b9256b380edf882abefff0ad944113cf0d98e4a951e8f4e31030` |
| `case2_exact_certificate.pkl` | `cfbc3c39d7a28013671144f43ef76f0498542eaf6d562dd624bba3311194e4aa` |
| `case1_residuals_exact.txt` | `f026228c422a213eed853f684b0e2fd98cbc51338ee11c3254a77bf053957c2b` |
| `case1_branch1_after_w_eqs.txt` | `dd6d122161388b8a1961c7e52e92b990e62d0e95c15acc5db56423a1cdf8ea44` |
| `case1_branch2_after_w_eqs.txt` | `0f7c9a8a01a18d725e2a9dc663eec0e4ea018ecb1b1959de1ae7be0214218692` |
| `case1_branch1_after_w.pkl` | `368a1dafdb6d26708b85d652a437c848a7676ba1419e4575ae74186f022621b9` |
| `h0_branch1_exact_certificate.pkl` | `664de005e99bc6a0e61ba479ba64ed57ddbfa9c5399b955cbb12237ac70f8186` |
| `h0_exact_certificate.pkl` | `d7844c2f8edea62be4e3a7fe8a160dc4b7b70efd1262993ec6d2ed7e78722a1a` |
| `hne0_deg35.pkl` | `082471d05a2a7ceebca9fd3a615d8fb6fddaee8ff80afae161a043f71edcb575` |
| `hne0_polred.pkl` | `5a6e423d74ef09fc9c7a7282c500bda566018d7e56a93124665796bbe417cedf` |
| `fixed_matrix_p71.npz` | `3fc0a958672d361a343e8fbeae77c1012f8a5bb4b8e1aa3fb9acdb670c8726dd` |
| `pivot_scalar_rows_p71.npy` | `be4215d7f0303e54890494f7ee03936eff0f6f3e416fe274d7c8cf579938190d` |
| `hard/h_certificate_exact.txt` | `0e48ffab32469ef8405a6945b16cf1521ddeb3c592ae4e5051968110a4dc656a` |
| `verify_all.sh` | `f40bacfa84d915b375b4158109d5120b2d964a8d8012be770750d310abfb5837` |

The exact base rings, variables, orders, localizations, and identity forms are
listed in [the reproduction note](../PAIR_72_108_REPRODUCTION.md).
