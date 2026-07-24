# Exact replay commands

The log-boundary and Poisson-square prefilters are:

```bash
.venv/bin/python plane-jc/cas/test_intrinsic_a2_boundary.py
.venv/bin/python plane-jc/cas/test_plane_boundary_exclusion.py
.venv/bin/python plane-jc/cas/test_finite_normalization_signatures.py
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
Its fast regression is:

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
the incomplete `(75,125)` row are rejected rather than assigned invented
Laurent bands.

The exact forced F2 `j=1` skeleton and its machine-readable residual
obligations are tested separately:

```bash
python3 plane-jc/cas/test_f2_75_125_frontend.py
python3 plane-jc/cas/f2_75_125_frontend.py
```

The second command emits JSON.  Its `frontend_complete` field is intentionally
false until the lower Laurent boundary has been classified exhaustively.

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
