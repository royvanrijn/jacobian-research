# Boundary-package compiler: stage-one prototype

## Status and outcome

This is a research prototype, not a realization theorem.  It reverses the
usual repository workflow by taking finite-normalization boundary data as
input and checking whether that data survives exact necessary conditions for
an affine-source polynomial Keller map.

The executable prototype is
[`boundary_package_compiler.py`](../scripts/boundary_package_compiler.py).
Its first benchmark suite establishes the following computational facts.

1. A degree-four \(A_4\) package with two selected index-two boundary primes
   and a three-punctured rational selected curve passes every implemented
   stage-one condition.
2. The degree-seven Fano action of
   \(\operatorname{GL}_3(\mathbb F_2)\), with two selected index-two primes
   over its involution branch, also passes.
3. A rational normalization \(\mathbb G_m\) with a nontrivial node pairing,
   finite-index unit lattice, and conductor-preserving involution passes.
4. A degree-three genus-one selected boundary with six simple branch values
   passes the exact monodromy, Riemann--Hurwitz, unit, and adjunction gates.
5. Conductor non-equivariance, ramification colored as lying in the affine
   Keller source, a coefficientwise determinant-ledger mismatch, and filling
   a hole of the actual value semigroup by saturation are rejected as
   intrinsic obstructions.
6. Feeding the exact symbolic-cone branch datum back into stage one gives one
   global \(L\)-prime with \((e,f)=(2,2)\).  It passes the degree,
   cycle-profile, and Riemann--Hurwitz checks and is rejected solely because
   that ramified prime is affine-colored.
7. The compiler now has a second affine-source lattice gate for any certified
   affine factorial core with a based unit lattice.  It reproduces the balanced
   wild-boundary Smith factors in `N=2,3,5,6,7`, rejects the stabilized `N=5`
   block by class-group torsion `Z/4`, and computes the exact order four of
   its named class `[L1]`, despite constant units and a unimodular synthetic
   projective-boundary input.
8. The factoriality restriction can be replaced by a finite presentation of
   the core class group plus explicit lifts of its relations.  The resulting
   block presentation detects nonsplit extensions: the regression with
   `V=(2)`, `R=(2)`, and correction `A=(1)` has class group `Z/4` and a lifted
   core class of exact order four, whereas `A=(0)` gives `Z/2 + Z/2`.
9. A shared retained-root Euler gate now rejects a certified balanced
   squarefree retained polynomial of degree `r>1` before source
   reconstruction.  It records
   `[U]=L^2+(r-1)L` and `chi_c(U)=r`.  The regression rejects `r=4`,
   passes `r=1`, and leaves missing or nonsquarefree hypotheses explicitly
   `uncertified` or `not_applicable` rather than guessing them.
10. A shared conductor/contact-loss audit now turns a proposed lower-band
    truncation into exact branchwise inequalities.  The legacy scalar form is
    `available >= conductor + derivative loss + pole loss + other loss`; the
    preferred form records every named input-to-output path and requires only
    `available(input) >= conductor + path loss(input,output)`.  Passing data
    certify invariance of the conductor matching cokernel and distinguished
    residue class.  Short or incomplete data remain non-obstructing
    `insufficient` or `uncertified` rows.
11. An optional colored-fan front end compiles masks, target pullbacks,
    derivative and conductor divisors, boundary colors, unit lattices, class
    groups, and bounded affine-modification choices into one valuation
    matrix.  It exactly replays the \(A_4\), \(D_5\), and Davenport ledgers,
    including both exceptional \(D_5\) rows.  Its feasible output carries an
    explicit nonlinear residue and never asserts affine completion.  See
    [`TOROIDAL_BOUNDARY_FEASIBILITY.md`](TOROIDAL_BOUNDARY_FEASIBILITY.md).

Passing has status `unknown`, never `realized`.  The four surviving packages
still need a root equation, reconstruction functions, and a proof that the
regular-reconstruction open is affine space.

## 1. Why the package has two boundary dimensions

For a quasi-finite map

\[
 F:\mathbb A^3\longrightarrow\mathbb A^3
\]

the canonical Zariski--Main boundary primes are surfaces.  The rational,
nodal, or elliptic curves discussed in this experiment are selected strata
or canonically extracted curve quotients of those surfaces.  The input
therefore keeps separate:

- a `LocalPrime`, carrying \((e,f)\), affine/boundary color, and selection;
- a `SelectedCurve`, carrying normalization genus, punctures, conductor,
  units, and adjunction.

This avoids treating a curve model as though it were itself a height-one
prime of a threefold normalization.

## 2. Package schema

The current `BoundaryPackage` has the original stage-one proof blocks plus
optional retained-root, conductor-jet, and toroidal blocks.  A separate
optional stage-two block carries realization references.

### Toroidal feasibility

`ToroidalBoundaryDatum` records a certified fan skeleton, boundary colors
carried by its rays, and a valuation matrix with one row per color.  Fixed
divisor identities are checked coefficientwise.  A
`ColoredDivisorSpanProblem` asks whether a target colored divisor lies in the
integral column span of named generator functions.  The compiler computes
the rational-span rank jump, the exact order of the target class in the
cokernel, and primitive proportional-row witnesses.  A nonzero witness

\[
 b\tau_i-a\tau_j
 \quad\text{when}\quad A_i=a p,\ A_j=b p
\]

rules out rational, integral, and nonnegative generator exponents at once;
this is an unbounded obstruction, not a box search.  It is theorem-bearing
only when the input marks the declared generator architecture exhaustive and
supplies its scope certificate.  See
[`COLORED_DIVISOR_SPAN_OBSTRUCTION.md`](COLORED_DIVISOR_SPAN_OBSTRUCTION.md).

Separately, bounded integral variables can encode modification exponents,
deletion choices, or support-function inequalities; the compiler enumerates
the declared box exactly and returns its Pareto-minimal models.  An empty
bounded problem is an obstruction only when its input marks the search
conclusive and supplies an exhaustive-scope certificate.  Such a search can
serve as a regression for a span obstruction, but it cannot replace the
unbounded theorem.

Selected submatrices use the existing Smith implementation for unit and
class-group screens.  Primitive rays and smooth cones are checked, while the
fan-incidence certificate remains theorem-bearing input.  Outside a
separately certified normal affine toric model, trivial units and class group
are necessary rather than sufficient for affine space.  Every feasible audit
therefore carries its unresolved nonlinear residue forward.

### Cover data

`CoverDatum` records:

- the degree \(N\);
- a named transitive permutation group;
- one explicit inertia permutation for every branch divisor;
- the primes above it over the declared ground function field, with
  \((e,f)\); repeating each ramification index \(f\) times recovers the
  geometric inertia-cycle profile;
- product-one data for compact curve covers; and
- commutation relations for inertia around intersecting SNC divisors.

For a prime with arithmetic residue degree \(f\), the geometric cycle
profile contains \(f\) cycles of length \(e\).  The compiler checks

\[
 \sum_j e_jf_j=N
\]

and equality with the cycle profile.  It computes the generated group rather
than trusting the group name.

Every prime with \(e>1\) must be boundary-colored.  A polynomial Keller map
is etale on its affine source, so a ramified affine-colored prime is an
immediate contradiction.

### Selected curves

The curve block records

\[
 p_a(C)=g(\widetilde C)+\delta(C)
\]

and tests adjunction

\[
 C\cdot(C+K_S)=2p_a(C)-2
\]

in a declared surface intersection lattice.  In particular, a rational
curve with one node has arithmetic genus one and is tested against the same
adjunction number as a smooth elliptic curve.

Puncture units are represented by a free lattice.  Non-torsion conductor
characters add integral linear equations.  Torsion conductor characters add
congruences and may replace the unit lattice by a finite-index sublattice
without changing its rank.

The conductor also contains explicit point pairings and automorphism actions.
Every declared automorphism must preserve the set of unordered conductor
pairs.

### Determinant ledger

For every recorded valuation, the compiler checks

\[
 \operatorname{ord}(J_\alpha)
 r\,\operatorname{ord}(\Delta\circ\alpha)
 =
 \operatorname{ord}(J_\beta\circ F).
\]

This is the coefficientwise form of the boundary-cancelled ledger from
[`CONTROLLED_BOUNDARY_SUSPENSIONS.md`](../cancellation/CONTROLLED_BOUNDARY_SUSPENSIONS.md).
It is distinct from the automatic tame finite-cover identity

\[
 (e-1)+1=e.
\]

For the \(A_4\) fixture the three rows are the exact pure-target ledger

\[
\begin{array}{c|ccc}
 &\operatorname{ord}\det D\Phi&
   \text{auxiliary contribution}&
   \operatorname{ord}\mathcal B(\Phi)\\ \hline
 W&2&1&3\\
 K&3&0&3\\
 L&1&1&2.
\end{array}
\]

These are the identities proved in
[`A4_PURE_TARGET_LEDGER_LIFT.md`](A4_PURE_TARGET_LEDGER_LIFT.md).  They make
the lift log-Keller; they do not make it an ordinary Keller map.

### Retained-root Euler prefilter

`RetainedRootEulerDatum` is accepted only when the package explicitly names
certificates for:

1. the balanced chart
   `P=xu,T=x^2u,Q=A(x^2u)(u-x^(N-1))`;
2. different support on the declared omitted fierce boundary;
3. squarefree nonzero retained root fibres; and
4. exactly one omitted fierce boundary.

The shared evaluator then emits one of `not_applicable`, `uncertified`,
`passes`, or `obstructed`.  In the applicable case it records

\[
 \#U(\mathbb F_q)=q^2+(n_q(A)-1)q,\qquad
 [U_{\bar k}]=L^2+(\deg A-1)L,\qquad
 \chi_c(U_{\bar k})=\deg A.
\]

Thus `deg(A)>1` is an affine-plane obstruction, while degree one is neutral
at this gate and proceeds to the class-lattice tests.  The compiler does not
infer squarefreeness, the chart, or the omitted divisor from a coarse
monodromy ledger.  Nonsquarefree collision presentations remain live.

### Conductor/contact-loss truncation

An optional tuple of `ConductorBranchJetDatum` records, on every completed
normalization branch in the finite conductor support,

```text
(conductor exponent, derivative loss, pole loss, other contact loss,
 available jet order)
```

together with certificates for the conductor, expression tree, and
band-to-normal valuation.  The theorem in
[`CONDUCTOR_JET_TRUNCATION.md`](../plane-jc/CONDUCTOR_JET_TRUNCATION.md)
proves that the conductor quotient, the finite matching-map cokernel, and its
distinguished residue class are unchanged when

\[
 n_i\ge c_i+d_i+\ell_i+\epsilon_i
\]

on every branch.  The emitted `conductor_jet_truncation` status is
`not_declared`, `uncertified`, `passes`, or `insufficient`.  Neither
`uncertified` nor `insufficient` is a Keller obstruction; those statuses say
that the omitted bands must still be recovered or certified.

The preferred `ConductorBranchSensitivityDatum` ledger replaces the branch
maximum by named inputs, named matching or residue outputs, and certified
expression trees.  If `lambda_(i,alpha,j)` is the loss on paths from input
`j` to output `alpha`, it asks only for

```text
available(i,j) >= conductor(i) + lambda(i,alpha,j).
```

Inputs absent from a certified complete dependency graph are reported as
unused and require no jet.  Every failing pair carries its exact deficit, so
a search can reconstruct only that input band.  Available orders can be
compiled from a certified normal-valuation vector and a finite valuation
frontier of omitted Newton exponents.  Its certificate must prove that every
omitted monomial has order at least the displayed minimum; on a Laurent cone,
a coordinatewise antichain alone need not do so.  The scalar and
dependency-sensitive inputs are mutually exclusive; the scalar schema remains
supported for existing certificates.

### Pole box and value semigroup

`PolarLedger` records:

- every divisor on which a reconstruction coordinate may have a pole;
- whether that divisor is outside the affine reconstruction open;
- the exact valuation vector and a finite pole bound for every coordinate;
- generators for the actual affine value semigroup; and
- required membership or nonmembership tests.

The compiler rejects a pole exceeding its bound, a pole on an affine-colored
divisor, and a false semigroup-membership assertion.  Membership is solved as
an exact nonnegative integer-combination problem.

The regression semigroup

\[
 \langle(2,0),(1,1),(0,2)\rangle
\]

has saturation \(\mathbb N^2\), but \((1,0)\) is a hole.  A package requiring
\((1,0)\) to be an actual value is rejected.  This is the smallest executable
guard against silently replacing a chart algebra by its normalization.

The ledger also names polar-completeness, initial-presentation, and
Rees-strictness certificates.  The prototype currently checks that the
references are present; replaying their algebraic content belongs to the
next implementation layer.

### Affine-source necessary data

There are two distinct localization matrices.  For a normal finite
normalization \(\bar X\) with distinguished open \(U=\mathbb A^n\), the
projective-complement localization sequence forces the boundary-class map

\[
 \bigoplus_i\mathbb Z[E_i]\longrightarrow\operatorname{Cl}(\bar X)
\]

to have trivial kernel and cokernel.  The prototype checks a supplied square
boundary-class matrix has determinant \(\pm1\), and rejects declared
nonzero source unit or class-group rank.

When the reconstructed open itself contains a certified dense affine
factorial open `W=Spec(A)`, with `A` a UFD and hence `Cl(W)=0`, there is a
dual and often sharper test.  Put
`M=Gamma(W,O_W)^*/k^*`, require an explicit free basis, and let
`D_1,...,D_r` be all codimension-one primes in `U-W`.  The unit-valuation map

\[
 V:M\longrightarrow\mathbb Z^r
\]

fits into

\[
 0\longrightarrow\Gamma(U,\mathcal O_U)^*/k^*
 \longrightarrow M\xrightarrow{V}\mathbb Z^r
 \longrightarrow\operatorname{Cl}(U)\longrightarrow0.       \tag{1}
\]

The new `FactorialCoreDatum` therefore names the unit generators, every
boundary prime, their full valuation matrix, and proof references for
normality, factoriality, and codimension-one completeness.  A torus is the
important special case `W=G_m^n`.  The dependency-free compiler computes the
Smith factors as gcds of minors.  A free kernel is a unit obstruction; a free
or torsion cokernel is a Weil-class obstruction.  A `CoreClassQuery` also
computes the exact order of a represented reflexive class `delta`: infinite
if `[V|delta]` raises rank, otherwise the quotient of the top determinantal
divisors of `V` and `[V|delta]`.  The exact theorem and a second SymPy
implementation are in the
[boundary-lattice prefilter](../plane-jc/BOUNDARY_LATTICE_PREFILTER.md#dual-torus-core-localization).

For the balanced wild block

\[
 V_N=\begin{pmatrix}1&0\\N-2&N-1\end{pmatrix},
\]

the kernel is zero and the cokernel is `Z/(N-1)`.  Adding an identity torus
direction changes the Smith diagonal to `(1,1,N-1)` and preserves the
obstruction, which is the three-dimensional regression used here.  The
query vector for `L1` has exact order `N-1`; the vector for `div(x)` has
order one, providing a positive membership control.

The more general `PresentedCoreDatum` does not require `Cl(W)=0`.  It records
a finite presentation

\[
 \mathbb Z^s\xrightarrow{R}\mathbb Z^c
 \longrightarrow\operatorname{Cl}(W)\longrightarrow0,
\]

lifts of the `c` generators to `U`, and the boundary-valuation matrix `A` of
rational functions witnessing the `s` relations.  Divisor localization then
presents the reconstruction class group by

\[
 \operatorname{Cl}(U)=\operatorname{coker}
 \begin{pmatrix}V&A\\0&R\end{pmatrix}.             \tag{2}
\]

Here `V` is still the complete unit-valuation matrix, and
`Gamma(U,O_U)^*/k^*=ker(V)`.  The datum therefore requires separate
certificates for normality, the presentation of `Cl(W)`, the lifted relation
witnesses, and completeness of the codimension-one complement.  Its named
class vectors have `r+c` coordinates and use the same exact determinantal-
divisor order formula.  The correction matrix cannot be omitted: the two
one-relation examples `V=R=(2)` with `A=(0)` and `A=(1)` have the same
boundary and core quotients but total class groups `Z/2 + Z/2` and `Z/4`.
The schema also permits zero unit columns and zero relation columns, so cores
with only constant units or with free class-group summands are not excluded.

<!-- status-consumer: BL1 e86cdcd66993bccc -->
<!-- status-consumer: PWB7 19f4f4ffc96227a3 -->
<!-- status-consumer: CJT1 afb70f90ff10f3d7 -->

The qualifier “necessary” can be removed only in a narrower toric category:
if the dense-torus action extends and `U` is a normal affine toric variety,
then trivial class group gives `A^r x G_m^(n-r)`, and constant units force
`A^n`.  The present compiler does not infer an extended torus action from a
Laurent open, so its general verdict remains a filter.

This is only a necessary filter.  Trivial units, class group, and canonical
ledger do not characterize affine space; the motivic obstruction in
[`ORIENTED_CUBIC_COX_CHART.md`](ORIENTED_CUBIC_COX_CHART.md) is a concrete
counterexample.

### Stage-two certificate

The status `realized` is unavailable from this stage-one compiler.  A
`StageTwoRealizationCertificate` must name all of:

1. a root equation;
2. completed local-factorization certificates;
3. two-sided reconstruction identities;
4. a polynomial-ring isomorphism for the regular-reconstruction open;
5. a constant-Jacobian certificate; and
6. a monodromy certificate.

An incomplete certificate leaves the result `unknown`.  Complete references
also remain `unknown` until a future symbolic verifier has replayed their
contents; strings naming alleged certificates are not themselves a proof.

## 3. Exact monodromy benchmarks

### Tetrahedral degree four

Use

\[
 \sigma_0=(12)(34),\qquad
 \sigma_1=(234),\qquad
 \sigma_\infty=(142).
\]

Their product is one in the convention used by the checker, they generate
\(A_4\), and their ramification contributions are \(2,2,2\).  Hence

\[
 2g_X-2=4(-2)+6=-2,
\]

so the source curve has genus zero.  The first branch divisor has exactly
two index-two primes, both selected and boundary-colored.

Attach the selected curve

\[
 \mathbb P^1\setminus\{0,1,\infty\}.
\]

Its unit lattice has rank two.  This is the same logarithmic curve type as
the residual \(L=0\) normalization in the current \(A_4\) frontier.

### Fano degree seven

In the point action on the seven nonzero vectors of \(\mathbb F_2^3\), use

\[
\begin{aligned}
 \tau_0&=(14)(36),\\
 \tau_1&=(125)(374),\\
 \tau_\infty&=(1736452).
\end{aligned}
\]

They have product one, generate a group of order \(168\), and have
contributions \(2,4,6\).  Thus

\[
 2g_X-2=7(-2)+12=-2.
\]

The involution branch contains two index-two primes and three unramified
primes.  This is the smallest exact benchmark for combining two selected
ramified components with the desired
\(\operatorname{GL}_3(\mathbb F_2)\) monodromy.

## 4. Conductor and elliptic benchmarks

### Nodal rational curve

Let the normalization be \(\mathbb G_m\), and identify the points
\(1\) and \(-1\).  The affine curve ring is the equalizer

\[
 \{f\in k[t,t^{-1}]:f(1)=f(-1)\}.
\]

A unit \(ct^m\) descends exactly when \(m\) is even.  The normalized unit
lattice \(\mathbb Z\) is therefore replaced by \(2\mathbb Z\): rank one,
index two.  The node has \(\delta=1\), so the curve has arithmetic genus one.

The positive fixture uses an involution swapping the two node branches.  The
negative fixture sends their pair to a different pair and is rejected.  This
is the abstract form of the obstruction isolated in
[`DAVENPORT_BOUNDARY_INVOLUTION.md`](DAVENPORT_BOUNDARY_INVOLUTION.md).

### Elliptic selected boundary

The elliptic fixture is a degree-three cover of \(\mathbb P^1\) with six
simple branch values.  Its branch permutations are paired transpositions
which generate \(S_3\) and have product one.  Riemann--Hurwitz gives

\[
 2g_X-2=3(-2)+6=0,
\]

so \(g_X=1\).

The selected affine curve is an elliptic curve with one puncture.  Its
puncture-supported degree-zero divisor lattice is zero, hence it has only
scalar units.  A plane-cubic adjunction ledger gives

\[
 3(3-3)=0=2g-2.
\]

This proves only that the proposed elliptic package survives these coarse
filters.  Constructing an affine-space Keller realization remains open.

## 5. First symbolic replay: the \(A_4\) core

The first stage-two adapter is
[`verify_boundary_package_a4_stage_two.py`](../scripts/verify_boundary_package_a4_stage_two.py).
It does not copy the \(A_4\) formulas.  It replays and consumes the exact
globals from:

- [`verify_a4_affine_keller_frontier.py`](../scripts/verify_a4_affine_keller_frontier.py);
- [`verify_a4_ledger_reduction.py`](../scripts/verify_a4_ledger_reduction.py);
- [`verify_a4_pure_target_ledger.py`](../scripts/verify_a4_pure_target_ledger.py);
  and
- [`verify_a4_keller_inverse_cover.py`](../scripts/verify_a4_keller_inverse_cover.py).

SHA-256 hashes of all four replayed scripts are emitted in the certificate.
The adapter proves the following matches.

1. The oriented inverse equation has degree four and discriminant
   \(4096\delta\), with \(D^2=\delta\) on the oriented target.
2. The specialization
   \[
   T^4-7T^2-3T+1
   \]
   has discriminant \(183^2\), is irreducible, and has irreducible cubic
   resolvent.  Together with the square generic discriminant, this certifies
   the exact generic \(A_4\) group rather than only a subgroup of \(S_4\).
3. Exact polynomial division extracts the cone determinant orders
   \[
   (v_W,v_K,v_L)(\det D\Phi)=(2,3,1).
   \]
4. The fourfold lift and target pullback both have orders
   \[
   (3,3,2),
   \]
   so the auxiliary contribution is \((1,0,1)\).  These three extracted rows
   agree exactly with the abstract package ledger.
5. The pure-target lift satisfies
   \[
   \det D\widehat\Phi=\mathcal B(\widehat\Phi_1,
   \widehat\Phi_2,\widehat\Phi_3),
   \]
   while the original cone determinant is nonconstant and the naive
   multiplicative target factorization has denominator
   \(4W^2K^3L\).
6. The normalization parameter for \(L=0\) gives a rank-two restricted cone
   map at the exact point \(t=W=1\).  Since the target cubic \(\mathcal B\)
   is irreducible, \(L\) dominates its generic point.  The pullback order
   two and Jacobian order one therefore make \(L\) one affine ramification
   prime with \(e=2\).
7. The divisors \(W\) and \(K\) do not dominate \(\mathcal B=0\):
   \(W=0\) maps to the origin, while \(K=0\) forces \(P=0\) and maps into
   the codimension-two locus \(P=\mathcal B=0\).
8. If
   \[
   s=\frac{P-Q}{R}
   \]
   is the line-slope parameter on the normalization of the target cubic,
   then restriction to the normalization \(t=U-V\) of \(L\) gives
   \[
   s=\frac{t^2}{2t-3}.
   \]
   Thus the residue extension is generated by
   \[
   t^2-2st+3s=0.
   \]
   Its discriminant is \(4s(s-3)\), which has odd valuations at \(s=0\)
   and \(s=3\), so the polynomial is irreducible over
   \(\mathbb Q(s)\).  Consequently \(f(L/\mathcal B)=2\).  Over the exact
   splitting extension
   \[
   r^2=s(s-3)
   \]
   it factors as
   \[
   (t-s-r)(t-s+r),
   \]
   explicitly exhibiting the two geometric branches.
9. The projective source and target cubics have node
   \([1:1:0]\), whose two normalization branches are \(0\) and \(3\).
   The residue map sends
   \[
   0\longmapsto0,\qquad 3\longmapsto3,\qquad
   \infty\longmapsto\infty,
   \]
   so it preserves the conductor pair and the three-puncture set.

The replay therefore upgrades the \(A_4\) package from a purely discrete
fixture to an exact **symbolic-core match**.  Its overall status remains
`unknown`, but the selected generic double-transposition profile can now be
sharpened.

The residue calculation corrects the coarser degree count in the first
version of this replay.  The affine \(L\)-prime contributes

\[
 e(L/\mathcal B)f(L/\mathcal B)=2\cdot2=4.
\]

It therefore exhausts the degree-four extension over the generic target
branch: there is no second missing global height-one prime.  After extending
the residue field by \(\sqrt{s(s-3)}\), the quadratic residue polynomial
splits and \(L\) produces two geometric \(e=2\) primes.  These give the two
cycles of the \(A_4\) double transposition.  Hence the selected
double-transposition **geometric cycle profile is matched**, but the
package's global prime ledger is not: it requests two \((e,f)=(2,1)\)
primes, whereas this architecture has one \((2,2)\) prime.  This does not
yet match the other two branch generators in the abstract \((2,3,3)\)
triangle fixture to divisors of the higher-dimensional cone.

The coloring is also not matched: the package asks that both geometric
index-two primes lie on the normalization boundary, as a Keller realization
requires, whereas both arise from \(L=0\) inside the affine cone source.
The construction problem must split or replace this global residue-degree-two
prime and move the resulting ramification outside the affine source without
changing the \(A_4\) field.

Independently, a logarithmic target-Jacobian identity is not an ordinary
constant-Jacobian factorization through affine space.  The existing
denominator obstruction proves that the obvious target factorization fails.

This distinction is machine-recorded as

```text
package_ledger_matches = true
pure_target_log_keller = true
ordinary_cone_keller = false
dominant_affine_ramified_prime = "L"
ramification_index = 2
residue_degree = 2
local_degree_contribution = 4
remaining_degree_over_target_branch = 0
geometric_ramified_prime_count = 2
conductor_pair_preserved = true
contracted_source_divisors = ["W", "K"]
selected_double_transposition_profile_matched = true
abstract_global_prime_profile_matched = false
abstract_boundary_coloring_matched = false
affine_space_realization_matched = false
```

Thus the adapter advances the synthesis without treating a shared
monodromy group or a log-Keller identity as a boundary-package realization.

## 6. What stage one still lacks

The prototype now implements finite pole boxes and exact membership in a
finitely generated actual semigroup.  The next layer must make the remaining
certificate references from
[`FINITE_VALUATION_ALGEBRAIZATION.md`](../cancellation/FINITE_VALUATION_ALGEBRAIZATION.md):

1. derive denominator/polar completeness from supplied equations;
2. derive the pole bounds rather than trust the package;
3. compute generators and relevant holes of the actual value semigroup;
4. replay a finite Khovanskii/SAGBI presentation; and
5. verify Rees-module strictness by saturation.

It must also distinguish an actual boundary-class or normal-core valuation
matrix computed from a normalization from a synthetically supplied
candidate.  The compiler now refuses to treat the core matrix as conclusive
without named normality, complete-codimension-one-boundary data, and either
core factoriality or a class presentation with lifted relation witnesses; it
does not replay those certificates.  The benchmark
matrices in the current script are inputs to the necessary tests, not proofs
that the corresponding normalizations exist.

After that, the first symbolic target should be the \(A_4\) package.  A
first adapter now recovers the known root equation and residual \(WL\)
ledger, recognizes the pure-target log-Keller lift, and computes the complete
generic branch datum as one affine \((e,f)=(2,2)\) \(L\)-prime whose geometric
splitting gives the double transposition.  The next \(A_4\) task is to split
or replace that global prime with two residue-degree-one boundary primes and
solve the boundary-coloring mismatch.
The second monodromy target is the Fano package from
[`ABSOLUTE_SUNADA_KELLER_RESEARCH.md`](ABSOLUTE_SUNADA_KELLER_RESEARCH.md).

## 7. Reproduction

Run

```bash
python3 scripts/verify_boundary_package_compiler.py
python3 scripts/verify_boundary_package_a4_stage_two.py
```

The verifier checks four surviving `unknown` packages and eight intrinsic
`obstructed` packages, including the factorial-core, presented-core, and
retained-root Euler fixtures, then
emits their complete machine-readable reports.
The second command replays the four existing \(A_4\) symbolic certificates,
extracts their divisor orders, computes the residue degree and conductor map
over the target cubic, compares them to the abstract package, and emits the
remaining realization gates.
