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

The current `BoundaryPackage` has five proof-relevant blocks.

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

For a normal finite normalization \(\bar X\) with distinguished open
\(U=\mathbb A^n\), the localization sequence forces the boundary-class map

\[
 \bigoplus_i\mathbb Z[E_i]\longrightarrow\operatorname{Cl}(\bar X)
\]

to have trivial kernel and cokernel.  The prototype checks a supplied square
boundary-class matrix has determinant \(\pm1\), and rejects declared
nonzero source unit or class-group rank.

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

It must also distinguish an actual boundary-class matrix computed from a
normalization from a synthetically supplied unimodular candidate.  The
benchmark matrices in the current script are inputs to the necessary test,
not proofs that the corresponding normalizations exist.

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

The verifier checks four surviving `unknown` packages and five intrinsic
`obstructed` packages, then emits their complete machine-readable reports.
The second command replays the four existing \(A_4\) symbolic certificates,
extracts their divisor orders, computes the residue degree and conductor map
over the target cubic, compares them to the abstract package, and emits the
remaining realization gates.
